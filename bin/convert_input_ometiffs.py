#!/usr/bin/env python3
import importlib.resources
from argparse import ArgumentParser
from pathlib import Path
from math import ceil, log2
import anndata
import pandas as pd
import numpy as np
from ome_utils import find_ome_tiffs
from sklearn.preprocessing import StandardScaler
from sprm import modules
from xarray import DataArray
import tracemalloc
from bioio import BioImage

desired_pixel_size_for_pyramid=250

data_dir_possibilities = [
    Path("/data"),
    Path(__file__).parent / "data",
]

adata_paths = {
    "Stanford TMC": {"LI": Path("20260107_newSPRM_64CODEX_SLI_annotated.h5ad")},
    "General Electric RTI": {"SK": Path("20260324_skin_v5_12data_leiden15_noCD3_refine_labeled_with_spatial.h5ad")},
    # other tissues : other paths,
}

antibodies_dict = {
    "BCL2": "BCL-2",
    "CollagenIV": ["CollIV", "Collagen IV", "collagen IV", "COLIV"],
    "Cytokeratin": "cytokeratin",
    "eCAD": ["E-CAD", "ECAD"],
    "HLA-DR": "HLADR",
    "Hoechst1": "HOECHST1",
    "PanCK": "panCK",
    "Podoplanin": ["Podoplan", "podoplanin", "PDPN"],
    "Synaptophysin": ["Synapt", "Synapto"],
    "aDefensin5": ["aDef5", "aDefensin 5"],
    "MUC1": ["MUC1/EMA", "MUC-1"],
    "NKG2D (CD314)": ["NKG2D", "NKG2G"],
    "aSMA": ["SMActin", "a-SMA", "SMA"],
    "MUC2": "MUC-2",
    "Foxp3": "FoxP3",
}


def standardize_antb_df(antibodies_df: pd.DataFrame) -> pd.DataFrame:
    """
    Helper function to standardize antibody names.
    """
    for idx, row in antibodies_df.iterrows():
        new_name = find_antibody_key(idx)
        antibodies_df = antibodies_df.rename(index={idx: new_name})
    return antibodies_df


def find_antibody_key(value: str) -> str:
    """
    Helper function to standardize antibody names.
    """
    value_lower = value.strip().lower()
    for key, val in antibodies_dict.items():
        if isinstance(val, str) and val.strip().lower() == value_lower:
            return key
        elif isinstance(val, list) and value_lower in [v.strip().lower() for v in val]:
            return key
    return value


def find_data_file(tissue, provider):
    for path in data_dir_possibilities:
        if (f := path / adata_paths[provider][tissue]).is_file():
            print("Found training data file at", f)
            return f
    message_pieces = [f"Couldn't find data directory; tried:"]
    message_pieces.extend([f"\t{path}" for path in data_dir_possibilities])
    return None


def check_tissue(tissue):
    print(adata_paths.keys())
    if tissue in adata_paths.keys():
        print("found tissue adata path")
        return True
    else:
        return False


def find_expr_mask_dir(base_dir: Path) -> tuple[Path, Path]:
    if (d := base_dir / "pipeline_output").is_dir():
        return d / "expr", d / "mask"
    if (d := base_dir / "stitched").is_dir():
        return d / "expressions", d / "mask"
    raise ValueError("Couldn't find image and mask directories")


def align_skin_vars(test_adata, train_adata):
    # standardization commented out for skin
    test_lowercase_map = {v.lower(): v for v in test_adata.var_names}
    test_var_lowercase = list(test_lowercase_map.keys())

    print("Training data variables:", train_adata.var_names)
    train_var_lowercase = [v.lower() for v in train_adata.var_names]
    print("Test data variables before standardizing:", test_adata.var_names)
    common_vars_lowercase = [v for v in train_var_lowercase if v in test_var_lowercase]
    common_vars_original = [test_lowercase_map[v] for v in common_vars_lowercase]
    print("Common variables:", common_vars_original)
    # test_adata.var_names = list(test_adata.var_names)
    test_adata = test_adata[:, common_vars_original].copy()
    train_vars_lower_list = [v.lower() for v in train_adata.var_names]
    test_vars_lower_list = [v.lower() for v in test_adata.var_names]

    return test_adata, train_vars_lower_list, test_vars_lower_list


options_file = importlib.resources.files("sprm") / "options.txt"
output_dir_base = Path("sprm_features_outputs")


def convert(expr: Path, mask: Path):
    output_dir = output_dir_base / expr.stem
    output_dir.mkdir(exist_ok=True, parents=True)

    csv_base = expr.name.split('.', 1)[0]

    core = modules.preprocessing.run(
        img_file=expr,
        mask_file=mask,
        output_dir=output_dir,
        options=options_file,
    )
    features = modules.cell_features.run(
        core_data=core,
        output_dir=output_dir,
        compute_texture=False,
    )

    mean_expr = DataArray(
        features.mean_vector.squeeze(),
        coords=[
            core.mask.channel_labels,
            core.mask.cell_index,
            core.im.channel_labels,
        ],
        dims=["mask_channel", "cell_index", "expr_channel"],
    )

    expr_array = mean_expr.loc["cell", :, :].to_numpy()
    scaled_expr_array = StandardScaler().fit_transform(expr_array)
    # Don't assign obsm={"X_spatial": ...} here, since we want this stored
    # as a DataFrame with Y, X columns, but the index of such a DataFrame
    # and the index of the overall AnnData must match when instantiating
    # in that way.
    image_adata = anndata.AnnData(
        X=scaled_expr_array,
        obs=pd.DataFrame(index = [f"{csv_base}-"+str(i) for i in mean_expr.coords["cell_index"].to_series().tolist()]),
        var=pd.DataFrame(index=mean_expr.coords["expr_channel"]),
    )

    image_adata.var_names_make_unique()
    # So, create the DataFrame after the AnnData, using .obs_names as the
    # index, to make sure everything matches with minimal effort.
    cell_centers_df = pd.DataFrame(
        core.cell_centers[core.mask.interior_cells],
        index=image_adata.obs_names,
        columns=["X", "Y", "Z"],
    ).loc[:, ["Y", "X"]]
    image_adata.obsm["X_spatial"] = cell_centers_df
    image_adata.obs["unique_region"] = expr.stem

    print('SPRM conversion complete')

    return image_adata


def main(directory: Path, tissue: str, provider: str):
    tracemalloc.start()
    expr_dir, mask_dir = find_expr_mask_dir(directory)
    exprs = sorted(find_ome_tiffs(expr_dir))
    masks = sorted(find_ome_tiffs(mask_dir))

    adatas = []
    for expr, mask in zip(exprs, masks):
        adatas.append(convert(expr, mask))

    adata = anndata.concat(adatas, index_unique="-")
    # Check for model
    train_adata_path = find_data_file(tissue, provider)
    if not train_adata_path:
        print(f"There is no STELLAR model for {tissue}.")
        adata.write_h5ad("cell_data.h5ad")
        return
    train_adata = anndata.read_h5ad(train_adata_path)
    # standardize antibody names
    test_adata, train_vars, test_vars = align_skin_vars(adata, train_adata)
    if train_vars != test_vars:
        missing_vars = set(train_vars) - set(test_vars)
        print("The following variables are missing from the test data:", missing_vars)
        print("Exiting program, STELLAR will not run.")
    else:
        print("All required variables are present in the test data. Proceeding to STELLAR.")
    test_adata.write_h5ad("cell_data.h5ad")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("directory", type=Path)
    p.add_argument("tissue", type=str)
    p.add_argument("provider")
    args = p.parse_args()

    main(args.directory, args.tissue, args.provider)
