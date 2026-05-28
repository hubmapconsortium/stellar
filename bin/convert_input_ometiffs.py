#!/usr/bin/env python3
import importlib.resources
from argparse import ArgumentParser
from pathlib import Path
from math import ceil, log2
import anndata
import pandas as pd
import spatialdata as sd
import numpy as np
from ome_utils import find_ome_tiffs
from sklearn.preprocessing import StandardScaler
from sprm import modules
from xarray import DataArray
from spatialdata.models import Image2DModel, Labels2DModel, PointsModel, TableModel
import tracemalloc
from bioio import BioImage

desired_pixel_size_for_pyramid=250

data_dir_possibilities = [
    Path("/data"),
    Path(__file__).parent / "data",
]

adata_paths = {
    "intestine": Path("20260107_newSPRM_64CODEX_SLI_annotated.h5ad"),
    # other tissues : other paths,
}

possible_providers = [] # TODO: put providers here

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


def find_data_file(tissue):
    for path in data_dir_possibilities:
        if (f := path / adata_paths[tissue]).is_file():
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


def write_pseudo_adata():
    with open("results.txt", "w") as f:
        X = np.random.rand(3, 5)
        obs_df = pd.DataFrame(
            index=[f"cell_{i}" for i in range(3)],
            data = {"fake": np.random.randint(0,3, size=3), "stage": "mock"}
        )
        var_df = pd.DataFrame(
            index=[f"gene_{j}" for j in range(5)]
        )
        fake_adata = anndata.AnnData(X=X, obs=obs_df, var=var_df)
        fake_adata.write_h5ad("no_model.h5ad")


def find_expr_mask_dir(base_dir: Path) -> tuple[Path, Path]:
    if (d := base_dir / "pipeline_output").is_dir():
        return d / "expr", d / "mask"
    if (d := base_dir / "stitched").is_dir():
        return d / "expressions", d / "mask"
    raise ValueError("Couldn't find image and mask directories")


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
    expr_adata = image_adata.copy()
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

    print('Starting SpatialData conversion')
    print("Loading image data")
    image = BioImage(expr)
    image_data_squeezed = image.data.squeeze()
    print("... done. Original shape:", image.data.shape)

    image_scale_factors = (2,) * ceil(
        log2(max(image_data_squeezed.shape[1:]) / desired_pixel_size_for_pyramid)
    )

    img_for_sdata = Image2DModel.parse(
        data=image_data_squeezed,
        c_coords=image.channel_names,
        scale_factors=image_scale_factors,
    )

    print("Loading mask data")
    mask = BioImage(mask)
    mask_data = {ch: mask.data[0, i, 0, :, :] for i, ch in enumerate(mask.channel_names)}
    cell_indexes = sorted(set(mask_data["cells"].flat) - {0})
    print("... done.", len(cell_indexes), "cells")

    masks_for_sdata = {
        ch: Labels2DModel.parse(
            data=mask_array,
            scale_factors=image_scale_factors,
        )
        for ch, mask_array in mask_data.items()
    }

    cell_centers_copy = cell_centers_df.copy()
    cell_centers_copy.columns = ['y', 'x']
    shapes_for_sdata = PointsModel.parse(cell_centers_copy)

    tables = {
        "mean_expr": TableModel.parse(expr_adata),
    }

    for table in tables.values():
        table.obs["cell_id"] = pd.Series([int(i.split('-')[1]) for i in table.obs.index.values], index=table.obs.index)
        table.obs["region"] = pd.Categorical(["cells"] * len(table))
        table.uns["spatialdata_attrs"] = {
            "region": "cells",
            "region_key": "region",
            "instance_key": "cell_id",
        }

    sdata = sd.SpatialData(
        images={"expression": img_for_sdata},
        points={"centers": shapes_for_sdata},
        labels=masks_for_sdata,
        tables=tables,
    )

    sdata_name = f"{csv_base}_spatialdata.zarr"
    print("Saving SpatialData object to", sdata_name)
    print(sdata)
    sdata.write(sdata_name, overwrite=True)

    return image_adata


def main(directory: Path, tissue: str):
    # Check if a model is available before opening the image
    # TODO: check for provider as well when I receive that info
    adata_path = find_data_file(tissue)
    if not adata_path:
        print(f"There is no STELLAR model for {tissue}.")
        write_pseudo_adata()
        return

    train_adata = anndata.read_h5ad(adata_path)
    tracemalloc.start()
    expr_dir, mask_dir = find_expr_mask_dir(directory)
    exprs = sorted(find_ome_tiffs(expr_dir))
    masks = sorted(find_ome_tiffs(mask_dir))

    adatas = []
    for expr, mask in zip(exprs, masks):
        adatas.append(convert(expr, mask))

    adata = anndata.concat(adatas, index_unique="-")
    # Check if antibody names match
    test_var = standardize_antb_df(adata.var)
    print("Training data variables:", train_adata.var_names)
    print("Test data variables before standardizing:", adata.var_names)
    print("Test data variables after standardizing:", test_var)
    adata.var = test_var
    # Markers must all match and be in the same order
    common_vars = [v for v in train_adata.var_names if v in adata.var_names]
    print("Common variables (Training Order):", common_vars)
    test_adata = adata[:, common_vars].copy()
    if train_adata.var_names.to_list() != test_adata.var_names.to_list():
        missing_vars = train_adata.var_names.to_list().difference(test_adata.var_names.to_list())
        print("The following variables are missing from the test data:", missing_vars)
        write_pseudo_adata()
        return

    test_adata.write_h5ad("cell_data.h5ad")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("directory", type=Path)
    p.add_argument("tissue", type=str)
    args = p.parse_args()

    main(args.directory, args.tissue)
