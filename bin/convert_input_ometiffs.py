#!/usr/bin/env python3
import importlib.resources
from argparse import ArgumentParser
from pathlib import Path

import anndata
import pandas as pd
from ome_utils import find_ome_tiffs
from sklearn.preprocessing import StandardScaler
from sprm import modules
from xarray import DataArray


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
        obs=pd.DataFrame(index=mean_expr.coords["cell_index"]),
        var=pd.DataFrame(index=mean_expr.coords["expr_channel"]),
    )
    # So, create the DataFrame after the AnnData, using .obs_names as the
    # index, to make sure everything matches with minimal effort.
    cell_centers_df = pd.DataFrame(
        core.cell_centers[core.mask.interior_cells],
        index=image_adata.obs_names,
        columns=["X", "Y", "Z"],
    ).loc[:, ["Y", "X"]]
    image_adata.obsm["X_spatial"] = cell_centers_df
    image_adata.obs["unique_region"] = expr.stem
    return image_adata


def main(directory: Path):
    expr_dir, mask_dir = find_expr_mask_dir(directory)
    exprs = sorted(find_ome_tiffs(expr_dir))
    masks = sorted(find_ome_tiffs(mask_dir))

    adatas = []
    for expr, mask in zip(exprs, masks):
        adatas.append(convert(expr, mask))

    adata = anndata.concat(adatas)
    adata.write_h5ad("cell_data.h5ad")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("directory", type=Path)
    args = p.parse_args()

    main(args.directory)
