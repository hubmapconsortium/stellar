#!/usr/bin/env python3
import importlib.resources
from argparse import ArgumentParser
from pathlib import Path

import anndata
import numpy as np
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
        compute_texture=False,  # Set to True if you want texture features too
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
    spatial_coords = np.zeros((len(core.mask.cell_index), 2))

    expr_array = mean_expr.loc["cell", :, :].to_numpy()
    scaled_expr_array = StandardScaler().fit_transform(expr_array)
    image_adata = anndata.AnnData(
        X=scaled_expr_array,
        obs=pd.DataFrame(index=mean_expr.coords["cell_index"]),
        var=pd.DataFrame(index=mean_expr.coords["expr_channel"]),
        obsm={"X_spatial": spatial_coords},
    )
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
