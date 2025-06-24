#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from typing import Iterable

import anndata
import pandas as pd
from aicsimageio import AICSImage
from ome_utils import find_ome_tiffs
from skimage.measure import regionprops
from sklearn.preprocessing import StandardScaler


def get_ome_tiff_paths(input_dir: Path) -> Iterable[tuple[Path, Path]]:
    """
    Yields 2-tuples:
     [0] full Path to source file
     [1] output file Path (source file relative to input_dir)
    """
    for ome_tiff in find_ome_tiffs(input_dir):
        yield ome_tiff, ome_tiff.relative_to(input_dir)


def convert_image_data(expr_image_file: Path, mask_image_file: Path) -> anndata.AnnData:
    print("Loading image data")
    image = AICSImage(expr_image_file)
    image_data_squeezed = image.data.squeeze()
    print("... done. Original shape:", image.data.shape)

    print("Loading mask data")
    mask = AICSImage(mask_image_file)
    mask_data = {
        ch: mask.data[0, i, 0, :, :] for i, ch in enumerate(mask.channel_names)
    }
    cell_indexes = sorted(set(mask_data["cells"].flat) - {0})
    print("... done.", len(cell_indexes), "cells")

    total_expr = pd.DataFrame(0.0, index=cell_indexes, columns=image.channel_names)
    mean_expr = pd.DataFrame(0.0, index=cell_indexes, columns=image.channel_names)
    centers = pd.DataFrame(0.0, index=cell_indexes, columns=["y", "x"])

    for rp in regionprops(mask_data["cells"]):
        centers.loc[rp.label, :] = rp.centroid
        pixel_data = image_data_squeezed[:, *(rp.coords.T)]
        total_vec = pixel_data.sum(axis=1)
        total_expr.loc[rp.label, :] = total_vec
        mean_expr.loc[rp.label, :] = total_vec / rp.area

    scaled_expr_array = StandardScaler().fit_transform(mean_expr)
    image_adata = anndata.AnnData(
        X=scaled_expr_array,
        obs=pd.DataFrame(index=mean_expr.index),
        var=pd.DataFrame(index=mean_expr.columns),
        obsm={"X_spatial": centers.values},
    )

    image_adata.obs["unique_region"] = expr_image_file.stem
    print(image_adata)

    return image_adata


def find_expr_mask_dir(base_dir: Path) -> tuple[Path, Path]:
    if (d := base_dir / "pipeline_output").is_dir():
        return d / "expr", d / "mask"
    if (d := base_dir / "stitched").is_dir():
        return d / "expressions", d / "mask"
    raise ValueError("Couldn't find image and mask directories")


def main(directory: Path):
    expr_dir, mask_dir = find_expr_mask_dir(directory)
    exprs = sorted(find_ome_tiffs(expr_dir))
    masks = sorted(find_ome_tiffs(mask_dir))

    adatas = []
    for expr, mask in zip(exprs, masks):
        adatas.append(convert_image_data(expr, mask))

    adata = anndata.concat(adatas)
    adata.write_h5ad("cell_data.h5ad")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("directory", type=Path)
    args = p.parse_args()

    main(args.directory)
