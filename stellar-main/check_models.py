
import os
import sys
from argparse import ArgumentParser
from pathlib import Path

import anndata as ad
import numpy as np
import torch
from bioio import BioImage
from datasets import GraphDataset, load_hubmap_data, load_tonsilbe_data
from ome_utils import find_ome_tiffs
from STELLAR import STELLAR
from utils import prepare_save_dir

model_paths = {
    "intestine": Path("models/intestine_placeholder.pt"),
    # other tissues : other paths,
}

possible_providers = []


def find_expr_mask_dir(base_dir: Path) -> tuple[Path, Path]:
    if (d := base_dir / "pipeline_output").is_dir():
        return d / "expr", d / "mask"
    if (d := base_dir / "stitched").is_dir():
        return d / "expressions", d / "mask"
    raise ValueError("Couldn't find image and mask directories")


def write_result(val):
    with open("results.txt", "w") as f:
        f.write(val)


def main(directory, tissue, provider):
    # Check if the dataset is from an approved provider
    if provider not in possible_providers:
        write_result("False")
        return
    # Check if the there is a model that matches the tissue for the dataset
    if tissue not in model_paths.keys():
        write_result("False")
        return
    # Open provider image and then model to check for marker names
    # provider will be ometiff
    expr_dir, mask_dir = find_expr_mask_dir(directory)
    exprs = sorted(find_ome_tiffs(expr_dir))
    provider_imgs = [BioImage(expr) for expr in exprs]
    provider_ch_names = [str(i) for i in provider_imgs[0].channel_names]
    print("Provider channel names:")
    print(provider_ch_names)
    # model will be h5ad
    model = ad.read_h5ad(model_paths[tissue])
    model_ch_names =  model.var_keys()
    print("Model channel names:")
    print(model_ch_names)
    # Make sure the channel names are identical
    if provider_ch_names != model_ch_names:
        write_result("False")
        return
    # Everything passed if we get here, pass this to the rest of the pipeline
    write_result("True")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("tissue", type=Path)
    p.add_argument("provider", type=Path)
    p.add_argument("directory", type=Path)
    args = p.parse_args()

    main(args.directory, args.tissue, args.provider)