
import os
import sys
from argparse import ArgumentParser
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch
from bioio import BioImage
from datasets import GraphDataset, load_hubmap_data, load_tonsilbe_data
from ome_utils import find_ome_tiffs
from STELLAR import STELLAR
from utils import prepare_save_dir

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


def standardize_antb_names(antibodies):
    """
    Helper function to standardize antibody names.
    """
    for name in antibodies:
        new_name = find_antibody_key(name)
        antibodies = [new_name if x==name else x for x in antibodies]
    return antibodies


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

def find_data_file(tissue) -> Path:
    for path in data_dir_possibilities:
        if (f := path / adata_paths[tissue]).is_file():
            print("Found training data file at", f)
            return f
    message_pieces = [f"Couldn't find data directory; tried:"]
    message_pieces.extend([f"\t{path}" for path in data_dir_possibilities])
    raise FileNotFoundError("\n".join(message_pieces))


def find_expr_mask_dir(base_dir: Path) -> tuple[Path, Path]:
    if (d := base_dir / "pipeline_output").is_dir():
        return d / "expr", d / "mask"
    if (d := base_dir / "stitched").is_dir():
        return d / "expressions", d / "mask"
    raise ValueError("Couldn't find image and mask directories")


def write_result(val):
    with open("results.txt", "w") as f:
        f.write(val)


def main(directory, tissue):
    # TODO: Check if the dataset is from an approved provider
    # if provider not in possible_providers:
    #     write_result("False")
    #     return
    # Check if the there is a model that matches the tissue for the dataset
    if tissue not in adata_paths.keys():
        print(f"There is no STELLAR model for {tissue}.")
        write_result("False")
        return
    # Open provider image and then training data to check for marker names
    expr_dir, mask_dir = find_expr_mask_dir(directory)
    exprs = sorted(find_ome_tiffs(expr_dir))
    provider_imgs = [BioImage(expr) for expr in exprs]
    provider_ch_names = [str(i) for i in provider_imgs[0].channel_names]
    print("Provider channel names before standardizing:")
    print(provider_ch_names)
    print("Provider channel names after standardizing:")
    standardize_antb_names(provider_ch_names)
    print(provider_ch_names)
    adata = ad.read_h5ad(adata_paths[tissue])
    train_ch_names = adata.var_names.to_list()
    print("Training channel names:")
    print(train_ch_names)
    # Markers must all match and be in the same orer
    common_vars = [v for v in train_ch_names if v in provider_ch_names]
    print("Common variables (Training Order):", common_vars)
    # Make sure the channel names are identical
    if common_vars != train_ch_names:
        print("Provider channel names do not include the following required channels to run STELLAR:")
        print([v for v in train_ch_names if v not in provider_ch_names])
        write_result("False")
        return
    # Everything passed if we get here, pass this to the rest of the pipeline
    write_result(f"passed")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("tissue", type=Path)
    p.add_argument("directory", type=Path)
    args = p.parse_args()

    main(args.directory, args.tissue)