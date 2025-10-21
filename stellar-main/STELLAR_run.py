import argparse
import os
from pathlib import Path

import numpy as np
import torch
from datasets import GraphDataset, load_hubmap_data, load_tonsilbe_data
from STELLAR import STELLAR
from utils import prepare_save_dir

# TODO: generalize as appropriate with new or multiple references
#   and move the functionality to find this
data_filename = "20251007_cleaned_trainingdata_yang.h5ad"

data_dir_possibilities = [
    Path("/data"),
    Path(__file__).parent / "data",
]


def find_data_file() -> Path:
    for path in data_dir_possibilities:
        if (f := path / data_filename).is_file():
            print("Found training data file at", f)
            return f
    message_pieces = [f"Couldn't find data directory; tried:"]
    message_pieces.extend([f"\t{path}" for path in data_dir_possibilities])
    raise FileNotFoundError("\n".join(message_pieces))


def main():
    parser = argparse.ArgumentParser(description="STELLAR")
    parser.add_argument("cell_data_h5ad", type=Path)

    parser.add_argument(
        "--seed", type=int, default=1, metavar="S", help="random seed (default: 1),"
    )
    parser.add_argument(
        "--name",
        type=str,
        default="STELLAR",
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=20,
    )
    parser.add_argument(
        "--lr",
        type=float,
        default=1e-3,
    )
    parser.add_argument(
        "--wd",
        type=float,
        default=5e-2,
    )
    parser.add_argument(
        "--num-heads",
        type=int,
        default=22,
    )
    parser.add_argument(
        "--num-seed-class",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--sample-rate",
        type=float,
        default=0.5,
    )
    parser.add_argument(
        "-b",
        "--batch-size",
        default=1,
        type=int,
        metavar="N",
        help="mini-batch size",
    )
    parser.add_argument(
        "--distance_thres",
        default=50,
        type=int,
    )
    parser.add_argument(
        "--savedir",
        type=Path,
        default=Path(),
    )

    args = parser.parse_args()
    args.cuda = torch.cuda.is_available()
    args.device = torch.device("cuda" if args.cuda else "cpu")

    # Seed the run and create saving directory
    args.name = "STELLAR"
    args = prepare_save_dir(args, __file__)

    (
        labeled_X,
        labeled_y,
        unlabeled_X,
        labeled_edges,
        unlabeled_edges,
        inverse_dict,
        unlabeled_cell_indexes,
    ) = load_hubmap_data(
        find_data_file(),
        args.cell_data_h5ad,
        args.distance_thres,
        args.sample_rate,
    )
    dataset = GraphDataset(
        labeled_X,
        labeled_y,
        unlabeled_X,
        labeled_edges,
        unlabeled_edges,
    )

    stellar = STELLAR(args, dataset)
    stellar.train()
    _, results = stellar.pred()

    out_dir = Path("stellar")
    out_dir.mkdir(exist_ok=True, parents=True)
    with open(out_dir / f"{args.cell_data_h5ad.stem}.csv", "w") as f:
        print("ID,STELLAR_CellType", file=f)
        for cell_id, cell_type_id in zip(unlabeled_cell_indexes, results):
            cell_type = inverse_dict[cell_type_id]
            print(f"{cell_id},{cell_type}", file=f)

    print("done")


if __name__ == "__main__":
    main()
