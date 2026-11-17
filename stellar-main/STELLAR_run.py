import argparse
from pathlib import Path
from typing import Optional

import torch
from datasets import GraphDataset, load_hubmap_data
from STELLAR import STELLAR
from utils import prepare_save_dir

# TODO: generalize as appropriate with new or multiple references
#   and move the functionality to find this
data_filename = "20241108_HuBMAP_intestine_annotated.h5ad"

data_dir_possibilities = [
    Path("/data"),
    Path(__file__).parent / "data",
]


def find_data_file(tissue: str) -> Optional[Path]:
    for path in data_dir_possibilities:
        if (d := path / tissue).is_dir():
            h5ads = list(d.glob("*.h5ad"))
            assert (l := len(h5ads)) == 1, f"Need 1 H5AD file in {d}; found {l}"
            print("Found training data file at", h5ads[0])
            return h5ads[0]
    message_pieces = [f"Couldn't find data directory; tried:"]
    message_pieces.extend([f"\t{path}" for path in data_dir_possibilities])
    print("\n".join(message_pieces))
    return None


def main():
    parser = argparse.ArgumentParser(description="STELLAR")
    parser.add_argument("cell_data_h5ad", type=Path)
    parser.add_argument("model_dir_name")

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

    out_dir = Path("stellar")
    out_dir.mkdir(exist_ok=True, parents=True)

    maybe_training_data_path = find_data_file(args.model_dir_name)
    if maybe_training_data_path is None:
        # No training data for this tissue. Nothing to do; leave empty directory.
        return

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
        maybe_training_data_path,
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

    with open(out_dir / f"{args.cell_data_h5ad.stem}.csv", "w") as f:
        print("ID,STELLAR_CellType", file=f)
        for cell_id, cell_type_id in zip(unlabeled_cell_indexes, results):
            cell_type = inverse_dict[cell_type_id]
            cell_id_adj = cell_id.rsplit("-", 2)[1]
            print(f"{cell_id_adj},{cell_type}", file=f)

    print("done")


if __name__ == "__main__":
    main()
