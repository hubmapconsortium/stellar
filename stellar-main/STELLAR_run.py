import argparse
import os
import sys
from pathlib import Path

import anndata as ad
import json
import numpy as np
import pandas as pd
import torch
from datasets import GraphDataset, load_hubmap_data, load_tonsilbe_data
from STELLAR import STELLAR
from utils import prepare_save_dir

# TODO: generalize as appropriate with new or multiple references
#   and move the functionality to find this

data_dir_possibilities = [
    Path("/data"),
    Path(__file__).parent / "data",
]

pretrained_model_paths = {
    "Stanford TMC": {"LI": [Path("models/20260502_64CODEX_stellar_trained_model_origin—version.pt"),
                  Path(__file__).parent / "models/20260502_64CODEX_stellar_trained_model_origin—version.pt"],},
    "General Electric RTI": {"SK": [Path("models/20260612_12CODEX_stellar_trained_model_origin—cell_type_6_alldatasets.pt"),
                  Path(__file__).parent / "models/20260612_12CODEX_stellar_trained_model_origin—cell_type_6_alldatasets.pt"],},
    # other providers: tissues: paths,
}

reference_paths = {
    "Stanford TMC": {"LI": Path("20260107_newSPRM_64CODEX_SLI_annotated.h5ad")},
    "General Electric RTI": {"SK": Path("20260324_skin_v5_12data_leiden15_noCD3_refine_labeled_with_spatial.h5ad")},
    # other providers: tissues: paths,
}


def find_model_file(tissue, provider):
    if provider in pretrained_model_paths.keys():
        for path in pretrained_model_paths[provider][tissue]:
            if path.is_file():
                return path
    else:
        return None


def find_data_file(tissue, provider):
    for path in data_dir_possibilities:
        if (f := path / reference_paths[provider][tissue]).is_file():
            print("Found training data file at", f)
            return f
    message_pieces = [f"Couldn't find data directory; tried:"]
    message_pieces.extend([f"\t{path}" for path in data_dir_possibilities])
    return None


def create_cell_type_manifest(df, outdir):
    cell_type_manifest_dict = {}
    # There is not a CL mapping now, but this is coded to add it easily later.
    for column_header in ['STELLAR_CellType']:
        sub_dict = {
            val: int((df[column_header] == val).sum())
            for val in df[column_header].unique()
        }
        # Remove NaN key if it exists
        sub_dict = {k: v for k, v in sub_dict.items() if not pd.isna(k)}
        cell_type_manifest_dict[column_header] = sub_dict

    with open(f'{outdir}/cell_type_manifest.json', 'w') as f:
        json.dump(cell_type_manifest_dict, f)


def main():
    parser = argparse.ArgumentParser(description="STELLAR")
    parser.add_argument("cell_data_h5ad", type=Path)
    parser.add_argument("tissue", type=str)
    parser.add_argument("provider", type=str)

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
        find_data_file(args.tissue, args.provider),
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

    # Get model path if model exists, exit program if it doesn't
    model_path = find_model_file(args.tissue, args.provider)
    if model_path:
        checkpoint = torch.load(model_path, weights_only=True)
        # Get the number of cell types from the model
        saved_args = checkpoint.get('args', {})
        args.num_heads = saved_args['num_heads']
        stellar = STELLAR(args, dataset)
        stellar.model.load_state_dict(checkpoint['model_state'])
        stellar.model.eval()
        _, results = stellar.pred()

        out_dir = Path("stellar")
        out_dir.mkdir(exist_ok=True, parents=True)
        print(results)
        idxs, annotations = zip(unlabeled_cell_indexes, results)
        predictions_df = pd.DataFrame({'ID': idxs,
                                      'STELLAR_CellType': annotations})
        print(predictions_df)
        annotations_csv = out_dir / f"{args.cell_data_h5ad.stem}.csv"
        predictions_df.to_csv(annotations_csv, index=False)
        create_cell_type_manifest(predictions_df, out_dir)

        # Should I include the accuracy evaluation from Yang's notebook?
        print("done")
    else: # Exit the program nicely if there is no model or reference
        print(f"No pretrained model found for {args.tissue}.")
        # write a csv with only cell IDs and no columns
        out_dir = Path("stellar")
        out_dir.mkdir(exist_ok=True, parents=True)
        with open(out_dir / f"{args.cell_data_h5ad.stem}.csv", "w") as f:
            print("ID,STELLAR_CellType", file=f)
            for cell_id, cell_type_id in unlabeled_cell_indexes:
                print(f"{cell_id}", file=f)
        sys.exit("Exiting STELLAR...")


if __name__ == "__main__":
    main()
