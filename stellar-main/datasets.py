from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import pairwise_distances
from torch_geometric.data import Data, InMemoryDataset


def get_hubmap_edge_index(pos, regions, distance_thres, *, max_neighbors_per_node=None):
    """
    Build edges within each region using a cKDTree radius query (no dense NxN matrix).
    Optional: cap neighbors per node to further limit edges.
    """
    edge_list = []
    regions_unique = np.unique(regions)

    for reg in regions_unique:
        locs = np.where(regions == reg)[0]
        if len(locs) < 2:
            continue

        P = pos[locs, :]                 # (n_r, 2)
        tree = cKDTree(P)

        # vectorized: list of arrays of neighbor indices (including self)
        neigh = tree.query_ball_point(P, r=float(distance_thres))

        if max_neighbors_per_node is not None:
            # cap neighbors per node within radius using nearest queries
            for i, js in enumerate(neigh):
                # remove self if present
                js = [j for j in js if j != i]
                if max_neighbors_per_node and len(js) > max_neighbors_per_node:
                    # get k nearest within radius using a second query
                    # query returns fixed k; filter by radius + drop self
                    k = max_neighbors_per_node + 1
                    d, idx = tree.query(P[i], k=k, distance_upper_bound=float(distance_thres))
                    # idx/d may include invalid entries when fewer than k neighbors exist
                    valid = [(jj, dd) for jj, dd in zip(np.atleast_1d(idx), np.atleast_1d(d))
                             if np.isfinite(dd) and jj != i]
                    # take up to max_neighbors_per_node nearest
                    js = [jj for jj, _ in valid[:max_neighbors_per_node]]

                for j in js:
                    edge_list.append([locs[i], locs[j]])
        else:
            # no cap: add all neighbors within radius (except self)
            for i, js in enumerate(neigh):
                for j in js:
                    if j != i:
                        edge_list.append([locs[i], locs[j]])

    return edge_list


def get_tonsilbe_edge_index(pos, distance_thres):
    # construct edge indexes in one region
    edge_list = []
    dists = pairwise_distances(pos)
    dists_mask = dists < distance_thres
    np.fill_diagonal(dists_mask, 0)
    edge_list = np.transpose(np.nonzero(dists_mask)).tolist()
    return edge_list


def load_hubmap_data(
    labeled_file: Path, unlabeled_file: Path, distance_thres, sample_rate
):
    train_adata_full = anndata.read_h5ad(labeled_file)
    print("Training data variables:", train_adata_full.var_names)

    test_adata_full = anndata.read_h5ad(unlabeled_file)
    print("Test data variables:", test_adata_full.var_names)

    common_vars = sorted(
        set(train_adata_full.var_names) & set(test_adata_full.var_names)
    )
    print("Common variables:", common_vars)

    test_adata = test_adata_full[:, common_vars].copy()
    unlabeled_pos = test_adata.obsm["X_spatial"]
    unlabeled_regions = test_adata.obs["unique_region"]

    train_sel = np.random.choice(
        range(len(train_adata_full)),
        size=round(sample_rate * len(train_adata_full)),
        replace=False,
    )
    train_adata = train_adata_full[train_sel, common_vars].copy()

    labeled_pos = train_adata.obsm["X_spatial"]
    labeled_regions = train_adata.obs["File_ID"]

    train_y = train_adata.obs["cell_type"]
    cell_types = sorted(set(train_y))
    # we here map class in texts to categorical numbers and also save an inverse_dict to map the numbers back to texts
    cell_type_dict = {}
    inverse_dict = {}
    for i, cell_type in enumerate(cell_types):
        cell_type_dict[cell_type] = i
        inverse_dict[i] = cell_type
    train_y = np.array([cell_type_dict[x] for x in train_y])
    labeled_edges = get_hubmap_edge_index(labeled_pos, labeled_regions, distance_thres)
    unlabeled_edges = get_hubmap_edge_index(
        unlabeled_pos, unlabeled_regions, distance_thres
    )
    return (
        train_adata.X,
        train_y,
        test_adata.X,
        labeled_edges,
        unlabeled_edges,
        inverse_dict,
        test_adata.obs_names,
    )


def load_tonsilbe_data(filename, distance_thres, sample_rate):
    df = pd.read_csv(filename)
    train_df = df.loc[df["sample_name"] == "tonsil"]
    train_df = train_df.sample(n=round(sample_rate * len(train_df)), random_state=1)
    test_df = df.loc[df["sample_name"] == "Barretts Esophagus"]
    train_X = train_df.iloc[:, 1:-4].values
    test_X = test_df.iloc[:, 1:-4].values
    train_y = train_df["cell_type"].str.lower()
    labeled_pos = train_df.iloc[:, -4:-2].values
    unlabeled_pos = test_df.iloc[:, -4:-2].values
    cell_types = np.sort(list(set(train_y))).tolist()
    cell_type_dict = {}
    inverse_dict = {}
    for i, cell_type in enumerate(cell_types):
        cell_type_dict[cell_type] = i
        inverse_dict[i] = cell_type
    train_y = np.array([cell_type_dict[x] for x in train_y])
    labeled_edges = get_tonsilbe_edge_index(labeled_pos, distance_thres)
    unlabeled_edges = get_tonsilbe_edge_index(unlabeled_pos, distance_thres)
    return train_X, train_y, test_X, labeled_edges, unlabeled_edges, inverse_dict


class GraphDataset(InMemoryDataset):
    def __init__(self, labeled_X, labeled_y, unlabeled_X, labeled_edges, unlabeled_edges, transform=None,):
        self.root = '.'
        super(GraphDataset, self).__init__(self.root, transform)
        self.labeled_data = Data(x=torch.FloatTensor(labeled_X), edge_index=torch.LongTensor(labeled_edges).T, y=torch.LongTensor(labeled_y))
        self.unlabeled_data = Data(x=torch.FloatTensor(unlabeled_X), edge_index=torch.LongTensor(unlabeled_edges).T)
    def __len__(self):
        return 2
    def __getitem__(self, idx):
        return self.labeled_data, self.unlabeled_data