import os.path
from pathlib import Path

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import SAGEConv


def prepare_save_dir(args, filename):
    """Create saving directory."""
    runner_name = os.path.basename(filename).split(".")[0]
    model_dir = Path(f"experiments/{runner_name}/{args.name}")
    model_dir.mkdir(exist_ok=True, parents=True)
    args.savedir = model_dir
    return args


def entropy(x):
    """
    Helper function to compute the entropy over the batch
    input: batch w/ shape [b, num_classes]
    output: entropy value [is ideally -log(num_classes)]
    """
    EPS = 1e-8
    x_ = torch.clamp(x, min=EPS)
    b = x_ * torch.log(x_)

    if len(b.size()) == 2:  # Sample-wise entropy
        return -b.sum(dim=1).mean()
    elif len(b.size()) == 1:  # Distribution-wise entropy
        return -b.sum()
    else:
        raise ValueError(f"Input tensor is {len(b.size())}-Dimensional")


class MarginLoss(nn.Module):

    def __init__(self, m=0.2, weight=None, s=10):
        super(MarginLoss, self).__init__()
        self.m = m
        self.s = s
        self.weight = weight

    def forward(self, x, target):
        index = torch.zeros_like(x, dtype=torch.bool)
        index.scatter_(1, target.data.view(-1, 1), 1)
        x_m = x - self.m * self.s

        output = torch.where(index, x_m, x)
        return F.cross_entropy(output, target, weight=self.weight)


class MLPHead(nn.Module):
    """
    Your MLPClassifier, adapted to be a reusable head:
    - in_dim can be the GNN hidden size (e.g. 128)
    - returns logits; optionally returns the penultimate 256-d features
    """
    def __init__(self, in_dim, num_classes, dropout=0.1):
        super().__init__()
        self.fc1 = nn.Linear(in_dim, 256)
        self.bn1 = nn.BatchNorm1d(256)
        self.fc2 = nn.Linear(256, 1024)
        self.bn2 = nn.BatchNorm1d(1024)
        self.fc3 = nn.Linear(1024, 256)
        self.bn3 = nn.BatchNorm1d(256)
        self.drop = nn.Dropout(p=dropout)
        self.fc4 = nn.Linear(256, num_classes)

    def forward(self, x, return_feat=False):
        x = F.relu(self.bn1(self.fc1(x))); x = self.drop(x)
        x = F.relu(self.bn2(self.fc2(x))); x = self.drop(x)
        x = F.relu(self.bn3(self.fc3(x))); x = self.drop(x)
        penultimate = x
        logits = self.fc4(penultimate)
        if return_feat:
            return logits, penultimate
        return logits


class Encoder(nn.Module):
    """
    GNN encoder (Linear -> ReLU -> SAGEConv) + your MLP head.
    Returns:
      out       : [N, num_classes] logits
      feat      : penultimate 256-d features from the MLP head (useful for contrastive parts)
      out_feat  : graph features right after SAGEConv (pre-MLP)
    """
    def __init__(self, x_dim, num_cls, hid_dim=256, dropout=0.1):
        super(Encoder, self).__init__()
        self.x_dim = x_dim
        self.relu = nn.ReLU()
        self.conv1 = nn.Linear(x_dim, hid_dim)
        self.conv2 = SAGEConv(hid_dim, hid_dim)
        self.head  = MLPHead(in_dim=hid_dim, num_classes=num_cls, dropout=dropout)

    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        x = self.relu(self.conv1(x))
        x = self.conv2(x, edge_index)
        out_feat = x                          # graph feature before MLP
        out, penultimate = self.head(x, return_feat=True)
        feat = penultimate                    # 256-d features from your MLP
        return out, feat, out_feat


class FCNet(nn.Module):
    """
    Fully-connected baseline that uses your MLP head directly on node features.
    Keeps the (out, feat, out_feat) signature for compatibility.
    """
    def __init__(self, x_dim, num_cls, dropout=0.3):
        super(FCNet, self).__init__()
        self.head = MLPHead(in_dim=x_dim, num_classes=num_cls, dropout=dropout)

    def forward(self, data):
        x = data.x
        out, penultimate = self.head(x, return_feat=True)
        # We don't have a graph feature branch; reuse penultimate for both
        return out, penultimate, penultimate
