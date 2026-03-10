#!/usr/bin/env python3
import importlib.resources
from argparse import ArgumentParser
from pathlib import Path
from math import ceil, log2
import anndata
import pandas as pd
import spatialdata as sd
from ome_utils import find_ome_tiffs
from sklearn.preprocessing import StandardScaler
from sprm import modules
from xarray import DataArray
from spatialdata.models import Image2DModel, Labels2DModel, PointsModel, TableModel
import tracemalloc
from bioio import BioImage

desired_pixel_size_for_pyramid=250

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

    csv_base = expr.name.split('.', 1)[0]

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
        obs=pd.DataFrame(index = [f"{csv_base}-"+str(i) for i in mean_expr.coords["cell_index"].to_series().tolist()]),
        var=pd.DataFrame(index=mean_expr.coords["expr_channel"]),
    )
    print(image_adata.obs.index)
    expr_adata = image_adata.copy()
    # So, create the DataFrame after the AnnData, using .obs_names as the
    # index, to make sure everything matches with minimal effort.
    cell_centers_df = pd.DataFrame(
        core.cell_centers[core.mask.interior_cells],
        index=image_adata.obs_names,
        columns=["X", "Y", "Z"],
    ).loc[:, ["Y", "X"]]
    image_adata.obsm["X_spatial"] = cell_centers_df
    image_adata.obs["unique_region"] = expr.stem

    print('SPRM conversion complete')

    print('Starting SpatialData conversion')
    print("Loading image data")
    image = BioImage(expr)
    image_data_squeezed = image.data.squeeze()
    print("... done. Original shape:", image.data.shape)

    image_scale_factors = (2,) * ceil(
        log2(max(image_data_squeezed.shape[1:]) / desired_pixel_size_for_pyramid)
    )

    img_for_sdata = Image2DModel.parse(
        data=image_data_squeezed,
        c_coords=image.channel_names,
        scale_factors=image_scale_factors,
    )

    print("Loading mask data")
    mask = BioImage(mask)
    mask_data = {ch: mask.data[0, i, 0, :, :] for i, ch in enumerate(mask.channel_names)}
    cell_indexes = sorted(set(mask_data["cells"].flat) - {0})
    print("... done.", len(cell_indexes), "cells")

    masks_for_sdata = {
        ch: Labels2DModel.parse(
            data=mask_array,
            scale_factors=image_scale_factors,
        )
        for ch, mask_array in mask_data.items()
    }

    cell_centers_copy = cell_centers_df.copy()
    cell_centers_copy.columns = ['y', 'x']
    shapes_for_sdata = PointsModel.parse(cell_centers_copy)

    tables = {
        "mean_expr": TableModel.parse(expr_adata),
    }

    for table in tables.values():
        table.obs["cell_id"] = pd.Series([int(i.split('-')[1]) for i in table.obs.index.values], index=table.obs.index)
        table.obs["region"] = pd.Categorical(["cells"] * len(table))
        table.uns["spatialdata_attrs"] = {
            "region": "cells",
            "region_key": "region",
            "instance_key": "cell_id",
        }

    sdata = sd.SpatialData(
        images={"expression": img_for_sdata},
        points={"centers": shapes_for_sdata},
        labels=masks_for_sdata,
        tables=tables,
    )

    sdata_name = f"{csv_base}_spatialdata.zarr"
    print("Saving SpatialData object to", sdata_name)
    print(sdata)
    sdata.write(sdata_name, overwrite=True)

    return image_adata


def main(directory: Path):
    tracemalloc.start()
    expr_dir, mask_dir = find_expr_mask_dir(directory)
    exprs = sorted(find_ome_tiffs(expr_dir))
    masks = sorted(find_ome_tiffs(mask_dir))

    adatas = []
    for expr, mask in zip(exprs, masks):
        adatas.append(convert(expr, mask))

    adata = anndata.concat(adatas, index_unique="-")
    adata.write_h5ad("cell_data.h5ad")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("directory", type=Path)
    args = p.parse_args()

    main(args.directory)
