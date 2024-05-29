import os
import sys
import argparse
import scanpy as sc
import numpy as np

from scanpy import read_10x_mtx
from pathlib import Path
from typing import Union, Dict, Optional
import json
import numpy as np
import pandas as pd
from matplotlib.image import imread
import anndata
from anndata import AnnData, read_csv

from matplotlib import pyplot as plt
import scipy.stats
from scipy.sparse import csr_matrix

def resample_counts_inplace(adata, downsample=True, upsample=True, target_total=2*10**3, n_iter=10, seed=None):
    
    is_csr_matrix = type(adata.X) is csr_matrix

    if not seed is None:
        np.random.seed(seed)

    for i in range(n_iter):
        v = adata.to_df().values.astype(int)
        totals_per_obs = v.sum(axis=1)

        if downsample:
            wh = np.where(totals_per_obs > target_total)[0]
            for i in wh:
                wh_nonzero_genes = np.where(v[i] > 0)[0]
                diff = totals_per_obs[i] - target_total
                for j in np.random.choice(wh_nonzero_genes, diff, replace=True):
                    if v[i, j] > 0:
                        v[i, j] -= 1

        if upsample:
            wh = np.where(totals_per_obs < target_total)[0]
            for i in wh:
                wh_nonzero_genes = np.where(v[i] > 0)[0]
                diff = target_total - totals_per_obs[i]
                for j in np.random.choice(wh_nonzero_genes, diff, replace=True):
                    if v[i, j] > 0:
                        v[i, j] += 1

        adata.X = v
    
    if is_csr_matrix:
        adata.X = csr_matrix(adata.X)

    return


def histplotQC(se_data, bins, ax):
    try:
        ax.hist(se_data, density=True, bins=bins, color='navy', alpha=0.3)
        kde = scipy.stats.gaussian_kde(se_data)
        xx = np.linspace(min(se_data), max(se_data), 300)
        ax.set_xlabel(se_data.name)
        ax.set_ylabel('Density')
        ax.plot(xx, kde(xx), color='crimson')
        ax.set_xlim([0, ax.get_xlim()[1]])
    except:
        pass
    return


def read_visium_mtx(
    path: Union[str, Path],
    genome: Optional[str] = None,
    *,
    count_file: str = "filtered_feature_bc_matrix.h5",
    library_id: str = None,
    load_images: Optional[bool] = True,
    source_image_path: Optional[Union[str, Path]] = None,
) -> AnnData:
    """\
    Read 10x-Genomics-formatted visum dataset.
    In addition to reading regular 10x output,
    this looks for the `spatial` folder and loads images,
    coordinates and scale factors.
    Based on the `Space Ranger output docs`_.
    See :func:`~scanpy.pl.spatial` for a compatible plotting function.
    .. _Space Ranger output docs: https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/overview
    Parameters
    ----------
    path
        Path to directory for visium datafiles.
    genome
        Filter expression to genes within this genome.
    count_file
        Which file in the passed directory to use as the count file. Typically would be one of:
        'filtered_feature_bc_matrix.h5' or 'raw_feature_bc_matrix.h5'.
    library_id
        Identifier for the visium library. Can be modified when concatenating multiple adata objects.
    source_image_path
        Path to the high-resolution tissue image. Path will be included in
        `.uns["spatial"][library_id]["metadata"]["source_image_path"]`.
    Returns
    -------
    Annotated data matrix, where observations/cells are named by their
    barcode and variables/genes by gene name. Stores the following information:
    :attr:`~anndata.AnnData.X`
        The data matrix is stored
    :attr:`~anndata.AnnData.obs_names`
        Cell names
    :attr:`~anndata.AnnData.var_names`
        Gene names
    :attr:`~anndata.AnnData.var`\\ `['gene_ids']`
        Gene IDs
    :attr:`~anndata.AnnData.var`\\ `['feature_types']`
        Feature types
    :attr:`~anndata.AnnData.uns`\\ `['spatial']`
        Dict of spaceranger output files with 'library_id' as key
    :attr:`~anndata.AnnData.uns`\\ `['spatial'][library_id]['images']`
        Dict of images (`'hires'` and `'lowres'`)
    :attr:`~anndata.AnnData.uns`\\ `['spatial'][library_id]['scalefactors']`
        Scale factors for the spots
    :attr:`~anndata.AnnData.uns`\\ `['spatial'][library_id]['metadata']`
        Files metadata: 'chemistry_description', 'software_version', 'source_image_path'
    :attr:`~anndata.AnnData.obsm`\\ `['spatial']`
        Spatial spot coordinates, usable as `basis` by :func:`~scanpy.pl.embedding`.
    """
    path = Path(path)
    adata = read_10x_mtx(path / 'raw_feature_bc_matrix')
    
    adata.uns["spatial"] = dict()

    if library_id is None:
        library_id = 'library_id'

    adata.uns["spatial"][library_id] = dict()

    if load_images:
        files = dict(
            tissue_positions_file=path / 'spatial/tissue_positions_list.csv',
            tissue_positions_file_2=path / 'spatial/tissue_positions.csv',
            scalefactors_json_file=path / 'spatial/scalefactors_json.json',
            hires_image=path / 'spatial/tissue_hires_image.png',
            lowres_image=path / 'spatial/tissue_lowres_image.png',
        )

        adata.uns["spatial"][library_id]['images'] = dict()
        for res in ['hires', 'lowres']:
            try:
                adata.uns["spatial"][library_id]['images'][res] = imread(
                    str(files[f'{res}_image'])
                )
            except Exception:
                raise OSError(f"Could not find '{res}_image'")

        # read json scalefactors
        adata.uns["spatial"][library_id]['scalefactors'] = json.loads(
            files['scalefactors_json_file'].read_bytes()
        )
        
        adata.uns["spatial"][library_id]["metadata"] = {k: "NA" for k in ("chemistry_description", "software_version")}

        # read coordinates
        if os.path.isfile(files['tissue_positions_file']):
            positions = pd.read_csv(files['tissue_positions_file'], header=None)
            positions.columns = [
                'barcode',
                'in_tissue',
                'array_row',
                'array_col',
                'pxl_col_in_fullres',
                'pxl_row_in_fullres',
            ]
            positions.index = positions['barcode']
        elif os.path.isfile(files['tissue_positions_file_2']):
            positions = pd.read_csv(files['tissue_positions_file_2'], header=0, index_col=0)
            positions.columns = ['in_tissue', 'array_row', 'array_col', 'pxl_col_in_fullres', 'pxl_row_in_fullres']
            positions['barcode'] = positions.index.values
        else:
            raise NotImplementedError

        adata.obs = adata.obs.join(positions, how="left")

        adata.obsm['spatial'] = adata.obs[
            ['pxl_row_in_fullres', 'pxl_col_in_fullres']
        ].to_numpy()
        adata.obs.drop(
            columns=['barcode', 'pxl_row_in_fullres', 'pxl_col_in_fullres'],
            inplace=True,
        )

        # put image path in uns
        if source_image_path is not None:
            # get an absolute path
            source_image_path = str(Path(source_image_path).resolve())
            adata.uns["spatial"][library_id]["metadata"]["source_image_path"] = str(
                source_image_path
            )

    return adata
