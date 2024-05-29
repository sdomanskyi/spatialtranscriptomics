#!/opt/conda/bin/python

import os
os.environ["NUMBA_CACHE_DIR"] = "./tmp"

import sys
import argparse
import scanpy as sc
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import cm
   
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('True', 'yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('False', 'no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Preprocess single cell transcriptomics data.')

parser.add_argument('--fileNameST', metavar='name', type=str, default='st_data.norm.h5ad', help='')

parser.add_argument('--STdeconvolveSCclusterIds', metavar='name', type=str, default='STdeconvolve_sc_cluster_ids.csv', help='')
parser.add_argument('--STdeconvolvePropNormName', metavar='name', type=str, default='STdeconvolve_prop_norm.csv', help='')
parser.add_argument('--STdeconvolveBetaNormName', metavar='name', type=str, default='STdeconvolve_beta_norm.csv', help='')
parser.add_argument('--STdeconvolveSCloadings', metavar='name', type=str, default='STdeconvolve_sc_pca_feature_loadings.csv', help='')
parser.add_argument('--STdeconvolveSCclusterMarkers', metavar='name', type=str, default='STdeconvolve_sc_cluster_markers.csv', help='')
parser.add_argument('--SPOTlightPropNorm', metavar='name', type=str, default='SPOTlight_prop_norm.csv', help='')
parser.add_argument('--BayesSpaceClusters', metavar='name', type=str, default='bayes_spot_cluster.csv', help='')

parser.add_argument('--nHVGs', metavar='nHVGs', type=int, default=30, help='Number of principal components.')
parser.add_argument('--nPCs', metavar='PC', type=int, default=30, help='Number of principal components.')
parser.add_argument('--resolution', metavar='name', type=float, default=0.4, help='')

parser.add_argument('--scanpy_UMAP_st_sc', metavar='name', type=str, default='scanpy_UMAP_st_sc.png', help='')
parser.add_argument('--scanorama_UMAP_st_sc', metavar='name', type=str, default='scanorama_UMAP_st_sc.png', help='')

parser.add_argument('--Topics_LDA_spatial', metavar='name', type=str, default='Topics_LDA_spatial.png', help='')
parser.add_argument('--Topics_NMF_spatial', metavar='name', type=str, default='Topics_NMF_spatial.png', help='')

parser.add_argument('--UMAP_st_spots_clusters', metavar='name', type=str, default='UMAP_st_spots_clusters.png', help='')
parser.add_argument('--UMAP_clusters_embedding_density', metavar='name', type=str, default='UMAP_clusters_embedding_density.png', help='')
parser.add_argument('--UMAP_LDA_topics', metavar='name', type=str, default='UMAP_LDA_topics.png', help='')
parser.add_argument('--UMAP_NMF_topics', metavar='name', type=str, default='UMAP_NMF_topics.png', help='')

parser.add_argument('--Clusters_scanpy_spatial', metavar='name', type=str, default='Clusters_scanpy_spatial.png', help='')

parser.add_argument('--violin_topics_LDA', metavar='name', type=str, default='violin_topics_LDA.png', help='')
parser.add_argument('--violin_topics_NMF', metavar='name', type=str, default='violin_topics_NMF.png', help='')

parser.add_argument('--saveFileST', metavar='savefile', type=str, default='st_data.processed.h5ad', help='Path to a file to save h5ad data into.')

args = parser.parse_args()


# Main script
# See more settings at:
# https://scanpy.readthedocs.io/en/stable/generated/scanpy._settings.ScanpyConfig.html
sc.settings.figdir = ''
sc.set_figure_params(dpi_save=300, facecolor='white')

for dirName in ['show/', 'umap/', 'umap_density_clusters_/', 'umap_density_sclusters_/', 'violin/', 'heatmap/']:
    if not os.path.exists(dirName):
        os.makedirs(dirName)

# Read preprocessed data
st_adata = sc.read(args.fileNameST)
st_adata.X = st_adata.X.asfptype()

# Add LDA proportions to adata
df_theta = pd.read_csv(args.STdeconvolvePropNormName, index_col=0)
df_theta.columns = 'Topic LDA ' + df_theta.columns.astype(str)
for col in df_theta.columns:
    st_adata.obs[col] = df_theta[col]

# Add spatially aware clusters
BayesSpaceClusters = pd.read_csv(args.BayesSpaceClusters, index_col=0)['spatial.cluster']
st_adata.obs['sclusters'] = (BayesSpaceClusters.reindex(st_adata.obs.index)-1).astype(str).astype('category')

# Clustering and Selection
sc.pp.highly_variable_genes(st_adata, flavor="seurat", n_top_genes=args.nHVGs)
sc.pp.pca(st_adata, n_comps=args.nPCs, zero_center=True, use_highly_variable=True)
sc.pp.neighbors(st_adata)
sc.tl.umap(st_adata)
sc.tl.leiden(st_adata, key_added="clusters", resolution=args.resolution)

identity = "clusters"

if True:
    # Make plots with proportions
    plt.rcParams["figure.figsize"] = (5, 5)
    keys = st_adata.obs.columns[st_adata.obs.columns.str.contains('Topic LDA ')]
    sc.pl.spatial(st_adata, img_key="hires", color=keys, ncols=10, save='/' + args.Topics_LDA_spatial)

    # Make plots of spatial ST spots clusters
    plt.rcParams["figure.figsize"] = (10, 10)
    sc.pl.spatial(st_adata, img_key="hires", color=[identity], save='/' + args.Clusters_scanpy_spatial) # groups=['1', '2', '3']

    # Plot topics proportions in UMAP
    plt.rcParams["figure.figsize"] = (4, 4)
    keys = st_adata.obs.columns[st_adata.obs.columns.str.contains('Topic LDA ')]
    sc.pl.umap(st_adata, color=keys, wspace=0.4, ncols=5, save='/' + args.UMAP_LDA_topics)

    # Make violin plots of topic proportions of ST spots by clusters
    plt.rcParams["figure.figsize"] = (3.5, 3.5)
    keys = st_adata.obs.columns[st_adata.obs.columns.str.contains('Topic LDA ')]
    sc.pl.violin(st_adata, keys, jitter=0.4, groupby=identity, rotation=0, save='/' + args.violin_topics_LDA)

if os.path.isfile(args.SPOTlightPropNorm):
    df_theta = pd.read_csv(args.SPOTlightPropNorm, index_col=0).set_index('barcodes').drop('res_ss', axis=1)
    df_theta.columns = 'Topic NMF ' + df_theta.columns.astype(str)
    for col in df_theta.columns:
        st_adata.obs[col] = df_theta[col]

    # Make plots with proportions
    plt.rcParams["figure.figsize"] = (5, 5)
    keys = st_adata.obs.columns[st_adata.obs.columns.str.contains('Topic NMF ')]
    sc.pl.spatial(st_adata, img_key="hires", color=keys, ncols=10, save='/' + args.Topics_NMF_spatial)

    # Plot topics proportions in UMAP
    plt.rcParams["figure.figsize"] = (4, 4)
    keys = st_adata.obs.columns[st_adata.obs.columns.str.contains('Topic NMF ')]
    sc.pl.umap(st_adata, color=keys, wspace=0.4, ncols=5, save='/' + args.UMAP_NMF_topics)

    # Make violin plots of topic proportions of ST spots by clusters
    plt.rcParams["figure.figsize"] = (3.5, 3.5)
    keys = st_adata.obs.columns[st_adata.obs.columns.str.contains('Topic NMF ')]
    sc.pl.violin(st_adata, keys, jitter=0.4, groupby=identity, rotation=0, save='/' + args.violin_topics_NMF)

# Make plots of UMAP of ST spots clusters
plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(st_adata, color=["clusters", "sclusters", "total_counts", "n_genes_by_counts"], wspace=0.4, save='/' + args.UMAP_st_spots_clusters)

sc.tl.embedding_density(st_adata, basis='umap', groupby=identity)
sc.pl.embedding_density(st_adata, groupby=identity, ncols=10, save='/' + args.UMAP_clusters_embedding_density)

sc.tl.rank_genes_groups(st_adata, identity, method='wilcoxon', key_added="wilcoxon")
sc.pl.rank_genes_groups_heatmap(st_adata, n_genes=10, key="wilcoxon", groupby=identity, show_gene_labels=True, save='/' + 'DEG_heatmap.png')

st_adata.write(args.saveFileST)

exit(0)
