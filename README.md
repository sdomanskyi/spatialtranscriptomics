
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)

## Introduction

**ST-downstream-processing** is a bioinformatics best-practice analysis pipeline for 10x Visium spatial transcriptomics data downstream analysis. The pipeline for processing spatially-resolved gene counts with spatial coordinates, image data, and optionally single cell RNA-seq data, designed for 10x genomics Visium Spatial Geen Expression and Single Cell transcriptomics data. Specifically, input data can be output of 10x SpaceRanger and CellRanger.

The are numerous methods for ST data analysis, and this research area is rapidly developing. The pipeline may be useful to a community working in the area of ST. The pipeline that performs a set of analyses, not limited to but including quality control, normalization, deconvolution of spots into cell types and topics, resolution enhancement, clustering, selection, spatially-variable genes, etc. The software is made of R and Python with the software packages containerized.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies.

## Pipeline summary

1. Normalization and Quality Control ([`scran`](https://doi.org/doi:10.18129/B9.bioc.scran), [`scanpy`](https://github.com/theislab/scanpy))
2. Spots cell types and cell topics deconvolution ([`STdeconvolve`](https://jef.works/STdeconvolve/), [`SPOTlight`](https://github.com/MarcElosua/SPOTlight))
3. Spot resolution enhancement ([`BayesSpace`](https://github.com/edward130603/BayesSpace))
4. Dimensionality reduction and projection ([`scanpy`](https://github.com/theislab/scanpy), [`Seurat`](https://satijalab.org/seurat/))
5. Integration with scRNA-seq data ([`scanorama`](https://github.com/brianhie/scanorama), [`BBKNN`](https://github.com/Teichlab/bbknn)) 
6. Clustering of the spots ([`scanpy Leiden`](https://arxiv.org/abs/1810.08473), [`BayesSpace`](https://github.com/edward130603/BayesSpace))
7. Visualization of clusters and features in spatial coordinates and 2D projection layout ([`scanpy`](https://github.com/theislab/scanpy), [`Seurat`](https://satijalab.org/seurat/))
8. Identification of spatially variable features ([`SpatialDE`](https://github.com/Teichlab/SpatialDE))


The pipeline combines multiple tools, toolkits and platforms:

+ [`STdeconvolve`](https://jef.works/STdeconvolve/) - R implementation of LDA-based cell-topics deconvolution of spots.
+ [`SPOTlight`](https://github.com/MarcElosua/SPOTlight) - R implementation of NMF-based cell-types deconvolution of spots.
+ [`BayesSpace`](https://github.com/edward130603/BayesSpace) - R package for spatially-aware clustering and resolution enhancement.
+ [`SpatialDE`](https://github.com/Teichlab/SpatialDE) - Python package for identification of spatially variable genes.
+ [`scanpy`](https://github.com/theislab/scanpy) - scalable Python toolkit for analyzing single-cell gene expression data.
+ [`anndata`](https://github.com/theislab/anndata) - a Python package for handling annotated data matrices in memory and on disk.
+ [`Bioconductor`](https://www.bioconductor.org/) - software resource for the analysis of genomic data. Based primarily on the statistical R programming language.
+ [`Seurat`](https://satijalab.org/seurat/) - R toolkit for single cell genomics.
+ [`scran`](https://doi.org/doi:10.18129/B9.bioc.scran) - R package implements miscellaneous functions for analysis and interpretation of single-cell RNA-seq data.
+ [`SpatialExperiment`](https://doi.org/doi:10.18129/B9.bioc.SpatialExperiment) - R package for storing, retrieving spatial coordinates and gene expression.

## Quick Start

> **Note:** As a temporary measure the singularity containers necessary to run this pipeline were uploaded to https://doi.org/10.5281/zenodo.6266243. Download the two containers and edit the `nextflow.config` parameters `container_python` and `container_r`.

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/)

> **Note:** Test datasets and their description are located at [`nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/spatialtranscriptomics).

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```console
    git clone https://github.com/TheJacksonLaboratory/ST-downstream-processing.git --branch dev --single-branch
    ```

    Allocate resources to run an interactive job, run a quick test (~5-10 min):

    ```console
    srun -p dev -q dev -t 8:00:00 --cpus-per-task=1 --mem=2G -J ijob --pty /bin/bash
    ```

    Run a SLURM-based quick test of the pipeline (~10 min):
    ```console
    nextflow run main.nf -profile test,slurm,singularity
    ```

4. To run your own analysis edit file `run.sh` and run it:

    ```console
    ./run.sh
    ```

    Example samplesheet.csv:

    ```
    sample_id,species,st_data_dir,sc_data_dir,sc_annotation_data_dir,sc_annotation_counts,sc_annotation_labels
    sample1,Mouse,/path/to/spaceranger/outs1/,/path/to/filtered_feature_bc_matrix1/,,,
    sample2,Human,/path/to/spaceranger/outs2/,,,,
    sample3,Human,/path/to/spaceranger/outs3/,,/path/to/annotated/data/,counts.csv.gz,celltypes.csv.gz
    ```

## Documentation

> Note: The documentation is under development and will be updated as soon as possible.

The detailed documentation is split into the following pages:

* [Usage](docs/usage.md)
    * An overview of how the pipeline works, how to run it and a description of all of the different command-line flags.
* [Output](docs/output.md)
    * An overview of the different results produced by the pipeline and how to interpret them.


## Credits

The pipeline was originally developed by The Jackson Laboratory. This project has been supported by grants from the US National Institutes of Health [U24CA224067](https://reporter.nih.gov/project-details/10261367) and [U54AG075941](https://reporter.nih.gov/project-details/10376627). Original authors:

+ [Dr. Sergii Domanskyi](https://github.com/sdomanskyi)
+ Prof. Jeffrey Chuang
+ Dr. Anuj Srivastava

The pipeline is further refactored and simplified in collaboration with the [National Genomics Infastructure](https://ngisweden.scilifelab.se/) within [SciLifeLab](https://scilifelab.se/): https://github.com/nf-core/spatialvi


## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
