#!/usr/local/bin/Rscript

library(argparse)
library(Seurat)
library(SeuratDisk)

# Parse command-line arguments
parser <- ArgumentParser()
args <- parser$add_argument_group("Agruments", "required and optional arguments")
args$add_argument("--processedDataDir", help="Path to outputs", metavar="dir", required=TRUE)
args$add_argument("--rawDataDir", help="Path to image source data", metavar="dir", required=TRUE)
args <- parser$parse_args()

# Convert h5ad to h5seurat
SeuratDisk::Convert(paste0(args$processedDataDir, 'st_adata_raw.h5ad'), dest="h5seurat", overwrite=TRUE)
print("Created h5 data raw")

SeuratDisk::Convert(paste0(args$processedDataDir, 'st_adata_processed.h5ad'), dest="h5seurat", overwrite=TRUE)
print("Created h5 data processed")

# Create Seurat objects
se_st_raw <- SeuratDisk::LoadH5Seurat(paste0(args$processedDataDir, 'st_adata_raw.h5seurat'), misc=FALSE)
se_st <- SeuratDisk::LoadH5Seurat(paste0(args$processedDataDir, 'st_adata_processed.h5seurat'), misc=FALSE)

# Add raw counts
names <- dimnames(se_st@assays[["RNA"]]@counts)
se_st@assays[["RNA"]]@counts <- se_st_raw@assays[["RNA"]]@counts[names[[1]], names[[2]]]

# Remove temporary files and objects
rm(se_st_raw, names)
unlink(paste0(args$processedDataDir, 'st_adata_raw.h5seurat'))
unlink(paste0(args$processedDataDir, 'st_adata_processed.h5seurat'))

# Assign spatially-aware cluters as main ident in Seurat object
Idents(object=se_st) <- se_st@meta.data[["sclusters"]]

# Load and assign image to Seurat object
image <- Read10X_Image(image.dir=file.path(args$rawDataDir, 'spatial'), filter.matrix=TRUE)[Cells(x=se_st)]
DefaultAssay(object=image) <- "RNA"
se_st[["slice1"]] <- image
rm(image)

# Save RDS object containing Seurat object
saveRDS(object=se_st, file=paste0(args$processedDataDir, 'seurat_st.rds'))

quit(status=0)
