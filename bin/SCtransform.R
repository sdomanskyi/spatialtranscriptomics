#!/usr/local/bin/Rscript

# Load packages
library(argparse)
library(Seurat)
library(scater)


# Parse command-line arguments
parser <- ArgumentParser()

args <- parser$add_argument_group("Agruments", "required and optional arguments")

args$add_argument("--filePath", help = "Path to csv.gz counts file", metavar="dir", required=TRUE)
args$add_argument("--npCountsOutputName", help = "Name of csv.gz counts file", metavar="file", required=TRUE)
args$add_argument("--fileVarRaw", help = "Name of csv.gz var file", metavar="file", required=TRUE)
args$add_argument("--fileObsRaw", help = "Name of csv.gz obs file", metavar="file", required=TRUE)
args$add_argument("--npNormalizedOutputName", help = "Name of csv.gz normalized counts file", metavar="file", required=TRUE)

args <- parser$parse_args()


# Main script
matrix_ct <- read.csv(paste0(args$filePath, args$npCountsOutputName), row.names=1)
genes <- read.csv(paste0(args$filePath, args$fileVarRaw), row.names=1)
obs <- read.csv(paste0(args$filePath, args$fileObsRaw), row.names=1)
rownames(matrix_ct) <- rownames(genes)
colnames(matrix_ct) <- rownames(obs)
print(dim(matrix_ct))

se <- Seurat::CreateSeuratObject(counts=matrix_ct, min.cells=1, min.genes=1, project="data")
print(dim(se))

se <- Seurat::SCTransform(object=se, return.only.var.genes=FALSE, do.correct.umi=TRUE, verbose=FALSE, variable.features.n=2000)

write.csv(as.matrix(se@assays$SCT@data), file=gzfile(paste0(args$filePath, args$npNormalizedOutputName)))
#write.csv(se@assays$SCT@var.features, file=gzfile(paste0(args$filePath, 'var.features.', args$npNormalizedOutputName)))
#write.csv(se@assays$SCT@scale.data, file=gzfile(paste0(args$filePath, 'scale.data.', args$npNormalizedOutputName)))

quit(status=0)
