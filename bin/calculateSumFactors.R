#!/usr/local/bin/Rscript

# Load packages
library(argparse)
library(SpatialExperiment)
library(scran)


# Parse command-line arguments
parser <- ArgumentParser()

args <- parser$add_argument_group("Agruments", "required and optional arguments")

args$add_argument("--filePath", help = "Path to csv.gz counts file", metavar="dir", required=TRUE)
args$add_argument("--npCountsOutputName", help = "Name of csv.gz counts file", metavar="file", required=TRUE)
args$add_argument("--npFactorsOutputName", help = "Name of csv.gz factors file", metavar="file", required=TRUE)

args <- parser$parse_args()


# Main script
matrix_st <- read.csv(paste0(args$filePath, args$npCountsOutputName), row.names=1)
print(dim(matrix_st))

spe <- SpatialExperiment(list(counts=matrix_st))
sfs <- calculateSumFactors(spe, cluster=quickCluster(spe))
print(length(sfs))

write.csv(sfs, gzfile(paste0(args$filePath, args$npFactorsOutputName)), row.names=FALSE)

quit(status=0)
