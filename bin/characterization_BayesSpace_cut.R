#!/usr/local/bin/Rscript

# Load packages
library(argparse)
library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(Matrix)
#library(reticulate)


# Parse command-line arguments
parser <- ArgumentParser()

args <- parser$add_argument_group("Agruments", "required and optional arguments")

args$add_argument("--filePath", help="Path to csv.gz counts file", metavar="dir", required=TRUE)
args$add_argument("--nameX", default="st_adata_X.csv.gz", help="Path to X", metavar="file", required=FALSE)
args$add_argument("--nameVar", default="st_adata.var.csv", help="Path to features metadata", metavar="file", required=FALSE)
args$add_argument("--nameObs", default="st_adata.obs.csv", help="Path to observation metadata", metavar="file", required=FALSE)
args$add_argument("--countsFactor", default=100, help="factor", metavar="factor", required=FALSE)
args$add_argument("--numberHVG", default=2000, type="integer", help="factor", metavar="factor", required=FALSE)
args$add_argument("--numberPCs", default=30, type="integer", help="factor", metavar="factor", required=FALSE)
args$add_argument("--minClusters", default=2, type="integer", help="factor", metavar="factor", required=FALSE)
args$add_argument("--maxClusters", default=10, type="integer", help="factor", metavar="factor", required=FALSE)
args$add_argument("--optimalQ", default=5, type="integer", help="factor", metavar="factor", required=FALSE)
args$add_argument("--smoothing", default=2, type="double", help="gamma parameter", metavar="factor", required=FALSE)
args$add_argument("--STplatform", default="Visium", help="Technology grid", metavar="factor", required=FALSE)

args$add_argument("--qtuneSaveName", default="st_bayes_qtune.png", help="file name", metavar="file", required=FALSE)
args$add_argument("--bayesClustersName", default="st_bayes_clusters.png", help="file name", metavar="file", required=FALSE)
args$add_argument("--bayesClustersEnhancedName", default="st_bayes_clusters_enhanced.png", help="file name", metavar="file", required=FALSE)
args$add_argument("--bayesOriginalEnhancedFeatures", default="st_bayes_original_and_enhanced.png", help="file name", metavar="file", required=FALSE)
args$add_argument("--bayesEnhancedMarkers", default="bayes_enhanced_markers.csv", help="file name", metavar="file", required=FALSE)
args$add_argument("--bayesSubspotCoord", default="bayes_subspot_cluster_and_coord.csv", help="file name", metavar="file", required=FALSE)
args$add_argument("--bayesSpotCluster", default="bayes_spot_cluster.csv", help="file name", metavar="file", required=FALSE)
args$add_argument("--logNormalize", default=TRUE, help="file name", metavar="file", required=FALSE)
args <- parser$parse_args()

# Main script
set.seed(123)
#np <- import("numpy")
normDataDir <- args$filePath

# Load gene names
rowData <- read.csv(paste0(normDataDir, args$nameVar))$X
row_df <- as.data.frame(rowData)
row.names(row_df) <- rowData

# Load coordinates of spots
st_obs_all <- read.csv(paste0(normDataDir, args$nameObs))
col_df <- as.data.frame(st_obs_all$X)
row.names(col_df) <- st_obs_all$X
col_df$imagerow <- st_obs_all$array_row
col_df$imagecol <- st_obs_all$array_col
colnames(col_df) <- c('id', 'imagerow', 'imagecol')
col_df$row <- st_obs_all$array_row
col_df$col <- st_obs_all$array_col

# Load counts data
count.data <- exp(read.csv(paste0(normDataDir, args$nameX), row.names=1)) - 1
colnames(count.data) <- st_obs_all$X
rownames(count.data) <- rowData

dsp <- SingleCellExperiment(assays=list(counts=count.data), rowData=row_df, colData=col_df)
dsp <- spatialPreprocess(dsp, platform=args$STplatform, n.PCs=args$numberPCs, n.HVGs=args$numberHVG, log.normalize=args$logNormalize)
dsp <- spatialCluster(dsp, q=args$optimalQ, platform=args$STplatform, d=args$numberPCs, init.method="mclust", model="t", gamma=args$smoothing, nrep=1000, burn.in=100, save.chain=FALSE)

clusterPlot(dsp, palette=c("purple", "cyan", "blue", "yellow", "red"), color=NA) + theme_bw() + xlab("Column") + ylab("Row") + labs(fill="BayesSpace\ncluster", title="Spatial clustering")
ggsave(paste0(normDataDir, args$bayesClustersName), dpi=600, scale=0.75, width=8, height=8, units="in")

df_spot_cluster <- subset(as.data.frame(dsp@colData), select=-c(imagecol, imagerow))
write.csv(df_spot_cluster, file=paste0(args$filePath, args$bayesSpotCluster))

quit(status=0)
