#!/usr/local/bin/Rscript

# Load packages
library(argparse)
library(dplyr)
library(ggplot2)
library(Seurat)
library(STdeconvolve)
library(corrplot)

# Parse command-line arguments
parser <- ArgumentParser()

args <- parser$add_argument_group("Agruments", "required and optional arguments")

args$add_argument("--filePath", help="Path to csv.gz counts file", metavar="dir", required=TRUE)
args$add_argument("--outsPath", help="Path to data", metavar="dir", required=TRUE)

args$add_argument("--nameX", default="st_adata_X.csv.gz", help="Path to X", metavar="file", required=FALSE)
args$add_argument("--nameVar", default="st_adata.var.csv", help="Path to features metadata", metavar="file", required=FALSE)
args$add_argument("--nameObs", default="st_adata.obs.csv", help="Path to observation metadata", metavar="file", required=FALSE)

args$add_argument("--useSCdata", default=TRUE, help="Mode", metavar="flag", required=FALSE)

args$add_argument("--SCnameX", default="sc_adata_X.csv.gz", help="Path to X", metavar="file", required=FALSE)
args$add_argument("--SCnameVar", default="sc_adata.var.csv", help="Path to features metadata", metavar="file", required=FALSE)
args$add_argument("--SCnameObs", default="sc_adata.obs.csv", help="Path to observation metadata", metavar="file", required=FALSE)

args$add_argument("--fileh5", default="raw_feature_bc_matrix.h5", help="File HDF5", metavar="file", required=FALSE) # "filtered_feature_bc_matrix.h5"

args$add_argument("--outsSubDir", default="raw_feature_bc_matrix/", help="dir", metavar="file", required=FALSE)
args$add_argument("--mtxGeneColumn", default=2, type="integer", help="columns index", metavar="col", required=FALSE)
args$add_argument("--countsFactor", default=100, type="integer", help="factor", metavar="factor", required=FALSE)

args$add_argument("--corpusRemoveAbove", default=1.0, type="double", help="factor", metavar="factor", required=FALSE)
args$add_argument("--corpusRemoveBelow", default=0.05, type="double", help="factor", metavar="factor", required=FALSE)

args$add_argument("--LDAminTopics", default=8, type="integer", help="factor", metavar="factor", required=FALSE)
args$add_argument("--LDAmaxTopics", default=9, type="integer", help="factor", metavar="factor", required=FALSE)
args$add_argument("--LDAsaveFile", default="STdeconvolve_optLDA.rds", help="File to save LDA RDS", metavar="file", required=FALSE)

args$add_argument("--STdeconvolveScatterpiesName", default="STdeconvolve_st_scatterpies.png", help="dir", metavar="file", required=FALSE)
args$add_argument("--STdeconvolveScatterpiesSize", default=2.85, type="double", help="dir", metavar="file", required=FALSE)
args$add_argument("--STdeconvolveFeaturesSizeFactor", default=1.0, type="double", help="dir", metavar="file", required=FALSE)
args$add_argument("--STdeconvolveFeaturesName", default="STdeconvolve_st_prop.png", help="dir", metavar="file", required=FALSE)
args$add_argument("--STdeconvolveCorrName", default="STdeconvolve_st_prop_corr.png", help="dir", metavar="file", required=FALSE)
args$add_argument("--STdeconvolvePropNormName", default="STdeconvolve_prop_norm.csv", help="dir", metavar="file", required=FALSE)
args$add_argument("--STdeconvolveBetaNormName", default="STdeconvolve_beta_norm.csv", help="dir", metavar="file", required=FALSE)
args$add_argument("--STdeconvolveSCclustersName", default="STdeconvolve_sc_clusters.png", help="dir", metavar="file", required=FALSE)
args$add_argument("--STdeconvolveSCclusterIds", default="STdeconvolve_sc_cluster_ids.csv", help="dir", metavar="file", required=FALSE)
args$add_argument("--STdeconvolveSCpca", default="STdeconvolve_sc_pca.csv", help="dir", metavar="file", required=FALSE)
args$add_argument("--STdeconvolveSCloadings", default="STdeconvolve_sc_pca_feature_loadings.csv", help="dir", metavar="file", required=FALSE)
args$add_argument("--STdeconvolveSCclusterMarkers", default="STdeconvolve_sc_cluster_markers.csv", help="dir", metavar="file", required=FALSE)
args$add_argument("--trainedLDA", default="", help="dir", metavar="file", required=FALSE)

args <- parser$parse_args()

# Main script
set.seed(123)
#np <- import("numpy")
normDataDir <- args$filePath

filename <- list.files(path=args$outsPath, pattern=args$fileh5)[1]
print(args$outsPath)


library("jsonlite")
library("png")

Read10X_Image <- function(image.dir, image.name = "tissue_lowres_image.png", filter.matrix = TRUE, ...) {
  image <- readPNG(source = file.path(image.dir, image.name))
  scale.factors <- fromJSON(txt = file.path(image.dir, 'scalefactors_json.json'))
  tissue.positions.path <- Sys.glob(paths = file.path(image.dir, 'tissue_positions*'))
  tissue.positions <- read.csv(
    file = tissue.positions.path,
    col.names = c('barcodes', 'tissue', 'row', 'col', 'imagerow', 'imagecol'),
    header = ifelse(
      test = basename(tissue.positions.path) == "tissue_positions.csv",
      yes = TRUE,
      no = FALSE
    ),
    as.is = TRUE,
    row.names = 1
  )
  if (filter.matrix) {
    tissue.positions <- tissue.positions[which(x = tissue.positions$tissue == 1), , drop = FALSE]
  }
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
  spot.radius <-  unnormalized.radius / max(dim(x = image))
  return(new(
    Class = 'VisiumV1',
    image = image,
    scale.factors = scalefactors(
      spot = scale.factors$tissue_hires_scalef,
      fiducial = scale.factors$fiducial_diameter_fullres,
      hires = scale.factors$tissue_hires_scalef,
      scale.factors$tissue_lowres_scalef
    ),
    coordinates = tissue.positions,
    spot.radius = spot.radius
  ))
}


if (!is.na(filename)) {
    print(filename)
    se_st <- Seurat::Load10X_Spatial(data.dir = args$outsPath, filename = filename)
} else {
    image <- Read10X_Image(image.dir=file.path(args$outsPath, 'spatial'), filter.matrix=TRUE)
    m <- Read10X(paste0(args$outsPath, args$outsSubDir), gene.column=args$mtxGeneColumn)
    m <- m[,row.names(image@coordinates)]
    m <- m[,colSums(m)>0]
    se_st <- CreateSeuratObject(counts=m, assay="Spatial")
    image <- image[Cells(x=se_st)]
    DefaultAssay(object=image) <- "Spatial"
    se_st[["slice1"]] <- image
}

matrix_st <- data.matrix(read.csv(paste0(normDataDir, args$nameX), row.names=1))
st_genes <- read.csv(paste0(normDataDir, args$nameVar))$X
st_obs <- read.csv(paste0(normDataDir, args$nameObs))$X
rownames(matrix_st) <- st_genes
colnames(matrix_st) <- st_obs
se_st@assays$Spatial@counts <- as(matrix_st, "sparseMatrix")
se_st@assays$Spatial@data <- as(matrix_st, "sparseMatrix")

se_st@assays$Spatial@counts <- as(round(exp(se_st@assays$Spatial@counts) - 1), Class='dgCMatrix')

corpus <- restrictCorpus(se_st@assays$Spatial@counts, removeAbove=args$corpusRemoveAbove, removeBelow=args$corpusRemoveBelow)
corpus <- corpus + 1
print(dim(as.matrix(corpus)))
print(sum(colSums(as.matrix(corpus))==0))


if (args$trainedLDA=="") {
    ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(args$LDAminTopics, args$LDAmaxTopics, by = 1))
    optLDA <- optimalModel(models = ldas, opt = "min")
    saveRDS(object = optLDA, file = paste0(args$filePath, args$LDAsaveFile))
} else {
    print("Loading provided LDA")
    optLDA <- readRDS(file = args$trainedLDA)
}

results <- getBetaTheta(optLDA, t(as.matrix(corpus)))
deconProp <- results$theta
posk <- se_st@images[["slice1"]]@coordinates[c("imagerow", "imagecol")]
names(posk) <- c('y', 'x')
posk$x <- posk$x * se_st@images[["slice1"]]@scale.factors[["lowres"]]
posk$y <- dim(se_st@images[["slice1"]])[1] - posk$y * se_st@images[["slice1"]]@scale.factors[["lowres"]]
vizAllTopics(deconProp, posk, r=args$STdeconvolveScatterpiesSize, lwd=0, overlay=se_st@images[["slice1"]]@image)
ggsave(paste0(args$filePath, args$STdeconvolveScatterpiesName), dpi=600, scale=1.0, width=8, height=8, units='in')


# Individual topics proportions spatial overlay
decon_df <- deconProp %>% data.frame() %>% tibble::rownames_to_column("barcodes")

print(colnames(decon_df))
print(dim(decon_df))

# Cell type proportions correlation
decon_mtrx_sub <- deconProp[, colnames(deconProp)[which(colnames(deconProp) != "res_ss")]]
decon_mtrx_sub <- decon_mtrx_sub[, colSums(decon_mtrx_sub) > 0]
colnames(decon_mtrx_sub) <- paste("X", colnames(decon_mtrx_sub), sep = "")
decon_cor <- cor(decon_mtrx_sub)
p.mat <- corrplot::cor.mtest(mat = decon_mtrx_sub, conf.level = 0.95)
ggcorrplot::ggcorrplot(corr = decon_cor, p.mat = p.mat[[1]], hc.order = TRUE, type = "full", insig = "blank",
  lab = TRUE, outline.col = "lightgrey", method = "square", colors = c("#6D9EC1", "white", "#E46726"),
  title = "Cell type proportions correlation", legend.title = "Correlation\n(Pearson)") +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 22, hjust = 0.5, face = "bold"),
    legend.text = ggplot2::element_text(size = 12), legend.title = ggplot2::element_text(size = 15),
    axis.text.x = ggplot2::element_text(angle = 90), axis.text = ggplot2::element_text(size = 18, vjust = 0.5))
ggsave(paste0(args$filePath, args$STdeconvolveCorrName), dpi=600, scale=0.75, width=8, height=8, units='in')

write.csv(results$theta, file=paste0(args$filePath, args$STdeconvolvePropNormName))
write.csv(t(results$beta), file=paste0(args$filePath, args$STdeconvolveBetaNormName))

quit(status=0)
