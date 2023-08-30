if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biobase")
BiocManager::install("limma")
BiocManager::install("SingleCellExperiment")
BiocManager::install("TOAST")

library(Biobase)
library(limma)
library(SingleCellExperiment)
library(dplyr)
library(tidyr)
library(Seurat)
library(ComplexHeatmap)

# install devtools if necessary
install.packages('devtools')

# install the MuSiC package
devtools::install_github('xuranw/MuSiC')

# load
library(MuSiC)


#####
# import integrated and subclustered neuronal neuronal dataset
neurons.subset <- readRDS("./Data/neurons.subset.rds")

# set order of cell types:
level_list <- c()
for(j in 1:19){
  level_list <- c(level_list, paste("e", j, sep = ""))
}
for(j in 1:36){
  level_list <- c(level_list, paste("i", j, sep = ""))
}
neurons.subset@active.ident <- factor(neurons.subset@active.ident,
                                      levels = level_list)

setwd("./Data")

# I already merged the two spreadsheets into one .csv
bulk_df <- as.data.frame(read.csv("bulk_df.csv" ))
rownames(bulk_df) <- bulk_df[,1] # set gene ids as row names
bulk_df <- bulk_df[,-c(1)] # delete gene ids and width parameter

# create data frame with 47 rows (sequencing replicates) and 2 columns (features)
df <- data.frame(matrix(ncol = 2, nrow = 47))

# provide column names
colnames(df) <- c('fish', 'region')
rownames(df) <- colnames(bulk_df)

helpermatrix <- as.data.frame(rownames(df))

df[,1:2] <- helpermatrix %>%
  separate("rownames(df)", c("fish", "region"))

rm(helpermatrix) # clean-up

# make sure features are stored as factors
df$fish <- factor(df$fish)
df$region <- factor(df$region)

DM1_ids <- which(df$region == "DM1")
levels(df$region) <- c(levels(df$region),"DM")
df[c(DM1_ids),]$region <- rep("DM", length(DM1_ids))
df$region <- droplevels(df$region)
# sanity check should return 1
all(rownames(df) == colnames(bulk_df))

exprs_data <- as.matrix(bulk_df)
rownames(exprs_data) <- rownames(bulk_df)

eset <- ExpressionSet(assayData = exprs_data, phenoData = AnnotatedDataFrame(df))
bulk.mtx <- exprs(eset)
colnames(bulk.mtx) <- eset$region
# transform Seurat object into SingleCellExperiment:
sc_eset <- as.SingleCellExperiment(neurons.subset, assay = "RNA")
# set sample to sample identity
sc_eset$sample <- neurons.subset@meta.data$orig.ident


# music2 deconvolution
# music2 deconvolution with TOAST
set.seed(123)
est = music_prop(bulk.mtx = bulk.mtx, sc.sce = sc_eset, clusters = "ident", samples = "sample")
names(est)

# Jitter plot of estimated cell type proportions
jitter.fig = Jitter_Est(list(data.matrix(est$Est.prop.weighted),
                             data.matrix(est$Est.prop.allgene)),
                        method.name = c('MuSiC', 'NNLS'), title = 'Jitter plot of Est Proportions')

df <- as.data.frame(est$Est.prop.weighted)
df$rownames_col <- rownames(est$Est.prop.weighted)
df_avg <- aggregate(. ~ rownames_col, data = df, FUN = mean)
rownames(df_avg) <- df_avg$rownames_col
df_avg$rownames_col <- NULL

library(circlize)
col_fun = colorRamp2(c(0, 0.1, 0.4), c("white", "darkgrey", "black"))
col_fun(seq(-3, 3))

Heatmap(df_avg, col = col_fun, rect_gp = gpar(col = "grey", lwd = 1))
