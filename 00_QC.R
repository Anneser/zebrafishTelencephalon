library(Seurat)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(limma)


setwd("./20220125")
data_dir <- getwd()
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data <- Read10X(data.dir = data_dir)

nUMI = colSums(data)
names(nUMI) = colnames(data)
percent.mito <- colSums(data[grep("^mt-",rownames(data)),])/nUMI

SEO = CreateSeuratObject(counts = data , 
                         min.cells = 5 ,
                         min.features = 200 ,
                         project = "Anneser")

SEO <- AddMetaData(object = SEO, 
                   metadata = percent.mito, 
                   col.name = "percent.mito")

# subset dataset for QC
SEO <- subset(SEO , subset = nFeature_RNA > 500  & 
                percent.mito < 0.1 )


##### second dataset (22.03.2022)
setwd("./20220322")
data_dir <- getwd()
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data <- Read10X(data.dir = data_dir)

nUMI = colSums(data)
names(nUMI) = colnames(data)
percent.mito <- colSums(data[grep("^mt-",rownames(data)),])/nUMI


SEO2 = CreateSeuratObject(counts = data , 
                          min.cells = 5 ,
                          min.features = 200 ,
                          project = "Anneser")

SEO2 <- AddMetaData(object = SEO2, 
                    metadata = percent.mito, 
                    col.name = "percent.mito")

# subset dataset for QC
SEO2 <- subset(SEO2 , subset = nFeature_RNA > 500  & 
                 percent.mito < 0.1 )

#####
setwd("./Cosacak")

data_dir <- getwd()
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data <- Read10X(data.dir = data_dir)

nUMI = colSums(data)
names(nUMI) = colnames(data)
percent.mito <- colSums(data[grep("^mt-",rownames(data)),])/nUMI
percent.gfp = data[rownames(data) %in% "GFP",]/nUMI 


Cosacak = CreateSeuratObject(counts = data , 
                             min.cells = 5 ,
                             min.features = 200 ,
                             project = "Cosacak")

Cosacak <- AddMetaData(object = Cosacak, 
                       metadata = percent.mito, 
                       col.name = "percent.mito")
Cosacak <- AddMetaData(object = Cosacak,
                       metadata = percent.gfp,
                       col.name = "percent.gfp")


# subset dataset for QC
Cosacak <- subset(Cosacak , subset = nFeature_RNA > 500 & 
                    percent.mito < 0.1 )

# for adult fish 1
setwd("./Spanjaard/fish1")
data_dir <- getwd()
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data <- Read10X(data.dir = data_dir)

nUMI = colSums(data)
names(nUMI) = colnames(data)
percent.mito <- colSums(data[grep("^mt-",rownames(data)),])/nUMI


GSE1 = CreateSeuratObject(counts = data , 
                          min.cells = 5 ,
                          min.features = 200 ,
                          project = "Spanjaard")

GSE1 <- AddMetaData(object = GSE1, 
                    metadata = percent.mito, 
                    col.name = "percent.mito")

# for adult fish 2
setwd("./Spanjaard/fish2")
data_dir <- getwd()
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data <- Read10X(data.dir = data_dir)

nUMI = colSums(data)
names(nUMI) = colnames(data)
percent.mito <- colSums(data[grep("^mt-",rownames(data)),])/nUMI


GSE2 = CreateSeuratObject(counts = data , 
                          min.cells = 5 ,
                          min.features = 200 ,
                          project = "Spanjaard")

GSE2 <- AddMetaData(object = GSE2, 
                    metadata = percent.mito, 
                    col.name = "percent.mito")

# merge two objects
spanjaard <- merge(x = GSE1, y = GSE2)
# subset dataset for QC
spanjaard <- subset(spanjaard , subset = nFeature_RNA > 500 & 
                      percent.mito < 0.1 )


# import dataset from Raj et al. 2018
load('./Raj/GSE105010_fall.inDrops_v2.RData')

# update Seurat object (if necessary)
fall <- UpdateSeuratObject(fall)
fall@project.name = "Raj"
# according to Raj et al: 
telencephalon = c(12, 49, 32, 42, 2, 15)

tel_id <- fall@meta.data$Cluster2.5 %in% telencephalon
tel_Raj <- fall[,tel_id]
tel_Raj@meta.data$orig.ident = rep("Raj", length(tel_Raj@meta.data$orig.ident))
tel_Raj <- AddMetaData(object = tel_Raj, 
                       metadata = tel_Raj@meta.data$pt.mito, 
                       col.name = "percent.mito")


Anneser <- merge(x = SEO, y = SEO2)
all_df <- merge(x = Anneser, y = c(Cosacak, spanjaard, tel_Raj))
scData.list <- SplitObject(all_df, split.by = "orig.ident")

# clean up
rm(Cosacak, data, figure, GSE1, GSE2, nCount_plot, nFeat_plot, SEO, SEO2, spanjaard, Anneser)

#normalize and identify variable features for each dataset independently
scData.list <- lapply(X = scData.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

nFeat_plot <- VlnPlot(all_df,
                      features = "nFeature_RNA" ,  group.by = "orig.ident") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("number of detected genes") +
  ggtitle("")

nCount_plot <- VlnPlot(all_df,
                       features = "nCount_RNA" ,  group.by = "orig.ident") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("RNA counts") +
  ggtitle("")

pMito_plot <- VlnPlot(all_df,
                      features = "percent.mito" ,  group.by = "orig.ident") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("fraction of mitochondrial gene counts") +
  ggtitle("")


figure <- ggarrange(nFeat_plot , nCount_plot , pMito_plot,
                    labels = c("A", "B", "C") , 
                    ncol = 3)

figure


# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = scData.list)
scData.anchors <- FindIntegrationAnchors(object.list = scData.list, anchor.features = features)

# this command creates an integrated data assay
scData.combined <- IntegrateData(anchorset = scData.anchors)

# specify that we will perform downstream analysis on the corrected data, note that the original
# unmodified data still resides in the "RNA" assay
DefaultAssay(scData.combined) <- "integrated"

# Run standard workflow
scData.combined <- ScaleData(scData.combined, verbose = FALSE, vars.to.regress = c("nCount_RNA", "nFeature_RNA"))
scData.combined <- RunPCA(scData.combined, npcs = 40, verbose = FALSE)

# continue clustering with 33 PCs
scData.combined <- RunUMAP(scData.combined, reduction = "pca", dims = 1:33)
scData.combined <- FindNeighbors(scData.combined, reduction = "pca", dims = 1:33)

# figure out how to bootstrap resolution for solid clustering
scData.combined <- FindClusters(scData.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(scData.combined, reduction = "umap", group.by = "orig.ident", order = c("Cosacak", "Spanjaard", "Anneser", "Raj"))#, "Raj"))
p2 <- DimPlot(scData.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

# continue here after performing main analysis scripts 01 and 02. 

cells <- readRDS('./Data/telencephalon_cells.rds')
neurons.subset <- readRDS("./Data/neurons.subset.rds")

cells <- AddMetaData(object = cells, 
                   metadata = cells@active.ident, 
                   col.name = "active.ident")

neurons.subset <- AddMetaData(object = neurons.subset, 
                     metadata = neurons.subset@active.ident, 
                     col.name = "active.ident")

nFeat_plot_cells <- VlnPlot(cells,
                      features = "nFeature_RNA" ,  group.by = "active.ident") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("number of detected genes") +
  ggtitle("")

nCount_plot_cells <- VlnPlot(cells,
                       features = "nCount_RNA" ,  group.by = "active.ident") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("RNA counts") +
  ggtitle("")

pMito_plot_cells <- VlnPlot(cells,
                      features = "percent.mito" ,  group.by = "active.ident") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("fraction of mitochondrial gene counts") +
  ggtitle("")

figure <- ggarrange(nFeat_plot_cells , nCount_plot_cells , pMito_plot_cells,
                    labels = c("A", "B", "C") , 
                    ncol = 3)

figure

nFeat_plot_neurons <- VlnPlot(neurons.subset,
                            features = "nFeature_RNA" ,  group.by = "active.ident") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("number of detected genes") +
  ggtitle("")

nCount_plot_neurons <- VlnPlot(neurons.subset,
                             features = "nCount_RNA" ,  group.by = "active.ident") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("RNA counts") +
  ggtitle("")

pMito_plot_neurons <- VlnPlot(neurons.subset,
                            features = "percent.mito" ,  group.by = "active.ident") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("fraction of mitochondrial gene counts") +
  ggtitle("")

figure <- ggarrange(nFeat_plot , nCount_plot , pMito_plot, nFeat_plot_cells , nCount_plot_cells , pMito_plot_cells, nFeat_plot_neurons , nCount_plot_neurons , pMito_plot_neurons,
                    labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I") , 
                    ncol = 3, nrow = 3)

figure
