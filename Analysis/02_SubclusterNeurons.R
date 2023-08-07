# library import
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(clustree)
library(dendextend)
library(ComplexHeatmap)

# define functions
##### helper functions
# by Benjamin Doran, slightly adapted by LA

ChooseClusterResolutionDownsample <- function(
    input.srobj, n.pcs, sample.name =  format(paste(Sys.Date(), "SilhouetteScore", sep="_")),
    res.low = .01, res.high=3, res.n = 40, bias = "over", print_fig = F, figdir=F, coarseness = 3) {
  
  ######## step 1: save the input seurat object as a new temporary object, 
  ########         dont want to overwrite or change the original one with all of the parameter scans
  
  srobj.tmp = input.srobj 
  # in case there have been other things calculated in the metadata, just cut down to simplify/avoid errors
  srobj.tmp@meta.data = srobj.tmp@meta.data[,c(2:3)] # should just be the nUMI and nGene
  
  ######## step 2: calculate the FindClusters over a large range of resolutions
  print("Performing parameter scan over multiple resolutions...")
  
  set.res = round(exp(seq(log(res.low), log(res.high), length.out=res.n)), digits=3)
  srobj.tmp = FindClusters(srobj.tmp, dims.use = n.pcs, k.param=20,
                           resolution=set.res[1], save.SNN=T, 
                           verbose=FALSE)
  
  for(i in 2:length(set.res)){
    srobj.tmp = FindClusters(
      srobj.tmp, resolution=set.res[i], verbose=FALSE)
    print(paste("          ", round(100*i/length(set.res)), "% done with parameter scan", sep=""))
  }
  
  ######## step 3: output plot of how the resolution changes the number of clusters you get
  n.clusters = vector(mode="numeric", length=length(set.res))
  names(n.clusters) = set.res
  for(i in 1:length(n.clusters)){
    n.clusters[i] = length(table(as.vector(srobj.tmp@meta.data[,paste("integrated_snn_res.", names(n.clusters)[i], sep="")])))
  }
  
  ######## step 4: calculate the silhouette width for each resolution
  print("Computing a silhouette width for each cell, for each resolution...")
  require(cluster)
  
  dist.temp = cor(t(srobj.tmp@reductions$pca@cell.embeddings[,1:n.pcs]), method="pearson")
  random.cells.choose = sample(1:nrow(dist.temp), round(nrow(dist.temp)/10, digits=0))
  dist.temp.downsample = dist.temp[random.cells.choose, random.cells.choose]
  sil.all.matrix = matrix(data=NA, nrow=nrow(dist.temp.downsample), ncol=0)
  
  for(i in 1:length(set.res)){
    clusters.temp = as.numeric(as.vector(
      srobj.tmp@meta.data[random.cells.choose,paste("integrated_snn_res.", set.res[i], sep="")]))
    if(length(table(clusters.temp))>1){
      sil.out = silhouette(clusters.temp, as.dist(1-as.matrix(dist.temp.downsample)))
      sil.all.matrix = cbind(sil.all.matrix, sil.out[,3])
    }
    if(length(table(clusters.temp))==1){
      sil.all.matrix = cbind(sil.all.matrix, rep(0, length(clusters.temp)))
    }
    print(paste("          ", round(100*i/length(set.res)), "% done with silhouette calculation", sep=""))
    
  }
  
  ######## step 5: calculate summary metric to compare the silhouette distributions,
  ########         average has worked well so far... could get fancier
  
  print("Identifying a best resolution to maximize silhouette width")
  sil.average = colMeans(sil.all.matrix)
  names(sil.average) = set.res
  
  ######## step 6: automate choosing resolution that maximizes the silhouette 
  hist.out = hist(sil.average, length(sil.average)/1.2,  plot=FALSE)
  
  #  take the ones that fall into the top bin, 
  #  and the max OR MIN of those  ******* can change this to under vs over cluster
  if(bias=="over"){
    resolution.choice = as.numeric(max(
      names(sil.average[which(sil.average>hist.out$breaks[length(hist.out$breaks)-1])])))
  }
  if(bias=="under"){
    resolution.choice = as.numeric(min(
      names(sil.average[which(sil.average>hist.out$breaks[length(hist.out$breaks)-1])])))
  }
  
  # get the silhouette of the best resolution: 
  silhouette.best = as.numeric(sil.average[paste(resolution.choice)])
  
  print(paste("Best Resolution Choice: ", resolution.choice, ", with average silhouette score of: ",
              round(silhouette.best, digits=3), ", giving ", as.numeric(n.clusters[paste(resolution.choice)]),
              " clusters", sep=""))
  
  ######### step 7: output plot and data 
  
  if (print_fig) {
    setwd(figdir)
    
    print(paste0("Outputting summary statistics and returning seurat object... ",
                 "This will create a pdf in your output directory,",
                 " and will return your input seurat object amended with the best choice multiplied by 3",
                 " for clusters (found as Best.Clusters in the meta.data matrix, and set to your new ident)..."))
    
    pdf(paste(sample.name, ".pdf", sep=""),
        width=10, height=4, useDingbats=FALSE)
    par(mfrow=c(1,3))
    # Resolution vs # of Clusters
    plot(set.res, n.clusters, col="black", pch=19,
         type="p", xlab="Resolution", ylab="# Clusters",
         main="Resolution vs. # Clusters")
    # Resolution vs Average Silhouette
    plot(set.res, sil.average, col="black", pch=19,
         type="p", xlab="Resolution", ylab="Average Silhouette",
         main="Resolution vs. Average Silhouette")
    abline(h=hist.out$breaks[length(hist.out$breaks)-1], col="firebrick3", lty=2)
    abline(v=resolution.choice, col="dodgerblue2", lty=2)
    
    # N Clusters vs Average Silhouette
    plot(n.clusters, sil.average, col="black", pch=19,
         type="p", xlab="# Clusters", ylab="Average Silhouette",
         main="# Clusters vs. Average Silhouette")
    abline(h=hist.out$breaks[length(hist.out$breaks)-1], col="firebrick3", lty=2)
    abline(v=as.numeric(n.clusters[paste(resolution.choice)]), col="dodgerblue2", lty=2)
    dev.off()
  }
  
  ######## step 8: return the original seurat object, with the metadata containing a 
  ########         concatenated vector with the clusters defined by the best choice here,
  ########         as well as the ident set to this new vector
  
  resolution.choice <- resolution.choice * coarseness
  
  input.srobj <- FindClusters(input.srobj, dims.use = n.pcs,
                              resolution=resolution.choice, 
                              verbose=FALSE)

  return(input.srobj)
}

StandardWorkflow <- function(sc.object, scale = T){
  if(scale == T){
    sc.object <- ScaleData(sc.object, verbose = FALSE, vars.to.regress = c("nCount_RNA", "nFeature_RNA"))
  }
  
  sc.object <- RunPCA(sc.object, npcs = 40, verbose = FALSE)
  print("scaling and running PCA...")
  # calculate how much variance is explained by each PC
  pct <- sc.object[["pca"]]@stdev / sum(sc.object[["pca"]]@stdev) * 100
  # calculate cumulative percent for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # returns 21
  # take the higher value so that at least 90 % of the variance will be explained.
  pcs <- max(co1, co2)
  
  print(paste("number of PCs to use for subsequent analysis: ", pcs, sep = ""))
  
  # Visualize cutoff
  # Create a dataframe with values
  plot_df <- data.frame(pct = pct, 
                        cumu = cumu, 
                        rank = 1:length(pct))
  
  # Elbow plot to visualize 
  elbow_plot <- ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
    geom_text() + 
    geom_vline(xintercept = 90, color = "black", linetype= "dashed") + 
    geom_hline(yintercept = min(pct[pct > 5]), color = "black", linetype= "dashed") +
    xlab("cumulative variance explained [%]") +
    ylab("variance explained [%]") +
    ggtitle("Elbow Plot") +
    theme_classic()
  print(elbow_plot)
  
  # continue clustering with the determined number of PCs
  print("running UMAP...")
  sc.object <- RunUMAP(sc.object, reduction = "pca", dims = 1:pcs)
  print("finding neighbors...")
  sc.object <- FindNeighbors(sc.object, reduction = "pca", dims = 1:pcs)
  return(sc.object)
}

# import pre-processed single-cell sequencing dataset:
cells <- readRDS('./Data/telencephalon_cells.rds')
# export as h5ad file:
SaveH5Seurat(cells, "scdr.h5seurat")
Convert("scdr.h5seurat", dest="h5ad") # or set assay="RNA".

# select neurons:
neurons <- subset(cells, ident = c("GABA", "Glut"))
neurons <- subset(neurons, subset = nFeature_RNA > 1000)
saveRDS(neurons, file = './Data/neurons_raw.RDS')

# run normal workflow
neurons <- StandardWorkflow(neurons, scale = T)
# use 35 PCs for cluster resolution:
neurons <- ChooseClusterResolutionDownsample(neurons, 35, 
                                                    coarseness = 3, print_fig = F, figdir="./Plots")
# Visualization

pdf("./Plots/neurons_DimPlot.pdf")
DimPlot(neurons, reduction = "umap", group.by = "orig.ident", order = c("replicate 1", "replicate 2"))
DimPlot(neurons, reduction = "umap", label = TRUE, repel = TRUE)
FeaturePlot(neurons, features = c("slc17a6a", "slc32a1"), min.cutoff = "q9")
dev.off()

##### SUBSET GLUTAMATERGIC & GABAERGIC CELLS

glut <- c(1, 2, 5, 7, 8, 9, 11, 12, 13, 16, 17, 19, 20, 25) 
gaba <- c(0, 3, 4, 6, 10, 14, 15, 18, 21, 22, 23, 24, 26, 27, 28) 

neurons_gaba <- subset(neurons, ident = gaba)
neurons_glut <- subset(neurons, ident = glut)

# filter
slc32a1_expression = GetAssayData(object = neurons_glut, assay = "RNA", slot = "data")["slc32a1",]
neg_ids = names(which(slc32a1_expression==0))
neurons_glut_filtered = subset(neurons_glut, cells = neg_ids)

slc32a1_expression = GetAssayData(object = neurons_gaba, assay = "RNA", slot = "data")["slc32a1",]
pos_ids = names(which(slc32a1_expression>0))
neurons_gaba_filtered = subset(neurons_gaba,cells=pos_ids)

# check result
FeaturePlot(neurons_gaba,"slc32a1")
FeaturePlot(neurons_gaba_filtered,"slc32a1")
FeaturePlot(neurons_glut,"slc32a1")
FeaturePlot(neurons_glut_filtered,"slc32a1")

# Run standard workflow: 
neurons_glut_filtered <- StandardWorkflow(neurons_glut_filtered, scale = T)
neurons_gaba_filtered <- StandardWorkflow(neurons_gaba_filtered, scale = T)

# use 35 PCs for cluster resolution:
neurons_glut_filtered <- ChooseClusterResolutionDownsample(neurons_glut_filtered, 35, coarseness = 1,
                                                           print_fig = F, figdir="./Plots")
neurons_gaba_filtered <- ChooseClusterResolutionDownsample(neurons_gaba_filtered, 35, coarseness = 1,
                                                           print_fig = F, figdir="./Plots")


# Visualization
pdf("./Plots/neurons_filtered_DimPlot.pdf")
DimPlot(neurons_glut_filtered, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(neurons_gaba_filtered, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()

##### Glutamatergic subclustering

# Rename according to marker genes and proximity
new.cluster.ids <- c("e5",#0
                     "e11",#1
                     "e3",#2
                     "e2",#3
                     "e13",#4
                     "e10",#5
                     "e1",#6
                     "e7",#7
                     "e14",#8
                     "e9",#9
                     "e19",#10
                     "e16",#11
                     "e8",#12
                     "e6",#13
                     "e12",#14
                     "e18",#15
                     "e4",#16
                     "e17",#17
                     "e15"#18
)

names(new.cluster.ids) <- levels(neurons_glut_filtered)
neurons_glut_filtered <- RenameIdents(neurons_glut_filtered, new.cluster.ids)
levels(x = neurons_glut_filtered) <- c("e1","e2","e3","e4","e5",
                                              "e6","e7","e8","e9","e10",
                                              "e11","e12","e13","e14","e15",
                                              "e16","e17","e18","e19")                   

glut.markers <- FindAllMarkers(neurons_glut_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
glut.markers_top25 <- glut.markers %>% group_by(cluster) %>% top_n(25, avg_log2FC)
write.csv(glut.markers_top25, "./Analysis/glut.markers_top25.csv")

##### GABAergic subclustering

new.cluster.ids <- c("i19",#0
                     "i11",#1
                     "i17",#2
                     "i30",#3
                     "i22",#4
                     "i16",#5
                     "i18",#6
                     "i10",#7
                     "i13",#8
                     "i28",#9
                     "i8",#10
                     "i12",#11
                     "i9",#12 
                     "i36",#13
                     "i26",#14 
                     "i27",#15
                     "i29",#16 
                     "i4",#17
                     "i14",#18
                     "i23",#19
                     "i7",#20
                     "i25",#21
                     "i5",#22
                     "i6",#23
                     "i3",#24
                     "i2",#25
                     "i20",#26
                     "i21",#27
                     "i15",#28
                     "i1",#29
                     "i35",#30
                     "i32",#31 
                     "i31",#32
                     "i34",#33
                     "i24",#34 
                     "i33" #35
                 
)
names(new.cluster.ids) <- levels(neurons_gaba_filtered)
neurons_gaba_filtered <- RenameIdents(neurons_gaba_filtered, new.cluster.ids)

levels(x = neurons_gaba_filtered) <- c("i1", "i2", "i3", "i4", "i5", 
                                       "i6", "i7", "i8", "i9", "i10",
                                       "i11", "i12", "i13", "i14", "i15",
                                       "i16", "i17", "i18", "i19", "i20",
                                       "i21", "i22", "i23", "i24", "i25",
                                       "i26", "i27", "i28", "i29", "i30",
                                       "i31", "i32", "i33", "i34", "i35", 
                                       "i36")         

pdf("./Plots/neurons_subclustered_DimPlot.pdf")
DimPlot(neurons_glut_filtered, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE,  label.box = FALSE) + NoLegend()
DimPlot(neurons_gaba_filtered, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE,  label.box = FALSE) + NoLegend()
dev.off()

gaba.markers <- FindAllMarkers(neurons_gaba_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gaba.markers_top25 <- gaba.markers %>% group_by(cluster) %>% top_n(25, avg_log2FC)
write.csv(gaba.markers_top25, "./Analysis/gaba.markers_top25.csv")

saveRDS(neurons_glut_filtered, file="./Data/neurons.glut.reclustered.rds")
saveRDS(neurons_gaba_filtered, file="./Data/neurons.gaba.reclustered.rds")


##### Rename cells in all neurons dataset
gaba_idents <- Idents(neurons_gaba_filtered)
glut_idents <- Idents(neurons_glut_filtered)
neuron_idents <- c(gaba_idents, glut_idents)
# Stash cell identity classes
neurons[["old.ident"]] <- Idents(object = neurons)
# subset and rename
cell_idents <- names(neuron_idents)
neurons.subset <- neurons[,cell_idents]
neurons.subset <-SetIdent(neurons.subset, cells = cell_idents, neuron_idents)

# Visualization
pdf("./Plots/neurons_subset_DimPlot.pdf")
DimPlot(neurons.subset, reduction = "umap", group.by = "orig.ident", order = c("replicate 1", "replicate 2"))
DimPlot(neurons.subset,  reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE,  label.box = FALSE) + NoLegend()
dev.off()

# export
saveRDS(neurons.subset, file = "./Data/neurons.subset.rds")

# save updated cluster names in cells dataset and export as .csv
# Stash cell identity classes
cells[["old.ident"]] <- Idents(object = cells)
scdr_full <-as.data.frame(cells[["old.ident"]])
scdr_full$supercluster <- scdr_full$old.ident
levels(scdr_full$supercluster) <- c(levels(scdr_full$supercluster), levels(neurons))
scdr_full[names(Idents(neurons)),]$supercluster <- Idents(neurons)

scdr_full$subcluster <- scdr_full$old.ident
levels(scdr_full$subcluster) <- c(levels(scdr_full$subcluster), levels(neuron_idents))
scdr_full[cell_idents,]$subcluster <- neuron_idents
# udpate cells object
cells <- SetIdent(object = cells, cells = row.names(scdr_full), value = scdr_full$subcluster)
DimPlot(cells, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE,  label.box = FALSE) + NoLegend()

# save .csv
write.csv(scdr_full, "./Analysis/scdr_full.csv")

# export as h5ad file:

SaveH5Seurat(cells, "./Analysis/scdr_v2.h5seurat")
Convert("./Analysis/scdr_v2.h5seurat", dest="h5ad") # or set assay="RNA".

# to visualize individual genes:

plot <- FeaturePlot(neurons.subset, "zic2a", min.cutoff = 0, cols = c("lightgrey", "red"), order = FALSE) + FontSize(main = 0) + NoAxes() + NoLegend() 

# remove axes and ticks, titles
plot <- plot + theme_void() + labs(title = NULL) + NoLegend() 

# save the plot with transparent background as png
png(".Plots/zic2a.png", width = 450, height = 250,  bg = "transparent")
print(plot)
dev.off()

