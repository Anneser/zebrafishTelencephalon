library(Seurat)
library(dplyr)
library(cowplot)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(clustree)

#install.packages('BiocManager')
#BiocManager::install('SeuratDisk')

# define functions

loadscRNAdataset <- function(data_path, n_feature = 500, mito_cutoff = 0.1, project = "NA"){
  if (!require(Seurat)) library(Seurat)
  setwd(data_path)
  available_files <- list.files(data_path) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
  necessary_files <- c("barcodes.tsv.gz" , "features.tsv.gz", "matrix.mtx.gz")
  if (identical(sum(sapply(necessary_files, function(i) sum(grepl(i, available_files)))), as.integer(3))){
    print("all necessary files exist in directory...")
  }
  
  data <- Read10X(data.dir = data_path)
  print("loading files...")
  nUMI = colSums(data)
  
  print(paste("median of nUMIs per cell: ", median(nUMI), sep=""))
  names(nUMI) = colnames(data)
  
  percent.mito <- colSums(data[grep("^mt-",rownames(data)),])/nUMI
  print(paste("median of mitochondrial genes per cell (fraction): ", round(median(percent.mito), digits = 3), sep=""))
  
  SEO = CreateSeuratObject(counts = data , 
                           min.cells = 5 ,
                           min.features = 200 ,
                           project = project)
  SEO <- AddMetaData(object = SEO, 
                     metadata = nUMI, 
                     col.name = "nUMI")
  
  SEO <- AddMetaData(object = SEO, 
                     metadata = percent.mito, 
                     col.name = "percent.mito")
  
  
  # subset dataset for QC
  SEO <- subset(SEO , subset = nFeature_RNA > 500  & 
                  percent.mito < 0.1 )
  print(paste("Quality Control - filtering cells with less than", n_feature, " features and more than ",
              mito_cutoff, " mitochondrial genes.", sep=""))
  print(paste(length(SEO@assays$RNA@counts@p), " cells remain in dataset, corresponding to ", 
              round(length(SEO@assays$RNA@counts@p)/length(nUMI)*100, digits = 3), " %.", sep=""))
  
  return(SEO)
}

StandardWorkflow <- function(sc.object){
  
  sc.object <- ScaleData(sc.object, verbose = FALSE, vars.to.regress = c("nCount_RNA", "nFeature_RNA")) # vars.to.regress
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

FindSignificantPCs <- function(sc.object, n_shuffles = 10, PCs = 50){
  sc.temp <- sc.object
  DefaultAssay(sc.temp) <- "integrated"
  atts <- attributes(sc.temp@assays$integrated@scale.data)
  sc.temp <- RunPCA(sc.temp, verbose = FALSE)
  VarExpl <- sc.temp[["pca"]]@stdev
  print("begin shuffling...")
  
  # pre-allocate empty dataframe
  df_PCs <- data.frame(matrix(ncol = n_shuffles, nrow = PCs))
  
  
  for(i in 1:n_shuffles) {
    print(paste("shuffle #", i, sep = ""))
    for(j in 1:length(atts$dimnames[[1]])){
      sc.temp@assays$integrated@scale.data[j,] <- sample(sc.temp@assays$integrated@scale.data[j,])
    }
    sc.temp <- RunPCA(sc.temp, verbose = FALSE)
    df_PCs[,i] <- c(sc.temp[["pca"]]@stdev)
    
  }
  
  return(df_PCs)
}

ChooseClusterResolutionDownsample <- function(
    input.srobj, n.pcs, sample.name =  format(paste(Sys.Date(), "SilhouetteScore", sep="_")),
    res.low = .01, res.high=3, res.n = 40, bias = "over", print_fig = F, figdir=F) {
  
  ######## step 1: save the input seurat object as a new temporary object, 
  ########  we dont want to overwrite or change the original one with all of the parameter scans
  
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
  ########  average has worked well so far... could be different
  
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
                 " and will return your input seurat object amended with the best choice",
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
  
  Best.Clusters = srobj.tmp@meta.data[,paste("integrated_snn_res.", resolution.choice, sep="")]
  
  input.srobj$Best.Clusters = Best.Clusters
  Idents(input.srobj) = input.srobj$Best.Clusters
  input.srobj@misc$resolution.choice <- resolution.choice
  return(input.srobj)
}

# load Datasets

SEO <- loadscRNAdataset(".\20220125", project = "replicate 1")
SEO2 <- loadscRNAdataset(".\20220322", project = "replicate 2")

# merge datasets
all_df <- merge(x = SEO, y = SEO2)

# export as h5ad file:
SaveH5Seurat(all_df, "scdr_v2.h5seurat")
Convert("scdr_v2.h5seurat", dest="h5ad") # or set assay="RNA".

scData.list <- SplitObject(all_df, split.by = "orig.ident")

# clean up
rm(SEO, SEO2, all_df)

#normalize and identify variable features for each dataset independently
scData.list <- lapply(X = scData.list, FUN = function(x) {
  #x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 1e4)
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = scData.list)
scData.anchors <- FindIntegrationAnchors(object.list = scData.list, anchor.features = features)

# this command creates an integrated data assay
scData.combined <- IntegrateData(anchorset = scData.anchors)

# specify that we will perform downstream analysis on the corrected data, note that the original
# unmodified data still resides in the "RNA" assay
DefaultAssay(scData.combined) <- "integrated"

scData.combined <- StandardWorkflow(scData.combined)

scData.combined <- ChooseClusterResolutionDownsample(scData.combined, 33)

# Visualization
pdf("./Plots/cells_DimPlot1.pdf")
DimPlot(scData.combined, reduction = "umap", group.by = "orig.ident", order = c("replicate 1", "replicate 2"))
dev.off()
pdf("./Plots/cells_DimPlot2.pdf")
DimPlot(scData.combined, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()

cells.markers <- FindAllMarkers(object = scData.combined, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
cells.markers_top25 <- cells.markers %>% group_by(cluster) %>% top_n(25, avg_log2FC)
write.csv(cells.markers_top25, file = "./Analysis/cells.markers_top25.csv")

# Visualize marker genes
pdf("./Plots/cells_FeaturePlots.pdf")
FeaturePlot(object = scData.combined, features = c("snap25a","slc17a7a","slc17a6a","slc32a1"), min.cutoff = "q9")
FeaturePlot(object = scData.combined, features = c("olig1","gpr17", "gfap","mki67"), min.cutoff = "q9")
FeaturePlot(object = scData.combined, features = c("olig2","mbpa","sox10", "kcnj11l"), min.cutoff = "q9")
FeaturePlot(object = scData.combined, features = c("cldnk","mpz","cd59"), min.cutoff = "q9")
FeaturePlot(object = scData.combined, features = c("cspg4","her4.1","pcna", "mki67"), min.cutoff = "q9")
FeaturePlot(object = scData.combined, features = c("ccl36.1","il2rb", "c1qa", "c1qb"), min.cutoff = "q9")
FeaturePlot(object = scData.combined, features = c("c1qc","marco","tbx18", "pfn1"), min.cutoff = "q9")
FeaturePlot(object = scData.combined, features = c("kiss1","gng8", "dcn", "rbp4"), min.cutoff = "q9")
FeaturePlot(object = scData.combined, features = c("slc2a1a","klf2a", "angptl3", "exorh"), min.cutoff = "q9")
dev.off()

saveRDS(scData.combined, file="./Data/scData_combined.rds")

new.cluster.ids.l1 <- c("Glut", #0
                        "ependymoglia", # 1
                        "GABA", # 2
                        "Glut", # 3
                        "GABA", # 4
                        "Glut", # 5
                        "Glut", # 6
                        "GABA", # 7
                        "OPC",  # 8
                        "GABA", # 9
                        "Glut", # 10
                        "oligodendrocytes", # 11
                        "NPC", # 12
                        "leukocytes", # 13
                        "GABA", # 14
                        "microglia", # 15
                        "GABA", # 16
                        "cycling cells", # 17
                        "habenula", # 18
                        "differentiating OPC", # 19
                        "microglia", # 20
                        "End", # 21
                        "vEnd", # 22
                        "vSMC", # 23
                        "epiphysis", # 24
                        "ependymoglia", # 25
                        "habenula") # 26
                        
new.cluster.ids.l2 <- c("Glut1", #0
                        "ependymoglia1", # 1
                        "GABA1", # 2
                        "Glut2", # 3
                        "GABA2", # 4
                        "Glut3", # 5
                        "Glut4", # 6
                        "GABA3", # 7
                        "OPC",  # 8
                        "GABA4", # 9
                        "Glut5", # 10
                        "oligodendrocytes", # 11
                        "NPC", # 12
                        "leukocytes", # 13
                        "GABA", # 14
                        "microglia1", # 15
                        "GABA", # 16
                        "cycling cells", # 17
                        "habenula1", # 18
                        "differentiating OPC", # 19
                        "microglia2", # 20
                        "End", # 21 
                        "vEnd", # 22
                        "vSMC", # 23
                        "epiphysis", # 24
                        "ependymoglia2", # 25
                        "habenula2") # 26
                                      

names(new.cluster.ids.l1) <- levels(scData.combined)
names(new.cluster.ids.l2) <- levels(scData.combined)

scData.combined <- RenameIdents(scData.combined, new.cluster.ids.l1)

levels(x = scData.combined) <- c("Glut", 
                                 "GABA",
                                 "NPC",
                                 "ependymoglia", 
                                 "cycling cells",
                                 "OPC", 
                                 "differentiating OPC", 
                                 "oligodendrocytes", 
                                 "microglia",
                                 "leukocytes",
                                 "End", 
                                 "vEnd",
                                 "vSMC",
                                 "epiphysis",
                                 "habenula")

pdf("./Plots/cells_DimPlot_Level1.pdf")
DimPlot(scData.combined, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()

saveRDS(scData.combined, file="./Data/telencephalon_cells.rds")
scdr_v2 <- as.data.frame(cells@active.ident)

# DotPlot of marker gene expression
pdf("./Plots/cells_DotPlot_Level1.pdf")
DotPlot(scData.combined, assay = "RNA", features = c("snap25a",
                                      "slc17a6a",
                                      "slc32a1",
                                      "sox12",
                                      "gfap",
                                      "mki67",
                                      "olig1",
                                      "gpr17",
                                      "mbpa",
                                      "marco",
                                      "il2rb",
                                      "dcn",
                                      "slc2a1a",
                                      "angptl3",
                                      "exorh",
                                      "kiss1" )) + RotatedAxis()
dev.off()

# which clusters are neurons (snap25a expression) --> subcluster for better resolution
neurons = c(0,1,2,4,5,6,8,9,15)

scData.combined_neurons_ID <- scData.combined@meta.data$Best.Clusters %in% neurons
scData.neurons <- scData.combined[,scData.combined_neurons_ID]

save(scData.neurons, file = "./Data/scDataNeurons.RData")