library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(ggpubr)
#library(limma)
library(dendextend)
library(ComplexHeatmap)
library(genefilter)

library(circlize)
col_fun = colorRamp2(c(-0.3, 0, 0.3), c("blue", "white", "red"))
col_fun <- colorRamp2(c(-0.3, 0, 0.3),c('#6E377D',"white",'#EBB400'))
col_fun(seq(-3, 3))

# Load Datasets
neurons.glut <- readRDS('./Data/neurons.glut.reclustered.rds')
neurons.gaba <- readRDS('./Data/neurons.gaba.reclustered.rds')
bulk.DEG.avg <- readRDS('./Data/NGS_DEG.rds')

# Compute average count of informative genes in single cell data
glut.cluster.averages <- AverageExpression(neurons.glut)
gaba.cluster.averages <- AverageExpression(neurons.gaba)

# find genes that are differentially expressed between regions based on bulk sequencing data
common_geneset <- intersect(row.names(glut.cluster.averages$integrated), row.names(bulk.DEG.avg))
sc_corr_matrix <- glut.cluster.averages$integrated[common_geneset,]
bulk_corr_matrix <- bulk.DEG.avg[common_geneset,]
# compute row-wise z-score to get standardized expression level
bcm_zs <- (bulk_corr_matrix - rowMeans(bulk_corr_matrix)) / rowSds(bulk_corr_matrix)
scm_zs <- (sc_corr_matrix - rowMeans(sc_corr_matrix)) / rowSds(sc_corr_matrix)
# purge rows with NA:
na_id <- which(is.na(scm_zs), arr.ind = TRUE)
row_id <- unique(na_id[,1])
bcm_zs <- bcm_zs[-row_id,]
scm_zs <- scm_zs[-row_id,]

# compute p-values for all correlations:
p_glut <- as.data.frame(matrix(data = NA,
                               nrow = 4,
                               ncol = ncol(scm_zs)))
colnames(p_glut) <- colnames(scm_zs)
rownames(p_glut) <- colnames(bcm_zs)

# fill matrix step-wise
for(celltype in colnames(p_glut)){
  for(region in rownames(p_glut)){
    pv <- cor.test(bcm_zs[,region], scm_zs[, celltype])
    p_glut[region, celltype] <- pv$p.value
  }
}
# Bonferroni Correction for multiple testing
p_glut <- p_glut * length(rownames(p_glut)) * length(colnames(p_glut))

# order for cols and rows
Rowv <- bcm_zs %>% t %>% dist %>% hclust %>% as.dendrogram() %>% ladderize()
Colv <- scm_zs %>% t %>% dist %>% hclust %>% as.dendrogram() %>% ladderize()

column_order = paste0("e", seq(1, 19))

m_glut <- t(cor(scm_zs, bcm_zs))
ht_opt$TITLE_PADDING = unit(c(8, 8), "points")
hm_glut <- Heatmap(m_glut, name = "Pearson", column_title_side = "top",
                   column_title_gp = gpar(fill = "black", col = "white", border = "white"),
                   column_title = "GLUTAMATERGIC", col = col_fun, cluster_rows = Rowv,
                   #cluster_columns = Colv,
                   column_order = column_order, column_names_rot = 45, 
                   row_split = 4, row_gap = unit(2, "mm"), border = TRUE,
                   left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 2:4),
                                                                   labels = c("Dp", "Dl", "Dm", "Dc"), 
                                                         labels_gp = gpar(col = "white", fontsize = 10))),
                   # add significance levels
                   cell_fun = function(j, i, x, y, w, h, fill) {
                     if(p_glut[i, j] < 0.001 && m_glut[i, j]>0){
                       grid.text("***", x, y)
                     } else if(p_glut[i, j] < 0.01 && m_glut[i, j]>0) {
                       grid.text("**", x, y)
                     } else if(p_glut[i, j] < 0.05 && m_glut[i, j]>0) {
                       grid.text("*", x, y)
                     }}
                   )



# Same analysis with GABAergic cells should result in basically no correlation:
#gaba.cluster.averages <- AverageExpression(scData.neurons.gaba)
# find genes that are differentially expressed between regions based on bulk sequencing data
common_geneset <- intersect(row.names(gaba.cluster.averages$integrated), row.names(bulk.DEG.avg))
sc_corr_matrix <- gaba.cluster.averages$integrated[common_geneset,]
bulk_corr_matrix <- bulk.DEG.avg[common_geneset,]
# compute row-wise z-score to get standardized expression level
bcm_zs <- (bulk_corr_matrix - rowMeans(bulk_corr_matrix)) / rowSds(bulk_corr_matrix)
scm_zs <- (sc_corr_matrix - rowMeans(sc_corr_matrix)) / rowSds(sc_corr_matrix)
# purge rows with NA:
na_id <- which(is.na(scm_zs), arr.ind = TRUE)
row_id <- unique(na_id[,1])
bcm_zs <- bcm_zs[-row_id,]
scm_zs <- scm_zs[-row_id,]

# compute p-values for all correlations:
p_gaba <- as.data.frame(matrix(data = NA,
                               nrow = 6,
                               ncol = ncol(scm_zs)))
colnames(p_gaba) <- colnames(scm_zs)
rownames(p_gaba) <- colnames(bcm_zs)

# fill matrix step-wise
for(celltype in colnames(p_gaba)){
  for(region in rownames(p_gaba)){
    pv <- cor.test(bcm_zs[,region], scm_zs[, celltype])
    p_gaba[region, celltype] <- pv$p.value
  }
}

p_gaba <- p_gaba * length(rownames(p_gaba)) * length(colnames(p_gaba))

# compute dendrogram for the GABAergic clusters
Colv <- scm_zs %>% t %>% dist %>% hclust %>% as.dendrogram() %>% ladderize()
m_gaba <- t(cor(scm_zs, bcm_zs))

hm_gaba <- Heatmap(m_gaba, name = "GABAERGIC", column_title_side = "top",
                   column_title = "GABAERGIC", col = col_fun, cluster_rows = Rowv,
                   column_title_gp = gpar(fill = "black", col = "white", border = "white"),
                   cluster_columns = Colv, column_names_rot = 45, 
                   row_split = 3, row_gap = unit(2, "mm"), border = TRUE,
                   show_heatmap_legend = FALSE,
                   # add significance levels
                   cell_fun = function(j, i, x, y, w, h, fill) {
                     if(p_gaba[i, j] < 0.001 && m_gaba[i, j]>0) {
                       grid.text("***", x, y)
                     } else if(p_gaba[i, j] < 0.01 && m_gaba[i, j]>0) {
                       grid.text("**", x, y)
                     } else if(p_gaba[i, j] < 0.05 && m_gaba[i, j]>0) {
                       grid.text("*", x, y)
                     }})
                   
                  # left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 2:4),
                  #                                                  labels = c("DM", "Dp | Dl", "Dc"), 
                  #                                                  labels_gp = gpar(col = "white", fontsize = 10))))
pdf("./Plots/CorrelationBulkSingleCell.pdf")
hm_glut + hm_gaba
dev.off()
