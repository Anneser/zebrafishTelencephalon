# import libraries

library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(clustree)
library(dendextend)
library(randomForest)
library(biomaRt)
library(circlize)
library(rjson)
library(ComplexHeatmap)
library(readxl)
library(biomaRt)

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

# expression of TF, ribosomal genes, GPCRs, cell adhesion molecules, proteoglycans, zinc fingers, ion channels in neuron types.
# download KEGG lists from: https://www.kegg.jp/brite/dre03310.keg <-- adapt number as appropriate
# import locally compiled kegg list
kegg <- read.csv("./Data/KEGG.csv")
gene_names <- neurons.subset@assays$RNA@counts@Dimnames[1][[1]]

ensembl <- useEnsembl( "ensembl", dataset = "drerio_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", 
                                 "entrezgene_id", 
                                 "external_gene_name", 
                                 "refseq_mrna"
                                 ),
                  filters = "external_gene_name",
                  values = gene_names,
                  mart = ensembl )
idx <- match(gene_names , genemap$external_gene_name)
gene_list <- as.data.frame(gene_names)
gene_list$entrez <- genemap$entrezgene_id[ idx ]
gene_list$ensembl <- genemap$ensembl_gene_id[ idx ]

idx <- match(gene_list$entrez, kegg$entrezgene_id)
gene_list$KEGG_family <- kegg$KEGG_family[ idx ]
gene_list$KEGG_node <- kegg$KEGG_node[ idx]
gene_list$KEGG_subnode <- kegg$KEGG_subnode[ idx ]

annotated_gene_list <- gene_list[is.na(gene_list$KEGG_node) == 0,]

########### Cumulative Gene Expression Analysis

neuron_expression <- AverageExpression(neurons.subset)
neuron_expression_10 <- neuron_expression
cell_id <- neurons.subset@active.ident
for (i in 1:length(colnames(neuron_expression_10$RNA))) {
  active_id <- colnames(neuron_expression_10$RNA)[[i]]
  active_cell_id <- which(cell_id == active_id)
  active_cells <- neurons.subset@assays$RNA[,active_cell_id]
  neuron_expression_10$RNA[,active_id] <- rowSums(active_cells > 0)/length(active_cell_id) >= 0.1
}

# toggle for different modes of calculation:
neuron_expression <- neuron_expression_10

# transcription factors
tf_id = annotated_gene_list$KEGG_family == "dre03000"
tf_list = annotated_gene_list$gene_names[tf_id]
tf_expression <- neuron_expression$RNA[tf_list,]
# get rid of zero-rows
zero_id = which(rowSums(tf_expression)!= 0 )
tf_expression = tf_expression[zero_id,]
cumsum_tf = sort(colSums(tf_expression != 0))/max(colSums(tf_expression != 0))

# homeodomain transcription factors
homeodomain_list = subset(annotated_gene_list[annotated_gene_list$KEGG_family == "dre03000",]$KEGG_subnode, grepl('^Homeo domain', annotated_gene_list[annotated_gene_list$KEGG_family == "dre03000",]$KEGG_subnode))
homeodomain_id = annotated_gene_list$KEGG_subnode %in% homeodomain_list
homeodomain_list = annotated_gene_list$gene_names[homeodomain_id]
homeodomain_expression <- neuron_expression$RNA[homeodomain_list,]
zero_id = which(rowSums(homeodomain_expression)!= 0 )
homeodomain_expression = homeodomain_expression[zero_id,]
cumsum_homeodomain = sort(colSums(homeodomain_expression != 0))/max(colSums(homeodomain_expression != 0))

# Zinc fingers
zf_id = annotated_gene_list$KEGG_node == "Zinc finger"
zf_list <- annotated_gene_list$gene_names[zf_id]
zf_expression <- neuron_expression$RNA[zf_list,]
zero_id = which(rowSums(zf_expression)!=0)
zf_expression <- zf_expression[zero_id,]
cumsum_zf = sort(colSums(zf_expression != 0)) / max(colSums(zf_expression!=0))

# ribosomal proteins
ribosome_id = annotated_gene_list$KEGG_node == "Ribosomal proteins"
ribosome_list = annotated_gene_list$gene_names[ribosome_id]
ribosome_expression <- neuron_expression$RNA[ribosome_list,]
zero_id = which(rowSums(ribosome_expression)!= 0 )
ribosome_expression = ribosome_expression[zero_id,]
cumsum_ribosome = sort(colSums(ribosome_expression != 0))/max(colSums(ribosome_expression != 0))

# Proteoglycans
pg_id <- annotated_gene_list$KEGG_family == "dre00535"
pg_list <- annotated_gene_list$gene_names[pg_id]
pg_expression <- neuron_expression$RNA[pg_list,]
zero_id <- which(rowSums(pg_expression)!= 0)
pg_expression <- pg_expression[zero_id,]
cumsum_pg <- sort(colSums(pg_expression != 0)) / max(colSums(pg_expression != 0))

# Ion Channels
ic_id <- annotated_gene_list$KEGG_family == "dre04040"
ic_list <- annotated_gene_list$gene_names[ic_id]
ic_expression <- neuron_expression$RNA[ic_list,]
zero_id <- which(rowSums(ic_expression)!= 0)
ic_expression <- ic_expression[zero_id,]
cumsum_ic <- sort(colSums(ic_expression != 0)) / max(colSums(ic_expression != 0))

# Cell Adhesion Molecules 
cam_id <- annotated_gene_list$KEGG_family == "dre04515"
cam_list <- annotated_gene_list$gene_names[cam_id]
cam_expression <- neuron_expression$RNA[cam_list,]
zero_id <- which(rowSums(cam_expression)!= 0)
cam_expression <- cam_expression[zero_id,]
cumsum_cam <- sort(colSums(cam_expression != 0)) / max(colSums(cam_expression != 0))

# GPCRs
gpcr_id <- annotated_gene_list$KEGG_family == "dre04030"
gpcr_list <- annotated_gene_list$gene_names[gpcr_id]
gpcr_expression <- neuron_expression$RNA[gpcr_list,]
zero_id <- which(rowSums(gpcr_expression)!= 0)
gpcr_expression <- gpcr_expression[zero_id,]
cumsum_gpcr <- sort(colSums(gpcr_expression != 0)) / max(colSums(gpcr_expression != 0))

# plot cumulative expression across cell types of all marker gene categories
cumsum_df <- data.frame(matrix(ncol = 3, nrow = length(cumsum_tf)))
colnames(cumsum_df) <- c("celltype", "value", "gene_category")
cumsum_df$value <- cumsum_tf
cumsum_df$celltype <- c(1:55)
cumsum_df$gene_category <- rep("All TFs", 55)

# append zinc finger
helper_df <- data.frame(matrix(ncol = 3, nrow = length(cumsum_tf)))
colnames(helper_df) <- c("celltype", "value", "gene_category")
helper_df$value <- cumsum_zf
helper_df$celltype <- c(1:55)
helper_df$gene_category <- rep("Zinc finger", 55)
cumsum_df <- rbind(cumsum_df, helper_df)
# append homeodomain
helper_df <- data.frame(matrix(ncol = 3, nrow = length(cumsum_tf)))
colnames(helper_df) <- c("celltype", "value", "gene_category")
helper_df$value <- cumsum_homeodomain
helper_df$celltype <- c(1:55)
helper_df$gene_category <- rep("Homeodomain", 55)
cumsum_df <- rbind(cumsum_df, helper_df)
# append ribosome
helper_df <- data.frame(matrix(ncol = 3, nrow = length(cumsum_tf)))
colnames(helper_df) <- c("celltype", "value", "gene_category")
helper_df$value <- cumsum_ribosome
helper_df$celltype <- c(1:55)
helper_df$gene_category <- rep("Ribosomal proteins", 55)
cumsum_df <- rbind(cumsum_df, helper_df)
# append proteoglycans
helper_df <- data.frame(matrix(ncol = 3, nrow = length(cumsum_pg)))
colnames(helper_df) <- c("celltype", "value", "gene_category")
helper_df$value <- cumsum_pg
helper_df$celltype <- c(1:length(cumsum_pg))
helper_df$gene_category <- rep("Proteoglycans", length(cumsum_pg))
cumsum_df <- rbind(cumsum_df, helper_df)
# append ion channels 
helper_df <- data.frame(matrix(ncol = 3, nrow = length(cumsum_pg)))
colnames(helper_df) <- c("celltype", "value", "gene_category")
helper_df$value <- cumsum_ic
helper_df$celltype <- c(1:length(cumsum_ic))
helper_df$gene_category <- rep("Ion Channels", length(cumsum_ic))
cumsum_df <- rbind(cumsum_df, helper_df)
# append cell adhesion molecules 
helper_df <- data.frame(matrix(ncol = 3, nrow = length(cumsum_pg)))
colnames(helper_df) <- c("celltype", "value", "gene_category")
helper_df$value <- cumsum_cam
helper_df$celltype <- c(1:length(cumsum_cam))
helper_df$gene_category <- rep("Cell adhesion molecules", length(cumsum_cam))
cumsum_df <- rbind(cumsum_df, helper_df)
# append GPCRs 
helper_df <- data.frame(matrix(ncol = 3, nrow = length(cumsum_pg)))
colnames(helper_df) <- c("celltype", "value", "gene_category")
helper_df$value <- cumsum_gpcr
helper_df$celltype <- c(1:length(cumsum_gpcr))
helper_df$gene_category <- rep("GPCRs", length(cumsum_gpcr))
cumsum_df <- rbind(cumsum_df, helper_df)

cumsum_df$gene_category <- factor(cumsum_df$gene_category, levels=c("All TFs",
                                                                    "Homeodomain",
                                                                    "Zinc finger",
                                                                    "Cell adhesion molecules",
                                                                    "Ion Channels",
                                                                    "Proteoglycans",
                                                                    "GPCRs",
                                                                    "Ribosomal proteins"))


ggplot(cumsum_df, aes(celltype, value, colour = gene_category)) +
  geom_line(size = 1)+
  geom_point(size = 2)+
  ylim(c(0,1)) +
  ylab("cumulative fraction of expressed genes per family") +
  xlab("number of cell types expressing cumulative fraction") +
  scale_colour_brewer("Gene Family", type = "seq", palette = "Set2") +
  theme_bw() + 
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.background = element_rect(fill = "white", size = 4, colour = "white"),
    legend.justification = c(0, 1),
    #legend.position = c(0, 1),
    legend.position = "none",
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    axis.text=element_text(size=12), #change font size of axis text
    axis.title=element_text(size=12),
    panel.grid.major = element_line(colour = "grey70", size = 0.2),
    panel.grid.minor = element_blank()
  )
#########

############ Compute Cumulative Expression Levels per Cell Type (as implemented in Hain et al., 2022)
object.raw.data <- as.matrix(neurons.subset@assays$RNA@counts)
pct.matrix = matrix(data=NA, nrow=length(gene_names), ncol=length(level_list))
rownames(pct.matrix) <- gene_names
colnames(pct.matrix) <- level_list
thresh.min=0
for (i in level_list){
  cells.cluster <- WhichCells(object=neurons.subset, idents=i)
  data.cluster <- object.raw.data[,colnames(object.raw.data) %in% cells.cluster]
  pct.cluster <- round(apply(object.raw.data[gene_names, cells.cluster, drop = F],1,function(x)return(length(x[x>thresh.min])/length(x))),3)
  pct.matrix[,i] <- pct.cluster
}
pct.df <- as.data.frame(pct.matrix)
saveRDS(pct.df, "./Data/neurons_genes_percent_expressed.rds")

# transcription factors
tf_id = annotated_gene_list$KEGG_family == "dre03000"
tf_list = annotated_gene_list$gene_names[tf_id]
# homeodomain transcription factors
homeodomain_list = subset(annotated_gene_list[annotated_gene_list$KEGG_family == "dre03000",]$KEGG_subnode, grepl('^Homeo domain', annotated_gene_list[annotated_gene_list$KEGG_family == "dre03000",]$KEGG_subnode))
homeodomain_id = annotated_gene_list$KEGG_subnode %in% homeodomain_list
homeodomain_list = annotated_gene_list$gene_names[homeodomain_id]
# Zinc fingers
zf_id = annotated_gene_list$KEGG_node == "Zinc finger"
zf_list <- annotated_gene_list$gene_names[zf_id]
# ribosomal proteins
ribosome_id = annotated_gene_list$KEGG_node == "Ribosomal proteins"
ribosome_list = annotated_gene_list$gene_names[ribosome_id]
# Proteoglycans
pg_id <- annotated_gene_list$KEGG_family == "dre00535"
pg_list <- annotated_gene_list$gene_names[pg_id]
# Ion Channels
ic_id <- annotated_gene_list$KEGG_family == "dre04040"
ic_list <- annotated_gene_list$gene_names[ic_id]
# Cell Adhesion Molecules 
cam_id <- annotated_gene_list$KEGG_family == "dre04515"
cam_list <- annotated_gene_list$gene_names[cam_id]
# GPCRs
gpcr_id <- annotated_gene_list$KEGG_family == "dre04030"
gpcr_list <- annotated_gene_list$gene_names[gpcr_id]

#transcription factors
TF_intersection <- intersect(rownames(pct.df),tf_list)
pct.df_TFs <- pct.df[rownames(pct.df) %in% TF_intersection,]

TFs_cluster_counts <- as.data.frame(matrix(data=NA,nrow=length(rownames(pct.df_TFs)),ncol=2))
names(TFs_cluster_counts) <- c("gene","cluster_count")

for (i in 1:length(rownames(pct.df_TFs))){
  gene_row <- pct.df_TFs[i,]
  cluster_count <- length(gene_row[gene_row > 0.25])
  TFs_cluster_counts$gene[i] <- rownames(gene_row)
  TFs_cluster_counts$cluster_count[i] <- cluster_count
}

# homeodomain transcription factors
HD_intersection <- intersect(rownames(pct.df),homeodomain_list)
pct.df_HDs <- pct.df[rownames(pct.df) %in% HD_intersection,]

HDs_cluster_counts <- as.data.frame(matrix(data=NA,nrow=length(rownames(pct.df_HDs)),ncol=2))

names(HDs_cluster_counts) <- c("gene","cluster_count")

for (i in 1:length(rownames(pct.df_HDs))){
  gene_row <- pct.df_HDs[i,]
  cluster_count <- length(gene_row[gene_row > 0.25])
  HDs_cluster_counts$gene[i] <- rownames(gene_row)
  HDs_cluster_counts$cluster_count[i] <- cluster_count
}
# Zinc fingers
zf_intersection <- intersect(rownames(pct.df),zf_list)
pct.df_zfs <- pct.df[rownames(pct.df) %in% zf_intersection,]

zfs_cluster_counts <- as.data.frame(matrix(data=NA,nrow=length(rownames(pct.df_zfs)),ncol=2))

names(zfs_cluster_counts) <- c("gene","cluster_count")

for (i in 1:length(rownames(pct.df_zfs))){
  gene_row <- pct.df_zfs[i,]
  cluster_count <- length(gene_row[gene_row > 0.25])
  zfs_cluster_counts$gene[i] <- rownames(gene_row)
  zfs_cluster_counts$cluster_count[i] <- cluster_count
}
# ribosomal proteins
ribosome_intersection <- intersect(rownames(pct.df),ribosome_list)
pct.df_ribosome <- pct.df[rownames(pct.df) %in% ribosome_intersection,]

ribosome_cluster_counts <- as.data.frame(matrix(data=NA,nrow=length(rownames(pct.df_ribosome)),ncol=2))

names(ribosome_cluster_counts) <- c("gene","cluster_count")

for (i in 1:length(rownames(pct.df_ribosome))){
  gene_row <- pct.df_ribosome[i,]
  cluster_count <- length(gene_row[gene_row > 0.25])
  ribosome_cluster_counts$gene[i] <- rownames(gene_row)
  ribosome_cluster_counts$cluster_count[i] <- cluster_count
}
# Proteoglycans
pg_intersection <- intersect(rownames(pct.df),pg_list)
pct.df_pg <- pct.df[rownames(pct.df) %in% pg_intersection,]

pg_cluster_counts <- as.data.frame(matrix(data=NA,nrow=length(rownames(pct.df_pg)),ncol=2))

names(pg_cluster_counts) <- c("gene","cluster_count")

for (i in 1:length(rownames(pct.df_pg))){
  gene_row <- pct.df_pg[i,]
  cluster_count <- length(gene_row[gene_row > 0.25])
  pg_cluster_counts$gene[i] <- rownames(gene_row)
  pg_cluster_counts$cluster_count[i] <- cluster_count
}
# Ion Channels
ic_intersection <- intersect(rownames(pct.df),ic_list)
pct.df_ic <- pct.df[rownames(pct.df) %in% ic_intersection,]

ic_cluster_counts <- as.data.frame(matrix(data=NA,nrow=length(rownames(pct.df_ic)),ncol=2))

names(ic_cluster_counts) <- c("gene","cluster_count")

for (i in 1:length(rownames(pct.df_ic))){
  gene_row <- pct.df_ic[i,]
  cluster_count <- length(gene_row[gene_row > 0.25])
  ic_cluster_counts$gene[i] <- rownames(gene_row)
  ic_cluster_counts$cluster_count[i] <- cluster_count
}
# Cell Adhesion Molecules 
cam_intersection <- intersect(rownames(pct.df),cam_list)
pct.df_cam <- pct.df[rownames(pct.df) %in% cam_intersection,]

cam_cluster_counts <- as.data.frame(matrix(data=NA,nrow=length(rownames(pct.df_cam)),ncol=2))

names(cam_cluster_counts) <- c("gene","cluster_count")

for (i in 1:length(rownames(pct.df_cam))){
  gene_row <- pct.df_cam[i,]
  cluster_count <- length(gene_row[gene_row > 0.25])
  cam_cluster_counts$gene[i] <- rownames(gene_row)
  cam_cluster_counts$cluster_count[i] <- cluster_count
}
# GPCRs
gpcr_intersection <- intersect(rownames(pct.df),gpcr_list)
pct.df_gpcr <- pct.df[rownames(pct.df) %in% gpcr_intersection,]

gpcr_cluster_counts <- as.data.frame(matrix(data=NA,nrow=length(rownames(pct.df_gpcr)),ncol=2))

names(gpcr_cluster_counts) <- c("gene","cluster_count")

for (i in 1:length(rownames(pct.df_gpcr))){
  gene_row <- pct.df_gpcr[i,]
  cluster_count <- length(gene_row[gene_row > 0.25])
  gpcr_cluster_counts$gene[i] <- rownames(gene_row)
  gpcr_cluster_counts$cluster_count[i] <- cluster_count
}

# count data
countsZF <- zfs_cluster_counts[order(zfs_cluster_counts$cluster_count, decreasing = F),]
countsZF <- countsZF[countsZF$cluster_count>0,]
tZF <- table(countsZF$cluster_count)
xZF <- as.numeric(names(tZF))
yZF <- as.numeric(tZF)

countsHD <- HDs_cluster_counts[order(HDs_cluster_counts$cluster_count, decreasing = F),]
countsHD <- countsHD[countsHD$cluster_count>0,]
tHD <- table(countsHD$cluster_count)
xHD <- as.numeric(names(tHD))
yHD <- as.numeric(tHD)

counts <- TFs_cluster_counts[order(TFs_cluster_counts$cluster_count, decreasing = F),]
countsTFs <- counts[counts$cluster_count>0,]
tTFs <- table(countsTFs$cluster_count)
xTFs <- as.numeric(names(tTFs))
yTFs <- as.numeric(tTFs)

counts_ICs <- ic_cluster_counts[order(ic_cluster_counts$cluster_count, decreasing = F),]
countsICs <- counts_ICs[counts_ICs$cluster_count>0,]
tICs <- table(countsICs$cluster_count)
xICs <- as.numeric(names(tICs))
yICs <- as.numeric(tICs)


counts_GPCRs <- gpcr_cluster_counts[order(gpcr_cluster_counts$cluster_count, decreasing = F),]
countsGPCRs <- counts_GPCRs[counts_GPCRs$cluster_count>0,]
tGPCRs <- table(countsGPCRs$cluster_count)
xGPCRs <- as.numeric(names(tGPCRs))
yGPCRs <- as.numeric(tGPCRs)

counts_CAMs <- cam_cluster_counts[order(cam_cluster_counts$cluster_count, decreasing = F),]
countsCAMs <- counts_CAMs[counts_CAMs$cluster_count>0,]
tCAMs <- table(countsCAMs$cluster_count)
xCAMs <- as.numeric(names(tCAMs))
yCAMs <- as.numeric(tCAMs)

counts_PGs <- pg_cluster_counts[order(pg_cluster_counts$cluster_count, decreasing = F),]
countsPGs <- counts_PGs[counts_PGs$cluster_count>0,]
tPGs <- table(countsPGs$cluster_count)
xPGs <- as.numeric(names(tPGs))
yPGs <- as.numeric(tPGs)

counts_ribos <- ribosome_cluster_counts[order(ribosome_cluster_counts$cluster_count, decreasing = F),]
countsribos <- counts_ribos[counts_ribos$cluster_count>0,]
tRibos <- table(countsribos$cluster_count)
xRibos <- as.numeric(names(tRibos))
yRibos <- as.numeric(tRibos)

pdf("./Plots/TF_cluster_counts_by_family.pdf")
plot(xTFs,(cumsum(yTFs)/max(cumsum(yTFs))), cex=1, type="o", pch=19, xlab = "Number of neuron types expressing", ylab = "Cumulative fraction",xlim=c(0,55),ylim=c(0,1), col="#66C2A5")
lines(xZF,(cumsum(yZF)/max(cumsum(yZF))), cex=1, type="o", pch=19, col="#8DA0CB")
lines(xHD,(cumsum(yHD)/max(cumsum(yHD))), cex=1, type="o", pch=19, col="#FC8D62")
lines(xGPCRs,(cumsum(yGPCRs)/max(cumsum(yGPCRs))), cex=1, type="o", pch=19, col="#E5C494")
lines(xICs,(cumsum(yICs)/max(cumsum(yICs))), cex=1, type="o", pch=19, col="#A6D854")
lines(xCAMs,(cumsum(yCAMs)/max(cumsum(yCAMs))), cex=1, type="o", pch=19, col="#E78AC3")
lines(xPGs,(cumsum(yPGs)/max(cumsum(yPGs))), cex=1, type="o", pch=19, col="#FFD92F")
lines(xRibos,(cumsum(yRibos)/max(cumsum(yRibos))), cex=1, type="o", pch=19, col="#B3B3B3")
legend("bottomright", legend=c("all TFs", "zinc finger","homeodomain","GPCRs","ion channels","cell adhesion molecules","proteoglycans","ribosomal genes"),
       col=c("#66C2A5", "#8DA0CB","#FC8D62","#E5C494","#A6D854","#E78AC3","#FFD92F","#B3B3B3"), pch = 19, lty=1, cex=1)
dev.off()
