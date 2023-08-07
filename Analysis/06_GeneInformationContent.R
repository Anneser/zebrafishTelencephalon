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
library(factoextra)
library(cluster)

# set seed to reproduce random forest values
set.seed(2023)

# import integrated and subclustered neuronal neuronal dataset
neurons.subset <- readRDS("./Data/neurons.subset.rds")

# set order of cell types:
celltype_order <- c()
for(j in 1:19){
  celltype_order <- c(celltype_order, paste("e", j, sep = ""))
}
for(j in 1:36){
  celltype_order <- c(celltype_order, paste("i", j, sep = ""))
}
neurons.subset@active.ident <- factor(neurons.subset@active.ident,
                                      levels = celltype_order)

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


# compute expression levels of genes in three different degrees of reliability: 5, 10, or 25 % of cells per cluster have non-zero count
neuron_expression <- AverageExpression(neurons.subset)
cell_id <- neurons.subset@active.ident
for (i in 1:length(colnames(neuron_expression$RNA))) {
  active_id <- colnames(neuron_expression$RNA)[[i]]
  active_cell_id <- which(cell_id == active_id)
  active_cells <- neurons.subset@assays$RNA[,active_cell_id]
  neuron_expression$RNA[,active_id] <- rowSums(active_cells > 0)/length(active_cell_id) # adapt for being a percentage matrix
}

neuron_expression_1 <- neuron_expression
cell_id <- neurons.subset@active.ident
for (i in 1:length(colnames(neuron_expression_1$RNA))) {
  active_id <- colnames(neuron_expression_1$RNA)[[i]]
  active_cell_id <- which(cell_id == active_id)
  active_cells <- neurons.subset@assays$RNA[,active_cell_id]
  neuron_expression_1$RNA[,active_id] <- rowSums(active_cells > 0)/length(active_cell_id) >= 0.01 # adapt for being a percentage matrix
}

neuron_expression_5 <- neuron_expression
cell_id <- neurons.subset@active.ident
for (i in 1:length(colnames(neuron_expression_5$RNA))) {
  active_id <- colnames(neuron_expression_5$RNA)[[i]]
  active_cell_id <- which(cell_id == active_id)
  active_cells <- neurons.subset@assays$RNA[,active_cell_id]
  neuron_expression_5$RNA[,active_id] <- rowSums(active_cells > 0)/length(active_cell_id) >= 0.05 # adapt for being a percentage matrix
}

neuron_expression_10 <- neuron_expression
cell_id <- neurons.subset@active.ident
for (i in 1:length(colnames(neuron_expression_10$RNA))) {
  active_id <- colnames(neuron_expression_10$RNA)[[i]]
  active_cell_id <- which(cell_id == active_id)
  active_cells <- neurons.subset@assays$RNA[,active_cell_id]
  neuron_expression_10$RNA[,active_id] <- rowSums(active_cells > 0)/length(active_cell_id) >= 0.1 # adapt for being a percentage matrix
}

neuron_expression_25 <- neuron_expression
cell_id <- neurons.subset@active.ident
for (i in 1:length(colnames(neuron_expression_25$RNA))) {
  active_id <- colnames(neuron_expression_25$RNA)[[i]]
  active_cell_id <- which(cell_id == active_id)
  active_cells <- neurons.subset@assays$RNA[,active_cell_id]
  neuron_expression_25$RNA[,active_id] <- rowSums(active_cells > 0)/length(active_cell_id) >= 0.25 # adapt for being a percentage matrix
}

# toggle for different modes of calculation:
#neuron_expression <- neuron_expression_10
colors <- c('#6E377D',"white",'#EBB400')

# expression of "classical" neuromodulator GPCRs: dopamine, norepinephrin, ACh, serotonin
ach_id <- annotated_gene_list$KEGG_subnode == "Acetylcholine (muscarinic)" 
ach_list <- annotated_gene_list$gene_names[ach_id]

dop_id <- annotated_gene_list$KEGG_subnode == "Dopamine" 
dop_list <- annotated_gene_list$gene_names[dop_id]

ne_id <- annotated_gene_list$KEGG_subnode == "Adrenaline" 
ne_list <- annotated_gene_list$gene_names[ne_id]

se_id <- annotated_gene_list$KEGG_subnode == "Serotonin"
se_list <- annotated_gene_list$gene_names[se_id]

metabotropic_list <- c(ach_list, dop_list, ne_list, se_list)
metabotropic_expression <- neuron_expression$RNA[metabotropic_list,]

split = c(rep(1, 19), rep(2,36))
cha = HeatmapAnnotation(
  empty = anno_empty(border = FALSE),
  foo = anno_block(gp = gpar(fill = 2:3), labels = c("glutamatergic", "GABAergic"))
)

lha = rowAnnotation(
  empty = anno_empty(border = FALSE),
  foo = anno_block(gp = gpar(fill = 4:8), labels = c("ACh",
                                                     "dopamine",
                                                     "norepinephrine",
                                                     "serotonin")))

rsplit = c(rep(1, length(ach_list)),
           rep(2, length(dop_list)),
           rep(3, length(ne_list)),
           rep(4, length(se_list)))

metabotropic_expression <- neuron_expression_1$RNA[metabotropic_list,]

hm1 <- Heatmap(metabotropic_expression[,celltype_order], #col = col_fun,
        col = colorRamp2(c(0, 1), c("#6E377D", "#EBB400")),
        row_order = rownames(metabotropic_expression),
        rect_gp = gpar(col = "white", lwd = 0.1),
        cluster_columns = FALSE,
        #column_order = celltype_order,
        column_title = "expression of classical neuromodulators, cutoff = 0.01",
        column_split = split, top_annotation = cha,
        row_split = rsplit, left_annotation = lha)

metabotropic_expression <- neuron_expression_5$RNA[metabotropic_list,]

hm5 <- Heatmap(metabotropic_expression[,celltype_order], #col = col_fun,
        col = colorRamp2(c(0, 1), c("white", "red")),
        row_order = rownames(metabotropic_expression),
        rect_gp = gpar(col = "grey", lwd = 0.1),
        cluster_columns = FALSE,
        #column_order = celltype_order,
        column_title = "expression of classical neuromodulators, cutoff = 0.05",
        column_split = split, top_annotation = cha,
        row_split = rsplit, left_annotation = lha)

metabotropic_expression <- neuron_expression_10$RNA[metabotropic_list,]

hm10 <- Heatmap(metabotropic_expression[,celltype_order], #col = col_fun,
        col = colorRamp2(c(0, 1), c("darkblue", "red")),
        row_order = rownames(metabotropic_expression),
        rect_gp = gpar(col = "white", lwd = 0.1),
        cluster_columns = FALSE,
        #column_order = celltype_order,
        column_title = "expression of classical neuromodulators, cutoff = 0.1",
        column_split = split, top_annotation = cha,
        row_split = rsplit, left_annotation = lha)

metabotropic_expression <- neuron_expression_25$RNA[metabotropic_list,]

hm25 <- Heatmap(metabotropic_expression[,celltype_order], #col = col_fun,
        col = colorRamp2(c(0, 1), c("purple", "green")),
        row_order = rownames(metabotropic_expression),
        rect_gp = gpar(col = "white", lwd = 0.1),
        cluster_columns = FALSE,
        #column_order = celltype_order,
        column_title = "expression of classical neuromodulators, cutoff = 0.25",
        column_split = split, top_annotation = cha,
        row_split = rsplit, left_annotation = lha)

pdf("./Plots/BinaryExpressionNeuromodulators.pdf")
hm1
hm5
hm10
hm25
dev.off()

# average count of modulators per cell type
metabotropic_expression <- neuron_expression$RNA[metabotropic_list,]
exp.cutoff = 0.01
ACh_average = colSums((metabotropic_expression[,celltype_order] >= exp.cutoff)[1:length(ach_list),])
dop_average = colSums((metabotropic_expression[,celltype_order] >= exp.cutoff)[(1+length(ach_list)):(length(ach_list)+length(dop_list)),])
ne_average = colSums((metabotropic_expression[,celltype_order] >= exp.cutoff)[(1+length(ach_list)+length(dop_list)):(length(ach_list)+length(dop_list)+length(ne_list)),])
se_average = colSums((metabotropic_expression[,celltype_order] >= exp.cutoff)[(1+length(ach_list)+length(dop_list)+length(ne_list)):(length(ach_list)+length(dop_list)+length(ne_list)+ length(se_list)),])

metabotropic_average <- data.frame(expression = c(ACh_average, dop_average, ne_average, se_average),
                                   neuromodulator = c(rep("ACh", length(ACh_average)),
                                                      rep("dopamine", length(dop_average)),
                                                      rep("norepinephrine", length(ne_average)),
                                                      rep("serotonin", length(se_average))),
                                   class = c(rep(c(rep("glutamatergic", 19),
                                                   rep("GABAergic", 36)),4)),
                                   cutoff = rep(0.01, 55))

exp.cutoff = c(0.05, 0.10, 0.25)
for(i in exp.cutoff){
  ACh_average = colSums((metabotropic_expression[,celltype_order] >= i)[1:length(ach_list),])
  dop_average = colSums((metabotropic_expression[,celltype_order] >= i)[(1+length(ach_list)):(length(ach_list)+length(dop_list)),])
  ne_average = colSums((metabotropic_expression[,celltype_order] >= i)[(1+length(ach_list)+length(dop_list)):(length(ach_list)+length(dop_list)+length(ne_list)),])
  se_average = colSums((metabotropic_expression[,celltype_order] >= i)[(1+length(ach_list)+length(dop_list)+length(ne_list)):(length(ach_list)+length(dop_list)+length(ne_list)+ length(se_list)),])
  helper_df <- data.frame(expression = c(ACh_average, dop_average, ne_average, se_average),
                                     neuromodulator = c(rep("ACh", length(ACh_average)),
                                                        rep("dopamine", length(dop_average)),
                                                        rep("norepinephrine", length(ne_average)),
                                                        rep("serotonin", length(se_average))),
                                     class = c(rep(c(rep("glutamatergic", 19),
                                                     rep("GABAergic", 36)),4)),
                                     cutoff = rep(i, 55))
  metabotropic_average <- rbind(metabotropic_average, helper_df)
}

meta_plot <- ggplot(metabotropic_average, aes(x=neuromodulator, y = expression, fill = class)) + 
  geom_boxplot() +
  geom_point(position=position_jitterdodge(jitter.width = 0.1,)) +
  ylab("number of expressed receptor types") +
  xlab("neuromodulatory system") +
  ylim(c(0, 17)) +
  scale_colour_brewer("neuron class", type = "seq", palette = "Set2") +
  theme_bw() + 
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.background = element_rect(fill = "white", size = 4, colour = "white"),
    legend.justification = c(0, 1),
    # legend.position = c(0, 1),
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    axis.text=element_text(size=12), #change font size of axis text
    axis.title=element_text(size=12),
    panel.grid.major = element_line(colour = "grey70", size = 0.2),
    panel.grid.minor = element_blank()
  )

pdf("./Plots/NumberOfNmReceptors.pdf")
meta_plot + facet_wrap(. ~ cutoff, ncol = 2)
dev.off()

##### Across different gene families, which are the most informative for classifying cell types?
# Load gene lists
# transcription factors based on KEGG database
tf_id = annotated_gene_list$KEGG_family == "dre03000"
tf_list = annotated_gene_list$gene_names[tf_id] # n = 987
tf_list_filtered <- tf_list[rowSums(neuron_expression$RNA[tf_list,])>0] # n = 970

# homeodomain transcription factors
homeodomain_list = subset(annotated_gene_list[annotated_gene_list$KEGG_family == "dre03000",]$KEGG_subnode, grepl('^Homeo domain', annotated_gene_list[annotated_gene_list$KEGG_family == "dre03000",]$KEGG_subnode))
homeodomain_id = annotated_gene_list$KEGG_subnode %in% homeodomain_list
homeodomain_list = annotated_gene_list$gene_names[homeodomain_id] # n = 206
homeodomain_list <- homeodomain_list[rowSums(neuron_expression$RNA[homeodomain_list,])>0] # n = 200

# Zinc fingers
zf_id = annotated_gene_list$KEGG_node == "Zinc finger"
zf_list <- annotated_gene_list$gene_names[zf_id] # n = 330
zf_list <- zf_list[rowSums(neuron_expression$RNA[zf_list,])>0] # n = 328

# ribosomal proteins
ribosome_id = annotated_gene_list$KEGG_node == "Ribosomal proteins"
ribosome_list = annotated_gene_list$gene_names[ribosome_id] # n = 86
ribosome_list <- ribosome_list[rowSums(neuron_expression$RNA[ribosome_list,])>0] # n = 86

# Proteoglycans
pg_id <- annotated_gene_list$KEGG_family == "dre00535"
pg_list <- annotated_gene_list$gene_names[pg_id] # n = 47
pg_list <- pg_list[rowSums(neuron_expression$RNA[pg_list,])>0] # n = 47

# Ion Channels
ic_id <- annotated_gene_list$KEGG_family == "dre04040"
ic_list <- annotated_gene_list$gene_names[ic_id] # n = 354
ic_list <- ic_list[rowSums(neuron_expression$RNA[ic_list,])>0] # n = 354

# Cell Adhesion Molecules 
cam_id <- annotated_gene_list$KEGG_family == "dre04515"
cam_list <- annotated_gene_list$gene_names[cam_id] # n = 151
cam_list <- cam_list[rowSums(neuron_expression$RNA[cam_list,])>0] # n = 149

# GPCRs
gpcr_id <- annotated_gene_list$KEGG_family == "dre04030"
gpcr_list <- annotated_gene_list$gene_names[gpcr_id] # n = 367
gpcr_list <- gpcr_list[rowSums(neuron_expression$RNA[gpcr_list,])>0] # n = 361

##### PREDICTION MODEL
# TRANSCRIPTION FACTORS
# subset genes
# KEGG dataset
rf.df <- as.data.frame(t(as.data.frame(neurons.subset@assays$RNA[tf_list_filtered,])))
# scale & center
rf.df <- as.data.frame(scale(rf.df, center = TRUE, scale = TRUE))
rf.df$cell_ident <- as.data.frame(neurons.subset@active.ident)[[1]]
sc.rf <- randomForest(y = rf.df$cell_ident , x = rf.df[tf_list_filtered] , data=rf.df, importance=TRUE,
                      proximity=TRUE)

# OOB estimate of error rate: 21.45 % (500 trees), no. of variables tried at each split: 31
conf_matrix <- sc.rf$confusion[celltype_order, celltype_order]
hmtf1 <- Heatmap(conf_matrix/(apply(conf_matrix,1,sum)), 
        column_order = celltype_order, 
        row_order = celltype_order,
        col = colorRamp2(c(0, 0.5, 1), c("white", "grey", "black")))
# Plotting model
plot(sc.rf)
importance(sc.rf)
varImpPlot(sc.rf)

##### RIBOSOMAL GENES

rf.df <- as.data.frame(t(as.data.frame(neurons.subset@assays$RNA[ribosome_list,])))

# scale & center
rf.df <- as.data.frame(scale(rf.df, center = TRUE, scale = TRUE))
rf.df$cell_ident <- as.data.frame(neurons.subset@active.ident)[[1]]
sc.rf.ribo <- randomForest(y = rf.df$cell_ident , x = rf.df[ribosome_list] , data=rf.df, importance=TRUE,
                      proximity=TRUE)
print(sc.rf.ribo)
# OOB estimate of error rate: 81.8 % (500 trees), no. of variables tried at each split: 9
conf_matrix <- sc.rf.ribo$confusion[celltype_order, celltype_order]
hmribo <- Heatmap(conf_matrix/(apply(conf_matrix,1,sum)), 
        column_order = celltype_order, 
        row_order = celltype_order,
        col = colorRamp2(c(0, 0.5, 1), c("white", "grey", "black")))

# GPCRs
rf.df <- as.data.frame(t(as.data.frame(neurons.subset@assays$RNA[gpcr_list,])))
# scale & center
rf.df <- as.data.frame(scale(rf.df, center = TRUE, scale = TRUE))
rf.df$cell_ident <- as.data.frame(neurons.subset@active.ident)[[1]]
sc.rf.gpcr <- randomForest(y = rf.df$cell_ident , x = rf.df[gpcr_list] , data=rf.df, importance=TRUE,
                        proximity=TRUE)
print(sc.rf.gpcr)
# OOB estimate of error rate: 42.53 % (500 trees), no. of variables tried at each split: 9
conf_matrix <- sc.rf.gpcr$confusion[celltype_order, celltype_order]
hmgpcr <- Heatmap(conf_matrix/(apply(conf_matrix,1,sum)), 
        column_order = celltype_order, 
        row_order = celltype_order,
        col = colorRamp2(c(0, 0.5, 1), c("white", "grey", "black")))

# VISUALIZATION
pdf("./Plots/HeatmapRandomForestPrediction.pdf")
hmtf1
hmgpcr
hmribo
dev.off()

# Create random forest from all annotated genes
rf.df <- as.data.frame(t(as.data.frame(neurons.subset@assays$RNA[annotated_gene_list$gene_names,])))
# scale & center
rf.df <- as.data.frame(scale(rf.df, center = TRUE, scale = TRUE))
rf.df$cell_ident <- as.data.frame(neurons.subset@active.ident)[[1]]
na_id <- which(is.na(rf.df), arr.ind=TRUE)
rf.df <- rf.df[, -c(unique(na_id[,2]))]
sc.rf.annotatedGenes <- randomForest(y = rf.df$cell_ident , x = rf.df[colnames(rf.df[1:length(colnames(rf.df))-1])] , data=rf.df, importance=TRUE,
                           proximity=TRUE)
print(sc.rf.annotatedGenes)
# OOB estimate of error rate: 17.04 % (500 trees), no. of variables tried at each split: 45
conf_matrix <- sc.rf.annotatedGenes$confusion[celltype_order, celltype_order]
hm_annotGenes <- Heatmap(conf_matrix/(apply(conf_matrix,1,sum)), 
                  column_order = celltype_order, 
                  row_order = celltype_order,
                  col = colorRamp2(c(0, 0.5, 1), c("white", "grey", "black")))


top100genes <- names(importance(sc.rf.annotatedGenes,type = 2)[order(importance(sc.rf.annotatedGenes,type = 2), decreasing = TRUE),][1:500])
summary_table <- table(annotated_gene_list[annotated_gene_list$gene_names %in% intersect(annotated_gene_list$gene_names, top100genes),]$KEGG_family)
summary_table_all <- table(annotated_gene_list$KEGG_family)


##### Compare prediction errors of different gene families
family_list <- list(ic_list,
                    # rownames(pg_expression),
                    # rownames(ribosome_expression),
                    # rownames(tf_expression),
                    zf_list,
                    # rownames(homeodomain_expression),
                    gpcr_list)
# preallocate dataframe
err.rate.df2 <- data.frame(ngenes = rep(rep(c(50,100,150,200,250,300),each = 1),3),
                           gfamily = c(rep("ion channel", 6), # 353
                                       rep("zinc finger", 6), # 329
                                       rep("gpcr", 6)), # 361
                           error.rate = rep(NA, 18))
df_id = 1
           
for(item in family_list){
  for(i in c(50,100,150,200,250,300)){
    random_genes <- sample(item, i)
    rf.df <- as.data.frame(t(as.data.frame(neurons.subset@assays$RNA[random_genes,])))
    rf.df <- as.data.frame(scale(rf.df, center = TRUE, scale = TRUE))
    rf.df$cell_ident <- as.data.frame(neurons.subset@active.ident)[[1]]
    sc.rf <- randomForest(y = rf.df$cell_ident , x = rf.df[random_genes] , data=rf.df, importance=TRUE,
                          proximity=TRUE)
    err.rate.df2$error.rate[df_id] <- min(sc.rf$err.rate[,1])
    df_id <- df_id + 1
    print(df_id)
  }
}


ggplot(err.rate.df2, aes(ngenes, error.rate, color=gfamily)) + 
  geom_jitter(size=3, width = 2) + geom_smooth(method='loess') +
  geom_hline(yintercept = 0.2322, size = 2, linetype = 2, color = "grey") + 
  annotate("text", x = 250, y = 0.216, label = "all transcription factors") +
  geom_hline(yintercept = 0.8417, size = 2, linetype = 2, color = "grey") +
  annotate("text", x = 250, y = 0.818, label = "all ribosomal proteins") +
  
  
  xlab("number of randomly sampled genes") +
  xlim(c(0,310)) + 
  ylab("random forest error rate") + 
  ylim(c(0,1)) + 
  scale_colour_brewer("gene family", type = "seq", palette = "Set2") +
  theme_bw() + 
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.background = element_rect(fill = "white", size = 4, colour = "white"),
    legend.justification = c(0, 1),
    # legend.position = c(0, 1),
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    axis.text=element_text(size=12), #change font size of axis text
    axis.title=element_text(size=12),
    panel.grid.major = element_line(colour = "grey70", size = 0.2),
    panel.grid.minor = element_blank()
  )

###################################
# modulatory interaction networks #
###################################

GPCR_df <- read_xlsx('C:\\Users\\anneluka\\Documents\\Experiments\\scRNAseq\\analysis\\KEGG\\GPCR.xlsx')

# first determine how many receptors are actually found in the dataset (thresholds: 10 and 25 %)
neuron_expression_10 <- neuron_expression$RNA
for (i in 1:length(colnames(neuron_expression_10))) {
  active_id <- colnames(neuron_expression_10)[[i]]
  active_cell_id <- which(cell_id == active_id)
  active_cells <- scData.neurons.subset@assays$RNA[,active_cell_id]
  neuron_expression_10[,active_id] <- rowSums(active_cells > 0)/length(active_cell_id) >= 0.1
}

GPCR_df$receptor_expressed_10 <- c(rep(0, length(GPCR_df$GPCR_Family)))
GPCR_df$receptor_expressed_25 <- c(rep(0, length(GPCR_df$GPCR_Family)))
id = 1
for(entrez in GPCR_df$receptor_entrez){
  gene_id <- which(annotated_gene_list$entrez == entrez)
  gene <- annotated_gene_list$gene_names[gene_id]
  GPCR_df$receptor_expressed_10[id] <- sum(neuron_expression_10$RNA[gene , ])
  GPCR_df$receptor_expressed_25[id] <- sum(neuron_expression_25$RNA[gene , ])
  id = id + 1
}

# get neuropeptides from the GPCR_df object and mine for expression in scRNAseq dataset:
npp_list <- intersect(unique(GPCR_df$ligand), rownames(neuron_expression$RNA))
# train a model to predict cell class based on gene families to determine information content: transcription factors
npp.df <- as.data.frame(t(as.data.frame(neurons.subset@assays$RNA[npp_list,])))
npp.df <- as.data.frame(scale(rf.df, center = TRUE, scale = TRUE))
npp.df$cell_ident <- as.data.frame(neurons.subset@active.ident)[[1]]
npp.rf <- randomForest(y = npp.df$cell_ident , x = npp.df[npp_list] , data=npp.df, importance=TRUE,
                       proximity=TRUE)

# heatmap / dot plot of transcription factors determining cell type identity
tf_list_subset <- c("foxo1a",
                    "foxg1a",
                    "mef2aa",
                    "mef2cb",
                    "eomesa",
                    "emx1",
                    "emx2",
                    "emx3",
                    "nr2f2",
                    "nfia",
                    "sox4a",
                    "foxp2",
                    "pax6a",
                    "sox2",
                    "sox5",
                    "sox6",
                    "dlx5a",
                    "zic3",
                    "zbtb18",
                    "zeb2b",
                    "etv1")

DefaultAssay(neurons.subset) <- "RNA"
DotPlot(neurons.subset, features = tf_list_subset) + RotatedAxis()

gpcr_list <- c("adcyap1r1a",
               "adra2c",
               "tacr1a",
               "npy2rl",
               "adrb2a",
               "sstr3",
               "oprk1",
               "gpr78a",
               "npy8br",
               "kiss1ra")

DotPlot(neurons.subset, features = gpcr_list) + RotatedAxis()

npp_list <- c("cckb",
              "adcyap1b",
              "penkb",
              "pyya",
              "vgf",
              "tac1",
              "uts1",
              "trh",
              "apln",
              "pdyn")

DotPlot(neurons.subset, features = npp_list) + RotatedAxis()

dot_list <- c(tf_list_subset, gpcr_list, npp_list)

DotPlot(neurons.subset, features = dot_list) + RotatedAxis()
