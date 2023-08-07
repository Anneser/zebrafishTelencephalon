# import libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
#library(tidyverse)
library(Matrix)
library(ggplot2)
#library(ggpubr)
#library(clustree)
#library(dendextend)
#library(randomForest)
library(biomaRt)
#library(circlize)
#library(rjson)
library(ComplexHeatmap)
library(readxl)
#library(igraph)
library(stringr)
library(purrr)

# set seed for reproducibility
set.seed(2023)

# import integrated and subclustered neuronal neuronal dataset
neurons.subset <- readRDS("./Data/neurons.subset.rds")

neuron_marker <- FindAllMarkers(neurons.subset, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)

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
neuron_expression <- AverageExpression(neurons.subset)
neuron_expression_25 <- neuron_expression
cell_id <- neurons.subset@active.ident
for (i in 1:length(colnames(neuron_expression_25$RNA))) {
  active_id <- colnames(neuron_expression_25$RNA)[[i]]
  active_cell_id <- which(cell_id == active_id)
  active_cells <- neurons.subset@assays$RNA[,active_cell_id]
  neuron_expression_25$RNA[,active_id] <- rowSums(active_cells > 0)/length(active_cell_id) >= 0.25 # adapt for being a percentage matrix
}

kegg <- read.csv("./Data/KEGG.csv")
gene_ids <- kegg$entrezgene_id

ensembl <- useEnsembl( "ensembl", dataset = "drerio_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", 
                                 "entrezgene_id", 
                                 "external_gene_name", 
                                 "refseq_mrna"
),
filters = "entrezgene_id",
values = gene_ids,
mart = ensembl )
idx <- match(gene_ids , genemap$entrezgene_id)

gene_list <- as.data.frame(gene_ids)
gene_list$entrez <- genemap$entrezgene_id[ idx ]
gene_list$ensembl <- genemap$ensembl_gene_id[ idx ]
gene_list$external_gene_name <- genemap$external_gene_name[ idx]

idx <- match(gene_list$entrez, kegg$entrezgene_id)

gene_list$KEGG_family <- kegg$KEGG_family[ idx ]
gene_list$KEGG_node <- kegg$KEGG_node[ idx]
gene_list$KEGG_subnode <- kegg$KEGG_subnode[ idx ]


file_list <- list.files(path = ".Data/SAMap_informativeGenes/", pattern = ".csv", all.files = FALSE,
                        full.names = TRUE, recursive = FALSE,
                        ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
celltypes_fish <- str_extract(file_list, "(?<=zf_)[^_]+")
celltypes_lizard <- str_extract(file_list, "(?<=lz_)[^.]+")

# Create a Vector with Columns
columns = c("Var1","Freq","type", "celltype") 

#Create a Empty DataFrame with 0 rows and n columns
results_df = data.frame(matrix(nrow = 0, ncol = length(columns))) 

# Assign column names
colnames(results_df) = columns

conserved_genes_df <- data.frame(matrix(nrow = 0, ncol = 11))
colnames(conserved_genes_df) <- c("X", "index", 
                                  "pval_zf", 
                                  "pval_pv",
                                  "zebrafish",
                                  "lizard",
                                  "KEGG_family",
                                  "KEGG_node",
                                  "KEGG_subnode",
                                  "zf",
                                  "lz")

for(path in file_list){
  SAMap <- read.csv(path)
  celltype <- str_extract(path, "(?<=zf_)[^_]+")
  celltype_lizard <- str_extract(path, "(?<=lz_)[^.]+")
  idx <- match(SAMap$zebrafish, gene_list$external_gene_name)
  SAMap$KEGG_family <- kegg$KEGG_family[ idx ]
  SAMap$KEGG_node <- kegg$KEGG_node[ idx]
  SAMap$KEGG_subnode <- kegg$KEGG_subnode[ idx ]
  print(celltype)
  print(celltype_lizard)
  
  
  # Use the distinct function to only keep unique rows based on the zebrafish and KEGG_family columns
  df_unique <- SAMap %>% distinct(zebrafish, KEGG_family, .keep_all = TRUE)
  print(table(df_unique$KEGG_family))
  samap_df <- as.data.frame(table(df_unique$KEGG_family)/length(df_unique$index))
  samap_df$type <- rep("conserved", length(samap_df$Var1))
  samap_df$celltype <- rep(celltype, length(samap_df$Var1))
  
  # collect names of transcription factors (dre3000) and GPCRs (dre4030)
  tf_df <- drop_na(df_unique[df_unique$KEGG_family == "dre03000",])
  gpcr_df <- drop_na(df_unique[df_unique$KEGG_family == "dre04030",])
  conserved_df <- rbind(tf_df, gpcr_df)
  conserved_df$zf <- celltype
  conserved_df$lz <- celltype_lizard
  conserved_genes_df <- rbind(conserved_genes_df, conserved_df)
  
  # collect number of genes for each comparison
  results_df <- rbind(results_df, samap_df)
}

write.csv(conserved_genes_df, ".Data/SAMap_conserved_genes.csv")



# load the necessary package
install.packages("UpSetR")
library(UpSetR)
library(reshape2)
library(dplyr)
library(purrr)
library(dplyr)
library(tidyverse)

# Subset the dataframe based on the KEGG_family 
df_dre03000 <- conserved_genes_df[conserved_genes_df$KEGG_family == "dre03000", ]
df_dre04030 <- conserved_genes_df[conserved_genes_df$KEGG_family == "dre04030", ]

zf_celltypes <- c("e10", "e11", "e13", "e14")
lz_celltypes <- c("TEGLUT17", "TEGLUT19", "TEGLUT22", "TEGLUT25")

# Select only the comparisons you are interested in
comparisons <- c("e10_TEGLUT17", "e11_TEGLUT19", "e13_TEGLUT22", "e14_TEGLUT25")

# Install necessary packages
install.packages(c("reshape2", "pheatmap"))
library(reshape2)
library(pheatmap)
library(grDevices)

# Create a custom color palette
color_palette <- colorRampPalette(c("white", "black"))(2)


# Select only the comparisons you are interested in (toggle)
comparisons <- c("e10_TEGLUT17", "e11_TEGLUT19", "e13_TEGLUT22", "e14_TEGLUT25")
#comparisons <- c("i5_TEGABA6", "i7_TEGABA7", "i9_TEGABA12")
#comparisons <- c("i33_DIGABA7", "i34_DIGABA6", "i34_DIGABA7", "i34_DIGABA8")


# Create a new column that combines cell types of species A and B
df_dre03000$comparison <- paste(df_dre03000$zf, df_dre03000$lz, sep="_")

df_dre04030$comparison <- paste(df_dre04030$zf, df_dre04030$lz, sep="_")

# Subset for comparisons of interest
df_subset <- df_dre03000[df_dre03000$comparison %in% comparisons,]

# Create a new column to represent binary presence
df_subset$value <- 1

# Spread data to wide format
df_wide <- dcast(df_subset, comparison ~ zebrafish, value.var = "value", fill = 0)

# Remove comparison column before plotting
df_heatmap <- df_wide[,-1]

# Generate heatmap
pheatmap(df_heatmap, cluster_rows = FALSE, cluster_cols = TRUE,color = color_palette)

df_subset <- df_dre04030[df_dre04030$comparison %in% comparisons,]

# Create a new column to represent binary presence
df_subset$value <- 1

# Spread data to wide format
df_wide <- dcast(df_subset, comparison ~ zebrafish, value.var = "value", fill = 0)

# Remove comparison column before plotting
df_heatmap <- df_wide[,-1]

# Generate heatmap
pheatmap(df_heatmap, cluster_rows = FALSE, cluster_cols = TRUE,color = color_palette)
