# only run once for initial installation.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("biomaRt")


library(DESeq2)
library(UpSetR)
library(tidyr)
library(RColorBrewer)
library(pheatmap)
library(biomaRt)
library(ggplot2)
library(ggpubr)
library(readxl)


setwd("./Data/")

# the counts originating from the split fastq files have been merged
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

# create deseq2 dataset
dds <- DESeqDataSetFromMatrix(countData = bulk_df,
                              colData = df,
                              design = ~ region)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,] # cut low-count genes

# just go for differential expression analysis to find enriched genes in different regions
dds <- DESeq(dds)

# summarize dataset for correlation analysis
bulk_df2 = bulk_df
colnames(bulk_df2) <- df$region
bulk_SE <- sapply(levels(df$region), function(x) rowMeans(bulk_df2[names(bulk_df2) %in% x]))

##### show DEGs
res <- results(dds)
resultsNames(dds)

resOrdered <- res[order(res$pvalue),]
summary(resOrdered)
sum(res$padj < 0.1, na.rm=TRUE)
plotMA(res, ylim = c(-6, 6))

DEG_list <- rownames(resOrdered[1:sum(res$padj < 0.1, na.rm=TRUE),])

# visulization of samples in feature space
vsd <- vst(dds, blind=FALSE)
#rld <- rlog(dds, blind=FALSE) # boy, this sure takes a long time. 

sampleDists <- dist(t(assay(vsd))) 

sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(vsd$fish, vsd$region, sep="-")
rownames(sampleDistMatrix) <- vsd$region
colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

colors <- c("#ff6897ff", "#f362e2f3", "#429effff", "#b385ffff")

PCA_plot <- plotPCA(vsd, intgroup=c("region") , ntop = 500, returnData = TRUE) 

pca_panel <- ggplot(data = PCA_plot , aes(x = PC1, y = PC2, color = region)) +
  geom_point(size = 5) +
  labs(
    x = "PC 1" ,
    y = "PC 2"
  ) +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        text = element_text(size = 20)) +
  scale_color_manual(values=colors)

pca_panel

# annotate dataset
ensembl <- useEnsembl(biomart = "genes", dataset = "drerio_gene_ensembl")

annotations <- getBM(attributes = c('external_gene_name', 'entrezgene_id' , 'go_id' , 'name_1006' , "interpro" , "interpro_description"),
                     filters = 'external_gene_name',
                     values = row.names(res), 
                     mart = ensembl)
annotations <- annotations[!duplicated(annotations$external_gene_name),]

res$external_gene_name <- row.names(res)
res <- as.data.frame(res)
# Merge the deseq2 and biomart output
annot_res <- merge(res, annotations, by="external_gene_name" , all.x = TRUE)

# quantitatively assess how many GPCRs are expressed.
# First, get a list with annotations for all genes
gene_names <- res$external_gene_name
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

# now, load KEGG list with all GPCRs (search term "dre04030")
kegg <- read.csv("./KEGG.csv")
# match expressed genes with the KEGG list:
idx <- match(gene_list$entrez, kegg$entrezgene_id)
gene_list$KEGG_family <- kegg$KEGG_family[ idx ]
gene_list$KEGG_node <- kegg$KEGG_node[ idx]
gene_list$KEGG_subnode <- kegg$KEGG_subnode[ idx ]

annotated_gene_list <- gene_list[is.na(gene_list$KEGG_node) == 0,]
gpcr_id <- annotated_gene_list$KEGG_family == "dre04030"
gpcr_list <- annotated_gene_list$gene_names[gpcr_id] # 207 expressed GPCRs in total

# Here, load the manually curated list of neuromodulatory GPCRs:
GPCR_df <- read.csv("GPCR.csv", sep=";")
gpcr_nm_list <- unique(GPCR_df$receptor) # 342 neuromodulatory GPCRs
#check expression:
gpcr_id <- match(gpcr_nm_list, annotated_gene_list$gene_names)
gpcr_nm_list <- annotated_gene_list$gene_names[gpcr_id] # 342 neuromodulatory GPCRs (many NAs)

# compile complete list:
gpcr_complete_list <- union(gpcr_list, gpcr_nm_list) # 213 expressed GPCRs
gpcr_nm_list <- intersect(gpcr_list, gpcr_nm_list) # 129 expressed neuromodulatory GPCRs


# create a subset containing only differentially expressed receptors
# 1. subset to differentially expressed genes
DE_res <- annot_res[which(annot_res$padj<0.05),]
# 2. subset GPCRs
DE_GPCRs <- intersect(DE_res$external_gene_name, gpcr_complete_list) # gives 79
# now check how many GPCRs in general are differentially expressed:
DE_nmGPCRs <- intersect(DE_res$external_gene_name, gpcr_nm_list) # gives 51
# 3. get expression strength for these receptors in different brain areas
Expr_GPCRs <- dds[match(unique(DE_GPCRs) , res$external_gene_name),]
rld_GPCRs <- rlog(Expr_GPCRs, blind=FALSE)  # first option for clustering
assay_GPCRs <- assay(Expr_GPCRs)
assay_GPCRs <- t(apply(assay_GPCRs[,] , 1 , scale))

# plot GPCR: numbers in total vs differentially expressed
# step 1: prepare data as list of lists
GPCR_data_list <- list(GPCRs = gpcr_complete_list, nmGPCRs = gpcr_nm_list, DEGPCRs = DE_GPCRs, DEnmGPCRs = DE_nmGPCRs)
names(GPCR_data_list) <- c("all GPCRs", "neuromodulatory GPCRs", "differentially expressed GPCRs", "differentially expressed neuromodulatory GPCRs")
upset(fromList(GPCR_data_list), order.by = "freq", point.size = 3.5, line.size = 2, nsets = 4, main.bar.color = "black", matrix.color = "darkgrey",
      text.scale = c(1.3, 1.3, 1.3,1.3,2,1.3))

# create subset containing all differentially expressed genes
DEG_genes <- dds[DEG_list,]
rld_DEG_genes <- rlog(DEG_genes, blind=FALSE)
assay_DEG_rld <- assay(rld_DEG_genes) 
assay_DEG <- assay(DEG_genes)

# adapt column names to only reflect region
colnames(assay_DEG_rld) <- df$region
# average columns with the same name
regions <- unique(colnames(assay_DEG_rld))
assay_DEG_avg <- sapply(regions, function(x)  rowMeans(assay_DEG_rld[,colnames(assay_DEG_rld) %in% x], na.rm = TRUE))

saveRDS(assay_DEG_avg, file = "./NGS_DEG.rds")

######### CONTINUE HERE

df_gpcr <- as.data.frame(colData(dds)[,c("region","fish")])
pheatmap(assay(rld_GPCRs), cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df_gpcr)

# create dataframe for averaging over brain regions
GPCR_df <- as.data.frame(assay_GPCRs)
colnames(GPCR_df) <- df$region
anat_regions <- levels(df$region)
average_GPCR <- sapply(anat_regions, function(x) rowMeans(GPCR_df[, which(colnames(GPCR_df)==x)]))
rownames(average_GPCR) <- rownames(average_GPCR)
pheatmap(average_GPCR)

# plot distribution of expression strength 
hist(log(resOrdered$baseMean))

abline(v=mean(average_GPCR),col="red", lwd = 2)
abline(v=min(average_GPCR),col="red")
abline(v=max(average_GPCR),col="red")
