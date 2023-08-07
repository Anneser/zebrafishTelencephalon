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
library(igraph)


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

# compute expression levels of genes as 10 %, 25 %, or 50 % of cells per cluster have non-zero count
neuron_expression <- AverageExpression(neurons.subset)
cell_id <- neurons.subset@active.ident
neuron_expression_10 <- neuron_expression$RNA
neuron_expression_25 <- neuron_expression_10
neuron_expression_50 <- neuron_expression_10
for (i in 1:length(colnames(neuron_expression_10))) {
  active_id <- colnames(neuron_expression_10)[[i]]
  active_cell_id <- which(cell_id == active_id)
  active_cells <- neurons.subset@assays$RNA[,active_cell_id]
  neuron_expression_10[,active_id] <- rowSums(active_cells > 0)/length(active_cell_id) >= 0.1 # 10%
  neuron_expression_25[,active_id] <- rowSums(active_cells > 0)/length(active_cell_id) >= 0.25 # 25%
  neuron_expression_50[,active_id] <- rowSums(active_cells > 0)/length(active_cell_id) >= 0.5 # 50%
}

GPCR_df <- read_xlsx('./Data/GPCR.xlsx')
ligand_list <- unique(GPCR_df[GPCR_df$neuropeptide == 1 , ]$ligand)

##### Expression of neuropeptides
df_npp <- data.frame(matrix(0, ncol=length(ligand_list), nrow = length(celltype_order)))
colnames(df_npp) <- ligand_list
rownames(df_npp) <- celltype_order

for(ligand in ligand_list){
  gene_id <- which(rownames(neuron_expression_10) == ligand)
  celltype_id <- which(neuron_expression_10[gene_id,]>0)
  df_npp[celltype_id, ligand] = df_npp[celltype_id, ligand] + 1
  celltype_id <- which(neuron_expression_25[gene_id,]>0)
  df_npp[celltype_id, ligand] = df_npp[celltype_id, ligand] + 1
  celltype_id <- which(neuron_expression_50[gene_id,]>0)
  df_npp[celltype_id, ligand] = df_npp[celltype_id, ligand] + 1
}
zero_id <- which(colSums(df_npp)==0)
df_npp <- df_npp[,-zero_id]

##### Expression of GPCRs

# create empty dataframe
df <- data.frame(matrix(0, ncol = length(celltype_order), nrow = length(celltype_order)))

#provide column names
colnames(df) <- celltype_order
rownames(df) <- celltype_order

# compute adjacency matrix for purely telencephalic networks
for(ligand in ligand_list){
  print(ligand)
  gene_id <- which(rownames(neuron_expression_25) == ligand)
  if(sum(neuron_expression_25[gene_id,])>0){
    print(sum(neuron_expression_25[gene_id,]>0))
    celltype_id <- which(neuron_expression_10[gene_id,]>0)
    receptors <- GPCR_df[GPCR_df$ligand == ligand,]$receptor
    for(receptor in receptors){
      receptor_id <- which(rownames(neuron_expression_25) == receptor)
      celltype_id_receptor <- which(neuron_expression_25[receptor_id,]>0)
      print(sum(neuron_expression_25[receptor_id,]>0))
      df[celltype_id, celltype_id_receptor] <- df[celltype_id, celltype_id_receptor] + 1
    }
  }
}

# build the graph object
network <- graph_from_adjacency_matrix(as.matrix(df))
V(network)$type <- c(rep("excitatory", 19), rep("inhibitory", 36))
# plot it
plot(network)
plot(network, edge.arrow.size=.5, vertex.label.color="black", vertex.label.dist=0,
     vertex.color=c( "pink", "skyblue")[1+(V(network)$type=="inhibitory")] ) 

##### Gs subgroup
# create empty dataframe
df.Gs <- data.frame(matrix(0, ncol = length(celltype_order), nrow = length(celltype_order)))

#provide column names
colnames(df.Gs) <- celltype_order
rownames(df.Gs) <- celltype_order

GPCR_df_Gs <- GPCR_df[GPCR_df$GPCR_type == "Gs",]
ligand_list_Gs <- unique(GPCR_df_Gs[GPCR_df_Gs$neuropeptide == 1 , ]$ligand)

for(ligand in ligand_list_Gs){
  print(ligand)
  if(is.na(ligand) == 0){
    gene_id <- which(rownames(neuron_expression_25) == ligand)
    if(sum(neuron_expression_25[gene_id,])>0){
      celltype_id <- which(neuron_expression_25[gene_id,]>0)
      receptors <- GPCR_df[GPCR_df$ligand == ligand,]$receptor
      for(receptor in receptors){
        receptor_id <- which(rownames(neuron_expression_25) == receptor)
        celltype_id_receptor <- which(neuron_expression_25[receptor_id,]>0)
        df.Gs[celltype_id, celltype_id_receptor] <- df.Gs[celltype_id, celltype_id_receptor] + 1
      }
    }
  }
}


# build the graph object
network <- graph_from_adjacency_matrix(as.matrix(df.Gs[1:55, 1:55]), weighted = TRUE)
V(network)$type <- c(rep("excitatory", 19), rep("inhibitory", 36))

# Step 3: Compute edge widths based on the number of possible channels of information propagation
E(network)$width <- E(network)$weight

# Step 4: Plot the graph
plot(network, edge.width = E(network)$width)

# plot it
network_s <- simplify( network, remove.multiple = T, remove.loops = F, 
                 edge.attr.comb=c(weight="sum") )

plot(network_s, edge.arrow.size=.5, vertex.label.color="black", vertex.label.dist=0,edge.curved=.25,
     
     vertex.color=c( "pink", "skyblue")[1+(V(network)$type=="inhibitory")] , main = expression("G"[s]* " subfamily"))


##### Gq subgroup
# create empty dataframe
df.Gq <- data.frame(matrix(0, ncol = length(celltype_order), nrow = length(celltype_order)))

#provide column names
colnames(df.Gq) <- celltype_order
rownames(df.Gq) <- celltype_order

GPCR_df_Gq <- GPCR_df[GPCR_df$GPCR_type == "Gq",]
ligand_list_Gq <- unique(GPCR_df_Gq[GPCR_df_Gq$neuropeptide == 1 , ]$ligand)

for(ligand in ligand_list_Gq){
  print(ligand)
  if(is.na(ligand) == 0){
    gene_id <- which(rownames(neuron_expression_25) == ligand)
    if(sum(neuron_expression_25[gene_id,])>0){
      celltype_id <- which(neuron_expression_25[gene_id,]>0)
      receptors <- GPCR_df[GPCR_df$ligand == ligand,]$receptor
      for(receptor in receptors){
        receptor_id <- which(rownames(neuron_expression_25) == receptor)
        celltype_id_receptor <- which(neuron_expression_25[receptor_id,]>0)
        df.Gq[celltype_id_receptor, celltype_id] <- df.Gq[celltype_id_receptor, celltype_id] + 1
      }
    }
  }
}


# build the graph object
network <- graph_from_adjacency_matrix(as.matrix(df.Gq[1:55, 1:55]))
V(network)$type <- c(rep("excitatory", 19), rep("inhibitory", 36))
# plot it
plot(network)
network_s <- simplify( network, remove.multiple = T, remove.loops = F, 
                       
                       edge.attr.comb=c(weight="sum") )

plot(network_s, edge.arrow.size=.5, vertex.label.color="black", vertex.label.dist=0,edge.curved=.25,
     
     vertex.color=c( "pink", "skyblue")[1+(V(network)$type=="inhibitory")] , main = expression("G"[q]* " subfamily"))

##### Gi/o subgroup
# create empty dataframe
df.Gio <- data.frame(matrix(0, ncol = length(celltype_order), nrow = length(celltype_order)))

#provide column names
colnames(df.Gio) <- celltype_order
rownames(df.Gio) <- celltype_order

GPCR_df_Gio <- GPCR_df[GPCR_df$GPCR_type == "Gi/o",]
ligand_list_Gio <- unique(GPCR_df_Gio[GPCR_df_Gio$neuropeptide == 1 , ]$ligand)

for(ligand in ligand_list_Gio){
  print(ligand)
  if(is.na(ligand) == 0){
    gene_id <- which(rownames(neuron_expression_25) == ligand)
    if(sum(neuron_expression_25[gene_id,])>0){
      celltype_id <- which(neuron_expression_25[gene_id,]>0)
      receptors <- GPCR_df[GPCR_df$ligand == ligand,]$receptor
      for(receptor in receptors){
        receptor_id <- which(rownames(neuron_expression_25) == receptor)
        celltype_id_receptor <- which(neuron_expression_25[receptor_id,]>0)
        df.Gio[celltype_id, celltype_id_receptor] <- df.Gio[celltype_id, celltype_id_receptor] + 1
      }
    }
  }
}


# build the graph object
network <- graph_from_adjacency_matrix(as.matrix(df.Gio[1:55, 1:55]))

# plot it
plot(network)
V(network)$type <- c(rep("excitatory", 19), rep("inhibitory", 36))
# plot it
plot(network)
network_s <- simplify( network, remove.multiple = T, remove.loops = F)

plot(network_s, edge.arrow.size=.5, vertex.label.color="black", vertex.label.dist=0,edge.curved=.25,
     
     vertex.color=c( "pink", "skyblue")[1+(V(network)$type=="inhibitory")] , main = expression("G"[i/o]* " subfamily"))

##### check NPP expression

# Load Datasets
neurons.glut <- readRDS('./Data/neurons.glut.reclustered.rds')
neurons.gaba <- readRDS('./Data/neurons.gaba.reclustered.rds')
NPP_list <- c("npy", 
              "sst1.1",
              "sst1.2",
              "vip",
              "tac3a",
              "tac3b",
              "ccka",
              "cckb",
              "penka",
              "penkb",
              "crhb",
              "cort",
              "tac1",
              "pdyn",
              "pnoca",
              "pnocb",
              "trh",
              "adcyap1a",
              "adcyap1b",
              "grp",
              "nts",
              "nmba",
              "nmbb", 
              "galn")

GPCR_list <- c("npy1r", 
               "npy2r",
               "npy4r",
               "npy7r",
               "sstr1a",
               "sstr1b",
               "sstr2a",
               "sstr3",
               "oprd1a",
               "oprd1b",
               "oprm1",
               "oprk1",
               "oprl1",
               "vipr1b",
               "vipr2",
               "crhr1",
               "tacr1a",
               "tacr2",
               "tacr3a",
               "tacr3l",
               "adcyap1r1a",
               "adcyap1r1b",
               "grpr",
               "ntsr1",
               "nmbr",
               "galr1a",
               "galr1b",
               "galr2a",
               "galr2b")

glut.cluster.averages <- AverageExpression(neurons.glut)
gaba.cluster.averages <- AverageExpression(neurons.gaba)
neuron.averages <- AverageExpression(neurons.subset)

npp_matrix_neurons <- neuron.averages$integrated[NPP_list,]
gpcr_matrix_neurons <- neuron.averages$RNA[GPCR_list,]

npp_matrix_glut <- glut.cluster.averages$integrated[NPP_list,]
gpcr_matrix_glut <- glut.cluster.averages$RNA[GPCR_list,]
npp_matrix_gaba <- gaba.cluster.averages$integrated[NPP_list,]
gpcr_matrix_gaba <- gaba.cluster.averages$RNA[GPCR_list,]


row.logfracmax <- function(x){
  log10(x/max(x))
}

for (i in 1:length(NPP_list)){
  id <- (npp_matrix_neurons[i,]<max(npp_matrix_neurons[i,]/100000))
  npp_matrix_neurons[i,id] <- max(npp_matrix_neurons[i,])/100000
}
for (i in 1:length(GPCR_list)){
  id <- (gpcr_matrix_neurons[i,]<max(gpcr_matrix_neurons[i,]/1000))
  gpcr_matrix_neurons[i,id] <- max(gpcr_matrix_neurons[i,])/1000
}

library(circlize)

col_fun = colorRamp2(c(-5, -2.5, 0), c("white", "darkblue", "red"))
split = c(rep("excitatory", 19), rep("inhibitory", 36))
cell_types = c(paste("e", 1:19, sep = ""), paste("i", 1:36, sep = ""))

npp <- Heatmap(t(apply(npp_matrix_neurons, 1, FUN = row.logfracmax)),
               name = "NPP", rect_gp = gpar(col = "white", lwd = 1),
               column_title = "cell types", row_title = "genes",
               column_split = split,
               column_names_rot = 45,
               #row_order = NPP_list,
               #column_km = 3,
               col = col_fun,
               row_km = 2,
               column_order = cell_types)

col_fun = colorRamp2(c(-3, -1.5, 0), c("darkblue", "white", "red"))
gpcr <- Heatmap(t(apply(gpcr_matrix_neurons, 1, FUN = row.logfracmax)),
                column_title = "cell types", row_title = "genes",
                name = "GPCR", rect_gp = gpar(col = "white", lwd = 1),
                column_split = split,
                column_names_rot = 45,
                column_order = cell_types,
                col = col_fun,
                #column_km = 3,
                row_km = 4)

ht_list <- npp %v% gpcr
draw(ht_list)



##### following Ripoll-Sánchez et al., 2023: 
# https://www.biorxiv.org/content/10.1101/2022.10.30.514396v2
# Compute neuropeptidergic network properties for individual ligand-GPCR networks:
# integrative | broadcasting | local | assortative
# for this, we compute the adjacency matrices for each individual 
# peptide-receptor pair

# create empty dataframe
df <- data.frame(matrix(0, ncol = length(celltype_order), nrow = length(celltype_order)))

#provide column and row names
colnames(df) <- celltype_order
rownames(df) <- celltype_order

# compute adjacency matrix for individual P-R-pairs (purely telencephalic)
# also assemble dataframe listing number of neurons expressing peptide and 
# number of neurons expression receptors
# start by pre-allocating that dataframe:
mesoscale_df <- as.data.frame(matrix(0, ncol = 4, nrow = 1))
names(mesoscale_df) <- c("ligand", "peptide", "receptor", "family")

list_id = 0
for(family_list in c(ligand_list_Gs, ligand_list_Gq, ligand_list_Gio)){
  
  for(ligand in family_list){
    print(ligand)
    gene_id <- which(rownames(neuron_expression_25) == ligand)
    if(sum(neuron_expression_25[gene_id,])>0){
      # create empty dataframe
      df <- data.frame(matrix(0, ncol = length(celltype_order), nrow = length(celltype_order)))
      #provide column and row names
      colnames(df) <- celltype_order
      rownames(df) <- celltype_order
      print(sum(neuron_expression_25[gene_id,]>0))
      celltype_id <- which(neuron_expression_25[gene_id,]>0)
      receptors <- GPCR_df[GPCR_df$ligand == ligand,]$receptor
      receptor_max <- 0
      for(receptor in receptors){
        receptor_id <- which(rownames(neuron_expression_25) == receptor)
        celltype_id_receptor <- which(neuron_expression_25[receptor_id,]>0)
        print(sum(neuron_expression_25[receptor_id,]>0))
        if(sum(neuron_expression_25[receptor_id,]) > receptor_max){
          receptor_max <- sum(neuron_expression_25[receptor_id,])
        }
        df[celltype_id, celltype_id_receptor] <- df[celltype_id, celltype_id_receptor] + 1
      }
      new_entry <- c(ligand, sum(neuron_expression_25[gene_id,]), receptor_max, list_id)
      names(new_entry) <- c("ligand", "peptide", "receptor", "family")
      mesoscale_df <- rbind(mesoscale_df, new_entry)
      
      # build the graph object
      network <- graph_from_adjacency_matrix(as.matrix(df[1:55, 1:55]), weighted = TRUE, mode = "directed")
      
      # Step 3: Compute edge widths based on the number of possible channels of information propagation
      #E(network)$width <- E(network)$weight
      
      V(network)$type <- c(rep("excitatory", 19), rep("inhibitory", 36))
      Isolated = which(degree(network)==0)
      network = delete.vertices(network, Isolated)
      network2 = simplify(network)
      
      plot(network2, main = ligand,  vertex.label.color="black", vertex.label.dist=0, edge.arrow.size=.5,
           vertex.color=c( "pink", "skyblue")[1+(V(network)$type=="inhibitory")] ) 
    }
  }
}
# subset mesoscale_df to only include actual pairs:
mesoscale_df <- mesoscale_df[which(mesoscale_df$receptor!=0),]

for(ligand in mesoscale_df$ligand){
  if(ligand %in% ligand_list_Gio){
    mesoscale_df[mesoscale_df$ligand==ligand,]$family <- "Gio"
  }
  if(ligand %in% ligand_list_Gs){
    mesoscale_df[mesoscale_df$ligand==ligand,]$family <- "Gs"
  }
  else{
    mesoscale_df[mesoscale_df$ligand==ligand,]$family <- "Gq"
  }
}

plot(mesoscale_df$receptor, mesoscale_df$peptide, 
     ylab = "number of neurons expressing NPP", 
     xlab = "number of neurons expressing receptor", 
     xlim=c(0, 55),
     ylim=c(0, 32)
     )

points(mesoscale_df$receptor, mesoscale_df$peptide, pch = 16, col = factor(mesoscale_df$family), cex = 1.5) 
legend("topleft",
       legend = c("Gq", "Gs"),
       fill = 1:2,       # Color of the squares
       border = "black") # Color of the border of the squares

# number of ligand-receptor interactions per cell type
# create dataframe with 3 columns for monoaminergic, peptidergic and 
# peptidergic (intratelencephalic) receptors
LR_df <- as.data.frame(matrix(0, ncol = 3, nrow = 55))
names(LR_df) <- c("monoaminergic", "peptidergic", "peptidergic (intratelencephalic)")
rownames(LR_df) <- celltype_order
LR_df$monoaminergic <- colSums(metabotropic_expression) # take from 06_GeneInformationContent.R script
ligand_list <- unique(GPCR_df[GPCR_df$neuropeptide == 1 , ]$ligand)
df <- data.frame(matrix(0, ncol = 1, nrow = length(celltype_order)))
# provide column names
rownames(df) <- celltype_order

# compute adjacency matrix for purely telencephalic networks
for(ligand in ligand_list){
  print(ligand)
  gene_id <- which(rownames(neuron_expression_25) == ligand)
  if(sum(neuron_expression_25[gene_id,])>0){
    print(sum(neuron_expression_25[gene_id,]>0))
    celltype_id <- which(neuron_expression_25[gene_id,]>0)
    receptors <- GPCR_df[GPCR_df$ligand == ligand,]$receptor
    for(receptor in receptors){
      receptor_id <- which(rownames(neuron_expression_25) == receptor)
      celltype_id_receptor <- which(neuron_expression_25[receptor_id,]>0)
      print(sum(neuron_expression_25[receptor_id,]>0))
      df[celltype_id_receptor,] <- df[celltype_id_receptor,] + 1
    }
  }
}
LR_df$`peptidergic (intratelencephalic)` <- df[,1]

df <- data.frame(matrix(0, ncol = 1, nrow = length(celltype_order)))
# provide column names
rownames(df) <- celltype_order

# compute adjacency matrix for purely telencephalic networks
for(ligand in ligand_list){
  print(ligand)
  receptors <- GPCR_df[GPCR_df$ligand == ligand,]$receptor
  for(receptor in receptors){
    receptor_id <- which(rownames(neuron_expression_25) == receptor)
    celltype_id_receptor <- which(neuron_expression_25[receptor_id,]>0)
    print(sum(neuron_expression_25[receptor_id,]>0))
    df[celltype_id_receptor,] <- df[celltype_id_receptor,] + 1
    }
  }
LR_df$peptidergic<- df[,1]

heatmap = Heatmap(as.matrix(LR_df), name = "Receptor Count",
                  show_row_names = TRUE,  
                  show_column_names = TRUE,  
                  cluster_rows = FALSE,  
                  cluster_columns = FALSE, 
                  column_names_rot = 90,
                  col = colorRamp2(c(min(df), max(df)), c("white", "red"))  # change colors from white (min value) to red (max value)
)
draw(heatmap, heatmap_legend_side = "right")

library(tidyr)

# Convert data frame to long format for ggplot2
df_long <- gather(LR_df, key = "Receptor_Type", value = "Count")

# Create stacked histograms with different colors for each receptor type
ggplot(df_long, aes(x = Count, fill = Receptor_Type)) +
  geom_histogram(bins = 30, alpha = 0.5, position = "identity") +
  facet_wrap(~ Receptor_Type, ncol = 1) +
  theme_minimal() +
  labs(x = "Receptor Count", y = "Frequency")


#########################################################################
# mesoscale properties: Can cell types be sorted by network similarity? #
# #######################################################################

# With this approach we try to find underlying patterns in the gene 
# expression table of relevant neuromodulatory genes. 
# We include:
# - neuropeptide precursors
# - monoaminergic receptor
# - neuropeptidergic receptors
# - other neuromodulatory receptors


library(tsne)
library(umap)
#install.packages(c("FactoMineR", "factoextra", "cluster"))
library(FactoMineR)
library(factoextra)
library(cluster)

# get list of relevant genes
GPCR_df <- read_xlsx('./Data/GPCR.xlsx')
ligand_list <- unique(GPCR_df[GPCR_df$neuropeptide == 1 , ]$ligand)
npp_genes <- union(ligand_list, NPP_list) # list of neuropeptide precursors

# list of monoaminergic receptor
ach_id <- annotated_gene_list$KEGG_subnode == "Acetylcholine (muscarinic)" 
ach_list <- annotated_gene_list$gene_names[ach_id]
dop_id <- annotated_gene_list$KEGG_subnode == "Dopamine" 
dop_list <- annotated_gene_list$gene_names[dop_id]
ne_id <- annotated_gene_list$KEGG_subnode == "Adrenaline" 
ne_list <- annotated_gene_list$gene_names[ne_id]
se_id <- annotated_gene_list$KEGG_subnode == "Serotonin"
se_list <- annotated_gene_list$gene_names[se_id]
his_id <- annotated_gene_list$KEGG_subnode == "histamine"
his_list <- annotated_gene_list$gene_names[his_id]
monoamine_genes <- c(ach_list, dop_list, ne_list, se_list, his_list)

monoamine_df <- scale(t(neuron_expression$RNA[monoamine_genes,]), center = TRUE, scale = TRUE)
# perform tsne reduction
MA_tsne <- tsne(monoamine_df)
plot(MA_tsne[,1], MA_tsne[,2])
# perform umap embedding
MA_umap <- umap(monoamine_df)
# perform PCA
results <- prcomp(monoamine_df)
#reverse the signs of the scores
results$x <- -1*results$x
results$celltype <- c(rep(1, 19), rep(2, 36))
#calculate total variance explained by each principal component
plot(results$sdev^2 / sum(results$sdev^2))
plot(results$x[,1], results$x[,2], col = results$celltype, pch = 19)

# list of neuropeptidergic receptors
nppgpcr_genes <- union(GPCR_list, unique(GPCR_df[GPCR_df$neuropeptide == 1 , ]$receptor))

nppgpcr_df <- scale(t(neuron_expression$RNA[intersect(nppgpcr_genes, rownames(neuron_expression$RNA)),]), center = TRUE, scale = TRUE)
# perform tsne reduction
nppgpcr_tsne <- tsne(nppgpcr_df)
plot(nppgpcr_tsne[,1], nppgpcr_tsne[,2])
# perform umap embedding
nppgpcr_umap <- umap(nppgpcr_df)
plot(nppgpcr_umap$layout[,1], nppgpcr_umap$layout[,2])
# perform PCA
results <- prcomp(nppgpcr_df)
#reverse the signs of the scores
results$x <- -1*results$x
results$celltype <- c(rep(1, 19), rep(2, 36))
#calculate total variance explained by each principal component
plot(results$sdev^2 / sum(results$sdev^2))
plot(results$x[,1], results$x[,2], col = results$celltype, pch = 19)

# list of other neuromodulatory receptors
ta_list <- unique(GPCR_df[GPCR_df$ligand == "trace amine" , ]$receptor)
glut_list <- unique(GPCR_df[GPCR_df$ligand == "glutamate" , ]$receptor)
gaba_list <- unique(GPCR_df[GPCR_df$ligand == "GABA" , ]$receptor)
otherNM_genes <- c(ta_list, glut_list, gaba_list)

# join lists:
nm_list <- c(npp_genes, monoamine_genes, nppgpcr_genes, otherNM_genes)
idx <- match(nm_list, rownames(neuron_expression$RNA))

# create dataframe
df <- neuron_expression$RNA[na.omit(idx),]

# normalize for dimensionality reduction
df_scaled <- scale(t(df), center = TRUE, scale = TRUE)

# perform tsne reduction
df_tsne <- tsne(df_scaled)

# perform umap embedding
df_umap <- umap(df_scaled)

results <- prcomp(df_scaled)
#reverse the signs of the scores
results$x <- -1*results$x
results$celltype <- c(rep(1, 19), rep(2, 36))
#calculate total variance explained by each principal component
plot(results$sdev^2 / sum(results$sdev^2))

plot(results$x[,1], results$x[,2], col = results$celltype, pch = 19)


# set order of cell types:
celltype_order <- c()
for(j in 1:9){
  celltype_order <- c(celltype_order, paste("e0", j, sep = ""))
}
for(j in 10:19){
  celltype_order <- c(celltype_order, paste("e", j, sep = ""))
}
for(j in 1:9){
  celltype_order <- c(celltype_order, paste("i0", j, sep = ""))
}
for(j in 10:36){
  celltype_order <- c(celltype_order, paste("i", j, sep = ""))
}
rownames(df_scaled) <- celltype_order
# Find the cell type with the highest expression for each gene
highest_expr_celltypes <- apply(df_scaled, 2, function(col) rownames(df_scaled)[which.max(col)])

# Sort the columns (genes) based on the cell type with highest expression
sorted_expr_data <- df_scaled[, order(highest_expr_celltypes)]

# Create the heatmap
heatmap(sorted_expr_data, Colv = NA, Rowv = NA, scale = "none")
library(circlize)
col_fun = colorRamp2(c(-2, 2,  6), c("white","lightgrey", "black"))
col_fun(seq(-2, 6))

Heatmap(sorted_expr_data, cluster_rows = FALSE, 
        cluster_columns = FALSE, col = col_fun)

Heatmap(sorted_expr_data, col = col_fun)

ht1 = Heatmap(sorted_expr_data[1:19,1:63], cluster_rows = FALSE, cluster_columns = FALSE, col = col_fun)
ht2 = Heatmap(sorted_expr_data[20:55,64:224], cluster_rows = FALSE, cluster_columns = FALSE, col = col_fun)

ht1 


# Principal Component Analysis for Dimensionality Reduction
res.pca <- PCA(df_scaled, graph = FALSE)

# Get principal components
pca_data <- get_pca_ind(res.pca)
pca_data$clusters_cells <- c(rep(3, 19), rep(2, 36))

# You can visualize the result of PCA 
# Change the axes according to your needs
fviz_pca_ind(res.pca, 
             label = "none", 
             habillage = as.factor(pca_data$clusters_cells), 
             palette = "jco", 
             addEllipses = FALSE, 
             ellipse.level = 0.95)


# Create the data for the chart (number of genes per neuromodulatory family)
A <- c(31, 106, 68, 19)

# Plot the bar chart 
barplot(A, ylab = "number of genes")
