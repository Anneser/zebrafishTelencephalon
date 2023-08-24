# Analysis code for zebrafish telencephalon sequencing data
In this repository, we provide the analysis code for scRNA-seq and bulk sequencing data from zebrafish telencephalon.

## contents
- ./analysis/
  Analysis scripts.
- ./data/
  Additional files necessary for analysis.

## data availability
  Data can be downloaded as indicated in the manuscript and from the SRA depository linked with the BioProject PRJNA964500

## read mapping
We used the transcriptome annotation provided by the [Lawson lab](https://www.umassmed.edu/lawson-lab/reagents/zebrafish-transcriptome/) (V4.3.2), 
which can be freely downloaded from the linked page. For more information, see [their publication](https://elifesciences.org/articles/55792). 
For mapping, we used the [STARsolo package](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md). We copied the CellRanger whitelist 
for the 10x V3 chemistry (3M-february-2018.txt) to a local directory and mapped:

```
STAR \
--genomeDir /link/to/genome/STARindex.2.7.9a.BSgenome.Drerio.UCSC.danRer11.noAlt.Lawsone4.3.2.sjdb100/ \
--readFilesIn /link/to/fasta_files.fastq.gz \
--soloFeatures Gene Velocyto \
--soloType CB_UMI_Simple \
--soloCBwhitelist /link/to/CellRanger_3M_february_2018.txt \
--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloStrand Forward \
--outFileNamePrefix /STARsolo_2.7.9_3065F1/ \
--runThreadN 12 \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI CR UR CB UB GX GN
```

## session info
```
R version 4.2.1 (2022-06-23 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale:
[1] LC_COLLATE=German_Germany.utf8  LC_CTYPE=German_Germany.utf8    LC_MONETARY=German_Germany.utf8 LC_NUMERIC=C                   
[5] LC_TIME=German_Germany.utf8    

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ComplexHeatmap_2.12.1       randomForest_4.7-1.1        dendextend_1.16.0           clustree_0.5.0              ggraph_2.1.0               
 [6] ggpubr_0.4.0                ggplot2_3.4.3               Matrix_1.6-1                dplyr_1.0.10                SeuratObject_4.1.3         
[11] Seurat_4.3.0.1              pheatmap_1.0.12             RColorBrewer_1.1-3          tidyr_1.3.0                 DESeq2_1.36.0              
[16] SummarizedExperiment_1.26.1 Biobase_2.56.0              MatrixGenerics_1.8.1        matrixStats_0.62.0          GenomicRanges_1.48.0       
[21] GenomeInfoDb_1.32.4         IRanges_2.30.1              S4Vectors_0.34.0            BiocGenerics_0.42.0   
```

