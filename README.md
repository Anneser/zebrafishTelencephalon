<p align="center">
<img src="./Data/ColorfulTelencephalon_smoothed.jpg" width="280" height="280">
</p>

# Analysis code for zebrafish telencephalon sequencing data
In this repository, we provide the analysis code for scRNA-seq and bulk sequencing data from zebrafish telencephalon.


## contents
- ./analysis/
  Analysis scripts.
  - 00_QC.R computes basic quality control values for the scRNA-seq dataset.
  - 01_ClusterAllCells.R performs normalization and high-level clustering of the cells, according to general standards the Seurat pipeline.
  - 02_SubclusterNeurons.R selects neuronal clusters and subclusters GABAergic and glutamatergic cell types.
  - 03_BulkSequencingAnalysis.R performs the DeSEQ2 pipeline for our bulk sequencing dataset, so that we identify DEGs for different pallial regions.
  - 04_BulkCorrelationAnalysis.R correlates the pallial regions and the neuronal cell types to establish potential regions of origion for the single cell types.
  - 05_GeneFamilyAnalysis.R calculates the expression of genes belonging to different gene families within our identified cell types.
  - 06_GeneInformationContent.R performs a Random Forest analysis with which we assess how informative different gene families are to retrieve the cell type structure of the dataset.
  - 07_NeuromodulatoryNetworks.R computes expression thresholds for GPCRs, calculates graph structures of distinct neuromodulatory networks, and performs dimensionality reduction on the neuromodulation-associated gene expression dataset.
  - 08_SAMapConservedGenes.R takes the output of the SAMap pipeline to extract conserved GPCRs and transcription factors.
  - 09_MuSiC.R uses the [MuSiC](https://xuranw.github.io/MuSiC/articles/MuSiC.html) package to deconvolve contributions of single cell types to the four pallial regions from which we have obtained bulk sequencing data. 
- ./data/
  Additional files necessary for analysis.

## data availability
  Data can be downloaded as indicated in the manuscript and from the SRA depository linked with the BioProject PRJNA964500

## single cell read mapping
Reads were mapped to the zebrafish genome as provided by UCSC (danRer11, May 2017) excluding the altenative loci. We used the 
transcriptome annotation provided by the [Lawson lab](https://www.umassmed.edu/lawson-lab/reagents/zebrafish-transcriptome/) (V4.3.2), 
which can be freely downloaded from the linked page. For more information, see [their publication](https://elifesciences.org/articles/55792). 
For mapping, we used the [STARsolo package](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) (version 2.7.9a). We copied 
the CellRanger whitelist for the 10x V3 chemistry (3M-february-2018.txt) to a local directory and mapped:

```
STAR \
  --runThreadN 16 \
  --runMode genomeGenerate \
  --genomeDir /link/to/genome/STARindex.2.7.9a.BSgenome.Drerio.UCSC.danRer11.noAlt.Lawson4.3.2.sjdb100 \
  --genomeFastaFiles /link/to/genome/BSgenome.Drerio.UCSC.danRer11.noAlt.fa \
  --genomeSAindexNbases 14 \
  --sjdbGTFfile /link/to/annotation/V4.3.2.gtf \
  --sjdbOverhang 100

STAR \
--genomeDir /link/to/genome/STARindex.2.7.9a.BSgenome.Drerio.UCSC.danRer11.noAlt.Lawsone4.3.2.sjdb100/ \
--readFilesIn /link/to/fastq_file_1.gz /link/to/fastq_file_2.gz\
--soloFeatures Gene Velocyto \
--soloType CB_UMI_Simple \
--soloCBwhitelist /link/to/CellRanger_3M_february_2018.txt \
--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloStrand Forward \
--outFileNamePrefix /STARsolo_2.7.9_fastq_file/ \
--runThreadN 12 \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI CR UR CB UB GX GN
```

## bulk RNA-seq read mapping
We used the same STAR index files and used the "GeneCounts" option ('--quantMode') to count the number of reads per gene while mapping: 

```
STAR \
--genomeDir /link/to/genome/STARindex.2.7.9a.BSgenome.Drerio.UCSC.danRer11.noAlt.Lawsone4.3.2.sjdb100/ \
--readFilesIn /link/to/fastq_file.gz
--outFileNamePrefix /STARsolo_2.7.9_fastq_file/ \
--runThreadN 12 \
--readFilesCommand zcat \
--quantMode GeneCounts \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI
```

Counts from the split fastq files were merged and non-unique gene symbols were distinguished by adding a number (i.e.: the same way, it is done automatically in the STARsolo mode)

## SAMap
We closely followed the [vignette](https://github.com/atarashansky/SAMap/blob/main/SAMap_vignette.ipynb) provided in the SAMap repository.
After obtaining blast mappings between zebrafish and bearded dragon, we run the following python code:
```
from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder, transfer_annotations,
                            sankey_plot, chord_plot, CellTypeTriangles,
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd    
from samap.utils import save_samap   

fn1    = 'lizard/Pogona_vitticeps_withAnnotation.h5ad'
fn2    = 'zebrafish/scdr_v2_withAnnotation.h5ad'

filenames = {'lz':fn1,'zf':fn2}   
sm = SAMAP(
        filenames,
        f_maps = '../maps/',
        keys = {'zf':'cluster','lz':'newcluster'},
        save_processed=False
)

sm.run(neigh_from_keys = {'lz':True,'zf':True}, ncpus=24)

# should finish with:
#<samalg.SAM object at 0x7f8ec2400050>
```

