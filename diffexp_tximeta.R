## Tximeta text
#We first need to read our data into R. To do that, we will use a package called
#[tximeta](https://bioconductor.org/packages/release/bioc/html/tximeta.html).
#We use tximeta for two main reasons: first, it facilitates summarizing transcript-level counts
#from Salmon to the gene-level for differential expression analysis. Second, tximeta enhances
#incorporates metadata (e.g. transcript locations, transcript and genome source and version, 
#appropriate chromosome lengths, etc) for each transcriptome. This ensures computational reproducibility 
#by attaching critical annotation information to the data object, such that exact quantifications can be 
#reproduced from raw data (all software versions are also attached to the data object).

# [tximeta vignette](https://bioconductor.org/packages/release/bioc/vignettes/tximeta/inst/doc/tximeta.html)


### installations (not req'd for binder)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("tximeta")
#BiocManager::install("DESeq2")
#install.packages("tidyverse")
#install.packages("readr")
#suppressPackageStartupMessages(library(SummarizedExperiment)) # think this comes with tximeta already (I didn't need it)


library(tximeta)
library(DESeq2)
library(tidyverse)
library(readr)

# relevant paths should work if we're running this from within the 2020-ggg-201b-rnaseq folder
samples <- read_csv("rnaseq_samples.csv")
indexDir <- file.path("rnaseq", "quant", "sc_ensembl_index")
fastaFTP <- c("ftp://ftp.ensembl.org/pub/release-99/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz")
gtfPath <- "ftp://ftp.ensembl.org/pub/release-99/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.99.gtf.gz"

makeLinkedTxome(indexDir=indexDir,
                source="Ensembl", 
                organism="Saccharomyces cerevisiae", 
                release="99", 
                fasta=fastaFTP,
                gtf=gtfPath, 
                genome="GCA_000146045.2",
                write=FALSE)


se <- tximeta(samples)
colData(se)

# to look at the data:
assayNames(se)
rowRanges(se)
seqinfo(se)

# summarize to gene level
gse <- summarizeToGene(se)

# look at the data
rowRanges(gse)
mcols(gse)

# start deseq2 portion
dds <- DESeqDataSet(gse, ~condition)

## proceed as normal

# before making heatmap, get sample-condition from samples df
annot_col <- samples %>%
        column_to_rownames('names') %>%
        select(condition) %>%
        as.data.frame()

