# Setting a reasonable p-value threshold to be used throughout
p_cutoff <- 0.1

# this is a fold change cutoff
FC_cutoff <- log2(1.1)


# Don't run this in today's binder, it's just for reference!

#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("readr")
#install.packages("pheatmap")
# install.packages('RSQlite')
#install.packages("data.table")
#install.packages("kableExtra")
#install.packages("reshape2")
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", version = "3.11")
#BiocManager::install("tximeta")
#BiocManager::install("DESeq2")
#BiocManager::install("SummarizedExperiment")
#BiocManager::install("apeglm")

# Loading packages
library(SummarizedExperiment)
library(tximeta)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(readr)
library(pheatmap)
library(kableExtra)
library(reshape2)
library(apeglm)
# import file of sample info
samples <- read_csv("rnaseq_samples.csv")
# look at sample info
head(samples)

# the count files from last time and the index
indexDir <- file.path("rnaseq", "quant", "sc_ensembl_index")

# The same reference files from last time
gtfLocal <- "yeast_ref/Saccharomyces_cerevisiae.R64-1-1.99.gtf.gz"
fastaLocal <- "yeast_ref/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"

# assign variables to tximeta
makeLinkedTxome(indexDir = indexDir, source = "Ensembl", organism = "Saccharomyces cerevisiae",
                release = "99", genome = "GCA_000146045.2", fasta = fastaLocal, gtf = gtfLocal,
                write = FALSE)

# turn info into a summarized experiment object
se <- tximeta(samples)

# look at the new summarized experiment object
colData(se)
assayNames(se)
rowRanges(se)
seqinfo(se)
# summarize to gene
gse <- summarizeToGene(se)

# look at gse counts info
head(assays(gse)[["counts"]])
# melt the dataframe from wide to tall and save it as 'explore'
explore <- melt(assays(gse)[["counts"]])
# add column names to the new dataframe
colnames(explore) <- c("gene", "sample", "count")
# look at the first few rows of the new version
head(explore, 1000) %>%
  head()
# plot all counts by gene
ggplot(explore, aes(gene, count)) + geom_point()


# plot all counts by gene with colors
ggplot(explore, aes(gene, count, col = sample)) + geom_point()

# plot all counts by gene per sample
ggplot(explore, aes(gene, count)) + geom_point() + facet_wrap(~sample)


# calculate the principle components
yPRcomp <- prcomp(x = assays(gse)[["counts"]], center = TRUE, scale = TRUE)

# rotate the dataframe so it plots the right direction
yPCA <- as.data.frame(yPRcomp$rotation)

# Plot with ggplot, colored by sample
ggplot(yPCA, aes(PC1, PC2, col = rownames(yPCA))) + geom_point()

# add a column that describes the metadata
yPCA$condition <- c("wt", "wt", "wt", "snf2", "snf2", "snf2")

# plot again and colorize by the new column
ggplot(yPCA, aes(PC1, PC2, col = condition)) + geom_point()

# create our model
dds <- DESeqDataSet(se = gse, design = ~condition)
# run the differential expression
dds <- DESeq(dds)
# quick view all of the results
results(dds)

# save dynamically filtered results
res <- results(dds, lfcThreshold = FC_cutoff, alpha = p_cutoff, independentFiltering = TRUE)

# save results to a dataframe so we can view them more easily
resdf <- as.data.frame(results(dds))

# display the contents of the results variable
head(resdf)

# MA plot with unshrunken estimates
plotMA(res, alpha = p_cutoff)
# save shrunken results using the same contrast from the original results
shrunkres <- lfcShrink(dds = dds, coef = 2)

# MA plot with shrunken estimates
plotMA(shrunkres, alpha = p_cutoff)


# Plot your dispersion
plotDispEsts(dds)

# plot the gene with the lowest adjusted pvalue
plotCounts(dds, gene = which.min(res$padj), intgroup = "condition")

# Print session info
sessionInfo()

vsd <- vst(dds)
# calculate sample distances
sample_dists <- assay(vsd) %>%
  t() %>%
  dist() %>%
  as.matrix()

# calculate the MDS values from the distance matrix
mdsData <- data.frame(cmdscale(sample_dists))
mds <- cbind(mdsData, as.data.frame(colData(vsd)))  # combine with sample data

# And plot with ggplot2
ggplot(mds, aes(X1, X2, shape = condition)) + geom_point(size = 3) + theme_minimal()
# select genes
genes <- order(res$log2FoldChange, decreasing = TRUE)[1:20]

# use the samples dataframe to make a small dataframe of the metadata
annot_col <- samples %>%
  tibble::column_to_rownames("names") %>%
  dplyr::select(condition) %>%
  as.data.frame()

# plot the heatmap
pheatmap(assay(vsd)[genes, ])
