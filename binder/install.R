pkgs = c("BiocManager", "ggplot2", "dplyr", "readr", "pheatmap", "kableExtra" ,"knitr", "rmarkdown")
ncores = parallel::detectCores()
install.packages(pkgs, Ncpus = ncores)

if (!requireNamespace("BiocManager", quietly = TRUE))
 install.packages("BiocManager")

BiocManager::install("tximeta", version="3.10")
BiocManager::install("DESeq2")
