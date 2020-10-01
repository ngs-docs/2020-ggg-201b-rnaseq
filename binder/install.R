pkgs = c("BiocManager", "ggplot2", "dplyr", "readr", "pheatmap", "kableExtra" ,"knitr", "rmarkdown", "RSQLite", "data.table", "kableExtra")
ncores = parallel::detectCores()
install.packages(pkgs, Ncpus = ncores)

install.packages("BiocManager")
BiocManager::install("tximeta", version="3.10")
BiocManager::install("DESeq2")
BiocManager::install("SummarizedExperiment")

