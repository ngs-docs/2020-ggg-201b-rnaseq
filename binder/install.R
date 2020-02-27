pkgs = c("BiocManager", "ggplot2", "dplyr", "readr", "pheatmap", "kableExtra" ,"knitr", "rmarkdown")
ncores = parallel::detectCores()
install.packages(pkgs, Ncpus = ncores)

BiocManager::install("tximeta")
BiocManager::install("DESeq2")
