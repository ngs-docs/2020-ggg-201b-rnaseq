pkgs = c("BiocManager", "tidyverse", "readr" ,"knitr", "rmarkdown")
ncores = parallel::detectCores()
install.packages(pkgs, Ncpus = ncores)

BiocManager::install("tximeta")
BiocManager::install("DESeq2")
