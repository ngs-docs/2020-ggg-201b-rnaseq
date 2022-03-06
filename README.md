# RNA-Seq in the Cloud

June 21 & 23, 2021

CFDE Workshop

Open access
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/nih-cfde/rnaseq-in-the-cloud/stable?labpath=rstudio)

Requires github login
[![Binder](https://aws-uswest2-binder.pangeo.io/badge_logo.svg)](https://aws-uswest2-binder.pangeo.io/v2/gh/nih-cfde/rnaseq-in-the-cloud/stable?urlpath=rstudio)

Adapted for CFDE workshop from materials developed by Taylor Reiter and N. Tessa Pierce for 2020-ggg-201b course at UC Davis.


## To run snakemake in binder --

1) start the binder

2) at the command line, run

```
snakemake -j 4 --use-conda
```

3) open `rnaseq-workflow.Rmd`

4) knit to HTML
