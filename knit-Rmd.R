#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


output_format = "html_document"
if (length(args)==1) {
  output_format = args[1]
}

rmarkdown::render("rnaseq-workflow.Rmd", output_format=output_format)
