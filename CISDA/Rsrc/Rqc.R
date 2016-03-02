#! /usr/bin/env Rscript

## Rqc - Quality Control Tool for High-Throughput Sequencing Data.
## https://github.com/labbcb/Rqc
## args[1]: directory path that contains input file
## args[2]: output directory path
## run this script, such as:  Rscript  Rqc.R   2-FASTQ    2-FASTQ/Results/Rqc

library(Rqc)

args <- commandArgs(TRUE)
print("args: ")
print(args[1])
print(args[2])
print("#############")

inputPath = args[1];
outPath   = args[2];
if( ! file.exists(outPath) ) { dir.create(outPath) }
print(getwd())

rqc(  path = inputPath,  pattern=".fastq$",   sample = FALSE,  outdir = outPath,
      file = "Rqc_report",   openBrowser = FALSE,    workers = multicoreWorkers()  )







