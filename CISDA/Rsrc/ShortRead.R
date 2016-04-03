#! /usr/bin/env Rscript

## ShortRead: a Bioconductor package for input, quality assessment and exploration of high-throughput sequence data.
## https://bioconductor.org/packages/release/bioc/html/ShortRead.html
## args[1]: directory path that contains input file
## args[2]: output directory path
## run this script, such as:  Rscript  ShortRead.R   2-FASTQ    2-FASTQ/Results/ShortRead

library(ShortRead)

args <- commandArgs(TRUE)
print("args: ")
print(args[1])
print(args[2])
print("#############")

inputPath = args[1];
outPath   = args[2];
if( ! file.exists(outPath) ) { dir.create(outPath) }
print(getwd())
QA2_html <- paste(outPath, "/qa2_html", sep = "")
QA_html  <- paste(outPath, "/qa_html",  sep = "")
#if( ! file.exists(QA2_html) ) { dir.create(QA2_html) }
#if( ! file.exists(QA_html)  ) { dir.create(QA_html)  }

fls <- dir(inputPath,  ".fastq$",  full=TRUE)
coll <- QACollate(
        QAFastqSource(fls),   QAReadQuality(),           QAAdapterContamination(),  QANucleotideUse(),   QAQualityUse(), 
        QASequenceUse(),      QAFrequentSequence(n=10),  QANucleotideByCycle(),     QAQualityByCycle()
)
x1 <- qa2(coll, BPPARAM=SerialParam(), verbose=TRUE)
report(x1, dest=QA2_html,  type="html")


x2 <- qa(dirPath=inputPath,   pattern=".fastq$",   type="fastq",   BPPARAM=SerialParam())
report(x2, dest=QA_html,  type="html")

   


