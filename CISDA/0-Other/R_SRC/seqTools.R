#! /usr/bin/env Rscript

## Package ‘seqTools’: Analysis of nucleotide, sequence and quality content on fastq files.
## https://www.bioconductor.org/packages/release/bioc/html/seqTools.html
## args[1]: directory path that contains input file
## args[2]: output directory path
## run this script, such as:  Rscript  seqTools.R   2-FASTQ     2-FASTQ/Results/seqTools

library(seqTools)

args <- commandArgs(TRUE)
print("args: ")
print(args[1])
print(args[2])
print("#############")

inputPath = args[1];
outPath   = args[2];
if( ! file.exists(outPath) ) { dir.create(outPath) }
print(getwd())


{
outPath1 <- paste(outPath, "/5-mer", sep = "")
if( ! file.exists(outPath1) ) { dir.create(outPath1) }    

filenames <- dir(path=inputPath,   pattern="*.fastq$")
fq <- fastqq(file.path(inputPath, filenames), k=5)

file1 <- paste(outPath1, "/1-seqTools.txt", sep = "")
sink(file1)
print(filenames)
cat("\n\n\n\n\n")
print("All Files:")
print(fq)
num <- length(filenames)
for (i in seq(1:num) ) {
        cat("\n\n\n\n\n")
	print("#######################################################")
	print(filenames[i])
	print(fq[i])
}
sink() 
  

file2 <- paste(outPath1, "/2-plotKmerCount.pdf", sep = "")
pdf(file2)
plotKmerCount(fq)
for (i in seq(1:num) ) {
    plotKmerCount(fq[i])
}
dev.off()

 
file3 <- paste(outPath1, "/3-plotNucFreq.pdf", sep = "")
pdf(file3)
for (i in seq(1:num) ) {
    plotNucFreq(fq, i)
}
dev.off()


file4 <- paste(outPath1, "/4-plotNucCount.pdf", sep = "")
pdf(file4)
plotNucCount(fq)
for (i in seq(1:num) ) {
    plotNucCount(fq[i])
}
plotNucCount(fq, 1)
for (i in seq(1:num) ) {
    plotNucCount(fq[i], 1)
}
plotNucCount(fq, 2)
for (i in seq(1:num) ) {
    plotNucCount(fq[i], 2)
}
plotNucCount(fq, 3)
for (i in seq(1:num) ) {
    plotNucCount(fq[i], 3)
}
plotNucCount(fq, 4)
for (i in seq(1:num) ) {
    plotNucCount(fq[i], 4)
}
dev.off()



file5 <- paste(outPath1, "/5-plotGCcontent.pdf", sep = "")
pdf(file5)  
plotGCcontent(fq) 
for (i in seq(1:num) ) {
    plotGCcontent(fq[i])
}
dev.off()


file6 <- paste(outPath1, "/6-plotPhredQuant.pdf", sep = "")
pdf(file6)  
plotMergedPhredQuant(fq, main = "Phred quantiles for all files")
for (i in seq(1:num) ) {
    plotPhredQuant(fq, i, paste("Phred quantiles for the file", i,  sep = "")  )
}
dev.off()



file7 <- paste(outPath1, "/7-Phred-distribution.pdf", sep = "")
pdf(file7)  
plotPhredDist(fq, main = "Phred distribution for all files")
for (i in seq(1:num) ) {
    plotPhredDist(fq[i], main = paste("Phred distribution for the file", i,  sep = "")  )
}
dev.off()


mtx <- cbDistMatrix(fq)
file8 <- paste(outPath1, "/8-DistMatrix.txt", sep = "")
write.table(mtx, file = file8, append = FALSE, quote = TRUE, sep = "\t",
                 eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                 col.names = TRUE, qmethod = c("escape", "double"),
                 fileEncoding = "")



hc <- hclust(as.dist(mtx))
hcd <- as.dendrogram(hc, lty=1, lwd=2)
file9 <- paste(outPath1, "/9-hclust.pdf", sep = "")
pdf(file9)  
plot(hcd, horiz=TRUE,  las=1, edgePar=list(lwd=2, lty=1, col="blue"))
plot(hcd, horiz=FALSE, las=1, edgePar=list(lwd=2, lty=1, col="green4"))
dev.off()
file9 <- paste(outPath1, "/9-hclust.svg", sep = "")
svg(file9)  
plot(hcd, horiz=TRUE,  las=1, edgePar=list(lwd=2, lty=1, col="blue"))
dev.off()
}


