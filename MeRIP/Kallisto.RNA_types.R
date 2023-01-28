library(stringr)


df1 <- read.table(file="abundance.tsv", header = TRUE, sep = "\t", quote = "\"'",   dec = "." )
dim(df1)
df1[1:5,]

genes = as.vector( df1[ , 1] )
tpms = as.vector( df1[ , 5] )
length( genes )
length( tpms )

genes2 = str_extract(genes, "\\|[^|]+\\|$")
genes3 = str_extract(genes2, "[^|]+")
## table(genes3)




AllResults_g <- "Z-FinalFigures"
if( ! file.exists(AllResults_g) ) { dir.create(path=AllResults_g, recursive = TRUE) }



bool_1 <- (tpms >= 0.0) & (tpms < 0.1)
bool_2 <- (tpms >= 0.1) & (tpms < 1.0)
bool_3 <- (tpms >= 1.0) & (tpms < 5.0)
bool_4 <- (tpms >= 5.0) & (tpms < 10)
bool_5 <- (tpms >= 10)  & (tpms < 100)
bool_6 <- (tpms > 100)

length(bool_1[bool_1])
length(bool_2[bool_2])
length(bool_3[bool_3])
length(bool_4[bool_4])
length(bool_5[bool_5])
length(bool_6[bool_6])



data1 = sort( table(genes3[bool_1]), decreasing = TRUE) 
data1A = data1[1:9]
data1A[10] = sum(data1[ c(10:length(data1)) ])
names(data1A) = c( names(data1[1:9]), "others" )
data1B = as.vector(data1A)
names(data1B) = names(data1A)



data2 = sort( table(genes3[bool_2]), decreasing = TRUE) 
data2A = data2[1:7]
data2A[8] = sum(data2[ c(8:length(data2)) ])
names(data2A) = c( names(data2[1:7]), "others" )
data2B = as.vector(data2A)
names(data2B) = names(data2A)




data3 = sort( table(genes3[bool_3]), decreasing = TRUE) 
data3A = data3[1:7]
data3A[8] = sum(data3[ c(8:length(data3)) ])
names(data3A) = c( names(data3[1:7]), "others" )
data3B = as.vector(data3A)
names(data3B) = names(data3A)




data4 = sort( table(genes3[bool_4]), decreasing = TRUE) 
data4A = data4[1:7]
data4A[8] = sum(data4[ c(8:length(data4)) ])
names(data4A) = c( names(data4[1:7]), "others" )
data4B = as.vector(data4A)
names(data4B) = names(data4A)





data5 = sort( table(genes3[bool_5]), decreasing = TRUE) 
data5A = data5[1:7]
data5A[8] = sum(data5[ c(8:length(data5)) ])
names(data5A) = c( names(data5[1:7]), "others" )
data5B = as.vector(data5A)
names(data5B) = names(data5A)




data6 = sort( table(genes3[bool_6]), decreasing = TRUE) 
data6A = data6[1:7]
data6A[8] = sum(data6[ c(8:length(data6)) ])
names(data6A) = c( names(data6[1:7]), "others" )
data6B = as.vector(data6A)
names(data6B) = names(data6A)




sink(  paste(AllResults_g, "1A.RNA_types.txt", sep="/")  )
cat("\n\n##### 1: \n")
print( data1  )
cat("\n\n##### 2: \n")
print( data2  )
cat("\n\n##### 3: \n")
print( data3  )
cat("\n\n##### 4: \n")
print( data4  )
cat("\n\n##### 5: \n")
print( data5  )
cat("\n\n##### 6: \n")
print( data6  )
sink()




pdf( paste(AllResults_g, "1B.RNA_types.pdf", sep="/") )
par(mfrow=c(2,3))

p1 <- barplot(  data1B   ,  las=2, main = " 0.0 <= TPM < 0.1 " )
text(x = p1, y = data1B/1.1 ,   labels = data1B, srt=90, col="red" )

p2 <- barplot(  data2B   ,  las=2, main = " 0.1 <= TPM < 1.0 " )
text(x = p2, y = data2B/1.1 ,   labels = data2B, srt=90, col="red" )

p3 <- barplot(  data3B   ,  las=2, main = " 1.0 <= TPM < 5.0 " )
text(x = p3, y = data3B/1.1 ,   labels = data3B, srt=90, col="red" )

p4 <- barplot(  data4B   ,  las=2, main = " 5.0 <= TPM < 10 " )
text(x = p4, y = data4B/1.1 ,   labels = data4B, srt=90, col="red" )

p5 <- barplot(  data5B   ,  las=2, main = " 10 <= TPM < 100 " )
text(x = p5, y = data5B/1.1 ,   labels = data5B, srt=90, col="red" )

p6 <- barplot(  data6B   ,  las=2, main = " TPM >= 100 " )
text(x = p6, y = data6B/1.1 ,   labels = data6B, srt=90, col="red" )

dev.off()





index1 = order( tpms, decreasing=TRUE )
tpms[index1]
df2 = df1[index1, ]
df2[1:10,]


write.table(x=df2, file = paste(AllResults_g, "2.TPM.sorted.txt", sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )






top100_gene = as.vector( df2[1:100, 1] )
top100_TPM  = as.vector( df2[1:100, 5] )
top100_gene
top100_TPM

top100_gene2 = str_extract(top100_gene, "\\|[^|]+\\|$")
top100_gene3 = str_extract(top100_gene2, "[^|]+")
top100_gene3
names(top100_TPM) = top100_gene3

library(RColorBrewer)
myColors <- brewer.pal(12, "Set3")
myColors[ as.factor(top100_gene3) ]



pdf( paste(AllResults_g, "3B.RNA_top100.pdf", sep="/") , width=15, height=5,)
  barplot(  log2(top100_TPM)   ,  las=2 , cex.axis=0.5,   cex.names=0.5 , col=myColors[ as.factor(top100_gene3) ] )
  barplot(  top100_TPM   ,  las=2 , cex.axis=0.5,   cex.names=0.5 , col=myColors[ as.factor(top100_gene3) ] )
dev.off()


sink( paste(AllResults_g, "3A.RNA_top100.txt", sep="/") )
print( sort(table(top100_gene3), decreasing = TRUE) )
sink()
























#############
bool_1 <- (tpms >= 1.0) & (tpms < 5.0)
bool_2 <- (tpms >= 5.0) & (tpms < 50)
bool_3 <- (tpms >= 50.0) 

length(bool_1[bool_1])
length(bool_2[bool_2])
length(bool_3[bool_3])


data1 = sort( table(genes3[bool_1]), decreasing = TRUE) 
data1A = data1[1:9]
data1A[10] = sum(data1[ c(10:length(data1)) ])
names(data1A) = c( names(data1[1:9]), "others" )
data1B = as.vector(data1A)
names(data1B) = names(data1A)


data2 = sort( table(genes3[bool_2]), decreasing = TRUE) 
data2A = data2[1:7]
data2A[8] = sum(data2[ c(8:length(data2)) ])
names(data2A) = c( names(data2[1:7]), "others" )
data2B = as.vector(data2A)
names(data2B) = names(data2A)


data3 = sort( table(genes3[bool_3]), decreasing = TRUE) 
data3A = data3[1:7]
data3A[8] = sum(data3[ c(8:length(data3)) ])
names(data3A) = c( names(data3[1:7]), "others" )
data3B = as.vector(data3A)
names(data3B) = names(data3A)


sink( paste(AllResults_g, "4A.RNA_types.txt", sep="/") ) 
cat("\n\n##### 1: \n")
print( data1  )
cat("\n\n##### 2: \n")
print( data2  )
cat("\n\n##### 3: \n")
print( data3  )
cat("\n\n##### 4: \n")
sink()


pdf(  paste(AllResults_g, "4B.RNA_types.pdf", sep="/") ) 
par(mfrow=c(2,3))

p1 <- barplot(  data1B   ,  las=2 , main = " 1 <= TPM < 5 "  )
text(x = p1, y = data1B/1.1 ,   labels = data1B, srt=90, col="red" )

p2 <- barplot(  data2B   ,  las=2 , main = " 5 <= TPM < 50 "  )
text(x = p2, y = data2B/1.1 ,   labels = data2B, srt=90, col="red" )

p3 <- barplot(  data3B   ,  las=2 , main = " TPM >= 50  "  )
text(x = p3, y = data3B/1.1 ,   labels = data3B, srt=90, col="red" )

dev.off()











