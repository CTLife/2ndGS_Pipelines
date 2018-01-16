## example:  Rscript   BSDA9.R   2_P9D-H2KN3CCXY_Rep1   2_P9D-H2KN3CCXY_Rep1.bismark.cov    readsDis_mC   



args <- commandArgs(TRUE)
print("args: ")
print(args[1])         
print(args[2])
print(args[3])         
print("#############")

inputDir   = args[1];     ## input path
inputFile  = args[2];     ## input file
outPath    = args[3];     ## output file path

inputFILE  = paste(inputDir, inputFile , sep="/")  
outPath2   = outPath       
if( ! file.exists(outPath2) ) { dir.create(outPath2, recursive = TRUE)  }
sign1 = Sys.time()  
sign1 = gsub(" ", "_", sign1)

outFILE_1  = paste(outPath2, "1_Methylated_Frequency.txt" , sep="/")  
outFILE_2  = paste(outPath2, "2_Methylated_CumulativeFrequency.txt" , sep="/")  
outFILE_3  = paste(outPath2, "3_UnMethylated_Frequency.txt" , sep="/")  
outFILE_4  = paste(outPath2, "4_UnMethylated_CumulativeFrequency.txt" , sep="/")  
outFILE_5  = paste(outPath2, "5_Total_Frequency.txt" , sep="/")  
outFILE_6  = paste(outPath2, "6_Total_CumulativeFrequency.txt" , sep="/")  
outFILE_7  = paste(outPath2, "7_all_merge.txt" , sep="/")  




###########################################################################
MyFrequency <- function(vector1) {
    dis_vec <- rep(0, 102)
    for(i in vector1) {
        if(i>100) {dis_vec[102]=dis_vec[102]+1;}
        for(j in c(0:100)) {
            if(i==j)  {dis_vec[j+1]=dis_vec[j+1]+1;}            
        }
    }
    return(dis_vec)
}


MyCumulativeFrequency <- function(vector1) {
    dis_vec <- rep(0, 102)
    for(i in vector1) {
        for(j in c(0:101)) {
            if(i>=j) {dis_vec[j+1]=dis_vec[j+1]+1;}
        }
    }
    return(dis_vec)
}

###########################################################################





matrix_1 <- read.table(inputFILE , header=FALSE, sep="\t",  quote = "",  comment.char = "") 
dim(matrix_1)
matrix_1[1:10,]

methylated <- matrix_1[,5]
unmethylat <- matrix_1[,6] 
all_un_me  <- methylated + unmethylat
length(all_un_me)  
length(methylated) 
length(unmethylat)

meth_ratio <- methylated/all_un_me
length(meth_ratio)
meth_ratio[1:10]
length(meth_ratio[meth_ratio>1])

methylated_Fre    <- MyFrequency(methylated)
methylated_CumFre <- MyCumulativeFrequency(methylated)
unmethylat_Fre    <- MyFrequency(unmethylat)
unmethylat_CumFre <- MyCumulativeFrequency(unmethylat)
all_un_me_Fre     <- MyFrequency(all_un_me)
all_un_me_CumFre  <- MyCumulativeFrequency(all_un_me)
all_merge <- cbind(methylated_Fre, methylated_CumFre, unmethylat_Fre, unmethylat_CumFre, all_un_me_Fre, all_un_me_CumFre)
num_reads_cov <- c(0:101)

if( sum(methylated_Fre) != nrow(matrix_1) ) { print("Err-1\n") }
if( sum(unmethylat_Fre) != nrow(matrix_1) ) { print("Err-2\n") }
if( sum(all_un_me_Fre)  != nrow(matrix_1) ) { print("Err-3\n") }
 
write.table(cbind(num_reads_cov, methylated_Fre),     file = outFILE_1, append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )
write.table(cbind(num_reads_cov, methylated_CumFre),  file = outFILE_2, append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )
write.table(cbind(num_reads_cov, unmethylat_Fre),     file = outFILE_3, append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )
write.table(cbind(num_reads_cov, unmethylat_CumFre),  file = outFILE_4, append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )
write.table(cbind(num_reads_cov, all_un_me_Fre),      file = outFILE_5, append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )
write.table(cbind(num_reads_cov, all_un_me_CumFre),   file = outFILE_6, append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )
write.table(cbind(num_reads_cov, all_merge),          file = outFILE_7, append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )


###########################################################################################
library("ggplot2") 

MyTheme_1 <- function(textSize1=14, hjust1=NULL, vjust1=NULL,  angle1=NULL) {    # "hjust=1, vjust=1, angle=30" for some boxplots.
  theme(  
    line  = element_line(colour="black",  size=1.0,   linetype=1,      lineend=NULL),                                                                                        ## all line elements.          局部优先总体,下面3个也是,只对非局部设置有效.   所有线属性.
    rect  = element_rect(colour="black",  size=1.0,   linetype=1,      fill="transparent" ),                                                                                 ## all rectangluar elements.    hjust=1: 靠右对齐.   所有矩形区域属性.
    text  = element_text(family="serif",  face="plain",  colour="black",  size=textSize1, hjust=0.5, vjust=0.5,   angle=0, lineheight=1.0,  margin = NULL, debug = NULL),    ## all text elements.           "serif" for a serif font. 所有文本相关属性.
    title = element_text(family="serif",  face="plain",  colour="black",  size=textSize1, hjust=0.5, vjust=0.5,   angle=0, lineheight=1.0,  margin = NULL, debug = NULL),    ## all title elements: plot, axes, legends.    hjust:水平对齐的方向.  所有标题属性.
    ## aspect.ratio = 1,   ##高宽比
    
    axis.title    = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## label of axes (element_text; inherits from text).  horizontal: 水平的, 水平线 
    axis.title.x  = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## x axis label (element_text; inherits from axis.title)
    axis.title.y  = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=90,      lineheight=1.0,  margin = NULL, debug = NULL),       ## y axis label (element_text; inherits from axis.title)
    axis.text     = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## tick labels along axes (element_text; inherits from text). 坐标轴刻度的标签的属性.                                                         
    axis.text.x   = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=hjust1, vjust=vjust1, angle=angle1,  lineheight=1.0,  margin = NULL, debug = NULL),       ## x axis tick labels (element_text; inherits from axis.text)
    axis.text.y   = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## y axis tick labels (element_text; inherits from axis.text)
    
    axis.ticks        = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),          ## tick marks along axes (element_line; inherits from line). 坐标轴刻度线.
    axis.ticks.x      = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),          ## x axis tick marks (element_line; inherits from axis.ticks)
    axis.ticks.y      = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),          ## y axis tick marks (element_line; inherits from axis.ticks)
    axis.ticks.length = grid::unit(2.0,   "mm",   data=NULL),                                      ## length of tick marks (unit), ‘"mm"’ Millimetres.  10 mm = 1 cm.  刻度线长度
    axis.line         = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL), 	 ## lines along axes (element_line; inherits from line). 坐标轴线
    axis.line.x       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL), 	 ## line along x axis (element_line; inherits from axis.line)
    axis.line.y       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),	   ## line along y axis (element_line; inherits from axis.line)
    
    legend.background    = element_rect(colour="transparent", size=1, linetype=1, fill="transparent" ), 	      ## background of legend (element_rect; inherits from rect)
    legend.spacing       = grid::unit(1, "mm", data=NULL), 	                                                    ## extra space added around legend (unit). linetype=1指的是矩形边框的类型.
    legend.key           = element_rect(colour="transparent", size=2, linetype=1, fill="transparent" ), 	      ## background underneath legend keys. 图例符号. size=1指的是矩形边框的大小.
    legend.key.size      = grid::unit(6,   "mm", data=NULL) , 	                                                ## size of legend keys   (unit; inherits from legend.key.size)
    legend.key.height    = grid::unit(6.5, "mm", data=NULL) , 	                                                ## key background height (unit; inherits from legend.key.size)
    legend.key.width     = grid::unit(8,   "mm", data=NULL) ,                                                   ## key background width  (unit; inherits from legend.key.size)
    legend.text          = element_text(family="serif", face=NULL, colour="black", size=textSize1, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	##legend item labels. 图例文字标签.
    legend.text.align    = 0, 	                    ## alignment of legend labels (number from 0 (left) to 1 (right))
    legend.title         = element_blank(),   	    ## title of legend (element_text; inherits from title)
    legend.title.align   = 0, 	                    ## alignment of legend title (number from 0 (left) to 1 (right))
    legend.position      = "right", 	              ## the position of legends. ("left", "right", "bottom", "top", or two-element numeric vector)
    legend.direction     = "vertical",        	    ## layout of items in legends  ("horizontal" or "vertical")   图例排列方向
    legend.justification = "center",      	        ## anchor point for positioning legend inside plot ("center" or two-element numeric vector)  图例居中方式
    legend.box           = NULL, 	                  ## arrangement of multiple legends ("horizontal" or "vertical")  多图例的排列方式
    legend.box.just      = NULL, 	                  ## justification of each legend within the overall bounding box, when there are multiple legends ("top", "bottom", "left", or "right")  多图例的居中方式
    
    panel.background   = element_rect(colour="transparent", size=0.0, linetype=1, fill="transparent" ),   	## background of plotting area, drawn underneath plot (element_rect; inherits from rect)
    panel.border       = element_rect(colour="black", size=0.5, linetype=1, fill=NA ), 	                    ## border around plotting area, drawn on top of plot so that it covers tick marks and grid lines. This should be used with fill=NA (element_rect; inherits from rect)
    panel.spacing      = grid::unit(1, "mm", data=NULL) , 	                                                ## margin around facet panels (unit)  分面绘图区之间的边距
    panel.spacing.x    = grid::unit(1, "mm", data=NULL) ,
    panel.spacing.y    = grid::unit(1, "mm", data=NULL) ,
    panel.grid         = element_blank(), 	                                                                ## grid lines (element_line; inherits from line)  绘图区网格线
    panel.grid.major   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	    ## major grid lines (element_line; inherits from panel.grid)  主网格线
    panel.grid.minor   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## minor grid lines (element_line; inherits from panel.grid)  次网格线
    panel.grid.major.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	    ## vertical major grid lines (element_line; inherits from panel.grid.major)
    panel.grid.major.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## horizontal major grid lines (element_line; inherits from panel.grid.major)
    panel.grid.minor.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## vertical minor grid lines (element_line; inherits from panel.grid.minor)
    panel.grid.minor.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## horizontal minor grid lines (element_line; inherits from panel.grid.minor)
    
    plot.background	= element_rect(colour="transparent", size=NULL, linetype=NULL, fill="transparent" ),                                                ## background of the entire plot (element_rect; inherits from rect)  整个图形的背景
    plot.title      = element_text(family="serif", face=NULL, colour="black", size=textSize1, hjust=0.5, vjust=0.5,   angle=NULL, lineheight=NULL),     ## plot title (text appearance) (element_text; inherits from title)  图形标题
    plot.margin     = grid::unit(c(5, 5, 5, 5), "mm", data=NULL), 	                                                                                    ## margin around entire plot (unit with the sizes of the top, right, bottom, and left margins)
    
    strip.background = element_rect(colour=NULL,    size=NULL, linetype=NULL, fill=NULL ), 	                                                      ## background of facet labels (element_rect; inherits from rect)  分面标签背景
    strip.text       = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	    ## facet labels (element_text; inherits from text)
    strip.text.x     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	    ## facet labels along horizontal direction (element_text; inherits from strip.text)
    strip.text.y     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL)   	  ## facet labels along vertical direction (element_text; inherits from strip.text) 
  ) 
} 


MySaveGgplot2_1 <- function(ggplot2Figure1,  path1, fileName1,  height1, width1) {
  SVG1 <- paste(path1,  "/",  "SVG",  sep = "",  collapse = NULL)
  PNG1 <- paste(path1,  "/",  "PNG",  sep = "",  collapse = NULL)
  PDF1 <- paste(path1,  "/",  "PDF",  sep = "",  collapse = NULL)
  EPS1 <- paste(path1,  "/",  "EPS",  sep = "",  collapse = NULL)
  if( ! file.exists(SVG1) ) { dir.create(SVG1) }
  if( ! file.exists(PNG1) ) { dir.create(PNG1) }
  if( ! file.exists(PDF1) ) { dir.create(PDF1) }
  if( ! file.exists(EPS1) ) { dir.create(EPS1) }
  ggsave( filename = paste(SVG1,  "/",  fileName1,  ".svg",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200 )
  ggsave( filename = paste(PNG1,  "/",  fileName1,  ".png",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200 )
  ggsave( filename = paste(PDF1,  "/",  fileName1,  ".pdf",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200 )
  ggsave( filename = paste(EPS1,  "/",  fileName1,  ".eps",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200,   device=cairo_ps)         
}


## relative frequency histogram, number of sample type is 1.
MyHistogram_1A <- function(vector2,  path2,   fileName2,  title2,  xLab2, height2=4,  width2=4,   xMin2=0,  xMax2=1.2,    yMin2=0,  yMax2=0.4) {
  vector2[vector2>xMax2] <- xMax2
  vector2[vector2<xMin2] <- xMin2
  dataframeB  <- data.frame( xAxis = vector2 ) 
  
  FigureTemp1 <- ggplot(data=dataframeB, aes(x=xAxis) )  +  xlab(xLab2) + ylab("Relative Frequency") +  ggtitle(title2)  +  
    ylim(yMin2, yMax2) + geom_histogram( binwidth = 0.02, aes( y = (..count..)/sum(..count..)), fill="grey45"  )  + 
    scale_x_continuous(limits=c(xMin2-0.02, xMax2+0.02),  breaks=c(0, 0.2,  0.4, 0.6, 0.8, 1.0), labels=c("0", "0.2",  "0.4", "0.6", "0.8", "1.0") ) +  
    MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )  
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_frequency-limitY",    sep="",  collapse=NULL),  height1=height2,  width1=width2)
  
  FigureTemp2 <- ggplot(data=dataframeB, aes(x=xAxis) )  +  xlab(xLab2) + ylab("Relative Frequency") +  ggtitle(title2)  + 
    geom_histogram( binwidth = 0.02, aes( y = (..count..)/sum(..count..)), fill="grey45"  )  +  
    scale_x_continuous(limits=c(xMin2-0.02, xMax2+0.02),  breaks=c(0, 0.2,  0.4, 0.6, 0.8, 1.0), labels=c("0", "0.2",  "0.4", "0.6", "0.8", "1.0")  ) + 
    MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 ) 
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "_frequency",           sep="",  collapse=NULL),  height1=height2,  width1=width2)
}
###########################################################################################





num_reads_cov2 <- c(num_reads_cov, num_reads_cov, num_reads_cov)
my_sample_type <- c( rep("methylated", 102),  rep("unmethylated", 102),  rep("total", 102) )
list_fre <- data.frame(data.frame( xAxis = num_reads_cov2, yAxis = c(methylated_Fre, unmethylat_Fre, all_un_me_Fre)/1000000,   sampleType=my_sample_type ))
list_cum <- data.frame(data.frame( xAxis = num_reads_cov2, yAxis = c(methylated_CumFre, unmethylat_CumFre, all_un_me_CumFre)/1000000,   sampleType=my_sample_type ))




##################################################################################################
FigureTemp1 <- ggplot( data = list_fre, aes(x = xAxis, y = yAxis,   colour=sampleType ) ) +  
               scale_colour_manual( values= c("red", "blue", "yellow4" ), breaks=c("methylated", "unmethylated", "total") ) +
               scale_x_continuous(  breaks=c(0,  25, 50, 75, 100 ), labels=c("0",   "25", "50", "75", "100"   )  ) +  geom_line(size=1) +  xlim(0,101) + 
               geom_point(alpha=0.5, size=2 ) +  xlab("Coverage (number of reads)") +   ylab("Number of C sites (10^6)") +   ggtitle("Frequency distribution of reads on C sites (total>=1)") +  MyTheme_1()
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=outPath2, fileName1=paste("1a_Frequency",  sign1,   sep="_")        ,  height1=5,  width1=7)



FigureTemp1 <- ggplot( data = list_fre, aes(x = xAxis, y = yAxis,   colour=sampleType ) ) +  
               scale_colour_manual( values= c("red", "blue", "yellow4" ), breaks=c("methylated", "unmethylated", "total") ) +
               scale_x_continuous(  breaks=c(5,  25, 50, 75, 100 ), labels=c("5",   "25", "50", "75", "100"   )  ) +  geom_line(size=1) +  xlim(5,101) + ylim(0,0.4) +
               geom_point(alpha=0.5, size=2 ) +  xlab("Coverage (number of reads)") +   ylab("Number of C sites (10^6)") +   ggtitle("Frequency distribution of reads on C sites (total>=1)") +  MyTheme_1()
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=outPath2, fileName1=paste("1b_Frequency",  sign1,   sep="_")        ,  height1=5,  width1=7)




FigureTemp1 <- ggplot( data = list_fre, aes(x = xAxis, y = yAxis,   colour=sampleType ) ) +  
               scale_colour_manual( values= c("red", "blue", "yellow4" ), breaks=c("methylated", "unmethylated", "total") ) +
               scale_x_continuous(  breaks=c(0,  2, 4, 6  ), labels=c("0",   "2", "4", "6"    )  ) +  geom_line(size=1) +  xlim(0,6) + 
               geom_point(alpha=0.5, size=2 ) +  xlab("Coverage (number of reads)") +   ylab("Number of C sites (10^6)") +   ggtitle("Frequency distribution of reads on C sites (total>=1)") +  MyTheme_1()
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=outPath2, fileName1=paste("1c_Frequency",  sign1,   sep="_")        ,  height1=5,  width1=7)



FigureTemp1 <- ggplot( data = list_fre, aes(x = xAxis, y = yAxis,   colour=sampleType ) ) +  
  scale_colour_manual( values= c("red", "blue", "yellow4" ), breaks=c("methylated", "unmethylated", "total") ) +
  scale_x_continuous(  breaks=c(0,  5, 10, 15, 20 ), labels=c("0",   "5", "10", "15", "20"   )  ) +  geom_line(size=1) + xlim(0,20) + 
  geom_point(alpha=0.5, size=2) +  xlab("Coverage (number of reads)") +   ylab("Number of C sites (10^6)") +   ggtitle("Frequency distribution of reads on C sites (total>=1)") +  MyTheme_1()
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=outPath2, fileName1=paste("2a_Frequency",  sign1,   sep="_")        ,  height1=5,  width1=7)




FigureTemp1 <- ggplot( data = list_fre, aes(x = xAxis, y = yAxis,   colour=sampleType ) ) +  
  scale_colour_manual( values= c("red", "blue", "yellow4" ), breaks=c("methylated", "unmethylated", "total") ) +
  scale_x_continuous(  breaks=c(2,  5, 10, 15, 20 ), labels=c("2",   "5", "10", "15", "20"   )  ) +  geom_line(size=1) + xlim(2,20) + ylim(0,1) +
  geom_point(alpha=0.5, size=2) +  xlab("Coverage (number of reads)") +   ylab("Number of C sites (10^6)") +   ggtitle("Frequency distribution of reads on C sites (total>=1)") +  MyTheme_1()
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=outPath2, fileName1=paste("2b_Frequency",  sign1,   sep="_")        ,  height1=5,  width1=7)







#################
FigureTemp1 <- ggplot( data = list_cum, aes(x = xAxis, y = yAxis,   colour=sampleType ) ) +  
  scale_colour_manual( values= c("red", "blue", "yellow4" ), breaks=c("methylated", "unmethylated", "total") ) +
  scale_x_continuous(  breaks=c(0,  25, 50, 75, 100 ), labels=c("0",   "25", "50", "75", "100")  ) +  geom_line(size=1) + xlim(0,101) +
  geom_point(alpha=0.5, size=2 ) +  xlab("Coverage (>= number of reads)") +   ylab("Number of C sites (10^6)") +   ggtitle("Cumulative frequency distribution of reads on C sites (total>=1)") +  MyTheme_1()
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=outPath2, fileName1=paste("3_CumulativeFrequency",  sign1,   sep="_")        ,  height1=5,  width1=7)







FigureTemp1 <- ggplot( data = list_cum, aes(x = xAxis, y = yAxis,   colour=sampleType ) ) +  
  scale_colour_manual( values= c("red", "blue", "yellow4" ), breaks=c("methylated", "unmethylated", "total") ) +
  scale_x_continuous(  breaks=c(0,  5, 10, 15, 20 ), labels=c("0",   "5", "10", "15", "20"   )   ) + geom_line(size=1) + xlim(0,20) + 
  geom_point(alpha=0.5, size=2) +  xlab("Coverage (>= number of reads)") +   ylab("Number of C sites (10^6)") +   ggtitle("Cumulative frequency distribution of reads on C sites (total>=1)") +  MyTheme_1()
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=outPath2, fileName1=paste("4_CumulativeFrequency",  sign1,   sep="_")        ,  height1=5,  width1=7)








FigureTemp1 <- ggplot( data = list_cum, aes(x = xAxis, y = yAxis,   colour=sampleType ) ) +  
  scale_colour_manual( values= c("red", "blue", "yellow4" ), breaks=c("methylated", "unmethylated", "total") ) +
  scale_x_continuous(  breaks=c(0,  2, 4, 6, 8 ), labels=c("0",   "2", "4", "6", "8"   )   ) + geom_line(size=1) + xlim(0,8) + 
  geom_point(alpha=0.5, size=2) +  xlab("Coverage (>= number of reads)") +   ylab("Number of C sites (10^6)") +   ggtitle("Cumulative frequency distribution of reads on C sites (total>=1)") +  MyTheme_1()
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=outPath2, fileName1=paste("5_CumulativeFrequency",  sign1,   sep="_")        ,  height1=5,  width1=7)






FigureTemp1 <- ggplot( data = list_cum, aes(x = xAxis, y = yAxis,   colour=sampleType ) ) +  
  scale_colour_manual( values= c("red", "blue", "yellow4" ), breaks=c("methylated", "unmethylated", "total") ) +
  scale_x_continuous(  breaks=c(5,  25, 50, 75, 100 ), labels=c("5",   "25", "50", "75", "100")  ) +  geom_line(size=1) + xlim(5,101) + ylim(0,1) + 
  geom_point(alpha=0.5, size=2 ) +  xlab("Coverage (>= number of reads)") +   ylab("Number of C sites (10^6)") +   ggtitle("Cumulative frequency distribution of reads on C sites (total>=1)") +  MyTheme_1()
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=outPath2, fileName1=paste("6_CumulativeFrequency",  sign1,   sep="_")        ,  height1=5,  width1=7)




######################################################################################################




MyHistogram_1A(vector2=methylated/all_un_me,  path2=outPath2,   
               fileName2=paste("7_MethylatedRatio_distribtion",  sign1,   sep="_"),  
               title2="Distributions of Methylated ratio of Cs",  xLab2="Methylated Ratio", 
               height2=5,  width2=5,   xMin2=0,  xMax2=1.0,    yMin2=0,  yMax2=0.1)  





MyHistogram_1A(vector2=methylated/all_un_me,  path2=outPath2,   
               fileName2=paste("8_MethylatedRatio_distribtion",  sign1,   sep="_"),  
               title2="Distributions of Methylated ratio of Cs",  xLab2="Methylated Ratio", 
               height2=5,  width2=5,   xMin2=0,  xMax2=1.0,    yMin2=0,  yMax2=0.05)  

 











