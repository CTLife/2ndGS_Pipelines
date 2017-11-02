## example:  Rscript  BSseq-Saturation-2.R   BSseq-Saturation-1/1_Bismark   2_P9D-H2KN3CCXY_Rep1    BSseq-Saturation-2/1_Bismark



args <- commandArgs(TRUE)
print("args: ")
print(args[1])         
print(args[2])
print(args[3])         
print("#############")

inputDir1  = args[1];     ## input dir
inputDir2  = args[2];     ## input dir
outPath    = args[3];     ## output file path

CHG_inputFile_1 = paste(inputDir1, "/", inputDir2, "/", "100.",   inputDir2,  "/",  "CHG_context_100.", inputDir2, ".bismark.cov",  sep="")  
CHH_inputFile_1 = paste(inputDir1, "/", inputDir2, "/", "100.",   inputDir2,  "/",  "CHH_context_100.", inputDir2, ".bismark.cov",  sep="")  
CpG_inputFile_1 = paste(inputDir1, "/", inputDir2, "/", "100.",   inputDir2,  "/",  "CpG_context_100.", inputDir2, ".bismark.cov",  sep="")  

CHG_inputFile_2 = paste(inputDir1, "/", inputDir2, "/", "80.",   inputDir2,  "/",  "CHG_context_80.", inputDir2, ".bismark.cov",  sep="")  
CHH_inputFile_2 = paste(inputDir1, "/", inputDir2, "/", "80.",   inputDir2,  "/",  "CHH_context_80.", inputDir2, ".bismark.cov",  sep="")  
CpG_inputFile_2 = paste(inputDir1, "/", inputDir2, "/", "80.",   inputDir2,  "/",  "CpG_context_80.", inputDir2, ".bismark.cov",  sep="")  

CHG_inputFile_3 = paste(inputDir1, "/", inputDir2, "/", "60.",   inputDir2,  "/",  "CHG_context_60.", inputDir2, ".bismark.cov",  sep="")  
CHH_inputFile_3 = paste(inputDir1, "/", inputDir2, "/", "60.",   inputDir2,  "/",  "CHH_context_60.", inputDir2, ".bismark.cov",  sep="")  
CpG_inputFile_3 = paste(inputDir1, "/", inputDir2, "/", "60.",   inputDir2,  "/",  "CpG_context_60.", inputDir2, ".bismark.cov",  sep="")  

CHG_inputFile_4 = paste(inputDir1, "/", inputDir2, "/", "40.",   inputDir2,  "/",  "CHG_context_40.", inputDir2, ".bismark.cov",  sep="")  
CHH_inputFile_4 = paste(inputDir1, "/", inputDir2, "/", "40.",   inputDir2,  "/",  "CHH_context_40.", inputDir2, ".bismark.cov",  sep="")  
CpG_inputFile_4 = paste(inputDir1, "/", inputDir2, "/", "40.",   inputDir2,  "/",  "CpG_context_40.", inputDir2, ".bismark.cov",  sep="")  

CHG_inputFile_5 = paste(inputDir1, "/", inputDir2, "/", "20.",   inputDir2,  "/",  "CHG_context_20.", inputDir2, ".bismark.cov",  sep="")  
CHH_inputFile_5 = paste(inputDir1, "/", inputDir2, "/", "20.",   inputDir2,  "/",  "CHH_context_20.", inputDir2, ".bismark.cov",  sep="")  
CpG_inputFile_5 = paste(inputDir1, "/", inputDir2, "/", "20.",   inputDir2,  "/",  "CpG_context_20.", inputDir2, ".bismark.cov",  sep="")  

CHG_inputFile_6 = paste(inputDir1, "/", inputDir2, "/", "1.",   inputDir2,  "/",  "CHG_context_1.", inputDir2, ".bismark.cov",  sep="")  
CHH_inputFile_6 = paste(inputDir1, "/", inputDir2, "/", "1.",   inputDir2,  "/",  "CHH_context_1.", inputDir2, ".bismark.cov",  sep="")  
CpG_inputFile_6 = paste(inputDir1, "/", inputDir2, "/", "1.",   inputDir2,  "/",  "CpG_context_1.", inputDir2, ".bismark.cov",  sep="")  


outPath2   = paste(outPath,  inputDir2,   sep="/")        
if( ! file.exists(outPath2) ) { dir.create(outPath2, recursive = TRUE)  }

CHG_matrix_1 <- read.table(CHG_inputFile_1 , header=FALSE, sep="\t",  quote = "",  comment.char = "") 
CHH_matrix_1 <- read.table(CHH_inputFile_1 , header=FALSE, sep="\t",  quote = "",  comment.char = "") 
CpG_matrix_1 <- read.table(CpG_inputFile_1 , header=FALSE, sep="\t",  quote = "",  comment.char = "") 
dim(CHG_matrix_1)
dim(CHH_matrix_1)
dim(CpG_matrix_1)

CHG_matrix_2 <- read.table(CHG_inputFile_2 , header=FALSE, sep="\t",  quote = "",  comment.char = "") 
CHH_matrix_2 <- read.table(CHH_inputFile_2 , header=FALSE, sep="\t",  quote = "",  comment.char = "") 
CpG_matrix_2 <- read.table(CpG_inputFile_2 , header=FALSE, sep="\t",  quote = "",  comment.char = "") 
dim(CHG_matrix_2)
dim(CHH_matrix_2)
dim(CpG_matrix_2)

CHG_matrix_3 <- read.table(CHG_inputFile_3 , header=FALSE, sep="\t",  quote = "",  comment.char = "") 
CHH_matrix_3 <- read.table(CHH_inputFile_3 , header=FALSE, sep="\t",  quote = "",  comment.char = "") 
CpG_matrix_3 <- read.table(CpG_inputFile_3 , header=FALSE, sep="\t",  quote = "",  comment.char = "") 
dim(CHG_matrix_3)
dim(CHH_matrix_3)
dim(CpG_matrix_3)

CHG_matrix_4 <- read.table(CHG_inputFile_4 , header=FALSE, sep="\t",  quote = "",  comment.char = "") 
CHH_matrix_4 <- read.table(CHH_inputFile_4 , header=FALSE, sep="\t",  quote = "",  comment.char = "") 
CpG_matrix_4 <- read.table(CpG_inputFile_4 , header=FALSE, sep="\t",  quote = "",  comment.char = "") 
dim(CHG_matrix_4)
dim(CHH_matrix_4)
dim(CpG_matrix_4)


CHG_matrix_5 <- read.table(CHG_inputFile_5 , header=FALSE, sep="\t",  quote = "",  comment.char = "") 
CHH_matrix_5 <- read.table(CHH_inputFile_5 , header=FALSE, sep="\t",  quote = "",  comment.char = "") 
CpG_matrix_5 <- read.table(CpG_inputFile_5 , header=FALSE, sep="\t",  quote = "",  comment.char = "") 
dim(CHG_matrix_5)
dim(CHH_matrix_5)
dim(CpG_matrix_5)


CHG_matrix_6 <- read.table(CHG_inputFile_6 , header=FALSE, sep="\t",  quote = "",  comment.char = "") 
CHH_matrix_6 <- read.table(CHH_inputFile_6 , header=FALSE, sep="\t",  quote = "",  comment.char = "") 
CpG_matrix_6 <- read.table(CpG_inputFile_6 , header=FALSE, sep="\t",  quote = "",  comment.char = "") 
dim(CHG_matrix_6)
dim(CHH_matrix_6)
dim(CpG_matrix_6)
CpG_matrix_6[1:10,]


CHG_1  <- CHG_matrix_1[,5] + CHG_matrix_1[,6] 
CHH_1  <- CHH_matrix_1[,5] + CHH_matrix_1[,6] 
CpG_1  <- CpG_matrix_1[,5] + CpG_matrix_1[,6] 
length(CHG_1)
length(CHH_1)
length(CpG_1)

CHG_1  <- CHG_matrix_1[,5] + CHG_matrix_1[,6] 
CHH_1  <- CHH_matrix_1[,5] + CHH_matrix_1[,6] 
CpG_1  <- CpG_matrix_1[,5] + CpG_matrix_1[,6] 
length(CHG_1)
length(CHH_1)
length(CpG_1)

CHG_2  <- CHG_matrix_2[,5] + CHG_matrix_2[,6] 
CHH_2  <- CHH_matrix_2[,5] + CHH_matrix_2[,6] 
CpG_2  <- CpG_matrix_2[,5] + CpG_matrix_2[,6] 
length(CHG_2)
length(CHH_2)
length(CpG_2)

CHG_3  <- CHG_matrix_3[,5] + CHG_matrix_3[,6] 
CHH_3  <- CHH_matrix_3[,5] + CHH_matrix_3[,6] 
CpG_3  <- CpG_matrix_3[,5] + CpG_matrix_3[,6] 
length(CHG_3)
length(CHH_3)
length(CpG_3)

CHG_4  <- CHG_matrix_4[,5] + CHG_matrix_4[,6] 
CHH_4  <- CHH_matrix_4[,5] + CHH_matrix_4[,6] 
CpG_4  <- CpG_matrix_4[,5] + CpG_matrix_4[,6] 
length(CHG_4)
length(CHH_4)
length(CpG_4)

CHG_4  <- CHG_matrix_4[,5] + CHG_matrix_4[,6] 
CHH_4  <- CHH_matrix_4[,5] + CHH_matrix_4[,6] 
CpG_4  <- CpG_matrix_4[,5] + CpG_matrix_4[,6] 
length(CHG_4)
length(CHH_4)
length(CpG_4)

CHG_5  <- CHG_matrix_5[,5] + CHG_matrix_5[,6] 
CHH_5  <- CHH_matrix_5[,5] + CHH_matrix_5[,6] 
CpG_5  <- CpG_matrix_5[,5] + CpG_matrix_5[,6] 
length(CHG_5)
length(CHH_5)
length(CpG_5)


CHG_6  <- CHG_matrix_6[,5] + CHG_matrix_6[,6] 
CHH_6  <- CHH_matrix_6[,5] + CHH_matrix_6[,6] 
CpG_6  <- CpG_matrix_6[,5] + CpG_matrix_6[,6] 
length(CHG_6)
length(CHH_6)
length(CpG_6)


num_reads <- c( length(CHG_1), length(CHG_1[CHG_1>=5]),  length(CHH_1), length(CHH_1[CHH_1>=5]),  length(CpG_1), length(CpG_1[CpG_1>=5]), 
                length(CHG_2), length(CHG_2[CHG_2>=5]),  length(CHH_2), length(CHH_2[CHH_2>=5]),  length(CpG_2), length(CpG_2[CpG_2>=5]),
                length(CHG_3), length(CHG_3[CHG_3>=5]),  length(CHH_3), length(CHH_3[CHH_3>=5]),  length(CpG_3), length(CpG_3[CpG_3>=5]),
                length(CHG_4), length(CHG_4[CHG_4>=5]),  length(CHH_4), length(CHH_4[CHH_4>=5]),  length(CpG_4), length(CpG_4[CpG_4>=5]),
                length(CHG_5), length(CHG_5[CHG_5>=5]),  length(CHH_5), length(CHH_5[CHH_5>=5]),  length(CpG_5), length(CpG_5[CpG_5>=5]),
                length(CHG_6), length(CHG_6[CHG_6>=5]),  length(CHH_6), length(CHH_6[CHH_6>=5]),  length(CpG_6), length(CpG_6[CpG_6>=5])  )
dotKind <- c("CHG_1reads", "CHG_5reads", "CHH_1reads", "CHH_5reads", "CpG_1reads", "CpG_5reads" )
dotKind2 <- rep(dotKind, 6)
my_xAxis <- c(100, 80, 60, 40, 20, 1)
my_xAxis2 <- rep(my_xAxis, each=6) 



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




 
###############################################################################################
myFinalFrame <- data.frame( xAxis = my_xAxis2, yAxis = num_reads,   sampleType= dotKind2 )

write.table( myFinalFrame, file = paste(outPath2, "/", inputDir2, ".numberReads.txt", sep=""), append = FALSE, 
             quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )




FigureTemp1 <- ggplot( data = myFinalFrame, aes(x = xAxis, y = yAxis/1000000,   colour=sampleType ) ) +  
               scale_colour_manual( values= c("CHG_1reads"="yellow2", "CHG_5reads"="yellow4", "CHH_1reads"="blue", 
                                              "CHH_5reads"="blue4", "CpG_1reads"="red", "CpG_5reads"="red4" )  ) +
               scale_x_continuous(  breaks=c(1,  20, 40, 60, 80, 100 ), labels=c("1",  "20", "40", "60", "80", "100")  ) +  geom_line(size=1) +   
               geom_point(alpha=0.5, size=2 ) +  xlab("Sampling proportion in all mapped reads (%)") +   ylab("Number of C sits (10^6)") +   
               ggtitle( paste( "Saturation Analysis for Sample", inputDir2 , sep=" ") ) + MyTheme_1()   
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=outPath2, fileName1= paste(inputDir2, ".all", sep="") ,  height1=5,  width1=7)





cpg_num_reads <- c( length(CpG_1), length(CpG_1[CpG_1>=5]),  length(CpG_1[CpG_1>=10]), 
                    length(CpG_2), length(CpG_2[CpG_2>=5]),  length(CpG_2[CpG_2>=10]), 
                    length(CpG_3), length(CpG_3[CpG_3>=5]),  length(CpG_3[CpG_3>=10]), 
                    length(CpG_4), length(CpG_4[CpG_4>=5]),  length(CpG_3[CpG_4>=10]), 
                    length(CpG_5), length(CpG_5[CpG_5>=5]),  length(CpG_4[CpG_5>=10]), 
                    length(CpG_6), length(CpG_6[CpG_6>=5]),  length(CpG_5[CpG_6>=10])  )
cpg_dotKind <- c("CpG_1reads", "CpG_5reads", "CpG_10reads" )
cpg_dotKind2 <- rep(cpg_dotKind, 6)
cpg_my_xAxis <- c(100, 80, 60, 40, 20, 1)
cpg_my_xAxis2 <- rep(cpg_my_xAxis, each=3) 

cpg_myFinalFrame <- data.frame( xAxis = cpg_my_xAxis2, yAxis = cpg_num_reads,   sampleType= cpg_dotKind2 )

FigureTemp2 <- ggplot( data = cpg_myFinalFrame, aes(x = xAxis, y = yAxis/1000000,   colour=sampleType ) ) +  
               scale_colour_manual( values= c( "CpG_1reads"="pink", "CpG_5reads"="red", "CpG_10reads"="red4"  ) , breaks=c("CpG_1reads", "CpG_5reads", "CpG_10reads") ) +
              scale_x_continuous(  breaks=c(1,  20, 40, 60, 80, 100 ), labels=c("1",  "20", "40", "60", "80", "100")  ) +  geom_line(size=1) +   
              geom_point(alpha=0.5, size=2 ) +  xlab("Sampling proportion in all mapped reads (%)") +   ylab("Number of CpG sits (10^6)") +   
              ggtitle( paste( "Saturation Analysis for Sample", inputDir2 , sep=" ") ) + MyTheme_1()   
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=outPath2, fileName1= paste(inputDir2, ".CpG", sep="") ,  height1=5,  width1=7)




