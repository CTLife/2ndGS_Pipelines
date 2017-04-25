2ndGS_Pipelines
---------------     
Perl 5 and R source codes for next generation (the second generation) sequencing data analysis by integrating lots of open-source softwares.
---------------                                                                  

+ `CISDA`: ChIP-Seq Data Analyzer (also contains MNase-Seq, DNase-Seq, ChIP-exo and all ChIPseq-like data.)                                              
                                                                  
+ `MESDA`: MethylC-Seq Data Analyzer (also contains Bisulfite Sequencing and Single-cell BS-seq data.)     
                                                                     
+ `MeDIP`: MeDIP-seq and hMeDIP-seq data analysis.                       
                                                                   
+ `RASDA`: RNA-Seq Data Analyzer  
                                      
+ `small_RASDA` : small or micro RNA-seq data analysis.                          
                                                                                          
+ `long_RASDA` : long non-coding RNA-seq data analysis.                                     
                                                                                                                   
+ `HISDA`: Hi-C  Data Analyzer            
                                                               
                                                               
  All the scripts have been fully tested and passed my quality inspection in Ubuntu Linux.                  
                                               
[Download](https://github.com/CTLife/2ndGS_Pipelines/releases)                   
---------------------------------------------------------------------------------------------                                                                     
To use all the above 4 pipelines, the recent versions of these softwares must be available in your Linux OS.          
Common Tools:                                        
1. [Perl 5](https://www.perl.org/) , [R](https://www.r-project.org/) , [NCBI SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/)   
2. bzip2, gunzip, tar, unrar, xz and unzip.  (You can install them by using "sudo apt isntall" in Ubuntu.)      
3. [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) ,  [MultiQC](http://multiqc.info/) ,  [fastq tools](http://homes.cs.washington.edu/~dcjones/fastq-tools/) ,  [FaQCs](https://github.com/chienchi/FaQCs) ,  [prinseq](http://prinseq.sourceforge.net/) ,   [fastqp](https://github.com/mdshw5/fastqp) ,  [QC3](https://github.com/slzhao/QC3) , [NGS QC Toolkit](http://www.nipgr.res.in/ngsqctoolkit.html) , [HTQC](https://sourceforge.net/projects/htqc/files/) ,  [Rqc](http://bioconductor.org/packages/release/bioc/html/Rqc.html) , [ShortRead](https://bioconductor.org/packages/release/bioc/html/ShortRead.html) ,  [seqTools](https://www.bioconductor.org/packages/release/bioc/html/seqTools.html)            
4. [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)                       
5. [SAMtools](http://www.htslib.org/), [Subread utilities](http://subread.sourceforge.net/),  [SAMstat](http://samstat.sourceforge.net/), [Qualimap2](http://qualimap.bioinfo.cipf.es/), [PRESEQ](http://smithlabresearch.org/software/preseq/), [Picard](http://broadinstitute.github.io/picard/),  [BamQC](https://github.com/s-andrews/BamQC)
6. [BEDtools](https://github.com/arq5x/bedtools2/releases) , [BEDOPS](http://bedops.readthedocs.org/en/latest/) ,  [Deeptools](http://deeptools.github.io/)
 


####Only for CISDA: 
1. [BWA](http://bio-bwa.sourceforge.net/) , [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml) , [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) ,   [Subread](http://subread.sourceforge.net/)  
2. [MACS2](https://github.com/taoliu/MACS/) , [Peak Caller Q](http://charite.github.io/Q/) ,  [HOMER](http://homer.salk.edu/homer/) , [DANPOS](https://sites.google.com/site/danposdoc/)             

                                  
####Only for MESDA:                            
                               
####Only for RASDA:                         
                            
####Only for HISDA:                                                  
                                                                                           
                                                                                                        
---------------------------------------------------------------------------------
###License: The GNU General Public License v3.0                    
                                                                         
