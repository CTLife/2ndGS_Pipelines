
#RASDA
###RNA-Seq Data Analyzer.

You can run the scripts step by step or run one of them for a specific purpose.                  
Please run `perl RASDA-N.pl -help` to know that what can the N-th script do.                         


_____________________________________________________________________________________________________________________________
Some rules are necessary for automatic NGS data processing, such as the names of raw reads files (SRA, FASTQ or compressed FASTQ format) must be fixed for ease parsing. 

###All Rules for pipeline RASDA:
    1. Raw sequencing reads can only be stored in SRA, FASTQ or compressed FASTQ files, because only SRA format, FASTQ format and compressed FASTQ format are supported. 
       For compressed FASTQ format, there are only 7 compressed archive file formats are supported: bz2, gz, tar.gz, tar, rar, xz, zip. 
       Their suffixes must be: ".bz2",  ".gz",  ".tar.gz",  ".tar",  ".rar",  ".xz",  ".zip".
       RASDA can extract the compressed FASTQ files automatically by recognizing their suffixes.

    2. The names of raw reads files (SRA, FASTQ or compressed FASTQ format) must be: 
                                                                              groupName_fileDescription_Rep[1-9]_[1-2]_Lane[1-2].xxx 
       groupName: a positive integer number. 
                  One group only is one condition (one or more biological replicates) and has the same groupName.              
       fileDescription: a string only contains number, English letters, decimal point and hyphen.
                        That is to say,  fileDescription must match the regular expression pattern "[-.0-9A-Za-z]+" in Perl.
       Rep[1-9]:  note biological replicates, it can only be Rep1, Rep2, Rep3, Rep4, Rep5, Rep6, Rep7, Rep8 or Rep9.
       [1-2]:     note paired-end sequencing files, it can only be 1 or 2. And it is only for paired-end sequencing.  (Optional)
       Lane[1-2]: note technical replicates, it can only be Lane1 or Lane2. (Optional)
       xxx:       note file format, it can only be sra, fastq, fastq.bz2, fastq.gz, fastq.tar.gz, fastq.tar, fastq.rar, fastq.xz or fastq.zip.   

    3. For biological replicats, only "_Rep[1-9]"  is different, others are same.
       For Paired-end files,     only "_[1-2]"     is different, others are same.
       For technical replicates, only "_Lane[1-2]" is different, others are same.

    4.File name examples for 3 groups:
                 1_RNAseq-AdultMouseHeart_Rep1.sra
                 1_HRNAseq-AdultMouseHeart_Rep2.sra

                 2_RNAseq-E12.5-MouseHeart_Rep1_Lane1.fastq.bz2
                 2_RNAseq-E12.5-MouseHeart_Rep1_Lane2.fastq.bz2

                 3_RNAseq-E17.5GATA4KO-Mouse-Lab_Rep1_1.fastq.tar
                 3_RNAseq-E17.5GATA4KO-Mouse-Lab_Rep1_2.fastq.rar
                 3_RNAseq-E17.5GATA4KO-Mouse-Lab_Rep2_1.fastq.tar
                 3_RNAseq-E17.5GATA4KO-Mouse-Lab_Rep2_2.fastq.rar

    5. In your work directory, you should create a folder named "0-Other", and put folder "R_SRC" and "Adapters" into the folder "0-Other".  


###Usage:                       
     Step 1  by using RASDA1.pl, more details by "perl  RASDA1.pl  -help".                      
     Step 2  by using RASDA2.pl, more details by "perl  RASDA2.pl  -help".                       
     Step 3  by using RASDA3.pl, more details by "perl  RASDA3.pl  -help".                                  
     ......
