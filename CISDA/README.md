# CISDA                  
### ChIP-Seq Data Analyzer (also contains MNase-seq, DNase-seq, ChIP-exo and all ChIPseq-like data.)           
You can run the scripts step by step or run one of them for a specific purpose.                                           
Please run `perl CISDAn.pl -help` to know that what can the n-th script do.                    
                                   
# [Download](https://github.com/CTLife/2ndGS_Pipelines/releases)                                                                                                                                                
__________________________________________________________________________________________________________________      
                                                       
Some rules are necessary for automatic NGS data processing, such as the names of raw reads files (SRA or compressed FASTQ format) must be fixed for ease parsing. 

### All Rules for NGS pipeline CISDA:
    1. Raw sequencing reads can only be stored in SRA or compressed FASTQ files, because only SRA format and compressed FASTQ format are supported. 
       For compressed FASTQ format, there are only 7 compressed archive file formats are supported: bz2, gz, tar.gz, tar, rar, xz, zip. 
       Their suffixes must be: ".bz2",  ".gz",  ".tar.gz",  ".tar",  ".rar",  ".xz",  ".zip".
       CISDA can extract the compressed FASTQ files automatically by recognizing their suffixes.

    2. The names of raw reads files (SRA or compressed FASTQ format) must be: 
                                                                              groupName_fileDescription_Rep[1-9]_[1-2]_Lane[1-2].xxx 
       groupName: a positive integer number. 
                  One group only has one ChIP-seq input sample (one or more biological replicates) and has the same groupName.              
       fileDescription: a string only contains number, English letters, decimal point and hyphen.
                        That is to say,  fileDescription must match the regular expression pattern "[-.0-9A-Za-z]+" in Perl.
                        For ChIP-seq input file, fileDescription must be initiated with "Input-".
       Rep[1-9]:  note biological replicates, it can only be Rep1, Rep2, Rep3, Rep4, Rep5, Rep6, Rep7, Rep8 or Rep9.
       [1-2]:     note paired-end sequencing files, it can only be 1 or 2. And it is only for paired-end sequencing.  (Optional)
       Lane[1-2]: note technical replicates, it can only be Lane1 or Lane2. (Optional)
       xxx:       note file format, it can only be sra, fastq.bz2, fastq.gz, fastq.tar.gz, fastq.tar, fastq.rar, fastq.xz or fastq.zip. 

    3. For biological replicats, only "_Rep[1-9]"  is different, others are same.
       For Paired-end files,     only "_[1-2]"     is different, others are same.
       For technical replicates, only "_Lane[1-2]" is different, others are same.

    4.File name examples for 3 groups:
                 1_H3K27ac-AdultMouseHeart_Rep1.sra
                 1_H3K27ac-AdultMouseHeart_Rep2.sra
                 1_H3K27me3-AdultMouseHeart_Rep1.sra
                 1_H3K27me3-AdultMouseHeart_Rep2.sra
                 1_Input-AdultMouseHeart_Rep1.sra
                 1_Input-AdultMouseHeart_Rep2.sra

                 2_H3K27ac-E12.5-MouseHeart_Rep1_Lane1.fastq.bz2
                 2_H3K27ac-E12.5-MouseHeart_Rep1_Lane2.fastq.bz2
                 2_H3K27ac-E12.5-MouseHeart_Rep2_Lane1.fastq.gz
                 2_H3K27ac-E12.5-MouseHeart_Rep2_Lane2.fastq.gz
                 2_Input-E12.5-MouseHeart_Rep1_Lane1.fastq.tar.gz
                 2_Input-E12.5-MouseHeart_Rep1_Lane2.fastq.tar.gz

                 3_E17.5GATA4-Mouse-Lab_Rep1_1.fastq.tar
                 3_E17.5GATA4-Mouse-Lab_Rep1_2.fastq.rar
                 3_E17.5GATA4-Mouse-Lab_Rep2_1.fastq.tar
                 3_E17.5GATA4-Mouse-Lab_Rep2_2.fastq.rar
                 3_Input-Mouse-Lab_Rep1_1_Lane1.fastq.xz
                 3_Input-Mouse-Lab_Rep1_1_Lane2.fastq.zip
                 3_Input-Mouse-Lab_Rep1_2_Lane1.fastq.xz
                 3_Input-Mouse-Lab_Rep1_2_Lane2.fastq.zip

    5. In your work directory, you should create a folder named "0-Other",  and put folder "Adapters",  "R_SRC", "Shortcuts", "SRC_Log" and "TXT_Files" and  into the folder "0-Other".                     


### Usage:                                            
     Step 1  by using CISDA1.pl, more details by "perl  CISDA1.pl  -help".
     Step 2  by using CISDA2.pl, more details by "perl  CISDA2.pl  -help".
     Step 3  by using CISDA3.pl, more details by "perl  CISDA3.pl  -help".
     ......


