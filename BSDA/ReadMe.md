# BSDA                  
### BS-Seq Data Analyzer (Such as RRBS, MethylC-seq and PBAT Data.)           
You can run the scripts step by step or run one of them for a specific purpose.  
Each script is able to work independently of all others.                
Please run `perl BSDAn.pl -help` to know that what can the n-th script do.  
                                                          
Some parameters are optimized for RRBS, please change the parameters of BSDA3.pl, BSDA4.pl and BSDA5.pl for processing WGBS data. For instance, the parameters of TrimGalore, Bismark, BSseeker2 and Picard MarkDuplicates are different for analyzing RRBS, MethylC-seq and PBAT data.    
                                   
# [Download](https://github.com/CTLife/2ndGS_Pipelines/releases)                                                                                                                                                
__________________________________________________________________________________________________________________      
                                                       
Some rules are necessary for automatic NGS data processing, such as the names of raw reads files (SRA or compressed FASTQ format) must be fixed for ease of parsing. The fixed pattern of file name also gives you a clear mind for your lots of samples.

### All Rules for the Pipeline BSDA:   
    1. Raw sequencing reads can only be stored in SRA or compressed FASTQ files, because only SRA format and compressed FASTQ format are supported. 
       For compressed FASTQ format, there are only 7 compressed archive file formats are supported: bz2, gz, tar.gz, tar, rar, xz, zip. 
       Their suffixes must be: ".bz2",  ".gz",  ".tar.gz",  ".tar",  ".rar",  ".xz",  ".zip".
       BSDA can extract the compressed FASTQ files automatically by recognizing their suffixes.

    2. The names of raw reads files (SRA or compressed FASTQ format) must be: 
                                                                              groupName_fileDescription_Rep[1-9]_[1-2]_Lane[1-2].xxx 
       groupName: a positive integer number. (Required)
       fileDescription: a string only contains number, English letters, decimal point and hyphen.
                        That is to say,  fileDescription must match the regular expression pattern "[-.0-9A-Za-z]+" in Perl.
                        (Required)
       Rep[1-9]:  note biological replicates, it can only be Rep1, Rep2, Rep3, Rep4, Rep5, Rep6, Rep7, Rep8 or Rep9.
                  (Required)
       [1-2]:     note paired-end sequencing files, it can only be 1 or 2. And it is only for paired-end sequencing.  (Optional)
       Lane[1-2]: note technical replicates, it can only be Lane1 or Lane2. (Optional)
       xxx:       note file format, it can only be sra, fastq.bz2, fastq.gz, fastq.tar.gz, fastq.tar, fastq.rar, fastq.xz or fastq.zip. (Required)

    3. For biological replicats, "_Rep[1-9]"  is different.
       For Paired-end files,     only "_[1-2]" is different, others are same.
       For technical replicates, only "_Lane[1-2]" is different, others are same.

    4.File name examples:
                 1_ctrl-AdultMouseHeart_Rep1.sra
                 1_ctrl-AdultMouseHeart_Rep2.sra
                 1_KO-AdultMouseHeart_Rep1.sra
                 1_KO-AdultMouseHeart_Rep2.sra
                 2_RRBS-E12.5-MouseHeart_Rep1_Lane1.fastq.bz2
                 2_RRBS-E12.5-MouseHeart_Rep1_Lane2.fastq.bz2
                 2_WGBS-E12.5-MouseHeart_Rep2_Lane1.fastq.gz
                 2_WGBS-E12.5-MouseHeart_Rep2_Lane2.fastq.gz
                 3_male-Mouse-Lab_Rep1_1.fastq.tar
                 3_male-Mouse-Lab_Rep1_2.fastq.rar
                 3_pup-Mouse-Lab_Rep2_1.fastq.tar
                 3_pup-Mouse-Lab_Rep2_2.fastq.rar

### Usage:                                            
     Step 1  by using BSDA1.pl, more details by "perl  BSDA1.pl  -help".
     Step 2  by using BSDA2.pl, more details by "perl  BSDA2.pl  -help".
     Step 3  by using BSDA3.pl, more details by "perl  BSDA3.pl  -help".
     ......


