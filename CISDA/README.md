###You can run the scripts step by step or run one of them for a specific purpose. Please run "perl CISDA-N.pl -h" to know that what can the N-th script do.
__________________________________________________________________________________________________________________      
                                                       

###In order to realize automatic NGS data processing, the format of file names used by the pipeline must be in fixed format for ease parsing.                 
                                  
                                       
###Rules for NGS pipeline CISDA:                                         

    0. Only SRA format and compressed FASTQ files can be used as input files. 
       And the format of file name must be fixed, such as:
                    SRA file:  [01-99]_Target_Treat_Space_Time_Other_Rep[1-9].sra

       Compressed FASTQ file:  [01-99]_Target_Treat_Space_Time_Other_Rep[1-9]_[1-2]_Lane[1-2].fastq.bz2 
                               [01-99]_Target_Treat_Space_Time_Other_Rep[1-9]_[1-2]_Lane[1-2].fastq.gz
                               [01-99]_Target_Treat_Space_Time_Other_Rep[1-9]_[1-2]_Lane[1-2].fastq.tar.gz

      Only seven compressed formats are supported, their suffixes:  
                                      ".bz2",  ".gz",  ".tar.gz",  ".tar",  ".rar",  ".xz",  ".zip"
      RASDA can uncompress the all compressed FASTQ files automatically by recognizing their suffixes.


    2. File name: 
       Number_Target_Treat_Space_Time_Other_RepNum(_[1-2]).xxx
      [01-99]_Target_Treat_Space_Time_Other_Rep[1-9]_[1-2]_Lane[1-2].xxx.xxx 

       Number: order of samples, 01~99.                                                       (It is same for biological or technical replicates.)
       Target: such as H3K27ac, RNA_Pol_2, mRNA.                                              (It is same for biological or technical replicates.)
       Treat:  such as sham, banding, normal, EEDKO.  You also can use "NA" or "NULL".        (It is same for biological or technical replicates.)
       Space: organ, tissue or cell type.                                                     (It is same for biological or technical replicates.)
       Time: such as  adult, E12.5.                                                           (It is same for biological or technical replicates.)
       Other: other informations of this sample, You also can use "NA" or "NULL".             (It is same for biological or technical replicates.)
       RepNum: Rep[1-9], biological replicates.                                               (It is same for technical replicates.)
       [1-2]: only for paired-end sequencing. (Optional)   
       Lane[1-2]: Technical replicates. (Optional)

       All of them only can be:  "0-9",   "A-Z",   "a-z",    ".",     "-"

       For instance: 
       01_mRNA_Hand2as_CM-6_L8-I20006_shenLab_Rep1_1.fastq.gz,     
       01_mRNA_Hand2as_CM-6_L8-I20006_shenLab_Rep1_2.fastq.gz,    
       02_mm_mRNA_adult_cardiomyocyte_SRR1614298_Rep2.sra,    
       03_CM_2014NC_E12.5_RNAseq_Gata4KO_Rep1.sra, 

    3. The first step of all scripts is checking file name by using regular expression, such as: 
        $pattern = "[-.0-9A-Za-z]"
         m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])_?([1-2]?)\.fastq$/
         m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])(_[1-2])?_Lane[1-2]\.fastq\.bz2$/
         m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))(_[1-2])?(_Lane[1-2])?(\.fastq)?\.(\S+)$/

    4. Quality statistics: first by FastQC, after it finished, then run other  quality statistic softwares.

    5. For biological replicats, only "_Rep[1-9]"    is different, others are same.

    6. For techniche replicates, only "(_Lane[1-2])" is different, others are same.

    7. For Paired end files,     only "(_[1-2])"     is different, others are same.



                                                           
###Usage:                                           
     Step 1  by using RASDA-1.pl, more details by "perl  RASDA-1.pl  -h".                
     Step 2  by using RASDA-2.pl, more details by "perl  RASDA-2.pl  -h".                            
     ......                                               





