In order to realize automatic NGS data processing, 
the format of file names used by the pipeline must be in fixed format for ease parsing. 


    Rules for NGS pipeline CISDA:

    0. Only SRA format and compressed FASTQ files can be used as input files. 
       And the format of file name must be fixed, such as:
                    SRA file:  [01-99]_Target_Treat_Space_Time_Other_Rep[1-9].sra

       Compressed FASTQ file:  [01-99]_Target_Treat_Space_Time_Other_Rep[1-9]_[1-2]_Lane[1-2].fastq.bz2 
                               [01-99]_Target_Treat_Space_Time_Other_Rep[1-9]_[1-2]_Lane[1-2].fastq.gz
                               [01-99]_Target_Treat_Space_Time_Other_Rep[1-9]_[1-2]_Lane[1-2].fastq.tar.gz

      Only seven compressed formats are supported, their suffixes:  
                                      ".bz2",  ".gz",  ".tar.gz",  ".tar",  ".rar",  ".xz",  ".zip"
      CISDA can uncompress the all compressed FASTQ files automatically by recognizing their suffixes.

    1. In one folder, there is only one file as ChIP-seq input.
       If there are two or more biological replicate files as ChIP-seq input, you should merge them before peak calling.
       So if there are two or more groups, you should separate them, one folder is only coresponding to one group (one ChIP-seq input file).


    2. File name: 
       Number_Target_Treat_Space_Time_Other_RepNum(_[1-2]).xxx
      [01-99]_Target_Treat_Space_Time_Other_Rep[1-9]_[1-2]_Lane[1-2].xxx.xxx 

       Number: order of samples, 01~99. (It is same for biological or technical replicates.)
       Target: such as H3K27ac, RNA_Pol_2, mRNA. For input file, it must be "Input". 
       Treat:  such as sham, banding, normal, EEDKO.  You also can use "NA" or "NULL".
       Space: organ, tissue or cell type.
       Time: such as  adult, E12.5.
       Other: other informations of this sample, You also can use "NA" or "NULL".
       RepNum: Rep[1-9], biological replicates.
       [1-2]: only for paired-end sequencing. (Optional)   
       Lane[1-2]: Technical replicates. (Optional)

       All of them only can be:  "0-9",   "A-Z",   "a-z",    ".",     "-"

       For instance: 
       01_H3K27ac_EEDKO_CM_Adult_HeLab_Rep1.fastq,     
       01_Input_ctrlKO_CM_Adult_HeLab_Rep1.fastq,    
       18_MethylC_NULL_CM_Adult_SRR1040648_Rep1_2.fastq,    
       02_RNA_NO_CM_p10_Enhancer_Rep1_Lane1.fastq.bz2, 
       01_GATA4_normal_VA_Adult_SRR1025220-Ab_Rep1.sra 
       04_GATA4_sham_VA_Adult_SRR1025223-fb_Rep2.sra  

    3. The first step of all scripts is checking file name by using regular expression, such as: 
        $pattern = "[-.0-9A-Za-z]"
         m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])_?([1-2]?)\.fastq$/
         m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])(_[1-2])?_Lane[1-2]\.fastq\.bz2$/
         m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))(_[1-2])?(_Lane[1-2])?(\.fastq)?\.(\S+)$/

    4. Quality statistics: first by FastQC, after it finished, then run other  quality statistic softwares.

    5. For biological replicats, "Number_Target_Treat_Space_Time_Other" are same, only "_Rep[1-9]" is different.

    6. For techniche replicates, only "(_Lane[1-2])" is different.

    7. For Paired end files, only "(_[1-2])" is different.







Usage:
     Step 1  by using CISDA-1.pl, more details by "perl  CISDA-1.pl  -h".
     Step 2  by using CISDA-2.pl, more details by "perl  CISDA-2.pl  -h".
     ......




Suggesting Names of Folder (Optional):
        1-rawReads, 
        2-FASTQ, 
        3-Filtered, 
        4-Mapping,       
        5-SortMapped,    
        6-FinalBAM,      
        7-OtherForamts,  
        8-MACS14,        
        9-MACS2,
        10-Q,         
        11-HOMER,        
        12-DANPOS,        
        13-DiffPeaks,     
        14-BedTools,     
        15-BEDOPS,          
        16-DeepTools,    
        17-HTSeq,         
        18-featureCount, 
        19-GO,           
        20-PathWay,      
        21-SeqFeature,     ##split fragment k-mer, GC-skew, GC-profile, DNA-flexibility, PSSM, site convervative index vectors, etc.
        22-GeneList,
        23-Conservation 


