#!/usr/bin/env perl5
use  strict;
use  warnings;
use  v5.20;





###################################################################################################################################################################################################
###################################################################################################################################################################################################
########## Help Infromation ##########
my $HELP_g = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use RASDA (RNA-Seq Data Analyzer), version 0.62, 2016-01-13.      
        RASDA is a Pipeline for Single-end and Paired-end RNA-Seq Data Analysis by Integrating Lots of Softwares.

        Step 3: Mapping the reads to the reference genome by using HISAT2, STAR2 and Subjunc.  
                Assess the quality of RNA-Seq reads to identify possible sequencing errors 
                or biases by using FastQC and Subread utilities.
        Usage:  
               perl  RASDA-3.pl    [-v]    [-h]    [-i inputDir]    [-o outDir]   [-MM 3]
        For instance: 
                     perl  RASDA-3.pl    -i 3-Filtered          -o 4-Mapping         -MM 3       
                     perl  RASDA-3.pl    --input 3-Filtered     --output 4-Mapping   --MaxMismatch 3 
                     perl  RASDA-3.pl    --input 3-Filtered     --output 4-Mapping   --MaxMismatch 3   >> RASDA-3.runLog  2>&1
     
        -------------------------------------------------------------------------------------------------------------------
        Optional arguments:
        -v, --version        Show version number of this program and exit.

        -h, --help           Show this help message and exit.

        Required arguments:
        -i inputDir,  --input inputDir        inputDir is the name of your input folder that contains your FASTQ files,
                                              the suffix of the FASTQ files must be ".fastq".    (no default)

        -o outDir,  --output outDir           outDir is the name of your output folder that contains running 
                                              results (SAM format) of this step.      (no default)

        -MM N,      --MaxMismatch N           Specify the maximum number of mis-matched bases allowed in the alignment. 
                                              Only for Subjunc.  default: N is 3.
        ------------------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines such as CISDA, MESDA and HISDA,
        please see the web site: https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies 
        and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.    
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------  
';


########## Version Infromation ##########
my $version_g = "  The Third Step of RASDA (RNA-Seq Data Analyzer), version 0.62, 2016-01-13.";


########## Keys and Values ##########
if ($#ARGV   == -1) { print  "\n$HELP_g\n\n";  exit 0; }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0) { @ARGV = (@ARGV, "-h");           }       ## when the number of command argumants is odd. 
my %args = @ARGV;


########## Initialize  Variables ##########
my $input_g  = '3-Filtered';      ## This is only an initialization  value or suggesting value, not default value.
my $output_g = '4-Mapping';       ## This is only an initialization  value or suggesting value, not default value.
my $MM_g     = 3;                 ## default: $MM_g is 3


########## Available Arguments ##########
my $available = "  -v  --version    -h  --help    -i  --input    -o    --output  -MM  --MaxMismatch      ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if($available !~ m/\s$key\s/) {print  "    Cann't recognize $key !!\n";  $boole_g = 1; }
}
if($boole_g == 1) {
    print "\n    The Command Line Arguments are wrong!\n";
    print   '    Please see help message by using "perl  RASDA-3.pl  -h". ';
    print "\n\n";
    exit 0;
}


########## Get Arguments ##########
if ( ( exists $args{'-v' } )  or  ( exists $args{'--version'      } )  )     { print  "\n$version_g\n\n";    exit 0; }
if ( ( exists $args{'-h' } )  or  ( exists $args{'--help'         } )  )     { print  "\n$HELP_g\n\n";       exit 0; }
if ( ( exists $args{'-i' } )  or  ( exists $args{'--input'        } )  )     { ($input_g  = $args{'-i' })  or  ($input_g  = $args{'--input'      });  }else{print   "\n -i or --input  is required.\n\n";   print  "\n$HELP_g\n\n";       exit 0; }                                               
if ( ( exists $args{'-o' } )  or  ( exists $args{'--output'       } )  )     { ($output_g = $args{'-o' })  or  ($output_g = $args{'--output'     });  }else{print   "\n -o or --output is required.\n\n";   print  "\n$HELP_g\n\n";       exit 0; }      
if ( ( exists $args{'-MM'} )  or  ( exists $args{'--MaxMismatch'  } )  )     { ($MM_g = $args{'-MM' })     or  ($MM_g     = $args{'--MaxMismatch'});  }


########### Conditions #############
$input_g  =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";
$output_g =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";
$MM_g     =~ m/^\d+$/   ||  die   "\n$HELP_g\n\n";


######### Print Command Arguments to Standard Output ###########
print  "\n\n
        ################ Your Arguments ###############################
                Input  folder:  $input_g
                Output folder:  $output_g
                Maximum number of mis-matched bases: $MM_g
        ###############################################################  
\n\n";


###################################################################################################################################################################################################
###################################################################################################################################################################################################





print "\n\n\n\n\n##################################################################################################\n";
print   "\nRunning......\n";
my $input2_g  = "$input_g/Results";
my $output2_g = "$output_g/Results";
if ( !(-e $input2_g) )   { mkdir $input2_g    ||  die; }
if ( !(-e $output_g) )   { mkdir $output_g    ||  die; }
if ( !(-e $output2_g))   { mkdir $output2_g   ||  die; }
opendir(my $DH_input, $input_g)  ||  die;     
my @inputFiles = readdir($DH_input);
my $pattern = "[-.0-9A-Za-z]+";





print "\n\n\n\n\n##################################################################################################\n";
print("\nChecking all the necessary softwares in this step......\n");

system("hisat2  --version         >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'    >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");

system("propmapped                >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'    >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");

system("qualityScores             >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'    >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");

system("subjunc  -v               >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'    >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");

system("STAR    --version         >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'    >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");

system("fastqc    -v              >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'    >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");







print "\n\n\n\n\n##################################################################################################\n";
print("\nChecking all the input file names......\n");
my $fileNameBool = 1;
for ( my $i=0; $i<=$#inputFiles; $i++ ) {  
        next unless $inputFiles[$i] =~ m/\.fastq$/;   
        next unless $inputFiles[$i] !~ m/^[.]/;
        next unless $inputFiles[$i] !~ m/[~]$/;
        next unless $inputFiles[$i] !~ m/^unpaired/;
        print("$inputFiles[$i] ......\n");    
        my $temp = $inputFiles[$i]; 
        $temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])/   or  die   "wrong-1: ## $temp ##";
        $temp =~ m/_(Rep[1-9])\.fastq$/  or  $temp =~ m/_(Rep[1-9])_?([1-2]?)\.fastq$/           or  die   "wrong-2: ## $temp ##";
        if($temp !~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))(_[1-2])?\.fastq$/) {
             $fileNameBool = 0;
        }
}
if($fileNameBool == 1)  {print("    All the file names are passed.\n\n");}





print "\n\n\n\n\n##################################################################################################\n";
print "\nDetecting single-end and paired-end FASTQ files in input folder......\n";
my @singleEnd = ();
my @pairedEnd = ();
open(seqFiles_FH, ">", "$output2_g/z-singleEnd-pairedEnd-Files.txt")  or  die; 
for ( my $i=0; $i<=$#inputFiles; $i++ ) {     
    next unless $inputFiles[$i] =~ m/\.fastq$/;
    next unless $inputFiles[$i] !~ m/^[.]/;
    next unless $inputFiles[$i] !~ m/[~]$/;
    next unless $inputFiles[$i] !~ m/^unpaired/;
    if ($inputFiles[$i] =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])\.fastq$/) {   ## sinlge end sequencing files.
        $inputFiles[$i] =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])\.fastq$/  or  die;  
        $singleEnd[$#singleEnd+1] =  $inputFiles[$i];
        print  "\n\n        Single-end sequencing files:  $inputFiles[$i]\n";
        print seqFiles_FH  "Single-end sequencing files: $inputFiles[$i]\n";
    }else{     ## paired end sequencing files.
        $inputFiles[$i] =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])_([1-2])\.fastq$/  or  die; 
        if ($inputFiles[$i] =~ m/^(\S+_\S+_\S+_\S+_Rep[1-9])_1\.fastq$/) { ## The two files of one paired sequencing sample are always side by side. 
            my $temp = $1;
            my $end1 = $temp."_1.fastq";
            my $end2 = $temp."_2.fastq";
            (-e  "$input_g/$end1")  or die;  
            (-e  "$input_g/$end2")  or die;
            $pairedEnd[$#pairedEnd+1] =  $end1;
            $pairedEnd[$#pairedEnd+1] =  $end2;
            print  "\n\n        Paired-end sequencing files: $end1,  $end2\n";
            print seqFiles_FH  "Paired-end sequencing files: $end1,  $end2\n";
        }
    }
}
( ($#pairedEnd+1)%2 == 0 )  or die;
print   seqFiles_FH  "\n\n\n\n\n";
print   seqFiles_FH  "All single-end sequencing files:  @singleEnd\n\n\n\n\n\n";
print   seqFiles_FH  "All paired-end sequencing files:  @pairedEnd\n\n\n\n\n\n";
print    "\n\n";
print    "\n\n        All single-end sequencing files:  @singleEnd\n\n";
print    "\n\n        All paired-end sequencing files:  @pairedEnd\n\n";
my $numSingle = $#singleEnd + 1;
my $numPaired = $#pairedEnd + 1;
print seqFiles_FH   "\nThere are $numSingle single-end sequencing files.\n";
print seqFiles_FH   "\nThere are $numPaired paired-end sequencing files.\n";
print     "\n\n        There are $numSingle single-end sequencing files.\n";
print     "\n\n        There are $numPaired paired-end sequencing files.\n";







 





{ ########## Start HISAT2
print "\n\n\n\n\n##################################################################################################\n";
print "\n\nMapping reads to the reference genome by using HISAT2......\n";
my $HISAT2_index = "/home/yongp/MyProgramFiles/6-2G-HTS/6-RNAseq/hisat2-2.0.1/mm9";
my $HISAT2   = "$output_g/1-HISAT2";  
my $HISAT2_2 = "$output_g/1-HISAT2/Results";
if ( !(-e $HISAT2)   )  { mkdir $HISAT2   || die; }
if ( !(-e $HISAT2_2) )  { mkdir $HISAT2_2 || die; }

for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {  
        print("$pairedEnd[$i] ......\n");    
        print("$pairedEnd[$i+1] ......\n");    
        $pairedEnd[$i]   =~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq"  eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$HISAT2_2/z-paired-end-files.txt")  or  die;
        print  tempFH  "$end1,  $end2\n";
        system("hisat2   --threads 16   -q  --phred33  -x $HISAT2_index    -1 $input_g/$end1.fastq   -2 $input_g/$end2.fastq   -S  $HISAT2/$temp.sam    >>$HISAT2_2/$temp.runLog   2>&1");       
}

for (my $i=0; $i<=$#singleEnd; $i++) { 
        print("$singleEnd[$i] ......\n");      
        $singleEnd[$i] =~ m/^(\S+_\S+_\S+_\S+_Rep[1-9])\.fastq$/   or  die; 
        my $temp = $1; 
        system("hisat2   --threads 16   -q  --phred33  -x $HISAT2_index     -U $input_g/$temp.fastq                            -S  $HISAT2/$temp.sam    >>$HISAT2_2/$temp.runLog   2>&1");                    
}


my $FastQC  = "$HISAT2_2/FastQC";
my $QCstat  = "$HISAT2_2/QCstat";
if ( !( -e $FastQC)  )   { mkdir $FastQC    ||  die; }
if ( !( -e $QCstat)  )   { mkdir $QCstat    ||  die; }
opendir(my $DH_map, $HISAT2) || die;     
my @mapFiles = readdir($DH_map);


print "\n\n\n\n\n##################################################################################################\n";
print "\n\nDetecting the quality of sam files by using FastQC......\n";
for (my $i=0; $i<=$#mapFiles; $i++) {
    next unless $mapFiles[$i] =~ m/\.sam$/;
    next unless $mapFiles[$i] !~ m/^[.]/;
    next unless $mapFiles[$i] !~ m/[~]$/;
    next unless $mapFiles[$i] !~ m/^unpaired/;
    my $temp = $mapFiles[$i]; 
    $temp =~ s/\.sam$//  ||  die; 
    system("fastqc    --outdir $FastQC     --threads 16    --format sam    --kmers 7     $HISAT2/$temp.sam   >>$FastQC/$temp.runLog        2>&1"); 
    system("wc  -l   $HISAT2/$temp.sam        >> $QCstat/numberOfLines.runLog   2>&1");  
    system("echo    '\n'                      >> $QCstat/numberOfLines.runLog   2>&1");
    system("propmapped   -i $HISAT2/$temp.sam                    -o $QCstat/$temp.prommapped      >> $QCstat/$temp.prommapped      2>&1"); 
    system("qualityScores   --SAMinput   -i $HISAT2/$temp.sam    -o $QCstat/$temp.qualityScores   >> $QCstat/$temp.qualityScores   2>&1"); 
}

} ########## End HISAT2














{ ########## Start Subjunc
print "\n\n\n\n\n##################################################################################################\n";
print "\nMapping reads to the reference genome by using subjunc......\n";
my $index_name = "/home/yongp/MyProgramFiles/6-2G-HTS/3-Mapping/subread-1.5.0/bin/mm9";
my $subjunc   = "$output_g/2-Subjunc";  
my $subjunc_2 = "$output_g/2-Subjunc/Results"; 
if ( !(-e $subjunc)   )  { mkdir $subjunc   || die; }
if ( !(-e $subjunc_2) )  { mkdir $subjunc_2 || die; }

for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {  
        print("$pairedEnd[$i] ......\n");    
        print("$pairedEnd[$i+1] ......\n");  
        $pairedEnd[$i]   =~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$subjunc_2/z-paired-end-files.txt")  or  die;
        print  tempFH  "$end1,  $end2\n";
        system("subjunc  -T 16  -M $MM_g       --SAMoutput   -i $index_name    -r $input_g/$end1.fastq   -R  $input_g/$end2.fastq   -o  $subjunc/$temp.sam        >>$subjunc_2/$temp.runLog   2>&1");       
}

for (my $i=0; $i<=$#singleEnd; $i++) {  
        print("$singleEnd[$i] ......\n");       
        $singleEnd[$i] =~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))\.fastq$/   or  die; 
        my $temp = $1; 
        system("subjunc  -T 16  -M $MM_g      --SAMoutput   -i $index_name     -r $input_g/$temp.fastq                               -o $subjunc/$temp.sam        >>$subjunc_2/$temp.runLog   2>&1");       
}


my $FastQC  = "$subjunc_2/FastQC";
my $QCstat  = "$subjunc_2/QCstat";
if ( !( -e $FastQC)  )   { mkdir $FastQC    ||  die; }
if ( !( -e $QCstat)  )   { mkdir $QCstat    ||  die; }
opendir(my $DH_map, $subjunc) || die;     
my @mapFiles = readdir($DH_map);  


print "\n\n\n\n\n##################################################################################################\n";
print "\nDetecting the quality of sam files by using FastQC......\n";
for (my $i=0; $i<=$#mapFiles; $i++) {
    next unless $mapFiles[$i] =~ m/\.sam$/;
    next unless $mapFiles[$i] !~ m/^[.]/;
    next unless $mapFiles[$i] !~ m/[~]$/;
    next unless $mapFiles[$i] !~ m/^unpaired/;
    my $temp = $mapFiles[$i]; 
    $temp =~ s/\.sam$//  ||  die; 
    system("fastqc    --outdir $FastQC     --threads 16    --format sam    --kmers 7     $subjunc/$temp.sam   >>$FastQC/$temp.runLog        2>&1"); 
    system("wc  -l   $subjunc/$temp.sam       >> $QCstat/numberOfLines.runLog   2>&1");  
    system("echo    '\n'                      >> $QCstat/numberOfLines.runLog   2>&1");
    system("propmapped   -i $subjunc/$temp.sam                    -o $QCstat/$temp.prommapped      >> $QCstat/$temp.prommapped      2>&1"); 
    system("qualityScores   --SAMinput   -i $subjunc/$temp.sam    -o $QCstat/$temp.qualityScores   >> $QCstat/$temp.qualityScores   2>&1"); 
}

} ########## Start Subjunc












{ ########## Start STAR
print "\n\n\n\n\n##################################################################################################\n";
print "\n\nMapping reads to the reference genome by using STAR......\n";
my $STAR_run   = ".STAR-2.5.0b/bin/Linux_x86_64/STAR";
my $STAR_index = ".STAR-2.5.0b/bin/Linux_x86_64/mm9";
my $STAR   = "$output_g/3-STAR"; 
my $STAR_2 = "$output_g/3-STAR/Results"; 
if ( !(-e $STAR)   )  { mkdir $STAR   || die; }
if ( !(-e $STAR_2) )  { mkdir $STAR_2 || die; }

for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {  
        print("$pairedEnd[$i] ......\n");    
        print("$pairedEnd[$i+1] ......\n");  
        $pairedEnd[$i]   =~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq"  eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$STAR_2/z-paired-end-files.txt")  or  die;
        print  tempFH  "$end1,  $end2\n";
        system("./$STAR_run   --runThreadN 16     --outFileNamePrefix  $STAR/$temp      --genomeDir $STAR_index      --readFilesIn $input_g/$end1.fastq    $input_g/$end2.fastq        >>$STAR_2/$temp.runLog   2>&1  ");       
}

for (my $i=0; $i<=$#singleEnd; $i++) {  
        print("$singleEnd[$i] ......\n");       
        $singleEnd[$i] =~ m/^(\S+_\S+_\S+_\S+_Rep[1-9])\.fastq$/   or  die; 
        my $temp = $1; 
        system("./$STAR_run   --runThreadN 16     --outFileNamePrefix  $STAR/$temp     --genomeDir $STAR_index      --readFilesIn $input_g/$temp.fastq                                 >>$STAR_2/$temp.runLog   2>&1  ");                    
}


system("rename  s/Aligned.out.sam/.sam/   $STAR/*Aligned.out.sam"); 


my $FastQC  = "$STAR_2/FastQC";
my $QCstat  = "$STAR_2/QCstat";
if ( !( -e $FastQC)  )   { mkdir $FastQC    ||  die; }
if ( !( -e $QCstat)  )   { mkdir $QCstat    ||  die; }
opendir(my $DH_map, $STAR) || die;     
my @mapFiles = readdir($DH_map);


print "\n\n\n\n\n##################################################################################################\n";
print "\n\nDetecting the quality of sam files by using FastQC......\n";
for (my $i=0; $i<=$#mapFiles; $i++) {
    next unless $mapFiles[$i] =~ m/\.sam$/;
    next unless $mapFiles[$i] !~ m/^[.]/;
    next unless $mapFiles[$i] !~ m/[~]$/;
    next unless $mapFiles[$i] !~ m/^unpaired/;
    my $temp = $mapFiles[$i]; 
    $temp =~ s/\.sam$//  ||  die; 
    system("fastqc    --outdir $FastQC     --threads 16    --format sam    --kmers 7     $STAR/$temp.sam   >>$FastQC/$temp.runLog        2>&1"); 
    system("wc  -l   $STAR/$temp.sam          >> $QCstat/numberOfLines.runLog   2>&1");  
    system("echo    '\n'                      >> $QCstat/numberOfLines.runLog   2>&1");
    system("propmapped   -i $STAR/$temp.sam                    -o $QCstat/$temp.prommapped      >> $QCstat/$temp.prommapped      2>&1"); 
    system("qualityScores   --SAMinput   -i $STAR/$temp.sam    -o $QCstat/$temp.qualityScores   >> $QCstat/$temp.qualityScores   2>&1"); 
}

} ########## End STAR











print "\n\n\n\n\n##################################################################################################\n";
print "\n\n        Job Done! Cheers! \n\n";





























