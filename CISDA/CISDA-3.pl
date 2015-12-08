#!/usr/bin/env perl5
use  strict;
use  warnings;





###################################################################################################################################################################################################
###################################################################################################################################################################################################
########## Help Infromation ##########
my $HELP_g = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use CISDA (ChIP-seq Data Analyzer), version 0.51, 2016-01-01.      
        CISDA is a Pipeline for Single-end and Paired-end ChIP-seq Data Analysis by Integrating Lots of Softwares.

        Step 3: Mapping the reads to the reference genome by using Subread, BWA and Bowtie.  
                Assess the quality of ChIP-seq reads to identify possible sequencing errors 
                or biases by using FastQC and Subread utilities.
        Usage:  
               perl  CISDA-3.pl    [-v]    [-h]    [-i inputDir]    [-o outDir]   [-MM 2]
        For instance: 
                     perl  CISDA-3.pl    -i 3-Filtered          -o 4-Mapping         -MM 2       
                     perl  CISDA-3.pl    --input 3-Filtered     --output 4-Mapping   --MaxMismatch 2 
                     perl  CISDA-3.pl    --input 3-Filtered     --output 4-Mapping   --MaxMismatch 2   >> 3-runLog.txt  2>&1
     
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
                                              Only for Subread.  default: N is 2 (50bp * 4% = 2bp).
        ------------------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines such as RASDA, MESDA and HISDA,
        please see the web site: https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies 
        and Center for Life Sciences (CLS), Peking University, China.    
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------  
';


########## Version Infromation ##########
my $version_g = "  The Third Step of CISDA (ChIP-seq Data Analyzer), version 0.51, 2016-01-01.";


########## Keys and Values ##########
if ($#ARGV   == -1) { print  "\n$HELP_g\n\n";  exit 0; }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0) { @ARGV = (@ARGV, "-h");           }       ## when the number of command argumants is odd. 
my %args = @ARGV;


########## Initialize  Variables ##########
my $input_g  = '3-Filtered';      ## This is only an initialization  value or suggesting value, not default value.
my $output_g = '4-Mapping';       ## This is only an initialization  value or suggesting value, not default value.
my $MM_g     = 2;                 ## default: $MM_g is 2


########## Available Arguments ##########
my $available = "  -v  --version    -h  --help    -i  --input    -o    --output  -MM  --MaxMismatch      ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if($available !~ m/\s$key\s/) {print  "    Cann't recognize $key !!\n";  $boole_g = 1; }
}
if($boole_g == 1) {
    print "\n    The Command Line Arguments are wrong!\n";
    print   '    Please see help message by using "perl  CISDA-3.pl  -h". ';
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
if ( !(-e $output_g))   { mkdir $output_g    ||  die; }
opendir(my $DH_input, $input_g)  ||  die;     
my @inputFiles = readdir($DH_input);
my $pattern = "[-.0-9A-Za-z]+";





print "\n\n\n\n\n##################################################################################################\n";
print("\nChecking all the necessary softwares in this step......\n");

system("subread-align   -v        >> $output_g/version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'    >> $output_g/version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output_g/version_softwares.txt   2>&1");

system("propmapped                >> $output_g/version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'    >> $output_g/version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output_g/version_softwares.txt   2>&1");

system("qualityScores             >> $output_g/version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'    >> $output_g/version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output_g/version_softwares.txt   2>&1");

system("bwa                       >> $output_g/version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'    >> $output_g/version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output_g/version_softwares.txt   2>&1");

system("bowtie  --version         >> $output_g/version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'    >> $output_g/version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output_g/version_softwares.txt   2>&1");

system("bowtie2 --version         >> $output_g/version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'    >> $output_g/version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output_g/version_softwares.txt   2>&1");





print "\n\n\n\n\n##################################################################################################\n";
print("\nChecking all the input file names......\n");
my $fileNameBool = 1;
for ( my $i=0; $i<=$#inputFiles; $i++ ) {  
        next unless $inputFiles[$i] =~ m/\.fastq$/;   
        next unless $inputFiles[$i] !~ m/^[.]/;
        next unless $inputFiles[$i] !~ m/[~]$/;
        next unless $inputFiles[$i] !~ m/^unpaired/;
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
open(seqFiles_FH, ">", "$output_g/singleEnd-pairedEnd-Files.txt")  or  die; 
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










{ ########## Start subread
print "\n\n\n\n\n##################################################################################################\n";
print "\nMapping reads to the reference genome by using Subread......\n";
my $index_name = "/home/yongp/MyProgramFiles/6-2G-HTS/3-Mapping/subread-1.5.0/bin/mm9";
my $subread = "$output_g/1-Subread";  
if ( !(-e $subread) )  { mkdir $subread || die; }

for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {  
        $pairedEnd[$i]   =~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$subread/paired-end-files.txt")  or  die;
        print  tempFH  "$end1,  $end2\n";
        system("subread-align  -T 16  -M $MM_g       --SAMoutput   -i $index_name    -r $input_g/$end1.fastq   -R  $input_g/$end2.fastq   -o  $subread/$temp.sam     -t 1    >>$subread/$temp.runLog   2>&1");       
}

for (my $i=0; $i<=$#singleEnd; $i++) {   
        $singleEnd[$i] =~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))\.fastq$/   or  die; 
        my $temp = $1; 
        system("subread-align  -T 16  -M $MM_g      --SAMoutput   -i $index_name    -r $input_g/$temp.fastq                              -o $subread/$temp.sam      -t 1    >>$subread/$temp.runLog   2>&1");       
}


my $FastQC  = "$subread/FastQC";
my $QCstat  = "$subread/QCstat";
if ( !( -e $FastQC)  )   { mkdir $FastQC    ||  die; }
if ( !( -e $QCstat)  )   { mkdir $QCstat    ||  die; }
opendir(my $DH_map, $subread) || die;     
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
    system("fastqc    --outdir $FastQC     --threads 16    --format sam    --kmers 7     $subread/$temp.sam   >>$FastQC/$temp.runLog        2>&1"); 
    system("wc  -l   $subread/$temp.sam       >> $QCstat/numberOfLines.runLog   2>&1");  
    system("echo    '\n'                      >> $QCstat/numberOfLines.runLog   2>&1");
    system("propmapped   -i $subread/$temp.sam                    -o $QCstat/$temp.prommapped      >> $QCstat/$temp.prommapped      2>&1"); 
    system("qualityScores   --SAMinput   -i $subread/$temp.sam    -o $QCstat/$temp.qualityScores   >> $QCstat/$temp.qualityScores   2>&1"); 
}

} ########## End subread










{ ########## Start bwa aln
print "\n\n\n\n\n##################################################################################################\n";
print "\n\nMapping reads to the reference genome by using bwa aln......\n";
my $bwa_index = "/home/yongp/MyProgramFiles/6-2G-HTS/3-Mapping/bwa.kit/mm9";
my $BWAaln = "$output_g/2-BWAaln";  ## shorter than 70bp
if ( !(-e $BWAaln) )  { mkdir $BWAaln || die; }

for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {
        $pairedEnd[$i]   =~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq"   eq   $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$BWAaln/paired-end-files.txt")  or  die;
        print  tempFH  "$end1,  $end2\n";
        system("bwa aln       -t 16   -f $BWAaln/$end1.sai      $bwa_index     $input_g/$end1.fastq    >>$BWAaln/$end1.runLog   2>&1");
        system("bwa aln       -t 16   -f $BWAaln/$end2.sai      $bwa_index     $input_g/$end2.fastq    >>$BWAaln/$end2.runLog   2>&1");
        system("bwa sampe             -f $BWAaln/$temp.sam      $bwa_index     $BWAaln/$end1.sai  $BWAaln/$end2.sai     $input_g/$end1.fastq  $input_g/$end2.fastq    >>$BWAaln/$temp.runLog   2>&1"); 
}

for (my $i=0; $i<=$#singleEnd; $i++) {   
        $singleEnd[$i] =~ m/^(\S+_\S+_\S+_\S+_Rep[1-9])\.fastq$/   or  die; 
        my $temp = $1; 
        system("bwa aln      -t 16   -f $BWAaln/$temp.sai    $bwa_index                          $input_g/$temp.fastq         >>$BWAaln/$temp.runLog   2>&1");
        system("bwa samse            -f $BWAaln/$temp.sam    $bwa_index     $BWAaln/$temp.sai    $input_g/$temp.fastq         >>$BWAaln/$temp.runLog   2>&1"); 
}


my $FastQC  = "$BWAaln/FastQC";
my $QCstat  = "$BWAaln/QCstat";
if ( !( -e $FastQC)  )   { mkdir $FastQC    ||  die; }
if ( !( -e $QCstat)  )   { mkdir $QCstat    ||  die; }
opendir(my $DH_map, $BWAaln) || die;     
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
    system("fastqc    --outdir $FastQC     --threads 16    --format sam    --kmers 7     $BWAaln/$temp.sam   >>$FastQC/$temp.runLog        2>&1"); 
    system("wc  -l   $BWAaln/$temp.sam        >> $QCstat/numberOfLines.runLog   2>&1");  
    system("echo    '\n'                      >> $QCstat/numberOfLines.runLog   2>&1");
    system("propmapped   -i $BWAaln/$temp.sam                    -o $QCstat/$temp.prommapped      >> $QCstat/$temp.prommapped      2>&1"); 
    system("qualityScores   --SAMinput   -i $BWAaln/$temp.sam    -o $QCstat/$temp.qualityScores   >> $QCstat/$temp.qualityScores   2>&1");   
}

}  ########## End bwa aln












{ ########## Start bwa mem
print "\n\n\n\n\n##################################################################################################\n";
print "\n\nMapping reads to the reference genome by using BWA mem......\n\n";
my $bwa_index = "/home/yongp/MyProgramFiles/6-2G-HTS/3-Mapping/bwa.kit/mm9";
my $BWAmem = "$output_g/3-BWA-mem";  ## longer  than 70bp
if ( !(-e $BWAmem) )  { mkdir $BWAmem || die; }


for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {
        $pairedEnd[$i]   =~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq"  eq  $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$BWAmem/paired-end-files.txt")  or  die;
        print  tempFH  "$end1,  $end2\n";
        system("bwa mem  -t 16   $bwa_index     $input_g/$end1.fastq  $input_g/$end2.fastq    >$BWAmem/$temp.sam"); 
}

for (my $i=0; $i<=$#singleEnd; $i++) {   
        $singleEnd[$i] =~ m/^(\S+_\S+_\S+_\S+_Rep[1-9])\.fastq$/   or  die; 
        my $temp = $1; 
        system("bwa mem  -t 16   $bwa_index     $input_g/$temp.fastq   >$BWAmem/$temp.sam");
}


my $FastQC  = "$BWAmem/FastQC";
my $QCstat  = "$BWAmem/QCstat";
if ( !( -e $FastQC)  )   { mkdir $FastQC    ||  die; }
if ( !( -e $QCstat)  )   { mkdir $QCstat    ||  die; }
opendir(my $DH_map, $BWAmem) || die;     
my @mapFiles = readdir($DH_map);


print "\n\nDetecting the quality of sam files by using FastQC......\n";
for (my $i=0; $i<=$#mapFiles; $i++) {
    next unless $mapFiles[$i] =~ m/\.sam$/;
    next unless $mapFiles[$i] !~ m/^[.]/;
    next unless $mapFiles[$i] !~ m/[~]$/;
    next unless $mapFiles[$i] !~ m/^unpaired/;
    my $temp = $mapFiles[$i]; 
    $temp =~ s/\.sam$//  ||  die; 
    system("fastqc    --outdir $FastQC     --threads 16    --format sam    --kmers 7     $BWAmem/$temp.sam   >>$FastQC/$temp.runLog        2>&1"); 
    system("wc  -l   $BWAmem/$temp.sam        >> $QCstat/numberOfLines.runLog   2>&1");  
    system("echo    '\n'                      >> $QCstat/numberOfLines.runLog   2>&1");
    system("propmapped   -i $BWAmem/$temp.sam                    -o $QCstat/$temp.prommapped       >> $QCstat/$temp.prommapped      2>&1"); 
    system("qualityScores   --SAMinput   -i $BWAmem/$temp.sam    -o $QCstat/$temp.qualityScores    >> $QCstat/$temp.qualityScores   2>&1"); 
}

}  ########## End bwa mem

















{ ########## Start Bowtie1
print "\n\n\n\n\n##################################################################################################\n";
print "\n\nMapping reads to the reference genome by using Bowtie1......\n";
my $Bowtie1_index = "/home/yongp/MyProgramFiles/6-2G-HTS/3-Mapping/bowtie-1.1.2/mm9";
my $Bowtie1 = "$output_g/4-Bowtie1";  
if ( !(-e $Bowtie1) )  { mkdir $Bowtie1 || die; }

for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {  
        $pairedEnd[$i]   =~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq"  eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$Bowtie1/paired-end-files.txt")  or  die;
        print  tempFH  "$end1,  $end2\n";
        system("bowtie   --threads 16   -q  --chunkmbs 640    --sam      --phred33-quals   $Bowtie1_index    -1 $input_g/$end1.fastq   -2 $input_g/$end2.fastq     $Bowtie1/$temp.sam    >>$Bowtie1/$temp.runLog   2>&1");       
}

for (my $i=0; $i<=$#singleEnd; $i++) {   
        $singleEnd[$i] =~ m/^(\S+_\S+_\S+_\S+_Rep[1-9])\.fastq$/   or  die; 
        my $temp = $1; 
        system("bowtie   --threads 16   -q  --chunkmbs 640    --sam      --phred33-quals   $Bowtie1_index     $input_g/$temp.fastq                                 $Bowtie1/$temp.sam    >>$Bowtie1/$temp.runLog   2>&1");                    
}


my $FastQC  = "$Bowtie1/FastQC";
my $QCstat  = "$Bowtie1/QCstat";
if ( !( -e $FastQC)  )   { mkdir $FastQC    ||  die; }
if ( !( -e $QCstat)  )   { mkdir $QCstat    ||  die; }
opendir(my $DH_map, $Bowtie1) || die;     
my @mapFiles = readdir($DH_map);


print "\n\nDetecting the quality of sam files by using FastQC......\n";
for (my $i=0; $i<=$#mapFiles; $i++) {
    next unless $mapFiles[$i] =~ m/\.sam$/;
    next unless $mapFiles[$i] !~ m/^[.]/;
    next unless $mapFiles[$i] !~ m/[~]$/;
    next unless $mapFiles[$i] !~ m/^unpaired/;
    my $temp = $mapFiles[$i]; 
    $temp =~ s/\.sam$//  ||  die; 
    system("fastqc    --outdir $FastQC     --threads 16    --format sam    --kmers 7     $Bowtie1/$temp.sam   >>$FastQC/$temp.runLog        2>&1"); 
    system("wc  -l   $Bowtie1/$temp.sam       >> $QCstat/numberOfLines.runLog   2>&1");  
    system("echo    '\n'                      >> $QCstat/numberOfLines.runLog   2>&1");
    system("propmapped   -i $Bowtie1/$temp.sam                    -o $QCstat/$temp.prommapped      >> $QCstat/$temp.prommapped      2>&1"); 
    system("qualityScores   --SAMinput   -i $Bowtie1/$temp.sam    -o $QCstat/$temp.qualityScores   >> $QCstat/$temp.qualityScores   2>&1"); 
}

} ########## End Bowtie1















{ ########## Start Bowtie2
print "\n\n\n\n\n##################################################################################################\n";
print "\n\nMapping reads to the reference genome by using Bowtie2......\n";
my $Bowtie2_index = "/home/yongp/MyProgramFiles/6-2G-HTS/3-Mapping/bowtie2-2.2.6/mm9";  
my $Bowtie2 = "$output_g/5-Bowtie2";  
if ( !(-e $Bowtie2) )  { mkdir $Bowtie2 || die; }

for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {  
        $pairedEnd[$i]   =~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$Bowtie2/paired-end-files.txt")  or  die;
        print  tempFH  "$end1,  $end2\n";
        system("bowtie2   --threads 16   -q   --phred33   --end-to-end    -x $Bowtie2_index    -1 $input_g/$end1.fastq        -2 $input_g/$end1.fastq     -S $Bowtie2/$temp.sam    >>$Bowtie2/$temp.runLog  2>&1");                                                                                                                                 
}

for (my $i=0; $i<=$#singleEnd; $i++) {   
        $singleEnd[$i] =~ m/^(\S+_\S+_\S+_\S+_Rep[1-9])\.fastq$/   or  die; 
        my $temp = $1; 
        system("bowtie2   --threads 16   -q   --phred33   --end-to-end    -x $Bowtie2_index    -U $input_g/$temp.fastq                                    -S $Bowtie2/$temp.sam    >>$Bowtie2/$temp.runLog  2>&1");                                                                
}


my $FastQC  = "$Bowtie2/FastQC";
my $QCstat  = "$Bowtie2/QCstat";
if ( !( -e $FastQC)  )   { mkdir $FastQC    ||  die; }
if ( !( -e $QCstat)  )   { mkdir $QCstat    ||  die; }
opendir(my $DH_map, $Bowtie2) || die;     
my @mapFiles = readdir($DH_map);


print "\n\nDetecting the quality of sam files by using FastQC......\n";
for (my $i=0; $i<=$#mapFiles; $i++) {
    next unless $mapFiles[$i] =~ m/\.sam$/;
    next unless $mapFiles[$i] !~ m/^[.]/;
    next unless $mapFiles[$i] !~ m/[~]$/;
    next unless $mapFiles[$i] !~ m/^unpaired/;
    my $temp = $mapFiles[$i]; 
    $temp =~ s/\.sam$//  ||  die; 
    system("fastqc    --outdir $FastQC     --threads 16    --format sam    --kmers 7     $Bowtie2/$temp.sam   >>$FastQC/$temp.runLog        2>&1"); 
    system("wc  -l   $Bowtie2/$temp.sam       >> $QCstat/numberOfLines.runLog   2>&1");  
    system("echo    '\n'                      >> $QCstat/numberOfLines.runLog   2>&1");
    system("propmapped   -i $Bowtie2/$temp.sam                    -o $QCstat/$temp.prommapped      >> $QCstat/$temp.prommapped      2>&1"); 
    system("qualityScores   --SAMinput   -i $Bowtie2/$temp.sam    -o $QCstat/$temp.qualityScores   >> $QCstat/$temp.qualityScores   2>&1"); 
}

}  ########## End Bowtie2










print "\n\n\n\n\n##################################################################################################\n";
print "\n\n        Job Done! Cheers! \n\n";





























