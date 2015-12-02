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

        Step 2: Remove adaptors and PCR primers, trim and filter the raw reads by using Trimmomatic,  
                and  assess the quality of the raw reads to identify possible sequencing errors 
                or biases by using FastQC, NGS_QC_Toolkit and FASTX-toolkit.
        Usage: 
               perl  CISDA-2.pl    [-v]    [-h]    [-i inputDir]    [-o outDir]  
        For instance: 
                     perl  CISDA-2.pl    -i 2-FASTQ          -o 3-Filtered           
                     perl  CISDA-2.pl    --input 2-FASTQ     --output 3-Filtered    
                     perl  CISDA-2.pl    --input 2-FASTQ     --output 3-Filtered    >> 2-runLog.txt  2>&1
     
        -------------------------------------------------------------------------------------------------------------------
        Optional arguments:
        -v, --version        Show version number of this program and exit.

        -h, --help           Show this help message and exit.

        Required arguments:
        -i inputDir,  --input inputDir        inputDir is the name of your input folder that contains your FASTQ files,
                                              the suffix of the FASTQ files must be ".fastq".    (no default)

        -o outDir,  --output outDir           outDir is the name of your output folder that contains running 
                                              results (fastq format) of this step.      (no default)
        ------------------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines such as RASDA, MESDA and HISDA,
        please see the web site: https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies 
        and Center for Life Sciences (CLS), Peking University, China.    
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------  
';


########## Version Infromation ##########
my $version_g = "  The Second Step of CISDA (ChIP-seq Data Analyzer), version 0.51, 2016-01-01.";


########## Keys and Values ##########
if ($#ARGV   == -1) { print  "\n$HELP_g\n\n";  exit 0; }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0) { @ARGV = (@ARGV, "-h");           }       ## when the number of command argumants is odd. 
my %args = @ARGV;


########## Initialize  Variables ##########
my $input_g  = '2-FASTQ';      ## This is only an initialization  value or suggesting value, not default value.
my $output_g = '3-Filtered';   ## This is only an initialization  value or suggesting value, not default value.


########## Available Arguments ##########
my $available = "  -v  --version    -h  --help    -i  --input    -o    --output  ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if($available !~ m/\s$key\s/) {print  "    Cann't recognize $key !!\n";  $boole_g = 1; }
}
if($boole_g == 1) {
    print "\n    The Command Line Arguments are wrong!\n";
    print   '    Please see help message by using "perl  CISDA-2.pl  -h". ';
    print "\n\n";
    exit 0;
}


########## Get Arguments ##########
if ( ( exists $args{'-v' } )  or  ( exists $args{'--version' } )  )     { print  "\n$version_g\n\n";    exit 0; }
if ( ( exists $args{'-h' } )  or  ( exists $args{'--help'    } )  )     { print  "\n$HELP_g\n\n";       exit 0; }
if ( ( exists $args{'-i' } )  or  ( exists $args{'--input'   } )  )     { ($input_g  = $args{'-i' })  or  ($input_g  = $args{'--input' });  }else{print   "\n -i or --input  is required.\n\n";   print  "\n$HELP_g\n\n";       exit 0; }                                               
if ( ( exists $args{'-o' } )  or  ( exists $args{'--output'  } )  )     { ($output_g = $args{'-o' })  or  ($output_g = $args{'--output'});  }else{print   "\n -o or --output is required.\n\n";   print  "\n$HELP_g\n\n";       exit 0; }      


########### Conditions #############
$input_g  =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";
$output_g =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";


######### Print Command Arguments to Standard Output ###########
print  "\n\n
        ################ Your Arguments ###############################
                Input  folder:  $input_g
                Output folder:  $output_g
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
my $Trimmomatic = "/home/yongp/MyProgramFiles/6-2G-HTS/2-NGSquality/Trimmomatic-0.35/trimmomatic-0.35.jar";
my $IlluQC_PRLL = "/home/yongp/MyProgramFiles/6-2G-HTS/2-NGSquality/NGSQCToolkit_v2.3.3/QC/IlluQC_PRLL.pl";
system("java  -jar    $Trimmomatic    PE      >> $output_g/version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'                >> $output_g/version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output_g/version_softwares.txt   2>&1");

system("java  -jar    $Trimmomatic    SE     >> $output_g/version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'               >> $output_g/version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output_g/version_softwares.txt   2>&1");

system("fastqc    -v             >> $output_g/version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'   >> $output_g/version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output_g/version_softwares.txt   2>&1");

system("perl  $IlluQC_PRLL  -h   >> $output_g/version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'   >> $output_g/version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output_g/version_softwares.txt   2>&1");

system("fastx_quality_stats  -h  >> $output_g/version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'   >> $output_g/version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output_g/version_softwares.txt   2>&1");

system("fastq_quality_boxplot_graph.sh   -h  >> $output_g/version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'               >> $output_g/version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output_g/version_softwares.txt   2>&1");

system("fastx_nucleotide_distribution_graph.sh  -h  >> $output_g/version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'                      >> $output_g/version_softwares.txt   2>&1");
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
        $temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])/   or  die  "wrong-1: ## $temp ##";
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





print "\n\n\n\n\n##################################################################################################\n";
print "\nFiltering the reads by using Trimmomatic ......\n";
my $readsNumFile = "$output_g/1and2-readsNumber.txt";
system("echo  FASTQ files in Folder $input_g. Need to be divided by 4:              >> $readsNumFile  2>&1");
for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {
        $pairedEnd[$i] =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])_1\.fastq$/   or  die;
        $pairedEnd[$i] =~ m/^(\S+)_1.fastq$/ or die;
        my $temp = $1;   
        my $end1 = $temp."_1.fastq";
        my $end2 = $temp."_2.fastq";
        ($end2 eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$output_g/paired-end-files.txt")  or  die;
        print  tempFH  "$end1,  $end2\n";
        system("java  -jar  $Trimmomatic  PE   -threads 16     $input_g/$end1  $input_g/$end2    $output_g/$end1  $output_g/unpaired-$end1    $output_g/$end2  $output_g/unpaired-$end2      ILLUMINACLIP:0-Other/TruSeqAdapter/All.fasta:2:30:10   LEADING:3   TRAILING:3   SLIDINGWINDOW:4:10   MINLEN:30   TOPHRED33    >>$output_g/$temp.runLog  2>&1");                       
        system("wc  -l   $input_g/$end1     >> $readsNumFile  2>&1");
        system("wc  -l   $input_g/$end2     >> $readsNumFile  2>&1");
        ##system("rm   $input_g/$end1");  
        ##system("rm   $input_g/$end2");          
}
for (my $i=0; $i<=$#singleEnd; $i++) {   
        $singleEnd[$i] =~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))\.fastq$/   or  die; 
        my $temp = $1; 
        system("java  -jar  $Trimmomatic  SE   -threads 16     $input_g/$temp.fastq  $output_g/$temp.fastq    ILLUMINACLIP:0-Other/TruSeqAdapter/All.fasta:2:30:10   LEADING:3   TRAILING:3   SLIDINGWINDOW:4:10   MINLEN:30    TOPHRED33    >>$output_g/$temp.runLog  2>&1");                    
        system("wc  -l   $input_g/$temp.fastq   >> $readsNumFile  2>&1");
        ##system("rm   $input_g/$temp.fastq");  
}
system("echo   '\n\n\n\n\n\n'     >> $readsNumFile  2>&1");





print "\n\n\n\n\n##################################################################################################\n";
print "\nDetecting single-end and paired-end FASTQ files in output folder ......\n";
opendir(my $DH_output, $output_g) || die;     
my @outputFiles = readdir($DH_output);
@singleEnd = ();
@pairedEnd = ();
open(seqFiles_FH, ">", "$output_g/singleEnd-pairedEnd-Files-thisFolder.txt")  or  die; 
for ( my $i=0; $i<=$#outputFiles; $i++ ) {     
    next unless $outputFiles[$i] =~ m/\.fastq$/;
    next unless $outputFiles[$i] !~ m/^[.]/;
    next unless $outputFiles[$i] !~ m/[~]$/;
    next unless $outputFiles[$i] !~ m/^unpaired/;
    if ($outputFiles[$i] =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])\.fastq$/) {   ## sinlge end sequencing files.
        $outputFiles[$i] =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])\.fastq$/  or  die;  
        $singleEnd[$#singleEnd+1] =  $outputFiles[$i];
        print  "\n\n        Single-end sequencing files: $outputFiles[$i]\n";
        print seqFiles_FH  "Single-end sequencing files: $outputFiles[$i]\n";
    }else{     ## paired end sequencing files.
        $outputFiles[$i] =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])_([1-2])\.fastq$/  or  die; 
        if ($outputFiles[$i] =~ m/^(\S+_\S+_\S+_\S+_Rep[1-9])_1\.fastq$/) { ## The two files of one paired sequencing sample are always side by side. 
            my $temp = $1;
            my $end1 = $temp."_1.fastq";
            my $end2 = $temp."_2.fastq";
            (-e  "$output_g/$end1")  or die;  
            (-e  "$output_g/$end2")  or die;
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
$numSingle = $#singleEnd + 1;
$numPaired = $#pairedEnd + 1;
print seqFiles_FH   "\nThere are $numSingle single-end sequencing files.\n";
print seqFiles_FH   "\nThere are $numPaired paired-end sequencing files.\n";
print     "\n\n        There are $numSingle single-end sequencing files.\n";
print     "\n\n        There are $numPaired paired-end sequencing files.\n";





my $FastQCdir       = "$output_g/FastQC";
my $FastQCdir_10mer = "$output_g/FastQC_10mer";
my $NGSQCToolkit    = "$output_g/NGSQCToolkit";
my $NGSQCToolPaired = "$output_g/NGSQCToolkit_PairedEnd";
my $FASTXtoolkit    = "$output_g/FASTXtoolkit";
if ( !( -e $FastQCdir)       )   { mkdir $FastQCdir        ||  die; }
if ( !( -e $FastQCdir_10mer) )   { mkdir $FastQCdir_10mer  ||  die; }
if ( !( -e $NGSQCToolkit)    )   { mkdir $NGSQCToolkit     ||  die; }
if ( !( -e $NGSQCToolPaired) )   { mkdir $NGSQCToolPaired  ||  die; }
if ( !( -e $FASTXtoolkit)    )   { mkdir $FASTXtoolkit     ||  die; }
  




print "\n\n\n\n\n##################################################################################################\n";
print "\nDetecting the quality of single-end FASTQ files by using FastQC......\n";
for ( my $i=0; $i<=$#singleEnd; $i++ ) {     
    my $temp = $singleEnd[$i]; 
    $temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])\.fastq$/   or  die;
    $temp =~ s/\.fastq$//  ||  die;
    system( "fastqc    --outdir $FastQCdir   --threads 16    --kmers 7     $output_g/$temp.fastq       >> $FastQCdir/$temp.runLog   2>&1" );  
} 





print "\n\n\n\n\n##################################################################################################\n";
print "\nDetecting the quality of paired-end FASTQ files by using FastQC and NGS_QC_Toolkit......\n";
for ( my $j=0; $j<=$#pairedEnd; $j=$j+2 ) {     
    my $temp1 = $pairedEnd[$j]; 
    my $temp2 = $pairedEnd[$j+1]; 
    $temp1 =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])_1\.fastq$/   or  die;
    $temp2 =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])_2\.fastq$/   or  die;
    $temp1 =~ s/\.fastq$//  ||  die;
    $temp2 =~ s/\.fastq$//  ||  die;
    system( "fastqc    --outdir $FastQCdir   --threads 16    --kmers 7     $output_g/$temp1.fastq       >> $FastQCdir/$temp1.runLog   2>&1" );  
    system( "fastqc    --outdir $FastQCdir   --threads 16    --kmers 7     $output_g/$temp2.fastq       >> $FastQCdir/$temp2.runLog   2>&1" );  
    print  seqFiles_FH  "\n\nquality statistics: $temp1,  $temp2\n";
    my $temp = $temp1;
    $temp =~ s/_1$//  ||  die;
    if ( !(-e "$NGSQCToolPaired/$temp") )   { mkdir  "$NGSQCToolPaired/$temp"  ||  die; }
    system( "perl   $IlluQC_PRLL     -pe $output_g/$temp1.fastq   $output_g/$temp2.fastq   N  A     -cpus 16   -onlyStat    -outputFolder $NGSQCToolPaired/$temp    >> $NGSQCToolPaired/$temp/$temp.runLog  2>&1" ); 
} 

    



print "\n\n\n\n\n##################################################################################################\n";
print "\nDetecting the quality of all FASTQ files by using FastQC, NGS_QC_Toolkit and FASTX-Toolkit......\n";
system("echo   FASTQ files in Folder $output_g. Need to be divided by 4:              >> $readsNumFile  2>&1");
for ( my $i=0; $i<=$#outputFiles; $i++ ) {     
    next unless $outputFiles[$i] =~ m/\.fastq$/;
    next unless $outputFiles[$i] !~ m/^[.]/;
    next unless $outputFiles[$i] !~ m/[~]$/;
    next unless $outputFiles[$i] !~ m/^unpaired/;
    my $temp = $outputFiles[$i]; 
    $temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])(_?)([1-2]?)\.fastq$/   or  die;
    $temp =~ s/\.fastq$//  ||  die;
    system( "fastqc    --outdir $FastQCdir_10mer    --threads 16    --kmers 10    $output_g/$temp.fastq       >> $FastQCdir_10mer/$temp.runLog   2>&1" );  
    if ( !(-e "$NGSQCToolkit/$temp") )   { mkdir  "$NGSQCToolkit/$temp"  ||  die; }
    system( "perl   $IlluQC_PRLL    -se $output_g/$temp.fastq   N  A    -cpus 16   -onlyStat    -outputFolder $NGSQCToolkit/$temp                >> $NGSQCToolkit/$temp/$temp.runLog   2>&1" );      
    system( "fastx_quality_stats                        -i $output_g/$temp.fastq                    -o $FASTXtoolkit/$temp.txt            >> $FASTXtoolkit/$temp.runLog   2>&1" ); 
    system( "fastq_quality_boxplot_graph.sh             -i $FASTXtoolkit/$temp.txt     -t $temp     -o $FASTXtoolkit/$temp.quality.png    >> $FASTXtoolkit/$temp.runLog   2>&1" ); 
    system( "fastx_nucleotide_distribution_graph.sh     -i $FASTXtoolkit/$temp.txt     -t $temp     -o $FASTXtoolkit/$temp.nucDis.png     >> $FASTXtoolkit/$temp.runLog   2>&1" );
    system("wc  -l   $output_g/$temp.fastq   >> $readsNumFile  2>&1");  
}
system("echo   '\n\n\n\n\n\n'     >> $readsNumFile  2>&1");





print "\n\n\n\n\n##################################################################################################\n";
print "\n\n        Job Done! Cheers! \n\n";




