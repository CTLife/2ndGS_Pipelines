#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.18;
## perl5 version >= 5.18,   you can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.





###################################################################################################################################################################################################
###################################################################################################################################################################################################


########## Help Infromation ##########
my $HELP_g = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use RASDA (RNA-Seq Data Analyzer), version 0.7.1, 2016-04-11.
        RASDA is a Pipeline for Single-end and Paired-end RNA-Seq Data Analysis by Integrating Lots of Softwares.

        Step 2: Remove adaptors and PCR primers, trim and filter the raw reads by using Skewer, AdapterRemoval, Cutadapt and Trimmomatic.
                And  assess the quality of the raw reads to identify possible sequencing errors or biases by using 20 softwares:
                Rqc, FastQC, kPAL, NGS_QC_Toolkit, FASTX-toolkit, FaQCs, prinseq, SolexaQA, stats in fastqutils of NGSUtils, htseq-qa, MultiQC,
                fqtools(alastair-droop), fastq-stats in EA-Utils, Raspberry, ht2-stat in HTQC, fastqp, QC3, fastq-tools, ShortRead and seqTools.

        Usage:
               perl  RASDA-2.pl    [-v]    [-h]    [-i inputDir]    [-o outDir]
        For instance:
                     perl  RASDA-2.pl    -i 2-FASTQ          -o 3-Filtered
                     perl  RASDA-2.pl    --input 2-FASTQ     --output 3-Filtered
                     perl  RASDA-2.pl    --input 2-FASTQ     --output 3-Filtered    >> RASDA-2.runLog  2>&1
                     RASDA-2.pl    --input 2-FASTQ     --output 3-Filtered    >> RASDA-2.runLog  2>&1

        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -v, --version        Show version number of this program and exit.

        -h, --help           Show this help message and exit.

        Required arguments:
        -i inputDir,  --input inputDir        "inputDir" is the name of input folder that contains your compressed
                                              fastq files or SRA files.
                                              The suffix of the compressed fastq files will be recognized. (no default)

        -o outDir,    --output outDir         "outDir" is the name of output folder that contains your running
                                              results (fastq format) of this step.      (no default)
        -----------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines such as CISDA, MESDA and HISDA,
        please see the web site: https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies
        and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';


########## Version Infromation ##########
my $version_g = "    The Second Step of RASDA (RNA-Seq Data Analyzer), version 0.7.1, 2016-04-11.";


########## Keys and Values ##########
if ($#ARGV   == -1)   { say  "\n$HELP_g\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-h") ;       }       ## when the number of command argumants is odd.
my %args = @ARGV;


########## Initialize  Variables ##########
my $input_g  = '2-FASTQ';           ## This is only an initialization value or suggesting value, not default value.
my $output_g = '3-Filtered';        ## This is only an initialization value or suggesting value, not default value.


########## Available Arguments ##########
my $available = "      -v  --version      -h  --help      -i  --input       -o    --output      ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole_g = 1; }
}
if($boole_g == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  RASDA-2.pl  -h' \n";
    exit 0;
}


########## Get Arguments ##########
if ( ( exists $args{'-v' } )  or  ( exists $args{'--version' } )  )     { say  "\n$version_g\n";    exit 0; }
if ( ( exists $args{'-h' } )  or  ( exists $args{'--help'    } )  )     { say  "\n$HELP_g\n";       exit 0; }
if ( ( exists $args{'-i' } )  or  ( exists $args{'--input'   } )  )     { ($input_g  = $args{'-i' })  or  ($input_g  = $args{'--input' });  }else{say   "\n -i or --input  is required.\n";   say  "\n$HELP_g\n";    exit 0; }
if ( ( exists $args{'-o' } )  or  ( exists $args{'--output'  } )  )     { ($output_g = $args{'-o' })  or  ($output_g = $args{'--output'});  }else{say   "\n -o or --output is required.\n";   say  "\n$HELP_g\n";    exit 0; }


########### Conditions #############
$input_g  =~ m/^\S+$/    ||  die   "\n\n$HELP_g\n\n";
$output_g =~ m/^\S+$/    ||  die   "\n\n$HELP_g\n\n";


######### Print Command Arguments to Standard Output ###########
say  "\n
        ################ Your Arguments ###############################
                Input  Folder:  $input_g
                Output Folder:  $output_g
        ###############################################################
\n";


###################################################################################################################################################################################################
###################################################################################################################################################################################################










say   "\n\n\n\n\n\n##################################################################################################";
say   "Running ......";
sub myMakeDir  {
    my $path = $_[0];
    if ( !( -e $path) )  { mkdir $path  ||  die; }
}
my $output2_g = "$output_g/QC_Results";
&myMakeDir("$output_g");
&myMakeDir("$output2_g");
opendir(my $DH_input, $input_g)  ||  die;
my @inputFiles = readdir($DH_input);
my $pattern = "[-.0-9A-Za-z]+";
my $numCores = 4;





## These commands must be available:
say   "\n\n\n\n\n\n##################################################################################################";
say   "Checking all the necessary softwares in this step ......";
sub printVersion  {
    my $software = $_[0];
    system("echo    '##############################################################################'  >> $output2_g/VersionsOfSoftwares.txt   2>&1");
    system("echo    '#########$software'                                                              >> $output2_g/VersionsOfSoftwares.txt   2>&1");
    system("$software                                                                                 >> $output2_g/VersionsOfSoftwares.txt   2>&1");
    system("echo    '\n\n\n\n\n\n'                                                                    >> $output2_g/VersionsOfSoftwares.txt   2>&1");
}
my $Trimmomatic = "/home/yp/.MyProgramFiles/3_HTS-2G/3-OtherCommonTools/Trim/Trimmomatic-0.36/trimmomatic-0.36.jar"; 
&printVersion(" skewer  -v ");
&printVersion(" AdapterRemoval  -v ");   
&printVersion(" cutadapt  --version ");
&printVersion(" java  -jar  $Trimmomatic  -version");
&printVersion(" fastqc    -v ");
&printVersion(" kpal  -v ");
&printVersion(" IlluQC_PRLL.pl  -h ");
&printVersion(" fastx_quality_stats  -h ");
&printVersion(" FaQCs.pl   --version ");
&printVersion(" prinseq-lite.pl  -version ");
&printVersion(" fastq-stats -h "); ## in EA-Utils.
&printVersion(" SolexaQA ");
&printVersion(" raspberry ");
&printVersion(" ht2-stat  -v ");
&printVersion(" fastqp    -h ");
&printVersion(" fastqutils ");
&printVersion(" htseq-qa  -h ");
&printVersion(" multiqc  --version ");
&printVersion(" fqtools  -v ");
&printVersion(" qc3.pl  -h ");
&printVersion(" fastq-uniq  --version ");  ## in fastq-tools.





say   "\n\n\n\n\n\n##################################################################################################";
say   "Checking all the input file names ......";
my @groupFiles = ();
my $fileNameBool = 1;
for ( my $i=0; $i<=$#inputFiles; $i++ ) {
        next unless $inputFiles[$i] !~ m/^[.]/;
        next unless $inputFiles[$i] !~ m/[~]$/;
        next unless $inputFiles[$i] !~ m/^QC_Results$/;
        next unless $inputFiles[$i] =~ m/\.fastq$/;
        say   "\t......$inputFiles[$i]" ;
        my $temp = $inputFiles[$i]; 
        $groupFiles[++$#groupFiles] = $inputFiles[$i];  
        $temp =~ m/^(\d+)_($pattern)_(Rep[1-9])/   or  die   "wrong-1: ## $temp ##";
        $temp =~ m/_(Rep[1-9])\.fastq$/  or  $temp =~ m/_(Rep[1-9])_?([1-2]?)\.fastq$/   or  die   "wrong-2: ## $temp ##";
        if($temp !~ m/^((\d+)_($pattern)_(Rep[1-9]))(_[1-2])?\.fastq$/) {
             $fileNameBool = 0;
        }
}
if($fileNameBool == 1)  { say    "\n\t\tAll the file names are passed.\n";  }
@groupFiles = sort(@groupFiles);
my $numGroup = 0;
my $noteGroup = 0;
for ( my $i=0; $i<=$#groupFiles; $i++ ) {
    $groupFiles[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])/  or  die;
    my $n1 = $1;
    $n1>=1  or  die;
    if($noteGroup != $n1) {say "\n\t\tGroup $n1:";  $numGroup++; }
    say  "\t\t\t$groupFiles[$i]";
    $noteGroup = $n1;
}
say  "\n\t\tThere are $numGroup groups.";





say   "\n\n\n\n\n\n##################################################################################################";
say   "Detecting single-end and paired-end FASTQ files in input folder ......";
my @singleEnd = ();
my @pairedEnd = ();
open(seqFiles_FH, ">", "$output2_g/singleEnd-pairedEnd-Files.txt")  or  die; 
for ( my $i=0; $i<=$#inputFiles; $i++ ) {     
    next unless $inputFiles[$i] =~ m/\.fastq$/;
    next unless $inputFiles[$i] !~ m/^[.]/;
    next unless $inputFiles[$i] !~ m/[~]$/;
    next unless $inputFiles[$i] !~ m/^unpaired/;
    say    "\t......$inputFiles[$i]"; 
    if ($inputFiles[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])\.fastq$/) {   ## sinlge end sequencing files.
        $inputFiles[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])\.fastq$/  or  die;  
        $singleEnd[$#singleEnd+1] =  $inputFiles[$i];
        say  "\t\t\t\tSingle-end sequencing files:  $inputFiles[$i]\n";
        say seqFiles_FH  "Single-end sequencing files: $inputFiles[$i]";
    }else{     ## paired end sequencing files.
        $inputFiles[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])_([1-2])\.fastq$/  or  die; 
        if ($inputFiles[$i] =~ m/((^\d+)_($pattern)_Rep[1-9])_1\.fastq$/) { ## The two files of one paired sequencing sample are always side by side. 
            my $temp = $1;
            my $end1 = $temp."_1.fastq";
            my $end2 = $temp."_2.fastq";
            (-e  "$input_g/$end1")  or die;  
            (-e  "$input_g/$end2")  or die;
            $pairedEnd[$#pairedEnd+1] =  $end1;
            $pairedEnd[$#pairedEnd+1] =  $end2;
            say  "\t\t\t\tPaired-end sequencing files: $end1,  $end2\n";
            say seqFiles_FH  "Paired-end sequencing files: $end1,  $end2";
        }
    }
}
( ($#pairedEnd+1)%2 == 0 )  or die;
say   seqFiles_FH  "\n\n\n\n\n";
say   seqFiles_FH  "All single-end sequencing files:@singleEnd\n\n\n";
say   seqFiles_FH  "All paired-end sequencing files:@pairedEnd\n\n\n";
say    "\t\t\t\tAll single-end sequencing files:@singleEnd\n\n";
say    "\t\t\t\tAll paired-end sequencing files:@pairedEnd\n\n";
my $numSingle = $#singleEnd + 1;
my $numPaired = $#pairedEnd + 1;
say seqFiles_FH   "\nThere are $numSingle single-end sequencing files.\n";
say seqFiles_FH   "\nThere are $numPaired paired-end sequencing files.\n";
say     "\t\t\t\tThere are $numSingle single-end sequencing files.\n";
say     "\t\t\t\tThere are $numPaired paired-end sequencing files.\n";





say    "\n\n\n\n\n\n##################################################################################################";
say    "Filtering the reads by using Trimmomatic ......";
my $Trimmomatic_dir  = "$output_g/1-Trimmomatic";
&myMakeDir("$Trimmomatic_dir");   
for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {
        $pairedEnd[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])_1\.fastq$/   or  die;
        $pairedEnd[$i] =~ m/^(\S+)_1.fastq$/ or die;
        my $temp = $1;   
        my $end1 = $temp."_1.fastq";
        my $end2 = $temp."_2.fastq";
        say   "\t......$end1";
        say   "\t......$end2";
        ($end2 eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$Trimmomatic_dir/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("java  -jar  $Trimmomatic  PE   -threads $numCores     $input_g/$end1  $input_g/$end2    $Trimmomatic_dir/$end1  $Trimmomatic_dir/unpaired-$end1    $Trimmomatic_dir/$end2  $Trimmomatic_dir/unpaired-$end2      ILLUMINACLIP:0-Other/TruSeqAdapter.fasta:2:30:10   LEADING:3   TRAILING:3   SLIDINGWINDOW:4:10   MINLEN:30   TOPHRED33    >>$Trimmomatic_dir/$temp.runLog  2>&1");                                                   
}
for (my $i=0; $i<=$#singleEnd; $i++) {   
        say   "\t......$singleEnd[$i]" ;
        $singleEnd[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.fastq$/   or  die; 
        my $temp = $1; 
        system("java  -jar  $Trimmomatic  SE   -threads $numCores     $input_g/$temp.fastq  $Trimmomatic_dir/$temp.fastq    ILLUMINACLIP:0-Other/TruSeqAdapter.fasta:2:30:10   LEADING:3   TRAILING:3   SLIDINGWINDOW:4:10   MINLEN:30    TOPHRED33    >>$Trimmomatic_dir/$temp.runLog  2>&1");                    
}





say    "\n\n\n\n\n\n##################################################################################################";
say    "Filtering the reads by using Skewer ......";
my $Skewer_dir  = "$output_g/2-Skewer";
&myMakeDir("$Skewer_dir");   
for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {
        $pairedEnd[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])_1\.fastq$/   or  die;
        $pairedEnd[$i] =~ m/^(\S+)_1.fastq$/ or die;
        my $temp = $1;   
        my $end1 = $temp."_1.fastq";
        my $end2 = $temp."_2.fastq";
        say   "\t......$end1";
        say   "\t......$end2";
        ($end2 eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$Skewer_dir/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("skewer   -x 0-Other/TruSeqAdapter.fasta   -y 0-Other/TruSeqAdapter.fasta   -q 10  -l 30  -n  -u  -o $Skewer_dir/$temp --threads $numCores   -m pe     $Trimmomatic_dir/$end1   $Trimmomatic_dir/$end2  >>$Skewer_dir/$temp.runLog  2>&1");                                                   
}
for (my $i=0; $i<=$#singleEnd; $i++) {   
        say   "\t......$singleEnd[$i]" ;
        $singleEnd[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.fastq$/   or  die; 
        my $temp = $1; 
        system("skewer   -x 0-Other/TruSeqAdapter.fasta    -q 10  -l 30  -n    -o $Skewer_dir/$temp   --threads $numCores     $Trimmomatic_dir/$temp.fastq    >>$Skewer_dir/$temp.runLog  2>&1");                                                   
}
system( "rm  $Trimmomatic_dir/*.fastq" );





say    "\n\n\n\n\n\n##################################################################################################";
say    "Filtering the reads by using AdapterRemoval ......";
my $AdapterRemoval_dir  = "$output_g/3-AdapterRemoval";
&myMakeDir("$AdapterRemoval_dir");   
for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {
        $pairedEnd[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])_1\.fastq$/   or  die;
        $pairedEnd[$i] =~ m/^(\S+)_1.fastq$/ or die;
        my $temp = $1;   
        my $end1 = $temp."_1.fastq";
        my $end2 = $temp."_2.fastq";
        say   "\t......$end1";
        say   "\t......$end2";
        ($end2 eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$AdapterRemoval_dir/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("AdapterRemoval  --file1 $Skewer_dir/$temp-trimmed-pair1.fastq  --file2 $Skewer_dir/$temp-trimmed-pair2.fastq    --basename $AdapterRemoval_dir/$temp    --trimns    --trimqualities    --minquality 10   --minlength 30    --threads $numCores    >>$AdapterRemoval_dir/$temp.runLog  2>&1");    
}
for (my $i=0; $i<=$#singleEnd; $i++) {   
        say   "\t......$singleEnd[$i]" ;
        $singleEnd[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.fastq$/   or  die; 
        my $temp = $1; 
        system("AdapterRemoval  --file1 $Skewer_dir/$temp-trimmed.fastq    --basename $AdapterRemoval_dir/$temp    --trimns    --trimqualities    --minquality 10   --minlength 30    --threads $numCores    >>$AdapterRemoval_dir/$temp.runLog  2>&1");    
}
system( "rename  -v   s/\.truncated/.fastq/    $AdapterRemoval_dir/*.truncated  " );
system( "rm  $Skewer_dir/*.fastq" );





say    "\n\n\n\n\n\n##################################################################################################";
say    "Filtering the reads by using Cutadapt ......";
my $Cutadapt_dir  = "$output_g/4-Cutadapt";
&myMakeDir("$Cutadapt_dir");   
for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {
        $pairedEnd[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])_1\.fastq$/   or  die;
        $pairedEnd[$i] =~ m/^(\S+)_1.fastq$/ or die;
        my $temp = $1;   
        my $end1 = $temp."_1.fastq";
        my $end2 = $temp."_2.fastq";
        say   "\t......$end1";
        say   "\t......$end2";
        ($end2 eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$Cutadapt_dir/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("cutadapt  --adapter=0-Other/TruSeqAdapter.fasta  --front=0-Other/TruSeqAdapter.fasta  --anywhere=0-Other/TruSeqAdapter.fasta   -A 0-Other/TruSeqAdapter.fasta  -G 0-Other/TruSeqAdapter.fasta   -B 0-Other/TruSeqAdapter.fasta    --quality-cutoff=10,10  --minimum-length=30   --max-n=0.1     -o $Cutadapt_dir/$end1 -p $Cutadapt_dir/$end2   $AdapterRemoval_dir/$temp.pair1.fastq  $AdapterRemoval_dir/$temp.pair2.fastq  >>$Cutadapt_dir/$temp.runLog  2>&1");                                                            
}
for (my $i=0; $i<=$#singleEnd; $i++) {   
        say   "\t......$singleEnd[$i]" ;
        $singleEnd[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.fastq$/   or  die; 
        my $temp = $1; 
        system("cutadapt  --adapter=0-Other/TruSeqAdapter.fasta  --front=0-Other/TruSeqAdapter.fasta  --anywhere=0-Other/TruSeqAdapter.fasta       --quality-cutoff=10,10  --minimum-length=30   --max-n=0.1     -o $Cutadapt_dir/$temp.fastq   $AdapterRemoval_dir/$temp.fastq   >>$Cutadapt_dir/$temp.runLog  2>&1");                                                            
}
system( "rm  $AdapterRemoval_dir/*.fastq" );
system( "rm  $AdapterRemoval_dir/*.discarded" );





say    "\n\n\n\n\n\n##################################################################################################";
say    "Filtering the reads by using Trimmomatic ......"; 
for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {
        $pairedEnd[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])_1\.fastq$/   or  die;
        $pairedEnd[$i] =~ m/^(\S+)_1.fastq$/ or die;
        my $temp = $1;   
        my $end1 = $temp."_1.fastq";
        my $end2 = $temp."_2.fastq";
        say   "\t......$end1";
        say   "\t......$end2";
        ($end2 eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$Trimmomatic_dir/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("java  -jar  $Trimmomatic  PE   -threads $numCores     $Cutadapt_dir/$end1  $Cutadapt_dir/$end2    $output_g/$end1  $output_g/unpaired-$end1    $output_g/$end2  $output_g/unpaired-$end2      ILLUMINACLIP:0-Other/TruSeqAdapter.fasta:2:30:10   LEADING:3   TRAILING:3   SLIDINGWINDOW:4:10   MINLEN:30   TOPHRED33    >>$output2_g/$temp.runLog  2>&1");                                                   
}
for (my $i=0; $i<=$#singleEnd; $i++) {   
        say   "\t......$singleEnd[$i]" ;
        $singleEnd[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.fastq$/   or  die; 
        my $temp = $1; 
        system("java  -jar  $Trimmomatic  SE   -threads $numCores     $Cutadapt_dir/$temp.fastq  $output_g/$temp.fastq    ILLUMINACLIP:0-Other/TruSeqAdapter.fasta:2:30:10   LEADING:3   TRAILING:3   SLIDINGWINDOW:4:10   MINLEN:30    TOPHRED33    >>$output2_g/$temp.runLog  2>&1");                    
}
system( "rm  $Cutadapt_dir/*.fastq" );   
system( "rm  $output_g/unpaired-*.fastq" );





say   "\n\n\n\n\n\n##################################################################################################";
say   "Detecting single-end and paired-end FASTQ files in output folder ......";
opendir(my $DH_output, $output_g) || die;     
my @outputFiles = readdir($DH_output);
@singleEnd = ();
@pairedEnd = ();
open(seqFiles_FH, ">", "$output2_g/singleEnd-pairedEnd-Files-thisFolder.txt")  or  die; 
for ( my $i=0; $i<=$#outputFiles; $i++ ) {     
    next unless $outputFiles[$i] =~ m/\.fastq$/;
    next unless $outputFiles[$i] !~ m/^[.]/;
    next unless $outputFiles[$i] !~ m/[~]$/;
    next unless $outputFiles[$i] !~ m/^unpaired/;
    say   "\t......$outputFiles[$i]";
    if ($outputFiles[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])\.fastq$/) {   ## sinlge end sequencing files.
        $outputFiles[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])\.fastq$/  or  die;  
        $singleEnd[$#singleEnd+1] =  $outputFiles[$i];
        say  "\t\t\t\tSingle-end sequencing files: $outputFiles[$i]\n";
        say seqFiles_FH  "Single-end sequencing files: $outputFiles[$i]\n";
    }else{     ## paired end sequencing files.
        $outputFiles[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])_([1-2])\.fastq$/  or  die; 
        if ($outputFiles[$i] =~ m/^((^\d+)_($pattern)_Rep[1-9])_1\.fastq$/) { ## The two files of one paired sequencing sample are always side by side. 
            my $temp = $1;
            my $end1 = $temp."_1.fastq";
            my $end2 = $temp."_2.fastq";
            (-e  "$output_g/$end1")  or die;  
            (-e  "$output_g/$end2")  or die;
            $pairedEnd[$#pairedEnd+1] =  $end1;
            $pairedEnd[$#pairedEnd+1] =  $end2;
            say  "\t\t\t\tPaired-end sequencing files: $end1,  $end2\n";
            say seqFiles_FH  "Paired-end sequencing files: $end1,  $end2\n";
        }
    }
}
( ($#pairedEnd+1)%2 == 0 )  or die;
say   seqFiles_FH  "\n\n\n\n\n";
say   seqFiles_FH  "All single-end sequencing files:@singleEnd\n\n\n";
say   seqFiles_FH  "All paired-end sequencing files:@pairedEnd\n\n\n";
say    "\t\t\t\tAll single-end sequencing files:@singleEnd\n\n";
say    "\t\t\t\tAll paired-end sequencing files:@pairedEnd\n\n";
$numSingle = $#singleEnd + 1;
$numPaired = $#pairedEnd + 1;
say seqFiles_FH   "\nThere are $numSingle single-end sequencing files.\n";
say seqFiles_FH   "\nThere are $numPaired paired-end sequencing files.\n";
say     "\t\t\t\tThere are $numSingle single-end sequencing files.\n";
say     "\t\t\t\tThere are $numPaired paired-end sequencing files.\n";





my $R_QC          = "$output2_g/1_R_QC";
my $FastQC        = "$output2_g/2_FastQC";
my $kPAL          = "$output2_g/3_kPAL";
my $NGSQC_Toolkit = "$output2_g/4_NGSQC_Toolkit";
my $FASTX_Toolkit = "$output2_g/5_FASTX_Toolkit";
my $FaQCs         = "$output2_g/6_FaQCs";
my $PRINSEQ       = "$output2_g/7_PRINSEQ";
my $EA_Utils      = "$output2_g/8_EA_Utils";
my $SolexaQA      = "$output2_g/9_SolexaQA";
my $Raspberry     = "$output2_g/10_Raspberry";
my $HTQC          = "$output2_g/11_HTQC";
my $Fastqp        = "$output2_g/12_Fastqp";
my $ShortRead     = "$output2_g/13_ShortRead";
my $seqTools      = "$output2_g/14_seqTools";
my $FastqUtils    = "$output2_g/15_FastqUtils";
my $HTSeqQA       = "$output2_g/16_HTSeqQA";
my $MultiQC       = "$output2_g/17_MultiQC";
my $FqTools       = "$output2_g/18_FqTools";
my $QC3           = "$output2_g/19_QC3";
my $FastqTools    = "$output2_g/20_FastqTools";
&myMakeDir("$R_QC");
&myMakeDir("$FastQC");
&myMakeDir("$kPAL");
&myMakeDir("$NGSQC_Toolkit");
&myMakeDir("$FASTX_Toolkit");
&myMakeDir("$FaQCs");
&myMakeDir("$PRINSEQ");
&myMakeDir("$EA_Utils");
&myMakeDir("$SolexaQA");
&myMakeDir("$Raspberry");
&myMakeDir("$HTQC");
&myMakeDir("$Fastqp");
&myMakeDir("$ShortRead");
&myMakeDir("$seqTools");
&myMakeDir("$FastqUtils");
&myMakeDir("$HTSeqQA");
&myMakeDir("$MultiQC");
&myMakeDir("$FqTools");
&myMakeDir("$QC3");
&myMakeDir("$FastqTools");





say   "\n\n\n\n\n\n##################################################################################################";
say   "Detecting the quality of all FASTQ files by using FastQC and multiqc ......";
for ( my $i=0; $i<=$#outputFiles; $i++ ) {
    next unless $outputFiles[$i] =~ m/\.fastq$/;
    next unless $outputFiles[$i] !~ m/^[.]/;
    next unless $outputFiles[$i] !~ m/[~]$/;
    my $temp = $outputFiles[$i];
    say    "\t......$temp";
    $temp =~ m/^(\d+)_($pattern)_(Rep[1-9])_?([1-2]?)\.fastq$/   or  die;
    $temp =~ s/\.fastq$//  ||  die;
    system( "fastqc    --outdir $FastQC    --threads $numCores  --format fastq   --kmers 7    $output_g/$temp.fastq       >> $FastQC/$temp.runLog     2>&1" );
}
system( "multiqc  --outdir $MultiQC          $FastQC/*_fastqc.zip      >> $MultiQC/multiqc.fastqc.runLog   2>&1" );





say   "\n\n\n\n\n\n##################################################################################################";
say   "Detecting the quality of all FASTQ files by using fastq-stats in EA-Utils, SolexaQA, raspberry, fastqp, FASTX-toolkit, fastqutils, htseq-qa, fqtools, fastq-tools ......";
for ( my $i=0; $i<=$#outputFiles; $i++ ) {
    next unless $outputFiles[$i] =~ m/\.fastq$/;
    next unless $outputFiles[$i] !~ m/^[.]/;
    next unless $outputFiles[$i] !~ m/[~]$/;
    my $temp = $outputFiles[$i];
    say    "\t......$temp";
    $temp =~ m/^(\d+)_($pattern)_(Rep[1-9])_?([1-2]?)\.fastq$/   or  die;
    $temp =~ s/\.fastq$//  ||  die;
    &myMakeDir("$FASTX_Toolkit/$temp");
    &myMakeDir("$EA_Utils/$temp");
    &myMakeDir("$SolexaQA/$temp");
    &myMakeDir("$FastqTools/$temp");
    system( "fastx_quality_stats                        -i $output_g/$temp.fastq                           -o $FASTX_Toolkit/$temp/$temp.txt            >> $FASTX_Toolkit/$temp.runLog   2>&1" );
    system( "fastq_quality_boxplot_graph.sh             -i $FASTX_Toolkit/$temp/$temp.txt     -t $temp     -o $FASTX_Toolkit/$temp/$temp.quality.png    >> $FASTX_Toolkit/$temp.runLog   2>&1" );
    system( "fastx_nucleotide_distribution_graph.sh     -i $FASTX_Toolkit/$temp/$temp.txt     -t $temp     -o $FASTX_Toolkit/$temp/$temp.nucDis.png     >> $FASTX_Toolkit/$temp.runLog   2>&1" );
    my $temp2 = $temp;
    $temp2 =~ s/\.//g;
    system( "fastq_to_fasta   -v -n   -i $output_g/$temp.fastq    -o $output_g/$temp2.fasta                                                         >> $FASTX_Toolkit/$temp.runLog   2>&1" );
    system( "echo     -e    '\n\nfastx_clipper:\n'   >> $FASTX_Toolkit/$temp.runLog  ");
    system( "fastx_clipper    -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA   -l 30   -n  -v    -i $output_g/$temp2.fasta   -o $output_g/$temp.clipped.fa   >> $FASTX_Toolkit/$temp.runLog   2>&1" );
    system( "echo     -e    '\n\nfastx_trimmer:\n'   >> $FASTX_Toolkit/$temp.runLog  ");
    system( "fastx_trimmer   -v   -m 30    -i $output_g/$temp.clipped.fa   -o $output_g/$temp.trimmed.fa                                            >> $FASTX_Toolkit/$temp.runLog   2>&1" );
    system( "echo     -e    '\n\nfastx_collapser:\n'   >> $FASTX_Toolkit/$temp.runLog  ");
    system( "fastx_collapser -v   -i $output_g/$temp.trimmed.fa    -o $output_g/$temp.collapsed.fa                                                  >> $FASTX_Toolkit/$temp.runLog   2>&1" );
    system( "fasta_clipping_histogram.pl    $output_g/$temp.collapsed.fa     $FASTX_Toolkit/$temp.collapsed.png                                     >> $FASTX_Toolkit/$temp.runLog   2>&1" );
    system( "rm  $output_g/*.fa" );

    system( "fastq-stats  -d  -s 100   -x $EA_Utils/$temp/$temp.fastxStatistics    -b $EA_Utils/$temp/$temp.baseBreakdown    -L $EA_Utils/$temp/$temp.lengthCounts  $output_g/$temp.fastq   >> $EA_Utils/$temp.runLog   2>&1" );
    system( "SolexaQA  analysis   $output_g/$temp.fastq    --minmax    --directory $SolexaQA/$temp       >> $SolexaQA/$temp.runLog    2>&1" );
    system( "raspberry  -t $numCores     $output_g/$temp.fastq                                           >> $Raspberry/$temp.runLog   2>&1" );
    system( "mv  $output_g/*.rlen   $Raspberry" );
    system( "fastqp    --nreads 20000000   --kmer 5    --output $Fastqp/$temp  --type fastq   --median-qual 30     $output_g/$temp.fastq     >> $Fastqp/$temp.runLog     2>&1 " );
    system( "fastqutils   stats  $output_g/$temp.fastq     >> $FastqUtils/$temp.runLog     2>&1 " );
    system( "htseq-qa   --type=fastq    --outfile=$HTSeqQA/$temp.pdf     $output_g/$temp.fastq    >> $HTSeqQA/$temp.runLog     2>&1 " );

    system( "fqtools   count        $output_g/$temp.fastq    >> $FqTools/$temp.runLog     2>&1 " );
    system( "fqtools   basetab  -a  $output_g/$temp.fastq    >> $FqTools/$temp.runLog     2>&1 " );
    system( "fqtools   qualtab      $output_g/$temp.fastq    >> $FqTools/$temp.runLog     2>&1 " );
    system( "fqtools   lengthtab    $output_g/$temp.fastq    >> $FqTools/$temp.runLog     2>&1 " );

    system( "fastq-kmers -k 1  $output_g/$temp.fastq    >> $FastqTools/$temp/$temp.1mer     2>&1 " );
    system( "fastq-kmers -k 2  $output_g/$temp.fastq    >> $FastqTools/$temp/$temp.2mer     2>&1 " );
    system( "fastq-kmers -k 3  $output_g/$temp.fastq    >> $FastqTools/$temp/$temp.3mer     2>&1 " );
    system( "fastq-kmers -k 4  $output_g/$temp.fastq    >> $FastqTools/$temp/$temp.4mer     2>&1 " );
    system( "fastq-qual        $output_g/$temp.fastq    >> $FastqTools/$temp/$temp.qual     2>&1 " );
    system( "fastq-uniq  --verbose  $output_g/$temp.fastq          >> $FastqTools/$temp/$temp.uniq     2>&1 " );
}





say   "\n\n\n\n\n\n##################################################################################################";
say   "Detecting the quality of single-end FASTQ files by using NGS_QC_Toolkit, FaQCs, prinseq, ht2-stat in HTQC ......";
for ( my $i=0; $i<=$#singleEnd; $i++ ) {
    my $temp = $singleEnd[$i];
    say   "\t......$singleEnd[$i]";
    $temp =~ m/^(\d+)_($pattern)_(Rep[1-9])\.fastq$/   or  die;
    $temp =~ s/\.fastq$//  ||  die;
    &myMakeDir("$NGSQC_Toolkit/$temp");
    &myMakeDir("$FaQCs/$temp");
    &myMakeDir("$PRINSEQ/$temp");
    &myMakeDir("$HTQC/$temp");
    system( "echo         '$output_g/$temp.fastq'   >> $QC3/fileLists  ");
    system( "IlluQC_PRLL.pl      -se $output_g/$temp.fastq     N    A      -cpus $numCores  -onlyStat     -outputFolder $NGSQC_Toolkit/$temp          >> $NGSQC_Toolkit/$temp.runLog   2>&1" );
    system( "FaQCs.pl     -prefix $temp     -t $numCores      -qc_only  -min_L 30    -d $FaQCs/$temp          -u $output_g/$temp.fastq                         >> $FaQCs/$temp.runLog          2>&1" );
    system( "ht2-stat  -i $output_g/$temp.fastq    -o  $HTQC/$temp   --se     --encode sanger   --threads $numCores  >> $HTQC/$temp.runLog   2>&1" );

    system( "prinseq-lite.pl  -out_format 1 -verbose  -fastq $output_g/$temp.fastq    -graph_data $PRINSEQ/$temp/$temp.gd      >> $PRINSEQ/$temp.runLog      2>&1" );
    system( "prinseq-graphs.pl   -i $PRINSEQ/$temp/$temp.gd   -png_all    -o $PRINSEQ/$temp/$temp                              >> $PRINSEQ/$temp.runLog      2>&1" );
    system( "prinseq-graphs.pl   -i $PRINSEQ/$temp/$temp.gd   -html_all   -o $PRINSEQ/$temp/$temp                              >> $PRINSEQ/$temp.runLog      2>&1" );
    system( "prinseq-lite.pl  -out_format 1  -verbose  -fastq $output_g/$temp.fastq    -stats_all                              >> $PRINSEQ/$temp.stats_all   2>&1" );
    system( "rm   $output_g/*_prinseq_* " );
}





say   "\n\n\n\n\n\n##################################################################################################";
say   "Detecting the quality of paired-end FASTQ files by using NGS_QC_Toolkit, FaQCs, prinseq, ht2-stat in HTQC ......";
for ( my $j=0; $j<=$#pairedEnd; $j=$j+2 ) {
    my $temp1 = $pairedEnd[$j];
    my $temp2 = $pairedEnd[$j+1];
    my $temp  = $temp2;
    say   "\t......$temp1";
    say   "\t......$temp2";
    $temp1 =~ m/^(\d+)_($pattern)_(Rep[1-9])_1\.fastq$/   or  die;
    $temp2 =~ m/^(\d+)_($pattern)_(Rep[1-9])_2\.fastq$/   or  die;
    $temp1 =~ s/\.fastq$//    ||  die;
    $temp2 =~ s/\.fastq$//    ||  die;
    $temp  =~ s/_2\.fastq$//  ||  die;
    &myMakeDir("$NGSQC_Toolkit/$temp");
    &myMakeDir("$FaQCs/$temp");
    &myMakeDir("$PRINSEQ/$temp");
    &myMakeDir("$HTQC/$temp");
    system( "echo         '$output_g/$temp1.fastq'  >> $QC3/fileLists  ");
    system( "echo         '$output_g/$temp2.fastq'  >> $QC3/fileLists  ");
    system( "IlluQC_PRLL.pl      -pe $output_g/$temp1.fastq  $output_g/$temp2.fastq   N    A      -cpus $numCores     -onlyStat          -outputFolder $NGSQC_Toolkit/$temp           >> $NGSQC_Toolkit/$temp.runLog   2>&1" );
    system( "FaQCs.pl     -prefix $temp     -t $numCores      -qc_only  -min_L 30    -d $FaQCs/$temp          -p $output_g/$temp1.fastq  $output_g/$temp2.fastq            >> $FaQCs/$temp.runLog          2>&1" );
    system( "ht2-stat  -i $output_g/$temp1.fastq  $output_g/$temp2.fastq   -o  $HTQC/$temp    --pe     --encode sanger   --threads 6  >> $HTQC/$temp.runLog   2>&1" );

    system( "prinseq-lite.pl   -out_format 1   -verbose  -fastq $output_g/$temp1.fastq  -fastq2 $output_g/$temp2.fastq   -graph_data $PRINSEQ/$temp/$temp.gd      >> $PRINSEQ/$temp.runLog   2>&1" );
    system( "prinseq-graphs.pl   -i $PRINSEQ/$temp/$temp.gd   -png_all    -o $PRINSEQ/$temp/$temp                 >> $PRINSEQ/$temp.runLog   2>&1" );
    system( "prinseq-graphs.pl   -i $PRINSEQ/$temp/$temp.gd   -html_all   -o $PRINSEQ/$temp/$temp                 >> $PRINSEQ/$temp.runLog   2>&1" );
    system( "prinseq-lite.pl   -out_format 1  -verbose  -fastq $output_g/$temp.fastq    -stats_all                                  >> $PRINSEQ/$temp.stats_all   2>&1" );
    system( "rm   $output_g/*_prinseq_* " );
}





say   "\n\n\n\n\n\n##################################################################################################";
say   "Detecting the quality of all FASTQ files by using kPAL, Rqc, ShortRead, seqTools and QC3 ......";
system( "kpal   count  -k 6   $output_g/*.fasta    $kPAL/rawReads.k6   >> $kPAL/rawReads.k6.runLog   2>&1");
system( "kpal   info                               $kPAL/rawReads.k6   >> $kPAL/rawReads.k6.runLog   2>&1");
system( "kpal   matrix   $kPAL/rawReads.k6   $kPAL/rawReads.k6.matrix  >> $kPAL/rawReads.k6.runLog   2>&1");
system( "rm   $output_g/*.fasta " );
system( "qc3.pl   -m f    -i $QC3/fileLists   -o $QC3   -t $numCores  >> $QC3/qc3.runLog   2>&1" );
system( "Rscript  0-Other/Rsrc/Rqc.R         $output_g    $R_QC        >> $R_QC/Rqc.runLog              2>&1" );
system( "Rscript  0-Other/Rsrc/ShortRead.R   $output_g    $ShortRead   >> $ShortRead/ShortRead.runLog   2>&1" );
system( "Rscript  0-Other/Rsrc/seqTools.R    $output_g    $seqTools    >> $seqTools/seqTools.runLog     2>&1" );





say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n";





## END
