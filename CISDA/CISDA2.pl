#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.18;
## Perl5 version >= 5.18,   you can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
## Help Infromation
my $HELP_g = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use CISDA (ChIP-Seq Data Analyzer), version 0.7.2, 2016-04-17.
        CISDA is a Pipeline for Single-end and Paired-end ChIP-Seq Data Analysis by Integrating Lots of Softwares.

        Step 2: Remove adaptors and PCR primers, trim and filter the raw reads by using Trimmomatic.
                And  assess the quality of the raw reads to identify possible sequencing errors or biases by using 12 softwares:
                FastQC, MultiQC, fastq-tools, FaQCs, prinseq, fastqp, QC3, NGS_QC_Toolkit, ht2-stat in HTQC, Rqc, ShortRead and seqTools.  

        Usage:
               perl  CISDA2.pl    [-version]    [-help]    [-in inputDir]    [-out outDir]
        For instance:
               perl  CISDA2.pl  -in 2-FASTQ   -out 3-Filtered    >> CISDA2.runLog  2>&1

        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.

        -help           Show this help message and exit.

        Required arguments:
        -input inputDir        "inputDir" is the name of input folder that contains your fastq files files.  (no default)

        -output outDir         "outDir" is the name of output folder that contains your running results (fastq format) of this step.  (no default)
        -----------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines such as RASDA, MESDA and HISDA,
        please see the web site: https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies
        and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';

## Version Infromation
my $version_g = "    The Second Step of CISDA (ChIP-Seq Data Analyzer), version 0.7.2, 2016-04-17.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP_g\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;    }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
my $input_g  = '2-FASTQ';           ## This is only an initialization value or suggesting value, not default value.
my $output_g = '3-Filtered';        ## This is only an initialization value or suggesting value, not default value.

## Available Arguments 
my $available_g = "   -version    -help    -in   -out  ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available_g !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole_g = 1; }
}
if($boole_g == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  CISDA2.pl  -help' \n";
    exit 0;
}

## Get Arguments 
if ( exists $args{'-version' }   )     { say  "\n$version_g\n";    exit 0; }
if ( exists $args{'-help'    }   )     { say  "\n$HELP_g\n";       exit 0; }
if ( exists $args{'-in'      }   )     { $input_g  = $args{'-in' };  }else{say   "\n -in  is required.\n";   say  "\n$HELP_g\n";    exit 0; }
if ( exists $args{'-out'     }   )     { $output_g = $args{'-out'};  }else{say   "\n -out is required.\n";   say  "\n$HELP_g\n";    exit 0; }

## Conditions
$input_g  =~ m/^\S+$/    ||  die   "\n\n$HELP_g\n\n";
$output_g =~ m/^\S+$/    ||  die   "\n\n$HELP_g\n\n";

## Print Command Arguments to Standard Output
say  "\n
        ################ Arguments ###############################
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
    if ( !( -e $path) )  { system("mkdir  -p  $path"); }
    if ( !( -e $path) )  { mkdir $path  ||  die; }
}
my $output2_g = "$output_g/QC_Results";
&myMakeDir("$output_g");
&myMakeDir("$output2_g");
opendir(my $DH_input_g, $input_g)  ||  die;
my @inputFiles_g = readdir($DH_input_g);
my $pattern_g    = "[-.0-9A-Za-z]+";
my $numCores_g   = 4;
###################################################################################################################################################################################################





###################################################################################################################################################################################################
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
my $Trimmomatic_g = "/home/yp/.MyProgramFiles/3_HTS2G/1_common/5_Trim/Trimmomatic-0.36/trimmomatic-0.36.jar";
&printVersion(" java  -jar  $Trimmomatic_g  -version");
&printVersion(" fastqc    -v ");
&printVersion(" multiqc           --version ");
&printVersion(" fastq-uniq        --version ");  ## in fastq-tools.
&printVersion(" FaQCs.pl          --version ");
&printVersion(" prinseq-lite.pl    -version ");
&printVersion(" fastqp    -h ");
&printVersion(" qc3.pl    -h ");
&printVersion(" IlluQC_PRLL.pl  -h ");
&printVersion(" ht2-stat  -v ");
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
say   "\n\n\n\n\n\n##################################################################################################";
say   "Checking all the input file names ......";
my @groupFiles = ();
my $fileNameBool = 1;
for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {
        next unless $inputFiles_g[$i] !~ m/^[.]/;
        next unless $inputFiles_g[$i] !~ m/[~]$/;
        next unless $inputFiles_g[$i] !~ m/^QC_Results$/;
        next unless $inputFiles_g[$i] =~ m/\.fastq$/;
        say   "\t......$inputFiles_g[$i]" ;
        my $temp = $inputFiles_g[$i];
        $groupFiles[++$#groupFiles] = $temp;
        $temp =~ m/^(\d+)_($pattern_g)_(Rep[1-9])_?([1-2]?)\.fastq$/    or  die   "wrong-1: ## $temp ##";
        if($temp !~ m/^((\d+)_($pattern_g)_(Rep[1-9]))(_[1-2])?(\.fastq)?/) {
             $fileNameBool = 0;
        }
}
if($fileNameBool == 1)  { say    "\n\t\tAll the file names are passed.\n";  }
@groupFiles   = sort(@groupFiles);
my $numGroup  = 0;
my $noteGroup = 0;
for ( my $i=0; $i<=$#groupFiles; $i++ ) {
    $groupFiles[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])/  or  die;
    my $n1 = $1;
    $n1>=1  or  die;
    if($noteGroup != $n1) {say "\n\t\tGroup $n1:";  $numGroup++; }
    say  "\t\t\t$groupFiles[$i]";
    $noteGroup = $n1;
}
say  "\n\t\tThere are $numGroup groups.";
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Detecting single-end and paired-end FASTQ files in input folder ......";     ## The fastq files are same between input folder and ouput folder.
my @singleEnd_g   = ();
my @pairedEnd_g   = ();
open(seqFiles_FH_g, ">", "$output2_g/singleEnd-pairedEnd-Files.txt")  or  die;
for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {
    next unless $inputFiles_g[$i] =~ m/\.fastq$/;
    next unless $inputFiles_g[$i] !~ m/^[.]/;
    next unless $inputFiles_g[$i] !~ m/[~]$/;
    next unless $inputFiles_g[$i] !~ m/^unpaired/;
    next unless $inputFiles_g[$i] !~ m/^QC_Results$/;
    say    "\t......$inputFiles_g[$i]";
    $inputFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])_?([1-2]?)\.fastq$/   or  die;
    if ($inputFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])\.fastq$/) {   ## sinlge end sequencing files.
        $inputFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])\.fastq$/  or  die;
        $singleEnd_g[$#singleEnd_g+1] =  $inputFiles_g[$i];
        say         "\t\t\t\tSingle-end sequencing files: $inputFiles_g[$i]\n";
        say  seqFiles_FH_g  "Single-end sequencing files: $inputFiles_g[$i]\n";
    }else{     ## paired end sequencing files.
        $inputFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])_([1-2])\.fastq$/  or  die;
        if ($inputFiles_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/) { ## The two files of one paired sequencing sample are always side by side.
            my $temp = $1;
            my $end1 = $temp."_1.fastq";
            my $end2 = $temp."_2.fastq";
            (-e  "$input_g/$end1")  or die;
            (-e  "$input_g/$end2")  or die;
            $pairedEnd_g[$#pairedEnd_g+1] =  $end1;
            $pairedEnd_g[$#pairedEnd_g+1] =  $end2;
            say        "\t\t\t\tPaired-end sequencing files: $end1,  $end2\n";
            say seqFiles_FH_g  "Paired-end sequencing files: $end1,  $end2\n";
        }
    }
}
( ($#pairedEnd_g+1)%2 == 0 )  or die;
say   seqFiles_FH_g  "\n\n\n\n\n";
say   seqFiles_FH_g  "All single-end sequencing files:@singleEnd_g\n\n\n";
say   seqFiles_FH_g  "All paired-end sequencing files:@pairedEnd_g\n\n\n";
say          "\t\t\t\tAll single-end sequencing files:@singleEnd_g\n\n";
say          "\t\t\t\tAll paired-end sequencing files:@pairedEnd_g\n\n";
my $numSingle_g = $#singleEnd_g + 1;
my $numPaired_g = $#pairedEnd_g + 1;
say seqFiles_FH_g   "\nThere are $numSingle_g single-end sequencing files.\n";
say seqFiles_FH_g   "\nThere are $numPaired_g paired-end sequencing files.\n";
say           "\t\t\t\tThere are $numSingle_g single-end sequencing files.\n";
say           "\t\t\t\tThere are $numPaired_g paired-end sequencing files.\n";
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
say    "\n\n\n\n\n\n##################################################################################################";
say    "Filtering the reads by using Trimmomatic ......"; 
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {
        $pairedEnd_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])_1\.fastq$/   or  die;
        $pairedEnd_g[$i] =~ m/^(\S+)_1.fastq$/ or die;
        my $temp = $1;   
        my $end1 = $temp."_1.fastq";
        my $end2 = $temp."_2.fastq";
        say   "\t......$end1";
        say   "\t......$end2";
        ($end2 eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$output2_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("java  -jar  $Trimmomatic_g  PE   -threads $numCores_g     $input_g/$end1  $input_g/$end2    $output_g/$end1  $output_g/unpaired-$end1    $output_g/$end2  $output_g/unpaired-$end2      ILLUMINACLIP:0-Other/Adapters/PE.fasta:1:30:10   LEADING:3   TRAILING:3   SLIDINGWINDOW:4:10   MINLEN:30   TOPHRED33    >>$output2_g/$temp.runLog  2>&1");                                                   
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {   
        say   "\t......$singleEnd_g[$i]" ;
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die; 
        my $temp = $1; 
        system("java  -jar  $Trimmomatic_g  SE   -threads $numCores_g     $input_g/$temp.fastq  $output_g/$temp.fastq    ILLUMINACLIP:0-Other/Adapters/SE.fasta:1:30:10   LEADING:3   TRAILING:3   SLIDINGWINDOW:4:10   MINLEN:30    TOPHRED33    >>$output2_g/$temp.runLog  2>&1");                    
}
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_FASTQ_1  {
    my $dir1      =  $_[0];   ## All the fastq files must be in this folder.
    my $QCresults = "$dir1/QC_Results";
    my $FastQC    = "$QCresults/1_FastQC";
    my $MultiQC   = "$QCresults/2_MultiQC";
    &myMakeDir($QCresults);
    &myMakeDir($FastQC);
    &myMakeDir($MultiQC);
    opendir(my $FH_Files, $dir1) || die;     
    my @Files = readdir($FH_Files);
    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of all FASTQ files by using FastQC and MultiQC ......";
    for ( my $i=0; $i<=$#Files; $i++ ) {
        next unless $Files[$i] =~ m/\.fastq$/;
        next unless $Files[$i] !~ m/^[.]/;
        next unless $Files[$i] !~ m/[~]$/;
        my $temp = $Files[$i];
        say    "\t......$temp";
        $temp =~ s/\.fastq$//  ||  die;
        system( "fastqc    --outdir $FastQC    --threads $numCores_g  --format fastq   --kmers 7    $dir1/$temp.fastq       >> $FastQC/$temp.runLog     2>&1" );
    }
    system( "multiqc  --verbose  --outdir $MultiQC          $FastQC/*_fastqc.zip      >> $MultiQC/MultiQC.FastQC.runLog   2>&1" );
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_FASTQ_2  {
    my $dir1      =  $_[0];   ## All the fastq files must be in this folder.
    my $QCresults = "$dir1/QC_Results";
    my $FastqTools= "$QCresults/3_FastqTools";
    &myMakeDir($QCresults);
    &myMakeDir($FastqTools);
    opendir(my $FH_Files, $dir1) || die;     
    my @Files = readdir($FH_Files);
    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of all FASTQ files by using fastq-tools ......";
    for ( my $i=0; $i<=$#Files; $i++ ) {
        next unless $Files[$i] =~ m/\.fastq$/;
        next unless $Files[$i] !~ m/^[.]/;
        next unless $Files[$i] !~ m/[~]$/;
        my $temp = $Files[$i];
        say    "\t......$temp";
        $temp =~ s/\.fastq$//  ||  die;
        &myMakeDir("$FastqTools/$temp");
        system( "fastq-kmers   -k 1     $dir1/$temp.fastq    >> $FastqTools/$temp/$temp.1mer     2>&1 " );
        system( "fastq-kmers   -k 2     $dir1/$temp.fastq    >> $FastqTools/$temp/$temp.2mer     2>&1 " );
        system( "fastq-kmers   -k 3     $dir1/$temp.fastq    >> $FastqTools/$temp/$temp.3mer     2>&1 " );
        system( "fastq-qual             $dir1/$temp.fastq    >> $FastqTools/$temp/$temp.qual     2>&1 " );
        system( "fastq-uniq  --verbose  $dir1/$temp.fastq    >> $FastqTools/$temp/$temp.uniq     2>&1 " );
    }
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_FASTQ_3  {
    my $dir1      =  $_[0];   ## All the fastq files must be in this folder.
    my $QCresults = "$dir1/QC_Results";
    my $Fastqp    = "$QCresults/4_Fastqp";
    my $FaQCs     = "$QCresults/5_FaQCs";
    &myMakeDir($QCresults);
    &myMakeDir($Fastqp);
    &myMakeDir($FaQCs);
    opendir(my $FH_Files, $dir1) || die;     
    my @Files = readdir($FH_Files);
    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of all FASTQ files by using fastqp, FaQCs ......";
    for ( my $i=0; $i<=$#Files; $i++ ) {
        next unless $Files[$i] =~ m/\.fastq$/;
        next unless $Files[$i] !~ m/^[.]/;
        next unless $Files[$i] !~ m/[~]$/;
        my $temp = $Files[$i];
        say    "\t......$temp";
        $temp =~ s/\.fastq$//  ||  die;
        system( "fastqp    --nreads 20000000   --kmer 4    --output $Fastqp/$temp  --type fastq   --median-qual 30     $dir1/$temp.fastq     >> $Fastqp/$temp.runLog     2>&1 " );
        system( "FaQCs.pl     -prefix $temp     -t $numCores_g      -qc_only  -min_L 30    -d $FaQCs/$temp          -u $dir1/$temp.fastq     >> $FaQCs/$temp.runLog          2>&1" );
    }
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_FASTQ_4  {
    my $dir1      =  $_[0];   ## All the fastq files must be in this folder.
    my $QCresults = "$dir1/QC_Results";
    my $PRINSEQ   = "$QCresults/6_PRINSEQ";
    &myMakeDir($QCresults);
    &myMakeDir($PRINSEQ);
    opendir(my $FH_Files, $dir1) || die;     
    my @Files = readdir($FH_Files);
    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of all FASTQ files by using prinseq ......";
    for ( my $i=0; $i<=$#Files; $i++ ) {
        next unless $Files[$i] =~ m/\.fastq$/;
        next unless $Files[$i] !~ m/^[.]/;
        next unless $Files[$i] !~ m/[~]$/;
        my $temp = $Files[$i];
        say    "\t......$temp";
        $temp =~ s/\.fastq$//  ||  die;
        &myMakeDir("$PRINSEQ/$temp");
        system( "prinseq-lite.pl  -out_format 1   -verbose    -fastq $dir1/$temp.fastq    -graph_data $PRINSEQ/$temp/$temp.gd      >> $PRINSEQ/$temp.runLog      2>&1" );
        system( "prinseq-graphs.pl   -i $PRINSEQ/$temp/$temp.gd   -png_all    -o $PRINSEQ/$temp/$temp                              >> $PRINSEQ/$temp.runLog      2>&1" );
        system( "prinseq-graphs.pl   -i $PRINSEQ/$temp/$temp.gd   -html_all   -o $PRINSEQ/$temp/$temp                              >> $PRINSEQ/$temp.runLog      2>&1" );
        system( "prinseq-lite.pl  -out_format 1  -verbose  -fastq  $dir1/$temp.fastq    -stats_all                                 >> $PRINSEQ/$temp.stats_all   2>&1" );
        system( "rm   $dir1/*_prinseq_* " );
    }
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
&myQC_FASTQ_1($output_g);
&myQC_FASTQ_2($output_g);
&myQC_FASTQ_3($output_g);
&myQC_FASTQ_4($output_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
say   "\n\n\n\n\n\n##################################################################################################";
say   "Detecting the quality of paired-end FASTQ files by using NGS_QC_Toolkit, FaQCs, prinseq, ht2-stat in HTQC ......";
my $NGSQC_Toolkit = "$output2_g/7_NGSQC_Toolkit_PairedEnd";
my $FaQCs         = "$output2_g/8_FaQCs_PairedEnd";
my $HTQC          = "$output2_g/9_HTQC_PairedEnd";
my $PRINSEQ       = "$output2_g/10_PRINSEQ_PairedEnd"; 
&myMakeDir($NGSQC_Toolkit); 
&myMakeDir($FaQCs); 
&myMakeDir($HTQC); 
&myMakeDir($PRINSEQ); 
for ( my $j=0; $j<=$#pairedEnd_g; $j=$j+2 ) {
    my $temp1 = $pairedEnd_g[$j];
    my $temp2 = $pairedEnd_g[$j+1];
    my $temp  = $temp2;
    say   "\t......$temp1";
    say   "\t......$temp2";
    $temp1 =~ m/^(\d+)_($pattern_g)_(Rep[1-9])_1\.fastq$/   or  die;
    $temp2 =~ m/^(\d+)_($pattern_g)_(Rep[1-9])_2\.fastq$/   or  die;
    $temp1 =~ s/\.fastq$//    ||  die;
    $temp2 =~ s/\.fastq$//    ||  die;
    $temp  =~ s/_2\.fastq$//  ||  die;
    &myMakeDir("$NGSQC_Toolkit/$temp"); 
    &myMakeDir("$FaQCs/$temp"); 
    &myMakeDir("$HTQC/$temp"); 
    &myMakeDir("$PRINSEQ/$temp"); 
    system( "IlluQC_PRLL.pl      -pe $output_g/$temp1.fastq  $output_g/$temp2.fastq   N    A      -cpus $numCores_g     -onlyStat          -outputFolder $NGSQC_Toolkit/$temp           >> $NGSQC_Toolkit/$temp.runLog   2>&1" );
    system( "FaQCs.pl     -prefix $temp     -t $numCores_g      -qc_only  -min_L 30    -d $FaQCs/$temp          -p $output_g/$temp1.fastq  $output_g/$temp2.fastq            >> $FaQCs/$temp.runLog          2>&1" );
    system( "ht2-stat  -i $output_g/$temp1.fastq  $output_g/$temp2.fastq   -o  $HTQC/$temp    --pe     --encode sanger   --threads $numCores_g  >> $HTQC/$temp.runLog   2>&1" );
    system( "prinseq-lite.pl   -out_format 1   -verbose  -fastq $output_g/$temp1.fastq  -fastq2 $output_g/$temp2.fastq   -graph_data $PRINSEQ/$temp/$temp.gd      >> $PRINSEQ/$temp.runLog   2>&1" );
    system( "prinseq-graphs.pl   -i $PRINSEQ/$temp/$temp.gd   -png_all    -o $PRINSEQ/$temp/$temp                 >> $PRINSEQ/$temp.runLog   2>&1" );
    system( "prinseq-graphs.pl   -i $PRINSEQ/$temp/$temp.gd   -html_all   -o $PRINSEQ/$temp/$temp                 >> $PRINSEQ/$temp.runLog   2>&1" );
    system( "prinseq-lite.pl   -out_format 1  -verbose  -fastq $output_g/$temp.fastq    -stats_all                                  >> $PRINSEQ/$temp.stats_all   2>&1" );
    system( "rm   $output_g/*_prinseq_* " );
}
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
say   "\n\n\n\n\n\n##################################################################################################";
say   "Detecting the quality of all FASTQ files by using QC3, Rqc, ShortRead and seqTools ......";
my $QC3        = "$output2_g/11_QC3";
my $R_QC       = "$output2_g/12_R_QC";
my $ShortRead  = "$output2_g/13_ShortRead";
my $seqTools   = "$output2_g/14_seqTools";  
&myMakeDir($QC3); 
&myMakeDir($R_QC); 
&myMakeDir($ShortRead); 
&myMakeDir($seqTools);
my @outputFiles_g = @inputFiles_g;
for ( my $i=0; $i<=$#outputFiles_g; $i++ ) {    
    next unless $outputFiles_g[$i] =~ m/\.fastq$/;
    next unless $outputFiles_g[$i] !~ m/^[.]/;
    next unless $outputFiles_g[$i] !~ m/[~]$/;
    next unless $outputFiles_g[$i] !~ m/^unpaired/;
    next unless $outputFiles_g[$i] !~ m/^QC_Results$/;
    system( "echo         '$output_g/$outputFiles_g[$i]'   >> $QC3/fileLists  ");
}
system( "qc3.pl   -m f    -i $QC3/fileLists   -o $QC3   -t $numCores_g   >> $QC3/qc3.runLog               2>&1" );
system( "Rscript  0-Other/R_SRC/Rqc.R         $output_g    $R_QC         >> $R_QC/Rqc.runLog              2>&1" );
system( "Rscript  0-Other/R_SRC/ShortRead.R   $output_g    $ShortRead    >> $ShortRead/ShortRead.runLog   2>&1" );
system( "Rscript  0-Other/R_SRC/seqTools.R    $output_g    $seqTools     >> $seqTools/seqTools.runLog     2>&1" );
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n";





## END
