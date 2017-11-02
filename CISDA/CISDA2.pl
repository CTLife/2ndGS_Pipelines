#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.22;

## Perl5 version >= 5.22
## You can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
## This step is not neccesary.
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $genome_g = '';  ## such as "hg38", "ce11", "mm10".    
my $input_g  = '';  ## such as "2-mergedFASTQ"   
my $output_g = '';  ## such as "3-rmDupFASTQ"

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use CISDA (ChIP-Seq Data Analyzer), version 0.9.0, 2017-10-01.
        CISDA is a Pipeline for Single-end and Paired-end ChIP-Seq Data Analysis by Integrating Lots of Softwares.

        Step 2: Discard optical and PCR duplicates by using ParDRe and Clumpify in BBMap.
                And  assess the quality of the raw reads to identify possible sequencing errors or biases by using 12 softwares:
                FastQC, fastq-tools, bbcountunique.sh in BBMap, AfterQC, FaQCs, fastqp, FastQ_Screen, PRINSEQ, QC3, seqTools, ShortRead and Rqc. 
                And aggregate the results from FastQC or FastQ_Screen analyses across many samples into a single report by using MultiQC.         

        Usage:
               perl  CISDA2.pl    [-version]    [-help]   [-genome RefGenome]    [-in inputDir]    [-out outDir]
        For instance:
               perl  CISDA2.pl   -genome hg38   -in 2-mergedFASTQ   -out 3-rmDupFASTQ    > CISDA2.runLog  

        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.

        -help           Show this help message and exit.

        Required arguments:
        -genome RefGenome   "RefGenome" is the short name of your reference genome, such as "hg38", "ce11", "mm10".    (no default)

        -in inputDir        "inputDir" is the name of input path that contains your FASTQ files.  (no default)

        -out outDir         "outDir" is the name of output path that contains your running results (fastq files) of this step.  (no default)
        -----------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines, please visit https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ Jie Qiao Lab, yongp@outlook.com, Key Laboratory of Assisted Reproduction at Third Hospital,  
        Academy for Advanced Interdisciplinary Studies, and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';

## Version Infromation
my $version = "    The Second Step of CISDA (ChIP-Seq Data Analyzer), version 0.9.0, 2017-10-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$genome_g = 'hg38';           ## This is only an initialization value or suggesting value, not default value.
$input_g  = '2-mergedFASTQ';  ## This is only an initialization value or suggesting value, not default value.
$output_g = '3-rmDupFASTQ';   ## This is only an initialization value or suggesting value, not default value.

## Available Arguments 
my $available = "   -version    -help   -genome   -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  CISDA2.pl  -help' \n";
    exit 0;
}

## Get Arguments 
if ( exists $args{'-version' }   )     { say  "\n$version\n";    exit 0; }
if ( exists $args{'-help'    }   )     { say  "\n$HELP\n";       exit 0; }
if ( exists $args{'-genome'  }   )     { $genome_g = $args{'-genome'  }; }else{say   "\n -genome is required.\n";   say  "\n$HELP\n";    exit 0; }
if ( exists $args{'-in'      }   )     { $input_g  = $args{'-in'      }; }else{say   "\n -in     is required.\n";   say  "\n$HELP\n";    exit 0; }
if ( exists $args{'-out'     }   )     { $output_g = $args{'-out'     }; }else{say   "\n -out    is required.\n";   say  "\n$HELP\n";    exit 0; }

## Conditions
$genome_g =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";
$input_g  =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";
$output_g =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";

## Print Command Arguments to Standard Output
say  "\n
        ################ Arguments ###############################
                Reference Genome:  $genome_g
                Input       Path:  $input_g
                Output      Path:  $output_g
        ###############################################################
\n";
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Running ......";

## These files must be available:
my $Rqc_g        = "0-Other/R_SRC/Rqc.R";
my $seqTools_g   = "0-Other/R_SRC/seqTools.R";
my $ShortRead_g  = "0-Other/R_SRC/ShortRead.R";

sub myMakeDir  {
    my $path = $_[0];
    if ( !( -e $path) )  { system("mkdir  -p  $path"); }
    if ( !( -e $path) )  { mkdir $path  ||  die;       }
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
{
## These commands must be available:
say   "\n\n\n\n\n\n##################################################################################################";
say   "Checking all the necessary softwares in this step ......";

sub printVersion  {
    my $software = $_[0];
    system("echo    '##############################################################################'  >> $output2_g/VersionsOfSoftwares.txt   2>&1");
    system("echo    '#########$software : '                                                           >> $output2_g/VersionsOfSoftwares.txt   2>&1");
    system("$software                                                                                 >> $output2_g/VersionsOfSoftwares.txt   2>&1");
    system("echo    '\n\n\n\n\n\n'                                                                    >> $output2_g/VersionsOfSoftwares.txt   2>&1");
}

&printVersion(" ParDRe  -h ");
&printVersion(" clumpify.sh  -h ");

&printVersion(" fastqc        --version ");
&printVersion(" multiqc       --version ");
&printVersion(" fastq-uniq    --version ");
&printVersion(" fastq_screen  --version ");
&printVersion(" after.py      --version ");
&printVersion(" FaQCs.pl         -h ");  
&printVersion(" fastqp           -h ");
&printVersion(" qc3.pl           -h ");
&printVersion(" prinseq-lite.pl  -h ");
&printVersion(" bbcountunique.sh -h ");

}
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
say   seqFiles_FH_g  "All single-end sequencing files: @singleEnd_g\n\n\n";
say   seqFiles_FH_g  "All paired-end sequencing files: @pairedEnd_g\n\n\n";
say          "\t\t\t\tAll single-end sequencing files: @singleEnd_g\n\n";
say          "\t\t\t\tAll paired-end sequencing files: @pairedEnd_g\n\n";
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
say    "Remove duplicates with clumpify.sh in BBMap and ParDRe ......"; 
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
        system("clumpify.sh  in=$input_g/$end1  in2=$input_g/$end2    out=$output_g/clumpify.$end1  out2=$output_g/clumpify.$end2     dedupe=t   optical=f  dupedist=40  -Xmx20g  >> $output2_g/$temp.clumpify.runLog  2>&1");                                                   
        system("ParDRe  -i $output_g/clumpify.$end1  -p $output_g/clumpify.$end2    -o  $output_g/$end1  -r $output_g/$end2     -t $numCores_g       >> $output2_g/$temp.ParDRe.runLog  2>&1");                                                   
        system( "rm   $output_g/clumpify.* " );
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {   
        say   "\t......$singleEnd_g[$i]" ;
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die; 
        my $temp = $1; 
        system("clumpify.sh  in=$input_g/$temp.fastq    out=$output_g/clumpify.$temp.fastq      dedupe=t   optical=f   dupedist=40  -Xmx20g >> $output2_g/$temp.clumpify.runLog  2>&1");  
        system("ParDRe  -i $output_g/clumpify.$temp.fastq    -o  $output_g/$temp.fastq   -t $numCores_g      >> $output2_g/$temp.ParDRe.runLog  2>&1");                                                   
        system( "rm   $output_g/clumpify.* " );                  
}
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_FASTQ_1  {
    my $dir1      =  $_[0];   ## All the fastq files must be in this folder.
    my $QCresults = "$dir1/QC_Results";
    my $FastQC    = "$QCresults/1_FastQC";
    my $MultiQC   = "$QCresults/MultiQC/FastQC";
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
    system( "multiqc  --title FastQC   --verbose  --export  --outdir $MultiQC          $FastQC      >> $MultiQC/MultiQC.FastQC.runLog   2>&1" );
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_FASTQ_2  {
    my $dir1      =  $_[0];   ## All the fastq files must be in this folder.
    my $QCresults = "$dir1/QC_Results";
    my $FastQ_Screen   = "$QCresults/2_FastQ_Screen";
    my $MultiQC        = "$QCresults/MultiQC/FastQ_Screen";
    &myMakeDir($QCresults);
    &myMakeDir($FastQ_Screen);
    &myMakeDir($MultiQC);
    opendir(my $FH_Files, $dir1) || die;     
    my @Files = readdir($FH_Files);
    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of all FASTQ files by using FastQ_Screen and MultiQC ......";
    for ( my $i=0; $i<=$#Files; $i++ ) {
        next unless $Files[$i] =~ m/\.fastq$/;
        next unless $Files[$i] !~ m/^[.]/;
        next unless $Files[$i] !~ m/[~]$/;
        my $temp = $Files[$i];
        say    "\t......$temp";
        $temp =~ s/\.fastq$//  ||  die;
        system( "fastq_screen  --aligner  bowtie2  --outdir $FastQ_Screen/$temp  --threads $numCores_g    --subset 1000000   $dir1/$temp.fastq      >> $FastQ_Screen/$temp.runLog      2>&1" );
    }
    system( "multiqc  --title FastQ_Screen   --verbose  --export  --outdir $MultiQC     $FastQ_Screen     >> $MultiQC/MultiQC.FastQ_Screen.runLog   2>&1" );
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_FASTQ_3  {
    my $dir1      =  $_[0];   ## All the fastq files must be in this folder.
    my $QCresults = "$dir1/QC_Results";
    my $PRINSEQ   = "$QCresults/3_PRINSEQ";
    my $BBUnique  = "$QCresults/4_BBCountUnique";
    &myMakeDir($QCresults);
    &myMakeDir($PRINSEQ);
    &myMakeDir($BBUnique);
    opendir(my $FH_Files, $dir1) || die;     
    my @Files = readdir($FH_Files);
    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of all FASTQ files by using prinseq and bbcountunique.sh ......";
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
        system( "rm   $dir1/*_prinseq_* " );
        system( "bbcountunique.sh  in=$dir1/$temp.fastq       out=$BBUnique/$temp.txt   samplerate=0.8  -Xmx20g   interval=100000    >> $BBUnique/$temp.runLog      2>&1" );
    } 
}

###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_FASTQ_4  {
    my $dir1        =  $_[0];   ## All the fastq files must be in this folder.
    my $QCresults   = "$dir1/QC_Results";
    my $FastqTools  = "$QCresults/5_fastq_tools";
    my $FaQCs       = "$QCresults/6_FaQCs";
    &myMakeDir($QCresults);
    &myMakeDir($FastqTools);
    &myMakeDir($FaQCs);
    opendir(my $FH_Files, $dir1) || die;     
    my @Files = readdir($FH_Files);
    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of all FASTQ files by using fastq_tools and FaQCs ......";
    for ( my $i=0; $i<=$#Files; $i++ ) {
        next unless $Files[$i] =~ m/\.fastq$/;
        next unless $Files[$i] !~ m/^[.]/;
        next unless $Files[$i] !~ m/[~]$/;
        my $temp = $Files[$i];
        say    "\t......$temp";
        $temp =~ s/\.fastq$//  ||  die;
        &myMakeDir("$FaQCs/$temp");
        system( "fastq-uniq        $dir1/$temp.fastq          >  $FastqTools/$temp.uniq2 ");
        system( "head  -n 100000   $FastqTools/$temp.uniq2    >  $FastqTools/$temp.uniq " );
        system( "rm  $FastqTools/$temp.uniq2" );
        system( "FaQCs.pl   -qc_only   -u $dir1/$temp.fastq    -prefix $temp  -t $numCores_g    -subset 1  -debug 1   -d $FaQCs/$temp   >> $FaQCs/$temp.runLog   2>&1" );
    }
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_FASTQ_5  {
    my $dir1      =  $_[0];   ## All the fastq files must be in this folder.
    my $QCresults = "$dir1/QC_Results";
    my $Fastqp    = "$QCresults/7_Fastqp";
    my $AfterQC   = "$QCresults/8_AfterQC";
    &myMakeDir($QCresults);
    &myMakeDir($Fastqp);
    &myMakeDir($AfterQC);
    opendir(my $FH_Files, $dir1) || die;     
    my @Files = readdir($FH_Files);
    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of all FASTQ files by using fastqp and AfterQC ......";
    for ( my $i=0; $i<=$#Files; $i++ ) {
        next unless $Files[$i] =~ m/\.fastq$/;
        next unless $Files[$i] !~ m/^[.]/;
        next unless $Files[$i] !~ m/[~]$/;
        my $temp = $Files[$i];
        say    "\t......$temp";
        $temp =~ s/\.fastq$//  ||  die;
        system( "fastqp    --nreads 200000   --kmer 4    --output $Fastqp/$temp  --type fastq   --median-qual 30   --count-duplicates   $dir1/$temp.fastq     >> $Fastqp/$temp.runLog     2>&1 " );
        system( "after.py  -1 $dir1/$temp.fastq  --qc_only  --qc_sample=200000  -g $AfterQC/$temp -b $AfterQC/$temp  -r $AfterQC/$temp  >> $AfterQC/$temp.runLog     2>&1 " );
    }
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
&myQC_FASTQ_1($output_g);
&myQC_FASTQ_2($output_g);
&myQC_FASTQ_3($output_g);
&myQC_FASTQ_4($output_g);
&myQC_FASTQ_5($output_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
say   "\n\n\n\n\n\n##################################################################################################";
say   "Detecting the quality of paired-end FASTQ files by using FaQCs, bbcountunique.sh, prinseq and AfterQC ......";

my $FaQCs    = "$output2_g/9_FaQCs_PairedEnd";    
my $bbunique = "$output2_g/10_bbcountunique_PairedEnd";    
my $prinseq  = "$output2_g/11_prinseq_PairedEnd";  
my $AfterQC  = "$output2_g/12_AfterQC_PairedEnd";  
&myMakeDir($FaQCs); 
&myMakeDir($bbunique); 
&myMakeDir($prinseq); 
&myMakeDir($AfterQC); 

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
    &myMakeDir("$FaQCs/$temp"); 
    &myMakeDir("$prinseq/$temp"); 
    system( "FaQCs.pl   -qc_only   -p $output_g/$temp1.fastq  $output_g/$temp2.fastq      -subset 1    -prefix $temp  -t $numCores_g   -debug 1   -d $FaQCs/$temp   >> $FaQCs/$temp.runLog    2>&1" );
    system( "bbcountunique.sh  in=$output_g/$temp1.fastq   in2=$output_g/$temp2.fastq     out=$bbunique/$temp.txt    samplerate=0.8  -Xmx20g   interval=100000      >> $bbunique/$temp.runLog 2>&1" );
    system( "prinseq-lite.pl  -out_format 1   -verbose    -fastq $output_g/$temp1.fastq   -fastq2  $output_g/$temp2.fastq  -graph_data $prinseq/$temp/$temp.gd      >> $prinseq/$temp.runLog  2>&1" );
    system( "prinseq-graphs.pl   -i $prinseq/$temp/$temp.gd   -png_all    -o $prinseq/$temp/$temp                              >> $prinseq/$temp.runLog      2>&1" );
    system( "prinseq-graphs.pl   -i $prinseq/$temp/$temp.gd   -html_all   -o $prinseq/$temp/$temp                              >> $prinseq/$temp.runLog      2>&1" );
    system( "rm   $output_g/*_prinseq_* " );
    system( "after.py  -1 $output_g/$temp1.fastq  -2 $output_g/$temp2.fastq   --qc_only  --qc_sample=200000  -g $AfterQC/$temp -b $AfterQC/$temp  -r $AfterQC/$temp >> $AfterQC/$temp.runLog  2>&1" );
}

}
###################################################################################################################################################################################################
 




say   "\n\n\n\n\n\n##################################################################################################";
say   "\n\n\n\n\n \t Now, you can stop this program!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n\n\n ";





###################################################################################################################################################################################################
{
say   "\n\n\n\n\n\n##################################################################################################";
say   "Detecting the quality of all FASTQ files by using QC3, Rqc, ShortRead and seqTools ......";

my $R_QC       = "$output2_g/13_R_QC";
my $ShortRead  = "$output2_g/14_ShortRead";
my $seqTools   = "$output2_g/15_seqTools"; 
my $QC3        = "$output2_g/16_QC3";
 
&myMakeDir($R_QC); 
&myMakeDir($ShortRead); 
&myMakeDir($seqTools);
&myMakeDir($QC3); 

system( "Rscript  $Rqc_g         $output_g    $R_QC         >> $R_QC/Rqc.runLog              2>&1" );
system( "Rscript  $ShortRead_g   $output_g    $ShortRead    >> $ShortRead/ShortRead.runLog   2>&1" );
system( "Rscript  $seqTools_g    $output_g    $seqTools     >> $seqTools/seqTools.runLog     2>&1" );

my @outputFiles_g = @inputFiles_g;
for ( my $i=0; $i<=$#outputFiles_g; $i++ ) {
    next unless $outputFiles_g[$i] =~ m/\.fastq$/;
    next unless $outputFiles_g[$i] !~ m/^[.]/;
    next unless $outputFiles_g[$i] !~ m/[~]$/;
    next unless $outputFiles_g[$i] !~ m/^unpaired/;
    next unless $outputFiles_g[$i] !~ m/^QC_Results$/;
    system( "echo         '$output_g/$outputFiles_g[$i]'   >> $QC3/fileLists  ");
}

system( "qc3.pl   -m f    -i $QC3/fileLists   -o $QC3   -t $numCores_g   >> $QC3/qc3.runLog  2>&1" );

}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";





## END
