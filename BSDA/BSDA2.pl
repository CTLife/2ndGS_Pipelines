#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  File::Copy;  
use  v5.22;
## Perl5 version >= 5.22
## You can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "2-rawFASTQ"   
my $output_g = '';  ## such as "3-finalFASTQ"

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use BSDA (BS-Seq Data Analyzer), version 0.9.4,  2018-02-01.
        BSDA is a Pipeline for Single-end and Paired-end BS-Seq Data Analysis by Integrating Lots of Softwares.

        Step 2: Remove adaptors and PCR primers, trim and filter the reads by using TrimGalore and fastp.

                And  assess the quality of the raw reads to identify possible sequencing errors or biases by using 8 softwares:
                FastQC, fastp, bbcountunique.sh in BBMap, fastq-tools, FastQ_Screen, fastqp and PRINSEQ. 
                And aggregate the results from FastQC and FastQ_Screen analyses across many samples into a single report 
                by using MultiQC.
     
                If this script works well, you do not need to check the the versions of the softwares or packages whcih are used in this script. 
                And you do not need to exactly match the versions of the softwares or packages.
                If some errors or warnings are reported, please check the versions of softwares or packages.

                The versions of softwares or packages are used in this script:  
                        Perl,   5.22.1 
                        cutadapt, 1.15
                        TrimGalore, 0.4.5
                        FastQC, 0.11.6
                        fastp,  0.12.2
                        bbcountunique.sh in BBMap, 37.77
                        fastq-tools,   0.8
                        FastQ_Screen,  0.11.4 
                        fastqp,        0.3.2
                        PRINSEQ,       0.20.4   
                        MultiQC,       1.3        

        Usage:
               perl  BSDA2.pl    [-version]    [-help]    [-in inputDir]    [-out outDir]
        For instance:
               perl  BSDA2.pl   -in 2-rawFASTQ   -out 3-finalFASTQ    > BSDA2.runLog  

        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.

        -help           Show this help message and exit.

        Required arguments:
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
my $version = "    The 2nd Step of BSDA (BS-Seq Data Analyzer), version 0.9.4,  2018-02-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$input_g  = '2-rawFASTQ';      ## This is only an initialization value or suggesting value, not default value.
$output_g = '3-finalFASTQ';    ## This is only an initialization value or suggesting value, not default value.

## Available Arguments 
my $available = "   -version    -help   -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  BSDA2.pl  -help' \n";
    exit 0;
}

## Get Arguments 
if ( exists $args{'-version' }   )     { say  "\n$version\n";    exit 0; }
if ( exists $args{'-help'    }   )     { say  "\n$HELP\n";       exit 0; }
if ( exists $args{'-in'      }   )     { $input_g  = $args{'-in'      }; }else{say   "\n -in     is required.\n";   say  "\n$HELP\n";    exit 0; }
if ( exists $args{'-out'     }   )     { $output_g = $args{'-out'     }; }else{say   "\n -out    is required.\n";   say  "\n$HELP\n";    exit 0; }

## Conditions
$input_g  =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";
$output_g =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";

## Print Command Arguments to Standard Output
say  "\n
        ################ Arguments ###############################
                Input   Path:  $input_g
                Output  Path:  $output_g
        ###############################################################
\n";
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Running ......";

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
my $numCores_g   = 8;
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

&printVersion(" cutadapt      --version ");
&printVersion(" trim_galore   --version ");

&printVersion(" fastqc           --version ");
&printVersion(" fastp            --help ");
&printVersion(" bbcountunique.sh -h ");
&printVersion(" bbmap.sh         --version ");
&printVersion(" fastq-uniq       --version ");
&printVersion(" fastq_screen     --version ");
&printVersion(" fastqp           -h ");
&printVersion(" prinseq-lite.pl  -version ");
&printVersion(" multiqc          --version ");
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
say    "Filtering the reads by using trim_galore ......"; 
&myMakeDir("$output2_g/trim_galore");          

for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {
        $pairedEnd_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])_1\.fastq$/   or  die;
        $pairedEnd_g[$i] =~ m/^(\S+)_1.fastq$/ or die;
        my $temp = $1;   
        my $end1 = $temp."_1.fastq";
        my $end2 = $temp."_2.fastq";
        say   "\t......$end1";
        say   "\t......$end2";
        ($end2 eq $pairedEnd_g[$i+1])  or  die;
        system("trim_galore  --quality 20   --phred33  --length 20    --rrbs     --output_dir  $output2_g/trim_galore    --paired  $input_g/$end1  $input_g/$end2          >> $output2_g/trim_galore/$temp.runLog  2>&1"); 
        my $end1a = $temp."_1_val_1.fq";
        my $end2a = $temp."_2_val_2.fq";        
        move("$output2_g/trim_galore/$end1a",  "$output_g/$end1");
        move("$output2_g/trim_galore/$end2a",  "$output_g/$end2");
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {   
        say   "\t......$singleEnd_g[$i]" ;
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die; 
        my $temp = $1; 
        system("trim_galore  --quality 20   --phred33  --length 20     --rrbs    --output_dir  $output2_g/trim_galore       $input_g/$temp.fastq      >> $output2_g/trim_galore/$temp.runLog  2>&1");   
        my $end1a = $temp."_trimmed.fq";
        move("$output2_g/trim_galore/$end1a",  "$output_g/$temp.fastq");
}
 

say    "\n\n\n\n\n\n##################################################################################################";
say    "Filtering the reads by using fastp ......"; 
&myMakeDir("$output_g/fastp");          
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {
        $pairedEnd_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])_1\.fastq$/   or  die;
        $pairedEnd_g[$i] =~ m/^(\S+)_1.fastq$/ or die;
        my $temp = $1;   
        my $end1 = $temp."_1.fastq";
        my $end2 = $temp."_2.fastq";
        say   "\t......$end1";
        say   "\t......$end2";
        ($end2 eq $pairedEnd_g[$i+1])  or  die;
        &myMakeDir("$output_g/fastp/$temp");          
        system( "fastp    --in1 $output_g/$end1  --in2 $output_g/$end2  --out1 $output_g/fastp/$end1  --out2  $output_g/fastp/$end2   -p  --report_title $temp  --thread $numCores_g     --json $output_g/fastp/$temp/$temp.json   --html $output_g/fastp/$temp/$temp.html   >> $output_g/fastp/$temp/$temp.runLog   2>&1" );                                           
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {   
        say   "\t......$singleEnd_g[$i]" ;
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die; 
        my $temp = $1; 
        &myMakeDir("$output_g/fastp/$temp");          
        system( "fastp    --in1 $output_g/$temp.fastq  --out1 $output_g/fastp/$temp.fastq   -p  --report_title $temp  --thread $numCores_g     --json $output_g/fastp/$temp/$temp.json   --html $output_g/fastp/$temp/$temp.html   >> $output_g/fastp/$temp/$temp.runLog   2>&1" );                                           
}

}

print("\n\n\n\n\n#########################################\n");
print("Now, you can start the next step!\n\n\n\n");
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_FASTQ_1  {
    my $dir1      =  $_[0];   ## All the fastq files must be in this folder.
    my $QCresults = "$dir1/QC_Results";
    my $FastQC    = "$QCresults/1_FastQC";
    my $fastp     = "$QCresults/2_fastp";
    my $MultiQC1  = "$QCresults/MultiQC/FastQC";
    my $MultiQC2  = "$QCresults/MultiQC/fastp";
    &myMakeDir($QCresults);
    &myMakeDir($FastQC);
    &myMakeDir($fastp);
    &myMakeDir($MultiQC1);
    &myMakeDir($MultiQC2);  
    opendir(my $FH_Files, $dir1) || die;       
    my @Files = readdir($FH_Files);
    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of all FASTQ files by using FastQC, fastp and MultiQC ......";
    for ( my $i=0; $i<=$#Files; $i++ ) {
        next unless $Files[$i] =~ m/\.fastq$/;
        next unless $Files[$i] !~ m/^[.]/;
        next unless $Files[$i] !~ m/[~]$/;
        my $temp = $Files[$i];
        say    "\t......$temp";
        $temp =~ s/\.fastq$//  ||  die;
        system( "fastqc    --outdir $FastQC    --threads $numCores_g  --format fastq   --kmers 7    $dir1/$temp.fastq       >> $FastQC/$temp.runLog     2>&1" );
        &myMakeDir("$fastp/$temp");        
        system( "fastp    --in1 $dir1/$temp.fastq  -p  --report_title $temp  --thread $numCores_g     --json $fastp/$temp/$temp.json   --html $fastp/$temp/$temp.html   >> $fastp/$temp/$temp.runLog   2>&1" );                                           
    }
    system( "multiqc  --title FastQC   --filename  FastQC  --module  fastqc  --verbose  --export  --outdir $MultiQC1    $FastQC  >> $MultiQC1/MultiQC.FastQC.runLog   2>&1" );
    system( "multiqc  --title fastp    --filename  fastp                     --verbose  --export  --outdir $MultiQC2    $fastp   >> $MultiQC2/MultiQC.fastp.runLog    2>&1" );
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_FASTQ_2  {
    my $dir1      =  $_[0];   ## All the fastq files must be in this folder.
    my $QCresults = "$dir1/QC_Results";
    my $FastQ_Screen   = "$QCresults/3_FastQ_Screen";
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
        system( "fastq_screen  --bisulfite    --aligner  bowtie2  --outdir $FastQ_Screen/$temp  --threads $numCores_g    --subset 100000   $dir1/$temp.fastq      >> $FastQ_Screen/$temp.runLog      2>&1" );
    }
    system( "multiqc  --title FastQ_Screen  --filename  FastQ_Screen  --module  fastq_screen   --verbose  --export  --outdir $MultiQC     $FastQ_Screen     >> $MultiQC/MultiQC.FastQ_Screen.runLog   2>&1" );
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_FASTQ_3  {
    my $dir1        =  $_[0];   ## All the fastq files must be in this folder.
    my $QCresults   = "$dir1/QC_Results";
    my $FastqTools  = "$QCresults/4_fastq_tools";
    &myMakeDir($QCresults);
    &myMakeDir($FastqTools);
    opendir(my $FH_Files, $dir1) || die;     
    my @Files = readdir($FH_Files);
    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of all FASTQ files by using fastq_tools and fastp ......";
    for ( my $i=0; $i<=$#Files; $i++ ) {
        next unless $Files[$i] =~ m/\.fastq$/;
        next unless $Files[$i] !~ m/^[.]/;
        next unless $Files[$i] !~ m/[~]$/;
        my $temp = $Files[$i];
        say    "\t......$temp";
        $temp =~ s/\.fastq$//  ||  die;
        system( "fastq-uniq        $dir1/$temp.fastq          >  $FastqTools/$temp.uniq2 ");
        system( "head  -n 100000   $FastqTools/$temp.uniq2    >  $FastqTools/$temp.uniq " );
        system( "rm  $FastqTools/$temp.uniq2" );
    }
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_FASTQ_4  {
    my $dir1      =  $_[0];   ## All the fastq files must be in this folder.
    my $QCresults = "$dir1/QC_Results";
    my $Fastqp    = "$QCresults/5_Fastqp";
    &myMakeDir($QCresults);
    &myMakeDir($Fastqp);
    opendir(my $FH_Files, $dir1) || die;     
    my @Files = readdir($FH_Files);
    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of all FASTQ files by using fastqp ......";
    for ( my $i=0; $i<=$#Files; $i++ ) {
        next unless $Files[$i] =~ m/\.fastq$/;
        next unless $Files[$i] !~ m/^[.]/;
        next unless $Files[$i] !~ m/[~]$/;
        my $temp = $Files[$i];
        say    "\t......$temp";
        $temp =~ s/\.fastq$//  ||  die;
        system( "fastqp    --nreads 200000   --kmer 4    --output $Fastqp/$temp  --type fastq   --median-qual 30   --count-duplicates   $dir1/$temp.fastq     >> $Fastqp/$temp.runLog     2>&1 " );
    }
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_FASTQ_5  {
    my $dir1      =  $_[0];   ## All the fastq files must be in this folder.
    my $QCresults = "$dir1/QC_Results";
    my $PRINSEQ   = "$QCresults/6_PRINSEQ";
    my $BBUnique  = "$QCresults/7_BBCountUnique";
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
&myQC_FASTQ_1($output_g);
&myQC_FASTQ_1("$output_g/fastp");
&myQC_FASTQ_2($output_g);
&myQC_FASTQ_3($output_g);
&myQC_FASTQ_4($output_g);
&myQC_FASTQ_5($output_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
say   "\n\n\n\n\n\n##################################################################################################";
say   "Detecting the quality of paired-end FASTQ files by using bbcountunique.sh and fastp ......";

my $bbunique = "$output2_g/8_bbcountunique_PairedEnd"; 
my $fastp = "$output2_g/9_fastp_PairedEnd";    
&myMakeDir($bbunique); 
&myMakeDir($fastp); 

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
    &myMakeDir("$fastp/$temp");    
    system( "bbcountunique.sh  in=$output_g/$temp1.fastq   in2=$output_g/$temp2.fastq     out=$bbunique/$temp.txt    samplerate=0.8  -Xmx20g   interval=100000      >> $bbunique/$temp.runLog 2>&1" );
    system( "fastp          --in1 $output_g/$temp1.fastq  --in2 $output_g/$temp2.fastq  -p  --report_title $temp  --thread $numCores_g     --json $fastp/$temp/$temp.json   --html $fastp/$temp/$temp.html   >> $fastp/$temp/$temp.runLog   2>&1" );                                           
}

}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";





## END
