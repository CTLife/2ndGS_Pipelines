#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.22;

## Perl5 version >= 5.22
## You can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $genome_g = '';  ## such as "mm10", "ce11", "hg38".    
my $input_g  = '';  ## such as "1-rawReads"   
my $output_g = '';  ## such as "2-mergedFASTQ"

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use RASDA (RNA-Seq Data Analyzer), version 0.9.0, 2017-10-01.
        RASDA is a Pipeline for Single-end and Paired-end RNA-Seq Data Analysis by Integrating Lots of Softwares.

        Step 1: Extract all compressed FASTQ files, or convert SRA to FASTQ format by using SRA_Toolkit, 
                and merge the two lanes of the same sample (merge technical replicates).

                And  assess the quality of the raw reads to identify possible sequencing errors or biases by using 12 softwares:
                FastQC, fastq-tools, bbcountunique.sh in BBMap, AfterQC, FaQCs, fastqp, FastQ_Screen, PRINSEQ, QC3, seqTools, ShortRead and Rqc. 
                And aggregate the results from FastQC or FastQ_Screen analyses across many samples into a single report by using MultiQC.         

        Usage:
               perl  RASDA1.pl    [-version]    [-help]   [-genome RefGenome]    [-in inputDir]    [-out outDir]
        For instance:
               perl  RASDA1.pl   -genome hg38   -in 1-rawReads   -out 2-mergedFASTQ    > RASDA1.runLog  

        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.

        -help           Show this help message and exit.

        Required arguments:
        -genome RefGenome   "RefGenome" is the short name of your reference genome, such as "mm10", "ce11", "hg38".    (no default)

        -in inputDir        "inputDir" is the name of input path that contains your compressed FASTQ files or SRA files.  (no default)

        -out outDir         "outDir" is the name of output path that contains your running results (fastq files) of this step.  (no default)
        -----------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines, please visit https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ Jie Qiao Lab, yongp@outlook.com, Key Laboratory of Assisted Reproduction at Third Hospital,  
        Academy for Advanced Interdisciplinary Studies, and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';

## Version Infromation
my $version = "    The First Step of RASDA (RNA-Seq Data Analyzer), version 0.9.0, 2017-10-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$genome_g = 'hg38';           ## This is only an initialization value or suggesting value, not default value.
$input_g  = '1-rawReads';     ## This is only an initialization value or suggesting value, not default value.
$output_g = '2-mergedFASTQ';  ## This is only an initialization value or suggesting value, not default value.

## Available Arguments 
my $available = "   -version    -help   -genome   -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  RASDA1.pl  -help' \n";
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

my $input2_g  = "$input_g/QC_Results";
my $output2_g = "$output_g/QC_Results";
&myMakeDir("$input2_g");
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

&printVersion(" fastq-dump  --version ");
&printVersion(" bzip2       --version ");
&printVersion(" gunzip      --version ");
&printVersion(" tar         --version ");
&printVersion(" xz          --version ");
&printVersion(" unzip       -hh       ");
&printVersion(" unrar   ");

&printVersion(" fastqc           --version ");
&printVersion(" multiqc          --version ");
&printVersion(" fastq-uniq       --version ");
&printVersion(" fastq_screen     --version ");
&printVersion(" after.py         --version ");
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
        say   "\t......$inputFiles_g[$i]" ;
        my $temp = $inputFiles_g[$i];
        $groupFiles[++$#groupFiles] = $temp;
        $temp =~ m/^(\d+)_($pattern_g)_(Rep[1-9])/    or  die   "wrong-1: ## $temp ##";
        $temp =~ m/_(Rep[1-9])\.sra$/  or  $temp =~ m/_(Rep[1-9])_?([1-2]?)(_Lane[1-2])?\.fastq/  or  die   "wrong-2: ## $temp ##";
        if($temp !~ m/^((\d+)_($pattern_g)_(Rep[1-9]))(_[1-2])?(_Lane[1-2])?(\.fastq)?/) {
             $fileNameBool = 0;
        }
}
if($fileNameBool == 1)  { say    "\n\t\tAll the file names are passed.\n";  }

@groupFiles   = sort(@groupFiles);
my $numGroup  = 0;
my $noteGroup = 0;
my $noteGroup2= 0;
for ( my $i=0; $i<=$#groupFiles; $i++ ) {
    $groupFiles[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])/  or  die;
    my $n1 = $1;
    $n1>=1  or  die;
    if($noteGroup2 != $n1) {say "\n\t\tGroup $n1:";  $numGroup++; }
    say  "\t\t\t$groupFiles[$i]";
    $noteGroup2 = $n1;  
    if($noteGroup < $n1) {$noteGroup = $n1;}
}
($noteGroup == $numGroup)  or  die "##($noteGroup == $numGroup)##\n\n";  ## so the first number of FASTQ file name must be from small to large, 1,2,3,4,5,6,7,8,9 ......   
say  "\n\t\tThere are $numGroup groups.";
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
say   "\n\n\n\n\n\n##################################################################################################";
say   "Converting SRA files into FASTQ files or extracting the compressed fastq files ......";
for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {
        next unless $inputFiles_g[$i] !~ m/^[.]/;
        next unless $inputFiles_g[$i] !~ m/[~]$/;
        next unless $inputFiles_g[$i] !~ m/^QC_Results$/;
        say  "\t......$inputFiles_g[$i]" ;
        my $temp = $inputFiles_g[$i];
        if ($temp =~ m/\.sra$/) {
                $temp =~ m/^(\d+)_($pattern_g)_(Rep[1-9])\.sra$/   or  die;
                system("fastq-dump   --split-3   --dumpbase   --outdir $input_g    $input_g/$temp   >> $input2_g/$temp.runLog  2>&1");  ## Formats sequence using base space.
        }elsif($temp =~ m/^((\d+)_($pattern_g)_(Rep[1-9])_?([1-2]?)(_Lane[1-2])?\.fastq)\.(\S+)$/) {
                $temp =~ m/^((\d+)_($pattern_g)_(Rep[1-9])_?([1-2]?)(_Lane[1-2])?\.fastq)\.(\S+)$/   or  die;
                my  $tempFastq = $1;
                my  $suffix1   = $7;    ## Only Nine compressed formats are supported, their suffixes:  ".bz2",  ".gz",  ".tar.gz",  ".tar",  ".rar",  ".xz",  ".zip", ".tar.bz2", ".tar.xz".
                my  $tempBool  = 0;
                if($suffix1 eq "bz2"    )  { $tempBool++;  system("bzip2   -cd         $input_g/$temp   >  $input_g/$tempFastq");  }
                if($suffix1 eq "gz"     )  { $tempBool++;  system("gunzip  --stdout    $input_g/$temp   >  $input_g/$tempFastq");  }
                if($suffix1 eq "tar.gz" )  { $tempBool++;  system("tar     -xzvf       $input_g/$temp  -C  $input_g");             }
                if($suffix1 eq "tar"    )  { $tempBool++;  system("tar     -xvf        $input_g/$temp  -C  $input_g");             }
                if($suffix1 eq "rar"    )  { $tempBool++;  system("unrar    e          $input_g/$temp      $input_g");             }
                if($suffix1 eq "xz"     )  { $tempBool++;  system("xz      -cd         $input_g/$temp   >  $input_g/$tempFastq");  }
                if($suffix1 eq "zip"    )  { $tempBool++;  system("unzip   -np         $input_g/$temp   >  $input_g/$tempFastq");  }
                if($suffix1 eq "tar.bz2")  { $tempBool++;  system("tar     -jxvf       $input_g/$temp  -C  $input_g");             }
                if($suffix1 eq "tar.xz" )  { $tempBool++;  system("tar     -xvJf       $input_g/$temp  -C  $input_g");             }
                if($tempBool  != 1)        { say("\n$temp is wrong!!\n");  die; }
        }else{
                $temp =~ m/^(\d+)_($pattern_g)_(Rep[1-9])_?([1-2]?)(_Lane[1-2])?\.fastq$/   or  die;
        }
}
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
say   "\n\n\n\n\n\n##################################################################################################";
say   "Merge the two lanes of the same sample or copy files ......";
opendir(my $DH_input, $input_g)  ||  die;
my @twoLanes = readdir($DH_input);
for (my $i=0; $i<=$#twoLanes; $i++) {
    next unless $twoLanes[$i] !~ m/^QC_Results$/;
    next unless $twoLanes[$i] =~ m/\.fastq$/;
    next unless $twoLanes[$i] !~ m/^[.]/;
    next unless $twoLanes[$i] !~ m/[~]$/;
    say    "\t......$twoLanes[$i]";
    $twoLanes[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])_?([1-2]?)(_Lane[1-2])?\.fastq$/   or  die;
    if( $twoLanes[$i] =~ m/^(\S+)_Lane1.fastq$/ ) {
        my $temp = $1;
        $temp =~ m/^(\d+)_($pattern_g)_(Rep[1-9])_?([1-2]?)$/   or  die;
        open(tempFH, ">>", "$output2_g/Merge-Two-Lanes.log")  or  die;
        my $lane1 = $temp."_Lane1.fastq";
        my $lane2 = $temp."_Lane2.fastq";
        system(       "cat  $input_g/$lane1  $input_g/$lane2   > $output_g/$temp.fastq" );
        say  tempFH   "cat  $input_g/$lane1  $input_g/$lane2   > $output_g/$temp.fastq";
    }
    if ($twoLanes[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])_?([1-2]?).fastq$/) {
        system(  "cp   $input_g/$twoLanes[$i]    $output_g/$twoLanes[$i]" );
    }
}
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Detecting single-end and paired-end FASTQ files ......";

opendir(my $DH_output_g, $output_g) || die;
my @outputFiles_g = readdir($DH_output_g);
my @singleEnd_g   = ();
my @pairedEnd_g   = ();
open(seqFiles_FH_g, ">", "$output2_g/singleEnd-pairedEnd-Files.txt")  or  die;

for ( my $i=0; $i<=$#outputFiles_g; $i++ ) {
    next unless $outputFiles_g[$i] =~ m/\.fastq$/;
    next unless $outputFiles_g[$i] !~ m/^[.]/;
    next unless $outputFiles_g[$i] !~ m/[~]$/;
    next unless $outputFiles_g[$i] !~ m/^unpaired/;
    next unless $outputFiles_g[$i] !~ m/^QC_Results$/;
    say    "\t......$outputFiles_g[$i]";
    $outputFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])_?([1-2]?)\.fastq$/   or  die;
    if ($outputFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])\.fastq$/) {   ## sinlge end sequencing files.
        $outputFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])\.fastq$/  or  die;
        $singleEnd_g[$#singleEnd_g+1] =  $outputFiles_g[$i];
        say         "\t\t\t\tSingle-end sequencing files: $outputFiles_g[$i]\n";
        say  seqFiles_FH_g  "Single-end sequencing files: $outputFiles_g[$i]\n";
    }else{     ## paired end sequencing files.
        $outputFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])_([1-2])\.fastq$/  or  die;
        if ($outputFiles_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/) { ## The two files of one paired sequencing sample are always side by side.
            my $temp = $1;
            my $end1 = $temp."_1.fastq";
            my $end2 = $temp."_2.fastq";
            (-e  "$output_g/$end1")  or die;
            (-e  "$output_g/$end2")  or die;
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
&myQC_FASTQ_1($input_g);
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
