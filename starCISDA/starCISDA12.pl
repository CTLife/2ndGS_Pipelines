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
my $input_g  = '';  ## such as "12A-MACS2/1_Trim_BWAmem/1-Peaks/run10"
my $output_g = '';  ## such as "13A-IDR/1_Trim_BWAmem/1-Peaks/run10"                

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use starCISDA (small-scale TELP-assisted rapid ChIP–seq Data Analyzer), version 0.9.0, 2017-10-01.
        starCISDA is a Pipeline for Single-end and Paired-end small-scale TELP-assisted rapid ChIP–seq Data Analysis by Integrating Lots of Softwares.

        Step 12: Caculate Irreproducible Discovery Rate (IDR) based on the ouptuts of MACS2 (narrowPeak or broadPeak format).
                 https://github.com/nboley/idr

        Usage:
               perl  starCISDA12.pl    [-version]    [-help]   [-genome RefGenome]    [-in inputDir]    [-out outDir]
        For instance:
               perl  starCISDA12.pl   -genome hg38   -in 12A-MACS2/1_Trim_BWAmem/1-Peaks/run10   -out 13A-IDR/1_Trim_BWAmem/1-Peaks/run10    > starCISDA12.runLog

        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.

        -help           Show this help message and exit.

        Required arguments:
        -genome RefGenome   "RefGenome" is the short name of your reference genome, such as "mm10", "ce11", "hg38".    (no default)

        -in inputDir        "inputDir" is the name of input path that contains your BAM files.  (no default)

        -out outDir         "outDir" is the name of output path that contains your running results of this step.  (no default)
        -----------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines, please visit https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ Jie Qiao Lab, yongp@outlook.com, Key Laboratory of Assisted Reproduction at Third Hospital,
        Academy for Advanced Interdisciplinary Studies, and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';

## Version Infromation
my $version = "    The Twelfth Step of starCISDA (small-scale TELP-assisted rapid ChIP–seq Data Analyzer), version 0.9.0, 2017-10-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$genome_g = 'hg38';                                    ## This is only an initialization value or suggesting value, not default value.
$input_g  = '12A-MACS2/1_Trim_BWAmem/1-Peaks/run10';   ## This is only an initialization value or suggesting value, not default value.
$output_g = '13A-IDR/1_Trim_BWAmem/1-Peaks/run10';     ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help   -genome   -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  starCISDA12.pl  -help' \n";
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
say    "\n\n\n\n\n\n##################################################################################################";
say    "Running......";

sub myMakeDir  {
    my $path = $_[0];
    if ( !( -e $path) )  { system("mkdir  -p  $path"); }
    if ( !( -e $path) )  { mkdir $path  ||  die; }
}

my $output2_g = "$output_g/QC_Results";
&myMakeDir($output_g);
&myMakeDir($output2_g);

opendir(my $DH_input_g, $input_g)  ||  die;
my @inputFiles_g = readdir($DH_input_g);
my $pattern_g    = "[-.0-9A-Za-z]+";
my $numCores_g   = 4;
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Checking all the necessary softwares in this step......" ;

sub printVersion  {
    my $software = $_[0];
    system("echo    '##############################################################################'  >> $output2_g/VersionsOfSoftwares.txt   2>&1");
    system("echo    '#########$software'                                                              >> $output2_g/VersionsOfSoftwares.txt   2>&1");
    system("$software                                                                                 >> $output2_g/VersionsOfSoftwares.txt   2>&1");
    system("echo    '\n\n\n\n\n\n'                                                                    >> $output2_g/VersionsOfSoftwares.txt   2>&1");
}

sub fullPathApp  {
    my $software = $_[0];
    say($software);
    system("which   $software  > yp_my_temp_1.$software.txt");
    open(tempFH, "<", "yp_my_temp_1.$software.txt")  or  die;
    my @fullPath1 = <tempFH>;
    ($#fullPath1 == 0)  or  die;
    system("rm  yp_my_temp_1.$software.txt");
    $fullPath1[0] =~ s/\n$//  or  die;
    return($fullPath1[0]);
}

&printVersion("idr  -h");
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Caculating Irreproducible Discovery Rate (IDR) ......";

{

my $temp = "1_A-KD17031177-15-R1-PCGF1";
my $rep1 = $temp . "_Rep1/$temp" . "_Rep1_peaks.broadPeak";
my $rep2 = "1_A-KD17031177-16-R1-PCGF1_Rep2/1_A-KD17031177-16-R1-PCGF1_Rep2_peaks.broadPeak";
&myMakeDir("$output_g/$temp");
system("idr  --samples  $input_g/$rep1    $input_g/$rep2     --input-file-type broadPeak    --output-file  $output_g/$temp/idrValues.txt   --soft-idr-threshold  0.05      --plot       > $output_g/$temp/$temp.runLog    2>&1");


$temp = "1_20161206WY-3-3-Ctrl-F1-H3K27me3";
$rep1 = $temp . "_Rep1/$temp" . "_Rep1_peaks.broadPeak";
$rep2 = "1_20161206WY-3-4-KO-F9-H3K27me3_Rep1/1_20161206WY-3-4-KO-F9-H3K27me3_Rep1_peaks.broadPeak";
&myMakeDir("$output_g/$temp");
system("idr  --samples  $input_g/$rep1    $input_g/$rep2    --input-file-type broadPeak    --output-file  $output_g/$temp/idrValues.txt   --soft-idr-threshold  0.05      --plot       > $output_g/$temp/$temp.runLog    2>&1");


}
###################################################################################################################################################################################################
## system("idr  --samples  $input_g/$rep1    $input_g/$rep2    --input-file-type broadPeak    --output-file  $output_g/$temp/idrValues.txt   --soft-idr-threshold  0.05      --plot       > $output_g/$temp/$temp.runLog    2>&1");





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";





## END
