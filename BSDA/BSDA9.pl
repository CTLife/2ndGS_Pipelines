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
my $input_g  = '';  ## such as "7-callMethylation"
my $output_g = '';  ## such as "10-coveredSites-Bismark"

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use BSDA (BS-Seq Data Analyzer), version 0.9.4,  2018-02-01.
        BSDA is a Pipeline for Single-end and Paired-end BS-Seq Data Analysis by Integrating Lots of Softwares.
                                                            
        Step 9: This script generates some figures for the correlation of number of covered C sites and reads.

                If this script works well, you do not need to check the the versions of the softwares or packages whcih are used in this script.
                And you do not need to exactly match the versions of the softwares or packages.
                If some errors or warnings are reported, please check the versions of softwares or packages.

                The versions of softwares or packages are used in this script:  
                        Perl, 5.22.1
                        R,    3.4.3

        Usage:
               perl  BSDA9.pl    [-version]    [-help]   [-genome RefGenome]    [-in inputDir]    [-out outDir]
        For instance:
               perl  BSDA9.pl   -genome hg38   -in 7-callMethylation   -out 10-coveredSites-Bismark    > BSDA9.runLog

        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.

        -help           Show this help message and exit.

        Required arguments:
        -genome RefGenome   "RefGenome" is the short name of your reference genome, such as "mm10", "ce11", "hg38".    (no default)

        -in inputDir        "inputDir" is the name of input path that contains your input files.  (no default)

        -out outDir         "outDir" is the name of output path that contains your running results of this step.  (no default)
        -----------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines, please visit https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ Jie Qiao Lab, yongp@outlook.com, Key Laboratory of Assisted Reproduction at Third Hospital,
        Academy for Advanced Interdisciplinary Studies, and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';

## Version Infromation
my $version = "    The 9th Step of BSDA (BS-Seq Data Analyzer), version 0.9.4,  2018-02-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$genome_g = 'hg38';                      ## This is only an initialization value or suggesting value, not default value.
$input_g  = '7-callMethylation';         ## This is only an initialization value or suggesting value, not default value.
$output_g = '10-coveredSites-Bismark';   ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help   -genome   -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  BSDA9.pl  -help' \n";
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
   
my $genomeFolder_g = "/media/yp/biox1/.RefGenomes/Shortcuts2/$genome_g";  
my $chromSize_g    = "$genomeFolder_g/$genome_g.chrom.sizes";
my $BAMpath_g      = "5-finalBAM/2A_Bismark";
&myMakeDir($output_g);

opendir(my $DH_input_g, $BAMpath_g)  ||  die;
my @inputFiles_g = readdir($DH_input_g);
my $pattern_g    = "[-.0-9A-Za-z]+";
my $numCores_g   = 8;
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
say   "\n\n\n\n\n\n##################################################################################################";
say   "Checking all the input file names ......";
my @groupFiles = ();
my $fileNameBool = 1;
for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {
        next unless $inputFiles_g[$i] =~ m/\.bam$/;
        next unless $inputFiles_g[$i] !~ m/^[.]/;
        next unless $inputFiles_g[$i] !~ m/[~]$/;
        next unless $inputFiles_g[$i] !~ m/^QC_Results$/;
        next unless $inputFiles_g[$i] !~ m/^unpaired/;
        say   "\t......$inputFiles_g[$i]" ;
        my $temp = $inputFiles_g[$i];
        $groupFiles[++$#groupFiles] = $inputFiles_g[$i];
        $temp =~ m/^(\d+)_($pattern_g)_(Rep[1-9])/   or  die   "wrong-1: ## $temp ##";
        $temp =~ m/_(Rep[1-9])\.bam$/   or  die   "wrong-2: ## $temp ##";
        if($temp !~ m/^((\d+)_($pattern_g)_(Rep[1-9]))(_[1-2])?\.bam$/) {
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
say   "Detecting BAM files in input folder ......";
my @BAMfiles_g = ();
{
open(seqFiles_FH, ">", "$output_g/BAM-Files.txt")  or  die; 
for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {     
    next unless $inputFiles_g[$i] =~ m/\.bam$/;
    next unless $inputFiles_g[$i] !~ m/^[.]/;
    next unless $inputFiles_g[$i] !~ m/[~]$/;
    next unless $inputFiles_g[$i] !~ m/^unpaired/;
    say    "\t......$inputFiles_g[$i]"; 
    $inputFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])\.bam$/  or  die;  
    $BAMfiles_g[$#BAMfiles_g+1] =  $inputFiles_g[$i];
    say        "\t\t\t\tBAM file:  $inputFiles_g[$i]\n";
    say   seqFiles_FH  "BAM file:  $inputFiles_g[$i]\n";
}

say   seqFiles_FH  "\n\n\n\n\n";  
say   seqFiles_FH  "All BAM files:@BAMfiles_g\n\n\n";
say        "\t\t\t\tAll BAM files:@BAMfiles_g\n\n";
my $num1 = $#BAMfiles_g + 1;
say seqFiles_FH   "\nThere are $num1 BAM files.\n";
say         "\t\t\t\tThere are $num1 BAM files.\n";
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Generating bedGraph and figures ......";

my $dir1 = "$output_g/1_CHG_cov1Reads";
my $dir2 = "$output_g/2_CHH_cov1Reads";
my $dir3 = "$output_g/3_CpG_cov1Reads";
&myMakeDir($dir1); 
&myMakeDir($dir2); 
&myMakeDir($dir3); 

my $dir1a = "$output_g/1_CHG_readsDis_mC";
my $dir2a = "$output_g/2_CHH_readsDis_mC";
my $dir3a = "$output_g/3_CpG_readsDis_mC";
&myMakeDir($dir1a); 
&myMakeDir($dir2a); 
&myMakeDir($dir3a); 

for (my $i=0; $i<=$#BAMfiles_g; $i++) {
    my $temp = $BAMfiles_g[$i]; 
    $temp =~ s/\.bam$//  ||  die; 
    say   "\t......$temp";

    my $SE1 = $temp;
    my $SE2 = $temp;
    $SE1 =~ s/_Rep/-end1_Rep/ or die;
    $SE2 =~ s/_Rep/-end2_Rep/ or die;
    my $PE_dir  = "$input_g/2A_Bismark/$temp";
    my $SE1_dir = "$input_g/2B_Bismark_unmapped_SE/$SE1";
    my $SE2_dir = "$input_g/2B_Bismark_unmapped_SE/$SE2";
   
    system("bismark2bedGraph   --output  $temp.CHG.bedGraph       --dir $dir1   --cutoff 1    --CX_context    --buffer_size  10G      $PE_dir/CHG_*       $SE1_dir/CHG_*         $SE2_dir/CHG_*        >  $dir1/$temp.runLog  2>&1 ");    
    system("bismark2bedGraph   --output  $temp.CHH.bedGraph       --dir $dir2   --cutoff 1    --CX_context    --buffer_size  10G      $PE_dir/CHH_*       $SE1_dir/CHH_*         $SE2_dir/CHH_*        >  $dir2/$temp.runLog  2>&1 ");    
    system("bismark2bedGraph   --output  $temp.CpG.bedGraph       --dir $dir3   --cutoff 1    --CX_context    --buffer_size  10G      $PE_dir/CpG_*       $SE1_dir/CpG_*         $SE2_dir/CpG_*        >  $dir3/$temp.runLog  2>&1 ");                          

    system("gunzip  --stdout    $dir1/$temp.CHG.bismark.cov.gz        >  $dir1/$temp.CHG.bismark.cov");
    system("gunzip  --stdout    $dir2/$temp.CHH.bismark.cov.gz        >  $dir2/$temp.CHH.bismark.cov");
    system("gunzip  --stdout    $dir3/$temp.CpG.bismark.cov.gz        >  $dir3/$temp.CpG.bismark.cov");
         
    &myMakeDir("$dir1a/$temp" ); 
    &myMakeDir("$dir2a/$temp" ); 
    &myMakeDir("$dir3a/$temp" ); 
 
    system("Rscript   BSDA9.R   $dir1   $temp.CHG.bismark.cov    $dir1a/$temp   ");         
    system("Rscript   BSDA9.R   $dir2   $temp.CHH.bismark.cov    $dir2a/$temp   ");         
    system("Rscript   BSDA9.R   $dir3   $temp.CpG.bismark.cov    $dir3a/$temp   ");  

    system("rm  -rf  $dir1/$temp.CHG.bismark.cov");     
    system("rm  -rf  $dir2/$temp.CHH.bismark.cov");     
    system("rm  -rf  $dir3/$temp.CpG.bismark.cov");     
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";



  

## END
