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
my $input_g  = '';  ## such as "7-finalBAM/1_Trim_BWAmem"
my $output_g = '';  ## such as "10-BigWig/1_Trim_BWAmem"

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use starCISDA (small-scale TELP-assisted rapid ChIP–seq Data Analyzer), version 0.9.0, 2017-10-01.
        starCISDA is a Pipeline for Single-end and Paired-end small-scale TELP-assisted rapid ChIP–seq Data Analysis by Integrating Lots of Softwares.
                                                            
        Step 9: Convert bam to bigwig by using bedtools and bedGraphToBigWig.
                Only reads density will be computed, no extend and no shift.
                So single-end and paired-end reads will not be distinguished.   

        Usage:
               perl  starCISDA9.pl    [-version]    [-help]   [-genome RefGenome]    [-in inputDir]    [-out outDir]
        For instance:
               perl  starCISDA9.pl   -genome hg38   -in 7-finalBAM/1_Trim_BWAmem   -out 10-BigWig/1_Trim_BWAmem    > starCISDA9.runLog

        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.

        -help           Show this help message and exit.

        Required arguments:
        -genome RefGenome   "RefGenome" is the short name of your reference genome, such as "mm10", "ce11", "hg38".    (no default)

        -in inputDir        "inputDir" is the name of input path that contains your BAM files.  (no default)

        -out outDir         "outDir" is the name of output path that contains your running results (bigwig files) of this step.  (no default)
        -----------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines, please visit https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ Jie Qiao Lab, yongp@outlook.com, Key Laboratory of Assisted Reproduction at Third Hospital,
        Academy for Advanced Interdisciplinary Studies, and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';

## Version Infromation
my $version = "    The Ninth Step of starCISDA (small-scale TELP-assisted rapid ChIP–seq Data Analyzer), version 0.9.0, 2017-10-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$genome_g = 'hg38';                       ## This is only an initialization value or suggesting value, not default value.
$input_g  = '7-finalBAM/1_Trim_BWAmem';   ## This is only an initialization value or suggesting value, not default value.
$output_g = '10-BigWig/1_Trim_BWAmem';    ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help   -genome   -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  starCISDA9.pl  -help' \n";
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

&printVersion("bedtools   --version");
&printVersion("bedGraphToBigWig");
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
say   "Detecting bam files in input folder ......";
sub numberOfLines  {
    my $filename = $_[0];
    say($filename);
    system("samtools view  $filename  |  wc -l   > yp_my_temp_1.xxxxx.txt");  
    open(tempFH, "<", "yp_my_temp_1.xxxxx.txt")  or  die;
    my @fullPath1 = <tempFH>; 
    ($#fullPath1 == 0)  or  die;
    system("rm   -rf   yp_my_temp_1.xxxxx.txt");
    $fullPath1[0] =~ m/^(\d+)\s+/  or  die;
    my $numLines = $1; 
    return($numLines);
}

my @bamfiles_g = ();
my @readsNum_g = ();
open(seqFiles_FH, ">", "$output2_g/bam-Files.txt")  or  die; 
for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {     
    next unless $inputFiles_g[$i] =~ m/\.bam$/;
    next unless $inputFiles_g[$i] !~ m/^[.]/;
    next unless $inputFiles_g[$i] !~ m/[~]$/;
    next unless $inputFiles_g[$i] !~ m/^unpaired/;
    next unless $inputFiles_g[$i] !~ m/^removed_/;
    say    "\t......$inputFiles_g[$i]"; 
    $inputFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])\.bam$/  or  die;  
    $bamfiles_g[$#bamfiles_g+1] =  $inputFiles_g[$i];
    $readsNum_g[$#readsNum_g+1] =  &numberOfLines("$input_g/$inputFiles_g[$i]");
    say        "\t\t\t\tbam file:  $inputFiles_g[$i]\n";
    say   seqFiles_FH  "bam file:  $inputFiles_g[$i]\n";
}

say   seqFiles_FH  "\n\n\n\n\n";  
say   seqFiles_FH  "All bam files:@bamfiles_g\n\n\n";
say        "\t\t\t\tAll bam files:@bamfiles_g\n\n";
my $num1 = $#bamfiles_g + 1;
say seqFiles_FH   "\nThere are $num1 bam files.\n";
say         "\t\t\t\tThere are $num1 bam files.\n";

###################################################################################################################################################################################################











###################################################################################################################################################################################################
{
say   "\n\n\n\n\n\n##################################################################################################";
say   "Convert bam to BigWig ......";
for (my $i=0; $i<=$#bamfiles_g; $i++) {
    my $temp = $bamfiles_g[$i]; 
    $temp =~ s/\.bam$//  ||  die; 
    say   "\t......$bamfiles_g[$i]";
    my $Temp_scale =  10**7/$readsNum_g[$i];
    print  "\t\t\t\t$bamfiles_g[$i]:\t$readsNum_g[$i]\t$Temp_scale\n";
    system("\t\tbedtools   genomecov  -bg   -split   -scale $Temp_scale   -ibam $input_g/$temp.bam   -g 0-Other/Shortcuts/$genome_g/$genome_g.chrom.sizes    > $output_g/$temp.bedGraph");
    sleep(3);
    system("bedGraphToBigWig    $output_g/$temp.bedGraph     0-Other/Shortcuts/$genome_g/$genome_g.chrom.sizes      $output_g/$temp.bw");
}
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";





## END
