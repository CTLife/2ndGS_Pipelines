#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.22;

## Perl5 version >= 5.22
## You can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
## Lots of memories and time are required.
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $genome_g = '';  ## such as "mm10", "ce11", "hg38".
my $input_g  = '';  ## such as "6-BAMPE/1_Bismark"
my $output_g = '';  ## such as "BSseq-Saturation-1/1_Bismark"

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use BSDA (BS-Seq Data Analyzer), version 0.9.0, 2017-10-01.
        BSDA is a Pipeline for Single-end and Paired-end BS-Seq Data Analysis by Integrating Lots of Softwares.
                                                            
        BSseq Saturation Analysis for Mapped Reads.

        Usage:
               perl  BSseq-Saturation-1.pl    [-version]    [-help]     [-in inputDir]    [-out outDir]
        For instance:
               perl  BSseq-Saturation-1.pl   -in 6-BAMPE/1_Bismark   -out BSseq-Saturation-1/1_Bismark    > BSseq-Saturation-1.runLog

        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.

        -help           Show this help message and exit.

        Required arguments:
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
my $version = "    The Seventh Step of BSDA (BS-Seq Data Analyzer), version 0.9.0, 2017-10-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$input_g  = '6-BAMPE/1_Bismark';              ## This is only an initialization value or suggesting value, not default value.
$output_g = 'BSseq-Saturation-1/1_Bismark';      ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help    -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  BSseq-Saturation-1.pl  -help' \n";
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

&myMakeDir($output_g);
my $pattern_g    = "[-.0-9A-Za-z]+";
opendir(my $DH_input_g, $input_g)  ||  die;
my @inputFiles_g = readdir($DH_input_g);

my $numCores_g   = 4;
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

for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {     
    next unless $inputFiles_g[$i] =~ m/\.bam$/;
    next unless $inputFiles_g[$i] !~ m/^[.]/;
    next unless $inputFiles_g[$i] !~ m/[~]$/;
    next unless $inputFiles_g[$i] !~ m/^unpaired/;
    say    "\t......$inputFiles_g[$i]"; 
    $inputFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])\.bam$/  or  die;  
    $BAMfiles_g[$#BAMfiles_g+1] =  $inputFiles_g[$i];
    say        "\t\t\t\tBAM file:  $inputFiles_g[$i]\n";
}

say        "\t\t\t\tAll BAM files:@BAMfiles_g\n\n";
my $num1 = $#BAMfiles_g + 1;
say         "\t\t\t\tThere are $num1 BAM files.\n";
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Extract mC ......";
for (my $i=0; $i<=$#BAMfiles_g; $i++) {
    my $temp = $BAMfiles_g[$i]; 
    $temp =~ s/\.bam$//  ||  die; 
    say   "\t......$BAMfiles_g[$i]";
    &myMakeDir("$output_g/$temp");  
  
    my $numReads_1 = readpipe("samtools   view  $input_g/$temp.bam   |  wc -l");
    $numReads_1 =~ s/\s//g;
    my $numReads_2 = $numReads_1 * 0.8;
    my $numReads_3 = $numReads_1 * 0.6;
    my $numReads_4 = $numReads_1 * 0.4;
    my $numReads_5 = $numReads_1 * 0.2;
    my $numReads_6 = $numReads_1 * 0.01;
    say("\nNumber of reads of $temp : #$numReads_1#\n");
    
    system("bedtools  sample  -n $numReads_1  -i $input_g/$temp.bam             > $output_g/$temp/100.$temp.bam");
    system("bedtools  sample  -n $numReads_2  -i $output_g/$temp/100.$temp.bam  > $output_g/$temp/80.$temp.bam");
    system("bedtools  sample  -n $numReads_3  -i $output_g/$temp/80.$temp.bam   > $output_g/$temp/60.$temp.bam");
    system("bedtools  sample  -n $numReads_4  -i $output_g/$temp/60.$temp.bam   > $output_g/$temp/40.$temp.bam");
    system("bedtools  sample  -n $numReads_5  -i $output_g/$temp/40.$temp.bam   > $output_g/$temp/20.$temp.bam");
    system("bedtools  sample  -n $numReads_6  -i $output_g/$temp/20.$temp.bam   > $output_g/$temp/1.$temp.bam");

    &myMakeDir("$output_g/$temp/100.$temp");
    &myMakeDir("$output_g/$temp/80.$temp");
    &myMakeDir("$output_g/$temp/60.$temp");
    &myMakeDir("$output_g/$temp/40.$temp");
    &myMakeDir("$output_g/$temp/20.$temp");
    &myMakeDir("$output_g/$temp/1.$temp");

    system("bismark_methylation_extractor  --single-end   --output $output_g/$temp/100.$temp  --mbias_off  --gzip    --comprehensive  --buffer_size 16G     $output_g/$temp/100.$temp.bam   > $output_g/$temp/100.$temp.runLog   2>&1");                                       
    system("bismark_methylation_extractor  --single-end   --output $output_g/$temp/80.$temp   --mbias_off  --gzip    --comprehensive  --buffer_size 16G     $output_g/$temp/80.$temp.bam    > $output_g/$temp/80.$temp.runLog    2>&1");                                       
    system("bismark_methylation_extractor  --single-end   --output $output_g/$temp/60.$temp   --mbias_off  --gzip    --comprehensive  --buffer_size 16G     $output_g/$temp/60.$temp.bam    > $output_g/$temp/60.$temp.runLog    2>&1");                                       
    system("bismark_methylation_extractor  --single-end   --output $output_g/$temp/40.$temp   --mbias_off  --gzip    --comprehensive  --buffer_size 16G     $output_g/$temp/40.$temp.bam    > $output_g/$temp/40.$temp.runLog    2>&1");                                       
    system("bismark_methylation_extractor  --single-end   --output $output_g/$temp/20.$temp   --mbias_off  --gzip    --comprehensive  --buffer_size 16G     $output_g/$temp/20.$temp.bam    > $output_g/$temp/20.$temp.runLog    2>&1");                                       
    system("bismark_methylation_extractor  --single-end   --output $output_g/$temp/1.$temp    --mbias_off  --gzip    --comprehensive  --buffer_size 16G     $output_g/$temp/1.$temp.bam     > $output_g/$temp/1.$temp.runLog     2>&1");                                       

    system("rm -rf $output_g/$temp/*.bam");

    my $file1 = "CHG_context_100.".$temp;
    my $file2 = "CHH_context_100.".$temp;
    my $file3 = "CpG_context_100.".$temp;
    system("bismark2bedGraph  --buffer_size 16G   --CX_context  --dir $output_g/$temp/100.$temp   -o $file1.bedGraph    $output_g/$temp/100.$temp/$file1.txt.gz");
    system("bismark2bedGraph  --buffer_size 16G   --CX_context  --dir $output_g/$temp/100.$temp   -o $file2.bedGraph    $output_g/$temp/100.$temp/$file2.txt.gz");
    system("bismark2bedGraph  --buffer_size 16G                 --dir $output_g/$temp/100.$temp   -o $file3.bedGraph    $output_g/$temp/100.$temp/$file3.txt.gz");
    system("gunzip  --stdout    $output_g/$temp/100.$temp/$file1.bismark.cov.gz   >  $output_g/$temp/100.$temp/$file1.bismark.cov  ");  
    system("gunzip  --stdout    $output_g/$temp/100.$temp/$file2.bismark.cov.gz   >  $output_g/$temp/100.$temp/$file2.bismark.cov  ");  
    system("gunzip  --stdout    $output_g/$temp/100.$temp/$file3.bismark.cov.gz   >  $output_g/$temp/100.$temp/$file3.bismark.cov  ");  
    
    $file1 = "CHG_context_80.".$temp;
    $file2 = "CHH_context_80.".$temp;
    $file3 = "CpG_context_80.".$temp;
    system("bismark2bedGraph  --buffer_size 16G   --CX_context  --dir $output_g/$temp/80.$temp   -o  $file1.bedGraph    $output_g/$temp/80.$temp/$file1.txt.gz");
    system("bismark2bedGraph  --buffer_size 16G   --CX_context  --dir $output_g/$temp/80.$temp   -o  $file2.bedGraph    $output_g/$temp/80.$temp/$file2.txt.gz");
    system("bismark2bedGraph  --buffer_size 16G                 --dir $output_g/$temp/80.$temp   -o  $file3.bedGraph    $output_g/$temp/80.$temp/$file3.txt.gz");
    system("gunzip  --stdout    $output_g/$temp/80.$temp/$file1.bismark.cov.gz   >  $output_g/$temp/80.$temp/$file1.bismark.cov  ");  
    system("gunzip  --stdout    $output_g/$temp/80.$temp/$file2.bismark.cov.gz   >  $output_g/$temp/80.$temp/$file2.bismark.cov  ");  
    system("gunzip  --stdout    $output_g/$temp/80.$temp/$file3.bismark.cov.gz   >  $output_g/$temp/80.$temp/$file3.bismark.cov  ");  
    
    $file1 = "CHG_context_60.".$temp;
    $file2 = "CHH_context_60.".$temp;
    $file3 = "CpG_context_60.".$temp;
    system("bismark2bedGraph  --buffer_size 16G   --CX_context  --dir $output_g/$temp/60.$temp   -o  $file1.bedGraph    $output_g/$temp/60.$temp/$file1.txt.gz");
    system("bismark2bedGraph  --buffer_size 16G   --CX_context  --dir $output_g/$temp/60.$temp   -o  $file2.bedGraph    $output_g/$temp/60.$temp/$file2.txt.gz");
    system("bismark2bedGraph  --buffer_size 16G                 --dir $output_g/$temp/60.$temp   -o  $file3.bedGraph    $output_g/$temp/60.$temp/$file3.txt.gz");
    system("gunzip  --stdout    $output_g/$temp/60.$temp/$file1.bismark.cov.gz   >  $output_g/$temp/60.$temp/$file1.bismark.cov  ");  
    system("gunzip  --stdout    $output_g/$temp/60.$temp/$file2.bismark.cov.gz   >  $output_g/$temp/60.$temp/$file2.bismark.cov  ");  
    system("gunzip  --stdout    $output_g/$temp/60.$temp/$file3.bismark.cov.gz   >  $output_g/$temp/60.$temp/$file3.bismark.cov  ");  
    
    $file1 = "CHG_context_40.".$temp;
    $file2 = "CHH_context_40.".$temp;
    $file3 = "CpG_context_40.".$temp;
    system("bismark2bedGraph  --buffer_size 16G   --CX_context  --dir $output_g/$temp/40.$temp   -o  $file1.bedGraph    $output_g/$temp/40.$temp/$file1.txt.gz");
    system("bismark2bedGraph  --buffer_size 16G   --CX_context  --dir $output_g/$temp/40.$temp   -o  $file2.bedGraph    $output_g/$temp/40.$temp/$file2.txt.gz");
    system("bismark2bedGraph  --buffer_size 16G                 --dir $output_g/$temp/40.$temp   -o  $file3.bedGraph    $output_g/$temp/40.$temp/$file3.txt.gz");
    system("gunzip  --stdout    $output_g/$temp/40.$temp/$file1.bismark.cov.gz   >  $output_g/$temp/40.$temp/$file1.bismark.cov  ");  
    system("gunzip  --stdout    $output_g/$temp/40.$temp/$file2.bismark.cov.gz   >  $output_g/$temp/40.$temp/$file2.bismark.cov  ");  
    system("gunzip  --stdout    $output_g/$temp/40.$temp/$file3.bismark.cov.gz   >  $output_g/$temp/40.$temp/$file3.bismark.cov  ");  
    
    $file1 = "CHG_context_20.".$temp;
    $file2 = "CHH_context_20.".$temp;
    $file3 = "CpG_context_20.".$temp;
    system("bismark2bedGraph  --buffer_size 16G   --CX_context  --dir $output_g/$temp/20.$temp   -o  $file1.bedGraph    $output_g/$temp/20.$temp/$file1.txt.gz");
    system("bismark2bedGraph  --buffer_size 16G   --CX_context  --dir $output_g/$temp/20.$temp   -o  $file2.bedGraph    $output_g/$temp/20.$temp/$file2.txt.gz");
    system("bismark2bedGraph  --buffer_size 16G                 --dir $output_g/$temp/20.$temp   -o  $file3.bedGraph    $output_g/$temp/20.$temp/$file3.txt.gz");
    system("gunzip  --stdout    $output_g/$temp/20.$temp/$file1.bismark.cov.gz   >  $output_g/$temp/20.$temp/$file1.bismark.cov  ");  
    system("gunzip  --stdout    $output_g/$temp/20.$temp/$file2.bismark.cov.gz   >  $output_g/$temp/20.$temp/$file2.bismark.cov  ");  
    system("gunzip  --stdout    $output_g/$temp/20.$temp/$file3.bismark.cov.gz   >  $output_g/$temp/20.$temp/$file3.bismark.cov  ");  
    
    $file1 = "CHG_context_1.".$temp;
    $file2 = "CHH_context_1.".$temp;
    $file3 = "CpG_context_1.".$temp;
    system("bismark2bedGraph  --buffer_size 16G   --CX_context  --dir $output_g/$temp/1.$temp   -o  $file1.bedGraph    $output_g/$temp/1.$temp/$file1.txt.gz");
    system("bismark2bedGraph  --buffer_size 16G   --CX_context  --dir $output_g/$temp/1.$temp   -o  $file2.bedGraph    $output_g/$temp/1.$temp/$file2.txt.gz");
    system("bismark2bedGraph  --buffer_size 16G                 --dir $output_g/$temp/1.$temp   -o  $file3.bedGraph    $output_g/$temp/1.$temp/$file3.txt.gz");
    system("gunzip  --stdout    $output_g/$temp/1.$temp/$file1.bismark.cov.gz   >  $output_g/$temp/1.$temp/$file1.bismark.cov  ");  
    system("gunzip  --stdout    $output_g/$temp/1.$temp/$file2.bismark.cov.gz   >  $output_g/$temp/1.$temp/$file2.bismark.cov  ");  
    system("gunzip  --stdout    $output_g/$temp/1.$temp/$file3.bismark.cov.gz   >  $output_g/$temp/1.$temp/$file3.bismark.cov  ");  
      
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";





## END
