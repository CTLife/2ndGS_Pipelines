#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.20;       
## perl5 version >= 5.20,   you can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
## Help Infromation
my $HELP_g = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use RASDA (RNA-Seq Data Analyzer), version 0.7.3, 2016-06-01.      
        RASDA is a Pipeline for Single-end and Paired-end RNA-Seq Data Analysis by Integrating Lots of Softwares.

        Step 8: Convert bam to bigwig by using bedtools and bedGraphToBigWig.
                Only reads density will be computed, no extend and no shift.
                So single-end and paired-end reads will not be distinguished.
        Usage:  
               perl  RASDA8.pl    [-version]    [-help]    [-in inputDir]    [-out outDir]   
        For instance: 
               perl  RASDA8.pl    -in 6-FinalBAM/5_STAR     -out 9-BigWig/5_STAR   >> RASDA8.runLog  2>&1
     
        -------------------------------------------------------------------------------------------------------------------
        Optional arguments:
        -v, --version        Show version number of this program and exit.

        -h, --help           Show this help message and exit.

        Required arguments:
        -in  inputDir       inputDir is the name of your input folder that contains your bam files, the suffix of the bam files must be ".bam".  (no default)

        -out outDir         outDir is the name of your output folder that contains running results (bigwig format) of this step.  (no default)
        ------------------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines such as CISDA, MESDA and HISDA,
        please see the web site: https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies 
        and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.    
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------  
';

## Version Information
my $version_g = "  The Eighth Step of RASDA (RNA-Seq Data Analyzer), version 0.7.3, 2016-06-01.";

## Keys and Values
if ($#ARGV   == -1) { print  "\n$HELP_g\n\n";  exit 0; }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0) { @ARGV = (@ARGV, "help");         }       ## when the number of command argumants is odd. 
my %args = @ARGV;

## Initialize  Variables
my $input_g   = '6-FinalBAM/5_STAR';       ## This is only an initialization  value or suggesting value, not default value.
my $output_g  = '9-BigWig/5_STAR';         ## This is only an initialization  value or suggesting value, not default value.

## Available Arguments
my $available = "  -version    -help    -in    -out    ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if($available !~ m/\s$key\s/) {print  "    Cann't recognize $key !!\n";  $boole_g = 1; }
}
if($boole_g == 1) {
    print "\n    The Command Line Arguments are wrong!\n";
    print   '    Please see help message by using "perl  RASDA8.pl  -help". ';
    print "\n\n";
    exit 0;
}

## Get Arguments
if ( exists $args{'-version' }   )     { say  "\n$version_g\n";   exit 0; }
if ( exists $args{'-help'    }   )     { say  "\n$HELP_g\n";      exit 0; }
if ( exists $args{'-in'      }   )     { $input_g  = $args{'-in' };       }else{say   "\n -in  is required.\n";          say  "\n$HELP_g\n";    exit 0; }
if ( exists $args{'-out'     }   )     { $output_g = $args{'-out'};       }else{say   "\n -out is required.\n";          say  "\n$HELP_g\n";    exit 0; }

## Conditions
$input_g   =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";
$output_g  =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";

## Print Command Arguments to Standard Output
print  "\n\n
        ################ Arguments ###############################
                Input  folder:  $input_g
                Output folder:  $output_g
        ###############################################################  
\n\n";
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
my $pattern_g  = "[-.0-9A-Za-z]+";
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
        next unless $inputFiles_g[$i] !~ m/^removed_/;
        say   "\t......$inputFiles_g[$i]" ; 
        my $temp = $inputFiles_g[$i]; 
        $groupFiles[++$#groupFiles] = $inputFiles_g[$i];  
        $temp =~ m/^(\d+)_($pattern_g)_(Rep[1-9])/   or  die   "wrong-1: ## $temp ##";
        $temp =~ m/_(Rep[1-9])\.bam$/  or    die   "wrong-2: ## $temp ##";
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
    system("\t\tbedtools   genomecov  -bg   -split   -scale $Temp_scale   -ibam $input_g/$temp.bam   -g 0-Other/Shortcuts/mm10.chrom.sizes    > $output_g/$temp.bedGraph");
    sleep(3);
    system("bedGraphToBigWig    $output_g/$temp.bedGraph     0-Other/Shortcuts/mm10.chrom.sizes      $output_g/$temp.bw");
}
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n";





## END
