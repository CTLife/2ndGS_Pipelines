#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.22;
use  String::Similarity;

## Perl5 version >= 5.22
## You can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $genome_g = '';  ## such as "mm10", "ce11", "hg38".
my $input_g  = '';  ## such as "6-BAMPE/1_BWAmem"
my $output_g = '';  ## such as "9B-Q/1_BWAmem"

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use CISDA (ChIP-Seq Data Analyzer), version 0.9.4,  2018-02-01.
        CISDA is a Pipeline for Single-end and Paired-end ChIP-Seq Data Analysis by Integrating Lots of Softwares.
                                                            
        Step 8: Call peaks for BAM or BAMPE by using Q (http://charite.github.io/Q/).  

        Usage:
               perl  CISDA8B.pl    [-version]    [-help]   [-genome RefGenome]    [-in inputDir]    [-out outDir]
        For instance:
               perl  CISDA8B.pl   -genome hg38   -in 6-BAMPE/1_BWAmem   -out 9B-Q/1_BWAmem    > CISDA8B.runLog

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
my $version = "    The 8th Step of CISDA (ChIP-Seq Data Analyzer), version 0.9.4,  2018-02-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$genome_g = 'hg38';               ## This is only an initialization value or suggesting value, not default value.
$input_g  = '6-BAMPE/1_BWAmem';   ## This is only an initialization value or suggesting value, not default value.
$output_g = '9B-Q/1_BWAmem';      ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help   -genome   -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  CISDA8B.pl  -help' \n";
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

&printVersion("Q  --version");

###################################################################################################################################################################################################





###################################################################################################################################################################################################                 
my @groupFiles_g = ();
{
say   "\n\n\n\n\n\n##################################################################################################";
say   "Checking all the input file names ......";     
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
        $groupFiles_g[++$#groupFiles_g] = $inputFiles_g[$i];  
        $temp =~ m/^(\d+)_($pattern_g)_(Rep[1-9])/   or  die   "wrong-1: ## $temp ##";
        $temp =~ m/_(Rep[1-9])\.bam$/  or    die   "wrong-2: ## $temp ##";
        if($temp !~ m/^((\d+)_($pattern_g)_(Rep[1-9]))(_[1-2])?\.bam$/) {
             $fileNameBool = 0;
        }
}

if($fileNameBool == 1)  { say    "\n\t\tAll the file names are passed.\n";  }
@groupFiles_g   = sort(@groupFiles_g);
my $numGroup  = 0;
my $noteGroup = 0;
for ( my $i=0; $i<=$#groupFiles_g; $i++ ) { 
    $groupFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])/  or  die;
    my $n1 = $1;
    $n1>=1  or  die;
    if($noteGroup != $n1) {say "\n\t\tGroup $n1:";  $numGroup++; }
    say  "\t\t\t$groupFiles_g[$i]";
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
open(seqFiles_FH, ">", "$output2_g/BAM-Files.txt")  or  die; 
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

for ( my $i=0; $i<=$#BAMfiles_g; $i++ ) { 
   $BAMfiles_g[$i] = "$input_g/$BAMfiles_g[$i]";
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Detecting input files for each group ......";
my @input_group_g = ();
{
my $numGroup  = 0;
my $noteGroup = 0;
for ( my $i=0; $i<=$#groupFiles_g; $i++ ) { 
    $groupFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])/  or  die;
    my $n1 = $1;
    $n1>=1  or  die;
    if($noteGroup != $n1) {
        say "\n\t\tGroup $n1:";  
        $numGroup++; 
    }
    say  "\t\t\t$groupFiles_g[$i]";
    if($groupFiles_g[$i] =~ m/^(\d+)_Input-/) {
        $input_group_g[++$#input_group_g] = $groupFiles_g[$i];
    }
    $noteGroup = $n1; 
}

## ($#input_group_g == $numGroup-1)   or  die  "\n\n##$#input_group_g,\t$numGroup##\n\n";
open(runFH, ">", "$output2_g/chip_input.txt")  or die;

for ( my $i=0; $i<=$#input_group_g; $i++ ) { 
    my $group_n = $i+1;
    print  runFH   "chip input for $group_n: $input_group_g[$i]\n";
    print          "chip input for $group_n: $input_group_g[$i]\n";
}

}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $genome_size_g = '';  
if($genome_g eq 'hg38')  {$genome_size_g = 'hs'; }
if($genome_g eq 'dm6' )  {$genome_size_g = 'dm'; }
if($genome_g eq 'mm10')  {$genome_size_g = 'mm'; }
print("\n\ngenome_size: $genome_size_g\n\n");

my $format_g = "BAMPE";   ## --format {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}



sub whichNearestString  {
    my @list = @_;
    my $ref_string = $list[0];
    my @string_array = @list; 
    my $similarity_temp = 0;
    my $target_string;
    #say($ref_string);
    #say(@string_array); 
    #say( $#string_array );  
    for(my $i=1; $i<=$#string_array; $i++) {
        my $similarity_1 = String::Similarity::similarity($ref_string, $string_array[$i]);
        if($similarity_1 > $similarity_temp) {$target_string = $string_array[$i]; $similarity_temp = $similarity_1; }
    }

    return($target_string);   
}

say("\n\nInput of each sample:");
for ( my $i=0; $i<=$#groupFiles_g; $i++ ) { 
    $groupFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])/  or  die;
    my $input_1 = &whichNearestString($groupFiles_g[$i], @input_group_g ); 
    print  "\t$groupFiles_g[$i],\t$input_1 ......\n";  
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Peak calling ......";

{
&myMakeDir("$output_g/1-Peaks");
my $outputDir1 = "$output_g/1-Peaks/run1";
my $outputDir2 = "$output_g/1-Peaks/run2";
my $outputDir3 = "$output_g/1-Peaks/run3";
my $outputDir4 = "$output_g/1-Peaks/run4";
my $outputDir5 = "$output_g/1-Peaks/run5";
my $outputDir6 = "$output_g/1-Peaks/run6";
my $outputDir7 = "$output_g/1-Peaks/run7";
my $outputDir8 = "$output_g/1-Peaks/run8";
my $outputDir9 = "$output_g/1-Peaks/run9";
my $outputDir10 = "$output_g/1-Peaks/run10";

&myMakeDir($outputDir1);
&myMakeDir($outputDir2);
&myMakeDir($outputDir3);
&myMakeDir($outputDir4);
&myMakeDir($outputDir5);
&myMakeDir($outputDir6);
&myMakeDir($outputDir7);
&myMakeDir($outputDir8);
&myMakeDir($outputDir9);
&myMakeDir($outputDir10);

for ( my $i=0; $i<=$#groupFiles_g; $i++ ) { 
    $groupFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])/  or  die;
    my $input_1 = &whichNearestString($groupFiles_g[$i], @input_group_g ); 
    print  "\t$groupFiles_g[$i],\t$input_1 ......\n";  
    my $temp = $groupFiles_g[$i];
    $temp =~ s/\.bam$//  or  die;
    &myMakeDir("$outputDir1/$temp");
    &myMakeDir("$outputDir2/$temp");
    &myMakeDir("$outputDir3/$temp");
    &myMakeDir("$outputDir4/$temp");
    &myMakeDir("$outputDir5/$temp");
    &myMakeDir("$outputDir6/$temp");
    &myMakeDir("$outputDir7/$temp");
    &myMakeDir("$outputDir8/$temp");
    &myMakeDir("$outputDir9/$temp");
    &myMakeDir("$outputDir10/$temp");

    system(`Q   --treatment-sample  $input_g/$temp.bam     --control-sample $input_g/$input_1    --fragment-length-average -1     --fragment-length-deviation 50    --p-value-cutoff 1        --out-prefix  $outputDir1/$temp/$temp       --thread-num $numCores_g   --verbose  --step-num 1000       >> $outputDir1/$temp.runLog    2>&1 `);                    
    system(`Q   --treatment-sample  $input_g/$temp.bam     --control-sample $input_g/$input_1    --fragment-length-average -1     --fragment-length-deviation 50    --p-value-cutoff 2        --out-prefix  $outputDir2/$temp/$temp       --thread-num $numCores_g   --verbose  --step-num 1000       >> $outputDir2/$temp.runLog    2>&1 `);                    
    system(`Q   --treatment-sample  $input_g/$temp.bam     --control-sample $input_g/$input_1    --fragment-length-average -1     --fragment-length-deviation 50    --p-value-cutoff 3        --out-prefix  $outputDir3/$temp/$temp       --thread-num $numCores_g   --verbose  --step-num 1000       >> $outputDir3/$temp.runLog    2>&1 `);                    
    system(`Q   --treatment-sample  $input_g/$temp.bam     --control-sample $input_g/$input_1    --fragment-length-average -1     --fragment-length-deviation 50    --p-value-cutoff 4        --out-prefix  $outputDir4/$temp/$temp       --thread-num $numCores_g   --verbose  --step-num 1000       >> $outputDir4/$temp.runLog    2>&1 `);                    
    system(`Q   --treatment-sample  $input_g/$temp.bam     --control-sample $input_g/$input_1    --fragment-length-average -1     --fragment-length-deviation 50    --p-value-cutoff 5        --out-prefix  $outputDir5/$temp/$temp       --thread-num $numCores_g   --verbose  --step-num 1000       >> $outputDir5/$temp.runLog    2>&1 `);                    
    system(`Q   --treatment-sample  $input_g/$temp.bam     --control-sample $input_g/$input_1    --fragment-length-average -1     --fragment-length-deviation 50    --p-value-cutoff 6        --out-prefix  $outputDir6/$temp/$temp       --thread-num $numCores_g   --verbose  --step-num 1000       >> $outputDir6/$temp.runLog    2>&1 `);                    
    system(`Q   --treatment-sample  $input_g/$temp.bam     --control-sample $input_g/$input_1    --fragment-length-average -1     --fragment-length-deviation 50    --p-value-cutoff 7        --out-prefix  $outputDir7/$temp/$temp       --thread-num $numCores_g   --verbose  --step-num 1000       >> $outputDir7/$temp.runLog    2>&1 `);                    
    system(`Q   --treatment-sample  $input_g/$temp.bam     --control-sample $input_g/$input_1    --fragment-length-average -1     --fragment-length-deviation 50    --p-value-cutoff 8        --out-prefix  $outputDir8/$temp/$temp       --thread-num $numCores_g   --verbose  --step-num 1000       >> $outputDir8/$temp.runLog    2>&1 `);                    
    system(`Q   --treatment-sample  $input_g/$temp.bam     --control-sample $input_g/$input_1    --fragment-length-average -1     --fragment-length-deviation 50    --p-value-cutoff 9        --out-prefix  $outputDir9/$temp/$temp       --thread-num $numCores_g   --verbose  --step-num 1000       >> $outputDir9/$temp.runLog    2>&1 `);                    
    system(`Q   --treatment-sample  $input_g/$temp.bam     --control-sample $input_g/$input_1    --fragment-length-average -1     --fragment-length-deviation 50    --p-value-cutoff 10       --out-prefix  $outputDir10/$temp/$temp      --thread-num $numCores_g   --verbose  --step-num 1000       >> $outputDir10/$temp.runLog   2>&1 `);                    
     
    system("Rscript  $outputDir1/$temp/$temp-Q-binding-characteristics.R");
    system("Rscript  $outputDir2/$temp/$temp-Q-binding-characteristics.R");
    system("Rscript  $outputDir3/$temp/$temp-Q-binding-characteristics.R");
    system("Rscript  $outputDir4/$temp/$temp-Q-binding-characteristics.R");
    system("Rscript  $outputDir5/$temp/$temp-Q-binding-characteristics.R");
    system("Rscript  $outputDir6/$temp/$temp-Q-binding-characteristics.R");
    system("Rscript  $outputDir7/$temp/$temp-Q-binding-characteristics.R");
    system("Rscript  $outputDir8/$temp/$temp-Q-binding-characteristics.R");
    system("Rscript  $outputDir9/$temp/$temp-Q-binding-characteristics.R");
    system("Rscript  $outputDir10/$temp/$temp-Q-binding-characteristics.R");

 }
  
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";


 

  
## END
