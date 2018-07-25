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
my $output_g = '';  ## such as "9C-GEM/1_BWAmem"

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use CISDA (ChIP-Seq Data Analyzer), version 0.9.4,  2018-02-01.
        CISDA is a Pipeline for Single-end and Paired-end ChIP-Seq Data Analysis by Integrating Lots of Softwares.
                                                            
        Step 8: Call peaks for BAM or BAMPE by using GEM (https://groups.csail.mit.edu/cgs/gem/).   

        Usage:
               perl  CISDA8C.pl    [-version]    [-help]   [-genome RefGenome]    [-in inputDir]    [-out outDir]
        For instance:
               perl  CISDA8C.pl   -genome hg38   -in 6-BAMPE/1_BWAmem   -out 9C-GEM/1_BWAmem    > CISDA8C.runLog

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
$output_g = '9C-GEM/1_BWAmem';    ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help   -genome   -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  CISDA8C.pl  -help' \n";
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

my $GEM_g = &fullPathApp("gem.jar");  
&printVersion("java  -jar   $GEM_g   --help");

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

&myMakeDir($outputDir1);
&myMakeDir($outputDir2);

my $temp_dir_g = "$genome_g.yongpeng.forGEM";
&myMakeDir($temp_dir_g);
system("gt  splitfasta   -splitdesc $temp_dir_g    0-Other/Shortcuts/$genome_g/$genome_g.fasta");

for ( my $i=0; $i<=$#groupFiles_g; $i++ ) { 
    $groupFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])/  or  die;
    my $input_1 = &whichNearestString($groupFiles_g[$i], @input_group_g ); 
    print  "\t$groupFiles_g[$i],\t$input_1 ......\n";
    my $temp = $groupFiles_g[$i];
    $temp =~ s/\.bam$//  or  die;
    &myMakeDir("$outputDir1/$temp");
    &myMakeDir("$outputDir2/$temp");

    system(`java  -Xmx10G  -jar   $GEM_g   --d /home/yp/AnnotationBED/GEM_ReadDistributionFile/Read_Distribution_default_ChIP-seq.txt    --g 0-Other/Shortcuts/$genome_g/$genome_g.chrom.sizes     --genome  $temp_dir_g   --t $numCores_g    --expt $input_g/$temp.bam    --ctrl $input_g/$input_1   --f SAM   --out $outputDir1/$temp/$temp    --k_min 6    --k_max 13   --q 2    --k_seqs 5000     >> $outputDir1/$temp.runLog    2>&1 `);                    
    system(`java  -Xmx10G  -jar   $GEM_g   --d /home/yp/AnnotationBED/GEM_ReadDistributionFile/Read_Distribution_default_ChIP-seq.txt    --g 0-Other/Shortcuts/$genome_g/$genome_g.chrom.sizes     --genome  $temp_dir_g   --t $numCores_g    --expt $input_g/$temp.bam    --ctrl $input_g/$input_1   --f SAM   --out $outputDir2/$temp/$temp    --k_min 6    --k_max 13   --q 3    --k_seqs 2000     >> $outputDir2/$temp.runLog    2>&1 `);                    

 }

system("rm  -rf  $temp_dir_g");
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";


 

  
## END
