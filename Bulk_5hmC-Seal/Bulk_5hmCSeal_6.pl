#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.22;
## Perl5 version >= 5.22
## You can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $genome_g = '';  ## such as "mm39", "ce11", "hg38".
my $input_g  = '';  ## such as "5_rmDups/1_BWA"
my $output_g = '';  ## such as "6_MACS2/1_BWA"

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Call peaks for BAM or BAMPE by using MACS2 without input.  
        Usage:
               perl  Bulk_5hmCSeal_6.pl    [-version]    [-help]   [-genome RefGenome]    [-in inputDir]    [-out outDir]
        For instance:
               nohup time  perl  Bulk_5hmCSeal_6.pl   -genome hg38   -in 5_rmDups/1_BWA   -out 6_MACS2/1_BWA    > Bulk_5hmCSeal_6.runLog.1_BWA.txt   2>&1    &
        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.
        -help           Show this help message and exit.
        Required arguments:
        -genome RefGenome   "RefGenome" is the short name of your reference genome, such as "mm39", "ce11", "hg38".    (no default)
        -in inputDir        "inputDir" is the name of input path that contains your BAM files.  (no default)
        -out outDir         "outDir" is the name of output path that contains your running results of this step.  (no default)
        -----------------------------------------------------------------------------------------------------------
        For more details about this pipeline and other NGS data analysis piplines, please visit https://github.com/CTLife/2ndGS_Pipelines
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';

## Version Infromation
my $version = "     The 6th Step, version 1.0,  2022-09-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$genome_g = 'hg38';             ## This is only an initialization value or suggesting value, not default value.
$input_g  = '5_rmDups/1_BWA';   ## This is only an initialization value or suggesting value, not default value.
$output_g = '6_MACS2/1_BWA';    ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help   -genome   -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  Bulk_5hmCSeal_6.pl  -help' \n";
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
my $numCores_g   = 16;
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

&printVersion("macs2  --version");
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
my $genome_size_g = '';  
if($genome_g eq 'hg38')  {$genome_size_g = 'hs'; }
if($genome_g eq 'dm6' )  {$genome_size_g = 'dm'; }
if($genome_g eq 'mm39')  {$genome_size_g = 'mm'; }
print("\n\ngenome_size: $genome_size_g\n\n");
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Peak calling for narrow peaks ......";

{
&myMakeDir("$output_g");
my $outputDir1 = "$output_g/narrow_run1";
my $outputDir2 = "$output_g/narrow_run2";
my $outputDir3 = "$output_g/narrow_run3";
my $outputDir4 = "$output_g/narrow_run4";
my $outputDir5 = "$output_g/narrow_run5";
my $outputDir6 = "$output_g/narrow_run6";
my $outputDir7 = "$output_g/narrow_run7";
my $outputDir8 = "$output_g/narrow_run8";
my $outputDir9 = "$output_g/narrow_run9";
my $outputDir10 = "$output_g/narrow_run10";
my $outputDir11 = "$output_g/narrow_run11";
my $outputDir12 = "$output_g/narrow_run12";
my $outputDir13 = "$output_g/narrow_run13";
my $outputDir14 = "$output_g/narrow_run14";
my $outputDir15 = "$output_g/narrow_run15";

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
&myMakeDir($outputDir11);
&myMakeDir($outputDir12);
&myMakeDir($outputDir13);
&myMakeDir($outputDir14);
&myMakeDir($outputDir15);

for ( my $i=0; $i<=$#BAMfiles_g; $i++ ) { 
    my $temp = $BAMfiles_g[$i];
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
    &myMakeDir("$outputDir11/$temp");
    &myMakeDir("$outputDir12/$temp");
    &myMakeDir("$outputDir13/$temp");
    &myMakeDir("$outputDir14/$temp");
    &myMakeDir("$outputDir15/$temp");

    ## Duplicates  were removed by picard.
    system("macs2  callpeak   --treatment $input_g/$temp.bam       --format BAMPE   --gsize $genome_size_g   --keep-dup all     --outdir $outputDir1/$temp    --name $temp      --verbose 3         --qvalue 0.05                 >> $outputDir1/$temp.runLog    2>&1");                    
    system("macs2  callpeak   --treatment $input_g/$temp.bam       --format BAMPE   --gsize $genome_size_g   --keep-dup all     --outdir $outputDir2/$temp    --name $temp      --verbose 3         --qvalue 0.1                  >> $outputDir2/$temp.runLog    2>&1");                    
    system("macs2  callpeak   --treatment $input_g/$temp.bam       --format BAMPE   --gsize $genome_size_g   --keep-dup all     --outdir $outputDir3/$temp    --name $temp      --verbose 3         --pvalue 0.00001              >> $outputDir3/$temp.runLog    2>&1");                    
    system("macs2  callpeak   --treatment $input_g/$temp.bam       --format BAMPE   --gsize $genome_size_g   --keep-dup all     --outdir $outputDir4/$temp    --name $temp      --verbose 3         --pvalue 0.0001               >> $outputDir4/$temp.runLog    2>&1");    
    system("macs2  callpeak   --treatment $input_g/$temp.bam       --format BAMPE   --gsize $genome_size_g   --keep-dup all     --outdir $outputDir5/$temp    --name $temp      --verbose 3         --pvalue 0.001                >> $outputDir5/$temp.runLog    2>&1");   
    
    system("macs2  callpeak   --treatment $input_g/$temp.bam       --format BAMPE   --gsize $genome_size_g   --keep-dup all     --outdir $outputDir6/$temp    --name $temp      --verbose 3         --nomodel     --qvalue 0.05              >> $outputDir6/$temp.runLog    2>&1");                    
    system("macs2  callpeak   --treatment $input_g/$temp.bam       --format BAMPE   --gsize $genome_size_g   --keep-dup all     --outdir $outputDir7/$temp    --name $temp      --verbose 3         --nomodel     --qvalue 0.1               >> $outputDir7/$temp.runLog    2>&1");                    
    system("macs2  callpeak   --treatment $input_g/$temp.bam       --format BAMPE   --gsize $genome_size_g   --keep-dup all     --outdir $outputDir8/$temp    --name $temp      --verbose 3         --nomodel     --pvalue 0.00001           >> $outputDir8/$temp.runLog    2>&1");                    
    system("macs2  callpeak   --treatment $input_g/$temp.bam       --format BAMPE   --gsize $genome_size_g   --keep-dup all     --outdir $outputDir9/$temp    --name $temp      --verbose 3         --nomodel     --pvalue 0.0001            >> $outputDir9/$temp.runLog    2>&1");    
    system("macs2  callpeak   --treatment $input_g/$temp.bam       --format BAMPE   --gsize $genome_size_g   --keep-dup all     --outdir $outputDir10/$temp   --name $temp      --verbose 3         --nomodel     --pvalue 0.001             >> $outputDir10/$temp.runLog   2>&1");   
    
    system("macs2  callpeak   --treatment $input_g/$temp.bam       --format BAMPE   --gsize $genome_size_g   --keep-dup all     --outdir $outputDir11/$temp   --name $temp      --verbose 3         --nomodel  --nolambda    --qvalue 0.05                   >> $outputDir11/$temp.runLog    2>&1");                    
    system("macs2  callpeak   --treatment $input_g/$temp.bam       --format BAMPE   --gsize $genome_size_g   --keep-dup all     --outdir $outputDir12/$temp   --name $temp      --verbose 3         --nomodel  --nolambda    --qvalue 0.1                    >> $outputDir12/$temp.runLog    2>&1");                    
    system("macs2  callpeak   --treatment $input_g/$temp.bam       --format BAMPE   --gsize $genome_size_g   --keep-dup all     --outdir $outputDir13/$temp   --name $temp      --verbose 3         --nomodel  --nolambda    --pvalue 0.00001                >> $outputDir13/$temp.runLog    2>&1");                    
    system("macs2  callpeak   --treatment $input_g/$temp.bam       --format BAMPE   --gsize $genome_size_g   --keep-dup all     --outdir $outputDir14/$temp   --name $temp      --verbose 3         --nomodel  --nolambda    --pvalue 0.0001                 >> $outputDir14/$temp.runLog    2>&1");    
    system("macs2  callpeak   --treatment $input_g/$temp.bam       --format BAMPE   --gsize $genome_size_g   --keep-dup all     --outdir $outputDir15/$temp   --name $temp      --verbose 3         --nomodel  --nolambda    --pvalue 0.001                  >> $outputDir15/$temp.runLog    2>&1");   

 }

}
###################################################################################################################################################################################################








###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";


 

  
## END
