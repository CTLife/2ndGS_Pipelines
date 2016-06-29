#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.20;
## Perl5 version >= 5.20,   you can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
## Help Infromation
my $HELP_g = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use CISDA (ChIP-Seq Data Analyzer), version 0.7.3, 2016-06-01.      
        CISDA is a Pipeline for Single-end and Paired-end ChIP-Seq Data Analysis by Integrating Lots of Softwares.

        Step 9: Call peaks for BAM or BAMPE by using MACS2.
        Usage:  
               perl  CISDA9.pl    [-v]    [-h]    [-i inputDir]    [-o outDir]   [-f format]  
        For instance: 
                     perl  CISDA9.pl    -i 6-FinalBAM/1_BWA          -o 10-MACS2    -f BAM                    
                     perl  CISDA9.pl    --input 6-FinalBAM/1_BWA     --output 10-MACS2  --format BAM   
                     perl  CISDA9.pl    --input 6-FinalBAM/1_BWA     --output 10-MACS2  --format BAM     >> CISDA9.runLog  2>&1
     
        -------------------------------------------------------------------------------------------------------------------
        Optional arguments:
        -v, --version        Show version number of this program and exit.

        -h, --help           Show this help message and exit.

        Required arguments:
        -i inputDir,  --input inputDir        inputDir is the name of your input folder that contains your BAM files,
                                              the suffix of the BAM files must be ".bam".    (no default)

        -o outDir,  --output outDir           outDir is the name of your output folder that contains running 
                                              results of this step.      (no default) 

        -f format,  --format format           format should be BAM or BAMPE.   (no default) 
        ------------------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines such as RASDA, MESDA and HISDA,
        please see the web site: https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies 
        and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.    
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------  
';

## Version Information
my $version_g = "  The Ninth Step of CISDA (ChIP-Seq Data Analyzer), version 0.7.3, 2016-06-01.";

## Keys and Values
if ($#ARGV   == -1) { print  "\n$HELP_g\n\n";  exit 0; }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0) { @ARGV = (@ARGV, "-h");           }       ## when the number of command argumants is odd. 
my %args = @ARGV;

## Initialize  Variables
my $input_g   = '6-FinalBAM/1_BWA';      ## This is only an initialization  value or suggesting value, not default value.
my $output_g  = '10-MACS2';              ## This is only an initialization  value or suggesting value, not default value.
my $format_g  = 'BAM';                   ## This is only an initialization  value or suggesting value, not default value.

## Available Arguments
my $available = "  -v  --version    -h  --help    -i  --input    -o    --output     -f  --format     ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if($available !~ m/\s$key\s/) {print  "    Cann't recognize $key !!\n";  $boole_g = 1; }
}
if($boole_g == 1) {
    print "\n    The Command Line Arguments are wrong!\n";
    print   '    Please see help message by using "perl  CISDA9.pl  -h". ';
    print "\n\n";
    exit 0;
}

## Get Arguments
if ( ( exists $args{'-v' } )  or  ( exists $args{'--version'      } )  )     { print  "\n$version_g\n\n";    exit 0; }
if ( ( exists $args{'-h' } )  or  ( exists $args{'--help'         } )  )     { print  "\n$HELP_g\n\n";       exit 0; }
if ( ( exists $args{'-i' } )  or  ( exists $args{'--input'        } )  )     { ($input_g  = $args{'-i' })   or  ($input_g  = $args{'--input'  });  }else{print   "\n -i or --input   is required.\n\n";   print  "\n$HELP_g\n\n";       exit 0; }                                               
if ( ( exists $args{'-o' } )  or  ( exists $args{'--output'       } )  )     { ($output_g = $args{'-o' })   or  ($output_g = $args{'--output' });  }else{print   "\n -o or --output  is required.\n\n";   print  "\n$HELP_g\n\n";       exit 0; }      
if ( ( exists $args{'-f' } )  or  ( exists $args{'--format'       } )  )     { ($format_g = $args{'-f' })   or  ($format_g = $args{'--format' });  }else{print   "\n -f or --format  is required.\n\n";   print  "\n$HELP_g\n\n";       exit 0; }      

## Conditions
$input_g   =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";
$output_g  =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";
($format_g eq "BAM")  or  ($format_g  eq "BAMPE")   ||  die   "\n$HELP_g\n\n";

## Print Command Arguments to Standard Output
print  "\n\n
        ################ Your Arguments ###############################
                Input  folder:  $input_g
                Output folder:  $output_g
                BAM file format:  $format_g
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
my $pattern_g = "[-.0-9A-Za-z]+";
my $numCores_g = 4;
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
&printVersion("macs2  --version");
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
open(seqFiles_FH, ">", "$output2_g/BAM-Files.txt")  or  die; 
for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {     
    next unless $inputFiles_g[$i] =~ m/\.bam$/;
    next unless $inputFiles_g[$i] !~ m/^[.]/;
    next unless $inputFiles_g[$i] !~ m/[~]$/;
    next unless $inputFiles_g[$i] !~ m/^unpaired/;
    next unless $inputFiles_g[$i] !~ m/^removed_/;
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
($#input_group_g == $numGroup-1)   or  die  "\n\n##$#input_group_g,\t$numGroup##\n\n";
open(runFH, ">", "$output2_g/chip_input.txt")  or die;
for ( my $i=0; $i<=$#input_group_g; $i++ ) { 
    my $group_n = $i+1;
    print  runFH   "chip input for $group_n: $input_group_g[$i]\n";
    print          "chip input for $group_n: $input_group_g[$i]\n";
}
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
my $outputDir10= "$output_g/1-Peaks/run10";
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
    my $n1 = $1;
    $n1>=1  or  die;
    $n1--;
    print  "\t$groupFiles_g[$i],\t$input_group_g[$n1] ......\n";
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
    system(`macs2  callpeak   --treatment $input_g/$temp.bam  --control $input_g/$input_group_g[$n1]   --format $format_g   --gsize mm   --keep-dup 1     --outdir $outputDir1/$temp    --name $temp   --bw 300    --mfold 5 50   --qvalue 0.05                 >> $outputDir1/$temp.runLog    2>&1`);                    
    system(`macs2  callpeak   --treatment $input_g/$temp.bam  --control $input_g/$input_group_g[$n1]   --format $format_g   --gsize mm   --keep-dup 1     --outdir $outputDir2/$temp    --name $temp   --bw 300    --mfold 5 50   --qvalue 0.05   --broad       >> $outputDir2/$temp.runLog    2>&1`);                    
    system(`macs2  callpeak   --treatment $input_g/$temp.bam  --control $input_g/$input_group_g[$n1]   --format $format_g   --gsize mm   --keep-dup 1     --outdir $outputDir3/$temp    --name $temp   --bw 300    --mfold 5 50   --pvalue 0.001                >> $outputDir3/$temp.runLog    2>&1`);                    
    system(`macs2  callpeak   --treatment $input_g/$temp.bam  --control $input_g/$input_group_g[$n1]   --format $format_g   --gsize mm   --keep-dup 1     --outdir $outputDir4/$temp    --name $temp   --bw 300    --mfold 5 50   --pvalue 0.001  --broad       >> $outputDir4/$temp.runLog    2>&1`);                    
    system(`macs2  callpeak   --treatment $input_g/$temp.bam  --control $input_g/$input_group_g[$n1]   --format $format_g   --gsize mm   --keep-dup 1     --outdir $outputDir5/$temp    --name $temp   --bw 300    --mfold 5 50   --pvalue 0.01   --broad       >> $outputDir5/$temp.runLog    2>&1`);                    
    system(`macs2  callpeak   --treatment $input_g/$temp.bam  --control $input_g/$input_group_g[$n1]   --format $format_g   --gsize mm   --keep-dup 1     --outdir $outputDir6/$temp    --name $temp   --bw 300    --mfold 5 50   --pvalue 0.05   --broad       >> $outputDir6/$temp.runLog    2>&1`);                    
    system(`macs2  callpeak   --treatment $input_g/$temp.bam  --control $input_g/$input_group_g[$n1]   --format $format_g   --gsize mm   --keep-dup 1     --outdir $outputDir7/$temp    --name $temp   --bw 300    --mfold 5 50   --qvalue 0.05   --nolambda    >> $outputDir7/$temp.runLog    2>&1`); 
    system(`macs2  callpeak   --treatment $input_g/$temp.bam  --control $input_g/$input_group_g[$n1]   --format $format_g   --gsize mm   --keep-dup 1     --outdir $outputDir8/$temp    --name $temp   --bw 300    --mfold 5 50   --qvalue 0.05   --nomodel     >> $outputDir8/$temp.runLog    2>&1`);                                       
    system(`macs2  callpeak   --treatment $input_g/$temp.bam  --control $input_g/$input_group_g[$n1]   --format $format_g   --gsize mm   --keep-dup 1     --outdir $outputDir9/$temp    --name $temp   --bw 300    --mfold 5 50   --qvalue 0.05   --bdg         >> $outputDir9/$temp.runLog    2>&1`);                    
    system(`macs2  callpeak   --treatment $input_g/$temp.bam                                           --format $format_g   --gsize mm   --keep-dup 1     --outdir $outputDir10/$temp   --name $temp   --bw 300    --mfold 5 50   --qvalue 0.05                 >> $outputDir10/$temp.runLog   2>&1`);                    
}

}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Caculate fragment size by using macs2 predictd ......";
{
my $outputDir1 = "$output_g/2-fragmentSize/1-default";
my $outputDir2 = "$output_g/2-fragmentSize/2-run2";
&myMakeDir($outputDir1);
&myMakeDir($outputDir2);

for(my $i=0; $i<=$#BAMfiles_g; $i++) {
    next unless $BAMfiles_g[$i] =~ m/\.bam$/;
    next unless $BAMfiles_g[$i] !~ m/^[.]/;
    next unless $BAMfiles_g[$i] !~ m/[~]$/;
    my $temp = $BAMfiles_g[$i];
    $temp =~ s/\.bam$//  or  die;
    system(`macs2  predictd   --ifile $input_g/$temp.bam   --format BAM   --gsize mm    --bw 300    --mfold 5 50      --outdir $outputDir1    --rfile $temp.R     >> $outputDir1/$temp.runLog   2>&1`);                    
    system(`macs2  predictd   --ifile $input_g/$temp.bam   --format BAM   --gsize mm    --bw 500    --mfold 3 300     --outdir $outputDir2    --rfile $temp.R     >> $outputDir2/$temp.runLog   2>&1`);                    
}

}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n";





## END
