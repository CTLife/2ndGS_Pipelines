#!/usr/bin/env perl5
use  strict;
use  warnings;
use  v5.20;





###################################################################################################################################################################################################
###################################################################################################################################################################################################
########## Help Infromation ##########
my $HELP_g = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use CISDA (ChIP-Seq Data Analyzer), version 0.62, 2016-01-13.      
        CISDA is a Pipeline for Single-end and Paired-end ChIP-Seq Data Analysis by Integrating Lots of Softwares.

        Step 10: Call peak by using MACS1.4.2.
        Usage:  
               perl  CISDA-10.pl    [-v]    [-h]    [-i inputDir]    [-o outDir]   
        For instance: 
                     perl  CISDA-10.pl    -i 7-FinalBAM/1-Subread          -o 11-MACS1.4.2                      
                     perl  CISDA-10.pl    --input 7-FinalBAM/1-Subread     --output 11-MACS1.4.2     
                     perl  CISDA-10.pl    --input 7-FinalBAM/1-Subread     --output 11-MACS1.4.2    >> CISDA-10.runLog  2>&1
     
        -------------------------------------------------------------------------------------------------------------------
        Optional arguments:
        -v, --version        Show version number of this program and exit.

        -h, --help           Show this help message and exit.

        Required arguments:
        -i inputDir,  --input inputDir        inputDir is the name of your input folder that contains your BAM files,
                                              the suffix of the BAM files must be ".bam".    (no default)

        -o outDir,  --output outDir           outDir is the name of your output folder that contains running 
                                              results of this step.      (no default) 
        ------------------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines such as RASDA, MESDA and HISDA,
        please see the web site: https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies 
        and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.    
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------  
';


########## Version Infromation ##########
my $version_g = "  The Tenth Step of CISDA (ChIP-Seq Data Analyzer), version 0.62, 2016-01-13.";


########## Keys and Values ##########
if ($#ARGV   == -1) { print  "\n$HELP_g\n\n";  exit 0; }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0) { @ARGV = (@ARGV, "-h");           }       ## when the number of command argumants is odd. 
my %args = @ARGV;


########## Initialize  Variables ##########
my $input_g   = '7-FinalBAM/1-Subread';      ## This is only an initialization  value or suggesting value, not default value.
my $output_g  = '11-MACS1.4.2';                  ## This is only an initialization  value or suggesting value, not default value.


########## Available Arguments ##########
my $available = "  -v  --version    -h  --help    -i  --input    -o    --output     -L  --fragLen     ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if($available !~ m/\s$key\s/) {print  "    Cann't recognize $key !!\n";  $boole_g = 1; }
}
if($boole_g == 1) {
    print "\n    The Command Line Arguments are wrong!\n";
    print   '    Please see help message by using "perl  CISDA-10.pl  -h". ';
    print "\n\n";
    exit 0;
}


########## Get Arguments ##########
if ( ( exists $args{'-v' } )  or  ( exists $args{'--version'      } )  )     { print  "\n$version_g\n\n";    exit 0; }
if ( ( exists $args{'-h' } )  or  ( exists $args{'--help'         } )  )     { print  "\n$HELP_g\n\n";       exit 0; }
if ( ( exists $args{'-i' } )  or  ( exists $args{'--input'        } )  )     { ($input_g  = $args{'-i' })   or  ($input_g  = $args{'--input'  });  }else{print   "\n -i or --input   is required.\n\n";   print  "\n$HELP_g\n\n";       exit 0; }                                               
if ( ( exists $args{'-o' } )  or  ( exists $args{'--output'       } )  )     { ($output_g = $args{'-o' })   or  ($output_g = $args{'--output' });  }else{print   "\n -o or --output  is required.\n\n";   print  "\n$HELP_g\n\n";       exit 0; }      


########### Conditions #############
$input_g   =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";
$output_g  =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";


######### Print Command Arguments to Standard Output ###########
print  "\n\n
        ################ Your Arguments ###############################
                Input  folder:  $input_g
                Output folder:  $output_g
        ###############################################################  
\n\n";


###################################################################################################################################################################################################
###################################################################################################################################################################################################





print "\n\n\n\n\n##################################################################################################\n";
print   "\nRunning......\n";
my $output2_g = "$output_g/Results";
if ( !(-e $output_g) )   { mkdir $output_g    ||  die; }
if ( !(-e $output2_g))   { mkdir $output2_g   ||  die; }
(-e $output_g)   ||  die;
my $pattern = "[-.0-9A-Za-z]+";







print "\n\n\n\n\n##################################################################################################\n";
print("\nChecking all the necessary softwares in this step......\n");

system("macs14  --version          >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'     >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");






opendir(my $DH_input, $input_g) || die;     
my @inputFiles = readdir($DH_input);

print("\nChecking all the input file names......\n");
my $fileNameBool = 1;
for ( my $i=0; $i<=$#inputFiles; $i++ ) {  
        next unless $inputFiles[$i] =~ m/^(\d{2})_/;   
        next unless $inputFiles[$i] =~ m/\.bam$/;   
        next unless $inputFiles[$i] !~ m/^[.]/;
        next unless $inputFiles[$i] !~ m/[~]$/;
        next unless $inputFiles[$i] !~ m/^unpaired/;
        my $temp = $inputFiles[$i]; 
        $temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]).bam$/   or  die  "wrong-1: ## $temp ##";
        if($temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])\.bam$/) {
             print("$inputFiles[$i]......\n");
        }else{
             $fileNameBool = 0;
        }
}
if($fileNameBool == 1)  {print("    All the file names are passed.\n\n");}





my $NGSinput = "";
my $bool = 0;
for(my $i=0; $i<=$#inputFiles; $i++) {
    next unless $inputFiles[$i] =~ m/\.bam$/;
    next unless $inputFiles[$i] !~ m/^[.]/;
    next unless $inputFiles[$i] !~ m/[~]$/;    
    if($inputFiles[$i] =~ m/^\d{2}_Input_/i) {
        $NGSinput = "$input_g/$inputFiles[$i]";  
        $bool++;
    }  
}
if ($bool != 1) { print  "#$bool#\n";  die;}
open(runFH, ">", "$output2_g/runLog.txt")  or die;
print   runFH   "\nNGS input file: $NGSinput\n\n";
print   "        \nNGS input file: $NGSinput\n\n";







{####################

my $outputDir1 = "$output_g/run1";
my $outputDir2 = "$output_g/run2";
my $outputDir3 = "$output_g/run3";
my $outputDir4 = "$output_g/run4";
my $outputDir5 = "$output_g/run5";
my $outputDir6 = "$output_g/run6";
my $outputDir7 = "$output_g/run7";
my $outputDir8 = "$output_g/run8";
my $outputDir9 = "$output_g/run9";
my $outputDir10= "$output_g/run10";
if (!(-e $outputDir1))   {mkdir $outputDir1    or  die; } 
if (!(-e $outputDir2))   {mkdir $outputDir2    or  die; } 
if (!(-e $outputDir3))   {mkdir $outputDir3    or  die; } 
if (!(-e $outputDir4))   {mkdir $outputDir4    or  die; } 
if (!(-e $outputDir5))   {mkdir $outputDir5    or  die; } 
if (!(-e $outputDir6))   {mkdir $outputDir6    or  die; } 
if (!(-e $outputDir7))   {mkdir $outputDir7    or  die; } 
if (!(-e $outputDir8))   {mkdir $outputDir8    or  die; } 
if (!(-e $outputDir9))   {mkdir $outputDir9    or  die; } 
if (!(-e $outputDir10))  {mkdir $outputDir10   or  die; } 


print "\n\n        Peak Calling......";
for(my $i=0; $i<=$#inputFiles; $i++) {
    next unless $inputFiles[$i] =~ m/\.bam$/;
    next unless $inputFiles[$i] !~ m/^[.]/;
    next unless $inputFiles[$i] !~ m/[~]$/;
    my $temp = $inputFiles[$i];
    $temp =~ s/\.bam$//  or  die;
    print  runFH  "\n        $temp......\n";
    print         "\n        $temp......\n";
    if ( !(-e "$outputDir1/$temp") )   {mkdir   "$outputDir1/$temp"    or  die; } 
    if ( !(-e "$outputDir2/$temp") )   {mkdir   "$outputDir2/$temp"    or  die; } 
    if ( !(-e "$outputDir3/$temp") )   {mkdir   "$outputDir3/$temp"    or  die; } 
    if ( !(-e "$outputDir4/$temp") )   {mkdir   "$outputDir4/$temp"    or  die; } 
    if ( !(-e "$outputDir5/$temp") )   {mkdir   "$outputDir5/$temp"    or  die; } 
    if ( !(-e "$outputDir6/$temp") )   {mkdir   "$outputDir6/$temp"    or  die; } 
    if ( !(-e "$outputDir7/$temp") )   {mkdir   "$outputDir7/$temp"    or  die; } 
    if ( !(-e "$outputDir8/$temp") )   {mkdir   "$outputDir8/$temp"    or  die; } 
    if ( !(-e "$outputDir9/$temp") )   {mkdir   "$outputDir9/$temp"    or  die; } 
    if ( !(-e "$outputDir10/$temp") )  {mkdir   "$outputDir10/$temp"   or  die; } 

    system(`macs14   --treatment=$input_g/$temp.bam   --control=$NGSinput   --name=$outputDir1/$temp/$temp   --format=BAM   --gsize=mm  --bw=300  --pvalue=1e-5   --mfold=10,30    --diag             >> $outputDir1/$temp.runLog   2>&1`);                                                                                                           
    system(`macs14   --treatment=$input_g/$temp.bam   --control=$NGSinput   --name=$outputDir2/$temp/$temp   --format=BAM   --gsize=mm  --bw=300  --pvalue=1e-5   --mfold=10,30    --keep-dup=auto    >> $outputDir2/$temp.runLog   2>&1`);   
    system(`macs14   --treatment=$input_g/$temp.bam   --control=$NGSinput   --name=$outputDir3/$temp/$temp   --format=BAM   --gsize=mm  --bw=300  --pvalue=1e-5   --mfold=10,30    --keep-dup=all     >> $outputDir3/$temp.runLog   2>&1`); 
                                                                                                                                                                                                                                                  
    system(`macs14   --treatment=$input_g/$temp.bam   --control=$NGSinput   --name=$outputDir4/$temp/$temp   --format=BAM   --gsize=mm  --bw=300  --pvalue=0.01   --mfold=10,30    --keep-dup=1       >> $outputDir4/$temp.runLog   2>&1`);                                                                                                           
    system(`macs14   --treatment=$input_g/$temp.bam   --control=$NGSinput   --name=$outputDir5/$temp/$temp   --format=BAM   --gsize=mm  --bw=300  --pvalue=0.01   --mfold=10,30    --keep-dup=auto    >> $outputDir5/$temp.runLog   2>&1`);   
    system(`macs14   --treatment=$input_g/$temp.bam   --control=$NGSinput   --name=$outputDir6/$temp/$temp   --format=BAM   --gsize=mm  --bw=300  --pvalue=0.01   --mfold=10,30    --keep-dup=all     >> $outputDir6/$temp.runLog   2>&1`);                                                                                                                                                                                                                                                   

    system(`macs14   --treatment=$input_g/$temp.bam   --control=$NGSinput   --name=$outputDir7/$temp/$temp   --format=BAM   --gsize=mm  --bw=300  --pvalue=1e-3   --mfold=5,50     --diag             >> $outputDir7/$temp.runLog   2>&1`);                                                                                                           
    system(`macs14   --treatment=$input_g/$temp.bam   --control=$NGSinput   --name=$outputDir8/$temp/$temp   --format=BAM   --gsize=mm  --bw=300  --pvalue=0.05   --mfold=5,50     --keep-dup=all     >> $outputDir8/$temp.runLog   2>&1`);   
    system(`macs14   --treatment=$input_g/$temp.bam   --control=$NGSinput   --name=$outputDir9/$temp/$temp   --format=BAM   --gsize=mm  --bw=300  --pvalue=1e-3   --mfold=5,50     --keep-dup=1       >> $outputDir9/$temp.runLog   2>&1`);                                                                                                                                                                                                                                                   
    system(`macs14   --treatment=$input_g/$temp.bam                         --name=$outputDir10/$temp/$temp  --format=BAM   --gsize=mm  --bw=300  --pvalue=1e-3   --mfold=5,50     --keep-dup=1       >> $outputDir10/$temp.runLog  2>&1`);                                                                                                                                                                                                                                                   
}

}##########################







print "\n\n\n\n\n##################################################################################################\n";
print "\n\n        Job Done! Cheers! \n\n";





























