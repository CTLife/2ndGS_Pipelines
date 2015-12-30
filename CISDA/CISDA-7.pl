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

        Step 7: Convert BAM to BIGWIG.
        Usage:  
               perl  CISDA-7.pl    [-v]    [-h]    [-i inputDir]    [-o outDir]   [-L n]
        For instance: 
                     perl  CISDA-7.pl    -i 7-FinalBAM/1-Subread          -o 8-BigWig          -L 200            
                     perl  CISDA-7.pl    --input 7-FinalBAM/1-Subread     --output 8-BigWig    --fragLen 200 
                     perl  CISDA-7.pl    --input 7-FinalBAM/1-Subread     --output 8-BigWig    --fragLen 200  >> CISDA-7.runLog  2>&1
     
        -------------------------------------------------------------------------------------------------------------------
        Optional arguments:
        -v, --version        Show version number of this program and exit.

        -h, --help           Show this help message and exit.

        Required arguments:
        -i inputDir,  --input inputDir        inputDir is the name of your input folder that contains your BAM files,
                                              the suffix of the BAM files must be ".bam".    (no default)

        -o outDir,  --output outDir           outDir is the name of your output folder that contains running 
                                              results (Bigwig format) of this step.      (no default)

        -L n,  --fragLen n                    n is fragment length for single-end reads.    (no default)        
        ------------------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines such as RASDA, MESDA and HISDA,
        please see the web site: https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies 
        and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.    
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------  
';


########## Version Infromation ##########
my $version_g = "  The Seventh Step of CISDA (ChIP-Seq Data Analyzer), version 0.62, 2016-01-13.";


########## Keys and Values ##########
if ($#ARGV   == -1) { print  "\n$HELP_g\n\n";  exit 0; }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0) { @ARGV = (@ARGV, "-h");           }       ## when the number of command argumants is odd. 
my %args = @ARGV;


########## Initialize  Variables ##########
my $input_g   = '7-FinalBAM/1-Subread';      ## This is only an initialization  value or suggesting value, not default value.
my $output_g  = '8-BigWig';                  ## This is only an initialization  value or suggesting value, not default value.
my $fragLen_g = 200;                         ## This is only an initialization  value or suggesting value, not default value.


########## Available Arguments ##########
my $available = "  -v  --version    -h  --help    -i  --input    -o    --output     -L  --fragLen     ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if($available !~ m/\s$key\s/) {print  "    Cann't recognize $key !!\n";  $boole_g = 1; }
}
if($boole_g == 1) {
    print "\n    The Command Line Arguments are wrong!\n";
    print   '    Please see help message by using "perl  CISDA-7.pl  -h". ';
    print "\n\n";
    exit 0;
}


########## Get Arguments ##########
if ( ( exists $args{'-v' } )  or  ( exists $args{'--version'      } )  )     { print  "\n$version_g\n\n";    exit 0; }
if ( ( exists $args{'-h' } )  or  ( exists $args{'--help'         } )  )     { print  "\n$HELP_g\n\n";       exit 0; }
if ( ( exists $args{'-i' } )  or  ( exists $args{'--input'        } )  )     { ($input_g  = $args{'-i' })   or  ($input_g  = $args{'--input'  });  }else{print   "\n -i or --input   is required.\n\n";   print  "\n$HELP_g\n\n";       exit 0; }                                               
if ( ( exists $args{'-o' } )  or  ( exists $args{'--output'       } )  )     { ($output_g = $args{'-o' })   or  ($output_g = $args{'--output' });  }else{print   "\n -o or --output  is required.\n\n";   print  "\n$HELP_g\n\n";       exit 0; }      
if ( ( exists $args{'-L' } )  or  ( exists $args{'--fragLen'      } )  )     { ($fragLen_g= $args{'-L' })   or  ($fragLen_g= $args{'--fragLen'});  }else{print   "\n -L or --fragLen is required.\n\n";   print  "\n$HELP_g\n\n";       exit 0; }      


########### Conditions #############
$input_g   =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";
$output_g  =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";
$fragLen_g =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";


######### Print Command Arguments to Standard Output ###########
print  "\n\n
        ################ Your Arguments ###############################
                Input  folder:  $input_g
                Output folder:  $output_g
                fragment length for single-end reads:  $fragLen_g
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

system("bamCoverage    --version   >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'     >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");
system("bamCorrelate    --version  >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'     >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");
system("bamFingerprint    --version  >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'       >> $output2_g/z-version_softwares.txt   2>&1");
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


my $BigWig1  = "$output_g/1-Extend-E7";
my $BigWig2  = "$output_g/2-Extend-1X";
my $BigWig3  = "$output_g/3-Extend-RPKM";

if ( !( -e $BigWig1) )   { mkdir $BigWig1  ||  die; }
if ( !( -e $BigWig2) )   { mkdir $BigWig2  ||  die; }
if ( !( -e $BigWig3) )   { mkdir $BigWig3  ||  die; }


my $EffectiveGenomeSize = 2150570000;
my $norFact = (10**7)*$fragLen_g; 
## There are three normalization methods: 10^7, 1X, RPKM.  The second method (1X) is the most reasonable method.

print "\n\nConvert BAM to BigWig ......\n";
my @BAMFiles  = '';
my @BAMFiles2 = '';
for(my $i=0; $i<=$#inputFiles; $i++) {
    next unless $inputFiles[$i] =~ m/\.bam$/;
    next unless $inputFiles[$i] !~ m/^[.]/;
    next unless $inputFiles[$i] !~ m/[~]$/;
    my $temp = $inputFiles[$i];
    $temp =~ s/\.bam$//  or  die;
    $BAMFiles[++$#BAMFiles] = "$input_g/$temp.bam";
    $BAMFiles2[++$#BAMFiles2] = "$temp";

    system(" bamCoverage   --numberOfProcessors 4    --bam $input_g/$temp.bam     --outFileName $BigWig1/$temp.bw     --outFileFormat bigwig       --normalizeTo1x  $norFact                --fragmentLength $fragLen_g      --binSize 1    --ignoreDuplicates  ");                
    system(" bamCoverage   --numberOfProcessors 4    --bam $input_g/$temp.bam     --outFileName $BigWig2/$temp.bw     --outFileFormat bigwig       --normalizeTo1x  $EffectiveGenomeSize    --fragmentLength $fragLen_g      --binSize 20   --ignoreDuplicates  ");                    
    system(" bamCoverage   --numberOfProcessors 4    --bam $input_g/$temp.bam     --outFileName $BigWig3/$temp.bw     --outFileFormat bigwig       --normalizeUsingRPKM                     --fragmentLength $fragLen_g      --binSize 20   --ignoreDuplicates  ");                                             
}


say(@BAMFiles); 
say(@BAMFiles2);

my $spearman   = "$output_g/spearman";
my $pearson    = "$output_g/pearson";
my $bamFinger  = "$output_g/bamFinger";

if ( !( -e $spearman)  )   { mkdir $spearman   ||  die; }
if ( !( -e $pearson)   )   { mkdir $pearson    ||  die; }
if ( !( -e $bamFinger) )   { mkdir $bamFinger  ||  die; }

system(" bamCorrelate  bins  --bamfiles @BAMFiles   -o $spearman/heatmap_spearman.svg   --corMethod spearman    --labels @BAMFiles2    --binSize 10000   --plotTitle Spearman_rank_correlation_coefficient  --numberOfProcessors 4   --outFileCorMatrix $spearman/spearman_matrix.txt       --plotNumbers   --fragmentLength $fragLen_g  --ignoreDuplicates  ");                                                
system(" bamCorrelate  bins  --bamfiles @BAMFiles   -o $pearson/heatmap_pearson.svg     --corMethod pearson     --labels @BAMFiles2    --binSize 10000   --plotTitle Pearson_rank_correlation_coefficient   --numberOfProcessors 4   --outFileCorMatrix $pearson/pearson_matrix.txt         --plotNumbers   --fragmentLength $fragLen_g  --ignoreDuplicates  ");    
system(" bamFingerprint      --bamfiles @BAMFiles   --plotFile  $bamFinger/Finger.svg    --binSize 500   --labels @BAMFiles2    --fragmentLength $fragLen_g   --numberOfSamples 1000000   --plotTitle bamFingerprint_Figure  --numberOfProcessors 4    --ignoreDuplicates");





print "\n\n\n\n\n##################################################################################################\n";
print "\n\n        Job Done! Cheers! \n\n";





























