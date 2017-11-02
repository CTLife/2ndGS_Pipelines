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
my $input_g  = '';  ## such as "9-BigWig/1_Bismark/1-allCs"
my $output_g = '';  ## such as "10-DeepTools/1_Bismark/1-allCs"

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use BSDA (BS-Seq Data Analyzer), version 0.9.0, 2017-10-01.
        BSDA is a Pipeline for Single-end and Paired-end BS-Seq Data Analysis by Integrating Lots of Softwares.
                                                            
        Step 9: Analyze all the bigwig files by using deepTools.  

        Usage:
               perl  BSDA9.pl    [-version]    [-help]   [-genome RefGenome]    [-in inputDir]    [-out outDir]
        For instance:
               perl  BSDA9.pl   -genome hg38   -in 9-BigWig/1_Bismark/1-allCs   -out 10-DeepTools/1_Bismark/1-allCs    > BSDA9.runLog

        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.

        -help           Show this help message and exit.

        Required arguments:
        -genome RefGenome   "RefGenome" is the short name of your reference genome, such as "mm10", "ce11", "hg38".    (no default)

        -in inputDir        "inputDir" is the name of input path that contains your bigwig files.  (no default)

        -out outDir         "outDir" is the name of output path that contains your running results of this step.  (no default)
        -----------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines, please visit https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ Jie Qiao Lab, yongp@outlook.com, Key Laboratory of Assisted Reproduction at Third Hospital,
        Academy for Advanced Interdisciplinary Studies, and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';

## Version Infromation
my $version = "    The 9th Step of BSDA (BS-Seq Data Analyzer), version 0.9.0, 2017-10-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$genome_g = 'hg38';                            ## This is only an initialization value or suggesting value, not default value.
$input_g  = '9-BigWig/1_Bismark/1-allCs';      ## This is only an initialization value or suggesting value, not default value.
$output_g = '10-DeepTools/1_Bismark/1-allCs';  ## This is only an initialization value or suggesting value, not default value.

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

$output_g = "$output_g/1-before-computeMatrix";
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

&printVersion("deeptools --version");

###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
say   "\n\n\n\n\n\n##################################################################################################";
say   "Checking all the input file names ......";
my @groupFiles = ();
my $fileNameBool = 1;
for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {
        next unless $inputFiles_g[$i] =~ m/\.bw$/;
        next unless $inputFiles_g[$i] !~ m/^[.]/;
        next unless $inputFiles_g[$i] !~ m/[~]$/;
        next unless $inputFiles_g[$i] !~ m/^QC_Results$/;
        next unless $inputFiles_g[$i] !~ m/^unpaired/;
        say   "\t......$inputFiles_g[$i]" ;
        my $temp = $inputFiles_g[$i];
        $groupFiles[++$#groupFiles] = $inputFiles_g[$i];
        $temp =~ m/^(\d+)_($pattern_g)_(Rep[1-9])/   or  die   "wrong-1: ## $temp ##";
        if($temp !~ m/^((\d+)_($pattern_g)_(Rep[1-9]))(.+)\.bw$/) {
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
say   "Detecting bigwig files in input folder ......";
my @BWfiles_g = ();
{
open(seqFiles_FH, ">", "$output2_g/BW-Files.txt")  or  die; 
for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {     
    next unless $inputFiles_g[$i] =~ m/\.bw$/;
    next unless $inputFiles_g[$i] !~ m/^[.]/;
    next unless $inputFiles_g[$i] !~ m/[~]$/;
    next unless $inputFiles_g[$i] !~ m/^unpaired/;
    say    "\t......$inputFiles_g[$i]"; 
    $inputFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])/  or  die;  
    $BWfiles_g[$#BWfiles_g+1] =  $inputFiles_g[$i];
    say        "\t\t\t\tBW file:  $inputFiles_g[$i]\n";
    say   seqFiles_FH  "BW file:  $inputFiles_g[$i]\n";
}

say   seqFiles_FH  "\n\n\n\n\n";  
say   seqFiles_FH  "All BW files:@BWfiles_g\n\n\n";
say        "\t\t\t\tAll BW files:@BWfiles_g\n\n";
my $num1 = $#BWfiles_g + 1;
say seqFiles_FH   "\nThere are $num1 BW files.\n";
say         "\t\t\t\tThere are $num1 BW files.\n";
}

my @BWfiles_g2 =  @BWfiles_g;
for ( my $i=0; $i<=$#BWfiles_g; $i++ ) { 
   $BWfiles_g[$i] = "$input_g/$BWfiles_g[$i]";
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
say   "\n\n\n\n\n\n##################################################################################################";
say   "Using deepTools ......";
    system("multiBigwigSummary bins   --bwfiles @BWfiles_g    --binSize 10000   --numberOfProcessors $numCores_g    --outRawCounts  $output_g/1-RawCounts.txt   --outFileName $output_g/1-results.npz    >> $output_g/1-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/1-results.npz    --whatToPlot heatmap        --corMethod pearson     -o $output_g/2-heatmap_pearson.pdf       --outFileCorMatrix $output_g/2-heatmap_pearson.txt         --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_g/2-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/1-results.npz    --whatToPlot heatmap        --corMethod spearman    -o $output_g/3-heatmap_spearman.pdf      --outFileCorMatrix $output_g/3-heatmap_spearman.txt        --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_g/3-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/1-results.npz    --whatToPlot scatterplot    --corMethod pearson     -o $output_g/4-scatterplot_pearson.pdf   --outFileCorMatrix $output_g/4-scatterplot_pearson.txt     --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_g/4-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/1-results.npz    --whatToPlot scatterplot    --corMethod spearman    -o $output_g/5-scatterplot_spearman.pdf  --outFileCorMatrix $output_g/5-scatterplot_spearman.txt    --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_g/5-runLog.txt   2>&1");                               

    system("multiBigwigSummary bins   --bwfiles @BWfiles_g    --binSize 5000   --numberOfProcessors $numCores_g    --outRawCounts  $output_g/6-RawCounts.txt   --outFileName $output_g/6-results.npz    >> $output_g/6-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/6-results.npz    --whatToPlot heatmap        --corMethod pearson     -o $output_g/7-heatmap_pearson.pdf        --outFileCorMatrix $output_g/7-heatmap_pearson.txt          --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_g/7-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/6-results.npz    --whatToPlot heatmap        --corMethod spearman    -o $output_g/8-heatmap_spearman.pdf       --outFileCorMatrix $output_g/8-heatmap_spearman.txt         --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_g/8-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/6-results.npz    --whatToPlot scatterplot    --corMethod pearson     -o $output_g/9-scatterplot_pearson.pdf    --outFileCorMatrix $output_g/9-scatterplot_pearson.txt      --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_g/9-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/6-results.npz    --whatToPlot scatterplot    --corMethod spearman    -o $output_g/10-scatterplot_spearman.pdf  --outFileCorMatrix $output_g/10-scatterplot_spearman.txt    --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_g/10-runLog.txt  2>&1");                               

    system("multiBigwigSummary bins   --bwfiles @BWfiles_g    --binSize 1000   --numberOfProcessors $numCores_g    --outRawCounts  $output_g/11-RawCounts.txt   --outFileName $output_g/11-results.npz    >> $output_g/11-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/11-results.npz    --whatToPlot heatmap        --corMethod pearson     -o $output_g/12-heatmap_pearson.pdf        --outFileCorMatrix $output_g/12-heatmap_pearson.txt          --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_g/12-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/11-results.npz    --whatToPlot heatmap        --corMethod spearman    -o $output_g/13-heatmap_spearman.pdf       --outFileCorMatrix $output_g/13-heatmap_spearman.txt         --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_g/13-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/11-results.npz    --whatToPlot scatterplot    --corMethod pearson     -o $output_g/14-scatterplot_pearson.pdf    --outFileCorMatrix $output_g/14-scatterplot_pearson.txt      --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_g/14-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/11-results.npz    --whatToPlot scatterplot    --corMethod spearman    -o $output_g/15-scatterplot_spearman.pdf   --outFileCorMatrix $output_g/15-scatterplot_spearman.txt     --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_g/15-runLog.txt  2>&1");                               
}
###################################################################################################################################################################################################    





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";


 

  
## END
