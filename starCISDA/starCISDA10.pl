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
my $output_g = '';  ## such as "11-DeepTools/1_Trim_BWAmem"

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use starCISDA (small-scale TELP-assisted rapid ChIP–seq Data Analyzer), version 0.9.0, 2017-10-01.
        starCISDA is a Pipeline for Single-end and Paired-end small-scale TELP-assisted rapid ChIP–seq Data Analysis by Integrating Lots of Softwares.
                                                            
        Step 10: Analyze all the BAM files by using deepTools.  

        Usage:
               perl  starCISDA10.pl    [-version]    [-help]   [-genome RefGenome]    [-in inputDir]    [-out outDir]
        For instance:
               perl  starCISDA10.pl   -genome hg38   -in 7-finalBAM/1_Trim_BWAmem   -out 11-DeepTools/1_Trim_BWAmem    > starCISDA10.runLog

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
my $version = "    The Tenth Step of starCISDA (small-scale TELP-assisted rapid ChIP–seq Data Analyzer), version 0.9.0, 2017-10-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$genome_g = 'hg38';                       ## This is only an initialization value or suggesting value, not default value.
$input_g  = '7-finalBAM/1_Trim_BWAmem';   ## This is only an initialization value or suggesting value, not default value.
$output_g = '11-DeepTools/1_Trim_BWAmem'; ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help   -genome   -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  starCISDA10.pl  -help' \n";
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

my @BAMfiles_g2 =  @BAMfiles_g;
for ( my $i=0; $i<=$#BAMfiles_g; $i++ ) { 
   $BAMfiles_g[$i] = "$input_g/$BAMfiles_g[$i]";
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
say   "\n\n\n\n\n\n##################################################################################################";
say   "Using deepTools ......";
    system("multiBamSummary bins   --bamfiles @BAMfiles_g    --binSize 10000   --numberOfProcessors $numCores_g    --outRawCounts  $output_g/1-RawCounts.txt   --outFileName $output_g/1-results.npz    >> $output_g/1-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/1-results.npz    --whatToPlot heatmap        --corMethod pearson     -o $output_g/2-heatmap_pearson.pdf       --outFileCorMatrix $output_g/2-heatmap_pearson.txt         --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_g/2-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/1-results.npz    --whatToPlot heatmap        --corMethod spearman    -o $output_g/3-heatmap_spearman.pdf      --outFileCorMatrix $output_g/3-heatmap_spearman.txt        --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_g/3-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/1-results.npz    --whatToPlot scatterplot    --corMethod pearson     -o $output_g/4-scatterplot_pearson.pdf   --outFileCorMatrix $output_g/4-scatterplot_pearson.txt     --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_g/4-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/1-results.npz    --whatToPlot scatterplot    --corMethod spearman    -o $output_g/5-scatterplot_spearman.pdf  --outFileCorMatrix $output_g/5-scatterplot_spearman.txt    --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_g/5-runLog.txt   2>&1");                               

    system("multiBamSummary bins   --bamfiles @BAMfiles_g    --binSize 5000   --numberOfProcessors $numCores_g    --outRawCounts  $output_g/6-RawCounts.txt   --outFileName $output_g/6-results.npz    >> $output_g/6-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/6-results.npz    --whatToPlot heatmap        --corMethod pearson     -o $output_g/7-heatmap_pearson.pdf        --outFileCorMatrix $output_g/7-heatmap_pearson.txt          --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_g/7-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/6-results.npz    --whatToPlot heatmap        --corMethod spearman    -o $output_g/8-heatmap_spearman.pdf       --outFileCorMatrix $output_g/8-heatmap_spearman.txt         --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_g/8-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/6-results.npz    --whatToPlot scatterplot    --corMethod pearson     -o $output_g/9-scatterplot_pearson.pdf    --outFileCorMatrix $output_g/9-scatterplot_pearson.txt      --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_g/9-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/6-results.npz    --whatToPlot scatterplot    --corMethod spearman    -o $output_g/10-scatterplot_spearman.pdf  --outFileCorMatrix $output_g/10-scatterplot_spearman.txt    --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_g/10-runLog.txt  2>&1");                               

    system("multiBamSummary bins   --bamfiles @BAMfiles_g    --binSize 1000   --numberOfProcessors $numCores_g    --outRawCounts  $output_g/11-RawCounts.txt   --outFileName $output_g/11-results.npz    >> $output_g/11-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/11-results.npz    --whatToPlot heatmap        --corMethod pearson     -o $output_g/12-heatmap_pearson.pdf        --outFileCorMatrix $output_g/12-heatmap_pearson.txt          --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_g/12-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/11-results.npz    --whatToPlot heatmap        --corMethod spearman    -o $output_g/13-heatmap_spearman.pdf       --outFileCorMatrix $output_g/13-heatmap_spearman.txt         --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_g/13-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/11-results.npz    --whatToPlot scatterplot    --corMethod pearson     -o $output_g/14-scatterplot_pearson.pdf    --outFileCorMatrix $output_g/14-scatterplot_pearson.txt      --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_g/14-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_g/11-results.npz    --whatToPlot scatterplot    --corMethod spearman    -o $output_g/15-scatterplot_spearman.pdf   --outFileCorMatrix $output_g/15-scatterplot_spearman.txt     --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_g/15-runLog.txt  2>&1");                               

    system("plotCoverage  --bamfiles   @BAMfiles_g    --plotFile $output_g/16-Coverage-check.pdf   --outRawCounts $output_g/16-RawCounts.txt    --plotHeight 20   --plotWidth 20  --numberOfProcessors $numCores_g    --ignoreDuplicates   >> $output_g/16-runLog.txt  2>&1");   
    system("bamPEFragmentSize  --bamfiles   @BAMfiles_g    --histogram   $output_g/17-FragmentSize.png   --numberOfProcessors $numCores_g    >> $output_g/17-runLog.txt  2>&1");  
}
###################################################################################################################################################################################################    





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Using deepTools for more ......";

my $genome_size_g = 0;
if($genome_g eq 'hg38') {$genome_size_g = 2451960000; }
if($genome_g eq 'mm10') {$genome_size_g = 2150570000; }
if($genome_g eq 'dm6')  {$genome_size_g = 121400000;  }
print("\n\ngenome_size: $genome_g, $genome_size_g\n\n");

for (my $i=0; $i<=$#BAMfiles_g2; $i++) {
    my $temp = $BAMfiles_g2[$i]; 
    $temp =~ s/\.bam$//  ||  die; 
    say   "\t......$BAMfiles_g[$i]";
    &myMakeDir("$output_g/$temp");
    &myMakeDir("$output_g/bigwig");
    system("computeGCBias   --bamfile $input_g/$temp.bam     --effectiveGenomeSize $genome_size_g     --genome /media/yp/biox1/.RefGenomes/$genome_g.2bit    --numberOfProcessors $numCores_g    --GCbiasFrequenciesFile  $output_g/$temp/1-GCBias.txt    --biasPlot  $output_g/$temp/1-GCBias.png    --fragmentLength 220        >> $output_g/$temp/1-runLog.txt  2>&1 ");                 
    system("bamCoverage     --bam $input_g/$temp.bam     --normalizeTo1x $genome_size_g      --numberOfProcessors $numCores_g    --outFileName  $output_g/bigwig/$temp.1X.bw      --extendReads   --ignoreDuplicates       >> $output_g/$temp/2-runLog.txt  2>&1 ");
    system("bamCoverage     --bam $input_g/$temp.bam     --normalizeUsingRPKM                --numberOfProcessors $numCores_g    --outFileName  $output_g/bigwig/$temp.RPKM.bw    --extendReads   --ignoreDuplicates       >> $output_g/$temp/3-runLog.txt  2>&1 ");
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";


 

  
## END
