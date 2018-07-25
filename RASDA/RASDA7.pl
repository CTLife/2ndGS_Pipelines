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
my $input_g  = '';  ## such as "7-BigWig/3A_STAR"
my $output_g = '';  ## such as "8-DeepTools/3A_STAR"

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use RASDA (RNA-Seq Data Analyzer), version 0.9.4,  2018-02-01.
        RASDA is a Pipeline for Single-end and Paired-end RNA-Seq Data Analysis by Integrating Lots of Softwares.
                                                            
        Step 7: Analyze all the bigwig and BAM files by using deepTools.  

        Usage:
               perl  RASDA7.pl    [-version]    [-help]   [-genome RefGenome]    [-in inputDir]    [-out outDir]
        For instance:
               perl  RASDA7.pl   -genome hg38   -in 7-BigWig/3A_STAR   -out 8-DeepTools/3A_STAR    > RASDA7.runLog

        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.

        -help           Show this help message and exit.

        Required arguments:
        -genome RefGenome   "RefGenome" is the short name of your reference genome, such as "mm10", "ce11", "hg38".    (no default)

        -in inputDir        "inputDir" is the name of input path that contains your BigWig files.  (no default)

        -out outDir         "outDir" is the name of output path that contains your running results of this step.  (no default)
        -----------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines, please visit https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ Jie Qiao Lab, yongp@outlook.com, Key Laboratory of Assisted Reproduction at Third Hospital,
        Academy for Advanced Interdisciplinary Studies, and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';

## Version Infromation
my $version = "    The Seventh Step of RASDA (RNA-Seq Data Analyzer), version 0.9.4,  2018-02-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$genome_g = 'hg38';                 ## This is only an initialization value or suggesting value, not default value.
$input_g  = '7-BigWig/3A_STAR';     ## This is only an initialization value or suggesting value, not default value.
$output_g = '8-DeepTools/3A_STAR';  ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help   -genome   -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  RASDA7.pl  -help' \n";
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
        $temp =~ m/_(Rep[1-9])\.bw$/   or  die   "wrong-2: ## $temp ##";
        if($temp !~ m/^((\d+)_($pattern_g)_(Rep[1-9]))(_[1-2])?\.bw$/) {
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
say   "Detecting BigWig files in input folder ......";
my @BigWigfiles_g = ();
{
open(seqFiles_FH, ">", "$output2_g/BigWig-Files.txt")  or  die; 
for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {     
    next unless $inputFiles_g[$i] =~ m/\.bw$/;
    next unless $inputFiles_g[$i] !~ m/^[.]/;
    next unless $inputFiles_g[$i] !~ m/[~]$/;
    next unless $inputFiles_g[$i] !~ m/^unpaired/;
    say    "\t......$inputFiles_g[$i]"; 
    $inputFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])\.bw$/  or  die;  
    $BigWigfiles_g[$#BigWigfiles_g+1] =  $inputFiles_g[$i];
    say        "\t\t\t\tBigWig file:  $inputFiles_g[$i]\n";
    say   seqFiles_FH  "BigWig file:  $inputFiles_g[$i]\n";
}

say   seqFiles_FH  "\n\n\n\n\n";  
say   seqFiles_FH  "All BigWig files:@BigWigfiles_g\n\n\n";
say        "\t\t\t\tAll BigWig files:@BigWigfiles_g\n\n";
my $num1 = $#BigWigfiles_g + 1;
say seqFiles_FH   "\nThere are $num1 BigWig files.\n";
say         "\t\t\t\tThere are $num1 BigWig files.\n";
}

my @BigWigfiles_g2 =  @BigWigfiles_g;    
my @BigWigfiles_g3 =  @BigWigfiles_g;    ## for BAM files
for ( my $i=0; $i<=$#BigWigfiles_g2; $i++ ) { 
   $BigWigfiles_g2[$i] = "$input_g/$BigWigfiles_g2[$i]";   ## add path
   $BigWigfiles_g3[$i] = "$input_g/$BigWigfiles_g3[$i]";   ## add path
   $BigWigfiles_g3[$i] =~ s/7-BigWig/5-finalBAM/ or  die;  
   $BigWigfiles_g3[$i] =~ s/\.bw$/.bam/ or  die;  
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Using multiBigwigSummary bin, plotCorrelation and plotPCA ......";
my $output_sub1 = "$output_g/1_multiBigwigSummary_bin";   
&myMakeDir($output_sub1);
{
    system("multiBigwigSummary bins   --bwfiles @BigWigfiles_g2    --smartLabels    --binSize 10000   --numberOfProcessors $numCores_g    --outRawCounts  $output_sub1/1-RawCounts.10000bpBin.txt   --outFileName $output_sub1/1-results.npz    >> $output_sub1/1-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/1-results.npz    --whatToPlot heatmap        --corMethod pearson     -o $output_sub1/2-heatmap_pearson.pdf       --outFileCorMatrix $output_sub1/2-heatmap_pearson.txt         --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/2-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/1-results.npz    --whatToPlot heatmap        --corMethod spearman    -o $output_sub1/3-heatmap_spearman.pdf      --outFileCorMatrix $output_sub1/3-heatmap_spearman.txt        --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/3-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/1-results.npz    --whatToPlot scatterplot    --corMethod pearson     -o $output_sub1/4-scatterplot_pearson.pdf   --outFileCorMatrix $output_sub1/4-scatterplot_pearson.txt     --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/4-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/1-results.npz    --whatToPlot scatterplot    --corMethod spearman    -o $output_sub1/5-scatterplot_spearman.pdf  --outFileCorMatrix $output_sub1/5-scatterplot_spearman.txt    --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/5-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/1-results.npz    -o $output_sub1/6-plotPCA.pdf  --outFileNameData $output_sub1/6-plotPCA.txt    --plotHeight 20   --plotWidth 20      >> $output_sub1/6-plotPCA-runLog.txt   2>&1");                               

    system("multiBigwigSummary bins   --bwfiles @BigWigfiles_g2    --smartLabels    --binSize 5000   --numberOfProcessors $numCores_g    --outRawCounts  $output_sub1/6-RawCounts.5000bpBin.txt   --outFileName $output_sub1/6-results.npz    >> $output_sub1/6-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/6-results.npz    --whatToPlot heatmap        --corMethod pearson     -o $output_sub1/7-heatmap_pearson.pdf        --outFileCorMatrix $output_sub1/7-heatmap_pearson.txt          --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_sub1/7-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/6-results.npz    --whatToPlot heatmap        --corMethod spearman    -o $output_sub1/8-heatmap_spearman.pdf       --outFileCorMatrix $output_sub1/8-heatmap_spearman.txt         --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_sub1/8-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/6-results.npz    --whatToPlot scatterplot    --corMethod pearson     -o $output_sub1/9-scatterplot_pearson.pdf    --outFileCorMatrix $output_sub1/9-scatterplot_pearson.txt      --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_sub1/9-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/6-results.npz    --whatToPlot scatterplot    --corMethod spearman    -o $output_sub1/10-scatterplot_spearman.pdf  --outFileCorMatrix $output_sub1/10-scatterplot_spearman.txt    --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_sub1/10-runLog.txt  2>&1");                               
    system("plotPCA         -in $output_sub1/6-results.npz    -o $output_sub1/11-plotPCA.pdf  --outFileNameData $output_sub1/11-plotPCA.txt    --plotHeight 20   --plotWidth 20      >> $output_sub1/11-plotPCA-runLog.txt   2>&1");                               

    system("multiBigwigSummary bins   --bwfiles @BigWigfiles_g2    --smartLabels    --binSize 2000   --numberOfProcessors $numCores_g    --outRawCounts  $output_sub1/11-RawCounts.2000bpBin.txt   --outFileName $output_sub1/11-results.npz    >> $output_sub1/11-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/11-results.npz    --whatToPlot heatmap        --corMethod pearson     -o $output_sub1/12-heatmap_pearson.pdf        --outFileCorMatrix $output_sub1/12-heatmap_pearson.txt          --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_sub1/12-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/11-results.npz    --whatToPlot heatmap        --corMethod spearman    -o $output_sub1/13-heatmap_spearman.pdf       --outFileCorMatrix $output_sub1/13-heatmap_spearman.txt         --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_sub1/13-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/11-results.npz    --whatToPlot scatterplot    --corMethod pearson     -o $output_sub1/14-scatterplot_pearson.pdf    --outFileCorMatrix $output_sub1/14-scatterplot_pearson.txt      --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_sub1/14-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/11-results.npz    --whatToPlot scatterplot    --corMethod spearman    -o $output_sub1/15-scatterplot_spearman.pdf   --outFileCorMatrix $output_sub1/15-scatterplot_spearman.txt     --plotHeight 20   --plotWidth 20   --skipZeros  --removeOutliers     >> $output_sub1/15-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/11-results.npz    -o $output_sub1/17-plotPCA.pdf  --outFileNameData $output_sub1/17-plotPCA.txt    --plotHeight 20   --plotWidth 20      >> $output_sub1/17-plotPCA-runLog.txt   2>&1");                               

}
###################################################################################################################################################################################################    





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Using bigwigCompare ......";
my $output_sub2 = "$output_g/2_bigwigCompare";   
&myMakeDir($output_sub2);
{

for (my $i=0; $i<=$#BigWigfiles_g-1;  $i++) {
    for (my $j=$i+1; $j<=$#BigWigfiles_g;  $j++) {
        my $file1 = $BigWigfiles_g[$i];   
        my $file2 = $BigWigfiles_g[$j];
        system( "bigwigCompare  -b1 $input_g/$file1   -b2 $input_g/$file2  --operation log2      --binSize 20  --numberOfProcessors $numCores_g  -o $output_sub2/$file1....$file2.log2.bw          >> $output_sub2/1-runLog.txt   2>&1  " );   
        system( "bigwigCompare  -b1 $input_g/$file1   -b2 $input_g/$file2  --operation subtract  --binSize 20  --numberOfProcessors $numCores_g  -o $output_sub2/$file1....$file2.subtract.bw      >> $output_sub2/2-runLog.txt   2>&1  " ); 
   }
}

}
###################################################################################################################################################################################################    





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Using bamCoverage to generate bigwig files ......";
my $genome_size_g = 0;
if($genome_g eq 'hg38') {$genome_size_g = 2451960000; }
if($genome_g eq 'mm10') {$genome_size_g = 2150570000; }
if($genome_g eq 'dm6')  {$genome_size_g = 121400000;  }
print("\n\ngenome_size: $genome_g, $genome_size_g\n\n");
my  $input_g2 = $input_g; 
$input_g2 =~ s/7-BigWig/5-finalBAM/ or die;
my $output_sub3 = "$output_g/3_bigwig_byDeepTools";   
&myMakeDir($output_sub3);
{

for (my $i=0; $i<=$#BigWigfiles_g-1;  $i++) {
    my $temp = $BigWigfiles_g[$i]; 
    $temp =~ s/\.bw$//  ||  die;    
    system("bamCoverage     --bam $input_g2/$temp.bam  --binSize 20   --normalizeUsing RPKM   --effectiveGenomeSize $genome_size_g      --numberOfProcessors $numCores_g    --outFileName  $output_sub3/$temp.RPKM.bw          >> $output_sub3/1-runLog.txt  2>&1 ");
    system("bamCoverage     --bam $input_g2/$temp.bam  --binSize 20   --normalizeUsing BPM    --effectiveGenomeSize $genome_size_g      --numberOfProcessors $numCores_g    --outFileName  $output_sub3/$temp.BPM.bw           >> $output_sub3/2-runLog.txt  2>&1 ");
}

}
###################################################################################################################################################################################################    





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Using some other commands in deepTools ......";
my $output_sub4 = "$output_g/4_DeepTools_Other";   
&myMakeDir($output_sub4);

system("plotCoverage       --bamfiles   @BigWigfiles_g3    --plotFile    $output_sub4/1-Coverage-check.pdf   --outRawCounts $output_sub4/1-RawCounts.txt    --plotHeight 12   --plotWidth 30  --numberOfProcessors $numCores_g      >> $output_sub4/1-runLog.txt  2>&1");   
system("bamPEFragmentSize  --bamfiles   @BigWigfiles_g3    --histogram   $output_sub4/2-FragmentSize.pdf    --plotFileFormat pdf  --numberOfProcessors $numCores_g    >> $output_sub4/2-runLog.txt  2>&1");  

for (my $i=0; $i<=$#BigWigfiles_g; $i++) {
    my $temp = $BigWigfiles_g[$i]; 
    $temp =~ s/\.bw$//  ||  die; 
    say   "\t......$BigWigfiles_g[$i]";
    system("computeGCBias   --bamfile $input_g2/$temp.bam     --effectiveGenomeSize $genome_size_g     --genome /media/yp/biox1/.RefGenomes/3-FASTA-2bit/$genome_g.2bit    --numberOfProcessors $numCores_g    --GCbiasFrequenciesFile  $output_sub4/$temp.GCBias.txt    --biasPlot  $output_sub4/$temp.GCBias.png    --fragmentLength 220        >> $output_sub4/$temp.GCBias.runLog.txt  2>&1 ");                 
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";


 

  
## END






