#!/usr/bin/env  perl
use  strict;
use  warnings;
use  v5.12;  ## Perl5 version >= 5.12
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $genome_g   = '';  ## such as "mm39", "hg38".
my $input_g    = '';  ## such as "2_finalFASTQ"
my $output_g   = '';  ## such as "3_rawBAM"

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Step 3: Mapping reads to the reference genome by using BWA aln or Bowtie2.
                And assess the quality of BAM files to identify possible mapping errors or biases by using 12 softwares:
                SAMtools, Subread utilities, FASTQC, SAMstat, qualimap, deepTools, PRESEQ, Picard, goleft, Bamtools, Sambamba, and MultiQC.

                If this script works well, you do not need to check the the versions of the softwares or packages which are used in this script. 
                And you do not need to exactly match the versions of the softwares or packages.
                If some errors or warnings are reported, please check the versions of softwares or packages.

                The versions of softwares or packages are used in this script (This is shown for reference only. And do not need to exactly match them.):  
                        bwa,       0.7.17-r1198-dirty  (bwa)
                        bowtie2,   2.4.5               (bowtie2  --version)
                        SAMtools,  1.15.1              (samtools)   
                        Subread,   2.0.3               (propmapped)
                        FASTQC,    0.11.9              (fastqc  -v)    
                        SAMstat,   1.5.1               (samstat   -v)      
                        qualimap,  v.2.2.2-dev         (qualimap  -v) 
                        deepTools, 3.5.1               (plotFingerprint --version)
                        PRESEQ,    3.1.1               (preseq)    
                        Picard,    2.21.6              (which  picard.jar)
                        goleft,    0.2.4               (goleft  -v)
                        Bamtools,  2.5.2               (bamtools -v)  
                        sambamba,  0.8.2               (sambamba)
                        MultiQC,   1.12                (multiqc   --version)       

        Usage:
               perl Bulk_5hmCSeal_3.pl    [-version]    [-help]   [-genome RefGenome]    [-in inputDir]    [-out outDir]  
        For instance:
               nohup time  perl Bulk_5hmCSeal_3.pl   -genome hg38   -in 2_finalFASTQ   -out 3_rawBAM     > Bulk_5hmCSeal_3.runLog.txt  2>&1    &  

        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.

        -help           Show this help message and exit.

        Required arguments:
        -genome RefGenome   "RefGenome" is the short name of your reference genome, such as "mm39", "ce11", "hg38".    (no default)

        -in inputDir        "inputDir" is the name of input path that contains your FASTQ files.  (no default)

        -out outDir         "outDir" is the name of output path that contains your running results (BAM files) of this step.  (no default)
        -----------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines, please visit https://github.com/CTLife/2ndGS_Pipelines
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';

## Version Infromation
my $version = "    The 3rd Step, version 1.0,  2022-09-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$genome_g = 'hg38';                    ## This is only an initialization value or suggesting value, not default value.
$input_g  = '2_finalFASTQ';            ## This is only an initialization value or suggesting value, not default value.
$output_g = '3_rawBAM';                ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help   -genome   -in   -out   ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  Bulk_5hmCSeal_3.pl  -help' \n";
    exit 0;
}

## Get Arguments
if ( exists $args{'-version' }   )     { say  "\n$version\n";    exit 0; }
if ( exists $args{'-help'    }   )     { say  "\n$HELP\n";       exit 0; }
if ( exists $args{'-genome'  }   )     { $genome_g = $args{'-genome'  }; }else{say   "\n -genome is required.\n";   say  "\n$HELP\n";    exit 0; }
if ( exists $args{'-in'      }   )     { $input_g  = $args{'-in'      }; }else{say   "\n -in     is required.\n";   say  "\n$HELP\n";    exit 0; }
if ( exists $args{'-out'     }   )     { $output_g = $args{'-out'     }; }else{say   "\n -out    is required.\n";   say  "\n$HELP\n";    exit 0; }

## Conditions
$genome_g =~ m/^\S+$/       ||  die   "\n\n$HELP\n\n";
$input_g  =~ m/^\S+$/       ||  die   "\n\n$HELP\n\n";
$output_g =~ m/^\S+$/       ||  die   "\n\n$HELP\n\n";

## Print Command Arguments to Standard Output
say  "\n
        ################ Arguments ###############################
                Reference Genome:  $genome_g
                Input       Path:  $input_g
                Output      Path:  $output_g
        ##########################################################
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
## Context specific:
my  $commonPath_g     = "/media/yp/yongpeng-2TB/MyProgramFiles/2_Aligners";
my  $BWA_index_g      = "";
my  $Bowtie2_index_g  = "";

if($genome_g eq "hg38") {
    $BWA_index_g       = "$commonPath_g/bwa/RefGenomes/hg38/hg38";  
    $Bowtie2_index_g   = "$commonPath_g/bowtie2/RefGenomes/hg38/hg38";
}

if($genome_g eq "mm39") {
    $BWA_index_g       = "$commonPath_g/bwa/RefGenomes/mm39/mm39";  
    $Bowtie2_index_g   = "$commonPath_g/bowtie2/RefGenomes/mm39/mm39";
}

say   "\n\n\n\n\n\n##################################################################################################";
say   $BWA_index_g;
say   $Bowtie2_index_g;
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

my  $Picard_g   = &fullPathApp("picard.jar");
&printVersion("bwa ");
&printVersion("bowtie2    --version");
&printVersion("samtools");
&printVersion("fastqc    -v");
&printVersion("samstat   -v");
&printVersion("preseq");
&printVersion("qualimap  -v");
&printVersion("multiqc   --version");
&printVersion("propmapped");     ## in subread
&printVersion("qualityScores");  ## in subread
&printVersion("plotFingerprint --version");
&printVersion("goleft  -v");
&printVersion("bamtools -v");  
&printVersion("sambamba"); 

&printVersion("java  -jar  $Picard_g   CollectIndependentReplicateMetrics  --version");
&printVersion("java  -jar  $Picard_g   CollectAlignmentSummaryMetrics      --version");
&printVersion("java  -jar  $Picard_g   CollectBaseDistributionByCycle      --version");
&printVersion("java  -jar  $Picard_g   CollectGcBiasMetrics                --version");
&printVersion("java  -jar  $Picard_g   CollectInsertSizeMetrics            --version");
&printVersion("java  -jar  $Picard_g   CollectJumpingLibraryMetrics        --version");
&printVersion("java  -jar  $Picard_g   CollectMultipleMetrics              --version");
&printVersion("java  -jar  $Picard_g   CollectOxoGMetrics                  --version");
&printVersion("java  -jar  $Picard_g   CollectQualityYieldMetrics          --version");
&printVersion("java  -jar  $Picard_g   CollectSequencingArtifactMetrics    --version");
&printVersion("java  -jar  $Picard_g   CollectTargetedPcrMetrics           --version");
&printVersion("java  -jar  $Picard_g   CollectWgsMetrics                   --version");
&printVersion("java  -jar  $Picard_g   EstimateLibraryComplexity           --version");
&printVersion("java  -jar  $Picard_g   MeanQualityByCycle                  --version");
&printVersion("java  -jar  $Picard_g   QualityScoreDistribution            --version");
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Detecting single-end and paired-end FASTQ files ......";
my @singleEnd_g   = ();
my @pairedEnd_g   = ();
open(seqFiles_FH_g, ">", "$output2_g/singleEnd-pairedEnd-Files.txt")  or  die;

for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {
    next unless $inputFiles_g[$i] !~ m/^[.]/;
    next unless $inputFiles_g[$i] !~ m/[~]$/;
    next unless $inputFiles_g[$i] !~ m/\.sh$/;
    next unless $inputFiles_g[$i] !~ m/\.txt$/;
    next unless $inputFiles_g[$i] !~ m/^unpaired/;
    next unless $inputFiles_g[$i] !~ m/^QC_Results$/;
    next unless $inputFiles_g[$i] !~ m/^Trimmomatic$/;
    next unless $inputFiles_g[$i] =~ m/\.fastq/;
    next unless $inputFiles_g[$i] =~ m/\.fastq\.gz$/;
    my $temp =  $inputFiles_g[$i]; 
    say    "\t......$temp";
    if ($temp !~ m/^(\S+)\.R[12]\.fastq/) {   ## sinlge end sequencing files.
        $temp =~ m/^(\S+)\.fastq/   or  die;
        $singleEnd_g[$#singleEnd_g+1] =  $temp;
        say         "\t\t\t\tSingle-end sequencing files: $temp\n";
        say  seqFiles_FH_g  "Single-end sequencing files: $temp\n";
    }else{     ## paired end sequencing files.
        $temp =~ m/^(\S+)\.R[12]\.fastq/  or  die;
        $pairedEnd_g[$#pairedEnd_g+1] =  $temp;
        say        "\t\t\t\tPaired-end sequencing files: $temp\n";
        say seqFiles_FH_g  "Paired-end sequencing files: $temp\n";
    }
}

@singleEnd_g  = sort(@singleEnd_g);
@pairedEnd_g  = sort(@pairedEnd_g);

( ($#pairedEnd_g+1)%2 == 0 )  or die;
say   seqFiles_FH_g  "\n\n\n\n\n";
say   seqFiles_FH_g  "All single-end sequencing files:@singleEnd_g\n\n\n";
say   seqFiles_FH_g  "All paired-end sequencing files:@pairedEnd_g\n\n\n";
say          "\t\t\t\tAll single-end sequencing files:@singleEnd_g\n\n";
say          "\t\t\t\tAll paired-end sequencing files:@pairedEnd_g\n\n";

my $numSingle_g = $#singleEnd_g + 1;
my $numPaired_g = $#pairedEnd_g + 1;

say seqFiles_FH_g   "\nThere are $numSingle_g single-end sequencing files.\n";
say seqFiles_FH_g   "\nThere are $numPaired_g paired-end sequencing files.\n";
say           "\t\t\t\tThere are $numSingle_g single-end sequencing files.\n";
say           "\t\t\t\tThere are $numPaired_g paired-end sequencing files.\n";


for ( my $i=0; $i<$#pairedEnd_g; $i=$i+2 ) {
    my $temp = $pairedEnd_g[$i]; 
    $temp =~ s/\.R1\.fastq/.R2.fastq/  or die "\n##Error-1: $temp ##\n\n";
    ($pairedEnd_g[$i+1] eq $temp) or die "\n##Error-2: $temp ## $pairedEnd_g[$i+1] ##\n\n";
}

print("\n\n\n\n\n#########################################\n");
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_BAM_1  {
    my $dir1      =  $_[0];   ## All the SAM files must be in this folder.
    opendir(my $FH_Files, $dir1) || die;
    my @Files = readdir($FH_Files);

    my $QCresults = "$dir1/QC_Results";
    my $SAMtools  = "$QCresults/1_SAMtools";
    my $FastQC    = "$QCresults/2_FastQC";
    my $qualimap  = "$QCresults/3_qualimap";
    my $samstat   = "$QCresults/4_samstat";
    my $Bamtools  = "$QCresults/5_Bamtools";
    my $MultiQC1  = "$QCresults/6_MultiQC1_FastQC";
    my $MultiQC2  = "$QCresults/6_MultiQC2_qualimap";
    my $MultiQC3  = "$QCresults/6_MultiQC3_SAMtools";
    my $MultiQC4  = "$QCresults/6_MultiQC4_Bamtools";
    my $MultiQC5  = "$QCresults/6_MultiQC5_Aligner";
    &myMakeDir($QCresults);
    &myMakeDir($SAMtools);
    &myMakeDir($FastQC);
    &myMakeDir($qualimap);
    &myMakeDir($samstat);
    &myMakeDir($Bamtools);
    &myMakeDir($MultiQC1);
    &myMakeDir($MultiQC2);
    &myMakeDir($MultiQC3);
    &myMakeDir($MultiQC4);
    &myMakeDir($MultiQC5);

    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of all SAM files by using SAMtools, FastQC, qualimap, samstat, Bamtools and MultiQC ......";
    for ( my $i=0; $i<=$#Files; $i++ ) {
        next unless $Files[$i] =~ m/\.sam$/;
        next unless $Files[$i] !~ m/^[.]/;
        next unless $Files[$i] !~ m/[~]$/;
        my $temp = $Files[$i];
        say    "\t......$temp";
        $temp =~ s/\.sam$//  ||  die;
        system( "samtools  sort  -o $dir1/$temp.bam   --output-fmt bam  -T $dir1/Temp...$temp   --threads $numCores_g    $dir1/$temp.sam    >> $SAMtools/$temp.sort.runLog    2>&1");
        system( "samtools  index           $dir1/$temp.bam      >> $SAMtools/$temp.index.runLog         2>&1");
        system( "samtools  flagstat        $dir1/$temp.bam      >> $SAMtools/$temp.flagstat.runLog      2>&1");
        system( "samtools  idxstats        $dir1/$temp.bam      >> $SAMtools/$temp.idxstats.runLog      2>&1");
        system( "fastqc    --outdir $FastQC    --threads $numCores_g  --format bam   --kmers 7    $dir1/$temp.bam                   >> $FastQC/$temp.runLog      2>&1" );
        system( "qualimap  bamqc  -bam $dir1/$temp.bam   -c  -ip  -nt $numCores_g   -outdir $qualimap/$temp   --java-mem-size=16G   >> $qualimap/$temp.runLog    2>&1" );
        system( "samstat   $dir1/$temp.bam      >> $samstat/$temp.runLog         2>&1");
        system( "rm   $dir1/$temp.sam" );
        system( "bamtools   count    -in  $dir1/$temp.bam      > $Bamtools/bamtools_count.$temp.txt  ");
        system( "bamtools   stats    -in  $dir1/$temp.bam      > $Bamtools/bamtools_stats.$temp.txt  ");   
    }

    system( "multiqc    --title FastQC     --verbose  --export   --outdir $MultiQC1   --pdf         $FastQC            >> $MultiQC1/MultiQC.FastQC.runLog     2>&1" );
    system( "multiqc    --title qualimap   --verbose  --export   --outdir $MultiQC2   --pdf         $qualimap          >> $MultiQC2/MultiQC.qualimap.runLog   2>&1" );
    system( "multiqc    --title SAMtools   --verbose  --export   --outdir $MultiQC3   --pdf         $SAMtools          >> $MultiQC3/MultiQC.SAMtools.runLog   2>&1" );
    system( "multiqc    --title BAMtools   --verbose  --export   --outdir $MultiQC4   --pdf         $Bamtools          >> $MultiQC4/MultiQC.BAMtools.runLog   2>&1" );
    system( "multiqc    --title Aligner    --verbose  --export   --outdir $MultiQC5   --pdf    --ignore QC_Results   $dir1   >> $MultiQC5/MultiQC.Aligner.runLog    2>&1" );
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_BAM_2  {
    my $dir1      =  $_[0];   ## All the BAM files must be in this folder.
    opendir(my $FH_Files, $dir1) || die;
    my @Files = readdir($FH_Files);

    my $QCresults = "$dir1/QC_Results";
    my $Fingerprint    = "$QCresults/7_Fingerprint_100bpBin";
    my $Fingerprint2   = "$QCresults/8_Fingerprint_500bpBin";
    my $goleft         = "$QCresults/9_goleft";
    my $MultiQC1       = "$QCresults/10_MultiQC1_Fingerprint_100bpBin";
    my $MultiQC2       = "$QCresults/10_MultiQC2_Fingerprint_500bpBin";
    my $MultiQC3       = "$QCresults/10_MultiQC3_goleft";
    &myMakeDir($QCresults);
    &myMakeDir($Fingerprint);
    &myMakeDir($Fingerprint2);
    &myMakeDir($goleft);
    &myMakeDir($MultiQC1);
    &myMakeDir($MultiQC2);
    &myMakeDir($MultiQC3);

    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of all BAM files by using plotFingerprint in deepTools, goleft and MultiQC ......";
    for ( my $i=0; $i<=$#Files; $i++ ) {
        next unless $Files[$i] =~ m/\.bam$/;
        next unless $Files[$i] !~ m/^[.]/;
        next unless $Files[$i] !~ m/[~]$/;
        my $temp = $Files[$i];
        say    "\t......$temp";
        $temp =~ s/\.bam$//  ||  die;
        system("plotFingerprint  --bamfiles $dir1/$temp.bam     --numberOfSamples 1000000    --plotFile $Fingerprint/$temp.pdf    --plotTitle $temp   --outRawCounts  $Fingerprint/$temp.cov   --outQualityMetrics $Fingerprint/$temp.Metrics.txt   --numberOfProcessors $numCores_g   --binSize 100    >> $Fingerprint/$temp.runLog    2>&1");                           
        system("plotFingerprint  --bamfiles $dir1/$temp.bam     --numberOfSamples 1000000    --plotFile $Fingerprint2/$temp.pdf   --plotTitle $temp   --outRawCounts  $Fingerprint2/$temp.cov  --outQualityMetrics $Fingerprint2/$temp.Metrics.txt  --numberOfProcessors $numCores_g   --binSize 500    >> $Fingerprint2/$temp.runLog   2>&1");                                   
        system("goleft   covstats    $dir1/$temp.bam  > $goleft/$temp.covstats " );
        system("goleft   indexcov  --sex chrX,chrY  -d $goleft/$temp  $dir1/$temp.bam  > $goleft/$temp.indexcov.runLog      2>&1" );
    }
    system("sleep 5s");
    system( "multiqc    --title Fingerprint    --verbose  --export   --outdir $MultiQC1  --pdf    $Fingerprint      >> $MultiQC1/MultiQC.Fingerprint.runLog    2>&1" );
    system( "multiqc    --title Fingerprint    --verbose  --export   --outdir $MultiQC2  --pdf    $Fingerprint2     >> $MultiQC2/MultiQC.Fingerprint.runLog    2>&1" );
    system( "multiqc    --title goleft         --verbose  --export   --outdir $MultiQC3  --pdf    $goleft           >> $MultiQC3/MultiQC.goleft.runLog    2>&1" );

}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_BAM_3  {
    my $dir1      =  $_[0];   ## All the BAM files must be in this folder.
    opendir(my $FH_Files, $dir1) || die;
    my @Files = readdir($FH_Files);

    my $QCresults = "$dir1/QC_Results";
    my $PRESEQ    = "$QCresults/11_PRESEQ";
    my $PicardDir = "$QCresults/12_Picard";
    my $MultiQC1  = "$QCresults/13_MultiQC1_PRESEQ";
    my $MultiQC2  = "$QCresults/13_MultiQC2_Picard";
    &myMakeDir($QCresults);
    &myMakeDir($PRESEQ);
    &myMakeDir($PicardDir);
    &myMakeDir($MultiQC1);
    &myMakeDir($MultiQC2);

    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of all BAM files by using PRESEQ, Picard and MultiQC ......";
    for ( my $i=0; $i<=$#Files; $i++ ) {
        next unless $Files[$i] =~ m/\.bam$/;
        next unless $Files[$i] !~ m/^[.]/;
        next unless $Files[$i] !~ m/[~]$/;
        my $temp = $Files[$i];
        say    "\t......$temp";
        $temp =~ s/\.bam$//  ||  die;
        system("preseq  c_curve     -output  $PRESEQ/$temp.c_curve.pe.PRESEQ       -step 1000000    -verbose   -pe    $dir1/$temp.bam    >> $PRESEQ/$temp.c_curve.pe.runLog   2>&1");
        system("preseq  c_curve     -output  $PRESEQ/$temp.c_curve.se.PRESEQ       -step 1000000    -verbose          $dir1/$temp.bam    >> $PRESEQ/$temp.c_curve.se.runLog   2>&1");
        system("preseq  lc_extrap   -output  $PRESEQ/$temp.lc_extrap.pe.PRESEQ     -step 1000000    -verbose   -pe    $dir1/$temp.bam    >> $PRESEQ/$temp.lc_extrap.pe.runLog   2>&1");
        system("preseq  lc_extrap   -output  $PRESEQ/$temp.lc_extrap.se.PRESEQ     -step 1000000    -verbose          $dir1/$temp.bam    >> $PRESEQ/$temp.lc_extrap.se.runLog   2>&1");

        &myMakeDir("$PicardDir/$temp");
        system("java  -jar   $Picard_g   BamIndexStats                           INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/0_BamIndexStats.txt                                                                                  >> $PicardDir/$temp/0.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectAlignmentSummaryMetrics          INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/1_CollectAlignmentSummaryMetrics.txt                                                                 >> $PicardDir/$temp/1.runLog   2>&1" );
        system("java  -jar   $Picard_g   EstimateLibraryComplexity               INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/2_EstimateLibraryComplexity.txt                                                                      >> $PicardDir/$temp/2.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectInsertSizeMetrics                INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/3_CollectInsertSizeMetrics.txt          HISTOGRAM_FILE=$PicardDir/$temp/3.pdf  MINIMUM_PCT=0.05      >> $PicardDir/$temp/3.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectJumpingLibraryMetrics            INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/4_CollectJumpingLibraryMetrics                                                                       >> $PicardDir/$temp/4.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectMultipleMetrics                  INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/5_CollectMultipleMetrics                                                                             >> $PicardDir/$temp/5.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectBaseDistributionByCycle          INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/6_CollectBaseDistributionByCycle        CHART_OUTPUT=$PicardDir/$temp/6.pdf                          >> $PicardDir/$temp/6.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectQualityYieldMetrics              INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/7_CollectQualityYieldMetrics                                                                         >> $PicardDir/$temp/7.runLog   2>&1" );
        system("java  -jar   $Picard_g   MeanQualityByCycle                      INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/8_MeanQualityByCycle                    CHART_OUTPUT=$PicardDir/$temp/9.pdf                          >> $PicardDir/$temp/8.runLog   2>&1" );
        system("java  -jar   $Picard_g   QualityScoreDistribution                INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/9_QualityScoreDistribution              CHART_OUTPUT=$PicardDir/$temp/10.pdf                         >> $PicardDir/$temp/9.runLog  2>&1" );
    }
    system( "multiqc  --title PRESEQ    --verbose  --export  --outdir $MultiQC1  --pdf          $PRESEQ        >> $MultiQC1/MultiQC.PRESEQ.runLog   2>&1" );
    system( "multiqc  --title Picard    --verbose  --export  --outdir $MultiQC2  --pdf          $PicardDir     >> $MultiQC2/MultiQC.Picard.runLog   2>&1" );
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_BAM_4  {
    my $dir1      =  $_[0];   ## All the BAM files must be in this folder.
    opendir(my $DH_map, $dir1) || die;
    my @mapFiles = readdir($DH_map);

    my $QCresults = "$dir1/QC_Results";
    my $SubreadUti= "$QCresults/14_SubreadUti";
    my $Sambamba  = "$QCresults/15_sambamba_flagstat";
    &myMakeDir("$QCresults");
    &myMakeDir("$SubreadUti");
    &myMakeDir("$Sambamba");

    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of bam files by using Subreads utilities and sambamba ......";
    for (my $i=0; $i<=$#mapFiles; $i++) {
           next unless $mapFiles[$i] =~ m/\.bam$/;
           next unless $mapFiles[$i] !~ m/^[.]/;
           next unless $mapFiles[$i] !~ m/[~]$/;
           my $temp = $mapFiles[$i];
           $temp =~ s/\.bam$//  ||  die;
           say   "\t......$mapFiles[$i]";

           system("propmapped   -i $dir1/$temp.bam                    -o $SubreadUti/$temp.prommapped      >> $SubreadUti/$temp.prommapped.runLog.txt      2>&1");                                                     
           system("propmapped   -i $dir1/$temp.bam       -f           -o $SubreadUti/$temp.prommapped      >> $SubreadUti/$temp.prommapped.runLog.txt      2>&1");
           system("propmapped   -i $dir1/$temp.bam       -f   -p      -o $SubreadUti/$temp.prommapped      >> $SubreadUti/$temp.prommapped.runLog.txt      2>&1");
           system("sambamba  flagstat  -t $numCores_g  -p      $dir1/$temp.bam      >> $Sambamba/$temp.sambamba_flagstat.runLog.txt      2>&1");

     }
}
###################################################################################################################################################################################################

 



###################################################################################################################################################################################################
my $BWA_g   = "$output_g/1_BWA";
&myMakeDir($BWA_g);
{ ########## Start BWA
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using BWA ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {
        my $end1 = $pairedEnd_g[$i];
        my $end2 = $pairedEnd_g[$i+1];
        say    "\n\t......$end1";
        say    "  \t......$end2\n";
        my $temp = $end1; 
        $temp =~ s/\.R1\.fastq//  or die "\n##Error-a1: $temp ##\n\n";
        $temp =~ s/\.gz//;
        $temp =~ s/\.bz2//;
        open(tempFH,  ">>",  "$BWA_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        ## system("bwa mem  -t $numCores_g         $BWA_index_g   $input_g/$end1  $input_g/$end2    > $BWA_g/$temp.sam");
        system("bwa aln    -n 0.04   -t $numCores_g    -f $BWA_g/$end1.sai     $BWA_index_g      $input_g/$end1    >> $BWA_g/$end1.runLog   2>&1");
        system("bwa aln    -n 0.04   -t $numCores_g    -f $BWA_g/$end2.sai     $BWA_index_g      $input_g/$end2    >> $BWA_g/$end2.runLog   2>&1");
        system("bwa sampe  -n 3                        -f $BWA_g/$temp.sam     $BWA_index_g      $BWA_g/$end1.sai  $BWA_g/$end2.sai     $input_g/$end1  $input_g/$end2    >> $BWA_g/$temp.runLog   2>&1"); 
}

for (my $i=0; $i<=$#singleEnd_g; $i++) {
        my $temp = $singleEnd_g[$i];
        say   "\n\t......$temp\n";
        $temp =~ s/\.fastq//  or die "\n##Error-b1: $temp ##\n\n";
        $temp =~ s/\.gz//;
        $temp =~ s/\.bz2//;
        system("bwa mem  -t $numCores_g        $BWA_index_g   $input_g/$singleEnd_g[$i]   > $BWA_g/$temp.sam"); 
}

&myQC_BAM_1($BWA_g);
&myMakeDir("$output2_g/1_MultiQC_BWA");
system( "multiqc    --title BWA        --verbose  --export   --outdir $output2_g/1_MultiQC_BWA      --pdf      $BWA_g     >> $output2_g/1_MultiQC_BWA/MultiQC.BWA.runLog    2>&1" );

}  ########## End BWA
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $Bowtie2_g   = "$output_g/2_Bowtie2";
&myMakeDir($Bowtie2_g);
{ ########## Start Bowtie2
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using Bowtie2 ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {
        my $end1 = $pairedEnd_g[$i];
        my $end2 = $pairedEnd_g[$i+1];
        say    "\n\t......$end1";
        say    "  \t......$end2\n";
        my $temp = $end1; 
        $temp =~ s/\.R1\.fastq//  or die "\n##Error-a2: $temp ##\n\n";
        $temp =~ s/\.gz//;
        $temp =~ s/\.bz2//;
        open(tempFH, ">>", "$Bowtie2_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("bowtie2     --threads $numCores_g   -q   --phred33      -x $Bowtie2_index_g    -1 $input_g/$end1        -2 $input_g/$end2     -S $Bowtie2_g/$temp.sam    >> $Bowtie2_g/$temp.runLog  2>&1");

}

for (my $i=0; $i<=$#singleEnd_g; $i++) {
        my $temp = $singleEnd_g[$i];
        say   "\n\t......$temp\n";
        $temp =~ s/\.fastq//  or die "\n##Error-b2: $temp ##\n\n";
        $temp =~ s/\.gz//;
        $temp =~ s/\.bz2//;
        system("bowtie2    --threads $numCores_g   -q   --phred33     -x $Bowtie2_index_g    -U $input_g/$singleEnd_g[$i]       -S $Bowtie2_g/$temp.sam    >> $Bowtie2_g/$temp.runLog  2>&1");
}
&myQC_BAM_1($Bowtie2_g);

&myMakeDir("$output2_g/2_MultiQC_Bowtie2");
system( "multiqc    --title Bowtie2        --verbose  --export   --outdir $output2_g/2_MultiQC_Bowtie2     --ignore QC_Results   --pdf    $Bowtie2_g     >> $output2_g/2_MultiQC_Bowtie2/MultiQC.Bowtie2.runLog    2>&1" );

}  ########## End Bowtie2
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $Bowtie2_g2   = "$output_g/3_Bowtie2_Local";
&myMakeDir($Bowtie2_g2);
{ ########## Start Bowtie2
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using Bowtie2 ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {
        my $end1 = $pairedEnd_g[$i];
        my $end2 = $pairedEnd_g[$i+1];
        say    "\n\t......$end1";
        say    "  \t......$end2\n";
        my $temp = $end1; 
        $temp =~ s/\.R1\.fastq//  or die "\n##Error-a2: $temp ##\n\n";
        $temp =~ s/\.gz//;
        $temp =~ s/\.bz2//;
        open(tempFH, ">>", "$Bowtie2_g2/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("bowtie2     --threads $numCores_g   -q   --phred33   --local   -x $Bowtie2_index_g    -1 $input_g/$end1        -2 $input_g/$end2     -S $Bowtie2_g2/$temp.sam    >> $Bowtie2_g2/$temp.runLog  2>&1");

}

for (my $i=0; $i<=$#singleEnd_g; $i++) {
        my $temp = $singleEnd_g[$i];
        say   "\n\t......$temp\n";
        $temp =~ s/\.fastq//  or die "\n##Error-b2: $temp ##\n\n";
        $temp =~ s/\.gz//;
        $temp =~ s/\.bz2//;
        system("bowtie2    --threads $numCores_g   -q   --phred33  --local   -x $Bowtie2_index_g    -U $input_g/$singleEnd_g[$i]       -S $Bowtie2_g2/$temp.sam    >> $Bowtie2_g2/$temp.runLog  2>&1");
}
&myQC_BAM_1($Bowtie2_g2);

&myMakeDir("$output2_g/3_MultiQC_Bowtie2_Local");
system( "multiqc    --title Bowtie2        --verbose  --export   --outdir $output2_g/3_MultiQC_Bowtie2_Local     --ignore QC_Results   --pdf    $Bowtie2_g2     >> $output2_g/3_MultiQC_Bowtie2_Local/MultiQC.Bowtie2_Local.runLog    2>&1" );

}  ########## End Bowtie2
###################################################################################################################################################################################################





###################################################################################################################################################################################################
&myQC_BAM_2($BWA_g);
&myQC_BAM_2($Bowtie2_g);
&myQC_BAM_2($Bowtie2_g2);

&myQC_BAM_3($BWA_g);
&myQC_BAM_3($Bowtie2_g);
&myQC_BAM_3($Bowtie2_g2);

&myQC_BAM_4($BWA_g);
&myQC_BAM_4($Bowtie2_g);
&myQC_BAM_4($Bowtie2_g2);
###################################################################################################################################################################################################






###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";





## END
