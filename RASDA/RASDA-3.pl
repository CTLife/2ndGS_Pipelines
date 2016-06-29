#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.20;       
## perl5 version >= 5.20,   you can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
## Help Infromation
my $HELP_g = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use RASDA (RNA-Seq Data Analyzer), version 0.7.3, 2016-06-01.      
        RASDA is a Pipeline for Single-end and Paired-end RNA-Seq Data Analysis by Integrating Lots of Softwares.

        Step 3: Mapping reads to the reference genome or transcriptome by using 6 softwares (mappers or aligners): 
                    RSEM+Bowtie2, Kallisto, Sailfish, Salmon; STAR, Trim_STAR, HISAT2, Trim_HISAT2. 
                Assess the quality of BAM files to identify possible sequencing errors or biases by using 15 softwares:
                    SAMtools, Subread utilities, FASTQC, SAMstat, BAMStats, qualimap, MultiQC, ezBAMQC, PRESEQ, Picard, fastqp, BamQC, QoRTs, RNA-SeQC and RSeQC.

        Usage:  
               perl  RASDA3.pl    [-version]    [-help]    [-in inputDir]    [-out outDir]   [-genome reference_genome]    
        For instance: 
               perl  RASDA3.pl   -in 3-Filtered    -out 4-Mapping  -genome mm10   >> RASDA3.runLog  2>&1

        -------------------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.

        -help           Show this help message and exit.

        Required arguments:
        -in inputDir        inputDir is the name of your input folder that contains your FASTQ files, the suffix of the FASTQ files must be ".fastq".    (no default)

        -out outDir         outDir is the name of your output folder that contains running results (BAM format) of this step.   (no default)   
                                                  
        -genome reference_genome        NGS reads wil be mapped to "reference_genome". 
                                        10 UCSC reference genomes are available:  hg38, mm10, dm6, ce11, sacCer3, danRer10, oryCun2, rn6, rheMac3, panTro4.  
                                              		hg38     for  Homo sapiens.
                                             		mm10     for  Mus musculus.
                                              		dm6      for  Drosophila melanogaster.
                                              		ce11     for  Caenorhabditis elegans.
                                              		sacCer3  for  Saccharomyces cerevisiae.
                                              		danRer10 for  Danio rerio (Zebrafish).
                                              		oryCun2  for  Rabbit.
                                              		rn6      for  Rat.
                                              		rheMac3  for  Rhesus.
                                              		panTro4  for  Chimp.                                                  
        ------------------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines such as CISDA, MESDA and HISDA,
        please see the web site: https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies 
        and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.    
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------  
';

## Version Infromation 
my $version_g = "  The Third Step of RASDA (RNA-Seq Data Analyzer), version 0.7.3, 2016-06-01.";

## Keys and Values
if ($#ARGV   == -1) { say  "\n$HELP_g\n";  exit 0; }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0) { @ARGV = (@ARGV, "-h");       }       ## when the number of command argumants is odd. 
my %args = @ARGV;

## Initialize  Variables
my $input_g  = '3-Filtered';      ## This is only an initialization  value or suggesting value, not default value.
my $output_g = '4-Mapping';       ## This is only an initialization  value or suggesting value, not default value.
my $genome_g = 'mm10';            ## This is only an initialization  value or suggesting value, not default value.

## Available Arguments
my $available = "  -version    -help   -in   -out  -genome   ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if($available !~ m/\s$key\s/) {say    "\n\tCann't recognize $key !!";  $boole_g = 1; }
}
if($boole_g == 1) {
    say   "\tThe Command Line Arguments are wrong!";
    say   "\tPlease see help message by using 'perl  RASDA3.pl  -help' \n";
    exit 0;
}

## Get Arguments
if ( exists $args{'-version' }   )     { say  "\n$version_g\n";   exit 0; }
if ( exists $args{'-help'    }   )     { say  "\n$HELP_g\n";      exit 0; }
if ( exists $args{'-in'      }   )     { $input_g  = $args{'-in' };       }else{say   "\n -in  is required.\n";          say  "\n$HELP_g\n";    exit 0; }
if ( exists $args{'-out'     }   )     { $output_g = $args{'-out'};       }else{say   "\n -out is required.\n";          say  "\n$HELP_g\n";    exit 0; }
if ( exists $args{'-genome'  }   )     { $genome_g = $args{'-genome'  };  }else{say   "\n -genome is required.\n";       say  "\n$HELP_g\n";    exit 0; } 

## Conditions
$input_g   =~ m/^\S+$/   ||  die   "The Command Line Arguments are wrong!\n$HELP_g\n\n";
$output_g  =~ m/^\S+$/   ||  die   "The Command Line Arguments are wrong!\n$HELP_g\n\n";
$genome_g  =~ m/^\S+$/   ||  die   "The Command Line Arguments are wrong!\n$HELP_g\n\n";

## Print Command Arguments to Standard Output
say  "\n
        ################ Arguments ###############################
                Input  folder:  $input_g
                Output folder:  $output_g
                Reference genome: $genome_g
        ###############################################################  
\n";
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
my $numCores_g   = 6;
###################################################################################################################################################################################################





###################################################################################################################################################################################################
## Context specific:
my  $commonPath1_g      = "/home/yp/.MyProgramFiles/3_HTS2G/4_ChIPseq/3_Aligners";
my  $TopHat2_index_g    = "$commonPath1_g/bowtie2-2.2.9/RefGenomes/$genome_g/$genome_g";
my  $Subjunc_index_g    = "$commonPath1_g/subread-1.5.0-p2/RefGenomes/$genome_g/$genome_g";
my  $mRNA_index = $genome_g."_mRNA";
my  $commonPath2_g      = "/home/yp/.MyProgramFiles/3_HTS2G/5_RNAseq/3_Aligners";
my  $STAR_index_g       = "$commonPath2_g/STAR-2.5.2a/RefGenomes/$genome_g";
my  $HISAT2_index_g     = "$commonPath2_g/hisat2-2.0.4/RefGenomes/$genome_g/$genome_g";
my  $RSEM_index_g       = "$commonPath2_g/RSEM/RefGenomes/$genome_g/$mRNA_index";
my  $Kallisto_index_g   = "$commonPath2_g/kallisto/RefGenomes/$genome_g/$mRNA_index";
my  $Sailfish_index_g   = "$commonPath2_g/SailfishBeta-0.10.0/RefGenomes/$genome_g/$mRNA_index";
my  $Salmon_index_g     = "$commonPath2_g/salmon-0.6.0/RefGenomes/$genome_g/$mRNA_index";
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
my  $BAMStats_g = &fullPathApp("BAMStats.jar");
my  $QoRTs_g    = &fullPathApp("QoRTs.jar");
my  $RNASeQC_g  = &fullPathApp("RNA-SeQC.jar");
my  $trim5_g    = 3;   ##bp
my  $trim3_g    = 3;   ##bp
&printVersion("bowtie2   --version");                      
&printVersion("tophat   -v");
&printVersion("subjunc  -v");
&printVersion("STAR  --version");
&printVersion("hisat2   --version ");  
&printVersion("rsem-calculate-expression   --version");
&printVersion("kallisto");
&printVersion("sailfish  --help ");
&printVersion("salmon    --help ");
&printVersion("samtools");
&printVersion("fastqc   -v");
&printVersion("samstat   -v");
&printVersion("java  -jar  $BAMStats_g  -h");
&printVersion("bamqc  -v");
&printVersion("preseq");
&printVersion("qualimap  -v");
&printVersion("fastqp   -h");
&printVersion("multiqc   --version");
&printVersion("propmapped");
&printVersion("qualityScores");  
&printVersion("ezBAMQC  -h");
&printVersion("java  -jar  $Picard_g   CollectAlignmentSummaryMetrics      --version");
&printVersion("java  -jar  $Picard_g   EstimateLibraryComplexity           --version");
&printVersion("java  -jar  $Picard_g   CollectInsertSizeMetrics            --version");
&printVersion("java  -jar  $Picard_g   CollectJumpingLibraryMetrics        --version");
&printVersion("java  -jar  $Picard_g   CollectMultipleMetrics              --version");
&printVersion("java  -jar  $Picard_g   CollectBaseDistributionByCycle      --version");
&printVersion("java  -jar  $Picard_g   CollectQualityYieldMetrics          --version");
&printVersion("java  -jar  $Picard_g   CollectWgsMetricsFromQuerySorted    --version");
&printVersion("java  -jar  $Picard_g   MeanQualityByCycle                  --version");
&printVersion("java  -jar  $Picard_g   QualityScoreDistribution            --version");
&printVersion("java  -jar  $QoRTs_g      --man QC");
&printVersion("java  -jar  $RNASeQC_g    --version");
&printVersion("tin.py                  --version");     ##RSeQC
&printVersion("bam_stat.py             --version");     ##RSeQC
&printVersion("clipping_profile.py     --version");     ##RSeQC
&printVersion("deletion_profile.py     --version");     ##RSeQC
&printVersion("geneBody_coverage.py    --version");     ##RSeQC
&printVersion("inner_distance.py       --version");     ##RSeQC
&printVersion("insertion_profile.py    --version");     ##RSeQC
&printVersion("junction_annotation.py  --version");     ##RSeQC
&printVersion("junction_saturation.py  --version");     ##RSeQC
&printVersion("mismatch_profile.py     --version");     ##RSeQC
&printVersion("read_distribution.py    --version");     ##RSeQC
&printVersion("read_duplication.py     --version");     ##RSeQC
&printVersion("read_GC.py              --version");     ##RSeQC
&printVersion("read_NVC.py             --version");     ##RSeQC
&printVersion("read_quality.py         --version");     ##RSeQC
&printVersion("RPKM_saturation.py      --version");     ##RSeQC
###################################################################################################################################################################################################


 


###################################################################################################################################################################################################              
{
say   "\n\n\n\n\n\n##################################################################################################";
say   "Checking all the input file names ......";
my @groupFiles = ();     
my $fileNameBool = 1;
for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {   
        next unless $inputFiles_g[$i] =~ m/\.fastq$/;  
        next unless $inputFiles_g[$i] !~ m/^[.]/;
        next unless $inputFiles_g[$i] !~ m/[~]$/;
        next unless $inputFiles_g[$i] !~ m/^QC_Results$/;
        next unless $inputFiles_g[$i] !~ m/^unpaired/;
        say   "\t......$inputFiles_g[$i]" ; 
        my $temp = $inputFiles_g[$i]; 
        $groupFiles[++$#groupFiles] = $inputFiles_g[$i];  
        $temp =~ m/^(\d+)_($pattern_g)_(Rep[1-9])/   or  die   "wrong-1: ## $temp ##";
        $temp =~ m/_(Rep[1-9])\.fastq$/  or  $temp =~ m/_(Rep[1-9])_?([1-2]?)\.fastq$/   or  die   "wrong-2: ## $temp ##";
        if($temp !~ m/^((\d+)_($pattern_g)_(Rep[1-9]))(_[1-2])?\.fastq$/) {
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
say   "Detecting single-end and paired-end FASTQ files in input folder ......";     ## The fastq files are same between input folder and ouput folder.
my @singleEnd_g   = ();
my @pairedEnd_g   = ();
open(seqFiles_FH_g, ">", "$output2_g/singleEnd-pairedEnd-Files.txt")  or  die;
for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {
    next unless $inputFiles_g[$i] =~ m/\.fastq$/;
    next unless $inputFiles_g[$i] !~ m/^[.]/;
    next unless $inputFiles_g[$i] !~ m/[~]$/;
    next unless $inputFiles_g[$i] !~ m/^unpaired/;
    next unless $inputFiles_g[$i] !~ m/^QC_Results$/;
    say    "\t......$inputFiles_g[$i]";
    $inputFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])_?([1-2]?)\.fastq$/   or  die;
    if ($inputFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])\.fastq$/) {   ## sinlge end sequencing files.
        $inputFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])\.fastq$/  or  die;
        $singleEnd_g[$#singleEnd_g+1] =  $inputFiles_g[$i];
        say         "\t\t\t\tSingle-end sequencing files: $inputFiles_g[$i]\n";
        say  seqFiles_FH_g  "Single-end sequencing files: $inputFiles_g[$i]\n";
    }else{     ## paired end sequencing files.
        $inputFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])_([1-2])\.fastq$/  or  die;
        if ($inputFiles_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/) { ## The two files of one paired sequencing sample are always side by side.
            my $temp = $1;
            my $end1 = $temp."_1.fastq";
            my $end2 = $temp."_2.fastq";
            (-e  "$input_g/$end1")  or die;
            (-e  "$input_g/$end2")  or die;
            $pairedEnd_g[$#pairedEnd_g+1] =  $end1;
            $pairedEnd_g[$#pairedEnd_g+1] =  $end2;
            say        "\t\t\t\tPaired-end sequencing files: $end1,  $end2\n";
            say seqFiles_FH_g  "Paired-end sequencing files: $end1,  $end2\n";
        }
    }
}
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
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_BAM_1  {
    my $dir1      =  $_[0];   ## All the SAM files must be in this folder.
    my $QCresults = "$dir1/QC_Results";
    my $SAMtools  = "$QCresults/1_SAMtools";
    my $FastQC    = "$QCresults/2_FastQC";
    my $qualimap  = "$QCresults/3_qualimap";
    my $samstat   = "$QCresults/4_samstat";
    my $MultiQC1  = "$QCresults/5_MultiQC_FastQC";
    my $MultiQC2  = "$QCresults/5_MultiQC_qualimap";
    &myMakeDir($QCresults);
    &myMakeDir($SAMtools);
    &myMakeDir($FastQC);
    &myMakeDir($qualimap);
    &myMakeDir($samstat);
    &myMakeDir($MultiQC1);
    &myMakeDir($MultiQC2);
    opendir(my $FH_Files, $dir1) || die;     
    my @Files = readdir($FH_Files);
    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of all SAM files by using SAMtools, FastQC, qualimap, samstat and MultiQC ......";
    for ( my $i=0; $i<=$#Files; $i++ ) {
        next unless $Files[$i] =~ m/\.sam$/;
        next unless $Files[$i] !~ m/^[.]/;
        next unless $Files[$i] !~ m/[~]$/;
        my $temp = $Files[$i];
        say    "\t......$temp";
        $temp =~ s/\.sam$//  ||  die;
        system("samtools  sort  -m 2G  -o $dir1/$temp.bam   --output-fmt bam  -T $dir1/yp_$temp   --threads $numCores_g    $dir1/$temp.sam    >>$SAMtools/$temp.runLog    2>&1");
        system("samtools  index           $dir1/$temp.bam      >>$SAMtools/$temp.runLog  2>&1");
        system("samtools  flagstat        $dir1/$temp.bam      >>$SAMtools/$temp.runLog  2>&1");
        system(`samtools  idxstats        $dir1/$temp.bam      >>$SAMtools/$temp.runLog  2>&1`);
        system( "fastqc    --outdir $FastQC    --threads $numCores_g  --format bam   --kmers 7    $dir1/$temp.bam                     >> $FastQC/$temp.runLog      2>&1" );
        system( "qualimap  bamqc   -bam $dir1/$temp.bam   -c  -ip  -nt $numCores_g   -outdir $qualimap/$temp   --java-mem-size=12G    >> $qualimap/$temp.runLog    2>&1" );
        system( "qualimap  rnaseq  -bam $dir1/$temp.bam   -gtf 0-Other/Shortcuts/mm10_RefSeq_GTF   -oc $temp.counts   -outdir $qualimap/rnaseq_$temp   --java-mem-size=12G   >> $qualimap/rnaseq_$temp.runLog    2>&1" );
        system( "samstat   $dir1/$temp.bam      >> $samstat/$temp.runLog         2>&1");  
        system( "rm   $dir1/$temp.sam" );  
    }
    system( "multiqc  --verbose  --outdir $MultiQC1          $FastQC/*_fastqc.zip      >> $MultiQC1/MultiQC.FastQC.runLog     2>&1" );
    system( "multiqc  --verbose  --outdir $MultiQC2          $qualimap/*               >> $MultiQC2/MultiQC.qualimap.runLog   2>&1" );
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_BAM_2  {
    my $dir1      =  $_[0];   ## All the BAM files must be in this folder.   
    my $QCresults = "$dir1/QC_Results";
    my $PicardDir = "$QCresults/6_Picard";
    my $BAMStats  = "$QCresults/7_BAMStats";
    my $ezBAMQC   = "$QCresults/8_ezBAMQC";
    my $MultiQC   = "$QCresults/9_MultiQC_Picard";
    &myMakeDir($QCresults);
    &myMakeDir($PicardDir);
    &myMakeDir($MultiQC);
    &myMakeDir($BAMStats);
    &myMakeDir($ezBAMQC);
    opendir(my $FH_Files, $dir1) || die;     
    my @Files = readdir($FH_Files);
    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of all BAM files by using Picard, MultiQC, BAMStats, ezBAMQC ......";
    for ( my $i=0; $i<=$#Files; $i++ ) {
        next unless $Files[$i] =~ m/\.bam$/;
        next unless $Files[$i] !~ m/^[.]/;
        next unless $Files[$i] !~ m/[~]$/;
        my $temp = $Files[$i];
        say    "\t......$temp";
        $temp =~ s/\.bam$//  ||  die;
 
        &myMakeDir("$PicardDir/$temp"); 
        system("java  -jar   $Picard_g   CollectAlignmentSummaryMetrics      INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/1_CollectAlignmentSummaryMetrics     R=0-Other/Shortcuts/$genome_g.fa                            >> $PicardDir/$temp/1.runLog   2>&1" );
        system("java  -jar   $Picard_g   EstimateLibraryComplexity           INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/2_EstimateLibraryComplexity                                                                      >> $PicardDir/$temp/2.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectInsertSizeMetrics            INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/3_CollectInsertSizeMetrics          HISTOGRAM_FILE=$PicardDir/$temp/3.pdf  MINIMUM_PCT=0.01      >> $PicardDir/$temp/3.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectJumpingLibraryMetrics        INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/4_CollectJumpingLibraryMetrics                                                                   >> $PicardDir/$temp/4.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectMultipleMetrics              INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/5_CollectMultipleMetrics                                                                         >> $PicardDir/$temp/5.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectBaseDistributionByCycle      INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/6_CollectBaseDistributionByCycle     CHART_OUTPUT=$PicardDir/$temp/6.pdf                         >> $PicardDir/$temp/6.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectQualityYieldMetrics          INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/7_CollectQualityYieldMetrics                                                                     >> $PicardDir/$temp/7.runLog   2>&1" ); 
        system("java  -jar   $Picard_g   CollectWgsMetricsFromQuerySorted    INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/8_CollectWgsMetricsFromQuerySorted                                                               >> $PicardDir/$temp/8.runLog   2>&1" );
        system("java  -jar   $Picard_g   MeanQualityByCycle                  INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/9_MeanQualityByCycle                 CHART_OUTPUT=$PicardDir/$temp/9.pdf                         >> $PicardDir/$temp/9.runLog   2>&1" );
        system("java  -jar   $Picard_g   QualityScoreDistribution            INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/10_QualityScoreDistribution          CHART_OUTPUT=$PicardDir/$temp/10.pdf                        >> $PicardDir/$temp/10.runLog  2>&1" ); 

        system("java  -jar  $BAMStats_g    --distances  --infile $dir1/$temp.bam   --lengths   --mapped    --outfile $BAMStats/$temp   --qualities   --starts   --view html        >> $BAMStats/$temp.runLog   2>&1");  
        system("ezBAMQC    -i $dir1/$temp.bam   --refgene 0-Other/Shortcuts/mm10_RefSeq_GTF   --outputDir $ezBAMQC/$temp  --stranded no  --label $temp  --threads $numCores_g      >> $ezBAMQC/$temp.runLog   2>&1");  
    }
    system( "multiqc  --verbose  --outdir $MultiQC          $PicardDir/*                 >> $MultiQC/MultiQC.Picard.runLog   2>&1" );
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_BAM_3  {
    my $dir1      =  $_[0];   ## All the BAM files must be in this folder.   
    my $QCresults = "$dir1/QC_Results";
    my $PRESEQ    = "$QCresults/10_PRESEQ";
    my $RSeQC     = "$QCresults/11_RSeQC";
    my $MultiQC   = "$QCresults/12_MultiQC_PRESEQ";
    &myMakeDir($QCresults);
    &myMakeDir($PRESEQ);
    &myMakeDir($MultiQC);
    &myMakeDir($RSeQC);
    opendir(my $FH_Files, $dir1) || die;     
    my @Files = readdir($FH_Files);
    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of all BAM files by using PRESEQ, MultiQC and RSeQC ......";
    for ( my $i=0; $i<=$#Files; $i++ ) {
        next unless $Files[$i] =~ m/\.bam$/;
        next unless $Files[$i] !~ m/^[.]/;
        next unless $Files[$i] !~ m/[~]$/;
        my $temp = $Files[$i];
        say    "\t......$temp";
        $temp =~ s/\.bam$//  ||  die;
        system("preseq  c_curve   -output  $PRESEQ/$temp.pe.PRESEQ     -step 1000000    -verbose   -pe  -bam  $dir1/$temp.bam    >> $PRESEQ/$temp.pe.runLog   2>&1");   
        system("preseq  c_curve   -output  $PRESEQ/$temp.se.PRESEQ     -step 1000000    -verbose        -bam  $dir1/$temp.bam    >> $PRESEQ/$temp.se.runLog   2>&1"); 

        &myMakeDir("$RSeQC/$temp");
        system("tin.py                  --input=$dir1/$temp.bam                                                        --refgene=0-Other/Shortcuts/mm10_RefSeq.bed         >> $RSeQC/$temp/1-tin.runLog                   2>&1");     
        system("bam_stat.py             --input-file=$dir1/$temp.bam                                                                                                       >> $RSeQC/$temp/2-bam_stat.runLog              2>&1");     
        system("clipping_profile.py     --input-file=$dir1/$temp.bam    --out-prefix=$RSeQC/$temp/3-clipping_profile     --sequencing=PE                                   >> $RSeQC/$temp/3-clipping_profile.runLog      2>&1");      
        system("deletion_profile.py     --input=$dir1/$temp.bam         --out-prefix=$RSeQC/$temp/4-deletion_profile     --read-align-length=50                            >> $RSeQC/$temp/4-deletion_profile.runLog      2>&1");      
        system("geneBody_coverage.py    --input=$dir1/$temp.bam         --out-prefix=$RSeQC/$temp/5-geneBody_coverage    --refgene=0-Other/Shortcuts/mm10_RefSeq.bed       >> $RSeQC/$temp/5-geneBody_coverage.runLog     2>&1");      
        system("inner_distance.py       --input-file=$dir1/$temp.bam    --out-prefix=$RSeQC/$temp/6-inner_distance       --refgene=0-Other/Shortcuts/mm10_RefSeq.bed       >> $RSeQC/$temp/6-inner_distance.runLog        2>&1");      
        system("insertion_profile.py    --input-file=$dir1/$temp.bam    --out-prefix=$RSeQC/$temp/7-insertion_profile    --sequencing=PE                                   >> $RSeQC/$temp/7-insertion_profile.runLog     2>&1");      
        system("junction_annotation.py  --input-file=$dir1/$temp.bam    --out-prefix=$RSeQC/$temp/8-junction_annotation  --refgene=0-Other/Shortcuts/mm10_RefSeq.bed       >> $RSeQC/$temp/8-junction_annotation.runLog   2>&1");      
        system("junction_saturation.py  --input-file=$dir1/$temp.bam    --out-prefix=$RSeQC/$temp/9-junction_saturation  --refgene=0-Other/Shortcuts/mm10_RefSeq.bed       >> $RSeQC/$temp/9-junction_saturation.runLog   2>&1");      
        system("mismatch_profile.py     --input=$dir1/$temp.bam         --out-prefix=$RSeQC/$temp/10-mismatch_profile     --read-align-length=50                            >> $RSeQC/$temp/10-mismatch_profile.runLog     2>&1");      
        system("read_distribution.py    --input-file=$dir1/$temp.bam                                                     --refgene=0-Other/Shortcuts/mm10_RefSeq.bed       >> $RSeQC/$temp/11-read_distribution.runLog    2>&1");      
        system("read_duplication.py     --input-file=$dir1/$temp.bam    --out-prefix=$RSeQC/$temp/11-read_duplication                                                       >> $RSeQC/$temp/12-read_duplication.runLog     2>&1");      
        system("read_GC.py              --input-file=$dir1/$temp.bam    --out-prefix=$RSeQC/$temp/12-read_GC                                                               >> $RSeQC/$temp/13-read_GC.runLog              2>&1");      
        system("read_NVC.py             --input-file=$dir1/$temp.bam    --out-prefix=$RSeQC/$temp/13-read_NVC             --nx                                             >> $RSeQC/$temp/14-read_NVC.runLog             2>&1");      
        system("read_quality.py         --input-file=$dir1/$temp.bam    --out-prefix=$RSeQC/$temp/14-read_quality                                                          >> $RSeQC/$temp/15-read_quality.runLog         2>&1");      
        system("RPKM_saturation.py      --input-file=$dir1/$temp.bam    --out-prefix=$RSeQC/$temp/15-RPKM_saturation      --refgene=0-Other/Shortcuts/mm10_RefSeq.bed      >> $RSeQC/$temp/16-RPKM_saturation.runLog      2>&1");      
    }
    system( "multiqc  --verbose  --outdir $MultiQC          $PRESEQ/*                 >> $MultiQC/MultiQC.PRESEQ.runLog   2>&1" );
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_BAM_4  {
    my $dir1      =  $_[0];   ## All the BAM files must be in this folder.
    my $QCresults = "$dir1/QC_Results";
    my $SubreadUti= "$QCresults/12_SubreadUti";
    my $Fastqp    = "$QCresults/13_fastqp";
    my $QoRTs     = "$QCresults/14_QoRTs";
    my $RNASeQC   = "$QCresults/15_RNASeQC";
    my $BamQC     = "$QCresults/16_BamQC";
    &myMakeDir("$QCresults");
    &myMakeDir("$SubreadUti");
    &myMakeDir("$Fastqp");
    &myMakeDir("$QoRTs");
    &myMakeDir("$RNASeQC");
    &myMakeDir("$BamQC");
    opendir(my $DH_map, $dir1) || die;     
    my @mapFiles = readdir($DH_map);
    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of bam files by using Subreads utilities, fastqp, BamQC, QoRTs and RNASeQC ......";
    for (my $i=0; $i<=$#mapFiles; $i++) {
           next unless $mapFiles[$i] =~ m/\.bam$/;
           next unless $mapFiles[$i] !~ m/^[.]/;
           next unless $mapFiles[$i] !~ m/[~]$/;
           my $temp = $mapFiles[$i]; 
           $temp =~ s/\.bam$//  ||  die; 
           say   "\t......$mapFiles[$i]";  
           system("propmapped   -i $dir1/$temp.bam                    -o $SubreadUti/$temp.prommapped      >> $SubreadUti/$temp.prommapped      2>&1");
           system("echo      '\n\n\n\n\n'                                                                  >> $SubreadUti/$temp.prommapped      2>&1"); 
           system("propmapped   -i $dir1/$temp.bam       -f           -o $SubreadUti/$temp.prommapped      >> $SubreadUti/$temp.prommapped      2>&1"); 
           system("echo      '\n\n\n\n\n'                                                                  >> $SubreadUti/$temp.prommapped      2>&1"); 
           system("propmapped   -i $dir1/$temp.bam       -f   -p      -o $SubreadUti/$temp.prommapped      >> $SubreadUti/$temp.prommapped      2>&1"); 
           system("qualityScores   --BAMinput   -i $dir1/$temp.bam    -o $SubreadUti/$temp.qualityScores   >> $SubreadUti/$temp.qualityScores   2>&1");
           system("fastqp    --nreads 50000000   --kmer 4    --output $Fastqp/$temp  --type bam   --median-qual 30     $dir1/$temp.bam     >> $Fastqp/$temp.runLog     2>&1 " );
           system("java -Xmx18000M -Xms5000M -jar $QoRTs_g QC  --generatePlots  $dir1/$temp.bam   0-Other/Shortcuts/mm10_RefSeq_GTF    $QoRTs/$temp   >> $QoRTs/$temp.runLog     2>&1  " );
           ##system("java   -jar $RNASeQC_g       " );
     }
     system( "bamqc  $dir1   --outdir  $BamQC        >> $BamQC/bamqc.runLog     2>&1 " );
}  
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $RSEM_g   = "$output_g/1_RSEM";  
&myMakeDir($RSEM_g);
{ ########## Start RSEM
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference transcriptome by using RSEM ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {  
        say    "\n\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]\n";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$RSEM_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("rsem-calculate-expression  -p $numCores_g  --paired-end      --bowtie2    --estimate-rspd   --append-names    $input_g/$end1.fastq  $input_g/$end2.fastq   $RSEM_index_g    $RSEM_g/$temp    >>$RSEM_g/$temp.runLog  2>&1 ");    
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {   
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("rsem-calculate-expression  -p $numCores_g  --fragment-length-mean 250    --fragment-length-sd 100     --bowtie2    --estimate-rspd   --append-names    $input_g/$temp.fastq     $RSEM_index_g    $RSEM_g/$temp    >>$RSEM_g/$temp.runLog  2>&1 ");    
}
}  ########## End RSEM
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $Kallisto_g   = "$output_g/2_Kallisto";  
&myMakeDir($Kallisto_g);
{ ########## Start Kallisto
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference transcriptome by using Kallisto ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {  
        say    "\n\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]\n";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$Kallisto_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("kallisto quant  --threads=$numCores_g       --index=$Kallisto_index_g    --output-dir=$Kallisto_g/$temp    $input_g/$end1.fastq  $input_g/$end2.fastq    >>$Kallisto_g/$temp.runLog  2>&1 ");    
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {   
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("kallisto quant  --threads=$numCores_g   --single -l 250 -s 100    --index=$Kallisto_index_g    --output-dir=$Kallisto_g/$temp    $input_g/$temp.fastq    >>$Kallisto_g/$temp.runLog  2>&1");    
}
}  ########## End Kallisto
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $Sailfish_g   = "$output_g/3_Sailfish";  
&myMakeDir($Sailfish_g);
{ ########## Start Sailfish
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference transcriptome by using Sailfish ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {  
        say    "\n\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]\n";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$Sailfish_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("sailfish quant  --threads $numCores_g    --libType IU      --index $Sailfish_index_g    --output $Sailfish_g/$temp    --mates1 $input_g/$end1.fastq  --mates2 $input_g/$end2.fastq    >>$Sailfish_g/$temp.runLog  2>&1 ");    
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {   
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("sailfish quant  --threads $numCores_g    --libType U      --index $Sailfish_index_g    --output $Sailfish_g/$temp    --unmatedReads $input_g/$temp.fastq      >>$Sailfish_g/$temp.runLog  2>&1 ");                         
}
}  ########## End Sailfish
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $Salmon_g   = "$output_g/4_Salmon";  
&myMakeDir($Salmon_g);
{ ########## Start Salmon
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference transcriptome by using Salmon ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {  
        say    "\n\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]\n";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$Salmon_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("salmon quant  --threads $numCores_g   --libType IU      --index $Salmon_index_g    --output $Salmon_g/$temp    --mates1 $input_g/$end1.fastq  --mates2 $input_g/$end2.fastq    >>$Salmon_g/$temp.runLog  2>&1 ");    
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {   
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("salmon quant  --threads $numCores_g    --libType U      --index $Salmon_index_g    --output $Salmon_g/$temp    --unmatedReads $input_g/$temp.fastq      >>$Salmon_g/$temp.runLog  2>&1 ");                         
}
}  ########## End Salmon
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $STAR_g   = "$output_g/5_STAR";  
&myMakeDir($STAR_g);
{ ########## Start STAR
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using STAR ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {  
        say    "\n\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]\n";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$STAR_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("STAR  --runMode alignReads    --runThreadN $numCores_g    --sjdbGTFfile 0-Other/Shortcuts/mm10_RefSeq_GTF    --quantMode   GeneCounts    --outFileNamePrefix  $STAR_g/$temp.   --genomeDir $STAR_index_g  --readFilesIn $input_g/$end1.fastq  $input_g/$end2.fastq   >>$STAR_g/$temp.runLog  2>&1 ");    
        system("rename  s/Aligned.out.sam/sam/   $STAR_g/*Aligned.out.sam");                                                                                                                        
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {   
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("STAR  --runMode alignReads   --runThreadN $numCores_g    --sjdbGTFfile 0-Other/Shortcuts/mm10_RefSeq_GTF     --quantMode   GeneCounts     --outFileNamePrefix  $STAR_g/$temp.   --genomeDir $STAR_index_g  --readFilesIn $input_g/$temp.fastq   >>$STAR_g/$temp.runLog  2>&1 ");    
        system("rename  s/Aligned.out.sam/sam/   $STAR_g/*Aligned.out.sam");                                                            
}
}  ########## End STAR
&myQC_BAM_1($STAR_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $Trim_STAR_g   = "$output_g/6_Trim_STAR";  
&myMakeDir($Trim_STAR_g);
my $input2_g = "2-FASTQ";
{ ########## Start STAR
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using Trim_STAR ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {  
        say    "\n\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]\n";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$Trim_STAR_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("STAR  --runMode alignReads    --runThreadN $numCores_g    --sjdbGTFfile 0-Other/Shortcuts/mm10_RefSeq_GTF  --quantMode   GeneCounts   --clip3pNbases $trim3_g    --clip5pNbases $trim5_g     --outFileNamePrefix  $Trim_STAR_g/$temp.   --genomeDir $STAR_index_g  --readFilesIn $input2_g/$end1.fastq  $input2_g/$end2.fastq   >>$Trim_STAR_g/$temp.runLog  2>&1 ");    
        system("rename  s/Aligned.out.sam/sam/   $Trim_STAR_g/*Aligned.out.sam");                                                                                                                        
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {   
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("STAR  --runMode alignReads   --runThreadN $numCores_g    --sjdbGTFfile 0-Other/Shortcuts/mm10_RefSeq_GTF  --quantMode   GeneCounts    --clip3pNbases $trim3_g    --clip5pNbases $trim5_g      --outFileNamePrefix  $Trim_STAR_g/$temp.   --genomeDir $STAR_index_g  --readFilesIn $input2_g/$temp.fastq   >>$Trim_STAR_g/$temp.runLog  2>&1 ");    
        system("rename  s/Aligned.out.sam/sam/   $Trim_STAR_g/*Aligned.out.sam");                                                            
}
}  ########## End STAR
&myQC_BAM_1($Trim_STAR_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $HISAT2_g   = "$output_g/7_HISAT2";  
&myMakeDir($HISAT2_g);
{ ########## Start HISAT2
say   "\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using HISAT2 ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {  
        say    "\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$HISAT2_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("hisat2         --threads $numCores_g   -q   --phred33   --end-to-end    -x $HISAT2_index_g    -1 $input_g/$end1.fastq        -2 $input_g/$end2.fastq     -S $HISAT2_g/$temp.sam    >>$HISAT2_g/$temp.runLog  2>&1");                                                                                                                                 
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {   
        say   "\t......$singleEnd_g[$i]";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("hisat2        --threads $numCores_g   -q   --phred33   --end-to-end    -x $HISAT2_index_g    -U $input_g/$temp.fastq                                    -S $HISAT2_g/$temp.sam    >>$HISAT2_g/$temp.runLog  2>&1");                                                                
}
}  ########## End HISAT2
&myQC_BAM_1($HISAT2_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $Trim_HISAT2_g   = "$output_g/8_Trim_HISAT2";  
&myMakeDir($Trim_HISAT2_g);
{ ########## Start HISAT2
say   "\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using Trim_HISAT2 ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {  
        say    "\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$Trim_HISAT2_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("hisat2         --threads $numCores_g   -q   --phred33   --end-to-end    -x $HISAT2_index_g   --trim5 $trim5_g  --trim3 $trim3_g   -1 $input2_g/$end1.fastq        -2 $input2_g/$end2.fastq     -S $Trim_HISAT2_g/$temp.sam    >>$Trim_HISAT2_g/$temp.runLog  2>&1");                                                                                                                                 
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {   
        say   "\t......$singleEnd_g[$i]";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("hisat2        --threads $numCores_g   -q   --phred33   --end-to-end    -x $HISAT2_index_g   --trim5 $trim5_g  --trim3 $trim3_g     -U $input2_g/$temp.fastq                                    -S $Trim_HISAT2_g/$temp.sam    >>$Trim_HISAT2_g/$temp.runLog  2>&1");                                                                
}
}  ########## End HISAT2
&myQC_BAM_1($Trim_HISAT2_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $RSEM2_g   = "$output_g/9_RSEM_raw";  
&myMakeDir($RSEM2_g);
{ ########## Start RSEM
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference transcriptome by using RSEM ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {  
        say    "\n\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]\n";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$RSEM2_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("rsem-calculate-expression  -p $numCores_g  --paired-end      --bowtie2    --estimate-rspd   --append-names    $input2_g/$end1.fastq  $input2_g/$end2.fastq   $RSEM_index_g    $RSEM2_g/$temp    >>$RSEM2_g/$temp.runLog  2>&1 ");    
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {   
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("rsem-calculate-expression  -p $numCores_g  --fragment-length-mean 250    --fragment-length-sd 100     --bowtie2    --estimate-rspd   --append-names    $input2_g/$temp.fastq     $RSEM_index_g    $RSEM2_g/$temp    >>$RSEM2_g/$temp.runLog  2>&1 ");    
}
}  ########## End RSEM
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $Kallisto2_g   = "$output_g/10_Kallisto_raw";  
&myMakeDir($Kallisto2_g);
{ ########## Start Kallisto
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference transcriptome by using Kallisto ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {  
        say    "\n\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]\n";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$Kallisto2_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("kallisto quant  --threads=$numCores_g       --index=$Kallisto_index_g    --output-dir=$Kallisto2_g/$temp    $input2_g/$end1.fastq  $input2_g/$end2.fastq    >>$Kallisto2_g/$temp.runLog  2>&1 ");    
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {   
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("kallisto quant  --threads=$numCores_g   --single -l 250 -s 100    --index=$Kallisto_index_g    --output-dir=$Kallisto2_g/$temp    $input2_g/$temp.fastq    >>$Kallisto2_g/$temp.runLog  2>&1");    
}
}  ########## End Kallisto
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $Sailfish2_g   = "$output_g/11_Sailfish_raw";  
&myMakeDir($Sailfish2_g);
{ ########## Start Sailfish
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference transcriptome by using Sailfish ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {  
        say    "\n\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]\n";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$Sailfish2_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("sailfish quant  --threads $numCores_g    --libType IU      --index $Sailfish_index_g    --output $Sailfish2_g/$temp    --mates1 $input2_g/$end1.fastq  --mates2 $input2_g/$end2.fastq    >>$Sailfish2_g/$temp.runLog  2>&1 ");    
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {   
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("sailfish quant  --threads $numCores_g    --libType U      --index $Sailfish_index_g    --output $Sailfish2_g/$temp    --unmatedReads $input2_g/$temp.fastq      >>$Sailfish2_g/$temp.runLog  2>&1 ");                         
}
}  ########## End Sailfish
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $Salmon2_g   = "$output_g/12_Salmon_raw";  
&myMakeDir($Salmon2_g);
{ ########## Start Salmon
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference transcriptome by using Salmon ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {  
        say    "\n\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]\n";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$Salmon2_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("salmon quant  --threads $numCores_g   --libType IU      --index $Salmon_index_g    --output $Salmon2_g/$temp    --mates1 $input2_g/$end1.fastq  --mates2 $input2_g/$end2.fastq    >>$Salmon2_g/$temp.runLog  2>&1 ");    
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {   
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("salmon quant  --threads $numCores_g    --libType U      --index $Salmon_index_g    --output $Salmon2_g/$temp    --unmatedReads $input2_g/$temp.fastq      >>$Salmon2_g/$temp.runLog  2>&1 ");                         
}
}  ########## End Salmon
###################################################################################################################################################################################################






###################################################################################################################################################################################################
&myQC_BAM_2($STAR_g);
&myQC_BAM_2($Trim_STAR_g);
&myQC_BAM_2($HISAT2_g);
&myQC_BAM_2($Trim_HISAT2_g);
&myQC_BAM_3($STAR_g);
&myQC_BAM_3($Trim_STAR_g);
&myQC_BAM_3($HISAT2_g);
&myQC_BAM_3($Trim_HISAT2_g);
&myQC_BAM_4($STAR_g);
&myQC_BAM_4($Trim_STAR_g);
&myQC_BAM_4($HISAT2_g);
&myQC_BAM_4($Trim_HISAT2_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n";





## END
