#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.18;       
## perl5 version >= 5.18,   you can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
## Help Infromation
my $HELP_g = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use CISDA (ChIP-Seq Data Analyzer), version 0.7.2, 2016-04-17.      
        CISDA is a Pipeline for Single-end and Paired-end ChIP-Seq Data Analysis by Integrating Lots of Softwares.

        Step 3: Mapping reads to the reference genome by using 5 softwares (mappers or aligners): 
                    BWA-aln/BWA-mem, Bowtie1/Bowtie2, Trim_Bowtie2, Subread, Trim_Subread. 
                Assess the quality of BAM files to identify possible sequencing errors or biases by using 10 softwares:
                    SAMtools, Subread utilities, FASTQC, samstat, qualimap, MultiQC, PRESEQ, Picard, fastqp, BamQC.

        Usage:  
               perl  CISDA3.pl    [-version]    [-help]    [-in inputDir]    [-out outDir]   [-mismatch N]   [-60bp yes]  [-genome reference_genome]    
        For instance: 
               perl  CISDA3.pl   -in 3-Filtered    -out 4-Mapping   -mismatch 6     -60bp yes   -genome mm10   >> CISDA3.runLog  2>&1

        -------------------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.

        -help           Show this help message and exit.

        Required arguments:
        -in inputDir        inputDir is the name of your input folder that contains your FASTQ files, the suffix of the FASTQ files must be ".fastq".    (no default)

        -out outDir         outDir is the name of your output folder that contains running results (BAM format) of this step.      (no default)

        -mismatch N         Specify the maximum number of mis-matched bases allowed in the alignment.  (no default)

        -60bp booleValue    booleValue is "no"  if reads length <= 60 bp, then Bowtie1 and BWA-aln will be invoked. 
                            booleValue is "yes" if reads length >  60 bp, then Bowtie2 and BWA-mem will be invoked.   (no default)
                                                  
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

        For more details about this pipeline and other NGS data analysis piplines such as RASDA, MESDA and HISDA,
        please see the web site: https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies 
        and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.    
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------  
';

## Version Infromation 
my $version_g = "  The Third Step of CISDA (ChIP-Seq Data Analyzer), version 0.7.2, 2016-04-17.";

## Keys and Values
if ($#ARGV   == -1) { say  "\n$HELP_g\n";  exit 0; }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0) { @ARGV = (@ARGV, "-h");       }       ## when the number of command argumants is odd. 
my %args = @ARGV;

## Initialize  Variables
my $input_g  = '3-Filtered';      ## This is only an initialization  value or suggesting value, not default value.
my $output_g = '4-Mapping';       ## This is only an initialization  value or suggesting value, not default value.
my $MM_g     = 6;                 ## This is only an initialization  value or suggesting value, not default value.
my $readLen_g= 'yes';             ## This is only an initialization  value or suggesting value, not default value.
my $genome_g = 'mm10';            ## This is only an initialization  value or suggesting value, not default value.

## Available Arguments
my $available = "  -version    -help   -in   -out  -mismatch   -60bp   -genome   ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if($available !~ m/\s$key\s/) {say    "\n\tCann't recognize $key !!";  $boole_g = 1; }
}
if($boole_g == 1) {
    say   "\tThe Command Line Arguments are wrong!";
    say   "\tPlease see help message by using 'perl  CISDA3.pl  -help' \n";
    exit 0;
}

## Get Arguments
if ( exists $args{'-version' }   )     { say  "\n$version_g\n";   exit 0; }
if ( exists $args{'-help'    }   )     { say  "\n$HELP_g\n";      exit 0; }
if ( exists $args{'-in'      }   )     { $input_g  = $args{'-in' };       }else{say   "\n -in  is required.\n";          say  "\n$HELP_g\n";    exit 0; }
if ( exists $args{'-out'     }   )     { $output_g = $args{'-out'};       }else{say   "\n -out is required.\n";          say  "\n$HELP_g\n";    exit 0; }
if ( exists $args{'-mismatch'}   )     { $MM_g     = $args{'-mismatch'};  }else{say   "\n -mismatch is required.\n";     say  "\n$HELP_g\n";    exit 0; } 
if ( exists $args{'-60bp'    }   )     { $readLen_g= $args{'-60bp'    };  }else{say   "\n -60bp is required.\n";         say  "\n$HELP_g\n";    exit 0; } 
if ( exists $args{'-genome'  }   )     { $genome_g = $args{'-genome'  };  }else{say   "\n -genome is required.\n";       say  "\n$HELP_g\n";    exit 0; } 

## Conditions
$input_g   =~ m/^\S+$/   ||  die   "The Command Line Arguments are wrong!\n$HELP_g\n\n";
$output_g  =~ m/^\S+$/   ||  die   "The Command Line Arguments are wrong!\n$HELP_g\n\n";
$MM_g      =~ m/^\d+$/   ||  die   "The Command Line Arguments are wrong!\n$HELP_g\n\n";
$readLen_g =~ m/^\S+$/   ||  die   "The Command Line Arguments are wrong!\n$HELP_g\n\n";
$genome_g  =~ m/^\S+$/   ||  die   "The Command Line Arguments are wrong!\n$HELP_g\n\n";

## say Command Arguments to Standard Output
say  "\n
        ################ Arguments ###############################
                Input  folder:  $input_g
                Output folder:  $output_g
                Maximum number of mis-matched bases: $MM_g
                Read length is more than 60bp: $readLen_g
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
my  $commonPath_g         = "/home/yp/.MyProgramFiles/3_HTS2G/2_ChIPseq/2_Mapping";
my  $BWA_index_g          = "$commonPath_g/bwa-0.7.13/RefGenomes/$genome_g/$genome_g";
my  $Bowtie1_index_g      = "$commonPath_g/bowtie-1.1.2/RefGenomes/$genome_g/$genome_g";
my  $Bowtie2_index_g      = "$commonPath_g/bowtie2-2.2.8/RefGenomes/$genome_g/$genome_g";
my  $Subread_index_g      = "$commonPath_g/subread-1.5.0-p1/RefGenomes/$genome_g/$genome_g";
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
my  $Picard_g = "/home/yp/.MyProgramFiles/3_HTS2G/1_common/3_afterMap/picard-tools-2.1.1/picard.jar";
my  $trim5_g  = 2;   ##bp
my  $trim3_g  = 5;   ##bp
&printVersion("bwa  aln");
&printVersion("bwa  mem");
&printVersion("bowtie    --version");
&printVersion("bowtie2   --version");
&printVersion("subread-align  -v");
&printVersion("samtools");
&printVersion("fastqc   -v");
&printVersion("samstat   -v");
&printVersion("bamqc  -v");
&printVersion("preseq");
&printVersion("qualimap  -v");
&printVersion("fastqp   -h");
&printVersion("multiqc   --version");
&printVersion("propmapped");
&printVersion("qualityScores");
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
    say   "Detecting the quality of all BAM files by using SAMtools, FastQC, qualimap, samstat and MultiQC ......";
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
        ##system("rm   $dir1/$temp.sam"); 
        system( "fastqc    --outdir $FastQC    --threads $numCores_g  --format bam   --kmers 7    $dir1/$temp.bam              >> $FastQC/$temp.runLog      2>&1" );
        system( "qualimap  bamqc  -bam $dir1/$temp.bam   -c  -nt $numCores_g   -outdir $qualimap/$temp   --java-mem-size=12G   >> $qualimap/$temp.runLog    2>&1" );
        system( "samstat   $dir1/$temp.bam      >> $samstat/$temp.runLog         2>&1");   
    }
    system( "multiqc  --verbose  --outdir $MultiQC1          $FastQC/*_fastqc.zip      >> $MultiQC1/MultiQC.FastQC.runLog     2>&1" );
    system( "multiqc  --verbose  --outdir $MultiQC2          $qualimap/*               >> $MultiQC2/MultiQC.qualimap.runLog   2>&1" );
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_BAM_2  {
    my $dir1      =  $_[0];   ## All the BAM files must be in this folder.
    my $QCresults = "$dir1/QC_Results";
    my $PRESEQ    = "$QCresults/6_PRESEQ";
    my $PicardDir = "$QCresults/7_Picard";
    my $MultiQC1  = "$QCresults/8_MultiQC_PRESEQ";
    my $MultiQC2  = "$QCresults/8_MultiQC_Picard";
    &myMakeDir($QCresults);
    &myMakeDir($PRESEQ);
    &myMakeDir($PicardDir);
    &myMakeDir($MultiQC1);
    &myMakeDir($MultiQC2);
    opendir(my $FH_Files, $dir1) || die;     
    my @Files = readdir($FH_Files);
    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of all BAM files by using PRESEQ, Picard and MultiQC ......";
    for ( my $i=0; $i<=$#Files; $i++ ) {
        next unless $Files[$i] =~ m/\.bam$/;
        next unless $Files[$i] !~ m/^[.]/;
        next unless $Files[$i] !~ m/[~]$/;
        my $temp = $Files[$i];
        say    "\t......$temp";
        $temp =~ s/\.bam$//  ||  die;
        system("preseq  c_curve   -output  $PRESEQ/$temp.pe.PRESEQ     -step 1000000    -verbose   -pe  -bam  $dir1/$temp.bam    >> $PRESEQ/$temp.pe.runLog   2>&1");   
        system("preseq  c_curve   -output  $PRESEQ/$temp.se.PRESEQ     -step 1000000    -verbose        -bam  $dir1/$temp.bam    >> $PRESEQ/$temp.se.runLog   2>&1");  
        &myMakeDir("$PicardDir/$temp"); 
        system("java  -jar   $Picard_g   CollectAlignmentSummaryMetrics      INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/1_CollectAlignmentSummaryMetrics     R=0-Other/ShortCuts/$genome_g.fa                            >> $PicardDir/$temp/1.runLog   2>&1" );
        system("java  -jar   $Picard_g   EstimateLibraryComplexity           INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/2_EstimateLibraryComplexity                                                                      >> $PicardDir/$temp/2.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectInsertSizeMetrics            INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/3_CollectInsertSizeMetrics          HISTOGRAM_FILE=$PicardDir/$temp/3.pdf  MINIMUM_PCT=0.01      >> $PicardDir/$temp/3.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectJumpingLibraryMetrics        INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/4_CollectJumpingLibraryMetrics                                                                   >> $PicardDir/$temp/4.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectMultipleMetrics              INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/5_CollectMultipleMetrics                                                                         >> $PicardDir/$temp/5.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectBaseDistributionByCycle      INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/6_CollectBaseDistributionByCycle     CHART_OUTPUT=$PicardDir/$temp/6.pdf                         >> $PicardDir/$temp/6.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectQualityYieldMetrics          INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/7_CollectQualityYieldMetrics                                                                     >> $PicardDir/$temp/7.runLog   2>&1" ); 
        system("java  -jar   $Picard_g   CollectWgsMetricsFromQuerySorted    INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/8_CollectWgsMetricsFromQuerySorted                                                               >> $PicardDir/$temp/8.runLog   2>&1" );
        system("java  -jar   $Picard_g   MeanQualityByCycle                  INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/9_MeanQualityByCycle                 CHART_OUTPUT=$PicardDir/$temp/9.pdf                         >> $PicardDir/$temp/9.runLog   2>&1" );
        system("java  -jar   $Picard_g   QualityScoreDistribution            INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/10_QualityScoreDistribution          CHART_OUTPUT=$PicardDir/$temp/10.pdf                        >> $PicardDir/$temp/10.runLog  2>&1" ); 
    }
    system( "multiqc  --verbose  --outdir $MultiQC1          $PRESEQ/*                 >> $MultiQC1/MultiQC.PRESEQ.runLog   2>&1" );
    system( "multiqc  --verbose  --outdir $MultiQC2          $PicardDir/*              >> $MultiQC2/MultiQC.Picard.runLog   2>&1" );
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_BAM_3  {
    my $dir1      =  $_[0];   ## All the BAM files must be in this folder.
    my $QCresults = "$dir1/QC_Results";
    my $SubreadUti= "$QCresults/9_SubreadUti";
    my $Fastqp    = "$QCresults/10_fastqp";
    my $BamQC     = "$QCresults/11_BamQC";
    &myMakeDir("$QCresults");
    &myMakeDir("$SubreadUti");
    &myMakeDir("$Fastqp");
    &myMakeDir("$BamQC");
    opendir(my $DH_map, $dir1) || die;     
    my @mapFiles = readdir($DH_map);
    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of bam files by using Subreads utilities, fastqp and BamQC ......";
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
     }
     system( "bamqc  $dir1   --outdir  $BamQC        >> $BamQC/bamqc.runLog     2>&1 " );
}  
###################################################################################################################################################################################################


               	  
 

################################################################################################################################################################################################### 
my $BWA_g  = "$output_g/1_BWA";   
&myMakeDir($BWA_g);
{ ## Start BWA
if ($readLen_g  eq  "no")  {  ## Start bwa aln, shorter than 60bp
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using bwa aln ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {
        say    "\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq"   eq   $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$BWA_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("bwa aln   -n 0.04  -o 2   -t $numCores_g   -R 10   -f $BWA_g/$end1.sai      $BWA_index_g     $input_g/$end1.fastq    >>$BWA_g/$end1.runLog   2>&1");
        system("bwa aln   -n 0.04  -o 2   -t $numCores_g   -R 10   -f $BWA_g/$end2.sai      $BWA_index_g     $input_g/$end2.fastq    >>$BWA_g/$end2.runLog   2>&1");
        system("bwa sampe -a 600    -f $BWA_g/$temp.sam      $BWA_index_g     $BWA_g/$end1.sai  $BWA_g/$end2.sai     $input_g/$end1.fastq  $input_g/$end2.fastq    >>$BWA_g/$temp.runLog   2>&1"); 
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {   
        say   "\t......$singleEnd_g[$i]";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("bwa aln   -n 0.04   -o 2  -t $numCores_g    -f $BWA_g/$temp.sai    $BWA_index_g                         $input_g/$temp.fastq         >>$BWA_g/$temp.runLog   2>&1");
        system("bwa samse                                   -f $BWA_g/$temp.sam    $BWA_index_g     $BWA_g/$temp.sai    $input_g/$temp.fastq         >>$BWA_g/$temp.runLog   2>&1");
}
}else{  ## Start bwa mem, longer  than 60bp
($readLen_g  eq "yes")   or   die;
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using BWA mem ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {
        say    "\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]";   
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq"  eq  $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$BWA_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("bwa mem  -t $numCores_g  -T 0   $BWA_index_g   $input_g/$end1.fastq  $input_g/$end2.fastq    >$BWA_g/$temp.sam"); 
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {  
        say   "\t......$singleEnd_g[$i]";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("bwa mem  -t $numCores_g  -T 0   $BWA_index_g   $input_g/$temp.fastq   >$BWA_g/$temp.sam");
}
}
} ## End BWA
&myQC_BAM_1($BWA_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $Bowtie_g   = "$output_g/2_Bowtie";  
&myMakeDir($Bowtie_g); 
{ ## Start Bowtie
if ($readLen_g  eq  "no")  {
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using Bowtie ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {  
        say    "\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq"  eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$Bowtie_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("bowtie   --threads $numCores_g   -q  --chunkmbs 6400  --best  --maxins 800   --sam   --phred33-quals   $Bowtie1_index_g    -1 $input_g/$end1.fastq   -2 $input_g/$end2.fastq    $Bowtie_g/$temp.sam    >>$Bowtie_g/$temp.runLog   2>&1"); 
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {   
        say   "\t......$singleEnd_g[$i]";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("bowtie   --threads $numCores_g   -q  --chunkmbs 6400  --best   --sam   --phred33-quals  $Bowtie1_index_g    $input_g/$temp.fastq   $Bowtie_g/$temp.sam   >>$Bowtie_g/$temp.runLog   2>&1");                
}
}else{
($readLen_g  eq "yes")   or die;
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using Bowtie2 ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {  
        say    "\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$Bowtie_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("bowtie2   --maxins 800  --threads $numCores_g   -q   --phred33   --end-to-end    -x $Bowtie2_index_g    -1 $input_g/$end1.fastq        -2 $input_g/$end2.fastq     -S $Bowtie_g/$temp.sam    >>$Bowtie_g/$temp.runLog  2>&1");                                                                                                                                 
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {   
        say   "\t......$singleEnd_g[$i]";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("bowtie2    --threads $numCores_g   -q   --phred33   --end-to-end    -x $Bowtie2_index_g    -U $input_g/$temp.fastq       -S $Bowtie_g/$temp.sam    >>$Bowtie_g/$temp.runLog  2>&1");                                                                
}
}
}  ## End Bowtie
&myQC_BAM_1($Bowtie_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $Bowtie2_g   = "$output_g/3_Trim_Bowtie2";  
&myMakeDir($Bowtie2_g); 
{  ## Start Trim_Bowtie2
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using Trim_Bowtie2 ......";
my $inputDir2 = "2-FASTQ";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {  
        say    "\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$Bowtie2_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("bowtie2   --maxins 800  --threads $numCores_g   -q   --phred33   --end-to-end   --trim5 $trim5_g    --trim3 $trim3_g   -x $Bowtie2_index_g    -1 $inputDir2/$end1.fastq        -2 $inputDir2/$end2.fastq     -S $Bowtie2_g/$temp.sam    >>$Bowtie2_g/$temp.runLog  2>&1");                                                                                                                                 
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {   
        say   "\t......$singleEnd_g[$i]";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("bowtie2    --threads $numCores_g   -q   --phred33   --end-to-end   --trim5 $trim5_g    --trim3 $trim3_g      -x $Bowtie2_index_g    -U $inputDir2/$temp.fastq       -S $Bowtie2_g/$temp.sam    >>$Bowtie2_g/$temp.runLog  2>&1");                                                                
}
}  ## End Trim_Bowtie2
&myQC_BAM_1($Bowtie2_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $subread_g  = "$output_g/4_Subread";  
&myMakeDir($subread_g);
{ ## Start subread
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using Subread ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {  
        say   "\t......$pairedEnd_g[$i]";
        say   "\t......$pairedEnd_g[$i+1]";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$subread_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("subread-align  -T $numCores_g  -I 10  -B 1  -M $MM_g   --SAMoutput  -d 0  -D 800   -i $Subread_index_g   -r $input_g/$end1.fastq   -R  $input_g/$end2.fastq   -o  $subread_g/$temp.sam   -t 1  >>$subread_g/$temp.runLog  2>&1");               
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {  
        say   "\t......$singleEnd_g[$i]"; 
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die; 
        my $temp = $1; 
        system("subread-align  -T $numCores_g  -I 10  -B 1  -M $MM_g   --SAMoutput   -i $Subread_index_g    -r $input_g/$temp.fastq    -o $subread_g/$temp.sam   -t 1    >>$subread_g/$temp.runLog   2>&1");       
}
} ## End subread
&myQC_BAM_1($subread_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $subread2_g  = "$output_g/5_Trim_Subread";  
&myMakeDir($subread2_g);
{ ## Start Trim_subread
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using Trim_Subread ......";
my $inputDir2 = "2-FASTQ";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {  
        say   "\t......$pairedEnd_g[$i]";
        say   "\t......$pairedEnd_g[$i+1]";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$subread2_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("subread-align  -T $numCores_g  -I 10  -B 1  -M $MM_g   --SAMoutput  -d 0  -D 800      --trim5 $trim5_g  --trim3 $trim3_g          -i $Subread_index_g   -r $inputDir2/$end1.fastq   -R  $inputDir2/$end2.fastq   -o  $subread2_g/$temp.sam   -t 1  >>$subread2_g/$temp.runLog  2>&1");               
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {  
        say   "\t......$singleEnd_g[$i]"; 
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die; 
        my $temp = $1; 
        system("subread-align  -T $numCores_g  -I 10  -B 1  -M $MM_g   --SAMoutput                    --trim5 $trim5_g  --trim3 $trim3_g          -i $Subread_index_g    -r $inputDir2/$temp.fastq    -o $subread2_g/$temp.sam   -t 1    >>$subread2_g/$temp.runLog   2>&1");       
}
} ## End Trim_subread
&myQC_BAM_1($subread2_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
&myQC_BAM_2($BWA_g);
&myQC_BAM_2($Bowtie_g);
&myQC_BAM_2($Bowtie2_g);
&myQC_BAM_2($subread_g);
&myQC_BAM_2($subread2_g);
&myQC_BAM_3($BWA_g);
&myQC_BAM_3($Bowtie_g);
&myQC_BAM_3($Bowtie2_g);
&myQC_BAM_3($subread_g);
&myQC_BAM_3($subread2_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n";





## END
