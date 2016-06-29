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

        Step 6: Remove unpaired reads and sort the BAM based on read name by using SAMtools. (Only for paired-end reads)
        Usage:  
               perl  RASDA6.pl    [-version]    [-help]    [-in inputDir]    [-out outDir]   
        For instance: 
               perl  RASDA6.pl    -in 6-FinalBAM/5_STAR     -out 7-SortName/5_STAR   >> RASDA6.runLog  2>&1
     
        -------------------------------------------------------------------------------------------------------------------
        Optional arguments:
        -v, --version        Show version number of this program and exit.

        -h, --help           Show this help message and exit.

        Required arguments:
        -in  inputDir       inputDir is the name of your input folder that contains your BAM files, the suffix of the BAM files must be ".bam".  (no default)

        -out outDir         outDir is the name of your output folder that contains running results (BAM format) of this step.  (no default)
        ------------------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines such as CISDA, MESDA and HISDA,
        please see the web site: https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies 
        and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.    
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------  
';

## Version Information
my $version_g = "  The Sixth Step of RASDA (RNA-Seq Data Analyzer), version 0.7.3, 2016-06-01.";

## Keys and Values
if ($#ARGV   == -1) { print  "\n$HELP_g\n\n";  exit 0; }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0) { @ARGV = (@ARGV, "help");         }       ## when the number of command argumants is odd. 
my %args = @ARGV;

## Initialize  Variables
my $input_g   = '6-FinalBAM/5_STAR';      ## This is only an initialization  value or suggesting value, not default value.
my $output_g  = '7-SortName/5_STAR';      ## This is only an initialization  value or suggesting value, not default value.

## Available Arguments
my $available = "  -version    -help    -in    -out    ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if($available !~ m/\s$key\s/) {print  "    Cann't recognize $key !!\n";  $boole_g = 1; }
}
if($boole_g == 1) {
    print "\n    The Command Line Arguments are wrong!\n";
    print   '    Please see help message by using "perl  RASDA6.pl  -help". ';
    print "\n\n";
    exit 0;
}

## Get Arguments
if ( exists $args{'-version' }   )     { say  "\n$version_g\n";   exit 0; }
if ( exists $args{'-help'    }   )     { say  "\n$HELP_g\n";      exit 0; }
if ( exists $args{'-in'      }   )     { $input_g  = $args{'-in' };       }else{say   "\n -in  is required.\n";          say  "\n$HELP_g\n";    exit 0; }
if ( exists $args{'-out'     }   )     { $output_g = $args{'-out'};       }else{say   "\n -out is required.\n";          say  "\n$HELP_g\n";    exit 0; }

## Conditions
$input_g   =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";
$output_g  =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";

## Print Command Arguments to Standard Output
print  "\n\n
        ################ Arguments ###############################
                Input  folder:  $input_g
                Output folder:  $output_g
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
my $pattern_g  = "[-.0-9A-Za-z]+";
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
        next unless $inputFiles_g[$i] !~ m/^removed_/;
        say   "\t......$inputFiles_g[$i]" ; 
        my $temp = $inputFiles_g[$i]; 
        $groupFiles[++$#groupFiles] = $inputFiles_g[$i];  
        $temp =~ m/^(\d+)_($pattern_g)_(Rep[1-9])/   or  die   "wrong-1: ## $temp ##";
        $temp =~ m/_(Rep[1-9])\.bam$/  or    die   "wrong-2: ## $temp ##";
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
sub  myQC_BAM_1  {
    my $dir1      =  $_[0];   ## All the BAM files must be in this folder.
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
        next unless $Files[$i] =~ m/\.bam$/;
        next unless $Files[$i] !~ m/^[.]/;
        next unless $Files[$i] !~ m/[~]$/;
        my $temp = $Files[$i];
        say    "\t......$temp";
        $temp =~ s/\.bam$//  ||  die;
        system("samtools  index           $dir1/$temp.bam      >>$SAMtools/$temp.runLog  2>&1");
        system("samtools  flagstat        $dir1/$temp.bam      >>$SAMtools/$temp.runLog  2>&1");
        system(`samtools  idxstats        $dir1/$temp.bam      >>$SAMtools/$temp.runLog  2>&1`);
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
    my $BAMStats  = "$QCresults/9_BAMStats";
    my $ezBAMQC   = "$QCresults/10_ezBAMQC";
    &myMakeDir($QCresults);
    &myMakeDir($PRESEQ);
    &myMakeDir($PicardDir);
    &myMakeDir($MultiQC1);
    &myMakeDir($MultiQC2);
    &myMakeDir($BAMStats);
    &myMakeDir($ezBAMQC);
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

        system("java  -jar   $Picard_g   CollectAlignmentSummaryMetrics      INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/1_CollectAlignmentSummaryMetrics     R=0-Other/Shortcuts/mm10.fa                            >> $PicardDir/$temp/1.runLog   2>&1" );
        system("java  -jar   $Picard_g   EstimateLibraryComplexity           INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/2_EstimateLibraryComplexity                                                                      >> $PicardDir/$temp/2.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectInsertSizeMetrics            INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/3_CollectInsertSizeMetrics          HISTOGRAM_FILE=$PicardDir/$temp/3.pdf  MINIMUM_PCT=0.01      >> $PicardDir/$temp/3.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectJumpingLibraryMetrics        INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/4_CollectJumpingLibraryMetrics                                                                   >> $PicardDir/$temp/4.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectMultipleMetrics              INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/5_CollectMultipleMetrics                                                                         >> $PicardDir/$temp/5.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectBaseDistributionByCycle      INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/6_CollectBaseDistributionByCycle     CHART_OUTPUT=$PicardDir/$temp/6.pdf                         >> $PicardDir/$temp/6.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectQualityYieldMetrics          INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/7_CollectQualityYieldMetrics                                                                     >> $PicardDir/$temp/7.runLog   2>&1" ); 
        system("java  -jar   $Picard_g   CollectWgsMetricsFromQuerySorted    INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/8_CollectWgsMetricsFromQuerySorted                                                               >> $PicardDir/$temp/8.runLog   2>&1" );
        system("java  -jar   $Picard_g   MeanQualityByCycle                  INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/9_MeanQualityByCycle                 CHART_OUTPUT=$PicardDir/$temp/9.pdf                         >> $PicardDir/$temp/9.runLog   2>&1" );
        system("java  -jar   $Picard_g   QualityScoreDistribution            INPUT=$dir1/$temp.bam   OUTPUT=$PicardDir/$temp/10_QualityScoreDistribution          CHART_OUTPUT=$PicardDir/$temp/10.pdf                        >> $PicardDir/$temp/10.runLog  2>&1" ); 

        system("java  -jar  $BAMStats_g    --distances  --infile $dir1/$temp.bam   --lengths   --mapped    --outfile $BAMStats/$temp   --qualities   --starts   --view html   >> $BAMStats/$temp.runLog   2>&1");  
        system("ezBAMQC    -i $dir1/$temp.bam   --refgene 0-Other/Shortcuts/mm10_RefSeq_GTF   --outputDir $ezBAMQC/$temp  --stranded no  --label $temp  --threads $numCores_g      >> $ezBAMQC/$temp.runLog   2>&1");  
    }
    system( "multiqc  --verbose  --outdir $MultiQC1          $PRESEQ/*                 >> $MultiQC1/MultiQC.PRESEQ.runLog   2>&1" );
    system( "multiqc  --verbose  --outdir $MultiQC2          $PicardDir/*              >> $MultiQC2/MultiQC.Picard.runLog   2>&1" );
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_BAM_3  {
    my $dir1      =  $_[0];   ## All the BAM files must be in this folder.
    my $QCresults = "$dir1/QC_Results";
    my $SubreadUti= "$QCresults/11_SubreadUti";
    my $Fastqp    = "$QCresults/12_fastqp";
    my $BamQC     = "$QCresults/13_BamQC";
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
say   "\n\n\n\n\n\n##################################################################################################";
say   "Sort all the reads based on read name ......";
for (my $i=0; $i<=$#BAMfiles_g; $i++) {
    my $temp = $BAMfiles_g[$i]; 
    $temp =~ s/\.bam$//  ||  die; 
    say   "\t......$BAMfiles_g[$i]";
    system("samtools  view  -hb   -f 2   -o $output_g/$temp.t1.sam   --output-fmt sam                          --threads $numCores_g     $input_g/$temp.bam        >>$output2_g/$temp.runLog    2>&1");
    system("samtools  sort  -m 2G  -n    -o $output_g/$temp.t2.bam   --output-fmt bam  -T $output_g/yp_$temp   --threads $numCores_g     $output_g/$temp.t1.sam    >>$output2_g/$temp.runLog    2>&1");
    system("samtools fixmate  $output_g/$temp.t2.bam     $output_g/$temp.bam     >>$output2_g/$temp.runLog    2>&1");   
    system("rm  $output_g/$temp.t1.sam");
    system("rm  $output_g/$temp.t2.bam");
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
&myQC_BAM_1($output_g);
&myQC_BAM_2($output_g);
&myQC_BAM_3($output_g); 
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n";





## END
