#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.18;       
## perl5 version >= 5.18,   you can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.





###################################################################################################################################################################################################
###################################################################################################################################################################################################


########## Help Infromation ##########
my $HELP_g = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use CISDA (ChIP-Seq Data Analyzer), version 0.7.1, 2016-04-11.      
        CISDA is a Pipeline for Single-end and Paired-end ChIP-Seq Data Analysis by Integrating Lots of Softwares.

        Step 5: Potential PCR duplicates are removed using "MarkDuplicates" from "Picards".
                Only the UMNRR4 with MAPQ is more than 30 are kept.
                UMNRR4: Uniquely mapped non-redundant reads with mismatch is no more than 4%. 
                Assess the quality of ChIP-Seq reads (in BAM files) to identify possible sequencing errors or biases by using 20 softwares:
                	SAMtools, FastQC, bamtools stats, bam stats (BamUtil), samstat, BamQC (Simon Andrews), PRESEQ, NGSQC, qualimap, fastqp, htseq-qa,   
                	multiqc, BBMap, qc3.pl, bamutils stats in ngsutils, biobambam2, Subread utilities, 10 tools in Picard, QuasR and Rqc.                 

        Usage:  
               perl  CISDA-5.pl    [-v]    [-h]    [-i inputDir]    [-o outDir]   
        For instance: 
                     perl  CISDA-5.pl    -i 5-MAPQ30/1-BWAmem          -o 6-FinalBAM/1-BWAmem         
                     perl  CISDA-5.pl    --input 5-MAPQ30/1-BWAmem     --output 6-FinalBAM/1-BWAmem   
                     perl  CISDA-5.pl    --input 5-MAPQ30/1-BWAmem     --output 6-FinalBAM/1-BWAmem      >> CISDA-5.runLog  2>&1
                           CISDA-5.pl    --input 5-MAPQ30/1-BWAmem     --output 6-FinalBAM/1-BWAmem      >> CISDA-5.runLog  2>&1
     
        -------------------------------------------------------------------------------------------------------------------
        Optional arguments:
        -v, --version        Show version number of this program and exit.

        -h, --help           Show this help message and exit.

        Required arguments:
        -i inputDir,  --input inputDir        inputDir is the name of your input folder that contains your BAM files,
                                              the suffix of the BAM files must be ".bam".    (no default)

        -o outDir,  --output outDir           outDir is the name of your output folder that contains running results.                                                 
        ------------------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines such as RASDA, MESDA and HISDA,
        please see the web site: https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies 
        and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.    
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------  
';


########## Version Infromation ##########
my $version_g = "  The Fifth Step of CISDA (ChIP-Seq Data Analyzer), version 0.7.1, 2016-04-11.";


########## Keys and Values ##########
if ($#ARGV   == -1) { say  "\n$HELP_g\n";  exit 0; }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0) { @ARGV = (@ARGV, "-h");       }       ## when the number of command argumants is odd. 
my %args = @ARGV;


########## Initialize  Variables ##########
my $input_g  = '5-MAPQ30/1-BWAmem';         ## This is only an initialization  value or suggesting value, not default value.
my $output_g = '6-FinalBAM/1-BWAmem';       ## This is only an initialization  value or suggesting value, not default value.


########## Available Arguments ##########
my $available = "  -v  --version    -h  --help    -i  --input    -o    --output     ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if($available !~ m/\s$key\s/) {say    "\n\tCann't recognize $key !!";  $boole_g = 1; }
}
if($boole_g == 1) {
    say   "\tThe Command Line Arguments are wrong!";
    say   "\tPlease see help message by using 'perl  CISDA-5.pl  -h' \n";
    exit 0;
}


########## Get Arguments ##########
if ( ( exists $args{'-v' } )  or  ( exists $args{'--version'      } )  )     { say  "\n$version_g\n";    exit 0; }
if ( ( exists $args{'-h' } )  or  ( exists $args{'--help'         } )  )     { say  "\n$HELP_g\n";       exit 0; }
if ( ( exists $args{'-i' } )  or  ( exists $args{'--input'        } )  )     { ($input_g   = $args{'-i'  })  or  ($input_g   = $args{'--input'      });  }else{say   "\n -i  or --input  is required.\n";          say  "\n$HELP_g\n";       exit 0; }                                               
if ( ( exists $args{'-o' } )  or  ( exists $args{'--output'       } )  )     { ($output_g  = $args{'-o'  })  or  ($output_g  = $args{'--output'     });  }else{say   "\n -o  or --output is required.\n";          say  "\n$HELP_g\n";       exit 0; }      


########### Conditions #############
$input_g   =~ m/^\S+$/   ||  die   "The Command Line Arguments are wrong!\n$HELP_g\n\n";
$output_g  =~ m/^\S+$/   ||  die   "The Command Line Arguments are wrong!\n$HELP_g\n\n";


######### say Command Arguments to Standard Output ###########
say  "\n
        ################ Your Arguments ###############################
                Input  folder:  $input_g
                Output folder:  $output_g
        ###############################################################  
\n";


###################################################################################################################################################################################################
###################################################################################################################################################################################################










say    "\n\n\n\n\n\n##################################################################################################";
say    "Running......";
sub myMakeDir  {
    my $path = $_[0];
    if ( !( -e $path) )  { mkdir $path  ||  die; }
}
my $output2_g = "$output_g/QC_Results";
system("mkdir  -p  $output_g");
&myMakeDir($output_g);
&myMakeDir($output2_g);
opendir(my $DH_input, $input_g)  ||  die;     
my @inputFiles = readdir($DH_input);
my $pattern = "[-.0-9A-Za-z]+";
my $numCores = 4;





say   "\n\n\n\n\n\n##################################################################################################";
say   "Checking all the necessary softwares in this step......" ;
sub printVersion  {
    my $software = $_[0];
    system("echo    '##############################################################################'  >> $output2_g/VersionsOfSoftwares.txt   2>&1");
    system("echo    '#########$software'                                                              >> $output2_g/VersionsOfSoftwares.txt   2>&1");
    system("$software                                                                                 >> $output2_g/VersionsOfSoftwares.txt   2>&1");
    system("echo    '\n\n\n\n\n\n'                                                                    >> $output2_g/VersionsOfSoftwares.txt   2>&1");
}
my  $Picard = "/home/yp/.MyProgramFiles/3_HTS-2G/2-BAMtools/picard-tools-2.1.1/picard.jar";
&printVersion("samtools");
&printVersion("fastqc   -v");
&printVersion("bamtools stats  -h");
&printVersion("bam   stats");
&printVersion("samstat   -v");
&printVersion("bamqc  -v");
&printVersion("preseq");
&printVersion("NGSQC  -h");
&printVersion("qualimap  -v");
&printVersion("fastqp   -h");
&printVersion("htseq-qa   -h");
&printVersion("multiqc   --version");
&printVersion("qc3.pl -h");
&printVersion("bamutils stats");
&printVersion("bammapdist  -h");   ##biobambam2
&printVersion("propmapped");
&printVersion("qualityScores");
&printVersion("commonkmers.sh");   ##BBMap
&printVersion("countgc.sh");       ##BBMap
&printVersion("java  -jar  $Picard   MarkDuplicates                      --version");
&printVersion("java  -jar  $Picard   CollectAlignmentSummaryMetrics      --version");
&printVersion("java  -jar  $Picard   EstimateLibraryComplexity           --version");
&printVersion("java  -jar  $Picard   CollectInsertSizeMetrics            --version");
&printVersion("java  -jar  $Picard   CollectJumpingLibraryMetrics        --version");
&printVersion("java  -jar  $Picard   CollectMultipleMetrics              --version");
&printVersion("java  -jar  $Picard   CollectBaseDistributionByCycle      --version");
&printVersion("java  -jar  $Picard   CollectQualityYieldMetrics          --version");
&printVersion("java  -jar  $Picard   CollectWgsMetricsFromQuerySorted    --version");
&printVersion("java  -jar  $Picard   MeanQualityByCycle                  --version");
&printVersion("java  -jar  $Picard   QualityScoreDistribution            --version");



 
                 
say   "\n\n\n\n\n\n##################################################################################################";
say   "Checking all the input file names ......";
my @groupFiles = ();     
my $fileNameBool = 1;
for ( my $i=0; $i<=$#inputFiles; $i++ ) {   
        next unless $inputFiles[$i] =~ m/\.bam$/;  
        next unless $inputFiles[$i] !~ m/^[.]/;
        next unless $inputFiles[$i] !~ m/[~]$/;
        next unless $inputFiles[$i] !~ m/^QC_Results$/;
        next unless $inputFiles[$i] !~ m/^removed_/;
        say   "\t......$inputFiles[$i]" ; 
        my $temp = $inputFiles[$i]; 
        $groupFiles[++$#groupFiles] = $inputFiles[$i];  
        $temp =~ m/^(\d+)_($pattern)_(Rep[1-9])/   or  die   "wrong-1: ## $temp ##";
        $temp =~ m/_(Rep[1-9])\.bam$/  or    die   "wrong-2: ## $temp ##";
        if($temp !~ m/^((\d+)_($pattern)_(Rep[1-9]))(_[1-2])?\.bam$/) {
             $fileNameBool = 0;
        }
}
if($fileNameBool == 1)  { say    "\n\t\tAll the file names are passed.\n";  }
@groupFiles   = sort(@groupFiles);
my $numGroup  = 0;
my $noteGroup = 0;
for ( my $i=0; $i<=$#groupFiles; $i++ ) { 
    $groupFiles[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])/  or  die;
    my $n1 = $1;
    $n1>=1  or  die;
    if($noteGroup != $n1) {say "\n\t\tGroup $n1:";  $numGroup++; }
    say  "\t\t\t$groupFiles[$i]";
    $noteGroup = $n1; 
}
say  "\n\t\tThere are $numGroup groups.";





say   "\n\n\n\n\n\n##################################################################################################";
say   "Detecting BAM files in input folder ......";
my @BAMfiles = ();
open(seqFiles_FH, ">", "$output2_g/BAM-Files.txt")  or  die; 
for ( my $i=0; $i<=$#inputFiles; $i++ ) {     
    next unless $inputFiles[$i] =~ m/\.bam$/;
    next unless $inputFiles[$i] !~ m/^[.]/;
    next unless $inputFiles[$i] !~ m/[~]$/;
    next unless $inputFiles[$i] !~ m/^removed_/;
    say    "\t......$inputFiles[$i]"; 
    $inputFiles[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])\.bam$/  or  die;  
    $BAMfiles[$#BAMfiles+1] =  $inputFiles[$i];
    say   "\t\t\t\tBAM file:  $inputFiles[$i]\n";
    say   seqFiles_FH  "BAM file: $inputFiles[$i]\n";

}
say   seqFiles_FH  "\n\n\n\n\n";  
say   seqFiles_FH  "All BAM files:@BAMfiles\n\n\n";
say    "\t\t\t\tAll BAM files:@BAMfiles\n\n";
my $num1 = $#BAMfiles + 1;
say seqFiles_FH   "\nThere are $num1 BAM files.\n";
say     "\t\t\t\tThere are $num1 BAM files.\n";





sub  myQC1  {
       my $relativePath1= $_[0];
       my $relativePath = "$relativePath1/QC_Results";
       my $SAMtools     = "$relativePath/1_SAMtools";
       my $SubreadUti   = "$relativePath/2_SubreadUti";
       my $FastQC       = "$relativePath/3_FastQC";
       my $samstat      = "$relativePath/4_samstat";
       my $MultiQC      = "$relativePath/5_MultiQC";
       my $PRESEQ       = "$relativePath/6_PRESEQ";
       my $qualimap     = "$relativePath/7_qualimap";
       my $PicardDir    = "$relativePath/8_Picard";
       &myMakeDir("$relativePath");
       &myMakeDir("$SAMtools");
       &myMakeDir("$FastQC");
       &myMakeDir("$samstat");
       &myMakeDir("$PRESEQ");
       &myMakeDir("$qualimap");
       &myMakeDir("$MultiQC");
       &myMakeDir("$SubreadUti");
       &myMakeDir("$PicardDir");
       opendir(my $DH_map, $relativePath1) || die;     
       my @mapFiles = readdir($DH_map);

       say   "\n\n\n\n\n\nDetecting the quality of bam files by using SAMtools, FASTQC, samstat, PRESEQ, qualimap, Subread utilities, Picard and MultiQC ......";
       for (my $i=0; $i<=$#mapFiles; $i++) {
           next unless $mapFiles[$i] =~ m/\.bam$/;
           next unless $mapFiles[$i] !~ m/^[.]/;
           next unless $mapFiles[$i] !~ m/[~]$/;
           my $temp = $mapFiles[$i]; 
           $temp =~ s/\.t\.bam$//  ||  die; 
           say   "\t......$mapFiles[$i]";  
           system("samtools  sort  -m 2G  -o $relativePath1/$temp.bam   --output-fmt bam  -T $relativePath1/1_$temp   --threads $numCores    $relativePath1/$temp.t.bam    >>$SAMtools/$temp.runLog    2>&1");                                                                                                                                                    
           system("samtools  index           $relativePath1/$temp.bam      >>$SAMtools/$temp.runLog  2>&1");
           system("samtools  flagstat        $relativePath1/$temp.bam      >>$SAMtools/$temp.runLog  2>&1");
           system(`samtools  idxstats        $relativePath1/$temp.bam      >>$SAMtools/$temp.runLog  2>&1`);
           system("rm   $relativePath1/$temp.t.bam"); 
           system("fastqc    --outdir $FastQC     --threads $numCores    --format bam    --kmers 7     $relativePath1/$temp.bam   >>$FastQC/$temp.runLog        2>&1"); 
           system("samstat   $relativePath1/$temp.bam      >> $samstat/$temp.runLog         2>&1");   
           system("qualimap  bamqc  -bam $relativePath1/$temp.bam   -c  -nt $numCores  -outdir $qualimap/$temp   --java-mem-size=12G   >>$qualimap/$temp.runLog    2>&1");
           system("propmapped   -i $relativePath1/$temp.bam                    -o $SubreadUti/$temp.prommapped      >> $SubreadUti/$temp.prommapped      2>&1");
           system("echo      '\n\n\n\n\n'                                                                           >> $SubreadUti/$temp.prommapped      2>&1"); 
           system("propmapped   -i $relativePath1/$temp.bam       -f           -o $SubreadUti/$temp.prommapped      >> $SubreadUti/$temp.prommapped      2>&1"); 
           system("echo      '\n\n\n\n\n'                                                                           >> $SubreadUti/$temp.prommapped      2>&1"); 
           system("propmapped   -i $relativePath1/$temp.bam       -f   -p      -o $SubreadUti/$temp.prommapped      >> $SubreadUti/$temp.prommapped      2>&1"); 
           system("qualityScores   --BAMinput   -i $relativePath1/$temp.bam    -o $SubreadUti/$temp.qualityScores   >> $SubreadUti/$temp.qualityScores   2>&1");
           system("preseq  c_curve   -output  $PRESEQ/$temp.pe.txt     -step 1000000    -verbose   -pe  -bam  $relativePath1/$temp.bam    >> $PRESEQ/$temp.pe.runLog   2>&1");   
           system("preseq  c_curve   -output  $PRESEQ/$temp.se.txt     -step 1000000    -verbose        -bam  $relativePath1/$temp.bam    >> $PRESEQ/$temp.se.runLog   2>&1");  
           &myMakeDir("$PicardDir/$temp"); 
           #system("java  -jar   $Picard   CollectAlignmentSummaryMetrics      INPUT=$relativePath1/$temp.bam   OUTPUT=$PicardDir/$temp/1_CollectAlignmentSummaryMetrics     R=0-Other/ShortCuts/$genome_g.fa                            >> $PicardDir/$temp/1.runLog   2>&1" );
           system("java  -jar   $Picard   EstimateLibraryComplexity           INPUT=$relativePath1/$temp.bam   OUTPUT=$PicardDir/$temp/2_EstimateLibraryComplexity                                                                      >> $PicardDir/$temp/2.runLog   2>&1" );
           system("java  -jar   $Picard   CollectInsertSizeMetrics            INPUT=$relativePath1/$temp.bam   OUTPUT=$PicardDir/$temp/3_CollectInsertSizeMetrics          HISTOGRAM_FILE=$PicardDir/$temp/3.pdf  MINIMUM_PCT=0.01      >> $PicardDir/$temp/3.runLog   2>&1" );
           system("java  -jar   $Picard   CollectJumpingLibraryMetrics        INPUT=$relativePath1/$temp.bam   OUTPUT=$PicardDir/$temp/4_CollectJumpingLibraryMetrics                                                                   >> $PicardDir/$temp/4.runLog   2>&1" );
           system("java  -jar   $Picard   CollectMultipleMetrics              INPUT=$relativePath1/$temp.bam   OUTPUT=$PicardDir/$temp/5_CollectMultipleMetrics                                                                         >> $PicardDir/$temp/5.runLog   2>&1" );
           system("java  -jar   $Picard   CollectBaseDistributionByCycle      INPUT=$relativePath1/$temp.bam   OUTPUT=$PicardDir/$temp/6_CollectBaseDistributionByCycle     CHART_OUTPUT=$PicardDir/$temp/6.pdf                         >> $PicardDir/$temp/6.runLog   2>&1" );
           system("java  -jar   $Picard   CollectQualityYieldMetrics          INPUT=$relativePath1/$temp.bam   OUTPUT=$PicardDir/$temp/7_CollectQualityYieldMetrics                                                                     >> $PicardDir/$temp/7.runLog   2>&1" ); 
           system("java  -jar   $Picard   CollectWgsMetricsFromQuerySorted    INPUT=$relativePath1/$temp.bam   OUTPUT=$PicardDir/$temp/8_CollectWgsMetricsFromQuerySorted                                                               >> $PicardDir/$temp/8.runLog   2>&1" );
           system("java  -jar   $Picard   MeanQualityByCycle                  INPUT=$relativePath1/$temp.bam   OUTPUT=$PicardDir/$temp/9_MeanQualityByCycle                 CHART_OUTPUT=$PicardDir/$temp/9.pdf                         >> $PicardDir/$temp/9.runLog   2>&1" );
           system("java  -jar   $Picard   QualityScoreDistribution            INPUT=$relativePath1/$temp.bam   OUTPUT=$PicardDir/$temp/10_QualityScoreDistribution          CHART_OUTPUT=$PicardDir/$temp/10.pdf                        >> $PicardDir/$temp/10.runLog  2>&1" ); 
       }
       system( "multiqc  --outdir $MultiQC     $FastQC/*_fastqc.zip      >> $MultiQC/multiqc.fastqc.runLog   2>&1" );
}  ## End myQC1



               	  
                	 
sub  myQC2  {
       my $relativePath1= $_[0];
       my $relativePath = "$relativePath1/QC_Results";
       my $bamtools    = "$relativePath/9_bamtools";
       my $BamUtil     = "$relativePath/10_BamUtil";
       my $BamQC       = "$relativePath/11_BamQC";
       my $NGSQC       = "$relativePath/12_NGSQC";
       my $htseqQA     = "$relativePath/13_htseqQA";
       my $Fastqp      = "$relativePath/14_fastqp";
       my $QC3         = "$relativePath/15_QC3";
       my $NGSutils    = "$relativePath/16_NGSutils";
       my $BBMap       = "$relativePath/17_BBMap";
       my $biobambam2  = "$relativePath/18_biobambam2";
       my $Rqc         = "$relativePath/19_Rqc";
       my $QuasR       = "$relativePath/20_QuasR";
       &myMakeDir("$relativePath");
       &myMakeDir("$bamtools");
       &myMakeDir("$BamUtil");
       &myMakeDir("$BamQC");
       &myMakeDir("$NGSQC");
       &myMakeDir("$htseqQA");
       &myMakeDir("$Fastqp");
       &myMakeDir("$QC3");
       &myMakeDir("$NGSutils");
       &myMakeDir("$BBMap");
       &myMakeDir("$biobambam2");
       &myMakeDir("$Rqc");
       &myMakeDir("$QuasR");      
       opendir(my $DH_map, $relativePath1) || die;     
       my @mapFiles = readdir($DH_map);
       #bamtools stats, bam stats (BamUtil), BamQC (Simon Andrews), NGSQC, fastqp, htseq-qa, BBMap, qc3.pl, bamutils stats in ngsutils, biobambam2, QuasR and Rqc. 
       say   "\n\n\n\n\n\nDetecting the quality of bam files by using other 12 tools ......";
       for (my $i=0; $i<=$#mapFiles; $i++) {
           next unless $mapFiles[$i] =~ m/\.bam$/;
           next unless $mapFiles[$i] !~ m/^[.]/;
           next unless $mapFiles[$i] !~ m/[~]$/;
           my $temp = $mapFiles[$i]; 
           $temp =~ s/\.bam$//  ||  die; 
           say   "\t......$mapFiles[$i]";  
           system( "bamtools stats  -in $relativePath1/$temp.bam    -insert  >> $bamtools/$temp.runLog   2>&1");
           system( "bam      stats  -in $relativePath1/$temp.bam    --basic   --pBaseQC $BamUtil/$temp.pBaseQC  >> $BamUtil/$temp.runLog   2>&1");
           system( "fastqp    --nreads 20000000   --kmer 5    --output $Fastqp/$temp  --type fastq   --median-qual 30     $relativePath1/$temp.bam     >> $Fastqp/$temp.runLog     2>&1 " );
           system( "bamutils stats    --nreads 20000000   -all     $relativePath1/$temp.bam     >> $NGSutils/$temp.runLog     2>&1 " );
       }
}  ## End myQC2





say   "\n\n\n\n\n\n##################################################################################################";
say   "Removing some reads ......";
for (my $i=0; $i<=$#BAMfiles; $i++) {
    my $temp = $BAMfiles[$i]; 
    $temp =~ s/\.bam$//  ||  die; 
    say   "\t......$BAMfiles[$i]";
    system("bamutils cleancigar    $input_g/$temp.bam    $output_g/$temp.cleancigar.bam   >> $output2_g/$temp.cleancigar.runLog     2>&1 ");
    system("java  -jar   $Picard   MarkDuplicates    REMOVE_DUPLICATES=true   INPUT=$output_g/$temp.cleancigar.bam    OUTPUT=$output_g/$temp.t.bam    METRICS_FILE=$output2_g/$temp.marked_dup_metrics.txt   >> $output2_g/$temp.runLog     2>&1 ");
    system("rm  $output_g/$temp.cleancigar.bam");
}
&myQC1($output_g);





say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n";





## END
