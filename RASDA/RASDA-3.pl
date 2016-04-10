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
        Welcome to use RASDA (RNA-Seq Data Analyzer), version 0.7.1, 2016-04-11.      
        RASDA is a Pipeline for Single-end and Paired-end RNA-Seq Data Analysis by Integrating Lots of Softwares.

        Step 3: Mapping reads to the reference genome by using popular softwares. 
                15 aligners (mappers) are available: 
                            Subjunc,  GSNAP      and Kallisto.
                            HISAT2,   TopHat2    and MapSplice2.
                            STAR,     BBMAP      and mrsFAST.                            
                            RSEM,     CRAC       and ContexMap2. 
                            Salmon,   Sailfish   and RapMap. 

                Assess the quality of RNA-Seq reads (in BAM files) to identify possible sequencing errors or biases by using 20 softwares:
                	SAMtools, FastQC, RSeQC, RSeQC, samstat, BamQC (Simon Andrews), PRESEQ, NGSQC, qualimap, fastqp, QoRTs, multiqc,   
                	BBMap, qc3.pl, bamutils stats in ngsutils, biobambam2, Subread utilities, 10 tools in Picard, QuasR and Rqc.                         

        Usage:  
               perl  RASDA-3.pl    [-v]    [-h]    [-i inputDir]    [-o outDir]   [-MM N]   [-g reference_genome ]   [-60 boole]  [-PE INT]
        For instance: 
                     perl  RASDA-3.pl    -i 3-Filtered          -o 4-Mapping         -MM 6    -g mm10     
                     perl  RASDA-3.pl    --input 3-Filtered     --output 4-Mapping   --MaxMismatch 6  --genome mm10
                     perl  RASDA-3.pl    --input 3-Filtered     --output 4-Mapping   --MaxMismatch 6  --genome mm10   >> RASDA-3.runLog  2>&1
                           RASDA-3.pl    --input 3-Filtered     --output 4-Mapping   --MaxMismatch 6  --genome mm10   >> RASDA-3.runLog  2>&1
     
        -------------------------------------------------------------------------------------------------------------------
        Optional arguments:
        -v, --version        Show version number of this program and exit.

        -h, --help           Show this help message and exit.

        Required arguments:
        -i inputDir,  --input inputDir        inputDir is the name of your input folder that contains your FASTQ files,
                                              the suffix of the FASTQ files must be ".fastq".    (no default)

        -o outDir,  --output outDir           outDir is the name of your output folder that contains running 
                                              results (BAM format) of this step.      (no default)

        -MM N,      --MaxMismatch N           Specify the maximum number of mis-matched bases allowed in the alignment. 
                                              Only for Subread.  (no default)

        -g reference_genome,      --genome reference_genome           NGS reads wil be mapped to "reference_genome". 
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

                                              The full path of the reference genome must be given at line 187 of this script.  
                                              You can also add other reference genomes at line 187 of this script. (no default) 
        ------------------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines such as CISDA, MESDA and HISDA,
        please see the web site: https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies 
        and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.    
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------  
';


########## Version Infromation ##########
my $version_g = "  The Third Step of RASDA (RNA-Seq Data Analyzer), version 0.7.1, 2016-04-11.";


########## Keys and Values ##########
if ($#ARGV   == -1) { say  "\n$HELP_g\n";  exit 0; }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0) { @ARGV = (@ARGV, "-h");       }       ## when the number of command argumants is odd. 
my %args = @ARGV;


########## Initialize  Variables ##########
my $input_g  = '3-Filtered';      ## This is only an initialization  value or suggesting value, not default value.
my $output_g = '4-Mapping';       ## This is only an initialization  value or suggesting value, not default value.
my $MM_g     = 6;                 ## This is only an initialization  value or suggesting value, not default value.
my $genome_g = 'mm10';            ## This is only an initialization  value or suggesting value, not default value.


########## Available Arguments ##########
my $available = "  -v  --version    -h  --help    -i  --input    -o    --output  -MM  --MaxMismatch    -g  --genome    ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if($available !~ m/\s$key\s/) {say    "\n\tCann't recognize $key !!";  $boole_g = 1; }
}
if($boole_g == 1) {
    say   "\tThe Command Line Arguments are wrong!";
    say   "\tPlease see help message by using 'perl  RASDA-3.pl  -h' \n";
    exit 0;
}


########## Get Arguments ##########
if ( ( exists $args{'-v' } )  or  ( exists $args{'--version'      } )  )     { say  "\n$version_g\n";    exit 0; }
if ( ( exists $args{'-h' } )  or  ( exists $args{'--help'         } )  )     { say  "\n$HELP_g\n";       exit 0; }
if ( ( exists $args{'-i' } )  or  ( exists $args{'--input'        } )  )     { ($input_g   = $args{'-i'  })  or  ($input_g   = $args{'--input'      });  }else{say   "\n -i  or --input  is required.\n";          say  "\n$HELP_g\n";       exit 0; }                                               
if ( ( exists $args{'-o' } )  or  ( exists $args{'--output'       } )  )     { ($output_g  = $args{'-o'  })  or  ($output_g  = $args{'--output'     });  }else{say   "\n -o  or --output is required.\n";          say  "\n$HELP_g\n";       exit 0; }      
if ( ( exists $args{'-MM'} )  or  ( exists $args{'--MaxMismatch'  } )  )     { ($MM_g      = $args{'-MM' })  or  ($MM_g      = $args{'--MaxMismatch'});  }else{say   "\n -MM or --MaxMismatch is required.\n";     say  "\n$HELP_g\n";       exit 0; } 
if ( ( exists $args{'-g' } )  or  ( exists $args{'--genome'       } )  )     { ($genome_g  = $args{'-g'  })  or  ($genome_g  = $args{'--genome'     });  }else{say   "\n -g  or --genome is required.\n";          say  "\n$HELP_g\n";       exit 0; } 


########### Conditions #############
$input_g   =~ m/^\S+$/   ||  die   "The Command Line Arguments are wrong!\n$HELP_g\n\n";
$output_g  =~ m/^\S+$/   ||  die   "The Command Line Arguments are wrong!\n$HELP_g\n\n";
$MM_g      =~ m/^\d+$/   ||  die   "The Command Line Arguments are wrong!\n$HELP_g\n\n";
$genome_g  =~ m/^\S+$/   ||  die   "The Command Line Arguments are wrong!\n$HELP_g\n\n";


######### say Command Arguments to Standard Output ###########
say  "\n
        ################ Your Arguments ###############################
                Input  folder:  $input_g
                Output folder:  $output_g
                Maximum number of mis-matched bases: $MM_g
                Reference genome: $genome_g
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
&myMakeDir($output_g);
&myMakeDir($output2_g);
opendir(my $DH_input, $input_g)  ||  die;     
my @inputFiles = readdir($DH_input);
my $pattern = "[-.0-9A-Za-z]+";
my $numCores = 4;





## Context specific:
my  $commonPath1        = "/home/yp/.MyProgramFiles/3_HTS-2G/4_ChIPseq/2_Mapping";
my  $commonPath2        = "/home/yp/.MyProgramFiles/3_HTS-2G/7_RNAseq/2_Mapping";
my  $Subjunc_index      = "$commonPath1/subread-1.5.0-p1/RefGenomes/$genome_g/$genome_g";
my  $GSNAP_index        = "$commonPath1/GSNAP/RefGenomes/$genome_g/$genome_g";
my  $Kallisto_index     = "$commonPath2/Kallisto/RefGenomes/$genome_g/$genome_g";  
my  $HISAT2_index       = "$commonPath2/hisat2/RefGenomes/$genome_g/$genome_g";
my  $TopHat2_index      = "$commonPath2/TopHat2/RefGenomes/$genome_g/$genome_g";
my  $MapSplice2_index   = "$commonPath2/MapSplice2/RefGenomes/$genome_g/$genome_g";
my  $STAR_index         = "$commonPath2/STAR/RefGenomes/$genome_g/$genome_g";
my  $BBMap_index        = "$commonPath2/BBMAP/RefGenomes/$genome_g/$genome_g";
my  $mrsFAST_index      = "$commonPath2/mrsFAST/RefGenomes/$genome_g/$genome_g/$genome_g";
my  $RSEM_index         = "$commonPath2/RSEM/RefGenomes/$genome_g";
my  $CRAC_index         = "$commonPath2/CRAC/RefGenomes/$genome_g/$genome_g"; 
my  $ContexMap2_index   = "$commonPath2/ContexMap2/RefGenomes/$genome_g/$genome_g";
my  $Salmon_index       = "$commonPath2/Salmon/RefGenomes/$genome_g/$genome_g";
my  $Sailfish_index     = "$commonPath2/Sailfish/RefGenomes/Shortcuts/$genome_g.fa";
my  $RapMap_index       = "$commonPath2/RapMap/RefGenomes/$genome_g";





say   "\n\n\n\n\n\n##################################################################################################";
say   "Checking all the necessary softwares in this step......" ;
sub printVersion  {
    my $software = $_[0];
    system("echo    '##############################################################################'  >> $output2_g/VersionsOfSoftwares.txt   2>&1");
    system("echo    '#########$software'                                                              >> $output2_g/VersionsOfSoftwares.txt   2>&1");
    system("$software                                                                                 >> $output2_g/VersionsOfSoftwares.txt   2>&1");
    system("echo    '\n\n\n\n\n\n'                                                                    >> $output2_g/VersionsOfSoftwares.txt   2>&1");
}
my  $Picard   = "/home/yp/.MyProgramFiles/3_HTS-2G/2-BAMtools/picard-tools-2.1.1/picard.jar";
my  $QoRTs    = "/home/yp/.MyProgramFiles/3_HTS-2G/7_RNAseq/1_beforeMap/QoRTs_1.0.7/QoRTs.jar";
my  $RNASeQC  = "/home/yp/.MyProgramFiles/3_HTS-2G/7_RNAseq/1_beforeMap/RNA-SeQC_v1.1.8.jar";  
&printVersion("subjunc  -v");
&printVersion("gsnap  --version");
&printVersion("kallisto");
&printVersion("hisat2  --version");
&printVersion("tophat  --version");
##&printVersion("mapsplice.py    --version");    
&printVersion("STAR  --version");
#&printVersion("rsem");
&printVersion("crac    -v");
#&printVersion("ContextMap_v2.7.1.jar");
&printVersion("salmon -v");
&printVersion("sailfish   -v");
&printVersion("rapmap");
&printVersion("mrsfast  -v");
&printVersion("bbmap.sh");
&printVersion("novoalign  --version");
&printVersion("samtools");
&printVersion("fastqc   -v");
&printVersion("java  -jar   $QoRTs");
&printVersion("java  -jar   $RNASeQC");
&printVersion("samstat   -v");
&printVersion("bamqc  -v");
&printVersion("preseq");
&printVersion("NGSQC  -h");
&printVersion("qualimap  -v");
&printVersion("fastqp   -h");
&printVersion("bam_stat.py  -h");
&printVersion("RNA_fragment_size.py  -h");
&printVersion("read_duplication.py  -h");
&printVersion("multiqc   --version");
&printVersion("qc3.pl -h");
&printVersion("bamutils stats");
&printVersion("bammapdist  -h");   ##biobambam2
&printVersion("propmapped");
&printVersion("qualityScores");
&printVersion("commonkmers.sh");   ##BBMap
&printVersion("countgc.sh");       ##BBMap
&printVersion("java  -jar  $Picard   ");
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
        next unless $inputFiles[$i] =~ m/\.fastq$/;  
        next unless $inputFiles[$i] !~ m/^[.]/;
        next unless $inputFiles[$i] !~ m/[~]$/;
        next unless $inputFiles[$i] !~ m/^QC_Results$/;
        next unless $inputFiles[$i] !~ m/^unpaired/;
        say   "\t......$inputFiles[$i]" ; 
        my $temp = $inputFiles[$i]; 
        $groupFiles[++$#groupFiles] = $inputFiles[$i];  
        $temp =~ m/^(\d+)_($pattern)_(Rep[1-9])/   or  die   "wrong-1: ## $temp ##";
        $temp =~ m/_(Rep[1-9])\.fastq$/  or  $temp =~ m/_(Rep[1-9])_?([1-2]?)\.fastq$/   or  die   "wrong-2: ## $temp ##";
        if($temp !~ m/^((\d+)_($pattern)_(Rep[1-9]))(_[1-2])?\.fastq$/) {
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
say   "Detecting single-end and paired-end FASTQ files in input folder ......";
my @singleEnd = ();
my @pairedEnd = ();
open(seqFiles_FH, ">", "$output2_g/singleEnd-pairedEnd-Files.txt")  or  die; 
for ( my $i=0; $i<=$#inputFiles; $i++ ) {     
    next unless $inputFiles[$i] =~ m/\.fastq$/;
    next unless $inputFiles[$i] !~ m/^[.]/;
    next unless $inputFiles[$i] !~ m/[~]$/;
    next unless $inputFiles[$i] !~ m/^unpaired/;
    say    "\t......$inputFiles[$i]"; 
    if ($inputFiles[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])\.fastq$/) {   ## sinlge end sequencing files.
        $inputFiles[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])\.fastq$/  or  die;  
        $singleEnd[$#singleEnd+1] =  $inputFiles[$i];
        say   "\t\t\t\tSingle-end sequencing files:  $inputFiles[$i]\n";
        say   seqFiles_FH  "Single-end sequencing files: $inputFiles[$i]\n";
    }else{     ## paired end sequencing files.
        $inputFiles[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])_([1-2])\.fastq$/  or  die; 
        if ($inputFiles[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))_1\.fastq$/) { ## The two files of one paired sequencing sample are always side by side. 
            my $temp = $1;
            my $end1 = $temp."_1.fastq";
            my $end2 = $temp."_2.fastq";
            (-e  "$input_g/$end1")  or die;  
            (-e  "$input_g/$end2")  or die;
            $pairedEnd[$#pairedEnd+1] =  $end1;
            $pairedEnd[$#pairedEnd+1] =  $end2;
            say   "\t\t\t\tPaired-end sequencing files: $end1,  $end2\n";
            say   seqFiles_FH  "Paired-end sequencing files: $end1,  $end2\n";
        }
    }
}
( ($#pairedEnd+1)%2 == 0 )  or die;
say   seqFiles_FH  "\n\n\n\n\n";
say   seqFiles_FH  "All single-end sequencing files:@singleEnd\n\n\n";
say   seqFiles_FH  "All paired-end sequencing files:@pairedEnd\n\n\n";
say    "\t\t\t\tAll single-end sequencing files:@singleEnd\n\n";
say    "\t\t\t\tAll paired-end sequencing files:@pairedEnd\n\n";
my $numSingle = $#singleEnd + 1;
my $numPaired = $#pairedEnd + 1;
say seqFiles_FH   "\nThere are $numSingle single-end sequencing files.\n";
say seqFiles_FH   "\nThere are $numPaired paired-end sequencing files.\n";
say     "\t\t\t\tThere are $numSingle single-end sequencing files.\n";
say     "\t\t\t\tThere are $numPaired paired-end sequencing files.\n";





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
       my $QoRTsDir     = "$relativePath/9_QoRTs";
       my $RNASeQCDir   = "$relativePath/10_RNA-SeQC";
       my $RSeQCDir     = "$relativePath/11_RSeQC";
       &myMakeDir("$relativePath");
       &myMakeDir("$SAMtools");
       &myMakeDir("$FastQC");
       &myMakeDir("$samstat");
       &myMakeDir("$PRESEQ");
       &myMakeDir("$qualimap");
       &myMakeDir("$MultiQC");
       &myMakeDir("$SubreadUti");
       &myMakeDir("$PicardDir");
       &myMakeDir("$QoRTsDir");
       &myMakeDir("$RNASeQCDir");
       &myMakeDir("$RSeQCDir");
       opendir(my $DH_map, $relativePath1) || die;     
       my @mapFiles = readdir($DH_map);

       say   "\n\n\n\n\n\nDetecting the quality of bam files by using SAMtools, FASTQC, samstat, PRESEQ, qualimap, QoRTs, RNA-SeQC, RSeQC, Subread utilities, Picard and MultiQC ......";
       for (my $i=0; $i<=$#mapFiles; $i++) {
           next unless $mapFiles[$i] =~ m/\.sam$/;
           next unless $mapFiles[$i] !~ m/^[.]/;
           next unless $mapFiles[$i] !~ m/[~]$/;
           my $temp = $mapFiles[$i]; 
           $temp =~ s/\.sam$//  ||  die; 
           say   "\t......$mapFiles[$i]";  
           system("samtools  sort  -m 2G  -o $relativePath1/$temp.t.bam   --output-fmt bam  -T $relativePath1/1_$temp   --threads $numCores    $relativePath1/$temp.sam    >>$SAMtools/$temp.runLog    2>&1");
           system("bamutils  cleancigar      $relativePath1/$temp.t.bam    $relativePath1/$temp.bam   >> $relativePath/$temp.cleancigar.runLog     2>&1 ");
           system("samtools  index           $relativePath1/$temp.bam      >>$SAMtools/$temp.runLog  2>&1");
           system("samtools  flagstat        $relativePath1/$temp.bam      >>$SAMtools/$temp.runLog  2>&1");
           system(`samtools  idxstats        $relativePath1/$temp.bam      >>$SAMtools/$temp.runLog  2>&1`);
           system("rm   $relativePath1/$temp.sam"); 
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
           #system("java  -jar   $QoRTs  QC  $relativePath1/$temp.bam     >> $SubreadUti/$temp.prommapped      2>&1"); 
           #system("propmapped   -i $relativePath1/$temp.bam       -f   -p      -o $SubreadUti/$temp.prommapped      >> $SubreadUti/$temp.prommapped      2>&1"); 
           #system("propmapped   -i $relativePath1/$temp.bam       -f   -p      -o $SubreadUti/$temp.prommapped      >> $SubreadUti/$temp.prommapped      2>&1"); 
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
           next unless $mapFiles[$i] =~ m/\.sam$/;
           next unless $mapFiles[$i] !~ m/^[.]/;
           next unless $mapFiles[$i] !~ m/[~]$/;
           my $temp = $mapFiles[$i]; 
           $temp =~ s/\.sam$//  ||  die; 
           say   "\t......$mapFiles[$i]";  
           system( "bamtools stats  -in $relativePath1/$temp.bam    -insert  >> $bamtools/$temp.runLog   2>&1");
           system( "bam      stats  -in $relativePath1/$temp.bam    --basic   --pBaseQC $BamUtil/$temp.pBaseQC  >> $BamUtil/$temp.runLog   2>&1");
           system( "fastqp    --nreads 20000000   --kmer 5    --output $Fastqp/$temp  --type fastq   --median-qual 30     $relativePath1/$temp.bam     >> $Fastqp/$temp.runLog     2>&1 " );
           system( "bamutils stats    --nreads 20000000   -all     $relativePath1/$temp.bam     >> $NGSutils/$temp.runLog     2>&1 " );
       }
}  ## End myQC2

 






if(1==1) { ########## Start HISAT2
say   "\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using HISAT2 ......";
my $HISAT2     = "$output_g/1-HISAT2";  
my $HISAT2_2   = "$output_g/1-HISAT2/QC_Results";  
&myMakeDir("$HISAT2");
&myMakeDir("$HISAT2_2");

for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {  
        say    "\t......$pairedEnd[$i]";
        say    "\t......$pairedEnd[$i+1]";
        $pairedEnd[$i]   =~ m/^((\d+)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d+)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$HISAT2_2/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("hisat2         --threads $numCores   -q   --phred33   --end-to-end    -x $HISAT2_index    -1 $input_g/$end1.fastq        -2 $input_g/$end2.fastq     -S $HISAT2/$temp.sam    >>$HISAT2_2/$temp.runLog  2>&1");                                                                                                                                 
}

for (my $i=0; $i<=$#singleEnd; $i++) {   
        say   "\t......$singleEnd[$i]";
        $singleEnd[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("hisat2        --threads $numCores   -q   --phred33   --end-to-end    -x $HISAT2_index    -U $input_g/$temp.fastq                                    -S $HISAT2/$temp.sam    >>$HISAT2_2/$temp.runLog  2>&1");                                                                
}

&myQC1($HISAT2);
}  ########## End HISAT2










if(1==1) { ########## Start Subjunc
say   "\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using Subjunc ......";
my $subjunc   = "$output_g/2-Subjunc";  
my $subjunc2  = "$output_g/2-Subjunc/QC_Results";  
&myMakeDir("$subjunc");
&myMakeDir("$subjunc2");

for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {  
        say   "\t......$pairedEnd[$i]";
        say   "\t......$pairedEnd[$i+1]";
        $pairedEnd[$i]   =~ m/^((\d+)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d+)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$subjunc2/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("subjunc  -T $numCores  -I 20  -B 1  -M $MM_g    --SAMoutput  -d 0  -D 1000   -i $Subjunc_index   -r $input_g/$end1.fastq   -R  $input_g/$end2.fastq   -o  $subjunc/$temp.sam   >>$subjunc2/$temp.runLog  2>&1");               
}

for (my $i=0; $i<=$#singleEnd; $i++) {  
        say   "\t......$singleEnd[$i]"; 
        $singleEnd[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.fastq$/   or  die; 
        my $temp = $1; 
        system("subjunc  -T $numCores  -I 20  -B 1  -M $MM_g   --SAMoutput   -i $Subjunc_index    -r $input_g/$temp.fastq    -o $subjunc/$temp.sam    >>$subjunc2/$temp.runLog   2>&1");       
}

&myQC1($subjunc);
} ########## End subjunc
  








 



if(1==2) { ########## Start GSNAP
say   "\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using GSNAP ......";
my $GSNAP   = "$output_g/4-GSNAP";  
my $GSNAP_2 = "$output_g/4-GSNAP/Results";
if ( !(-e $GSNAP)   )  { mkdir $GSNAP   || die; }
if ( !(-e $GSNAP_2) )  { mkdir $GSNAP_2 || die; }

for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {  
        say    "\t......$pairedEnd[$i]";
        say    "\t......$pairedEnd[$i+1]";
        $pairedEnd[$i]   =~ m/^((\d+)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d+)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$GSNAP_2/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("gsnap    --db=$GSNAP_index  -o $GSNAP/$temp.sam  --nthreads=8    -M $input_g/$end1.fastq  $input_g/$end2.fastq    >>$GSNAP_2/$temp.runLog  2>&1");                                                                                                                                 
        system("samtools  sort  -m 2G  -o $GSNAP/$temp.bam   -O bam  -T $GSNAP/1_$temp   -@ 8    $GSNAP/$temp.sam    >>$GSNAP_2/samtools_$temp.runLog    2>&1");
        system("samtools  index           $GSNAP/$temp.bam      >>$GSNAP_2/samtools_$temp.runLog  2>&1");
        system("samtools  flagstat        $GSNAP/$temp.bam      >>$GSNAP_2/samtools_$temp.runLog  2>&1");
        system(`samtools  idxstats        $GSNAP/$temp.bam      >>$GSNAP_2/samtools_$temp.runLog  2>&1`);
        system("rm   $GSNAP/$temp.sam");
}

for (my $i=0; $i<=$#singleEnd; $i++) {   
        say   "\t......$singleEnd[$i]";
        $singleEnd[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("gsnap    --db=$GSNAP_index  -o $GSNAP/$temp.sam  --nthreads=8    -M $input_g/$temp.fastq    >>$GSNAP_2/$temp.runLog");                                                                
        system("samtools  sort  -m 2G  -o $GSNAP/$temp.bam   -O bam  -T $GSNAP/1_$temp   -@ 8    $GSNAP/$temp.sam    >>$GSNAP_2/samtools_$temp.runLog    2>&1");
        system("samtools  index           $GSNAP/$temp.bam      >>$GSNAP_2/samtools_$temp.runLog  2>&1");
        system("samtools  flagstat        $GSNAP/$temp.bam      >>$GSNAP_2/samtools_$temp.runLog  2>&1");
        system(`samtools  idxstats        $GSNAP/$temp.bam      >>$GSNAP_2/samtools_$temp.runLog  2>&1`);
        system("rm   $GSNAP/$temp.sam");   
}

my $FastQC  = "$GSNAP_2/FastQC";
my $QCstat  = "$GSNAP_2/QCstat";
if ( !( -e $FastQC)  )   { mkdir $FastQC    ||  die; }
if ( !( -e $QCstat)  )   { mkdir $QCstat    ||  die; }
opendir(my $DH_map, $GSNAP) || die;     
my @mapFiles = readdir($DH_map);

say   "\n\nDetecting the quality of bam files ......";
for (my $i=0; $i<=$#mapFiles; $i++) {
    next unless $mapFiles[$i] =~ m/\.bam$/;
    next unless $mapFiles[$i] !~ m/^[.]/;
    next unless $mapFiles[$i] !~ m/[~]$/;
    my $temp = $mapFiles[$i]; 
    $temp =~ s/\.bam$//  ||  die; 
    say   "\t......$mapFiles[$i]";   
    system("fastqc    --outdir $FastQC     --threads 8    --format bam    --kmers 7     $GSNAP/$temp.bam   >>$FastQC/$temp.runLog        2>&1"); 
    system("samstat   $GSNAP/$temp.bam      >> $GSNAP_2/$temp.runLog         2>&1");
    system("qualimap  bamqc  -bam $GSNAP/$temp.bam   -c  -nt 8  -outdir $QCstat/qualimap_$temp   --java-mem-size=8G   >>$GSNAP_2/$temp.runLog    2>&1");
    system("propmapped   -i $GSNAP/$temp.bam                    -o $QCstat/$temp.prommapped      >> $QCstat/$temp.prommapped      2>&1");
    system("echo      '\n\n\n\n\n'                                                               >> $QCstat/$temp.prommapped      2>&1"); 
    system("propmapped   -i $GSNAP/$temp.bam       -f           -o $QCstat/$temp.prommapped      >> $QCstat/$temp.prommapped      2>&1"); 
    system("echo      '\n\n\n\n\n'                                                               >> $QCstat/$temp.prommapped      2>&1"); 
    system("propmapped   -i $GSNAP/$temp.bam       -f   -p      -o $QCstat/$temp.prommapped      >> $QCstat/$temp.prommapped      2>&1"); 
    system("qualityScores   --SAMinput   -i $GSNAP/$temp.bam    -o $QCstat/$temp.qualityScores   >> $QCstat/$temp.qualityScores   2>&1"); 
}

}  ########## End GSNAP





if(1==2) { ########## Start STAR
say   "\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using STAR ......";
my $STAR   = "$output_g/5-STAR";  
my $STAR_2 = "$output_g/5-STAR/Results";
if ( !(-e $STAR)   )  { mkdir $STAR   || die; }
if ( !(-e $STAR_2) )  { mkdir $STAR_2 || die; }
my $HOME = '$HOME';  

for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {  
        say    "\n\t......$pairedEnd[$i]";
        say    "\t......$pairedEnd[$i+1]\n";
        $pairedEnd[$i]   =~ m/^((\d+)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d+)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$STAR_2/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $STAR\n";
        system("STAR  --runMode alignReads  --runThreadN 8    --sjdbGTFfile 0-Other/$genome_g-Ensembl  --outFilterType Normal  --outFilterMultimapNmax 20  --alignSJoverhangMin 5  --alignSJDBoverhangMin 1   --outFilterMismatchNmax $MM_g   --alignIntronMin 20   --alignIntronMax 1000000    --alignMatesGapMax 1000000   --outFileNamePrefix  $HOME/$temp   --genomeDir $STAR_index  --readFilesIn $input_g/$end1.fastq $input_g/$end2.fastq   >>$STAR_2/$temp.runLog  2>&1 ");    
        system("cp  -rf   $HOME/$temp*   $STAR");
        system("rename  s/Aligned.out.sam/.sam/   $STAR/*Aligned.out.sam"); 
        system("samtools  sort  -m 2G  -o $STAR/$temp.bam   -O bam  -T $STAR/1_$temp   -@ 8    $STAR/$temp.sam    >>$STAR_2/samtools_$temp.runLog    2>&1");
        system("samtools  index           $STAR/$temp.bam      >>$STAR_2/samtools_$temp.runLog  2>&1");
        system("samtools  flagstat        $STAR/$temp.bam      >>$STAR_2/samtools_$temp.runLog  2>&1");
        system(`samtools  idxstats        $STAR/$temp.bam      >>$STAR_2/samtools_$temp.runLog  2>&1`);  
        system("rm  -rf  $HOME/$temp*"); 
        system("rm   $STAR/$temp.sam");                                                                                                                              
}

for (my $i=0; $i<=$#singleEnd; $i++) {   
        say   "\n\t......$singleEnd[$i]\n";
        $singleEnd[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("STAR  --runMode alignReads  --runThreadN 8    --sjdbGTFfile 0-Other/$genome_g-Ensembl  --outFilterType Normal  --outFilterMultimapNmax 20  --alignSJoverhangMin 5  --alignSJDBoverhangMin 1   --outFilterMismatchNmax $MM_g   --alignIntronMin 20   --alignIntronMax 1000000    --alignMatesGapMax 1000000   --outFileNamePrefix  $HOME/$temp   --genomeDir $STAR_index  --readFilesIn $input_g/$temp.fastq   >>$STAR_2/$temp.runLog  2>&1 ");    
        system("cp  -rf   $HOME/$temp*   $STAR");
        system("rename  s/Aligned.out.sam/.sam/   $STAR/*Aligned.out.sam"); 
        system("samtools  sort  -m 2G  -o $STAR/$temp.bam   -O bam  -T $STAR/1_$temp   -@ 8    $STAR/$temp.sam    >>$STAR_2/samtools_$temp.runLog    2>&1");
        system("samtools  index           $STAR/$temp.bam      >>$STAR_2/samtools_$temp.runLog  2>&1");
        system("samtools  flagstat        $STAR/$temp.bam      >>$STAR_2/samtools_$temp.runLog  2>&1");
        system(`samtools  idxstats        $STAR/$temp.bam      >>$STAR_2/samtools_$temp.runLog  2>&1`); 
        system("rm  -rf  $HOME/$temp*");    
        system("rm   $STAR/$temp.sam");                                                            
}

my $FastQC  = "$STAR_2/FastQC";
my $QCstat  = "$STAR_2/QCstat";
if ( !( -e $FastQC)  )   { mkdir $FastQC    ||  die; }
if ( !( -e $QCstat)  )   { mkdir $QCstat    ||  die; }
opendir(my $DH_map, $STAR) || die;     
my @mapFiles = readdir($DH_map);

say   "\n\nDetecting the quality of bam files ......";
for (my $i=0; $i<=$#mapFiles; $i++) {
    next unless $mapFiles[$i] =~ m/\.bam$/;
    next unless $mapFiles[$i] !~ m/^[.]/;
    next unless $mapFiles[$i] !~ m/[~]$/;
    my $temp = $mapFiles[$i]; 
    $temp =~ s/\.bam$//  ||  die; 
    say   "\t......$mapFiles[$i]";   
    system("fastqc    --outdir $FastQC     --threads 8    --format bam    --kmers 7     $STAR/$temp.bam   >>$FastQC/$temp.runLog        2>&1"); 
    system("samstat   $STAR/$temp.bam      >> $STAR_2/$temp.runLog         2>&1");
    system("qualimap  bamqc  -bam $STAR/$temp.bam   -c  -nt 8  -outdir $QCstat/qualimap_$temp   --java-mem-size=8G   >>$STAR_2/$temp.runLog    2>&1");
    system("propmapped   -i $STAR/$temp.bam                    -o $QCstat/$temp.prommapped      >> $QCstat/$temp.prommapped      2>&1");
    system("echo      '\n\n\n\n\n'                                                              >> $QCstat/$temp.prommapped      2>&1"); 
    system("propmapped   -i $STAR/$temp.bam       -f           -o $QCstat/$temp.prommapped      >> $QCstat/$temp.prommapped      2>&1"); 
    system("echo      '\n\n\n\n\n'                                                              >> $QCstat/$temp.prommapped      2>&1"); 
    system("propmapped   -i $STAR/$temp.bam       -f   -p      -o $QCstat/$temp.prommapped      >> $QCstat/$temp.prommapped      2>&1"); 
    system("qualityScores   --SAMinput   -i $STAR/$temp.bam    -o $QCstat/$temp.qualityScores   >> $QCstat/$temp.qualityScores   2>&1"); 
}

}  ########## End STAR





if(1==2) { ########## Start CRAC
say   "\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using CRAC ......";
my $CRAC   = "$output_g/6-CRAC";  
my $CRAC_2 = "$output_g/6-CRAC/Results";
if ( !(-e $CRAC)   )  { mkdir $CRAC   || die; }
if ( !(-e $CRAC_2) )  { mkdir $CRAC_2 || die; }

for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {  
        say    "\t......$pairedEnd[$i]";
        say    "\t......$pairedEnd[$i+1]";
        $pairedEnd[$i]   =~ m/^((\d+)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d+)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$CRAC_2/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("crac  -i $CRAC_index    -r $input_g/$end1.fastq  $input_g/$end2.fastq   -k 22  -o $CRAC/$temp.sam  --nb-threads 8  --summary $CRAC/$temp.summary   >>$CRAC_2/$temp.runLog  2>&1");                                                                                                                                 
        system("samtools  sort  -m 2G  -o $CRAC/$temp.bam   -O bam  -T $CRAC/1_$temp   -@ 8    $CRAC/$temp.sam    >>$CRAC_2/samtools_$temp.runLog    2>&1");
        system("samtools  index           $CRAC/$temp.bam      >>$CRAC_2/samtools_$temp.runLog  2>&1");
        system("samtools  flagstat        $CRAC/$temp.bam      >>$CRAC_2/samtools_$temp.runLog  2>&1");
        system(`samtools  idxstats        $CRAC/$temp.bam      >>$CRAC_2/samtools_$temp.runLog  2>&1`);
        system("rm   $CRAC/$temp.sam");                                                            
}

for (my $i=0; $i<=$#singleEnd; $i++) {   
        say   "\t......$singleEnd[$i]";
        $singleEnd[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("crac  -i $CRAC_index    -r $input_g/$temp.fastq   -k 22  -o $CRAC/$temp.sam  --nb-threads 8  --summary $CRAC/$temp.summary   >>$CRAC_2/$temp.runLog  2>&1");                                                                
        system("samtools  sort  -m 2G  -o $CRAC/$temp.bam   -O bam  -T $CRAC/1_$temp   -@ 8    $CRAC/$temp.sam    >>$CRAC_2/samtools_$temp.runLog    2>&1");
        system("samtools  index           $CRAC/$temp.bam      >>$CRAC_2/samtools_$temp.runLog  2>&1");
        system("samtools  flagstat        $CRAC/$temp.bam      >>$CRAC_2/samtools_$temp.runLog  2>&1");
        system(`samtools  idxstats        $CRAC/$temp.bam      >>$CRAC_2/samtools_$temp.runLog  2>&1`);
        system("rm   $CRAC/$temp.sam");                                                            
}

my $FastQC  = "$CRAC_2/FastQC";
my $QCstat  = "$CRAC_2/QCstat";
if ( !( -e $FastQC)  )   { mkdir $FastQC    ||  die; }
if ( !( -e $QCstat)  )   { mkdir $QCstat    ||  die; }
opendir(my $DH_map, $CRAC) || die;     
my @mapFiles = readdir($DH_map);

say   "\n\nDetecting the quality of bam files ......";
for (my $i=0; $i<=$#mapFiles; $i++) {
    next unless $mapFiles[$i] =~ m/\.bam$/;
    next unless $mapFiles[$i] !~ m/^[.]/;
    next unless $mapFiles[$i] !~ m/[~]$/;
    my $temp = $mapFiles[$i]; 
    $temp =~ s/\.bam$//  ||  die; 
    say   "\t......$mapFiles[$i]";   
    system("fastqc    --outdir $FastQC     --threads 8    --format bam    --kmers 7     $CRAC/$temp.bam   >>$FastQC/$temp.runLog        2>&1"); 
    system("samstat   $CRAC/$temp.bam      >> $CRAC_2/$temp.runLog         2>&1");      
    system("qualimap  bamqc  -bam $CRAC/$temp.bam   -c  -nt 8  -outdir $QCstat/qualimap_$temp   --java-mem-size=8G   >>$CRAC_2/$temp.runLog    2>&1");    
    system("propmapped   -i $CRAC/$temp.bam                    -o $QCstat/$temp.prommapped      >> $QCstat/$temp.prommapped      2>&1");
    system("echo      '\n\n\n\n\n'                                                              >> $QCstat/$temp.prommapped      2>&1"); 
    system("propmapped   -i $CRAC/$temp.bam       -f           -o $QCstat/$temp.prommapped      >> $QCstat/$temp.prommapped      2>&1"); 
    system("echo      '\n\n\n\n\n'                                                              >> $QCstat/$temp.prommapped      2>&1"); 
    system("propmapped   -i $CRAC/$temp.bam       -f   -p      -o $QCstat/$temp.prommapped      >> $QCstat/$temp.prommapped      2>&1"); 
    system("qualityScores   --SAMinput   -i $CRAC/$temp.bam    -o $QCstat/$temp.qualityScores   >> $QCstat/$temp.qualityScores   2>&1"); 
}

}  ########## End CRAC





if(1==2) { ########## Start BBMap
say   "\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using BBMap ......";
my $BBMap  = "$output_g/7-BBMap";  
my $BBMap2 = "$output_g/7-BBMap/Results";
if ( !(-e $BBMap)  )  { mkdir $BBMap   || die; }
if ( !(-e $BBMap2) )  { mkdir $BBMap2  || die; }

for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {  
        say    "\t......$pairedEnd[$i]";
        say    "\t......$pairedEnd[$i+1]";
        $pairedEnd[$i]   =~ m/^((\d+)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d+)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq"  eq  $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$BBMap2/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("bbmap.sh  path=$BBMap_index  in=$input_g/$end1.fastq  in2=$input_g/$end2.fastq   out=$BBMap/$temp.sam      pairlen=32000  overwrite=f  threads=8  maxindel=16000  >>$BBMap2/$temp.runLog  2>&1");                                                                                                                                 
        system("samtools  sort  -m 2G  -o $BBMap/$temp.bam   -O bam  -T $BBMap/1_$temp   -@ 8    $BBMap/$temp.sam    >>$BBMap2/samtools_$temp.runLog    2>&1");
        system("samtools  index           $BBMap/$temp.bam      >>$BBMap2/samtools_$temp.runLog  2>&1");
        system("samtools  flagstat        $BBMap/$temp.bam      >>$BBMap2/samtools_$temp.runLog  2>&1");
        system(`samtools  idxstats        $BBMap/$temp.bam      >>$BBMap2/samtools_$temp.runLog  2>&1`);
        system("rm   $BBMap/$temp.sam");                                                            
}

for (my $i=0; $i<=$#singleEnd; $i++) {   
        say   "\t......$singleEnd[$i]";
        $singleEnd[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("bbmap.sh  path=$BBMap_index  in=$input_g/$temp.fastq    out=$BBMap/$temp.sam   overwrite=f  threads=8     maxindel=16000   >>$BBMap2/$temp.runLog  2>&1");                                                                
        system("samtools  sort  -m 2G  -o $BBMap/$temp.bam   -O bam  -T $BBMap/1_$temp   -@ 8    $BBMap/$temp.sam    >>$BBMap2/samtools_$temp.runLog    2>&1");
        system("samtools  index           $BBMap/$temp.bam      >>$BBMap2/samtools_$temp.runLog  2>&1");
        system("samtools  flagstat        $BBMap/$temp.bam      >>$BBMap2/samtools_$temp.runLog  2>&1");
        system(`samtools  idxstats        $BBMap/$temp.bam      >>$BBMap2/samtools_$temp.runLog  2>&1`);
        system("rm   $BBMap/$temp.sam");                                                            
}

my $FastQC  = "$BBMap2/FastQC";
my $QCstat  = "$BBMap2/QCstat";
if ( !( -e $FastQC)  )   { mkdir $FastQC    ||  die; }
if ( !( -e $QCstat)  )   { mkdir $QCstat    ||  die; }
opendir(my $DH_map, $BBMap) || die;     
my @mapFiles = readdir($DH_map);

say   "\n\nDetecting the quality of bam files ......";
for (my $i=0; $i<=$#mapFiles; $i++) {
    next unless $mapFiles[$i] =~ m/\.bam$/;
    next unless $mapFiles[$i] !~ m/^[.]/;
    next unless $mapFiles[$i] !~ m/[~]$/;
    my $temp = $mapFiles[$i]; 
    $temp =~ s/\.bam$//  ||  die; 
    say   "\t......$mapFiles[$i]";   
    system("fastqc    --outdir $FastQC     --threads 8    --format bam    --kmers 7     $BBMap/$temp.bam   >>$FastQC/$temp.runLog        2>&1"); 
    system("samstat   $BBMap/$temp.bam      >> $BBMap2/$temp.runLog         2>&1");      
    system("qualimap  bamqc  -bam $BBMap/$temp.bam   -c  -nt 8  -outdir $QCstat/qualimap_$temp   --java-mem-size=8G   >>$BBMap2/$temp.runLog    2>&1");    
    system("propmapped   -i $BBMap/$temp.bam                    -o $QCstat/$temp.prommapped      >> $QCstat/$temp.prommapped      2>&1");
    system("echo      '\n\n\n\n\n'                                                               >> $QCstat/$temp.prommapped      2>&1"); 
    system("propmapped   -i $BBMap/$temp.bam       -f           -o $QCstat/$temp.prommapped      >> $QCstat/$temp.prommapped      2>&1"); 
    system("echo      '\n\n\n\n\n'                                                               >> $QCstat/$temp.prommapped      2>&1"); 
    system("propmapped   -i $BBMap/$temp.bam       -f   -p      -o $QCstat/$temp.prommapped      >> $QCstat/$temp.prommapped      2>&1"); 
    system("qualityScores   --SAMinput   -i $BBMap/$temp.bam    -o $QCstat/$temp.qualityScores   >> $QCstat/$temp.qualityScores   2>&1"); 
}

}  ########## End BBMap





say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n";





## END
