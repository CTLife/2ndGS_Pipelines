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

        Step 3: Mapping reads to the reference genome by using popular softwares. 
                15 aligners (mappers) are available: 
                            BWA,     Stampy    and GEM.
                            Bowtie,  CUSHAW3   and Pash3.
                            Subread, GSNAP     and SNAP.                            
                            mrsFAST, Novoalign and BBMap. 
                            Yara,    RazerS3   and SplazerS. 

                Assess the quality of ChIP-Seq reads (in BAM files) to identify possible sequencing errors or biases by using 20 softwares:
                	SAMtools, FastQC, bamtools stats, bam stats (BamUtil), samstat, BamQC (Simon Andrews), PRESEQ, NGSQC, qualimap, fastqp, htseq-qa,   
                	multiqc, BBMap, qc3.pl, bamutils stats in ngsutils, biobambam2, Subread utilities, 10 tools in Picard, QuasR and Rqc.                 

        Usage:  
               perl  CISDA-3.pl    [-v]    [-h]    [-i inputDir]    [-o outDir]   [-MM N]   [-g reference_genome ]   [-60 boole]  [-PE INT]
        For instance: 
                     perl  CISDA-3.pl    -i 3-Filtered          -o 4-Mapping         -MM 6    -g mm10   -60 yes  -PE 600
                     perl  CISDA-3.pl    --input 3-Filtered     --output 4-Mapping   --MaxMismatch 6  --genome mm10  --60bp yes   --MaxPE 600
                     perl  CISDA-3.pl    --input 3-Filtered     --output 4-Mapping   --MaxMismatch 6  --genome mm10  --60bp yes   --MaxPE 600   >> CISDA-3.runLog  2>&1
                           CISDA-3.pl    --input 3-Filtered     --output 4-Mapping   --MaxMismatch 6  --genome mm10  --60bp yes   --MaxPE 600   >> CISDA-3.runLog  2>&1
     
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

        -60 booleValue, --60bp booleValue     booleValue is "no"  if reads length <= 60 bp, then Bowtie1 and BWA-aln will be invoked. 
                                              booleValue is "yes" if reads length >  60 bp, then Bowtie2 and BWA-mem will be invoked.
                                              (no default)     

        -PE INT, --MaxPE INT     Maximum fragment/insert length for paired-end reads.   (no default)                                               
        ------------------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines such as RASDA, MESDA and HISDA,
        please see the web site: https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies 
        and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.    
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------  
';


########## Version Infromation ##########
my $version_g = "  The Third Step of CISDA (ChIP-Seq Data Analyzer), version 0.7.1, 2016-04-11.";


########## Keys and Values ##########
if ($#ARGV   == -1) { say  "\n$HELP_g\n";  exit 0; }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0) { @ARGV = (@ARGV, "-h");       }       ## when the number of command argumants is odd. 
my %args = @ARGV;


########## Initialize  Variables ##########
my $input_g  = '3-Filtered';      ## This is only an initialization  value or suggesting value, not default value.
my $output_g = '4-Mapping';       ## This is only an initialization  value or suggesting value, not default value.
my $MM_g     = 6;                 ## This is only an initialization  value or suggesting value, not default value.
my $genome_g = 'mm10';            ## This is only an initialization  value or suggesting value, not default value.
my $readLen_g= 'yes';             ## This is only an initialization  value or suggesting value, not default value.
my $maxPE_g  = '600';             ## This is only an initialization  value or suggesting value, not default value.


########## Available Arguments ##########
my $available = "  -v  --version    -h  --help    -i  --input    -o    --output  -MM  --MaxMismatch    -g  --genome  -60  --60bp   -PE  --MaxPE   ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if($available !~ m/\s$key\s/) {say    "\n\tCann't recognize $key !!";  $boole_g = 1; }
}
if($boole_g == 1) {
    say   "\tThe Command Line Arguments are wrong!";
    say   "\tPlease see help message by using 'perl  CISDA-3.pl  -h' \n";
    exit 0;
}


########## Get Arguments ##########
if ( ( exists $args{'-v' } )  or  ( exists $args{'--version'      } )  )     { say  "\n$version_g\n";    exit 0; }
if ( ( exists $args{'-h' } )  or  ( exists $args{'--help'         } )  )     { say  "\n$HELP_g\n";       exit 0; }
if ( ( exists $args{'-i' } )  or  ( exists $args{'--input'        } )  )     { ($input_g   = $args{'-i'  })  or  ($input_g   = $args{'--input'      });  }else{say   "\n -i  or --input  is required.\n";          say  "\n$HELP_g\n";       exit 0; }                                               
if ( ( exists $args{'-o' } )  or  ( exists $args{'--output'       } )  )     { ($output_g  = $args{'-o'  })  or  ($output_g  = $args{'--output'     });  }else{say   "\n -o  or --output is required.\n";          say  "\n$HELP_g\n";       exit 0; }      
if ( ( exists $args{'-MM'} )  or  ( exists $args{'--MaxMismatch'  } )  )     { ($MM_g      = $args{'-MM' })  or  ($MM_g      = $args{'--MaxMismatch'});  }else{say   "\n -MM or --MaxMismatch is required.\n";     say  "\n$HELP_g\n";       exit 0; } 
if ( ( exists $args{'-g' } )  or  ( exists $args{'--genome'       } )  )     { ($genome_g  = $args{'-g'  })  or  ($genome_g  = $args{'--genome'     });  }else{say   "\n -g  or --genome is required.\n";          say  "\n$HELP_g\n";       exit 0; } 
if ( ( exists $args{'-60'} )  or  ( exists $args{'--60bp'         } )  )     { ($readLen_g = $args{'-60' })  or  ($readLen_g = $args{'--60bp'       });  }else{say   "\n -60 or --60bp is required.\n";            say  "\n$HELP_g\n";       exit 0; } 
if ( ( exists $args{'-PE'} )  or  ( exists $args{'--MaxPE'        } )  )     { ($maxPE_g   = $args{'-PE' })  or  ($maxPE_g   = $args{'--MaxPE'      });  }else{say   "\n -PE or --MaxPE is required.\n";           say  "\n$HELP_g\n";       exit 0; } 


########### Conditions #############
$input_g   =~ m/^\S+$/   ||  die   "The Command Line Arguments are wrong!\n$HELP_g\n\n";
$output_g  =~ m/^\S+$/   ||  die   "The Command Line Arguments are wrong!\n$HELP_g\n\n";
$MM_g      =~ m/^\d+$/   ||  die   "The Command Line Arguments are wrong!\n$HELP_g\n\n";
$genome_g  =~ m/^\S+$/   ||  die   "The Command Line Arguments are wrong!\n$HELP_g\n\n";
$readLen_g =~ m/^\S+$/   ||  die   "The Command Line Arguments are wrong!\n$HELP_g\n\n";
$maxPE_g   =~ m/^\d+$/   ||  die   "The Command Line Arguments are wrong!\n$HELP_g\n\n";
($readLen_g  eq  "no")   or  ($readLen_g  eq  "yes")   or  die;


######### say Command Arguments to Standard Output ###########
say  "\n
        ################ Your Arguments ###############################
                Input  folder:  $input_g
                Output folder:  $output_g
                Maximum number of mis-matched bases: $MM_g
                Reference genome: $genome_g
                > 60bp: $readLen_g
                Maximum fragment length for paired-end reads: $maxPE_g 
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
my  $commonPath         = "/home/yp/.MyProgramFiles/3_HTS-2G/4_ChIPseq/2_Mapping";
my  $BWA_index          = "$commonPath/bwa-0.7.13/RefGenomes/$genome_g/$genome_g";
my  $Stampy_index       = "$commonPath/stampy-1.0.28/RefGenomes/$genome_g/$genome_g";
my  $GEM_index          = "$commonPath/gemtools-1.7.1/RefGenomes/$genome_g/$genome_g";  
my  $Bowtie1_index      = "$commonPath/bowtie-1.1.2/RefGenomes/$genome_g/$genome_g";
my  $Bowtie2_index      = "$commonPath/bowtie2-2.2.8/RefGenomes/$genome_g/$genome_g";
my  $Pash3_index        = "0-Other/ShortCuts/$genome_g.fa";    ##reference genome is fasta file.
my  $CUSHAW3_index      = "$commonPath/cushaw3-v3.0.3/RefGenomes/$genome_g/$genome_g";
my  $Subread_index      = "$commonPath/subread-1.5.0-p1/RefGenomes/$genome_g/$genome_g";
my  $GSNAP_index        = "$commonPath/GMAP-GSNAP/RefGenomes/$genome_g/$genome_g/$genome_g";
my  $SNAP_index         = "$commonPath/SNAP/RefGenomes/$genome_g";
my  $Yara_index         = "$commonPath/yara-0.9.5/RefGenomes/$genome_g/$genome_g"; 
my  $RazerS3_index      = "0-Other/ShortCuts/$genome_g.fa";    ##reference genome is fasta file.
my  $SplazerS_index     = "0-Other/ShortCuts/$genome_g.fa";    ##reference genome is fasta file.
my  $mrsFAST_index      = "$commonPath/mrsfast-3.3.8/RefGenomes/Shortcuts/$genome_g.fa";
my  $BBMap_index        = "$commonPath/bbmap/RefGenomes/$genome_g";
my  $Novoalign_index    = "$commonPath/novocraft/RefGenomes/$genome_g/$genome_g";





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
&printVersion("bwa");
&printVersion("bwa  aln");
&printVersion("bwa  mem");
&printVersion("stampy.py  --help");
&printVersion("gem-mapper  -h");
&printVersion("bowtie    --version");
&printVersion("bowtie2   --version");
&printVersion("pash3");
&printVersion("cushaw3");
&printVersion("subread-align  -v");
&printVersion("gsnap  --version");
&printVersion("snap-aligner");
&printVersion("yara_mapper  --version");
&printVersion("razers3  --version");
&printVersion("splazers  --version");
&printVersion("mrsfast  -v");
&printVersion("bbmap.sh");
&printVersion("novoalign  --version");
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
           next unless $mapFiles[$i] =~ m/\.sam$/;
           next unless $mapFiles[$i] !~ m/^[.]/;
           next unless $mapFiles[$i] !~ m/[~]$/;
           my $temp = $mapFiles[$i]; 
           $temp =~ s/\.sam$//  ||  die; 
           say   "\t......$mapFiles[$i]";  
           system("samtools  sort  -m 2G  -o $relativePath1/$temp.bam   --output-fmt bam  -T $relativePath1/1_$temp   --threads $numCores    $relativePath1/$temp.sam    >>$SAMtools/$temp.runLog    2>&1");
           system("samtools  index           $relativePath1/$temp.bam      >>$SAMtools/$temp.runLog  2>&1");
           system("samtools  flagstat        $relativePath1/$temp.bam      >>$SAMtools/$temp.runLog  2>&1");
           system(`samtools  idxstats        $relativePath1/$temp.bam      >>$SAMtools/$temp.runLog  2>&1`);
           system("rm   $relativePath1/$temp.sam"); 
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
           system("java  -jar   $Picard   CollectAlignmentSummaryMetrics      INPUT=$relativePath1/$temp.bam   OUTPUT=$PicardDir/$temp/1_CollectAlignmentSummaryMetrics     R=0-Other/ShortCuts/$genome_g.fa                            >> $PicardDir/$temp/1.runLog   2>&1" );
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


 







## BWA, Stampy and GEM.
if(1==1) {
if ($readLen_g  eq  "no")  {

{ ########## Start bwa aln, shorter than 70bp
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using bwa aln ......";
my $BWAaln  = "$output_g/1-BWAaln";   
&myMakeDir($BWAaln);
for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {
        say    "\t......$pairedEnd[$i]";
        say    "\t......$pairedEnd[$i+1]";
        $pairedEnd[$i]   =~ m/^((\d+)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d+)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq"   eq   $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$BWAaln/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("bwa aln   -n 0.06  -o 2   -t $numCores   -R 10   -f $BWAaln/$end1.sai      $BWA_index     $input_g/$end1.fastq    >>$BWAaln/$end1.runLog   2>&1");
        system("bwa aln   -n 0.06  -o 2   -t $numCores   -R 10   -f $BWAaln/$end2.sai      $BWA_index     $input_g/$end2.fastq    >>$BWAaln/$end2.runLog   2>&1");
        system("bwa sampe -a $maxPE_g    -f $BWAaln/$temp.sam      $BWA_index     $BWAaln/$end1.sai  $BWAaln/$end2.sai     $input_g/$end1.fastq  $input_g/$end2.fastq    >>$BWAaln/$temp.runLog   2>&1"); 
}
for (my $i=0; $i<=$#singleEnd; $i++) {   
        say   "\t......$singleEnd[$i]";
        $singleEnd[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("bwa aln   -n 0.06   -o 2  -t $numCores    -f $BWAaln/$temp.sai    $BWA_index                          $input_g/$temp.fastq         >>$BWAaln/$temp.runLog   2>&1");
        system("bwa samse                                 -f $BWAaln/$temp.sam    $BWA_index     $BWAaln/$temp.sai    $input_g/$temp.fastq         >>$BWAaln/$temp.runLog   2>&1");
}
&myQC1($BWAaln);
}  ########## End bwa aln

}else{
($readLen_g  eq "yes")   or die;

{ ########## Start bwa mem, longer  than 70bp
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using BWA mem ......";
my $BWAmem  = "$output_g/1-BWAmem";   
&myMakeDir($BWAmem); 
for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {
        say    "\t......$pairedEnd[$i]";
        say    "\t......$pairedEnd[$i+1]";
        $pairedEnd[$i]   =~ m/^((\d+)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d+)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq"  eq  $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$BWAmem/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("bwa mem  -t $numCores  -T 0   $BWA_index   $input_g/$end1.fastq  $input_g/$end2.fastq    >$BWAmem/$temp.sam"); 
}
for (my $i=0; $i<=$#singleEnd; $i++) {  
        say   "\t......$singleEnd[$i]";
        $singleEnd[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("bwa mem  -t $numCores  -T 0   $BWA_index   $input_g/$temp.fastq   >$BWAmem/$temp.sam");
}
&myQC1($BWAmem);
}  ########## End bwa mem
}
}


if(1==2) { ########## Start Stampy
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using Stampy ......";
my $Stampy   = "$output_g/1-Stampy";  
&myMakeDir($Stampy); 
my $BWAfolder  = "$output_g/1-BWAmem"; 
if ( !( -e $BWAfolder) )  { $BWAfolder  = "$output_g/1-BWAaln";  }  
opendir(my $DH_BWA, $BWAfolder)  ||  die;     
my @BWAFiles = readdir($DH_BWA);
for (my $i=0; $i<=$#BWAFiles; $i++) {   
        next unless $BWAFiles[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.bam$/;  
        say   "\t......$BWAFiles[$i]";
        $BWAFiles[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.bam$/   or  die;  
        my $temp = $1; 
        system("stampy.py    -g $Stampy_index  -h $Stampy_index  -o $Stampy/$temp.sam  --threads=$numCores   --bamkeepgoodreads  -M $BWAfolder/$temp.bam    >>$Stampy/$temp.runLog  2>&1");                                                                                                                                 
}
&myQC1($Stampy);
} ########## End Stampy


if(1==2) { ########## Start GEM
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using GEM ......";
my $GEM = "$output_g/1-GEM";  
&myMakeDir($GEM); 
&myQC1($GEM);
}  ########## End GEM










## Bowtie,  CUSHAW3 and Pash3.
if(1==1) {
if ($readLen_g  eq  "no")  {
{ ########## Start Bowtie1
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using Bowtie1 ......";
my $Bowtie1   = "$output_g/2-Bowtie1";  
&myMakeDir($Bowtie1); 
for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {  
        say    "\t......$pairedEnd[$i]";
        say    "\t......$pairedEnd[$i+1]";
        $pairedEnd[$i]   =~ m/^((\d+)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d+)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq"  eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$Bowtie1/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("bowtie   --threads $numCores   -q  --chunkmbs 6400  --maxins $maxPE_g   --sam   --phred33-quals   $Bowtie1_index    -1 $input_g/$end1.fastq   -2 $input_g/$end2.fastq    $Bowtie1/$temp.sam    >>$Bowtie1/$temp.runLog   2>&1"); 
}
for (my $i=0; $i<=$#singleEnd; $i++) {   
        say   "\t......$singleEnd[$i]";
        $singleEnd[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("bowtie   --threads $numCores   -q  --chunkmbs 6400    --sam   --phred33-quals  $Bowtie1_index    $input_g/$temp.fastq   $Bowtie1/$temp.sam   >>$Bowtie1/$temp.runLog   2>&1");                
}
&myQC1($Bowtie1);
} ########## End Bowtie1

}else{
($readLen_g  eq "yes")   or die;

{ ########## Start Bowtie2
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using Bowtie2 ......";
my $Bowtie2   = "$output_g/2-Bowtie2";  
&myMakeDir($Bowtie2); 
for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {  
        say    "\t......$pairedEnd[$i]";
        say    "\t......$pairedEnd[$i+1]";
        $pairedEnd[$i]   =~ m/^((\d+)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d+)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$Bowtie2/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("bowtie2   --maxins $maxPE_g  --threads $numCores   -q   --phred33   --end-to-end    -x $Bowtie2_index    -1 $input_g/$end1.fastq        -2 $input_g/$end2.fastq     -S $Bowtie2/$temp.sam    >>$Bowtie2/$temp.runLog  2>&1");                                                                                                                                 
}
for (my $i=0; $i<=$#singleEnd; $i++) {   
        say   "\t......$singleEnd[$i]";
        $singleEnd[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("bowtie2    --threads $numCores   -q   --phred33   --end-to-end    -x $Bowtie2_index    -U $input_g/$temp.fastq       -S $Bowtie2/$temp.sam    >>$Bowtie2/$temp.runLog  2>&1");                                                                
}
&myQC1($Bowtie2);
}  ########## End Bowtie2

}
}


if(1==2) { ########## Start CUSHAW3
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using CUSHAW3 ......";
my $CUSHAW3   = "$output_g/2-CUSHAW3";  
&myMakeDir($CUSHAW3);
for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {  
        say    "\t......$pairedEnd[$i]";
        say    "\t......$pairedEnd[$i+1]";
        $pairedEnd[$i]   =~ m/^((\d+)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d+)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$CUSHAW3/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("cushaw3 align  -r $CUSHAW3_index    -q $input_g/$end1.fastq  $input_g/$end2.fastq     -o $CUSHAW3/$temp.sam   -multi 1   -t $numCores  >>$CUSHAW3/$temp.runLog  2>&1");                                                                                                                                 
                                                           
}
for (my $i=0; $i<=$#singleEnd; $i++) {   
        say   "\t......$singleEnd[$i]";
        $singleEnd[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("cushaw3 align  -r $CUSHAW3_index  -f $input_g/$temp.fastq   -o $CUSHAW3/$temp.sam    -multi 1    -t $numCores    >>$CUSHAW3/$temp.runLog  2>&1");                                                                
                                                          
}
&myQC1($CUSHAW3);
}  ########## End CUSHAW3


if(1==2) { ########## Start Pash3
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using Pash3 ......";
my $Pash3   = "$output_g/2-Pash3";  
&myMakeDir($Pash3); 

&myQC1($Pash3);
}  ########## End Pash3










## Subread, GSNAP and SNAP.
if(1==1) { ########## Start subread
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using Subread ......";
my $subread  = "$output_g/3-Subread";  
&myMakeDir($subread);
for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {  
        say   "\t......$pairedEnd[$i]";
        say   "\t......$pairedEnd[$i+1]";
        $pairedEnd[$i]   =~ m/^((\d+)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d+)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$subread/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("subread-align  -T $numCores  -I 10  -B 1  -M $MM_g   --SAMoutput  -d 0  -D $maxPE_g   -i $Subread_index   -r $input_g/$end1.fastq   -R  $input_g/$end2.fastq   -o  $subread/$temp.sam   -t 1  >>$subread/$temp.runLog  2>&1");               
}
for (my $i=0; $i<=$#singleEnd; $i++) {  
        say   "\t......$singleEnd[$i]"; 
        $singleEnd[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.fastq$/   or  die; 
        my $temp = $1; 
        system("subread-align  -T $numCores  -I 10  -B 1  -M $MM_g   --SAMoutput   -i $Subread_index    -r $input_g/$temp.fastq    -o $subread/$temp.sam   -t 1    >>$subread/$temp.runLog   2>&1");       
}
&myQC1($subread);
} ########## End subread


if(1==2) { ########## Start GSNAP
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using GSNAP ......";
my $GSNAP   = "$output_g/3-GSNAP";  
&myMakeDir($GSNAP); 

&myQC1($GSNAP);
}  ########## End GSNAP

    
if(1==2) { ########## Start SNAP
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using SNAP ......";
my $SNAP   = "$output_g/3-SNAP";  
&myMakeDir($SNAP); 

&myQC1($SNAP);
}  ########## End SNAP










##  mrsFAST, Novoalign and BBMap.
if(1==1) { ########## Start mrsFAST
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using mrsFAST ......";      
my $mrsFAST  = "$output_g/4-mrsFAST";  
&myMakeDir($mrsFAST);
my  $input_g2 = "2-FASTQ";
for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {  
        say   "\t......$pairedEnd[$i]";
        say   "\t......$pairedEnd[$i+1]";
        $pairedEnd[$i]   =~ m/^((\d+)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d+)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$mrsFAST/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("mrsfast  --search  $mrsFAST_index  --mem 8   --threads $numCores  --pe   --seq1 $input_g2/$end1.fastq  --seq2  $input_g2/$end2.fastq   -e $MM_g   --min 100  --max $maxPE_g  --best  -o $mrsFAST/$temp.sam  >>$mrsFAST/$temp.runLog  2>&1");               
}
for (my $i=0; $i<=$#singleEnd; $i++) {  
        say   "\t......$singleEnd[$i]"; 
        $singleEnd[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.fastq$/   or  die; 
        my $temp = $1; 
        system("mrsfast  --search  $mrsFAST_index  --mem 8  --threads $numCores  --seq $input_g2/$temp.fastq    -e $MM_g   --best   -o $mrsFAST/$temp.sam    >>$mrsFAST/$temp.runLog   2>&1");       
}
&myQC1($mrsFAST);
} ########## End mrsFAST


if(1==2) { ########## Start Novoalign
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using Novoalign ......";
my $Novoalign  = "$output_g/4-Novoalign";  
&myMakeDir($Novoalign);
for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {  
        say    "\t......$pairedEnd[$i]";
        say    "\t......$pairedEnd[$i+1]";
        $pairedEnd[$i]   =~ m/^((\d+)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d+)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq"  eq  $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$Novoalign/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("novoalign  -d $Novoalign_index      -f $input_g/$end1.fastq  $input_g/$end2.fastq   -i PE 250,150   -o SAM     >$Novoalign/$temp.sam ");                                                                                                                                                                                         
}
for (my $i=0; $i<=$#singleEnd; $i++) {   
        say   "\t......$singleEnd[$i]";
        $singleEnd[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("novoalign  -d $Novoalign_index      -f $input_g/$temp.fastq     -o SAM     >$Novoalign/$temp.sam ");                                                                                                                                                                                       
}
&myQC1($Novoalign);
}  ########## End Novoalign


if(1==2) { ########## Start BBMap
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using BBMap ......";
my $BBMap  = "$output_g/4-BBMap";  
&myMakeDir($BBMap);
for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {  
        say    "\t......$pairedEnd[$i]";
        say    "\t......$pairedEnd[$i+1]";
        $pairedEnd[$i]   =~ m/^((\d+)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d+)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq"  eq  $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$BBMap/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("bbmap.sh  path=$BBMap_index  in=$input_g/$end1.fastq  in2=$input_g/$end2.fastq   out=$BBMap/$temp.sam      pairlen=$maxPE_g  threads=$numCores  maxindel=10  >>$BBMap/$temp.runLog  2>&1");                                                                                                                                                                                     
}

for (my $i=0; $i<=$#singleEnd; $i++) {   
        say   "\t......$singleEnd[$i]";
        $singleEnd[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("bbmap.sh  path=$BBMap_index  in=$input_g/$temp.fastq    out=$BBMap/$temp.sam  threads=$numCores     maxindel=10   >>$BBMap/$temp.runLog  2>&1");                                                                                                                           
}
&myQC1($BBMap);
}  ########## End BBMap










## Yara, RazerS3 and SplazerS.
if(1==2) { ########## Start Yara
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using Yara ......";
my $Yara   = "$output_g/5-Yara";  
&myMakeDir($Yara);
my $pattern2 = '(\s[AGCTNagctn]+\s)|(^@)';
for (my $i=0; $i<=$#pairedEnd; $i=$i+2) {  
        say    "\n\t......$pairedEnd[$i]";
        say    "\t......$pairedEnd[$i+1]\n";
        $pairedEnd[$i]   =~ m/^((\d+)_($pattern)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd[$i+1] =~ m/^((\d+)_($pattern)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;   
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd[$i+1])  or  die;
        open(tempFH, ">>", "$Yara/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("yara_mapper  --error-rate 6  --output-secondary   --threads $numCores    $Yara_index  $input_g/$end1.fastq $input_g/$end2.fastq  | grep -P '$pattern2'   > $Yara/$temp.sam");                                                                                                                             
}
for (my $i=0; $i<=$#singleEnd; $i++) {   
        say   "\n\t......$singleEnd[$i]\n";
        $singleEnd[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))\.fastq$/   or  die;  
        my $temp = $1; 
        system("yara_mapper --error-rate 6     --output-secondary   --threads $numCores    $Yara_index  $input_g/$temp.fastq  | grep -P '$pattern2'   > $Yara/$temp.sam");                                                          
}
&myQC1($Yara);
}  ########## End Yara


if(1==2) { ########## Start RazerS3
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using RazerS3 ......";
my $RazerS3   = "$output_g/5-RazerS3";  
&myMakeDir($RazerS3); 

&myQC1($RazerS3);
}  ########## End RazerS3


if(1==2) { ########## Start SplazerS
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using SplazerS ......";
my $SplazerS   = "$output_g/5-SplazerS";  
&myMakeDir($SplazerS); 

&myQC1($SplazerS);
}  ########## End SplazerS





say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n";





## END
