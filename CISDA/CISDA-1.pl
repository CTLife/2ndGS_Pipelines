#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.20;    ## perl5 version >= 5.20


                      


###################################################################################################################################################################################################
###################################################################################################################################################################################################


########## Help Infromation ##########
my $HELP_g = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use CISDA (ChIP-Seq Data Analyzer), version 0.7.0, 2016-03-05.
        CISDA is a Pipeline for Single-end and Paired-end ChIP-Seq Data Analysis by Integrating Lots of Softwares.

        Step 1: Extract all compressed FASTQ files, or convert SRA to FASTQ by using SRA_Toolkit,
                and merge the two lanes of the same sample (merge technical replicates).
                And  assess the quality of the raw reads to identify possible sequencing errors or biases by using
                Rqc, FastQC, kPAL, NGS_QC_Toolkit, FASTX-toolkit, FaQCs, prinseq, fastq-stats in EA-Utils, SolexaQA,
                Raspberry, ht2-stat in HTQC, SGA, fastqp,  ShortRead and seqTools.

        Usage:
               perl  CISDA-1.pl    [-v]    [-h]    [-i inputDir]    [-o outDir]
        For instance:
                     perl  CISDA-1.pl    -i 1-rawReads          -o 2-FASTQ
                     perl  CISDA-1.pl    --input 1-rawReads     --output 2-FASTQ
                     perl  CISDA-1.pl    --input 1-rawReads     --output 2-FASTQ    >> CISDA-1.runLog  2>&1

        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -v, --version        Show version number of this program and exit.

        -h, --help           Show this help message and exit.

        Required arguments:
        -i inputDir,  --input inputDir        "inputDir" is the name of input folder that contains your compressed
                                              fastq files or SRA files.
                                              The suffix of the compressed fastq files will be recognized. (no default)

        -o outDir,    --output outDir         "outDir" is the name of output folder that contains your running
                                              results (fastq format) of this step.      (no default)
        -----------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines such as RASDA, MESDA and HISDA,
        please see the web site: https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies
        and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';


########## Version Infromation ##########
my $version_g = "    The First Step of CISDA (ChIP-Seq Data Analyzer), version 0.7.0, 2016-03-05.";


########## Keys and Values ##########
if ($#ARGV   == -1)   { say  "\n$HELP_g\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-h") ;       }       ## when the number of command argumants is odd.
my %args = @ARGV;


########## Initialize  Variables ##########
my $input_g  = '1-rawReads';     ## This is only an initialization value or suggesting value, not default value.
my $output_g = '2-FASTQ';        ## This is only an initialization value or suggesting value, not default value.


########## Available Arguments ##########
my $available = "      -v  --version      -h  --help      -i  --input       -o    --output      ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole_g = 1; }
}
if($boole_g == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  CISDA-1.pl  -h' \n";
    exit 0;
}


########## Get Arguments ##########
if ( ( exists $args{'-v' } )  or  ( exists $args{'--version' } )  )     { say  "\n$version_g\n";    exit 0; }
if ( ( exists $args{'-h' } )  or  ( exists $args{'--help'    } )  )     { say  "\n$HELP_g\n";       exit 0; }
if ( ( exists $args{'-i' } )  or  ( exists $args{'--input'   } )  )     { ($input_g  = $args{'-i' })  or  ($input_g  = $args{'--input' });  }else{say   "\n -i or --input  is required.\n";   say  "\n$HELP_g\n";    exit 0; }                         
if ( ( exists $args{'-o' } )  or  ( exists $args{'--output'  } )  )     { ($output_g = $args{'-o' })  or  ($output_g = $args{'--output'});  }else{say   "\n -o or --output is required.\n";   say  "\n$HELP_g\n";    exit 0; }                             


########### Conditions #############
$input_g  =~ m/^\S+$/    ||  die   "\n\n$HELP_g\n\n";
$output_g =~ m/^\S+$/    ||  die   "\n\n$HELP_g\n\n";


######### Print Command Arguments to Standard Output ###########
say  "\n
        ################ Your Arguments ###############################
                Input  Folder:  $input_g
                Output Folder:  $output_g
        ###############################################################
\n";


###################################################################################################################################################################################################
###################################################################################################################################################################################################










say   "\n\n##################################################################################################";
say   "Running ......";
sub myMakeDir  {
    my $path = $_[0];
    if ( !( -e $path) )  { mkdir $path  ||  die; }
}
my $input2_g  = "$input_g/Results";
my $output2_g = "$output_g/Results";
&myMakeDir("$input2_g");
&myMakeDir("$output_g");
&myMakeDir("$output2_g");
opendir(my $DH_input, $input_g)  ||  die;
my @inputFiles = readdir($DH_input);
my $pattern = "[-.0-9A-Za-z]+";
my $numCores = 8;





## These commands must be available:
say   "\n\n##################################################################################################";
say   "Checking all the necessary softwares in this step ......";
sub printVersion  {
    my $software = $_[0];
    system("echo    '##############################################################################'  >> $output2_g/VersionsOfSoftwares.txt   2>&1");
    system("echo    '#########$software'                                                              >> $output2_g/VersionsOfSoftwares.txt   2>&1");
    system("$software                                                                                 >> $output2_g/VersionsOfSoftwares.txt   2>&1");
    system("echo    '\n\n\n\n\n\n'                                                                    >> $output2_g/VersionsOfSoftwares.txt   2>&1");
}
&printVersion(" fastq-dump  --version ");
&printVersion(" bzip2  --version ");
&printVersion(" gunzip  --version ");
&printVersion(" tar  --version ");
&printVersion(" unrar   ");
&printVersion(" xz    --version ");
&printVersion(" unzip -v ");
&printVersion(" fastqc    -v ");
&printVersion(" kpal  -v ");
&printVersion(" IlluQC_PRLL.pl  -h ");
&printVersion(" fastx_quality_stats  -h ");
&printVersion(" fastq_quality_boxplot_graph.sh   -h ");
&printVersion(" fastx_nucleotide_distribution_graph.sh  -h ");
&printVersion(" FaQCs.pl   --version ");
&printVersion(" prinseq-lite.pl  -version ");
&printVersion(" fastq-stats -v ");
&printVersion(" SolexaQA ");
&printVersion(" raspberry ");
&printVersion(" sga ");
&printVersion(" fastools  -v ");
&printVersion(" ht2-stat  -v ");
&printVersion(" fastqp    -h ");





say   "\n\n##################################################################################################";
say   "Checking all the input file names ......";
my @groupFiles = ();
my $fileNameBool = 1;
for ( my $i=0; $i<=$#inputFiles; $i++ ) {
        next unless $inputFiles[$i] !~ m/^[.]/;
        next unless $inputFiles[$i] !~ m/[~]$/;
        next unless $inputFiles[$i] !~ m/^Results$/;
        say   "\t......$inputFiles[$i]" ;
        my $temp = $inputFiles[$i];
        $groupFiles[++$#groupFiles] = $inputFiles[$i];
        $temp =~ m/^(\d+)_($pattern)_(Rep[1-9])/    or  die   "wrong-1: ## $temp ##";
        $temp =~ m/_(Rep[1-9])\.sra$/  or  $temp =~ m/_(Rep[1-9])_?([1-2]?)(_Lane[1-2])?\.fastq\.(\S+)$/  or  die   "wrong-2: ## $temp ##";
        if($temp !~ m/^((\d+)_($pattern)_(Rep[1-9]))(_[1-2])?(_Lane[1-2])?(\.fastq)?\.(\S+)$/) {
             $fileNameBool = 0;
        }
}
if($fileNameBool == 1)  { say    "\n\t\tAll the file names are passed.\n";  }
@groupFiles = sort(@groupFiles);
my $numGroup = 0;
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





say   "\n\n##################################################################################################";
say   "Converting SRA files into FASTQ files or Extracting the compressed fastq files ......";
for ( my $i=0; $i<=$#inputFiles; $i++ ) {
        next unless $inputFiles[$i] !~ m/^[.]/;
        next unless $inputFiles[$i] !~ m/[~]$/;
        next unless $inputFiles[$i] !~ m/^Results$/;
        say  "\t......$inputFiles[$i]" ;
        my $temp = $inputFiles[$i];
        if ($temp =~ m/\.sra$/) {
                $temp =~ m/^(\d+)_($pattern)_(Rep[1-9])\.sra$/   or  die;
                system("fastq-dump   --split-3   --dumpbase   --outdir $input_g    $input_g/$temp   >> $input2_g/$temp.runLog  2>&1");
        }else{
                $temp =~ m/^((\d+)_($pattern)_(Rep[1-9])_?([1-2]?)(_Lane[1-2])?\.fastq)\.(\S+)$/   or  die;
                my  $tempFastq = $1;
                my  $suffix1   = $7;    ## Only seven compressed formats are supported, their suffixes:  ".bz2",  ".gz",  ".tar.gz",  ".tar",  ".rar",  ".xz",  ".zip"
                my  $tempBool  = 0;
                if($suffix1 eq "bz2"    )  { $tempBool++;  system("bzip2   -cd         $input_g/$temp   >  $input_g/$tempFastq");  }
                if($suffix1 eq "gz"     )  { $tempBool++;  system("gunzip  --stdout    $input_g/$temp   >  $input_g/$tempFastq");  }
                if($suffix1 eq "tar.gz" )  { $tempBool++;  system("tar     -xzvf       $input_g/$temp  -C  $input_g");             }
                if($suffix1 eq "tar"    )  { $tempBool++;  system("tar     -xf         $input_g/$temp  -C  $input_g");             }
                if($suffix1 eq "rar"    )  { $tempBool++;  system("unrar    e          $input_g/$temp      $input_g");             }
                if($suffix1 eq "xz"     )  { $tempBool++;  system("xz      -cd         $input_g/$temp   >  $input_g/$tempFastq");  }
                if($suffix1 eq "zip"    )  { $tempBool++;  system("unzip   -np         $input_g/$temp   >  $input_g/$tempFastq");  }
                if($tempBool  != 1)        { say("$temp is wrong!!");  die; }
        }
}





say   "\n\n##################################################################################################";
say   "Merge the two lanes of the same sample or copy files ......";
my $FastQCdir1       = "$input2_g/FastQC";  
&myMakeDir("$FastQCdir1");
opendir($DH_input, $input_g) || die;
my @twoLanes = readdir($DH_input);
for (my $i=0; $i<=$#twoLanes; $i++) {
    next unless $twoLanes[$i] !~ m/^Results$/;
    next unless $twoLanes[$i] =~ m/\.fastq$/;
    next unless $twoLanes[$i] !~ m/^[.]/;
    next unless $twoLanes[$i] !~ m/[~]$/;
    say    "\t......$twoLanes[$i]";
    $twoLanes[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])_?([1-2]?)(_Lane[1-2])?\.fastq$/   or  die;
    if( $twoLanes[$i] =~ m/^(\S+)_Lane1.fastq$/ ) {
        my $temp = $1;
        $temp =~ m/^(\d+)_($pattern)_(Rep[1-9])_?([1-2]?)$/   or  die;
        open(tempFH, ">>", "$output2_g/Merge-Two-Lanes.log")  or  die;
        my $lane1 = $temp."_Lane1.fastq";
        my $lane2 = $temp."_Lane2.fastq";
        system(       "cat  $input_g/$lane1  $input_g/$lane2   > $output_g/$temp.fastq" );
        say  tempFH   "cat  $input_g/$lane1  $input_g/$lane2   > $output_g/$temp.fastq";
    }
    if ($twoLanes[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])_?([1-2]?).fastq$/) {
        system(  "cp   $input_g/$twoLanes[$i]    $output_g/$twoLanes[$i]" );
    }
    system( "fastqc    --outdir $FastQCdir1   --threads $numCores    --kmers 7     $input_g/$twoLanes[$i]       >> $FastQCdir1/$twoLanes[$i].runLog   2>&1" );
}
system( "Rscript  0-Other/Rsrc/Rqc.R   $input_g    $input2_g   >> $input2_g/Rqc.runLog   2>&1" );
system( "rm    $input_g/*.fastq" );





say   "\n\n##################################################################################################";
say   "Detecting single-end and paired-end FASTQ files ......";
opendir(my $DH_output, $output_g) || die;
my @outputFiles = readdir($DH_output);
my @singleEnd = ();
my @pairedEnd = ();
open(seqFiles_FH, ">", "$output2_g/singleEnd-pairedEnd-Files.txt")  or  die;
for ( my $i=0; $i<=$#outputFiles; $i++ ) {
    next unless $outputFiles[$i] =~ m/\.fastq$/;
    next unless $outputFiles[$i] !~ m/^[.]/;
    next unless $outputFiles[$i] !~ m/[~]$/;
    next unless $outputFiles[$i] !~ m/^unpaired/;
    next unless $outputFiles[$i] !~ m/^Results$/;
    say    "\t......$outputFiles[$i]";
    $outputFiles[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])_?([1-2]?)\.fastq$/   or  die;
    if ($outputFiles[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])\.fastq$/) {   ## sinlge end sequencing files.
        $outputFiles[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])\.fastq$/  or  die;
        $singleEnd[$#singleEnd+1] =  $outputFiles[$i];
        say  "\t\t\t\tSingle-end sequencing files: $outputFiles[$i]\n";
        say seqFiles_FH  "Single-end sequencing files: $outputFiles[$i]\n";
    }else{     ## paired end sequencing files.
        $outputFiles[$i] =~ m/^(\d+)_($pattern)_(Rep[1-9])_([1-2])\.fastq$/  or  die;
        if ($outputFiles[$i] =~ m/^((\d+)_($pattern)_(Rep[1-9]))_1\.fastq$/) { ## The two files of one paired sequencing sample are always side by side.
            my $temp = $1;
            my $end1 = $temp."_1.fastq";
            my $end2 = $temp."_2.fastq";
            (-e  "$output_g/$end1")  or die;
            (-e  "$output_g/$end2")  or die;
            $pairedEnd[$#pairedEnd+1] =  $end1;
            $pairedEnd[$#pairedEnd+1] =  $end2;
            say  "\t\t\t\tPaired-end sequencing files: $end1,  $end2\n";
            say seqFiles_FH  "Paired-end sequencing files: $end1,  $end2\n";
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





my $R_QC         = "$output2_g/1-Rqc";
my $FastQC       = "$output2_g/2-FastQC";
my $kPAL         = "$output2_g/3-kPAL";
my $NGSQCToolkit = "$output2_g/4-NGSQCToolkit";
my $FASTXtoolkit = "$output2_g/5-FASTXtoolkit";
my $FaQCs        = "$output2_g/6-FaQCs";
my $PRINSEQ      = "$output2_g/7-PRINSEQ";
my $EA_Utils     = "$output2_g/8-EA-Utils";
my $SolexaQA     = "$output2_g/9-SolexaQA";
my $Raspberry    = "$output2_g/10-raspberry";
my $HTQC         = "$output2_g/11-HTQC";
my $SGA          = "$output2_g/12-SGA";
my $Fastqp       = "$output2_g/13-fastqp";
my $ShortRead    = "$output2_g/14-ShortRead";
my $seqTools     = "$output2_g/15-seqTools";
&myMakeDir("$R_QC");
&myMakeDir("$FastQC");
&myMakeDir("$kPAL");
&myMakeDir("$NGSQCToolkit");
&myMakeDir("$FASTXtoolkit");
&myMakeDir("$FaQCs");
&myMakeDir("$PRINSEQ");
&myMakeDir("$EA_Utils");
&myMakeDir("$SolexaQA");
&myMakeDir("$Raspberry");
&myMakeDir("$HTQC");
&myMakeDir("$SGA");
&myMakeDir("$Fastqp");
&myMakeDir("$ShortRead");
&myMakeDir("$seqTools");





say   "\n\n##################################################################################################";
say   "Detecting the quality of all FASTQ files by using Rqc, FastQC, fastq-stats in EA-Utils, SolexaQA, raspberry, fastqp, FASTX-toolkit, ShortRead and seqTools ......";
system( "Rscript  0-Other/Rsrc/Rqc.R         $output_g    $R_QC        >> $R_QC/Rqc.runLog              2>&1" );
system( "Rscript  0-Other/Rsrc/ShortRead.R   $output_g    $ShortRead   >> $ShortRead/ShortRead.runLog   2>&1" );
system( "Rscript  0-Other/Rsrc/seqTools.R    $output_g    $seqTools    >> $seqTools/seqTools.runLog     2>&1" );
for ( my $i=0; $i<=$#outputFiles; $i++ ) {
    next unless $outputFiles[$i] =~ m/\.fastq$/;
    next unless $outputFiles[$i] !~ m/^[.]/;
    next unless $outputFiles[$i] !~ m/[~]$/;
    my $temp = $outputFiles[$i];
    say    "\t......$temp";
    $temp =~ m/^(\d+)_($pattern)_(Rep[1-9])_?([1-2]?)\.fastq$/   or  die;
    $temp =~ s/\.fastq$//  ||  die;
    &myMakeDir("$FastQC/$temp");
    &myMakeDir("$FASTXtoolkit/$temp");
    &myMakeDir("$EA_Utils/$temp");
    &myMakeDir("$SolexaQA/$temp");
    system( "fastqc    --outdir $FastQC/$temp    --threads $numCores    --kmers 7    $output_g/$temp.fastq       >> $FastQC/$temp.runLog     2>&1" );
    system( "fastx_quality_stats                        -i $output_g/$temp.fastq                          -o $FASTXtoolkit/$temp/$temp.txt            >> $FASTXtoolkit/$temp.runLog   2>&1" );
    system( "fastq_quality_boxplot_graph.sh             -i $FASTXtoolkit/$temp/$temp.txt     -t $temp     -o $FASTXtoolkit/$temp/$temp.quality.png    >> $FASTXtoolkit/$temp.runLog   2>&1" );
    system( "fastx_nucleotide_distribution_graph.sh     -i $FASTXtoolkit/$temp/$temp.txt     -t $temp     -o $FASTXtoolkit/$temp/$temp.nucDis.png     >> $FASTXtoolkit/$temp.runLog   2>&1" );
    system( "fastq-stats  -d  -s 100   -x $EA_Utils/$temp/$temp.fastxStatistics    -b $EA_Utils/$temp/$temp.baseBreakdown    -L $EA_Utils/$temp/$temp.lengthCounts  $output_g/$temp.fastq   >> $EA_Utils/$temp.runLog   2>&1" );
    system( "SolexaQA  analysis   $output_g/$temp.fastq    --minmax    --directory $SolexaQA/$temp       >> $SolexaQA/$temp.runLog    2>&1" );
    system( "raspberry  -t $numCores   $output_g/$temp.fastq                                             >> $Raspberry/$temp.runLog   2>&1" );
    system( "fastqp    --nreads 20000000   --kmer 5    --output $Fastqp/$temp  --type fastq   --median-qual 30     $output_g/$temp.fastq     >> $Fastqp/$temp.runLog     2>&1 " );
    system( "mv  $output_g/*.rlen   $Raspberry" );
}





say   "\n\n##################################################################################################";
say   "Detecting the quality of single-end FASTQ files by using kPAL, NGS_QC_Toolkit, FaQCs, prinseq, ht2-stat in HTQC and SGA ......";
for ( my $i=0; $i<=$#singleEnd; $i++ ) {
    my $temp = $singleEnd[$i];
    say   "\t......$singleEnd[$i]";
    $temp =~ m/^(\d+)_($pattern)_(Rep[1-9])\.fastq$/   or  die;
    $temp =~ s/\.fastq$//  ||  die;
    &myMakeDir("$NGSQCToolkit/$temp");
    &myMakeDir("$FaQCs/$temp");
    &myMakeDir("$PRINSEQ/$temp");
    &myMakeDir("$HTQC/$temp");
    &myMakeDir("$SGA/$temp");
    system( "fastools  fq2fa         $output_g/$temp.fastq        $kPAL/$temp.fa               >> $kPAL/$temp.runLog   2>&1");
    system( "IlluQC_PRLL.pl      -se $output_g/$temp.fastq     N    A      -cpus $numCores  -onlyStat     -outputFolder $NGSQCToolkit/$temp           >> $NGSQCToolkit/$temp.runLog   2>&1" );
    system( "FaQCs.pl     -prefix $temp     -t $numCores      -qc_only      -d $FaQCs/$temp          -u $output_g/$temp.fastq                         >> $FaQCs/$temp.runLog          2>&1" );
    system( "prinseq-lite.pl  -out_format 1 -verbose  -fastq $output_g/$temp.fastq    -graph_data $PRINSEQ/$temp/$temp.gd      >> $PRINSEQ/$temp.runLog   2>&1" );
    system( "prinseq-graphs.pl   -i $PRINSEQ/$temp/$temp.gd   -png_all    -o $PRINSEQ/$temp/$temp                 >> $PRINSEQ/$temp.runLog   2>&1" );
    system( "prinseq-graphs.pl   -i $PRINSEQ/$temp/$temp.gd   -html_all   -o $PRINSEQ/$temp/$temp                 >> $PRINSEQ/$temp.runLog   2>&1" );
    system( "prinseq-lite.pl  -out_format 1  -verbose  -fastq $output_g/$temp.fastq    -stats_all                               >> $PRINSEQ/$temp.stats_all   2>&1" );
    system( "ht2-stat  -i $output_g/$temp.fastq    -o  $HTQC/$temp   --se     --encode sanger   --threads $numCores  >> $HTQC/$temp.runLog   2>&1" );
    system( "sga  preprocess  --verbose    --out=$SGA/$temp/$temp.preprocess.fastq   --pe-mode=0  $output_g/$temp.fastq      >> $SGA/$temp.runLog   2>&1" );
    system( "sga index        --verbose  -a ropebwt   -t $numCores   $SGA/$temp/$temp.preprocess.fastq                     >> $SGA/$temp.runLog   2>&1" );
    system( "cp  *.preprocess.bwt      $SGA/$temp " );
    system( "cp  *.preprocess.rbwt     $SGA/$temp " );
    system( "cp  *.preprocess.sai      $SGA/$temp " );
    system( "cp  *.preprocess.rsai     $SGA/$temp " );
    system( "sga preqc         -t $numCores   $SGA/$temp/$temp.preprocess.fastq     > $SGA/$temp/$temp.preqc                             " );
    system( "sga-preqc-report.py  -o $SGA/$temp/$temp   --page_per_plot  $SGA/$temp/$temp.preqc                                         >> $SGA/$temp.runLog   2>&1" );
    system( "sga stats   --verbose   --threads=$numCores   --run-lengths     --kmer-distribution    $SGA/$temp/$temp.preprocess.fastq   >> $SGA/$temp.stats    2>&1" );
    system( "rm  $SGA/$temp/*.preprocess.* " );
    system( "rm  *.preprocess.* " );
    system( "rm   $output_g/*_prinseq_* " );
}





say   "\n\n##################################################################################################";
say   "Detecting the quality of paired-end FASTQ files by using kPAL, NGS_QC_Toolkit, FaQCs, prinseq, ht2-stat in HTQC and SGA ......";
for ( my $j=0; $j<=$#pairedEnd; $j=$j+2 ) {
    my $temp1 = $pairedEnd[$j];
    my $temp2 = $pairedEnd[$j+1];
    my $temp  = $temp2;
    say   "\t......$temp1";
    say   "\t......$temp2";
    $temp1 =~ m/^(\d+)_($pattern)_(Rep[1-9])_1\.fastq$/   or  die;
    $temp2 =~ m/^(\d+)_($pattern)_(Rep[1-9])_2\.fastq$/   or  die;
    $temp1 =~ s/\.fastq$//    ||  die;
    $temp2 =~ s/\.fastq$//    ||  die;
    $temp  =~ s/_2\.fastq$//  ||  die;
    &myMakeDir("$NGSQCToolkit/$temp");  
    &myMakeDir("$FaQCs/$temp");
    &myMakeDir("$PRINSEQ/$temp");
    &myMakeDir("$HTQC/$temp");
    &myMakeDir("$SGA/$temp");
    system( "fastools  fq2fa         $output_g/$temp1.fastq        $kPAL/$temp1.fa               >> $kPAL/$temp.runLog   2>&1");
    system( "fastools  fq2fa         $output_g/$temp2.fastq        $kPAL/$temp2.fa               >> $kPAL/$temp.runLog   2>&1");
    system( "IlluQC_PRLL.pl      -pe $output_g/$temp1.fastq  $output_g/$temp2.fastq   N    A      -cpus $numCores     -onlyStat          -outputFolder $NGSQCToolkit/$temp           >> $NGSQCToolkit/$temp.runLog   2>&1" );
    system( "FaQCs.pl     -prefix $temp     -t $numCores      -qc_only      -d $FaQCs/$temp          -p $output_g/$temp1.fastq  $output_g/$temp2.fastq   >> $FaQCs/$temp.runLog          2>&1" );
    system( "prinseq-lite.pl   -out_format 1   -verbose  -fastq $output_g/$temp1.fastq  -fastq2 $output_g/$temp2.fastq   -graph_data $PRINSEQ/$temp/$temp.gd      >> $PRINSEQ/$temp.runLog   2>&1" );
    system( "prinseq-graphs.pl   -i $PRINSEQ/$temp/$temp.gd   -png_all    -o $PRINSEQ/$temp/$temp                 >> $PRINSEQ/$temp.runLog   2>&1" );
    system( "prinseq-graphs.pl   -i $PRINSEQ/$temp/$temp.gd   -html_all   -o $PRINSEQ/$temp/$temp                 >> $PRINSEQ/$temp.runLog   2>&1" );
    system( "prinseq-lite.pl   -out_format 1  -verbose  -fastq $output_g/$temp.fastq    -stats_all                               >> $PRINSEQ/$temp.stats_all   2>&1" );
    system( "ht2-stat  -i $output_g/$temp1.fastq $output_g/$temp2.fastq   -o  $HTQC/$temp    --pe     --encode sanger   --threads 6  >> $HTQC/$temp.runLog   2>&1" );
    system( "sga  preprocess  --verbose  --out=$SGA/$temp/$temp.preprocess.fastq   --pe-mode=1  $output_g/$temp1.fastq   $output_g/$temp2.fastq   >> $SGA/$temp.runLog   2>&1" );
    system( "sga index        --verbose  -a ropebwt   -t $numCores   $SGA/$temp/$temp.preprocess.fastq                             >> $SGA/$temp.runLog   2>&1" );
    system( "cp  *.preprocess.bwt      $SGA/$temp " );
    system( "cp  *.preprocess.rbwt     $SGA/$temp " );
    system( "cp  *.preprocess.sai      $SGA/$temp " );
    system( "cp  *.preprocess.rsai     $SGA/$temp " );
    system( "sga preqc       -t $numCores   $SGA/$temp/$temp.preprocess.fastq     > $SGA/$temp/$temp.preqc                             " );
    system( "sga-preqc-report.py  -o $SGA/$temp/$temp   --page_per_plot  $SGA/$temp/$temp.preqc                            >> $SGA/$temp.runLog   2>&1" );
    system( "sga stats        --verbose   --threads=$numCores   --run-lengths     --kmer-distribution    $SGA/$temp/$temp.preprocess.fastq             >> $SGA/$temp.stats   2>&1" );
    system( "rm  $SGA/$temp/*.preprocess.* " );
    system( "rm  *.preprocess.* " );
    system( "rm   $output_g/*_prinseq_* " );
}
system( "kpal   count  -k 8   $kPAL/*.fa    $kPAL/rawReads.k8   >> $kPAL/rawReads.k8.runLog   2>&1");
system( "kpal   count  -k 6   $kPAL/*.fa    $kPAL/rawReads.k6   >> $kPAL/rawReads.k6.runLog   2>&1");
system( "kpal   count  -k 4   $kPAL/*.fa    $kPAL/rawReads.k4   >> $kPAL/rawReads.k4.runLog   2>&1");
system( "kpal   count  -k 2   $kPAL/*.fa    $kPAL/rawReads.k2   >> $kPAL/rawReads.k2.runLog   2>&1");
system( "kpal   info   $kPAL/rawReads.k8   >> $kPAL/rawReads.k8.runLog   2>&1");
system( "kpal   info   $kPAL/rawReads.k6   >> $kPAL/rawReads.k6.runLog   2>&1");
system( "kpal   info   $kPAL/rawReads.k4   >> $kPAL/rawReads.k4.runLog   2>&1");
system( "kpal   info   $kPAL/rawReads.k2   >> $kPAL/rawReads.k2.runLog   2>&1");
system( "kpal   matrix   $kPAL/rawReads.k8   $kPAL/rawReads.k8.matrix   >> $kPAL/rawReads.k8.runLog   2>&1");
system( "kpal   matrix   $kPAL/rawReads.k6   $kPAL/rawReads.k6.matrix   >> $kPAL/rawReads.k6.runLog   2>&1");
system( "kpal   matrix   $kPAL/rawReads.k4   $kPAL/rawReads.k4.matrix   >> $kPAL/rawReads.k4.runLog   2>&1");
system( "kpal   matrix   $kPAL/rawReads.k2   $kPAL/rawReads.k2.matrix   >> $kPAL/rawReads.k2.runLog   2>&1");
##system( "rm   $kPAL/*.fa " );
my $R_QC2 = "$R_QC"."-Final";
&myMakeDir($R_QC2);
system( "Rscript     0-Other/Rsrc/Rqc.R         $output_g    $R_QC2        >> $R_QC2/Rqc2.runLog   2>&1" );





say   "\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n";




