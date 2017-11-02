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
my $input_g  = '';  ## such as "4-finalFASTQ"
my $output_g = '';  ## such as "5-rawBAM"

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use RASDA (RNA-Seq Data Analyzer), version 0.9.0, 2017-10-01.
        RASDA is a Pipeline for Single-end and Paired-end RNA-Seq Data Analysis by Integrating Lots of Softwares.

        Step 4: Mapping reads to the reference transcriptome or genome by using 10 softwares (mappers or aligners):
                    Kallisto, Salmon, STAR, HISAT2, RSEM; BBMap, RapMap, Subjunc, Novoalign, GSNAP.
                Assess the quality of BAM files to identify possible sequencing errors or biases by using 13 softwares:
                    SAMtools, Subread utilities, FASTQC, SAMstat, qualimap, PRESEQ, Picard, goleft, deepTools, phantompeakqualtools, QoRTs, RNA-SeQC and RSeQC.
                And aggregate the results from Kallisto, Salmon, STAR, FastQC, Picard, Samtools, Preseq, Qualimap, goleft, RNA-SeQC and RSeQC analyses
                across many samples into a single report by using MultiQC.

        Usage:
               perl  RASDA4.pl    [-version]    [-help]   [-genome RefGenome]    [-in inputDir]    [-out outDir]
        For instance:
               perl  RASDA4.pl   -genome hg38   -in 4-finalFASTQ   -out 5-rawBAM    > RASDA4.runLog  2>&1

        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.

        -help           Show this help message and exit.

        Required arguments:
        -genome RefGenome   "RefGenome" is the short name of your reference genome, such as "mm10", "ce11", "hg38".    (no default)

        -in inputDir        "inputDir" is the name of input path that contains your FASTQ files.  (no default)

        -out outDir         "outDir" is the name of output path that contains your running results (BAM files) of this step.  (no default)
        -----------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines, please visit https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ Jie Qiao Lab, yongp@outlook.com, Key Laboratory of Assisted Reproduction at Third Hospital,
        Academy for Advanced Interdisciplinary Studies, and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';

## Version Infromation
my $version = "    The Fourth Step of RASDA (RNA-Seq Data Analyzer), version 0.9.0, 2017-10-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$genome_g = 'hg38';           ## This is only an initialization value or suggesting value, not default value.
$input_g  = '4-finalFASTQ';   ## This is only an initialization value or suggesting value, not default value.
$output_g = '5-rawBAM';       ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help   -genome   -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  RASDA4.pl  -help' \n";
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

my  $trim5_g  = 30;   ## bp
my  $trim3_g  = 30;   ## bp

###################################################################################################################################################################################################





## Kallisto, Salmon, STAR, HISAT2, RSEM; BBMap, RapMap, Subjunc, Novoalign, GSNAP.
###################################################################################################################################################################################################
## Context specific:
my  $commonPath_g      = "/media/yp/ProgramFiles/.MyProgramFiles/4_ChIPseq/5-Mapping";

my  $BWA_index_g       = "$commonPath_g/bwa/RefGenomes/$genome_g/$genome_g";
my  $Bowtie2_index_g   = "$commonPath_g/bowtie2/RefGenomes/$genome_g/$genome_g";
my  $BWA_ensembl_index_g = "$commonPath_g/bwa/RefGenomes/$genome_g.ensembl/$genome_g.ensembl";
my  $Bowtie2_ensembl_index_g   = "$commonPath_g/bowtie2/RefGenomes/$genome_g.ensembl/$genome_g.ensembl";

my  $Novoalign_index_g = "$commonPath_g/novocraft/RefGenomes/$genome_g/$genome_g";
my  $Subread_index_g   = "$commonPath_g/subread/RefGenomes/$genome_g/$genome_g";
my  $GSNAP_index_g     = "RefGenomes/$genome_g/$genome_g/$genome_g";
my  $BBMap_index_g     = "/media/yp/ProgramFiles/.MyProgramFiles/4_ChIPseq/3-Remove-Correct/bbmap/RefGenomes/$genome_g";
my  $Stampy_index_g    = "$commonPath_g/stampy/RefGenomes/$genome_g/$genome_g";
my  $NGM_index_g       = "$commonPath_g/NextGenMap/RefGenomes/Shortcuts/$genome_g/$genome_g.fasta";



my  $commonPath2_g      = "/media/yp/ProgramFiles/.MyProgramFiles/5-RNAseq";

my  $Kallisto_index_g  = "$commonPath2_g/kallisto/RefGenomes/$genome_g/$genome_g.RefSeq";
my  $Salmon_index_g    = "$commonPath2_g/Salmon/RefGenomes/$genome_g/$genome_g.RefSeq";
my  $STAR_index_g      = "$commonPath2_g/STAR/RefGenomes/$genome_g";
my  $HISAT2_index_g    = "$commonPath2_g/hisat2/RefGenomes/$genome_g/$genome_g";
my  $RSEM_index_g      = "$commonPath2_g/RSEM/RefGenomes/$genome_g/$genome_g.RefSeq";
my  $RapMap_index_g    = "$commonPath2_g/RapMap/RefGenomes/$genome_g/$genome_g.RefSeq";

my  $Kallisto_ensembl_index_g  = "$commonPath2_g/kallisto/RefGenomes/$genome_g.ensembl.cDNA/$genome_g.ensembl.cDNA";
my  $Salmon_ensembl_index_g    = "$commonPath2_g/Salmon/RefGenomes/$genome_g.ensembl.cDNA/$genome_g.ensembl.cDNA";
my  $RSEM_ensembl_index_g      = "$commonPath2_g/RSEM/RefGenomes/$genome_g.ensembl.cDNA/$genome_g.ensembl.cDNA";
my  $RapMap_ensembl_index_g    = "$commonPath2_g/RapMap/RefGenomes/$genome_g.ensembl.cDNA/$genome_g.ensembl.cDNA";
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

my  $Picard_g = &fullPathApp("picard.jar");
my  $phantompeakqualtools_g = &fullPathApp("run_spp.R");
my  $QoRTs_g = &fullPathApp("QoRTs.jar");
my  $RNASeQC_g = &fullPathApp("RNA-SeQC.jar");


&printVersion("kallisto");
&printVersion("salmon");
&printVersion("STAR    --version");
&printVersion("hisat2  --version");
&printVersion("rsem-calculate-expression --version");
&printVersion("rapmap");
&printVersion("subjunc -v");
&printVersion("novoalign --version");
&printVersion("gsnap --version");
&printVersion("bbmap.sh -h");

&printVersion("samtools");
&printVersion("fastqc    -v");
&printVersion("samstat   -v");
&printVersion("Rscript  $phantompeakqualtools_g");
&printVersion("preseq");
&printVersion("qualimap  -v");
&printVersion("multiqc   --version");
&printVersion("propmapped");     ## in subread
&printVersion("qualityScores");  ## in subread
&printVersion("plotFingerprint --version");
&printVersion("goleft  -v");

&printVersion("bam_stat.py             --version");  ## in RSeQC
&printVersion("geneBody_coverage.py    --version");  ## in RSeQC
&printVersion("inner_distance.py       --version");  ## in RSeQC
&printVersion("junction_annotation.py  --version");  ## in RSeQC
&printVersion("junction_saturation.py  --version");  ## in RSeQC
&printVersion("read_distribution.py    --version");  ## in RSeQC
&printVersion("read_duplication.py     --version");  ## in RSeQC
&printVersion("RPKM_saturation.py      --version");  ## in RSeQC
&printVersion("tin.py                  --version");  ## in RSeQC

&printVersion("java  -jar  $QoRTs_g");
&printVersion("java  -jar  $RNASeQC_g");

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
    my $qualimap0  = "$QCresults/3_qualimap_BAM";
    my $qualimap1 = "$QCresults/3_qualimap_RNA1_PE";
    my $qualimap2 = "$QCresults/3_qualimap_RNA2_SE";
    my $qualimap3 = "$QCresults/3_qualimap_RNA3_PE_ensembl";
    my $qualimap4 = "$QCresults/3_qualimap_RNA4_SE_ensembl";
    my $samstat   = "$QCresults/4_samstat";

    my $MultiQC1  = "$QCresults/5_MultiQC_FastQC";
    my $MultiQC2  = "$QCresults/5_MultiQC_qualimap0";
    my $MultiQC3  = "$QCresults/5_MultiQC_SAMtools";
    my $MultiQC5A = "$QCresults/5_MultiQC_qualimap_RNA1_PE";
    my $MultiQC5B = "$QCresults/5_MultiQC_qualimap_RNA2_EE";
    my $MultiQC5C = "$QCresults/5_MultiQC_qualimap_RNA3_PE_ensembl";
    my $MultiQC5D = "$QCresults/5_MultiQC_qualimap_RNA4_SE_ensembl";

    &myMakeDir($QCresults);
    &myMakeDir($SAMtools);
    &myMakeDir($FastQC);
    &myMakeDir($qualimap0);
    &myMakeDir($qualimap1);
    &myMakeDir($qualimap2);
    #&myMakeDir($qualimap3);
    #&myMakeDir($qualimap4);
    &myMakeDir($samstat);

    &myMakeDir($MultiQC1);
    &myMakeDir($MultiQC2);
    &myMakeDir($MultiQC3);
    &myMakeDir($MultiQC5A);
    &myMakeDir($MultiQC5B);
    #&myMakeDir($MultiQC5C);
    #&myMakeDir($MultiQC5D);

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
        system("samtools  index           $dir1/$temp.bam      >>$SAMtools/$temp.index.runLog  2>&1");
        system("samtools  flagstat        $dir1/$temp.bam      >>$SAMtools/$temp.flagstat      2>&1");
        system(`samtools  idxstats        $dir1/$temp.bam      >>$SAMtools/$temp.idxstats      2>&1`);
        system( "fastqc    --outdir $FastQC    --threads $numCores_g  --format bam   --kmers 7    $dir1/$temp.bam                   >> $FastQC/$temp.runLog      2>&1" );
        system( "qualimap  bamqc   -bam $dir1/$temp.bam   -c  -ip  -nt $numCores_g   -outdir $qualimap0/$temp   --java-mem-size=16G   >> $qualimap0/$temp.runLog    2>&1" );
        system( "qualimap  rnaseq  -bam $dir1/$temp.bam   -gtf 0-Other/Shortcuts/$genome_g/$genome_g.RefSeq.GTF     -oc $temp.computed_counts   --paired   -outdir $qualimap1/$temp   --java-mem-size=16G   >> $qualimap1/$temp.runLog    2>&1" );
        system( "qualimap  rnaseq  -bam $dir1/$temp.bam   -gtf 0-Other/Shortcuts/$genome_g/$genome_g.RefSeq.GTF     -oc $temp.computed_counts              -outdir $qualimap2/$temp   --java-mem-size=16G   >> $qualimap2/$temp.runLog    2>&1" );
        #system( "qualimap  rnaseq  -bam $dir1/$temp.bam   -gtf 0-Other/Shortcuts/$genome_g/$genome_g.ensembl.cDNA.GTF     -oc $temp.computed_counts   --paired   -outdir $qualimap3/$temp   --java-mem-size=16G   >> $qualimap3/$temp.runLog    2>&1" );
        #system( "qualimap  rnaseq  -bam $dir1/$temp.bam   -gtf 0-Other/Shortcuts/$genome_g/$genome_g.ensembl.cDNA.GTF     -oc $temp.computed_counts              -outdir $qualimap4/$temp   --java-mem-size=16G   >> $qualimap4/$temp.runLog    2>&1" );
        system( "samstat   $dir1/$temp.bam      >> $samstat/$temp.runLog         2>&1");
        system( "rm   $dir1/$temp.sam" );
    }

    system( "multiqc    --title FastQC     --verbose  --export   --outdir $MultiQC1          $FastQC            >> $MultiQC1/MultiQC.FastQC.runLog     2>&1" );
    system( "multiqc    --title qualimap   --verbose  --export   --outdir $MultiQC2          $qualimap0         >> $MultiQC2/MultiQC.qualimap.runLog   2>&1" );
    system( "multiqc    --title SAMtools   --verbose  --export   --outdir $MultiQC3          $SAMtools          >> $MultiQC3/MultiQC.SAMtools.runLog   2>&1" );

    system( "multiqc    --title qualimap   --verbose  --export   --outdir $MultiQC5A         $qualimap1         >> $MultiQC5A/MultiQC.qualimap.runLog   2>&1" );
    system( "multiqc    --title qualimap   --verbose  --export   --outdir $MultiQC5B         $qualimap2         >> $MultiQC5B/MultiQC.qualimap.runLog   2>&1" );
    #system( "multiqc    --title qualimap   --verbose  --export   --outdir $MultiQC5C         $qualimap3         >> $MultiQC5C/MultiQC.qualimap.runLog   2>&1" );
    #system( "multiqc    --title qualimap   --verbose  --export   --outdir $MultiQC5D         $qualimap4         >> $MultiQC5D/MultiQC.qualimap.runLog   2>&1" );
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_BAM_2  {
    my $dir1      =  $_[0];   ## All the BAM files must be in this folder.
    my $QCresults = "$dir1/QC_Results";
    my $Fingerprint    = "$QCresults/6_Fingerprint";
    my $Fingerprint2   = "$QCresults/7_Fingerprint2";
    my $goleft         = "$QCresults/8_goleft";
    my $phantompeak    = "$QCresults/9_phantompeakqualtools";
    my $MultiQC1       = "$QCresults/10_MultiQC_goleft";

    &myMakeDir($QCresults);
    &myMakeDir($Fingerprint);
    &myMakeDir($Fingerprint2);
    &myMakeDir($goleft);
    &myMakeDir($phantompeak);
    &myMakeDir($MultiQC1);

    opendir(my $FH_Files, $dir1) || die;
    my @Files = readdir($FH_Files);

    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of all BAM files by using plotFingerprint in deepTools, goleft , phantompeakqualtools and MultiQC ......";
    for ( my $i=0; $i<=$#Files; $i++ ) {
        next unless $Files[$i] =~ m/\.bam$/;
        next unless $Files[$i] !~ m/^[.]/;
        next unless $Files[$i] !~ m/[~]$/;
        my $temp = $Files[$i];
        say    "\t......$temp";
        $temp =~ s/\.bam$//  ||  die;
        system("plotFingerprint --bamfiles $dir1/$temp.bam   --extendReads 220  --numberOfSamples 1000000    --plotFile $Fingerprint/$temp.pdf    --plotTitle $temp   --outRawCounts  $Fingerprint/$temp.cov   --outQualityMetrics $Fingerprint/$temp.Metrics.txt   --numberOfProcessors $numCores_g   --binSize 500    >> $Fingerprint/$temp.runLog    2>&1");                           
        system("plotFingerprint --bamfiles $dir1/$temp.bam   --extendReads 220  --numberOfSamples 1000000    --plotFile $Fingerprint2/$temp.pdf   --plotTitle $temp   --outRawCounts  $Fingerprint2/$temp.cov  --outQualityMetrics $Fingerprint2/$temp.Metrics.txt  --numberOfProcessors $numCores_g   --binSize 5000   >> $Fingerprint2/$temp.runLog   2>&1");                                   
        system("goleft   covstats    $dir1/$temp.bam  > $goleft/$temp.covstats " );
        system("goleft   indexcov  --sex chrX,chrY  -d $goleft/$temp  $dir1/$temp.bam  > $goleft/$temp.indexcov.runLog      2>&1" );
        &myMakeDir("$phantompeak/$temp");
        system("Rscript    $phantompeakqualtools_g    -c=$dir1/$temp.bam   -p=$numCores_g   -odir=$phantompeak/$temp    -savd=$phantompeak/$temp/rdatafile.RData     -savp=$phantompeak/$temp/plotdatafile.pdf   -out=$phantompeak/$temp/resultfile.txt   >> $phantompeak/$temp.runLog   2>&1");
    }
    system("sleep 5s");
    system( "multiqc    --title goleft    --verbose  --export   --outdir $MultiQC1          $goleft     >> $MultiQC1/MultiQC.goleft.runLog    2>&1" );

}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_BAM_3  {
    my $dir1      =  $_[0];   ## All the BAM files must be in this folder.
    my $QCresults = "$dir1/QC_Results";
    my $PRESEQ    = "$QCresults/11_PRESEQ";
    my $PicardDir = "$QCresults/12_Picard";
    my $MultiQC1  = "$QCresults/13_MultiQC_PRESEQ";
    my $MultiQC2  = "$QCresults/13_MultiQC_Picard";

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
        system("preseq  c_curve     -output  $PRESEQ/$temp.c_curve.pe.PRESEQ       -step 1000000    -verbose   -pe  -bam  $dir1/$temp.bam    >> $PRESEQ/$temp.c_curve.pe.runLog   2>&1");
        system("preseq  c_curve     -output  $PRESEQ/$temp.c_curve.se.PRESEQ       -step 1000000    -verbose        -bam  $dir1/$temp.bam    >> $PRESEQ/$temp.c_curve.se.runLog   2>&1");
        system("preseq  lc_extrap   -output  $PRESEQ/$temp.lc_extrap.pe.PRESEQ     -step 1000000    -verbose   -pe  -bam  $dir1/$temp.bam    >> $PRESEQ/$temp.lc_extrap.pe.runLog   2>&1");
        system("preseq  lc_extrap   -output  $PRESEQ/$temp.lc_extrap.se.PRESEQ     -step 1000000    -verbose        -bam  $dir1/$temp.bam    >> $PRESEQ/$temp.lc_extrap.se.runLog   2>&1");

        &myMakeDir("$PicardDir/$temp");
        #system("java  -jar   $Picard_g   CollectIndependentReplicateMetrics      INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/0_CollectIndependentReplicateMetrics     VCF=null    MINIMUM_MQ=20                                    >> $PicardDir/$temp/0.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectAlignmentSummaryMetrics          INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/1_CollectAlignmentSummaryMetrics                                                                      >> $PicardDir/$temp/1.runLog   2>&1" );
        system("java  -jar   $Picard_g   EstimateLibraryComplexity               INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/2_EstimateLibraryComplexity                                                                           >> $PicardDir/$temp/2.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectInsertSizeMetrics                INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/3_CollectInsertSizeMetrics               HISTOGRAM_FILE=$PicardDir/$temp/3.pdf  MINIMUM_PCT=0.05      >> $PicardDir/$temp/3.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectJumpingLibraryMetrics            INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/4_CollectJumpingLibraryMetrics                                                                        >> $PicardDir/$temp/4.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectMultipleMetrics                  INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/5_CollectMultipleMetrics                                                                              >> $PicardDir/$temp/5.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectBaseDistributionByCycle          INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/6_CollectBaseDistributionByCycle         CHART_OUTPUT=$PicardDir/$temp/6.pdf                          >> $PicardDir/$temp/6.runLog   2>&1" );
        system("java  -jar   $Picard_g   CollectQualityYieldMetrics              INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/7_CollectQualityYieldMetrics                                                                          >> $PicardDir/$temp/7.runLog   2>&1" );
        #system("java  -jar   $Picard_g   CollectWgsMetrics                       INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/8_CollectWgsMetricsFromQuerySorted       REFERENCE_SEQUENCE=null                                      >> $PicardDir/$temp/8.runLog   2>&1" );
        system("java  -jar   $Picard_g   MeanQualityByCycle                      INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/9_MeanQualityByCycle                     CHART_OUTPUT=$PicardDir/$temp/9.pdf                          >> $PicardDir/$temp/9.runLog   2>&1" );
        system("java  -jar   $Picard_g   QualityScoreDistribution                INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/10_QualityScoreDistribution              CHART_OUTPUT=$PicardDir/$temp/10.pdf                         >> $PicardDir/$temp/10.runLog  2>&1" );
        #system("java  -jar   $Picard_g   CollectGcBiasMetrics                    INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/11_CollectGcBiasMetrics                  CHART_OUTPUT=$PicardDir/$temp/11.pdf   SUMMARY_OUTPUT=$PicardDir/$temp/11.summary.output                  >> $PicardDir/$temp/11.runLog  2>&1" );
        #system("java  -jar   $Picard_g   CollectOxoGMetrics                      INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/12_CollectOxoGMetrics                    REFERENCE_SEQUENCE=null                                      >> $PicardDir/$temp/12.runLog  2>&1" );
        #system("java  -jar   $Picard_g   CollectSequencingArtifactMetrics        INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/13_CollectSequencingArtifactMetrics                                       >> $PicardDir/$temp/13.runLog  2>&1" );
        #system("java  -jar   $Picard_g   CollectTargetedPcrMetrics               INPUT=$dir1/$temp.bam     OUTPUT=$PicardDir/$temp/14_CollectTargetedPcrMetrics                                        >> $PicardDir/$temp/14.runLog  2>&1" );
    }
    system( "multiqc  --title PRESEQ    --verbose  --export  --outdir $MultiQC1          $PRESEQ                 >> $MultiQC1/MultiQC.PRESEQ.runLog   2>&1" );
    system( "multiqc  --title Picard    --verbose  --export  --outdir $MultiQC2          $PicardDir              >> $MultiQC2/MultiQC.Picard.runLog   2>&1" );
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_BAM_4  {
    my $dir1      =  $_[0];   ## All the BAM files must be in this folder.
    my $QCresults = "$dir1/QC_Results";
    my $SubreadUti= "$QCresults/14_SubreadUti";

    &myMakeDir("$QCresults");
    &myMakeDir("$SubreadUti");

    opendir(my $DH_map, $dir1) || die;
    my @mapFiles = readdir($DH_map);

    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of bam files by using Subreads utilities and goleft ......";
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
     }
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
sub  myQC_BAM_RNA  {
    my $dir1      =  $_[0];   ## All the BAM files must be in this folder.
    my $QCresults = "$dir1/QC_Results";
    my $QoRTs     = "$QCresults/RNA_1_QoRTs";
    my $RSeQC     = "$QCresults/RNA_2_RSeQC";
    my $RNA_SeQC  = "$QCresults/RNA_3_RNA-SeQC";
    my $MultiQC1  = "$QCresults/RNA_4_MultiQC1_RSeQC";
    my $MultiQC2  = "$QCresults/RNA_4_MultiQC2_RNA-SeQC";

    &myMakeDir($QCresults);
    &myMakeDir($QoRTs);
    &myMakeDir($RSeQC);
    &myMakeDir($RNA_SeQC);
    &myMakeDir($MultiQC1);
    &myMakeDir($MultiQC2);

    opendir(my $FH_Files, $dir1) || die;
    my @Files = readdir($FH_Files);

    say   "\n\n\n\n\n\n##################################################################################################";
    say   "Detecting the quality of all BAM files by using QoRTs, RSeQC, RNA-SeQC and MultiQC ......";
    for ( my $i=0; $i<=$#Files; $i++ ) {
        next unless $Files[$i] =~ m/\.bam$/;
        next unless $Files[$i] !~ m/^[.]/;
        next unless $Files[$i] !~ m/[~]$/;
        my $temp = $Files[$i];
        say    "\t......$temp";
        $temp =~ s/\.bam$//  ||  die;

        &myMakeDir("$QoRTs/$temp");
        system("java  -jar  $QoRTs_g  QC    --generatePlots   $dir1/$temp.bam   0-Other/Shortcuts/$genome_g/$genome_g.RefSeq.GTF   $QoRTs/$temp   >> $QoRTs/$temp.runLog   2>&1");

        &myMakeDir("$RSeQC/$temp");
        system("tin.py                  --input=$dir1/$temp.bam                                                          --refgene=0-Other/GenesBED/$genome_g/$genome_g.UCSC_knownGene.bed         >> $RSeQC/$temp/1-tin.runLog                   2>&1");
        system("bam_stat.py             --input-file=$dir1/$temp.bam                                                                                                                               >> $RSeQC/$temp/2-bam_stat.runLog              2>&1");
        system("clipping_profile.py     --input-file=$dir1/$temp.bam    --out-prefix=$RSeQC/$temp/3-clipping_profile     --sequencing=PE                                                           >> $RSeQC/$temp/3-clipping_profile.runLog      2>&1");
        system("deletion_profile.py     --input=$dir1/$temp.bam         --out-prefix=$RSeQC/$temp/4-deletion_profile     --read-align-length=150                                                   >> $RSeQC/$temp/4-deletion_profile.runLog      2>&1");
        system("geneBody_coverage.py    --input=$dir1/$temp.bam         --out-prefix=$RSeQC/$temp/5-geneBody_coverage    --refgene=0-Other/GenesBED/$genome_g/$genome_g.UCSC_knownGene.bed         >> $RSeQC/$temp/5-geneBody_coverage.runLog     2>&1");
        system("inner_distance.py       --input-file=$dir1/$temp.bam    --out-prefix=$RSeQC/$temp/6-inner_distance       --refgene=0-Other/GenesBED/$genome_g/$genome_g.UCSC_knownGene.bed         >> $RSeQC/$temp/6-inner_distance.runLog        2>&1");
        system("insertion_profile.py    --input-file=$dir1/$temp.bam    --out-prefix=$RSeQC/$temp/7-insertion_profile    --sequencing=PE                                                           >> $RSeQC/$temp/7-insertion_profile.runLog     2>&1");
        system("junction_annotation.py  --input-file=$dir1/$temp.bam    --out-prefix=$RSeQC/$temp/8-junction_annotation  --refgene=0-Other/GenesBED/$genome_g/$genome_g.UCSC_knownGene.bed         >> $RSeQC/$temp/8-junction_annotation.runLog   2>&1");
        system("junction_saturation.py  --input-file=$dir1/$temp.bam    --out-prefix=$RSeQC/$temp/9-junction_saturation  --refgene=0-Other/GenesBED/$genome_g/$genome_g.UCSC_knownGene.bed         >> $RSeQC/$temp/9-junction_saturation.runLog   2>&1");
        system("mismatch_profile.py     --input=$dir1/$temp.bam         --out-prefix=$RSeQC/$temp/10-mismatch_profile     --read-align-length=150                                                  >> $RSeQC/$temp/10-mismatch_profile.runLog     2>&1");
        system("read_distribution.py    --input-file=$dir1/$temp.bam                                                     --refgene=0-Other/GenesBED/$genome_g/$genome_g.UCSC_knownGene.bed         >> $RSeQC/$temp/11-read_distribution.runLog    2>&1");
        system("read_duplication.py     --input-file=$dir1/$temp.bam    --out-prefix=$RSeQC/$temp/11-read_duplication                                                                              >> $RSeQC/$temp/12-read_duplication.runLog     2>&1");
        system("read_GC.py              --input-file=$dir1/$temp.bam    --out-prefix=$RSeQC/$temp/12-read_GC                                                                                       >> $RSeQC/$temp/13-read_GC.runLog              2>&1");
        system("read_NVC.py             --input-file=$dir1/$temp.bam    --out-prefix=$RSeQC/$temp/13-read_NVC             --nx                                                                     >> $RSeQC/$temp/14-read_NVC.runLog             2>&1");
        system("read_quality.py         --input-file=$dir1/$temp.bam    --out-prefix=$RSeQC/$temp/14-read_quality                                                                                  >> $RSeQC/$temp/15-read_quality.runLog         2>&1");
        system("RPKM_saturation.py      --input-file=$dir1/$temp.bam    --out-prefix=$RSeQC/$temp/15-RPKM_saturation      --refgene=0-Other/GenesBED/$genome_g/$genome_g.UCSC_knownGene.bed        >> $RSeQC/$temp/16-RPKM_saturation.runLog      2>&1");

        #system("java  -jar  $RNASeQC_g  -o $RNA_SeQC/$temp  -r 0-Other/Shortcuts/$genome_g/$genome_g.fasta     -s \"$temp|$dir1/$temp.bam|$temp\"     -t  0-Other/Shortcuts/$genome_g/$genome_g.RefSeq.GTF     >> $RNA_SeQC/$temp.runLog   2>&1");
    }
    system( "multiqc  --title RSeQC      --verbose  --export  --outdir $MultiQC1       $RSeQC                 >> $MultiQC1/MultiQC.RSeQC.runLog      2>&1" );
    #system( "multiqc  --title RNA_SeQC   --verbose  --export  --outdir $MultiQC2       $RNA_SeQC              >> $MultiQC2/MultiQC.RNA_SeQC.runLog   2>&1" );
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $Kallisto_1A_g   = "$output_g/1A_Kallisto";
&myMakeDir($Kallisto_1A_g);
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
        open(tempFH,  ">>",  "$Kallisto_1A_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("kallisto quant  --threads=$numCores_g       --index=$Kallisto_index_g    --output-dir=$Kallisto_1A_g/$temp    $input_g/$end1.fastq  $input_g/$end2.fastq    >>$Kallisto_1A_g/$temp.runLog  2>&1 ");
}

for (my $i=0; $i<=$#singleEnd_g; $i++) {
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;
        my $temp = $1;
        system("kallisto quant  --threads=$numCores_g   --single    -l 250    -s 100    --index=$Kallisto_index_g    --output-dir=$Kallisto_1A_g/$temp    $input_g/$temp.fastq    >>$Kallisto_1A_g/$temp.runLog  2>&1");
}

&myMakeDir("$output2_g/1A_MultiQC_Kallisto");
system( "multiqc    --title Kallisto        --verbose  --export   --outdir $output2_g/1A_MultiQC_Kallisto        $Kallisto_1A_g     >> $output2_g/1A_MultiQC_Kallisto/MultiQC.Kallisto.runLog    2>&1" );

}  ########## End Kallisto
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $Kallisto_1B_g   = "$output_g/1B_Kallisto_ensembl";
&myMakeDir($Kallisto_1B_g);
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
        open(tempFH,  ">>",  "$Kallisto_1B_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("kallisto quant  --threads=$numCores_g       --index=$Kallisto_ensembl_index_g    --output-dir=$Kallisto_1B_g/$temp    $input_g/$end1.fastq  $input_g/$end2.fastq    >>$Kallisto_1B_g/$temp.runLog  2>&1 ");
}

for (my $i=0; $i<=$#singleEnd_g; $i++) {
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;
        my $temp = $1;
        system("kallisto quant  --threads=$numCores_g   --single    -l 250    -s 100    --index=$Kallisto_ensembl_index_g    --output-dir=$Kallisto_1B_g/$temp    $input_g/$temp.fastq    >>$Kallisto_1B_g/$temp.runLog  2>&1");
}

&myMakeDir("$output2_g/1B_MultiQC_Kallisto");
system( "multiqc    --title Kallisto        --verbose  --export   --outdir $output2_g/1B_MultiQC_Kallisto        $Kallisto_1B_g     >> $output2_g/1B_MultiQC_Kallisto/MultiQC.Kallisto.runLog    2>&1" );

}  ########## End Kallisto
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $Salmon_2A_g   = "$output_g/2A_Salmon";
&myMakeDir($Salmon_2A_g);
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
        open(tempFH, ">>", "$Salmon_2A_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("salmon quant  --threads $numCores_g   --libType IU      --index $Salmon_index_g    --output $Salmon_2A_g/$temp    --mates1 $input_g/$end1.fastq  --mates2 $input_g/$end2.fastq    >>$Salmon_2A_g/$temp.runLog  2>&1 ");
}

for (my $i=0; $i<=$#singleEnd_g; $i++) {
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;
        my $temp = $1;
        system("salmon quant  --threads $numCores_g    --libType U      --index $Salmon_index_g    --output $Salmon_2A_g/$temp    --unmatedReads $input_g/$temp.fastq      >>$Salmon_2A_g/$temp.runLog  2>&1 ");
}

&myMakeDir("$output2_g/2A_MultiQC_Salmon");
system( "multiqc    --title Salmon        --verbose  --export   --outdir $output2_g/2A_MultiQC_Salmon        $Salmon_2A_g     >> $output2_g/2A_MultiQC_Salmon/MultiQC.Salmon.runLog    2>&1" );

}  ########## End Salmon
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $Salmon_2B_g   = "$output_g/2B_Salmon_ensembl";
&myMakeDir($Salmon_2B_g);
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
        open(tempFH, ">>", "$Salmon_2B_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("salmon quant  --threads $numCores_g   --libType IU      --index $Salmon_ensembl_index_g    --output $Salmon_2B_g/$temp    --mates1 $input_g/$end1.fastq  --mates2 $input_g/$end2.fastq    >>$Salmon_2B_g/$temp.runLog  2>&1 ");
}

for (my $i=0; $i<=$#singleEnd_g; $i++) {
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;
        my $temp = $1;
        system("salmon quant  --threads $numCores_g    --libType U      --index $Salmon_ensembl_index_g    --output $Salmon_2B_g/$temp    --unmatedReads $input_g/$temp.fastq      >>$Salmon_2B_g/$temp.runLog  2>&1 ");
}

&myMakeDir("$output2_g/2B_MultiQC_Salmon");
system( "multiqc    --title Salmon        --verbose  --export   --outdir $output2_g/2B_MultiQC_Salmon        $Salmon_2B_g     >> $output2_g/2B_MultiQC_Salmon/MultiQC.Salmon.runLog    2>&1" );

}  ########## End Salmon
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $STAR_3A_g   = "$output_g/3A_STAR";
&myMakeDir($STAR_3A_g);

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
        open(tempFH, ">>", "$STAR_3A_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("STAR  --runMode alignReads    --runThreadN $numCores_g    --sjdbGTFfile 0-Other/Shortcuts/$genome_g/$genome_g.RefSeq.GTF    --quantMode   GeneCounts    --outFileNamePrefix  $STAR_3A_g/$temp.   --genomeDir $STAR_index_g  --readFilesIn $input_g/$end1.fastq  $input_g/$end2.fastq   >>$STAR_3A_g/$temp.runLog  2>&1 ");
        system("rename  s/Aligned.out.sam/sam/   $STAR_3A_g/*Aligned.out.sam");
}

for (my $i=0; $i<=$#singleEnd_g; $i++) {
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;
        my $temp = $1;
        system("STAR  --runMode alignReads   --runThreadN $numCores_g    --sjdbGTFfile 0-Other/Shortcuts/$genome_g/$genome_g.RefSeq.GTF     --quantMode   GeneCounts     --outFileNamePrefix  $STAR_3A_g/$temp.   --genomeDir $STAR_index_g  --readFilesIn $input_g/$temp.fastq   >>$STAR_3A_g/$temp.runLog  2>&1 ");
        system("rename  s/Aligned.out.sam/sam/   $STAR_3A_g/*Aligned.out.sam");
}

&myMakeDir("$output2_g/3A_MultiQC_STAR");
system( "multiqc    --title STAR        --verbose  --export   --outdir $output2_g/3A_MultiQC_STAR        $STAR_3A_g     >> $output2_g/3A_MultiQC_STAR/MultiQC.STAR.runLog    2>&1" );

}  ########## End STAR

&myQC_BAM_1($STAR_3A_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $STAR_3B_g   = "$output_g/3B_Trim_STAR";
&myMakeDir($STAR_3B_g);

{ ########## Start STAR
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using STAR ......";
my $inputDir2 = "2-mergedFASTQ";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {
        say    "\n\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]\n";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$STAR_3B_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("STAR  --runMode alignReads    --runThreadN $numCores_g    --sjdbGTFfile 0-Other/Shortcuts/$genome_g/$genome_g.RefSeq.GTF    --quantMode   GeneCounts       --clip3pNbases $trim3_g    --clip5pNbases $trim5_g    --outFileNamePrefix  $STAR_3B_g/$temp.   --genomeDir $STAR_index_g  --readFilesIn $inputDir2/$end1.fastq  $inputDir2/$end2.fastq   >>$STAR_3B_g/$temp.runLog  2>&1 ");
        system("rename  s/Aligned.out.sam/sam/   $STAR_3B_g/*Aligned.out.sam");
}

for (my $i=0; $i<=$#singleEnd_g; $i++) {
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;
        my $temp = $1;
        system("STAR  --runMode alignReads   --runThreadN $numCores_g    --sjdbGTFfile 0-Other/Shortcuts/$genome_g/$genome_g.RefSeq.GTF     --quantMode   GeneCounts        --clip3pNbases $trim3_g    --clip5pNbases $trim5_g    --outFileNamePrefix  $STAR_3B_g/$temp.   --genomeDir $STAR_index_g  --readFilesIn $inputDir2/$temp.fastq   >>$STAR_3B_g/$temp.runLog  2>&1 ");
        system("rename  s/Aligned.out.sam/sam/   $STAR_3B_g/*Aligned.out.sam");
}

&myMakeDir("$output2_g/3B_MultiQC_STAR");
system( "multiqc    --title STAR        --verbose  --export   --outdir $output2_g/3B_MultiQC_STAR        $STAR_3B_g     >> $output2_g/3B_MultiQC_STAR/MultiQC.STAR.runLog    2>&1" );

}  ########## End STAR

&myQC_BAM_1($STAR_3B_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $HISAT2_4A_g   = "$output_g/4A_HISAT2";
&myMakeDir($HISAT2_4A_g);

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
        open(tempFH, ">>", "$HISAT2_4A_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("hisat2         --threads $numCores_g   -q   --phred33   --end-to-end    -x $HISAT2_index_g    -1 $input_g/$end1.fastq        -2 $input_g/$end2.fastq     -S $HISAT2_4A_g/$temp.sam    >>$HISAT2_4A_g/$temp.runLog  2>&1");
}

for (my $i=0; $i<=$#singleEnd_g; $i++) {
        say   "\t......$singleEnd_g[$i]";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;
        my $temp = $1;
        system("hisat2        --threads $numCores_g   -q   --phred33   --end-to-end    -x $HISAT2_index_g    -U $input_g/$temp.fastq                                    -S $HISAT2_4A_g/$temp.sam    >>$HISAT2_4A_g/$temp.runLog  2>&1");
}

}  ########## End HISAT2

&myQC_BAM_1($HISAT2_4A_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $HISAT2_4B_g   = "$output_g/4B_Trim_HISAT2";
&myMakeDir($HISAT2_4B_g);

{ ########## Start HISAT2
say   "\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using HISAT2 ......";
my $inputDir2 = "2-mergedFASTQ";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {
        say    "\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$HISAT2_4B_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("hisat2         --threads $numCores_g   --sp 2,0   -q   --phred33   --end-to-end    -x $HISAT2_index_g   --trim5 $trim5_g  --trim3 $trim3_g    -1 $inputDir2/$end1.fastq        -2 $inputDir2/$end2.fastq     -S $HISAT2_4B_g/$temp.sam    >>$HISAT2_4B_g/$temp.runLog  2>&1");
}

for (my $i=0; $i<=$#singleEnd_g; $i++) {
        say   "\t......$singleEnd_g[$i]";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;
        my $temp = $1;
        system("hisat2        --threads $numCores_g   --sp 2,0   -q   --phred33   --end-to-end    -x $HISAT2_index_g   --trim5 $trim5_g  --trim3 $trim3_g    -U $inputDir2/$temp.fastq                                    -S $HISAT2_4B_g/$temp.sam    >>$HISAT2_4B_g/$temp.runLog  2>&1");
}

}  ########## End HISAT2

&myQC_BAM_1($HISAT2_4B_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $RSEM_5A_g   = "$output_g/5A_RSEM";
&myMakeDir($RSEM_5A_g);
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
        open(tempFH, ">>", "$RSEM_5A_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("rsem-calculate-expression  -p $numCores_g  --paired-end      --bowtie2    --estimate-rspd   --append-names    $input_g/$end1.fastq  $input_g/$end2.fastq   $RSEM_index_g    $RSEM_5A_g/$temp    >>$RSEM_5A_g/$temp.runLog  2>&1 ");
}

for (my $i=0; $i<=$#singleEnd_g; $i++) {
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;
        my $temp = $1;
        system("rsem-calculate-expression  -p $numCores_g  --fragment-length-mean 250    --fragment-length-sd 100     --bowtie2    --estimate-rspd   --append-names    $input_g/$temp.fastq     $RSEM_index_g    $RSEM_5A_g/$temp    >>$RSEM_5A_g/$temp.runLog  2>&1 ");
}

}  ########## End RSEM
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $RSEM_5B_g   = "$output_g/5B_RSEM_ensembl";
&myMakeDir($RSEM_5B_g);
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
        open(tempFH, ">>", "$RSEM_5B_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("rsem-calculate-expression  -p $numCores_g  --paired-end      --bowtie2    --estimate-rspd   --append-names    $input_g/$end1.fastq  $input_g/$end2.fastq   $RSEM_ensembl_index_g    $RSEM_5B_g/$temp    >>$RSEM_5B_g/$temp.runLog  2>&1 ");
}

for (my $i=0; $i<=$#singleEnd_g; $i++) {
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;
        my $temp = $1;
        system("rsem-calculate-expression  -p $numCores_g  --fragment-length-mean 250    --fragment-length-sd 100     --bowtie2    --estimate-rspd   --append-names    $input_g/$temp.fastq     $RSEM_ensembl_index_g    $RSEM_5B_g/$temp    >>$RSEM_5B_g/$temp.runLog  2>&1 ");
}

}  ########## End RSEM
###################################################################################################################################################################################################





###################################################################################################################################################################################################
&myQC_BAM_RNA($STAR_3A_g);
&myQC_BAM_RNA($STAR_3B_g);
&myQC_BAM_RNA($HISAT2_4A_g);
&myQC_BAM_RNA($HISAT2_4B_g);

&myQC_BAM_2($STAR_3A_g);
&myQC_BAM_2($STAR_3B_g);
&myQC_BAM_2($HISAT2_4A_g);
&myQC_BAM_2($HISAT2_4B_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $BBMap_g  = "$output_g/6_BBMap";
&myMakeDir($BBMap_g);
{ ########## Start BBMap
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using BBMap ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {
        say    "\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]\n";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$BBMap_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("bbmap.sh     path=$BBMap_index_g       out=$BBMap_g/$temp.sam  maxindel=200000  local=t           threads=$numCores_g   in=$input_g/$end1.fastq  in2=$input_g/$end2.fastq   -Xmx26g >>$BBMap_g/$temp.runLog   2>&1");
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;
        my $temp = $1;
        system("bbmap.sh      path=$BBMap_index_g        out=$BBMap_g/$temp.sam   maxindel=200000    local=t         threads=$numCores_g   in=$input_g/$temp.fastq   -Xmx26g  >>$BBMap_g/$temp.runLog   2>&1");
}
}  ########## End BBMap
&myQC_BAM_1($BBMap_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $RapMap_7A_g   = "$output_g/7A_RapMap";
&myMakeDir($RapMap_7A_g);
{ ########## Start RapMap
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference transcriptome by using RapMap ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {
        say    "\n\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]\n";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$RapMap_7A_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("rapmap   quasimap  --numThreads $numCores_g   --maxNumHits 5      --index $RapMap_index_g    --output $RapMap_7A_g/$temp.sam    -1 $input_g/$end1.fastq  -2 $input_g/$end2.fastq    >>$RapMap_7A_g/$temp.runLog  2>&1 ");
}

for (my $i=0; $i<=$#singleEnd_g; $i++) {
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;
        my $temp = $1;
        system("rapmap   quasimap  --numThreads $numCores_g    --maxNumHits 5      --index $RapMap_index_g    --output $RapMap_7A_g/$temp.sam    -r $input_g/$temp.fastq      >>$RapMap_7A_g/$temp.runLog  2>&1 ");
}

}  ########## End RapMap
&myQC_BAM_1($RapMap_7A_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $RapMap_7B_g   = "$output_g/7B_RapMap_ensembl";
&myMakeDir($RapMap_7B_g);
{ ########## Start RapMap
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference transcriptome by using RapMap ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {
        say    "\n\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]\n";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$RapMap_7B_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,   $end2\n";
        system("rapmap   quasimap  --numThreads $numCores_g   --maxNumHits 5      --index $RapMap_ensembl_index_g    --output $RapMap_7B_g/$temp.sam    -1 $input_g/$end1.fastq  -2 $input_g/$end2.fastq    >>$RapMap_7B_g/$temp.runLog  2>&1 ");
}

for (my $i=0; $i<=$#singleEnd_g; $i++) {
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;
        my $temp = $1;
        system("rapmap   quasimap  --numThreads $numCores_g    --maxNumHits 5      --index $RapMap_ensembl_index_g    --output $RapMap_7B_g/$temp.sam    -r $input_g/$temp.fastq      >>$RapMap_7B_g/$temp.runLog  2>&1 ");
}

}  ########## End RapMap
&myQC_BAM_1($RapMap_7B_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $subread_g  = "$output_g/8_Subjunc";
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
        system("subjunc  -T $numCores_g  -I 15  --multiMapping  -B 1  -M 6   --SAMoutput  -d 10  -D 800   -i $Subread_index_g   -r $input_g/$end1.fastq   -R  $input_g/$end2.fastq   -o  $subread_g/$temp.sam     >>$subread_g/$temp.runLog  2>&1");
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {
        say   "\t......$singleEnd_g[$i]";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;
        my $temp = $1;
        system("subjunc  -T $numCores_g  -I 15  --multiMapping  -B 1  -M 6   --SAMoutput   -i $Subread_index_g    -r $input_g/$temp.fastq    -o $subread_g/$temp.sam        >>$subread_g/$temp.runLog   2>&1");
}
} ## End subread
&myQC_BAM_1($subread_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $GSNAP_g  = "$output_g/9_GSNAP";
&myMakeDir($GSNAP_g);
{ ########## Start GSNAP
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using GSNAP ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {
        say    "\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]\n";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq" eq $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$GSNAP_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("gsnap    --db=$GSNAP_index_g    --nthreads=$numCores_g   --novelsplicing=1  --output-file=$GSNAP_g/$temp.sam   $input_g/$end1.fastq  $input_g/$end2.fastq    >>$GSNAP_g/$temp.runLog   2>&1");
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {
        say   "\n\t......$singleEnd_g[$i]\n";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;
        my $temp = $1;
        system("gsnap    --db=$GSNAP_index_g    --nthreads=$numCores_g   --novelsplicing=1  --output-file=$GSNAP_g/$temp.sam       $input_g/$temp.fastq    >>$GSNAP_g/$temp.runLog   2>&1");
}
}  ########## End GSNAP
&myQC_BAM_1($GSNAP_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $Novoalign_g  = "$output_g/10_Novoalign";
&myMakeDir($Novoalign_g);
{ ########## Start Novoalign
say   "\n\n\n\n\n\n##################################################################################################";
say   "Mapping reads to the reference genome by using Novoalign ......";
for (my $i=0; $i<=$#pairedEnd_g; $i=$i+2) {
        say    "\t......$pairedEnd_g[$i]";
        say    "\t......$pairedEnd_g[$i+1]";
        $pairedEnd_g[$i]   =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_1\.fastq$/   or  die;
        $pairedEnd_g[$i+1] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))_2\.fastq$/   or  die;
        my $temp = $1;
        my $end1 = $temp."_1";
        my $end2 = $temp."_2";
        ("$end2.fastq"  eq  $pairedEnd_g[$i+1])  or  die;
        open(tempFH, ">>", "$Novoalign_g/paired-end-files.txt")  or  die;
        say  tempFH  "$end1,  $end2\n";
        system("novoalign  -a  -r Random  -d $Novoalign_index_g      -f $input_g/$end1.fastq  $input_g/$end2.fastq    -o SAM         > $Novoalign_g/$temp.sam   2> $Novoalign_g/$temp.alignment_stats.txt ");
}
for (my $i=0; $i<=$#singleEnd_g; $i++) {
        say   "\t......$singleEnd_g[$i]";
        $singleEnd_g[$i] =~ m/^((\d+)_($pattern_g)_(Rep[1-9]))\.fastq$/   or  die;
        my $temp = $1;
        system("novoalign  -a  -r Random   -d $Novoalign_index_g      -f $input_g/$temp.fastq     -o SAM         > $Novoalign_g/$temp.sam   2> $Novoalign_g/$temp.alignment_stats.txt ");
}
}  ########## End Novoalign
&myQC_BAM_1($Novoalign_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
&myQC_BAM_RNA($BBMap_g);
&myQC_BAM_RNA($RapMap_7A_g);
&myQC_BAM_RNA($RapMap_7B_g);
&myQC_BAM_RNA($subread_g);
&myQC_BAM_RNA($GSNAP_g);
&myQC_BAM_RNA($Novoalign_g);

&myQC_BAM_2($BBMap_g);
&myQC_BAM_2($RapMap_7A_g);
&myQC_BAM_2($RapMap_7B_g);
&myQC_BAM_2($subread_g);
&myQC_BAM_2($GSNAP_g);
&myQC_BAM_2($Novoalign_g);


&myQC_BAM_3($STAR_3A_g);
&myQC_BAM_3($STAR_3B_g);
&myQC_BAM_3($HISAT2_4A_g);
&myQC_BAM_3($HISAT2_4B_g);
&myQC_BAM_3($BBMap_g);
&myQC_BAM_3($RapMap_7A_g);
&myQC_BAM_3($RapMap_7B_g);
&myQC_BAM_3($subread_g);
&myQC_BAM_3($GSNAP_g);
&myQC_BAM_3($Novoalign_g);

&myQC_BAM_4($STAR_3A_g);
&myQC_BAM_4($STAR_3B_g);
&myQC_BAM_4($HISAT2_4A_g);
&myQC_BAM_4($HISAT2_4B_g);
&myQC_BAM_4($BBMap_g);
&myQC_BAM_4($RapMap_7A_g);
&myQC_BAM_4($RapMap_7B_g);
&myQC_BAM_4($subread_g);
&myQC_BAM_4($GSNAP_g);
&myQC_BAM_4($Novoalign_g);

###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";





## END
