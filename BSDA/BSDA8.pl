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
my $input_g  = '';  ## such as "8-bedGraph-Bismark"
my $output_g = '';  ## such as "9-Me-BigWig-Bismark"

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use BSDA (BS-Seq Data Analyzer), version 0.9.4,  2018-02-01.
        BSDA is a Pipeline for Single-end and Paired-end BS-Seq Data Analysis by Integrating Lots of Softwares.
                                                            
        Step 8: This script generates a cytosine methylation report for a genome of interest.
                Convert bedGraph to bigwig by using bedGraphToBigWig.

                If this script works well, you do not need to check the the versions of the softwares or packages whcih are used in this script.
                And you do not need to exactly match the versions of the softwares or packages.
                If some errors or warnings are reported, please check the versions of softwares or packages.

                The versions of softwares or packages are used in this script:  
                        Perl,      5.22.1
                        coverage2cytosine, 0.19.0   
                        bedGraphToBigWig, v4  

        Usage:
               perl  BSDA8.pl    [-version]    [-help]   [-genome RefGenome]    [-in inputDir]    [-out outDir]
        For instance:
               perl  BSDA8.pl   -genome hg38   -in 8-bedGraph-Bismark   -out 9-Me-BigWig-Bismark    > BSDA8.runLog

        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.

        -help           Show this help message and exit.

        Required arguments:
        -genome RefGenome   "RefGenome" is the short name of your reference genome, such as "mm10", "ce11", "hg38".    (no default)

        -in inputDir        "inputDir" is the name of input path that contains your BAM files.  (no default)

        -out outDir         "outDir" is the name of output path that contains your running results of this step.  (no default)
        -----------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines, please visit https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ Jie Qiao Lab, yongp@outlook.com, Key Laboratory of Assisted Reproduction at Third Hospital,
        Academy for Advanced Interdisciplinary Studies, and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';

## Version Infromation
my $version = "    The 8th Step of BSDA (BS-Seq Data Analyzer), version 0.9.4,  2018-02-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$genome_g = 'hg38';                  ## This is only an initialization value or suggesting value, not default value.
$input_g  = '8-bedGraph-Bismark';    ## This is only an initialization value or suggesting value, not default value.
$output_g = '9-Me-BigWig-Bismark';   ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help   -genome   -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  BSDA8.pl  -help' \n";
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
   
my $genomeFolder_g = "/media/yp/biox1/.RefGenomes/Shortcuts2/$genome_g";  
my $chromSize_g    = "$genomeFolder_g/$genome_g.chrom.sizes";
my $BAMpath_g      = "5-finalBAM/2A_Bismark";
&myMakeDir($output_g);

opendir(my $DH_input_g, $BAMpath_g)  ||  die;
my @inputFiles_g = readdir($DH_input_g);
my $pattern_g    = "[-.0-9A-Za-z]+";
my $numCores_g   = 8;
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
        say   "\t......$inputFiles_g[$i]" ;
        my $temp = $inputFiles_g[$i];
        $groupFiles[++$#groupFiles] = $inputFiles_g[$i];
        $temp =~ m/^(\d+)_($pattern_g)_(Rep[1-9])/   or  die   "wrong-1: ## $temp ##";
        $temp =~ m/_(Rep[1-9])\.bam$/   or  die   "wrong-2: ## $temp ##";
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
{
open(seqFiles_FH, ">", "$output_g/BAM-Files.txt")  or  die; 
for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {     
    next unless $inputFiles_g[$i] =~ m/\.bam$/;
    next unless $inputFiles_g[$i] !~ m/^[.]/;
    next unless $inputFiles_g[$i] !~ m/[~]$/;
    next unless $inputFiles_g[$i] !~ m/^unpaired/;
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
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Start bedGraphToBigWig ......";

my $dir1  = "$output_g/1_CHG_OB_bigwig";
my $dir2  = "$output_g/2_CHG_OT_bigwig";
my $dir3  = "$output_g/3_CHH_OB_bigwig";
my $dir4  = "$output_g/4_CHH_OT_bigwig";
my $dir5  = "$output_g/5_CpG_OB_bigwig";
my $dir6  = "$output_g/6_CpG_OT_bigwig";
my $dir7  = "$output_g/7_CH_OB_bigwig";
my $dir8  = "$output_g/8_CH_OT_bigwig";
my $dir9  = "$output_g/9_CHG_bigwig";
my $dir10 = "$output_g/10_CHH_bigwig";
my $dir11 = "$output_g/11_CpG_bigwig";
my $dir12 = "$output_g/12_CH_bigwig";
&myMakeDir($dir1); 
&myMakeDir($dir2); 
&myMakeDir($dir3); 
&myMakeDir($dir4); 
&myMakeDir($dir5); 
&myMakeDir($dir6); 
&myMakeDir($dir7); 
&myMakeDir($dir8); 
&myMakeDir($dir9); 
&myMakeDir($dir10); 
&myMakeDir($dir11); 
&myMakeDir($dir12); 

for (my $i=0; $i<=$#BAMfiles_g; $i++) {
    my $temp = $BAMfiles_g[$i]; 
    $temp =~ s/\.bam$//  ||  die; 
    say   "\t......$temp";

    system("gunzip  --stdout    $input_g/1_CHG_OB/$temp.CHG_OB.bedGraph.gz   >  $input_g/1_CHG_OB/$temp.CHG_OB.bedGraph");
    system("cat  $input_g/1_CHG_OB/$temp.CHG_OB.bedGraph    | sed '1d' |  sort  -k1,1  -k2,2n     >  $input_g/1_CHG_OB/$temp.CHG_OB.sorted.bedGraph");
    system("bedGraphToBigWig  $input_g/1_CHG_OB/$temp.CHG_OB.sorted.bedGraph   $chromSize_g      $dir1/$temp.CHG_OB.bw ");    
    system("rm   $input_g/1_CHG_OB/$temp.CHG_OB.bedGraph");

    system("gunzip  --stdout    $input_g/2_CHG_OT/$temp.CHG_OT.bedGraph.gz   >  $input_g/2_CHG_OT/$temp.CHG_OT.bedGraph");
    system("cat  $input_g/2_CHG_OT/$temp.CHG_OT.bedGraph    | sed '1d' |  sort  -k1,1  -k2,2n     >  $input_g/2_CHG_OT/$temp.CHG_OT.sorted.bedGraph");
    system("bedGraphToBigWig  $input_g/2_CHG_OT/$temp.CHG_OT.sorted.bedGraph   $chromSize_g      $dir2/$temp.CHG_OT.bw ");    
    system("rm   $input_g/2_CHG_OT/$temp.CHG_OT.bedGraph");

    system("gunzip  --stdout    $input_g/3_CHH_OB/$temp.CHH_OB.bedGraph.gz   >  $input_g/3_CHH_OB/$temp.CHH_OB.bedGraph");
    system("cat  $input_g/3_CHH_OB/$temp.CHH_OB.bedGraph    | sed '1d' |  sort  -k1,1  -k2,2n     >  $input_g/3_CHH_OB/$temp.CHH_OB.sorted.bedGraph");
    system("bedGraphToBigWig  $input_g/3_CHH_OB/$temp.CHH_OB.sorted.bedGraph   $chromSize_g      $dir3/$temp.CHH_OB.bw ");    
    system("rm   $input_g/3_CHH_OB/$temp.CHH_OB.bedGraph");

    system("gunzip  --stdout    $input_g/4_CHH_OT/$temp.CHH_OT.bedGraph.gz   >  $input_g/4_CHH_OT/$temp.CHH_OT.bedGraph");
    system("cat  $input_g/4_CHH_OT/$temp.CHH_OT.bedGraph    | sed '1d' |  sort  -k1,1  -k2,2n     >  $input_g/4_CHH_OT/$temp.CHH_OT.sorted.bedGraph");
    system("bedGraphToBigWig  $input_g/4_CHH_OT/$temp.CHH_OT.sorted.bedGraph   $chromSize_g      $dir4/$temp.CHH_OT.bw ");    
    system("rm   $input_g/4_CHH_OT/$temp.CHH_OT.bedGraph");

    system("gunzip  --stdout    $input_g/5_CpG_OB/$temp.CpG_OB.bedGraph.gz   >  $input_g/5_CpG_OB/$temp.CpG_OB.bedGraph");
    system("cat  $input_g/5_CpG_OB/$temp.CpG_OB.bedGraph    | sed '1d' |  sort  -k1,1  -k2,2n     >  $input_g/5_CpG_OB/$temp.CpG_OB.sorted.bedGraph");
    system("bedGraphToBigWig  $input_g/5_CpG_OB/$temp.CpG_OB.sorted.bedGraph   $chromSize_g      $dir5/$temp.CpG_OB.bw ");    
    system("rm   $input_g/5_CpG_OB/$temp.CpG_OB.bedGraph");

    system("gunzip  --stdout    $input_g/6_CpG_OT/$temp.CpG_OT.bedGraph.gz   >  $input_g/6_CpG_OT/$temp.CpG_OT.bedGraph");
    system("cat  $input_g/6_CpG_OT/$temp.CpG_OT.bedGraph    | sed '1d' |  sort  -k1,1  -k2,2n     >  $input_g/6_CpG_OT/$temp.CpG_OT.sorted.bedGraph");
    system("bedGraphToBigWig  $input_g/6_CpG_OT/$temp.CpG_OT.sorted.bedGraph   $chromSize_g      $dir6/$temp.CpG_OT.bw ");    
    system("rm   $input_g/6_CpG_OT/$temp.CpG_OT.bedGraph");

    system("gunzip  --stdout    $input_g/7_CH_OB/$temp.CH_OB.bedGraph.gz   >  $input_g/7_CH_OB/$temp.CH_OB.bedGraph");
    system("cat  $input_g/7_CH_OB/$temp.CH_OB.bedGraph    | sed '1d' |  sort  -k1,1  -k2,2n     >  $input_g/7_CH_OB/$temp.CH_OB.sorted.bedGraph");
    system("bedGraphToBigWig  $input_g/7_CH_OB/$temp.CH_OB.sorted.bedGraph   $chromSize_g      $dir7/$temp.CH_OB.bw ");    
    system("rm   $input_g/7_CH_OB/$temp.CH_OB.bedGraph");

    system("gunzip  --stdout    $input_g/8_CH_OT/$temp.CH_OT.bedGraph.gz   >  $input_g/8_CH_OT/$temp.CH_OT.bedGraph");
    system("cat  $input_g/8_CH_OT/$temp.CH_OT.bedGraph    | sed '1d' |  sort  -k1,1  -k2,2n     >  $input_g/8_CH_OT/$temp.CH_OT.sorted.bedGraph");
    system("bedGraphToBigWig  $input_g/8_CH_OT/$temp.CH_OT.sorted.bedGraph   $chromSize_g      $dir8/$temp.CH_OT.bw ");    
    system("rm   $input_g/8_CH_OT/$temp.CH_OT.bedGraph");

    system("gunzip  --stdout    $input_g/9_CHG/$temp.CHG.bedGraph.gz   >  $input_g/9_CHG/$temp.CHG.bedGraph");
    system("cat  $input_g/9_CHG/$temp.CHG.bedGraph    | sed '1d' |  sort  -k1,1  -k2,2n     >  $input_g/9_CHG/$temp.CHG.sorted.bedGraph");
    system("bedGraphToBigWig  $input_g/9_CHG/$temp.CHG.sorted.bedGraph   $chromSize_g      $dir9/$temp.CHG.bw ");    
    system("rm   $input_g/9_CHG/$temp.CHG.bedGraph");

    system("gunzip  --stdout    $input_g/10_CHH/$temp.CHH.bedGraph.gz   >  $input_g/10_CHH/$temp.CHH.bedGraph");
    system("cat  $input_g/10_CHH/$temp.CHH.bedGraph    | sed '1d' |  sort  -k1,1  -k2,2n     >  $input_g/10_CHH/$temp.CHH.sorted.bedGraph");
    system("bedGraphToBigWig  $input_g/10_CHH/$temp.CHH.sorted.bedGraph   $chromSize_g      $dir10/$temp.CHH.bw ");    
    system("rm   $input_g/10_CHH/$temp.CHH.bedGraph");

    system("gunzip  --stdout    $input_g/11_CpG/$temp.CpG.bedGraph.gz   >  $input_g/11_CpG/$temp.CpG.bedGraph");
    system("cat  $input_g/11_CpG/$temp.CpG.bedGraph    | sed '1d' |  sort  -k1,1  -k2,2n     >  $input_g/11_CpG/$temp.CpG.sorted.bedGraph");
    system("bedGraphToBigWig  $input_g/11_CpG/$temp.CpG.sorted.bedGraph   $chromSize_g      $dir11/$temp.CpG.bw ");    
    system("rm   $input_g/11_CpG/$temp.CpG.bedGraph");

    system("gunzip  --stdout    $input_g/12_CH/$temp.CH.bedGraph.gz   >  $input_g/12_CH/$temp.CH.bedGraph");
    system("cat  $input_g/12_CH/$temp.CH.bedGraph    | sed '1d' |  sort  -k1,1  -k2,2n     >  $input_g/12_CH/$temp.CH.sorted.bedGraph");
    system("bedGraphToBigWig  $input_g/12_CH/$temp.CH.sorted.bedGraph   $chromSize_g      $dir12/$temp.CH.bw ");    
    system("rm   $input_g/12_CH/$temp.CH.bedGraph");
             
}
###################################################################################################################################################################################################


 


###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Start coverage2cytosine ......";

my $dir13 = "$output_g/13_CX_MeReport";   
&myMakeDir($dir13);   

for (my $i=0; $i<=$#BAMfiles_g; $i++) {
    my $temp = $BAMfiles_g[$i]; 
    $temp =~ s/\.bam$//  ||  die; 
    say   "\t......$temp";
    system("coverage2cytosine   --gzip   --genome_folder $genomeFolder_g  --CX_context    --output  $dir13/$temp    $input_g/13_CX/$temp.CX.bismark.cov.gz   >  $dir13/$temp.runLog  2>&1 ");                                                    
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";



  

## END
