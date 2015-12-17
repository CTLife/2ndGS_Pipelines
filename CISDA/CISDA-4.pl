#!/usr/bin/env perl5
use  strict;
use  warnings;
use  v5.20;





###################################################################################################################################################################################################
###################################################################################################################################################################################################
########## Help Infromation ##########
my $HELP_g = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use CISDA (ChIP-Seq Data Analyzer), version 0.62, 2016-01-13.      
        CISDA is a Pipeline for Single-end and Paired-end ChIP-Seq Data Analysis by Integrating Lots of Softwares.

        Step 4: Remove unmapped reads, and mapped reads on chrUn, chrRandom and chrM, and sort the filtered reads.  
                Quality statistics by using FastQC, BamUtil, SAMtools, QualiMap and samstat.
        Usage:  
               perl  CISDA-4.pl    [-v]    [-h]    [-i inputDir]    [-o outDir]   
        For instance: 
                     perl  CISDA-4.pl    -i 4-Mapping          -o 5-SortMapped               
                     perl  CISDA-4.pl    --input 4-Mapping     --output 5-SortMapped    
                     perl  CISDA-4.pl    --input 4-Mapping     --output 5-SortMapped      >> CISDA-4.runLog  2>&1
     
        -------------------------------------------------------------------------------------------------------------------
        Optional arguments:
        -v, --version        Show version number of this program and exit.

        -h, --help           Show this help message and exit.

        Required arguments:
        -i inputDir,  --input inputDir        inputDir is the name of your input folder that contains your SAM files,
                                              the suffix of the SAM files must be ".sam".    (no default)

        -o outDir,  --output outDir           outDir is the name of your output folder that contains running 
                                              results (BAM format) of this step.      (no default)
        ------------------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines such as RASDA, MESDA and HISDA,
        please see the web site: https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies 
        and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.    
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------  
';


########## Version Infromation ##########
my $version_g = "  The Fourth Step of CISDA (ChIP-Seq Data Analyzer), version 0.62, 2016-01-13.";


########## Keys and Values ##########
if ($#ARGV   == -1) { print  "\n$HELP_g\n\n";  exit 0; }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0) { @ARGV = (@ARGV, "-h");           }       ## when the number of command argumants is odd. 
my %args = @ARGV;


########## Initialize  Variables ##########
my $input_g  = '5-SortMapped';      ## This is only an initialization  value or suggesting value, not default value.
my $output_g = '5-SortMapped';      ## This is only an initialization  value or suggesting value, not default value.


########## Available Arguments ##########
my $available = "  -v  --version    -h  --help    -i  --input    -o    --output   ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if($available !~ m/\s$key\s/) {print  "    Cann't recognize $key !!\n";  $boole_g = 1; }
}
if($boole_g == 1) {
    print "\n    The Command Line Arguments are wrong!\n";
    print   '    Please see help message by using "perl  CISDA-4.pl  -h". ';
    print "\n\n";
    exit 0;
}


########## Get Arguments ##########
if ( ( exists $args{'-v' } )  or  ( exists $args{'--version'      } )  )     { print  "\n$version_g\n\n";    exit 0; }
if ( ( exists $args{'-h' } )  or  ( exists $args{'--help'         } )  )     { print  "\n$HELP_g\n\n";       exit 0; }
if ( ( exists $args{'-i' } )  or  ( exists $args{'--input'        } )  )     { ($input_g  = $args{'-i' })  or  ($input_g  = $args{'--input'      });  }else{print   "\n -i or --input  is required.\n\n";   print  "\n$HELP_g\n\n";       exit 0; }                                               
if ( ( exists $args{'-o' } )  or  ( exists $args{'--output'       } )  )     { ($output_g = $args{'-o' })  or  ($output_g = $args{'--output'     });  }else{print   "\n -o or --output is required.\n\n";   print  "\n$HELP_g\n\n";       exit 0; }      


########### Conditions #############
$input_g  =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";
$output_g =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";


######### Print Command Arguments to Standard Output ###########
print  "\n\n
        ################ Your Arguments ###############################
                Input  folder:  $input_g
                Output folder:  $output_g
        ###############################################################  
\n\n";


###################################################################################################################################################################################################
###################################################################################################################################################################################################





print "\n\n\n\n\n##################################################################################################\n";
print   "\nRunning......\n";
my $output2_g = "$output_g/Results";
if ( !(-e $output_g) )   { mkdir $output_g    ||  die; }
if ( !(-e $output2_g))   { mkdir $output2_g   ||  die; }
(-e $output_g)   ||  die;
my $pattern = "[-.0-9A-Za-z]+";







print "\n\n\n\n\n##################################################################################################\n";
print("\nChecking all the necessary softwares in this step......\n");

system("samtools                  >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'    >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");

system("bam                       >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'    >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");

system("fastqc    -v              >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'    >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");

system("qualimap  bamqc           >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'    >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");

system("samstat                   >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'    >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");











##############################################################################################################
sub myFilterSAM  
##############################################################################################################
{
my $name1=$_[0];  ## input dir
my $name2=$_[1];  ## name of input file
my $name3=$_[2];  ## output dir

my $name4="All_$name2";      ## all reads
my $name5="$name2";          ## kept reads
my $name6="Removed_$name2";  ## removed reads

my $n4 = 0; ## all reads
my $n5 = 0; ## kept reads
my $n6 = 0; ## removed reads
  
my $num_chr1  = 0; ## number of reads on chr1
my $num_chr2  = 0; ## number of reads on chr2
my $num_chr3  = 0; ## number of reads on chr3
my $num_chr4  = 0; ## number of reads on chr4
my $num_chr5  = 0; ## number of reads on chr5
my $num_chr6  = 0; ## number of reads on chr6
my $num_chr7  = 0; ## number of reads on chr7
my $num_chr8  = 0; ## number of reads on chr8
my $num_chr9  = 0; ## number of reads on chr9
my $num_chr10 = 0; ## number of reads on chr10
my $num_chr11 = 0; ## number of reads on chr11
my $num_chr12 = 0; ## number of reads on chr12
my $num_chr13 = 0; ## number of reads on chr13
my $num_chr14 = 0; ## number of reads on chr14
my $num_chr15 = 0; ## number of reads on chr15
my $num_chr16 = 0; ## number of reads on chr16
my $num_chr17 = 0; ## number of reads on chr17
my $num_chr18 = 0; ## number of reads on chr18
my $num_chr19 = 0; ## number of reads on chr19
my $num_chrX  = 0; ## number of reads on chrX
my $num_chrY  = 0; ## number of reads on chrY


open(FILE1, "<", "$name1/$name2")  or die "$!";                    
open(FILE4, ">", "$name3/$name4")  or die "$!";  
open(FILE5, ">", "$name3/$name5")  or die "$!";  
open(FILE6, ">", "$name3/$name6")  or die "$!";  
open(FILE7, ">", "$name3/Results/z-$name2.numberOfReads")  or die "$!";

my  $numOfHeader = 0;
while (my $line1=<FILE1>) {
    if ($line1 =~ m/^@/) {
        print  FILE4  $line1   ;  
        print  FILE5  $line1   ;  
        print  FILE6  $line1   ;  
        $numOfHeader++;  
    }else{
        $line1 =~ m/^(\S+)\s+(\S+)\s+(\S+)\s+/  or die;
        my $chr = $3;
        $n4++;
        print  FILE4   $line1;
        given($chr) {
             when("chr1")  {print  FILE5   $line1; $num_chr1++;  $n5++; }       
             when("chr2")  {print  FILE5   $line1; $num_chr2++;  $n5++; }       
             when("chr3")  {print  FILE5   $line1; $num_chr3++;  $n5++; }       
             when("chr4")  {print  FILE5   $line1; $num_chr4++;  $n5++; }       
             when("chr5")  {print  FILE5   $line1; $num_chr5++;  $n5++; }       
             when("chr6")  {print  FILE5   $line1; $num_chr6++;  $n5++; }       
             when("chr7")  {print  FILE5   $line1; $num_chr7++;  $n5++; }       
             when("chr8")  {print  FILE5   $line1; $num_chr8++;  $n5++; }       
             when("chr9")  {print  FILE5   $line1; $num_chr9++;  $n5++; }       
             when("chr10") {print  FILE5   $line1; $num_chr10++; $n5++; }       
             when("chr11") {print  FILE5   $line1; $num_chr11++; $n5++; }       
             when("chr12") {print  FILE5   $line1; $num_chr12++; $n5++; }       
             when("chr13") {print  FILE5   $line1; $num_chr13++; $n5++; }       
             when("chr14") {print  FILE5   $line1; $num_chr14++; $n5++; }       
             when("chr15") {print  FILE5   $line1; $num_chr15++; $n5++; }       
             when("chr16") {print  FILE5   $line1; $num_chr16++; $n5++; }       
             when("chr17") {print  FILE5   $line1; $num_chr17++; $n5++; }       
             when("chr18") {print  FILE5   $line1; $num_chr18++; $n5++; }       
             when("chr19") {print  FILE5   $line1; $num_chr19++; $n5++; }       
             when("chrX")  {print  FILE5   $line1; $num_chrX++;  $n5++; }       
             when("chrY")  {print  FILE5   $line1; $num_chrY++;  $n5++; }   
             default       {print  FILE6   $line1; $n6++;               }    
        }
    }
}
print("numOfHeader:$numOfHeader\n\n");


print  FILE7  "All     Reads: $n4\n";
print  FILE7  "Kept    Reads: $n5\n";
print  FILE7  "Removed Reads: $n6\n\n\n";

print  FILE7  "chr1:  $num_chr1\n";
print  FILE7  "chr2:  $num_chr2\n";
print  FILE7  "chr3:  $num_chr3\n";
print  FILE7  "chr4:  $num_chr4\n";
print  FILE7  "chr5:  $num_chr5\n";
print  FILE7  "chr6:  $num_chr6\n";
print  FILE7  "chr7:  $num_chr7\n";
print  FILE7  "chr8:  $num_chr8\n";
print  FILE7  "chr9:  $num_chr9\n";
print  FILE7  "chr10: $num_chr10\n";
print  FILE7  "chr11: $num_chr11\n";
print  FILE7  "chr12: $num_chr12\n";
print  FILE7  "chr13: $num_chr13\n";
print  FILE7  "chr14: $num_chr14\n";
print  FILE7  "chr15: $num_chr15\n";
print  FILE7  "chr16: $num_chr16\n";
print  FILE7  "chr17: $num_chr17\n";
print  FILE7  "chr18: $num_chr18\n";
print  FILE7  "chr19: $num_chr19\n";
print  FILE7  "chrX:  $num_chrX\n";
print  FILE7  "chrY:  $num_chrY\n";

close FILE1;
close FILE4;
close FILE5;
close FILE6;
close FILE7;

}########################################### END 








my $input_subread  = "$input_g/1-Subread";
my $input_BWAaln   = "$input_g/2-BWAaln";
my $input_BWAmem   = "$input_g/3-BWA-mem";
my $input_Bowtie1  = "$input_g/4-Bowtie1";
my $input_Bowtie2  = "$input_g/5-Bowtie2";

my $output_subread = "$output_g/1-Subread";
my $output_BWAaln  = "$output_g/2-BWAaln";
my $output_BWAmem  = "$output_g/3-BWA-mem";
my $output_Bowtie1 = "$output_g/4-Bowtie1";
my $output_Bowtie2 = "$output_g/5-Bowtie2";

my $output2_subread = "$output_g/1-Subread/Results";
my $output2_BWAaln  = "$output_g/2-BWAaln/Results";
my $output2_BWAmem  = "$output_g/3-BWA-mem/Results";
my $output2_Bowtie1 = "$output_g/4-Bowtie1/Results";
my $output2_Bowtie2 = "$output_g/5-Bowtie2/Results";

if ( !(-e $output_subread) )   { mkdir $output_subread    ||  die; }
if ( !(-e $output_BWAaln ) )   { mkdir $output_BWAaln     ||  die; }
if ( !(-e $output_BWAmem ) )   { mkdir $output_BWAmem     ||  die; }
if ( !(-e $output_Bowtie1) )   { mkdir $output_Bowtie1    ||  die; }
if ( !(-e $output_Bowtie2) )   { mkdir $output_Bowtie2    ||  die; }

if ( !(-e $output2_subread) )   { mkdir $output2_subread    ||  die; }
if ( !(-e $output2_BWAaln ) )   { mkdir $output2_BWAaln     ||  die; }
if ( !(-e $output2_BWAmem ) )   { mkdir $output2_BWAmem     ||  die; }
if ( !(-e $output2_Bowtie1) )   { mkdir $output2_Bowtie1    ||  die; }
if ( !(-e $output2_Bowtie2) )   { mkdir $output2_Bowtie2    ||  die; }








{ ########## Start subread
print "\n\n\n\n\n##################################################################################################\n";
print("\nChecking all the input file names......\n");
opendir(my $DH_input, "$input_subread")  ||  die;     
my @inputFiles = readdir($DH_input);

my $fileNameBool = 1;
for ( my $i=0; $i<=$#inputFiles; $i++ ) {  
        next unless $inputFiles[$i] =~ m/\.sam$/;   
        next unless $inputFiles[$i] !~ m/^[.]/;
        next unless $inputFiles[$i] !~ m/[~]$/;
        next unless $inputFiles[$i] !~ m/^unpaired/;
        my $temp = $inputFiles[$i]; 
        $temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]).sam$/   or  die  "wrong-1: ## $temp ##";
        if($temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])\.sam$/) {
             print("$inputFiles[$i]......\n");
        }else{
             $fileNameBool = 0;
        }
}
if($fileNameBool == 1)  {print("    All the file names are passed.\n\n");}
 
my $qualimapDir = "$output2_subread/qualimap"; 
my $FastQCdir   = "$output2_subread/FastQC";  
my $outDirQC    = "$output2_subread/QCstatistics";

if ( !( -e $qualimapDir ) )   { mkdir $qualimapDir     ||  die; }
if ( !( -e $FastQCdir)    )   { mkdir $FastQCdir       ||  die; }
if ( !( -e $outDirQC)     )   { mkdir $outDirQC        ||  die; }

print "\n\n\n\n\n##################################################################################################\n";
print("\nAnalysis all the sam files ......\n");
for(my $i=0; $i<=$#inputFiles; $i++) {
    next unless $inputFiles[$i] =~ m/\.sam$/;
    next unless $inputFiles[$i] !~ m/^[.]/;
    next unless $inputFiles[$i] !~ m/[~]$/;
    print("$inputFiles[$i] ......\n");      
    my $temp = $inputFiles[$i];
    $temp =~ s/\.sam$//  or  die;

    &myFilterSAM("$input_subread",  "$temp.sam",  "$output_subread");   
    (-e "$output_subread/All_$temp.sam")       ||  die;
    (-e "$output_subread/$temp.sam")           ||  die;
    (-e "$output_subread/Removed_$temp.sam")   ||  die;
    system("samtools  sort   -O bam       -o $output_subread/All_$temp.bam        -T $output_subread/1_$temp      $output_subread/All_$temp.sam         >>$output2_subread/$temp.runLog    2>&1");
    system("samtools  sort   -O bam       -o $output_subread/$temp.bam            -T $output_subread/2_$temp      $output_subread/$temp.sam             >>$output2_subread/$temp.runLog    2>&1");
    system("samtools  sort   -O bam       -o $output_subread/Removed_$temp.bam    -T $output_subread/3_$temp      $output_subread/Removed_$temp.sam     >>$output2_subread/$temp.runLog    2>&1");

    if ( !( -e "$qualimapDir/All_$temp")     )   { mkdir  "$qualimapDir/All_$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_subread/All_$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/All_$temp      >>$output2_subread/$temp.runLog    2>&1");
    system("samstat   $output_subread/All_$temp.bam       >> $output2_subread/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                        >> $output2_subread/$temp.runLog    2>&1");

    if ( !( -e "$qualimapDir/$temp")     )   { mkdir  "$qualimapDir/$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_subread/$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/$temp      >>$output2_subread/$temp.runLog    2>&1");
    system("samstat   $output_subread/$temp.bam       >> $output2_subread/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                        >> $output2_subread/$temp.runLog    2>&1");   

    if ( !( -e "$qualimapDir/Removed_$temp")     )   { mkdir  "$qualimapDir/Removed_$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_subread/Removed_$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/Removed_$temp      >>$output2_subread/$temp.runLog    2>&1");
    system("samstat   $output_subread/Removed_$temp.bam       >> $output2_subread/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                            >> $output2_subread/$temp.runLog    2>&1");

    system("samtools  index       $output_subread/All_$temp.bam      >>$outDirQC/All_$temp.runLog 2>&1");
    system("samtools  flagstat    $output_subread/All_$temp.bam      >>$outDirQC/All_$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_subread/All_$temp.bam      >>$outDirQC/All_$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_subread/All_$temp.bam      >>$outDirQC/All_$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_subread/All_$temp.bam  --basic  --qual  --phred    >>$outDirQC/All_$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_subread/All_$temp.bam          >>$FastQCdir/All_$temp.runLog        2>&1");  

    system("samtools  index       $output_subread/$temp.bam      >>$outDirQC/$temp.runLog 2>&1");
    system("samtools  flagstat    $output_subread/$temp.bam      >>$outDirQC/$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_subread/$temp.bam      >>$outDirQC/$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_subread/$temp.bam      >>$outDirQC/$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_subread/$temp.bam  --basic  --qual  --phred    >>$outDirQC/$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_subread/$temp.bam          >>$FastQCdir/$temp.runLog        2>&1");  

    system("samtools  index       $output_subread/Removed_$temp.bam      >>$outDirQC/Removed_$temp.runLog 2>&1");
    system("samtools  flagstat    $output_subread/Removed_$temp.bam      >>$outDirQC/Removed_$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_subread/Removed_$temp.bam      >>$outDirQC/Removed_$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_subread/Removed_$temp.bam      >>$outDirQC/Removed_$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_subread/Removed_$temp.bam  --basic  --qual  --phred    >>$outDirQC/Removed_$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_subread/Removed_$temp.bam          >>$FastQCdir/Removed_$temp.runLog        2>&1");  

    system("rm   $output_subread/All_$temp.sam");
    system("rm   $output_subread/$temp.sam");
    system("rm   $output_subread/Removed_$temp.sam");
}
} ########## End subread










{ ########## Start BWAaln
print "\n\n\n\n\n##################################################################################################\n";
print("\nChecking all the input file names......\n");
opendir(my $DH_input, "$input_BWAaln")  ||  die;     
my @inputFiles = readdir($DH_input);

my $fileNameBool = 1;
for ( my $i=0; $i<=$#inputFiles; $i++ ) {  
        next unless $inputFiles[$i] =~ m/\.sam$/;   
        next unless $inputFiles[$i] !~ m/^[.]/;
        next unless $inputFiles[$i] !~ m/[~]$/;
        next unless $inputFiles[$i] !~ m/^unpaired/;
        my $temp = $inputFiles[$i]; 
        $temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]).sam$/   or  die  "wrong-1: ## $temp ##";
        if($temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])\.sam$/) {
             print("$inputFiles[$i]......\n");
        }else{
             $fileNameBool = 0;
        }
}
if($fileNameBool == 1)  {print("    All the file names are passed.\n\n");}

my $qualimapDir = "$output2_BWAaln/qualimap"; 
my $FastQCdir   = "$output2_BWAaln/FastQC";  
my $outDirQC    = "$output2_BWAaln/QCstatistics";

if ( !( -e $qualimapDir ) )   { mkdir $qualimapDir     ||  die; }
if ( !( -e $FastQCdir)    )   { mkdir $FastQCdir       ||  die; }
if ( !( -e $outDirQC)     )   { mkdir $outDirQC        ||  die; }

print "\n\n\n\n\n##################################################################################################\n";
print("\nAnalysis all the sam files ......\n");
for(my $i=0; $i<=$#inputFiles; $i++) {
    next unless $inputFiles[$i] =~ m/\.sam$/;
    next unless $inputFiles[$i] !~ m/^[.]/;
    next unless $inputFiles[$i] !~ m/[~]$/;
    print("$inputFiles[$i] ......\n");      
    my $temp = $inputFiles[$i];
    $temp =~ s/\.sam$//  or  die;

    &myFilterSAM("$input_BWAaln",  "$temp.sam",  "$output_BWAaln");   
    (-e "$output_BWAaln/All_$temp.sam")       ||  die;
    (-e "$output_BWAaln/$temp.sam")           ||  die;
    (-e "$output_BWAaln/Removed_$temp.sam")   ||  die;
    system("samtools  sort   -O bam       -o $output_BWAaln/All_$temp.bam        -T $output_BWAaln/1_$temp      $output_BWAaln/All_$temp.sam         >>$output2_BWAaln/$temp.runLog    2>&1");
    system("samtools  sort   -O bam       -o $output_BWAaln/$temp.bam            -T $output_BWAaln/2_$temp      $output_BWAaln/$temp.sam             >>$output2_BWAaln/$temp.runLog    2>&1");
    system("samtools  sort   -O bam       -o $output_BWAaln/Removed_$temp.bam    -T $output_BWAaln/3_$temp      $output_BWAaln/Removed_$temp.sam     >>$output2_BWAaln/$temp.runLog    2>&1");

    if ( !( -e "$qualimapDir/All_$temp")     )   { mkdir  "$qualimapDir/All_$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_BWAaln/All_$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/All_$temp      >>$output2_BWAaln/$temp.runLog    2>&1");
    system("samstat   $output_BWAaln/All_$temp.bam       >> $output2_BWAaln/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                        >> $output2_BWAaln/$temp.runLog    2>&1");

    if ( !( -e "$qualimapDir/$temp")     )   { mkdir  "$qualimapDir/$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_BWAaln/$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/$temp      >>$output2_BWAaln/$temp.runLog    2>&1");
    system("samstat   $output_BWAaln/$temp.bam       >> $output2_BWAaln/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                        >> $output2_BWAaln/$temp.runLog    2>&1");   

    if ( !( -e "$qualimapDir/Removed_$temp")     )   { mkdir  "$qualimapDir/Removed_$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_BWAaln/Removed_$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/Removed_$temp      >>$output2_BWAaln/$temp.runLog    2>&1");
    system("samstat   $output_BWAaln/Removed_$temp.bam       >> $output2_BWAaln/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                            >> $output2_BWAaln/$temp.runLog    2>&1");

    system("samtools  index       $output_BWAaln/All_$temp.bam      >>$outDirQC/All_$temp.runLog 2>&1");
    system("samtools  flagstat    $output_BWAaln/All_$temp.bam      >>$outDirQC/All_$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_BWAaln/All_$temp.bam      >>$outDirQC/All_$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_BWAaln/All_$temp.bam      >>$outDirQC/All_$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_BWAaln/All_$temp.bam  --basic  --qual  --phred    >>$outDirQC/All_$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_BWAaln/All_$temp.bam          >>$FastQCdir/All_$temp.runLog        2>&1");  

    system("samtools  index       $output_BWAaln/$temp.bam      >>$outDirQC/$temp.runLog 2>&1");
    system("samtools  flagstat    $output_BWAaln/$temp.bam      >>$outDirQC/$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_BWAaln/$temp.bam      >>$outDirQC/$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_BWAaln/$temp.bam      >>$outDirQC/$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_BWAaln/$temp.bam  --basic  --qual  --phred    >>$outDirQC/$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_BWAaln/$temp.bam          >>$FastQCdir/$temp.runLog        2>&1");  

    system("samtools  index       $output_BWAaln/Removed_$temp.bam      >>$outDirQC/Removed_$temp.runLog 2>&1");
    system("samtools  flagstat    $output_BWAaln/Removed_$temp.bam      >>$outDirQC/Removed_$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_BWAaln/Removed_$temp.bam      >>$outDirQC/Removed_$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_BWAaln/Removed_$temp.bam      >>$outDirQC/Removed_$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_BWAaln/Removed_$temp.bam  --basic  --qual  --phred    >>$outDirQC/Removed_$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_BWAaln/Removed_$temp.bam          >>$FastQCdir/Removed_$temp.runLog        2>&1");  

    system("rm   $output_BWAaln/All_$temp.sam");
    system("rm   $output_BWAaln/$temp.sam");
    system("rm   $output_BWAaln/Removed_$temp.sam");
}
} ########## End BWAaln










{ ########## Start BWAmem
print "\n\n\n\n\n##################################################################################################\n";
print("\nChecking all the input file names......\n");
opendir(my $DH_input, "$input_BWAmem")  ||  die;     
my @inputFiles = readdir($DH_input);

my $fileNameBool = 1;
for ( my $i=0; $i<=$#inputFiles; $i++ ) {  
        next unless $inputFiles[$i] =~ m/\.sam$/;   
        next unless $inputFiles[$i] !~ m/^[.]/;
        next unless $inputFiles[$i] !~ m/[~]$/;
        next unless $inputFiles[$i] !~ m/^unpaired/;
        my $temp = $inputFiles[$i]; 
        $temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]).sam$/   or  die  "wrong-1: ## $temp ##";
        if($temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])\.sam$/) {
             print("$inputFiles[$i]......\n");
        }else{
             $fileNameBool = 0;
        }
}
if($fileNameBool == 1)  {print("    All the file names are passed.\n\n");}

my $qualimapDir = "$output2_BWAmem/qualimap"; 
my $FastQCdir   = "$output2_BWAmem/FastQC";  
my $outDirQC    = "$output2_BWAmem/QCstatistics";

if ( !( -e $qualimapDir ) )   { mkdir $qualimapDir     ||  die; }
if ( !( -e $FastQCdir)    )   { mkdir $FastQCdir       ||  die; }
if ( !( -e $outDirQC)     )   { mkdir $outDirQC        ||  die; }

print "\n\n\n\n\n##################################################################################################\n";
print("\nAnalysis all the sam files ......\n");
for(my $i=0; $i<=$#inputFiles; $i++) {
    next unless $inputFiles[$i] =~ m/\.sam$/;
    next unless $inputFiles[$i] !~ m/^[.]/;
    next unless $inputFiles[$i] !~ m/[~]$/;
    print("$inputFiles[$i] ......\n");      
    my $temp = $inputFiles[$i];
    $temp =~ s/\.sam$//  or  die;

    &myFilterSAM("$input_BWAmem",  "$temp.sam",  "$output_BWAmem");   
    (-e "$output_BWAmem/All_$temp.sam")       ||  die;
    (-e "$output_BWAmem/$temp.sam")           ||  die;
    (-e "$output_BWAmem/Removed_$temp.sam")   ||  die;
    system("samtools  sort   -O bam       -o $output_BWAmem/All_$temp.bam        -T $output_BWAmem/1_$temp      $output_BWAmem/All_$temp.sam         >>$output2_BWAmem/$temp.runLog    2>&1");
    system("samtools  sort   -O bam       -o $output_BWAmem/$temp.bam            -T $output_BWAmem/2_$temp      $output_BWAmem/$temp.sam             >>$output2_BWAmem/$temp.runLog    2>&1");
    system("samtools  sort   -O bam       -o $output_BWAmem/Removed_$temp.bam    -T $output_BWAmem/3_$temp      $output_BWAmem/Removed_$temp.sam     >>$output2_BWAmem/$temp.runLog    2>&1");

    if ( !( -e "$qualimapDir/All_$temp")     )   { mkdir  "$qualimapDir/All_$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_BWAmem/All_$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/All_$temp      >>$output2_BWAmem/$temp.runLog    2>&1");
    system("samstat   $output_BWAmem/All_$temp.bam       >> $output2_BWAmem/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                        >> $output2_BWAmem/$temp.runLog    2>&1");

    if ( !( -e "$qualimapDir/$temp")     )   { mkdir  "$qualimapDir/$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_BWAmem/$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/$temp      >>$output2_BWAmem/$temp.runLog    2>&1");
    system("samstat   $output_BWAmem/$temp.bam       >> $output2_BWAmem/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                        >> $output2_BWAmem/$temp.runLog    2>&1");   

    if ( !( -e "$qualimapDir/Removed_$temp")     )   { mkdir  "$qualimapDir/Removed_$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_BWAmem/Removed_$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/Removed_$temp      >>$output2_BWAmem/$temp.runLog    2>&1");
    system("samstat   $output_BWAmem/Removed_$temp.bam       >> $output2_BWAmem/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                            >> $output2_BWAmem/$temp.runLog    2>&1");

    system("samtools  index       $output_BWAmem/All_$temp.bam      >>$outDirQC/All_$temp.runLog 2>&1");
    system("samtools  flagstat    $output_BWAmem/All_$temp.bam      >>$outDirQC/All_$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_BWAmem/All_$temp.bam      >>$outDirQC/All_$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_BWAmem/All_$temp.bam      >>$outDirQC/All_$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_BWAmem/All_$temp.bam  --basic  --qual  --phred    >>$outDirQC/All_$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_BWAmem/All_$temp.bam          >>$FastQCdir/All_$temp.runLog        2>&1");  

    system("samtools  index       $output_BWAmem/$temp.bam      >>$outDirQC/$temp.runLog 2>&1");
    system("samtools  flagstat    $output_BWAmem/$temp.bam      >>$outDirQC/$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_BWAmem/$temp.bam      >>$outDirQC/$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_BWAmem/$temp.bam      >>$outDirQC/$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_BWAmem/$temp.bam  --basic  --qual  --phred    >>$outDirQC/$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_BWAmem/$temp.bam          >>$FastQCdir/$temp.runLog        2>&1");  

    system("samtools  index       $output_BWAmem/Removed_$temp.bam      >>$outDirQC/Removed_$temp.runLog 2>&1");
    system("samtools  flagstat    $output_BWAmem/Removed_$temp.bam      >>$outDirQC/Removed_$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_BWAmem/Removed_$temp.bam      >>$outDirQC/Removed_$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_BWAmem/Removed_$temp.bam      >>$outDirQC/Removed_$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_BWAmem/Removed_$temp.bam  --basic  --qual  --phred    >>$outDirQC/Removed_$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_BWAmem/Removed_$temp.bam          >>$FastQCdir/Removed_$temp.runLog        2>&1");  

    system("rm   $output_BWAmem/All_$temp.sam");
    system("rm   $output_BWAmem/$temp.sam");
    system("rm   $output_BWAmem/Removed_$temp.sam");
}
} ########## End BWAmem










{ ########## Start Bowtie1
print "\n\n\n\n\n##################################################################################################\n";
print("\nChecking all the input file names......\n");
opendir(my $DH_input, "$input_Bowtie1")  ||  die;     
my @inputFiles = readdir($DH_input);

my $fileNameBool = 1;
for ( my $i=0; $i<=$#inputFiles; $i++ ) {  
        next unless $inputFiles[$i] =~ m/\.sam$/;   
        next unless $inputFiles[$i] !~ m/^[.]/;
        next unless $inputFiles[$i] !~ m/[~]$/;
        next unless $inputFiles[$i] !~ m/^unpaired/;
        my $temp = $inputFiles[$i]; 
        $temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]).sam$/   or  die  "wrong-1: ## $temp ##";
        if($temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])\.sam$/) {
             print("$inputFiles[$i]......\n");
        }else{
             $fileNameBool = 0;
        }
}
if($fileNameBool == 1)  {print("    All the file names are passed.\n\n");}

my $qualimapDir = "$output2_Bowtie1/qualimap"; 
my $FastQCdir   = "$output2_Bowtie1/FastQC";  
my $outDirQC    = "$output2_Bowtie1/QCstatistics";

if ( !( -e $qualimapDir ) )   { mkdir $qualimapDir     ||  die; }
if ( !( -e $FastQCdir)    )   { mkdir $FastQCdir       ||  die; }
if ( !( -e $outDirQC)     )   { mkdir $outDirQC        ||  die; }

print "\n\n\n\n\n##################################################################################################\n";
print("\nAnalysis all the sam files ......\n");
for(my $i=0; $i<=$#inputFiles; $i++) {
    next unless $inputFiles[$i] =~ m/\.sam$/;
    next unless $inputFiles[$i] !~ m/^[.]/;
    next unless $inputFiles[$i] !~ m/[~]$/;
    print("$inputFiles[$i] ......\n");      
    my $temp = $inputFiles[$i];
    $temp =~ s/\.sam$//  or  die;

    &myFilterSAM("$input_Bowtie1",  "$temp.sam",  "$output_Bowtie1");   
    (-e "$output_Bowtie1/All_$temp.sam")       ||  die;
    (-e "$output_Bowtie1/$temp.sam")           ||  die;
    (-e "$output_Bowtie1/Removed_$temp.sam")   ||  die;
    system("samtools  sort   -O bam       -o $output_Bowtie1/All_$temp.bam        -T $output_Bowtie1/1_$temp      $output_Bowtie1/All_$temp.sam         >>$output2_Bowtie1/$temp.runLog    2>&1");
    system("samtools  sort   -O bam       -o $output_Bowtie1/$temp.bam            -T $output_Bowtie1/2_$temp      $output_Bowtie1/$temp.sam             >>$output2_Bowtie1/$temp.runLog    2>&1");
    system("samtools  sort   -O bam       -o $output_Bowtie1/Removed_$temp.bam    -T $output_Bowtie1/3_$temp      $output_Bowtie1/Removed_$temp.sam     >>$output2_Bowtie1/$temp.runLog    2>&1");

    if ( !( -e "$qualimapDir/All_$temp")     )   { mkdir  "$qualimapDir/All_$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_Bowtie1/All_$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/All_$temp      >>$output2_Bowtie1/$temp.runLog    2>&1");
    system("samstat   $output_Bowtie1/All_$temp.bam       >> $output2_Bowtie1/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                        >> $output2_Bowtie1/$temp.runLog    2>&1");

    if ( !( -e "$qualimapDir/$temp")     )   { mkdir  "$qualimapDir/$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_Bowtie1/$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/$temp      >>$output2_Bowtie1/$temp.runLog    2>&1");
    system("samstat   $output_Bowtie1/$temp.bam       >> $output2_Bowtie1/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                        >> $output2_Bowtie1/$temp.runLog    2>&1");   

    if ( !( -e "$qualimapDir/Removed_$temp")     )   { mkdir  "$qualimapDir/Removed_$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_Bowtie1/Removed_$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/Removed_$temp      >>$output2_Bowtie1/$temp.runLog    2>&1");
    system("samstat   $output_Bowtie1/Removed_$temp.bam       >> $output2_Bowtie1/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                            >> $output2_Bowtie1/$temp.runLog    2>&1");

    system("samtools  index       $output_Bowtie1/All_$temp.bam      >>$outDirQC/All_$temp.runLog 2>&1");
    system("samtools  flagstat    $output_Bowtie1/All_$temp.bam      >>$outDirQC/All_$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_Bowtie1/All_$temp.bam      >>$outDirQC/All_$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_Bowtie1/All_$temp.bam      >>$outDirQC/All_$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_Bowtie1/All_$temp.bam  --basic  --qual  --phred    >>$outDirQC/All_$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_Bowtie1/All_$temp.bam          >>$FastQCdir/All_$temp.runLog        2>&1");  

    system("samtools  index       $output_Bowtie1/$temp.bam      >>$outDirQC/$temp.runLog 2>&1");
    system("samtools  flagstat    $output_Bowtie1/$temp.bam      >>$outDirQC/$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_Bowtie1/$temp.bam      >>$outDirQC/$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_Bowtie1/$temp.bam      >>$outDirQC/$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_Bowtie1/$temp.bam  --basic  --qual  --phred    >>$outDirQC/$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_Bowtie1/$temp.bam          >>$FastQCdir/$temp.runLog        2>&1");  

    system("samtools  index       $output_Bowtie1/Removed_$temp.bam      >>$outDirQC/Removed_$temp.runLog 2>&1");
    system("samtools  flagstat    $output_Bowtie1/Removed_$temp.bam      >>$outDirQC/Removed_$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_Bowtie1/Removed_$temp.bam      >>$outDirQC/Removed_$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_Bowtie1/Removed_$temp.bam      >>$outDirQC/Removed_$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_Bowtie1/Removed_$temp.bam  --basic  --qual  --phred    >>$outDirQC/Removed_$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_Bowtie1/Removed_$temp.bam          >>$FastQCdir/Removed_$temp.runLog        2>&1");  

    system("rm   $output_Bowtie1/All_$temp.sam");
    system("rm   $output_Bowtie1/$temp.sam");
    system("rm   $output_Bowtie1/Removed_$temp.sam");
}
} ########## End Bowtie1











{ ########## Start Bowtie2
print "\n\n\n\n\n##################################################################################################\n";
print("\nChecking all the input file names......\n");
opendir(my $DH_input, "$input_Bowtie2")  ||  die;     
my @inputFiles = readdir($DH_input);

my $fileNameBool = 1;
for ( my $i=0; $i<=$#inputFiles; $i++ ) {  
        next unless $inputFiles[$i] =~ m/\.sam$/;   
        next unless $inputFiles[$i] !~ m/^[.]/;
        next unless $inputFiles[$i] !~ m/[~]$/;
        next unless $inputFiles[$i] !~ m/^unpaired/;
        my $temp = $inputFiles[$i]; 
        $temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]).sam$/   or  die  "wrong-1: ## $temp ##";
        if($temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])\.sam$/) {
             print("$inputFiles[$i]......\n");
        }else{
             $fileNameBool = 0;
        }
}
if($fileNameBool == 1)  {print("    All the file names are passed.\n\n");}

my $qualimapDir = "$output2_Bowtie2/qualimap"; 
my $FastQCdir   = "$output2_Bowtie2/FastQC";  
my $outDirQC    = "$output2_Bowtie2/QCstatistics";

if ( !( -e $qualimapDir ) )   { mkdir $qualimapDir     ||  die; }
if ( !( -e $FastQCdir)    )   { mkdir $FastQCdir       ||  die; }
if ( !( -e $outDirQC)     )   { mkdir $outDirQC        ||  die; }

print "\n\n\n\n\n##################################################################################################\n";
print("\nAnalysis all the sam files ......\n");
for(my $i=0; $i<=$#inputFiles; $i++) {
    next unless $inputFiles[$i] =~ m/\.sam$/;
    next unless $inputFiles[$i] !~ m/^[.]/;
    next unless $inputFiles[$i] !~ m/[~]$/;
    print("$inputFiles[$i] ......\n");      
    my $temp = $inputFiles[$i];
    $temp =~ s/\.sam$//  or  die;

    &myFilterSAM("$input_Bowtie2",  "$temp.sam",  "$output_Bowtie2");   
    (-e "$output_Bowtie2/All_$temp.sam")       ||  die;
    (-e "$output_Bowtie2/$temp.sam")           ||  die;
    (-e "$output_Bowtie2/Removed_$temp.sam")   ||  die;
    system("samtools  sort   -O bam       -o $output_Bowtie2/All_$temp.bam        -T $output_Bowtie2/1_$temp      $output_Bowtie2/All_$temp.sam         >>$output2_Bowtie2/$temp.runLog    2>&1");
    system("samtools  sort   -O bam       -o $output_Bowtie2/$temp.bam            -T $output_Bowtie2/2_$temp      $output_Bowtie2/$temp.sam             >>$output2_Bowtie2/$temp.runLog    2>&1");
    system("samtools  sort   -O bam       -o $output_Bowtie2/Removed_$temp.bam    -T $output_Bowtie2/3_$temp      $output_Bowtie2/Removed_$temp.sam     >>$output2_Bowtie2/$temp.runLog    2>&1");

    if ( !( -e "$qualimapDir/All_$temp")     )   { mkdir  "$qualimapDir/All_$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_Bowtie2/All_$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/All_$temp      >>$output2_Bowtie2/$temp.runLog    2>&1");
    system("samstat   $output_Bowtie2/All_$temp.bam       >> $output2_Bowtie2/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                        >> $output2_Bowtie2/$temp.runLog    2>&1");

    if ( !( -e "$qualimapDir/$temp")     )   { mkdir  "$qualimapDir/$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_Bowtie2/$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/$temp      >>$output2_Bowtie2/$temp.runLog    2>&1");
    system("samstat   $output_Bowtie2/$temp.bam       >> $output2_Bowtie2/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                        >> $output2_Bowtie2/$temp.runLog    2>&1");   

    if ( !( -e "$qualimapDir/Removed_$temp")     )   { mkdir  "$qualimapDir/Removed_$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_Bowtie2/Removed_$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/Removed_$temp      >>$output2_Bowtie2/$temp.runLog    2>&1");
    system("samstat   $output_Bowtie2/Removed_$temp.bam       >> $output2_Bowtie2/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                            >> $output2_Bowtie2/$temp.runLog    2>&1");

    system("samtools  index       $output_Bowtie2/All_$temp.bam      >>$outDirQC/All_$temp.runLog 2>&1");
    system("samtools  flagstat    $output_Bowtie2/All_$temp.bam      >>$outDirQC/All_$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_Bowtie2/All_$temp.bam      >>$outDirQC/All_$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_Bowtie2/All_$temp.bam      >>$outDirQC/All_$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_Bowtie2/All_$temp.bam  --basic  --qual  --phred    >>$outDirQC/All_$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_Bowtie2/All_$temp.bam          >>$FastQCdir/All_$temp.runLog        2>&1");  

    system("samtools  index       $output_Bowtie2/$temp.bam      >>$outDirQC/$temp.runLog 2>&1");
    system("samtools  flagstat    $output_Bowtie2/$temp.bam      >>$outDirQC/$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_Bowtie2/$temp.bam      >>$outDirQC/$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_Bowtie2/$temp.bam      >>$outDirQC/$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_Bowtie2/$temp.bam  --basic  --qual  --phred    >>$outDirQC/$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_Bowtie2/$temp.bam          >>$FastQCdir/$temp.runLog        2>&1");  

    system("samtools  index       $output_Bowtie2/Removed_$temp.bam      >>$outDirQC/Removed_$temp.runLog 2>&1");
    system("samtools  flagstat    $output_Bowtie2/Removed_$temp.bam      >>$outDirQC/Removed_$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_Bowtie2/Removed_$temp.bam      >>$outDirQC/Removed_$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_Bowtie2/Removed_$temp.bam      >>$outDirQC/Removed_$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_Bowtie2/Removed_$temp.bam  --basic  --qual  --phred    >>$outDirQC/Removed_$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_Bowtie2/Removed_$temp.bam          >>$FastQCdir/Removed_$temp.runLog        2>&1");  

    system("rm   $output_Bowtie2/All_$temp.sam");
    system("rm   $output_Bowtie2/$temp.sam");
    system("rm   $output_Bowtie2/Removed_$temp.sam");
}
} ########## End Bowtie2









print "\n\n\n\n\n##################################################################################################\n";
print "\n\n        Job Done! Cheers! \n\n";





























