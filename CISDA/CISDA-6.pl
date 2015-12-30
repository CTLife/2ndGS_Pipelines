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

        Step 6: Only  retain the mapped reads with MAPQ>=30 and no PCR duplicates.  
                And repair the unpaired reads in paired-end file.
                Quality statistics by using FastQC, BamUtil, SAMtools, QualiMap and samstat.
        Usage:  
               perl  CISDA-6.pl    [-v]    [-h]    [-i inputDir]    [-o outDir]   
        For instance: 
                     perl  CISDA-6.pl    -i 6-MisMatch          -o 7-FinalBAM              
                     perl  CISDA-6.pl    --input 6-MisMatch     --output 7-FinalBAM     
                     perl  CISDA-6.pl    --input 6-MisMatch     --output 7-FinalBAM    >> CISDA-6.runLog  2>&1
     
        -------------------------------------------------------------------------------------------------------------------
        Optional arguments:
        -v, --version        Show version number of this program and exit.

        -h, --help           Show this help message and exit.

        Required arguments:
        -i inputDir,  --input inputDir        inputDir is the name of your input folder that contains your BAM files,
                                              the suffix of the BAM files must be ".bam".    (no default)

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
my $version_g = "  The Sixth Step of CISDA (ChIP-Seq Data Analyzer), version 0.62, 2016-01-13.";


########## Keys and Values ##########
if ($#ARGV   == -1) { print  "\n$HELP_g\n\n";  exit 0; }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0) { @ARGV = (@ARGV, "-h");           }       ## when the number of command argumants is odd. 
my %args = @ARGV;


########## Initialize  Variables ##########
my $input_g  = '7-FinalBAM';      ## This is only an initialization  value or suggesting value, not default value.
my $output_g = '7-FinalBAM';      ## This is only an initialization  value or suggesting value, not default value.


########## Available Arguments ##########
my $available = "  -v  --version    -h  --help    -i  --input    -o    --output     -MM  --MaxMismatch     ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if($available !~ m/\s$key\s/) {print  "    Cann't recognize $key !!\n";  $boole_g = 1; }
}
if($boole_g == 1) {
    print "\n    The Command Line Arguments are wrong!\n";
    print   '    Please see help message by using "perl  CISDA-6.pl  -h". ';
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

system("samstat                   >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'    >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");

system("removeDup                 >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'    >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");

system("repair                    >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'    >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");









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
        next unless $inputFiles[$i] =~ m/^(\d{2})_/;   
        next unless $inputFiles[$i] =~ m/\.bam$/;   
        next unless $inputFiles[$i] !~ m/^[.]/;
        next unless $inputFiles[$i] !~ m/[~]$/;
        next unless $inputFiles[$i] !~ m/^unpaired/;
        my $temp = $inputFiles[$i]; 
        $temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]).bam$/   or  die  "wrong-1: ## $temp ##";
        if($temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])\.bam$/) {
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
print("\nAnalysis all the bam files ......\n");
for(my $i=0; $i<=$#inputFiles; $i++) {
    next unless $inputFiles[$i] =~ m/^(\d{2})_/;   
    next unless $inputFiles[$i] =~ m/\.bam$/;
    next unless $inputFiles[$i] !~ m/^[.]/;
    next unless $inputFiles[$i] !~ m/[~]$/;
    print("$inputFiles[$i] ......\n");      
    my $temp = $inputFiles[$i];
    $temp =~ s/\.bam$//  or  die;

    system(`samtools  view  -h  $input_subread/$temp.bam  | wc  -l     >>$output2_subread/$temp.number    2>&1`);
    system(`samtools  view   -h   -q 30   -o $output_subread/$temp.q30.sam       $input_subread/$temp.bam      >>$output2_subread/$temp.runLog    2>&1`);    ## must use `,  cann't use "
    system(`wc  -l  $output_subread/$temp.q30.sam  >>$output2_subread/$temp.number    2>&1`);  
    if ( !( -e "$output_subread/1_$temp")     )   { mkdir "$output_subread/1_$temp"        ||  die; }
    system("removeDup    -r 1000    -t $output_subread/1_$temp    -c 999999999     -i  $output_subread/$temp.q30.sam   -o $output_subread/$temp.noDup.sam            >>$output2_subread/$temp.runLog    2>&1");  ## SAM file must contain header. 
    system("rm   -rf  $output_subread/1_$temp");
    system(`wc  -l  $output_subread/$temp.noDup.sam  >>$output2_subread/$temp.number    2>&1`);  
    system("repair    -S  -c  -T 16    -i  $output_subread/$temp.noDup.sam   -o $output_subread/$temp.unsort.bam            >>$output2_subread/$temp.runLog    2>&1");  ## SAM file must contain header.         
    system(`samtools  view  -h  $output_subread/$temp.unsort.bam  | wc  -l     >>$output2_subread/$temp.number    2>&1`);
    system("rm   $output_subread/$temp.q30.sam");
    system("rm   $output_subread/$temp.noDup.sam");
    system("samtools  sort   -O bam       -o $output_subread/$temp.bam        -T $output_subread/1_$temp      $output_subread/$temp.unsort.bam         >>$output2_subread/$temp.runLog    2>&1");
    system("rm   $output_subread/$temp.unsort.bam");

    if ( !( -e "$qualimapDir/$temp")     )   { mkdir  "$qualimapDir/$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_subread/$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/$temp    --java-mem-size=5G    >>$output2_subread/$temp.runLog    2>&1");
    system("samstat   $output_subread/$temp.bam           >> $output2_subread/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                        >> $output2_subread/$temp.runLog    2>&1");   

    system("samtools  index       $output_subread/$temp.bam      >>$outDirQC/$temp.runLog  2>&1");
    system("samtools  flagstat    $output_subread/$temp.bam      >>$outDirQC/$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_subread/$temp.bam      >>$outDirQC/$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_subread/$temp.bam      >>$outDirQC/$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_subread/$temp.bam  --basic  --qual  --phred    >>$outDirQC/$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_subread/$temp.bam          >>$FastQCdir/$temp.runLog        2>&1");  

    
}
} ########## End subread

















{ ########## Start BWAaln
print "\n\n\n\n\n##################################################################################################\n";
print("\nChecking all the input file names......\n");
opendir(my $DH_input, "$input_BWAaln")  ||  die;     
my @inputFiles = readdir($DH_input);

my $fileNameBool = 1;
for ( my $i=0; $i<=$#inputFiles; $i++ ) {  
        next unless $inputFiles[$i] =~ m/^(\d{2})_/;   
        next unless $inputFiles[$i] =~ m/\.bam$/;   
        next unless $inputFiles[$i] !~ m/^[.]/;
        next unless $inputFiles[$i] !~ m/[~]$/;
        next unless $inputFiles[$i] !~ m/^unpaired/;
        my $temp = $inputFiles[$i]; 
        $temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]).bam$/   or  die  "wrong-1: ## $temp ##";
        if($temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])\.bam$/) {
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
print("\nAnalysis all the bam files ......\n");
for(my $i=0; $i<=$#inputFiles; $i++) {
    next unless $inputFiles[$i] =~ m/^(\d{2})_/;   
    next unless $inputFiles[$i] =~ m/\.bam$/;
    next unless $inputFiles[$i] !~ m/^[.]/;
    next unless $inputFiles[$i] !~ m/[~]$/;
    print("$inputFiles[$i] ......\n");      
    my $temp = $inputFiles[$i];
    $temp =~ s/\.bam$//  or  die;

    system(`samtools  view  -h  $input_BWAaln/$temp.bam  | wc  -l     >>$output2_BWAaln/$temp.number    2>&1`);
    system(`samtools  view   -h   -q 30   -o $output_BWAaln/$temp.q30.sam       $input_BWAaln/$temp.bam      >>$output2_BWAaln/$temp.runLog    2>&1`);    ## must use `,  cann't use "
    system(`wc  -l  $output_BWAaln/$temp.q30.sam  >>$output2_BWAaln/$temp.number    2>&1`);  
    if ( !( -e "$output_BWAaln/1_$temp")     )   { mkdir "$output_BWAaln/1_$temp"        ||  die; }
    system("removeDup    -r 1000    -t $output_BWAaln/1_$temp    -c 999999999     -i  $output_BWAaln/$temp.q30.sam   -o $output_BWAaln/$temp.noDup.sam            >>$output2_BWAaln/$temp.runLog    2>&1");  ## SAM file must contain header. 
    system("rm   -rf  $output_BWAaln/1_$temp");
    system(`wc  -l  $output_BWAaln/$temp.noDup.sam  >>$output2_BWAaln/$temp.number    2>&1`);  
    system("repair    -S  -c  -T 16    -i  $output_BWAaln/$temp.noDup.sam   -o $output_BWAaln/$temp.unsort.bam            >>$output2_BWAaln/$temp.runLog    2>&1");  ## SAM file must contain header.         
    system(`samtools  view  -h  $output_BWAaln/$temp.unsort.bam  | wc  -l     >>$output2_BWAaln/$temp.number    2>&1`);
    system("rm   $output_BWAaln/$temp.q30.sam");
    system("rm   $output_BWAaln/$temp.noDup.sam");
    system("samtools  sort   -O bam       -o $output_BWAaln/$temp.bam        -T $output_BWAaln/1_$temp      $output_BWAaln/$temp.unsort.bam         >>$output2_BWAaln/$temp.runLog    2>&1");
    system("rm   $output_BWAaln/$temp.unsort.bam");

    if ( !( -e "$qualimapDir/$temp")     )   { mkdir  "$qualimapDir/$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_BWAaln/$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/$temp    --java-mem-size=5G    >>$output2_BWAaln/$temp.runLog    2>&1");
    system("samstat   $output_BWAaln/$temp.bam           >> $output2_BWAaln/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                        >> $output2_BWAaln/$temp.runLog    2>&1");   

    system("samtools  index       $output_BWAaln/$temp.bam      >>$outDirQC/$temp.runLog  2>&1");
    system("samtools  flagstat    $output_BWAaln/$temp.bam      >>$outDirQC/$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_BWAaln/$temp.bam      >>$outDirQC/$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_BWAaln/$temp.bam      >>$outDirQC/$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_BWAaln/$temp.bam  --basic  --qual  --phred    >>$outDirQC/$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_BWAaln/$temp.bam          >>$FastQCdir/$temp.runLog        2>&1");  

    
}
} ########## End BWAaln
















{ ########## Start BWAmem
print "\n\n\n\n\n##################################################################################################\n";
print("\nChecking all the input file names......\n");
opendir(my $DH_input, "$input_BWAmem")  ||  die;     
my @inputFiles = readdir($DH_input);

my $fileNameBool = 1;
for ( my $i=0; $i<=$#inputFiles; $i++ ) {  
        next unless $inputFiles[$i] =~ m/^(\d{2})_/;   
        next unless $inputFiles[$i] =~ m/\.bam$/;   
        next unless $inputFiles[$i] !~ m/^[.]/;
        next unless $inputFiles[$i] !~ m/[~]$/;
        next unless $inputFiles[$i] !~ m/^unpaired/;
        my $temp = $inputFiles[$i]; 
        $temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]).bam$/   or  die  "wrong-1: ## $temp ##";
        if($temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])\.bam$/) {
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
print("\nAnalysis all the bam files ......\n");
for(my $i=0; $i<=$#inputFiles; $i++) {
    next unless $inputFiles[$i] =~ m/^(\d{2})_/;   
    next unless $inputFiles[$i] =~ m/\.bam$/;
    next unless $inputFiles[$i] !~ m/^[.]/;
    next unless $inputFiles[$i] !~ m/[~]$/;
    print("$inputFiles[$i] ......\n");      
    my $temp = $inputFiles[$i];
    $temp =~ s/\.bam$//  or  die;

    system(`samtools  view  -h  $input_BWAmem/$temp.bam  | wc  -l     >>$output2_BWAmem/$temp.number    2>&1`);
    system(`samtools  view   -h   -q 30   -o $output_BWAmem/$temp.q30.sam       $input_BWAmem/$temp.bam      >>$output2_BWAmem/$temp.runLog    2>&1`);    ## must use `,  cann't use "
    system(`wc  -l  $output_BWAmem/$temp.q30.sam  >>$output2_BWAmem/$temp.number    2>&1`);  
    if ( !( -e "$output_BWAmem/1_$temp")     )   { mkdir "$output_BWAmem/1_$temp"        ||  die; }
    system("removeDup    -r 1000    -t $output_BWAmem/1_$temp    -c 999999999     -i  $output_BWAmem/$temp.q30.sam   -o $output_BWAmem/$temp.noDup.sam            >>$output2_BWAmem/$temp.runLog    2>&1");  ## SAM file must contain header. 
    system("rm   -rf  $output_BWAmem/1_$temp");
    system(`wc  -l  $output_BWAmem/$temp.noDup.sam  >>$output2_BWAmem/$temp.number    2>&1`);  
    system("repair    -S  -c  -T 16    -i  $output_BWAmem/$temp.noDup.sam   -o $output_BWAmem/$temp.unsort.bam            >>$output2_BWAmem/$temp.runLog    2>&1");  ## SAM file must contain header.         
    system(`samtools  view  -h  $output_BWAmem/$temp.unsort.bam  | wc  -l     >>$output2_BWAmem/$temp.number    2>&1`);
    system("rm   $output_BWAmem/$temp.q30.sam");
    system("rm   $output_BWAmem/$temp.noDup.sam");
    system("samtools  sort   -O bam       -o $output_BWAmem/$temp.bam        -T $output_BWAmem/1_$temp      $output_BWAmem/$temp.unsort.bam         >>$output2_BWAmem/$temp.runLog    2>&1");
    system("rm   $output_BWAmem/$temp.unsort.bam");

    if ( !( -e "$qualimapDir/$temp")     )   { mkdir  "$qualimapDir/$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_BWAmem/$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/$temp    --java-mem-size=5G    >>$output2_BWAmem/$temp.runLog    2>&1");
    system("samstat   $output_BWAmem/$temp.bam           >> $output2_BWAmem/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                        >> $output2_BWAmem/$temp.runLog    2>&1");   

    system("samtools  index       $output_BWAmem/$temp.bam      >>$outDirQC/$temp.runLog  2>&1");
    system("samtools  flagstat    $output_BWAmem/$temp.bam      >>$outDirQC/$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_BWAmem/$temp.bam      >>$outDirQC/$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_BWAmem/$temp.bam      >>$outDirQC/$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_BWAmem/$temp.bam  --basic  --qual  --phred    >>$outDirQC/$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_BWAmem/$temp.bam          >>$FastQCdir/$temp.runLog        2>&1");  

    
}
} ########## End BWAmem












{ ########## Start Bowtie1
print "\n\n\n\n\n##################################################################################################\n";
print("\nChecking all the input file names......\n");
opendir(my $DH_input, "$input_Bowtie1")  ||  die;     
my @inputFiles = readdir($DH_input);

my $fileNameBool = 1;
for ( my $i=0; $i<=$#inputFiles; $i++ ) {  
        next unless $inputFiles[$i] =~ m/^(\d{2})_/;   
        next unless $inputFiles[$i] =~ m/\.bam$/;   
        next unless $inputFiles[$i] !~ m/^[.]/;
        next unless $inputFiles[$i] !~ m/[~]$/;
        next unless $inputFiles[$i] !~ m/^unpaired/;
        my $temp = $inputFiles[$i]; 
        $temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]).bam$/   or  die  "wrong-1: ## $temp ##";
        if($temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])\.bam$/) {
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
print("\nAnalysis all the bam files ......\n");
for(my $i=0; $i<=$#inputFiles; $i++) {
    next unless $inputFiles[$i] =~ m/^(\d{2})_/;   
    next unless $inputFiles[$i] =~ m/\.bam$/;
    next unless $inputFiles[$i] !~ m/^[.]/;
    next unless $inputFiles[$i] !~ m/[~]$/;
    print("$inputFiles[$i] ......\n");      
    my $temp = $inputFiles[$i];
    $temp =~ s/\.bam$//  or  die;

    system(`samtools  view  -h  $input_Bowtie1/$temp.bam  | wc  -l     >>$output2_Bowtie1/$temp.number    2>&1`);
    system(`samtools  view   -h   -q 30   -o $output_Bowtie1/$temp.q30.sam       $input_Bowtie1/$temp.bam      >>$output2_Bowtie1/$temp.runLog    2>&1`);    ## must use `,  cann't use "
    system(`wc  -l  $output_Bowtie1/$temp.q30.sam  >>$output2_Bowtie1/$temp.number    2>&1`);  
    if ( !( -e "$output_Bowtie1/1_$temp")     )   { mkdir "$output_Bowtie1/1_$temp"        ||  die; }
    system("removeDup    -r 1000    -t $output_Bowtie1/1_$temp    -c 999999999     -i  $output_Bowtie1/$temp.q30.sam   -o $output_Bowtie1/$temp.noDup.sam            >>$output2_Bowtie1/$temp.runLog    2>&1");  ## SAM file must contain header. 
    system("rm   -rf  $output_Bowtie1/1_$temp");
    system(`wc  -l  $output_Bowtie1/$temp.noDup.sam  >>$output2_Bowtie1/$temp.number    2>&1`);  
    system("repair    -S  -c  -T 16    -i  $output_Bowtie1/$temp.noDup.sam   -o $output_Bowtie1/$temp.unsort.bam            >>$output2_Bowtie1/$temp.runLog    2>&1");  ## SAM file must contain header.         
    system(`samtools  view  -h  $output_Bowtie1/$temp.unsort.bam  | wc  -l     >>$output2_Bowtie1/$temp.number    2>&1`);
    system("rm   $output_Bowtie1/$temp.q30.sam");
    system("rm   $output_Bowtie1/$temp.noDup.sam");
    system("samtools  sort   -O bam       -o $output_Bowtie1/$temp.bam        -T $output_Bowtie1/1_$temp      $output_Bowtie1/$temp.unsort.bam         >>$output2_Bowtie1/$temp.runLog    2>&1");
    system("rm   $output_Bowtie1/$temp.unsort.bam");

    if ( !( -e "$qualimapDir/$temp")     )   { mkdir  "$qualimapDir/$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_Bowtie1/$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/$temp    --java-mem-size=5G    >>$output2_Bowtie1/$temp.runLog    2>&1");
    system("samstat   $output_Bowtie1/$temp.bam           >> $output2_Bowtie1/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                        >> $output2_Bowtie1/$temp.runLog    2>&1");   

    system("samtools  index       $output_Bowtie1/$temp.bam      >>$outDirQC/$temp.runLog  2>&1");
    system("samtools  flagstat    $output_Bowtie1/$temp.bam      >>$outDirQC/$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_Bowtie1/$temp.bam      >>$outDirQC/$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_Bowtie1/$temp.bam      >>$outDirQC/$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_Bowtie1/$temp.bam  --basic  --qual  --phred    >>$outDirQC/$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_Bowtie1/$temp.bam          >>$FastQCdir/$temp.runLog        2>&1");  

    
}
} ########## End Bowtie1













{ ########## Start Bowtie2
print "\n\n\n\n\n##################################################################################################\n";
print("\nChecking all the input file names......\n");
opendir(my $DH_input, "$input_Bowtie2")  ||  die;     
my @inputFiles = readdir($DH_input);

my $fileNameBool = 1;
for ( my $i=0; $i<=$#inputFiles; $i++ ) {  
        next unless $inputFiles[$i] =~ m/^(\d{2})_/;   
        next unless $inputFiles[$i] =~ m/\.bam$/;   
        next unless $inputFiles[$i] !~ m/^[.]/;
        next unless $inputFiles[$i] !~ m/[~]$/;
        next unless $inputFiles[$i] !~ m/^unpaired/;
        my $temp = $inputFiles[$i]; 
        $temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]).bam$/   or  die  "wrong-1: ## $temp ##";
        if($temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])\.bam$/) {
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
print("\nAnalysis all the bam files ......\n");
for(my $i=0; $i<=$#inputFiles; $i++) {
    next unless $inputFiles[$i] =~ m/^(\d{2})_/;   
    next unless $inputFiles[$i] =~ m/\.bam$/;
    next unless $inputFiles[$i] !~ m/^[.]/;
    next unless $inputFiles[$i] !~ m/[~]$/;
    print("$inputFiles[$i] ......\n");      
    my $temp = $inputFiles[$i];
    $temp =~ s/\.bam$//  or  die;

    system(`samtools  view  -h  $input_Bowtie2/$temp.bam  | wc  -l     >>$output2_Bowtie2/$temp.number    2>&1`);
    system(`samtools  view   -h   -q 30   -o $output_Bowtie2/$temp.q30.sam       $input_Bowtie2/$temp.bam      >>$output2_Bowtie2/$temp.runLog    2>&1`);    ## must use `,  cann't use "
    system(`wc  -l  $output_Bowtie2/$temp.q30.sam  >>$output2_Bowtie2/$temp.number    2>&1`);  
    if ( !( -e "$output_Bowtie2/1_$temp")     )   { mkdir "$output_Bowtie2/1_$temp"        ||  die; }
    system("removeDup    -r 1000    -t $output_Bowtie2/1_$temp    -c 999999999     -i  $output_Bowtie2/$temp.q30.sam   -o $output_Bowtie2/$temp.noDup.sam            >>$output2_Bowtie2/$temp.runLog    2>&1");  ## SAM file must contain header. 
    system("rm   -rf  $output_Bowtie2/1_$temp");
    system(`wc  -l  $output_Bowtie2/$temp.noDup.sam  >>$output2_Bowtie2/$temp.number    2>&1`);  
    system("repair    -S  -c  -T 16    -i  $output_Bowtie2/$temp.noDup.sam   -o $output_Bowtie2/$temp.unsort.bam            >>$output2_Bowtie2/$temp.runLog    2>&1");  ## SAM file must contain header.         
    system(`samtools  view  -h  $output_Bowtie2/$temp.unsort.bam  | wc  -l     >>$output2_Bowtie2/$temp.number    2>&1`);
    system("rm   $output_Bowtie2/$temp.q30.sam");
    system("rm   $output_Bowtie2/$temp.noDup.sam");
    system("samtools  sort   -O bam       -o $output_Bowtie2/$temp.bam        -T $output_Bowtie2/1_$temp      $output_Bowtie2/$temp.unsort.bam         >>$output2_Bowtie2/$temp.runLog    2>&1");
    system("rm   $output_Bowtie2/$temp.unsort.bam");

    if ( !( -e "$qualimapDir/$temp")     )   { mkdir  "$qualimapDir/$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_Bowtie2/$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/$temp    --java-mem-size=5G    >>$output2_Bowtie2/$temp.runLog    2>&1");
    system("samstat   $output_Bowtie2/$temp.bam           >> $output2_Bowtie2/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                        >> $output2_Bowtie2/$temp.runLog    2>&1");   

    system("samtools  index       $output_Bowtie2/$temp.bam      >>$outDirQC/$temp.runLog  2>&1");
    system("samtools  flagstat    $output_Bowtie2/$temp.bam      >>$outDirQC/$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_Bowtie2/$temp.bam      >>$outDirQC/$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_Bowtie2/$temp.bam      >>$outDirQC/$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_Bowtie2/$temp.bam  --basic  --qual  --phred    >>$outDirQC/$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_Bowtie2/$temp.bam          >>$FastQCdir/$temp.runLog        2>&1");  

    
}
} ########## End Bowtie2









print "\n\n\n\n\n##################################################################################################\n";
print "\n\n        Job Done! Cheers! \n\n";





























