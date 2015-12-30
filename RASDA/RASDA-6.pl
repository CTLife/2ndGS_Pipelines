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
        Welcome to use RASDA (RNA-Seq Data Analyzer), version 0.62, 2016-01-13.      
        RASDA is a Pipeline for Single-end and Paired-end RNA-Seq Data Analysis by Integrating Lots of Softwares.

        Step 6: Only  retain the mapped reads with MAPQ>=30 and no PCR duplicates.  
                And repair the unpaired reads in paired-end file.
                Quality statistics by using FastQC, BamUtil, SAMtools, QualiMap and samstat.
        Usage:  
               perl  RASDA-6.pl    [-v]    [-h]    [-i inputDir]    [-o outDir]   
        For instance: 
                     perl  RASDA-6.pl    -i 6-MisMatch          -o 7-FinalBAM              
                     perl  RASDA-6.pl    --input 6-MisMatch     --output 7-FinalBAM     
                     perl  RASDA-6.pl    --input 6-MisMatch     --output 7-FinalBAM    >> RASDA-6.runLog  2>&1
     
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

        For more details about this pipeline and other NGS data analysis piplines such as CISDA, MESDA and HISDA,
        please see the web site: https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies 
        and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.    
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------  
';


########## Version Infromation ##########
my $version_g = "  The Sixth Step of RASDA (RNA-Seq Data Analyzer), version 0.62, 2016-01-13.";


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
    print   '    Please see help message by using "perl  RASDA-6.pl  -h". ';
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









my $input_HISAT2   = "$input_g/1-HISAT2";
my $input_Subjunc  = "$input_g/2-Subjunc";
my $input_STAR     = "$input_g/3-STAR";

my $output_HISAT2  = "$output_g/1-HISAT2";
my $output_Subjunc = "$output_g/2-Subjunc";
my $output_STAR    = "$output_g/3-STAR";

my $output2_HISAT2  = "$output_g/1-HISAT2/Results";
my $output2_Subjunc = "$output_g/2-Subjunc/Results";
my $output2_STAR    = "$output_g/3-STAR/Results";

if ( !(-e $output_HISAT2)   )   { mkdir $output_HISAT2     ||  die; }
if ( !(-e $output_Subjunc)  )   { mkdir $output_Subjunc    ||  die; }
if ( !(-e $output_STAR)     )   { mkdir $output_STAR       ||  die; }

if ( !(-e $output2_HISAT2)  )   { mkdir $output2_HISAT2    ||  die; }
if ( !(-e $output2_Subjunc) )   { mkdir $output2_Subjunc   ||  die; }
if ( !(-e $output2_STAR)    )   { mkdir $output2_STAR      ||  die; }










{ ########## Start HISAT2
print "\n\n\n\n\n##################################################################################################\n";
print("\nChecking all the input file names......\n");
opendir(my $DH_input, "$input_HISAT2")  ||  die;     
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
 
my $qualimapDir = "$output2_HISAT2/qualimap"; 
my $FastQCdir   = "$output2_HISAT2/FastQC";  
my $outDirQC    = "$output2_HISAT2/QCstatistics";

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

    system(`samtools  view  -h  $input_HISAT2/$temp.bam  | wc  -l     >>$output2_HISAT2/$temp.number    2>&1`);
    system(`samtools  view   -h   -q 30   -o $output_HISAT2/$temp.q30.sam       $input_HISAT2/$temp.bam      >>$output2_HISAT2/$temp.runLog    2>&1`);    ## must use `,  cann't use "
    system(`wc  -l  $output_HISAT2/$temp.q30.sam  >>$output2_HISAT2/$temp.number    2>&1`);  
    if ( !( -e "$output_HISAT2/1_$temp")     )   { mkdir "$output_HISAT2/1_$temp"        ||  die; }
    system("removeDup    -r 100000    -t $output_HISAT2/1_$temp    -c 999999999     -i  $output_HISAT2/$temp.q30.sam   -o $output_HISAT2/$temp.noDup.sam            >>$output2_HISAT2/$temp.runLog    2>&1");  ## SAM file must contain header. 
    system("rm   -rf  $output_HISAT2/1_$temp");
    system(`wc  -l  $output_HISAT2/$temp.noDup.sam  >>$output2_HISAT2/$temp.number    2>&1`);  
    system("repair    -S  -c  -T 16    -i  $output_HISAT2/$temp.noDup.sam   -o $output_HISAT2/$temp.unsort.bam            >>$output2_HISAT2/$temp.runLog    2>&1");  ## SAM file must contain header.         
    system(`samtools  view  -h  $output_HISAT2/$temp.unsort.bam  | wc  -l     >>$output2_HISAT2/$temp.number    2>&1`);
    system("rm   $output_HISAT2/$temp.q30.sam");
    system("rm   $output_HISAT2/$temp.noDup.sam");
    system("samtools  sort   -O bam       -o $output_HISAT2/$temp.bam        -T $output_HISAT2/1_$temp      $output_HISAT2/$temp.unsort.bam         >>$output2_HISAT2/$temp.runLog    2>&1");
    system("rm   $output_HISAT2/$temp.unsort.bam");

    if ( !( -e "$qualimapDir/$temp")     )   { mkdir  "$qualimapDir/$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_HISAT2/$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/$temp    --java-mem-size=5G    >>$output2_HISAT2/$temp.runLog    2>&1");
    system("samstat   $output_HISAT2/$temp.bam           >> $output2_HISAT2/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                        >> $output2_HISAT2/$temp.runLog    2>&1");   

    system("samtools  index       $output_HISAT2/$temp.bam      >>$outDirQC/$temp.runLog  2>&1");
    system("samtools  flagstat    $output_HISAT2/$temp.bam      >>$outDirQC/$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_HISAT2/$temp.bam      >>$outDirQC/$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_HISAT2/$temp.bam      >>$outDirQC/$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_HISAT2/$temp.bam  --basic  --qual  --phred    >>$outDirQC/$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_HISAT2/$temp.bam          >>$FastQCdir/$temp.runLog        2>&1");  

    
}
} ########## End HISAT2











{ ########## Start Subjunc
print "\n\n\n\n\n##################################################################################################\n";
print("\nChecking all the input file names......\n");
opendir(my $DH_input, "$input_Subjunc")  ||  die;     
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
 
my $qualimapDir = "$output2_Subjunc/qualimap"; 
my $FastQCdir   = "$output2_Subjunc/FastQC";  
my $outDirQC    = "$output2_Subjunc/QCstatistics";

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

    system(`samtools  view  -h  $input_Subjunc/$temp.bam  | wc  -l     >>$output2_Subjunc/$temp.number    2>&1`);
    system(`samtools  view   -h   -q 30   -o $output_Subjunc/$temp.q30.sam       $input_Subjunc/$temp.bam      >>$output2_Subjunc/$temp.runLog    2>&1`);    ## must use `,  cann't use "
    system(`wc  -l  $output_Subjunc/$temp.q30.sam  >>$output2_Subjunc/$temp.number    2>&1`);  
    if ( !( -e "$output_Subjunc/1_$temp")     )   { mkdir "$output_Subjunc/1_$temp"        ||  die; }
    system("removeDup    -r 100000    -t $output_Subjunc/1_$temp    -c 999999999     -i  $output_Subjunc/$temp.q30.sam   -o $output_Subjunc/$temp.noDup.sam            >>$output2_Subjunc/$temp.runLog    2>&1");  ## SAM file must contain header. 
    system("rm   -rf  $output_Subjunc/1_$temp");
    system(`wc  -l  $output_Subjunc/$temp.noDup.sam  >>$output2_Subjunc/$temp.number    2>&1`);  
    system("repair    -S  -c  -T 16    -i  $output_Subjunc/$temp.noDup.sam   -o $output_Subjunc/$temp.unsort.bam            >>$output2_Subjunc/$temp.runLog    2>&1");  ## SAM file must contain header.         
    system(`samtools  view  -h  $output_Subjunc/$temp.unsort.bam  | wc  -l     >>$output2_Subjunc/$temp.number    2>&1`);
    system("rm   $output_Subjunc/$temp.q30.sam");
    system("rm   $output_Subjunc/$temp.noDup.sam");
    system("samtools  sort   -O bam       -o $output_Subjunc/$temp.bam        -T $output_Subjunc/1_$temp      $output_Subjunc/$temp.unsort.bam         >>$output2_Subjunc/$temp.runLog    2>&1");
    system("rm   $output_Subjunc/$temp.unsort.bam");


    if ( !( -e "$qualimapDir/$temp")     )   { mkdir  "$qualimapDir/$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_Subjunc/$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/$temp    --java-mem-size=5G    >>$output2_Subjunc/$temp.runLog    2>&1");
    system("samstat   $output_Subjunc/$temp.bam           >> $output2_Subjunc/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                        >> $output2_Subjunc/$temp.runLog    2>&1");   

    system("samtools  index       $output_Subjunc/$temp.bam      >>$outDirQC/$temp.runLog  2>&1");
    system("samtools  flagstat    $output_Subjunc/$temp.bam      >>$outDirQC/$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_Subjunc/$temp.bam      >>$outDirQC/$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_Subjunc/$temp.bam      >>$outDirQC/$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_Subjunc/$temp.bam  --basic  --qual  --phred    >>$outDirQC/$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_Subjunc/$temp.bam          >>$FastQCdir/$temp.runLog        2>&1");  

    
}
} ########## End Subjunc













{ ########## Start STAR
print "\n\n\n\n\n##################################################################################################\n";
print("\nChecking all the input file names......\n");
opendir(my $DH_input, "$input_STAR")  ||  die;     
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
 
my $qualimapDir = "$output2_STAR/qualimap"; 
my $FastQCdir   = "$output2_STAR/FastQC";  
my $outDirQC    = "$output2_STAR/QCstatistics";

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

    system(`samtools  view  -h  $input_STAR/$temp.bam  | wc  -l     >>$output2_STAR/$temp.number    2>&1`);
    system(`samtools  view   -h   -q 30   -o $output_STAR/$temp.q30.sam       $input_STAR/$temp.bam      >>$output2_STAR/$temp.runLog    2>&1`);    ## must use `,  cann't use "
    system(`wc  -l  $output_STAR/$temp.q30.sam  >>$output2_STAR/$temp.number    2>&1`);  
    if ( !( -e "$output_STAR/1_$temp")     )   { mkdir "$output_STAR/1_$temp"        ||  die; }
    system("removeDup    -r 100000    -t $output_STAR/1_$temp    -c 999999999     -i  $output_STAR/$temp.q30.sam   -o $output_STAR/$temp.noDup.sam            >>$output2_STAR/$temp.runLog    2>&1");  ## SAM file must contain header. 
    system("rm   -rf  $output_STAR/1_$temp");
    system(`wc  -l  $output_STAR/$temp.noDup.sam  >>$output2_STAR/$temp.number    2>&1`);  
    system("repair    -S  -c  -T 16    -i  $output_STAR/$temp.noDup.sam   -o $output_STAR/$temp.unsort.bam            >>$output2_STAR/$temp.runLog    2>&1");  ## SAM file must contain header.         
    system(`samtools  view  -h  $output_STAR/$temp.unsort.bam  | wc  -l     >>$output2_STAR/$temp.number    2>&1`);
    system("rm   $output_STAR/$temp.q30.sam");
    system("rm   $output_STAR/$temp.noDup.sam");
    system("samtools  sort   -O bam       -o $output_STAR/$temp.bam        -T $output_STAR/1_$temp      $output_STAR/$temp.unsort.bam         >>$output2_STAR/$temp.runLog    2>&1");
    system("rm   $output_STAR/$temp.unsort.bam");


    if ( !( -e "$qualimapDir/$temp")     )   { mkdir  "$qualimapDir/$temp"   ||  die; }       
    system("qualimap  bamqc  -bam $output_STAR/$temp.bam   -c    -gd MOUSE  -outdir $qualimapDir/$temp    --java-mem-size=5G    >>$output2_STAR/$temp.runLog    2>&1");
    system("samstat   $output_STAR/$temp.bam           >> $output2_STAR/$temp.runLog    2>&1");
    system("echo    '\n\n\n\n\n\n'                        >> $output2_STAR/$temp.runLog    2>&1");   

    system("samtools  index       $output_STAR/$temp.bam      >>$outDirQC/$temp.runLog  2>&1");
    system("samtools  flagstat    $output_STAR/$temp.bam      >>$outDirQC/$temp.runLog  2>&1");
    system(`samtools  idxstats    $output_STAR/$temp.bam      >>$outDirQC/$temp.runLog  2>&1`);
    system(`bam  validate  --in   $output_STAR/$temp.bam      >>$outDirQC/$temp.run=bam.runLog       2>&1`);
    system(`bam  stats     --in   $output_STAR/$temp.bam  --basic  --qual  --phred    >>$outDirQC/$temp.run=bam.runLog 2>&1`);
    system("fastqc    --outdir $FastQCdir          --threads 16    --format bam    --kmers 7     $output_STAR/$temp.bam          >>$FastQCdir/$temp.runLog        2>&1");  

    
}
} ########## End STAR









print "\n\n\n\n\n##################################################################################################\n";
print "\n\n        Job Done! Cheers! \n\n";





























