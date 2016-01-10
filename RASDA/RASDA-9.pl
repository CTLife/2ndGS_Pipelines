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

        Step 9: For each RNA-Seq sample, run cufflinks to assemble the read alignments obtained in the previous step.
        Usage:  
               perl  RASDA-9.pl    [-v]    [-h]    [-i inputDir]    [-o outDir]   
        For instance: 
                     perl  RASDA-9.pl    -i 7-FinalBAM/2-Subjunc          -o 10-Cufflinks                      
                     perl  RASDA-9.pl    --input 7-FinalBAM/2-Subjunc     --output 10-Cufflinks     
                     perl  RASDA-9.pl    --input 7-FinalBAM/2-Subjunc     --output 10-Cufflinks    >> RASDA-9.runLog  2>&1
     
        -------------------------------------------------------------------------------------------------------------------
        Optional arguments:
        -v, --version        Show version number of this program and exit.

        -h, --help           Show this help message and exit.

        Required arguments:
        -i inputDir,  --input inputDir        inputDir is the name of your input folder that contains your BAM files,
                                              the suffix of the BAM files must be ".bam".    (no default)

        -o outDir,  --output outDir           outDir is the name of your output folder that contains running 
                                              results of this step.      (no default) 
        ------------------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines such as CISDA, MESDA and HISDA,
        please see the web site: https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies 
        and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.    
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------  
';


########## Version Infromation ##########
my $version_g = "  The Nineth Step of RASDA (RNA-Seq Data Analyzer), version 0.62, 2016-01-13.";


########## Keys and Values ##########
if ($#ARGV   == -1) { print  "\n$HELP_g\n\n";  exit 0; }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0) { @ARGV = (@ARGV, "-h");           }       ## when the number of command argumants is odd. 
my %args = @ARGV;


########## Initialize  Variables ##########
my $input_g   = '7-FinalBAM/2-Subjunc';      ## This is only an initialization  value or suggesting value, not default value.
my $output_g  = '10-Cufflinks';                  ## This is only an initialization  value or suggesting value, not default value.


########## Available Arguments ##########
my $available = "  -v  --version    -h  --help    -i  --input    -o    --output     -L  --fragLen     ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if($available !~ m/\s$key\s/) {print  "    Cann't recognize $key !!\n";  $boole_g = 1; }
}
if($boole_g == 1) {
    print "\n    The Command Line Arguments are wrong!\n";
    print   '    Please see help message by using "perl  RASDA-9.pl  -h". ';
    print "\n\n";
    exit 0;
}


########## Get Arguments ##########
if ( ( exists $args{'-v' } )  or  ( exists $args{'--version'      } )  )     { print  "\n$version_g\n\n";    exit 0; }
if ( ( exists $args{'-h' } )  or  ( exists $args{'--help'         } )  )     { print  "\n$HELP_g\n\n";       exit 0; }
if ( ( exists $args{'-i' } )  or  ( exists $args{'--input'        } )  )     { ($input_g  = $args{'-i' })   or  ($input_g  = $args{'--input'  });  }else{print   "\n -i or --input   is required.\n\n";   print  "\n$HELP_g\n\n";       exit 0; }                                               
if ( ( exists $args{'-o' } )  or  ( exists $args{'--output'       } )  )     { ($output_g = $args{'-o' })   or  ($output_g = $args{'--output' });  }else{print   "\n -o or --output  is required.\n\n";   print  "\n$HELP_g\n\n";       exit 0; }      


########### Conditions #############
$input_g   =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";
$output_g  =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";


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

system("Cufflinks                  >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'     >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");
system("cuffmerge    -h            >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'     >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");
system("cuffquant                  >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'     >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");
system("cuffdiff                   >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'     >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");
system("cuffnorm                   >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'     >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");
system("cuffcompare                >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'     >> $output2_g/z-version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $output2_g/z-version_softwares.txt   2>&1");







opendir(my $DH_input, $input_g) || die;     
my @inputFiles = readdir($DH_input);

print("\nChecking all the input file names......\n");
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




my $outDir1 = "$output_g/1-cufflinks";
my $outDir2 = "$output_g/2-cuffmerge";
my $outDir3 = "$output_g/3-cuffquant";
my $outDir4 = "$output_g/4-cuffdiff";
my $outDir5 = "$output_g/5-CummeRbund";
my $outDir6 = "$output_g/6-cuffnorm";
my $outDir7 = "$output_g/7-cuffcompare";
if ( !(-e $outDir1) )  {mkdir $outDir1 or die; }
if ( !(-e $outDir2) )  {mkdir $outDir2 or die; }
if ( !(-e $outDir3) )  {mkdir $outDir3 or die; }
if ( !(-e $outDir4) )  {mkdir $outDir4 or die; }
if ( !(-e $outDir5) )  {mkdir $outDir5 or die; }
if ( !(-e $outDir6) )  {mkdir $outDir6 or die; }
if ( !(-e $outDir7) )  {mkdir $outDir7 or die; }





if(1==1) {   ###########################################################

print "\n\n(1) Assemble the read alignments ......\n";
## step1, cufflinks (with novel transcripts)
open(FH1, ">", "$outDir1/assemblies.txt"    )  or  die; 
for(my $i=0; $i<=$#inputFiles; $i++) {
    next unless $inputFiles[$i] =~ m/\.bam$/;
    next unless $inputFiles[$i] !~ m/^[.]/;
    next unless $inputFiles[$i] !~ m/[~]$/;
    my $temp = $inputFiles[$i];
    $temp =~ s/\.bam$//  or  die;
    system("cufflinks    --output-dir $outDir1/$temp        --num-threads 8     --GTF-guide 0-Other/mm9-RefSeq.gtf     --frag-bias-correct 0-Other/mm9.fa     --multi-read-correct    --library-type fr-unstranded          --label $temp     $input_g/$temp.bam     >>$outDir1/$temp.runLog      2>&1"   );                                     
    print  FH1    "$outDir1/$temp/transcripts.gtf\n";
}



print "\n\n(2) cuffmerge ......\n";
## step2, cuffmerge
system("cuffmerge      -o  $outDir2       --ref-gtf  0-Other/mm9-RefSeq.gtf     --ref-sequence 0-Other/mm9.fa     --num-threads 8   --keep-tmp  $outDir1/assemblies.txt      >>$outDir2/runLog.txt        2>&1");




print "\n\n(3) cuffquant ......\n";
## step3,    cuffquant 
for(my $i=0; $i<=$#inputFiles; $i++) {   
    next unless $inputFiles[$i] =~ m/\.bam$/; 
    next unless $inputFiles[$i] !~ m/^[.]/;
    next unless $inputFiles[$i] !~ m/[~]$/; 
    my $temp = $inputFiles[$i];
    $temp =~ s/\.bam$//  or  die;
    system("cuffquant     --output-dir $outDir3/$temp      --frag-bias-correct 0-Other/mm9.fa      --multi-read-correct    --num-threads 8       --library-type fr-unstranded         $outDir2/merged.gtf    $input_g/$temp.bam     >>$outDir3/$temp.runLog      2>&1"); 
}





} ###########################################################





my $BAM_Condition1 = "$input_g/01_mRNA_Hand2as_CM-6_L8-I20006_shenLab_Rep1.bam,$input_g/02_mm9_mRNA_adult_cardiomyocyte_SRR1614298_Rep2.bam";
my $BAM_Condition2 = "$input_g/03_CM_2014NC_E12.5_RNAseq_Gata4KO_Rep1.bam";
my $BAM_Condition3 = "";

my $CXB_Condition1 = "$outDir3/01_mRNA_Hand2as_CM-6_L8-I20006_shenLab_Rep1/abundances.cxb,$outDir3/02_mm9_mRNA_adult_cardiomyocyte_SRR1614298_Rep2/abundances.cxb";
my $CXB_Condition2 = "$outDir3/03_CM_2014NC_E12.5_RNAseq_Gata4KO_Rep1/abundances.cxb";
my $CXB_Condition3 = "";


print "\n\n(4) cuffdiff ......\n";
## step4,    cuffdiff 
system("cuffdiff   --output-dir $outDir4/BAM    --labels ctrl,KO     --num-threads 8    --library-type  fr-unstranded   $outDir2/merged.gtf    $BAM_Condition1      $BAM_Condition2          >>$outDir4/BAM-log.txt        2>&1");
system("cuffdiff   --output-dir $outDir4/CXB    --labels ctrl,KO     --num-threads 8    --library-type  fr-unstranded   $outDir2/merged.gtf    $CXB_Condition1      $CXB_Condition2          >>$outDir4/CXB-log.txt        2>&1");





print "\n\n(5) CummeRbund ......\n";
## step5,    CummeRbund




print "\n\n(6) cuffnorm ......\n";
## step6,    cuffnorm 
system("cuffnorm   --output-dir $outDir6/BAM    --labels ctrl,KO     --frag-bias-correct 0-Other/mm9.fa     --multi-read-correct   --num-threads 6         $outDir2/merged.gtf    $BAM_Condition1      $BAM_Condition2      >>$outDir6/BAM-log.txt        2>&1");
system("cuffnorm   --output-dir $outDir6/CXB    --labels ctrl,KO     --frag-bias-correct 0-Other/mm9.fa     --multi-read-correct   --num-threads 6         $outDir2/merged.gtf    $CXB_Condition1      $CXB_Condition2      >>$outDir6/CXB-log.txt        2>&1");                              




print "\n\n(7) cuffcompare ......\n";
## step7,    cuffcompare 
system("cuffcompare     -r 0-Other/mm9-RefSeq.gtf     -s 0-Other/mm9.fa       -o $outDir7     -i $outDir1/assemblies.txt   >>$outDir7/log.txt        2>&1");  





print "\n\n\n\n\n##################################################################################################\n";
print "\n\n        Job Done! Cheers! \n\n";












