#!/usr/bin/env perl5
use  strict;
use  warnings;






###################################################################################################################################################################################################
###################################################################################################################################################################################################

########## Help Infromation ##########
my $HELP_g = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use CISDA (ChIP-Seq Data Analyzer), version 0.51, 2016-01-01.    
        CISDA is a Pipeline for Single-end and Paired-end ChIP-Seq Data Analysis by Integrating Lots of Softwares.

        Step 1: Decompress all compressed FASTQ files, or convert SRA to FASTQ using SRA_Toolkit,  
                and merge the two lanes of the same sample (merge technical replicates), 
                and quality statistics by using FastQC, NGS_QC_Toolkit and FASTX-Toolkit.
  
        Usage: 
               perl  CISDA-1.pl    [-v]    [-h]    [-i inputDir]    [-o outDir]  

        For instance: 
                     perl  CISDA-1.pl    -i 1-rawReads          -o 2-FASTQ           
                     perl  CISDA-1.pl    --input 1-rawReads     --output 2-FASTQ    
                     perl  CISDA-1.pl    --input 1-rawReads     --output 2-FASTQ    >> 1-runLog.txt  2>&1
 
        
        ----------------------------------------------------------------------------------------------------------
        Optional arguments:

        -v, --version        Show version number of this program and exit.

        -h, --help           Show this help message and exit.


        Required arguments:

        -i inputDir,  --input inputDir        "inputDir" is the name of input folder that contains your compressed 
                                              fastq files or SRA files.  
                                              The suffix of the compressed fastq files will be recognized. (no default)

        -o outDir,  --output outDir           "outDir" is the name of output folder that contains your running 
                                              results (fastq format) of this step.      (no default)               
        -----------------------------------------------------------------------------------------------------------


        For more details about this pipeline and other NGS data analysis piplines such as RASDA, MESDA and HISDA,
        please see the web site:   
                                 https://github.com/CTLife/2ndGS_Pipelines


        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies 
        and Center for Life Sciences (CLS), Peking University, China.  
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------   
';



########## Version Infromation ##########
my $version_g = "    CISDA (ChIP-Seq Data Analyzer), version 0.51, 2016-01-01.";



########## Keys and Values ##########
if ($#ARGV   == -1)   { print  "\n$HELP_g\n\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-h") ;           }       ## when the number of command argumants is odd. 
my %args = @ARGV;



########## Initialize  Variables ##########
my $input_g  = '1-rawReads';     ## This is only an initialization  value or suggesting value, not default value.
my $output_g = '2-FASTQ';        ## This is only an initialization  value or suggesting value, not default value.



########## Available Arguments ##########
my $available = "      -v  --version      -h  --help      -i  --input       -o    --output      ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {print  "    Cann't recognize $key\n";  $boole_g = 1; }
}
if($boole_g == 1) {
    print "\n    The Command Line Arguments are wrong!\n";
    print '    Please see help message by using "perl  CISDA-1.pl  -h" ';
    print "\n\n";
    exit 0;
}



########## Get Arguments ##########
if ( ( exists $args{'-v' } )  or  ( exists $args{'--version' } )  )     { print  "\n$version_g\n\n";    exit 0; }
if ( ( exists $args{'-h' } )  or  ( exists $args{'--help'    } )  )     { print  "\n$HELP_g\n\n";       exit 0; }
if ( ( exists $args{'-i' } )  or  ( exists $args{'--input'   } )  )     { ($input_g  = $args{'-i' })  or  ($input_g  = $args{'--input' });  }else{print   "\n -i or --input  is required.\n\n";   print  "\n$HELP_g\n\n";       exit 0; }                                               
if ( ( exists $args{'-o' } )  or  ( exists $args{'--output'  } )  )     { ($output_g = $args{'-o' })  or  ($output_g = $args{'--output'});  }else{print   "\n -o or --output is required.\n\n";   print  "\n$HELP_g\n\n";       exit 0; }      



########### Conditions #############
$input_g  =~ m/^\S+$/    ||  die   "\n\n$HELP_g\n\n\n";
$output_g =~ m/^\S+$/    ||  die   "\n\n$HELP_g\n\n\n";




######### Print Command Arguments to Standard Output ###########
print  "\n\n
        ################ Your Arguments ###############################
                Input  Folder:  $input_g
                Output Folder:  $output_g
        ###############################################################  
\n\n";



###################################################################################################################################################################################################
###################################################################################################################################################################################################
















print   "\n\nRunning......\n\n";
opendir(my $DH_input, $input_g)  ||  die;     
my @inputFiles = readdir($DH_input);
my $pattern = "[-.0-9A-Za-z]+";






my $IlluQC_PRLL = "/home/yongp/MyProgramFiles/6-2G-HTS/2-NGSquality/NGSQCToolkit_v2.3.3/QC/IlluQC_PRLL.pl";
print("\n\nChecking all the necessary softwares in this step......\n\n");

system("fastq-dump  -h           >> $input_g/version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $input_g/version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'  >> $input_g/version_softwares.txt   2>&1");

system("fastqc    -v             >> $input_g/version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $input_g/version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'  >> $input_g/version_softwares.txt   2>&1");

system("perl  $IlluQC_PRLL  -h   >> $input_g/version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $input_g/version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'  >> $input_g/version_softwares.txt   2>&1");

system("fastx_quality_stats  -h  >> $input_g/version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $input_g/version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'  >> $input_g/version_softwares.txt   2>&1");

system("fastq_quality_boxplot_graph.sh   -h  >> $input_g/version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $input_g/version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'              >> $input_g/version_softwares.txt   2>&1");

system("fastx_nucleotide_distribution_graph.sh  -h  >> $input_g/version_softwares.txt   2>&1");
system("echo    '##############################################################################'  >> $input_g/version_softwares.txt   2>&1");
system("echo    '\n\n\n\n\n\n'                     >> $input_g/version_softwares.txt   2>&1");








print("\n\nChecking all the input file names......\n");
my $fileNameBool = 1;
for ( my $i=0; $i<=$#inputFiles; $i++ ) {     
        next unless $inputFiles[$i] !~ m/^[.]/;
        next unless $inputFiles[$i] !~ m/[~]$/;
        my $temp = $inputFiles[$i]; 
        $temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])/   or  die  "wrong-1: ## $temp ##";
        $temp =~ m/_(Rep[1-9])\.sra$/  or  $temp =~ m/_(Rep[1-9])_?([1-2]?)(_Lane[1-2])?\.fastq\.(\S+)$/  or  die   "wrong-2: ## $temp ##";
        if($temp !~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))(_[1-2])?(_Lane[1-2])?(\.fastq)?\.(\S+)$/) {
             $fileNameBool = 0;
        }
}
if($fileNameBool == 1)  {print("    All the file names are passed.\n\n");}









print "\n\nConverting SRA files into FASTQ files or Decompressing the compressed fastq files ......\n";
for ( my $i=0; $i<=$#inputFiles; $i++ ) {     
        next unless $inputFiles[$i] !~ m/^[.]/;
        next unless $inputFiles[$i] !~ m/[~]$/;
        my $temp = $inputFiles[$i]; 
        if ($temp =~ m/\.sra$/) {
                $temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])\.sra$/   or  die;
                system("fastq-dump   --split-3   --dumpbase   $input_g/$temp  --outdir $input_g   >> $input_g/$temp.runLog  2>&1");
        }else{
                $temp =~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])_?([1-2]?)(_Lane[1-2])?\.fastq)\.(\S+)$/   or  die;
                my  $tempFastq = $1;
                my  $suffix1   = $11;    ## Only seven compressed formats are supported, their suffixes:  ".bz2",  ".gz",  ".tar.gz",  ".tar",  ".rar",  ".xz",  ".zip"
                my  $tempBool  = 0;
                if($suffix1 eq "bz2"    )  { $tempBool++;  system("bzip2   -cd         $input_g/$temp   >  $input_g/$tempFastq");  }       
                if($suffix1 eq "gz"     )  { $tempBool++;  system("gunzip  --stdout    $input_g/$temp   >  $input_g/$tempFastq");  }  
                if($suffix1 eq "tar.gz" )  { $tempBool++;  system("tar     -xzvf       $input_g/$temp  -C  $input_g");             }  
                if($suffix1 eq "tar"    )  { $tempBool++;  system("tar     -xf         $input_g/$temp  -C  $input_g");             }  
                if($suffix1 eq "rar"    )  { $tempBool++;  system("unrar    e          $input_g/$temp      $input_g");             }  
                if($suffix1 eq "xz"     )  { $tempBool++;  system("xz      -cd         $input_g/$temp   >  $input_g/$tempFastq");  }  
                if($suffix1 eq "zip"    )  { $tempBool++;  system("unzip   -np         $input_g/$temp   >  $input_g/$tempFastq");  }
                if($tempBool  != 1) { print("$temp is wrong!!");  die; } 
        }
}







my $FastQCdir1       = "$input_g/FastQC";
if ( !( -e $FastQCdir1))   { mkdir $FastQCdir1  ||  die; }
if ( !(-e $output_g)   )   { mkdir $output_g    ||  die; }
print "\n\n\n\n\n##################################################################################################\n";
print "\n\n\n\n\nMerging the two lanes of the same sample or copy files ......\n";
opendir($DH_input, $input_g) || die;
my @twoLanes = readdir($DH_input);
for (my $i=0; $i<=$#twoLanes; $i++) {
    next unless $twoLanes[$i] =~ m/\.fastq$/; 
    next unless $twoLanes[$i] !~ m/^[.]/;
    next unless $twoLanes[$i] !~ m/[~]$/; 
    print  "\n        $twoLanes[$i]\n";
    $twoLanes[$i] =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])_?([1-2]?)(_Lane[1-2])?\.fastq$/   or  die;
    if( $twoLanes[$i] =~ m/^(\S+)_Lane1.fastq$/ ) {
        my $temp = $1;   
        $temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])_?([1-2]?)$/   or  die;
        open(tempFH, ">>", "$output_g/merge-two-lanes.log.txt")  or  die;
        my $lane1 = $temp."_Lane1.fastq";
        my $lane2 = $temp."_Lane2.fastq";
        system(       "cat  $input_g/$lane1  $input_g/$lane2   > $output_g/$temp.fastq" ); 
        print tempFH  "cat  $input_g/$lane1  $input_g/$lane2   > $output_g/$temp.fastq\n";
    }
    if ($twoLanes[$i] =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])_?([1-2]?).fastq$/) {
        system(  "cp   $input_g/$twoLanes[$i]    $output_g/$twoLanes[$i]" );
    }
    system( "fastqc    --outdir $FastQCdir1   --threads 16    --kmers 7     $input_g/$twoLanes[$i]       >> $FastQCdir1/$twoLanes[$i].runLog   2>&1" );  
}

##system( "rm    $input_g/*.fastq" );  






print "\n\n\n\n\n##################################################################################################\n";
print "\n\nDetecting single-end and paired-end FASTQ files......\n\n";
opendir(my $DH_output, $output_g) || die;     
my @outputFiles = readdir($DH_output);
my @singleEnd = ();
my @pairedEnd = ();
open(seqFiles_FH, ">", "$output_g/singleEnd-pairedEnd-Files.txt")  or  die; 
for ( my $i=0; $i<=$#outputFiles; $i++ ) {     
    next unless $outputFiles[$i] =~ m/\.fastq$/;
    next unless $outputFiles[$i] !~ m/^[.]/;
    next unless $outputFiles[$i] !~ m/[~]$/;
    next unless $outputFiles[$i] !~ m/^unpaired/; 
    $outputFiles[$i] =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])_?([1-2]?)(_Lane[1-2])?\.fastq$/   or  die;
    if ($outputFiles[$i] =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])\.fastq$/) {   ## sinlge end sequencing files.
        $outputFiles[$i] =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])\.fastq$/  or  die;  
        $singleEnd[$#singleEnd+1] =  $outputFiles[$i];
        print  "\n\n        Single-end sequencing files: $outputFiles[$i]\n";
        print seqFiles_FH  "Single-end sequencing files: $outputFiles[$i]\n";
    }else{     ## paired end sequencing files.
        $outputFiles[$i] =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])_([1-2])\.fastq$/  or  die; 
        if ($outputFiles[$i] =~ m/^((\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9]))_1\.fastq$/) { ## The two files of one paired sequencing sample are always side by side. 
            my $temp = $1;
            my $end1 = $temp."_1.fastq";
            my $end2 = $temp."_2.fastq";
            (-e  "$output_g/$end1")  or die;  
            (-e  "$output_g/$end2")  or die;
            $pairedEnd[$#pairedEnd+1] =  $end1;
            $pairedEnd[$#pairedEnd+1] =  $end2;
            print  "\n\n        Paired-end sequencing files: $end1,  $end2\n";
            print seqFiles_FH  "Paired-end sequencing files: $end1,  $end2\n";
        }
    }
}
( ($#pairedEnd+1)%2 == 0 )  or die;
print   seqFiles_FH  "\n\n\n\n\n";
print   seqFiles_FH  "All single-end sequencing files:@singleEnd\n\n\n\n\n\n";
print   seqFiles_FH  "All paired-end sequencing files:@pairedEnd\n\n\n\n\n\n";
print    "\n\n";
print    "\n\n        All single-end sequencing files:@singleEnd\n\n";
print    "\n\n        All paired-end sequencing files:@pairedEnd\n\n";
my $numSingle = $#singleEnd + 1;
my $numPaired = $#pairedEnd + 1;
print seqFiles_FH   "\nThere are $numSingle single-end sequencing files.\n";
print seqFiles_FH   "\nThere are $numPaired paired-end sequencing files.\n";
print     "\n\n        There are $numSingle single-end sequencing files.\n";
print     "\n\n        There are $numPaired paired-end sequencing files.\n";















my $FastQCdir       = "$output_g/FastQC";
my $FastQCdir_10mer = "$output_g/FastQC_10mer";
my $NGSQCToolkit    = "$output_g/NGSQCToolkit";
my $NGSQCToolPaired = "$output_g/NGSQCToolkit_PairedEnd";
my $FASTXtoolkit    = "$output_g/FASTXtoolkit";
if ( !( -e $FastQCdir)       )   { mkdir $FastQCdir        ||  die; }
if ( !( -e $FastQCdir_10mer) )   { mkdir $FastQCdir_10mer  ||  die; }
if ( !( -e $NGSQCToolkit)    )   { mkdir $NGSQCToolkit     ||  die; }
if ( !( -e $NGSQCToolPaired) )   { mkdir $NGSQCToolPaired  ||  die; }
if ( !( -e $FASTXtoolkit)    )   { mkdir $FASTXtoolkit     ||  die; }
  

print "\n\n\n\n\n##################################################################################################\n";
print "\nDetecting the quality of single-end FASTQ files by using FastQC......\n\n";
for ( my $i=0; $i<=$#singleEnd; $i++ ) {     
    my $temp = $singleEnd[$i]; 
    $temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])\.fastq$/   or  die;
    $temp =~ s/\.fastq$//  ||  die;
    system( "fastqc    --outdir $FastQCdir   --threads 16    --kmers 7     $output_g/$temp.fastq       >> $FastQCdir/$temp.runLog   2>&1" );  
} 





print "\n\nDetecting the quality of paired-end FASTQ files by using FastQC and NGS_QC_Toolkit......\n\n";
for ( my $j=0; $j<=$#pairedEnd; $j=$j+2 ) {     
    my $temp1 = $pairedEnd[$j]; 
    my $temp2 = $pairedEnd[$j+1]; 
    $temp1 =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])_1\.fastq$/   or  die;
    $temp2 =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])_2\.fastq$/   or  die;
    $temp1 =~ s/\.fastq$//  ||  die;
    $temp2 =~ s/\.fastq$//  ||  die;
    system( "fastqc    --outdir $FastQCdir   --threads 16    --kmers 7     $output_g/$temp1.fastq       >> $FastQCdir/$temp1.runLog   2>&1" );  
    system( "fastqc    --outdir $FastQCdir   --threads 16    --kmers 7     $output_g/$temp2.fastq       >> $FastQCdir/$temp2.runLog   2>&1" );  
    print  seqFiles_FH  "\n\nquality statistics: $temp1,  $temp2\n";
    my $temp = $temp1;
    $temp =~ s/_1$//  ||  die;
    if ( !(-e "$NGSQCToolPaired/$temp") )   { mkdir  "$NGSQCToolPaired/$temp"  ||  die; }
    system( "perl  $IlluQC_PRLL     -pe $output_g/$temp1.fastq   $output_g/$temp2.fastq   N  A    -cpus 16     -onlyStat    -outputFolder $NGSQCToolPaired/$temp    >> $NGSQCToolPaired/$temp/$temp.runLog  2>&1" ); 
} 

    


print "\n\nDetecting the quality of all FASTQ files by using FastQC, NGS_QC_Toolkit and FASTX-Toolkit......\n\n";
for ( my $i=0; $i<=$#outputFiles; $i++ ) {     
    next unless $outputFiles[$i] =~ m/\.fastq$/;
    next unless $outputFiles[$i] !~ m/^[.]/;
    next unless $outputFiles[$i] !~ m/[~]$/;
    my $temp = $outputFiles[$i]; 
    $temp =~ m/^(\d{2})_($pattern)_($pattern)_($pattern)_($pattern)_($pattern)_(Rep[1-9])_?([1-2]?)\.fastq$/   or  die;
    $temp =~ s/\.fastq$//  ||  die;
    system( "fastqc    --outdir $FastQCdir_10mer    --threads 16    --kmers 10    $output_g/$temp.fastq       >> $FastQCdir_10mer/$temp.runLog   2>&1" );  
    if ( !(-e "$NGSQCToolkit/$temp") )   { mkdir  "$NGSQCToolkit/$temp"  ||  die; }
    system( "perl   $IlluQC_PRLL   -se $output_g/$temp.fastq   N  A      -cpus 16   -onlyStat    -outputFolder $NGSQCToolkit/$temp         >> $NGSQCToolkit/$temp/$temp.runLog   2>&1" );      
    system( "fastx_quality_stats                        -i $output_g/$temp.fastq                    -o $FASTXtoolkit/$temp.txt            >> $FASTXtoolkit/$temp.runLog   2>&1" ); 
    system( "fastq_quality_boxplot_graph.sh             -i $FASTXtoolkit/$temp.txt     -t $temp     -o $FASTXtoolkit/$temp.quality.png    >> $FASTXtoolkit/$temp.runLog   2>&1" ); 
    system( "fastx_nucleotide_distribution_graph.sh     -i $FASTXtoolkit/$temp.txt     -t $temp     -o $FASTXtoolkit/$temp.nucDis.png     >> $FASTXtoolkit/$temp.runLog   2>&1" ); 
}





print "\n\n        Job Done! Cheers! \n\n\n\n\n";































