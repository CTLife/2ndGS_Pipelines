#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.22;
## Perl5 version >= 5.22
## You can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "2B_Bismark_unmapped_SE", global variable.   
my $output_g = '';  ## such as "2B_Bismark_unmapped_SE_Summary", global variable. 

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------                     
        Usage:
               perl bismark_summary_SE.pl    [-version]    [-help]     [-in inputDir]    [-out outDir]
        For instance:
               perl bismark_summary_SE.pl    -in 2B_Bismark_unmapped_SE   -out 2B_Bismark_unmapped_SE_Summary    >bismark_summary_SE.runLog  

        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.

        -help           Show this help message and exit.

        Required arguments:
        -in inputDir        "inputDir" is the name of input path.  (no default)

        -out outDir         "outDir" is the name of output path that contains your running results of this step.  (no default)
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';

## Version Infromation
my $version = "    bismark_summary_SE, version 0.9.4,  2018-02-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$input_g  = '2B_Bismark_unmapped_SE';     ## This is only an initialization value or suggesting value, not default value.
$output_g = '2B_Bismark_unmapped_SE_Summary';     ## This is only an initialization value or suggesting value, not default value.

## Available Arguments 
my $available = "   -version    -help    -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl bismark_summary_SE.pl  -help' \n";
    exit 0;
}

## Get Arguments 
if ( exists $args{'-version' }   )     { say  "\n$version\n";    exit 0; }
if ( exists $args{'-help'    }   )     { say  "\n$HELP\n";       exit 0; }
if ( exists $args{'-in'      }   )     { $input_g  = $args{'-in'      }; }else{say   "\n -in     is required.\n";   say  "\n$HELP\n";    exit 0; }
if ( exists $args{'-out'     }   )     { $output_g = $args{'-out'     }; }else{say   "\n -out    is required.\n";   say  "\n$HELP\n";    exit 0; }

## Conditions
$input_g  =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";
$output_g =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";

## Print Command Arguments to Standard Output
say  "\n
        ################ Arguments ###############################
                Input   Path:  $input_g
                Output  Path:  $output_g
        ##########################################################
\n";
}
###################################################################################################################################################################################################




sub myMakeDir  {
    my $path = $_[0];
    if ( !( -e $path) )  { system("mkdir  -p  $path"); }
    if ( !( -e $path) )  { mkdir $path  ||  die;       }
}

my $pattern_g  = "[-.0-9A-Za-z]+";
my $numCores_g = 4;
&myMakeDir($output_g);
opendir(my $DH_input_g, $input_g)  ||  die;
my @inputFiles_g = readdir($DH_input_g);

open(FH1, ">", "$output_g/bismark_summary_SE.txt" ) or die "$!"; 
print  FH1  "Sample_Name\tTotal_Sequence_Pairs \tUnique_Best_Hit\tUnique_Best_Hit_Ratio \tNo_Alignments\tNo_Alignments_Ratio\t" .      ##6
            "Multi_Mapped\tMulti_Mapped_Ratio\t not_Extracted\tnot_Extracted_Ratio\t Total_Fragments  \tTotal_Cs\t" .                  ##6                   
            "mCpG\tmCHG\tmCHH\tmC_unknown\tmCH\tmC_all \tun_mCpG\tun_mCHG\tun_mCHH\tun_mC_unknown\tun_mCH\tun_mC_all\t" .              ##12
            "mCpG_Ratio\tmCHG_Ratio\tmCHH_Ratio\tmC_unknown_Ratio\tmCH_Ratio\tmC_all_Ratio\t" .                                        ##6
            "methylated_C_in_CpG_context\tmethylated_C_in_CHG_context\tmethylated_C_in_CHH_context\t" .                                ##3
            "methylated_C_in_unknown_context \tSum_mC_unmC \tSum_mC_unmC(CG+CH+unknown)\tTotal_mC_Ratio(CG+CH+unknown)\t" .            ##4
            "Total_unmC_Ratio(CG+CH+unknown)\n";

my $numFiles = 0;
for ( my $i=0; $i<=$#inputFiles_g; $i++ ) { 
        next unless $inputFiles_g[$i] =~ m/^(\d+)_($pattern_g)_Rep\d+_[12]+_SE_report.txt$/;     
        next unless $inputFiles_g[$i] !~ m/^[.]/;
        next unless $inputFiles_g[$i] !~ m/[~]$/;
        next unless $inputFiles_g[$i] !~ m/\.runLog$/;
        next unless $inputFiles_g[$i] !~ m/\.bam/;
        $numFiles = $numFiles+1;
        my $file1 = "$input_g/$inputFiles_g[$i]";    
        open(tempFH1, "<",   $file1)   or die;
        my @lines1 = <tempFH1>; 

        my $Sample_1               = $inputFiles_g[$i];
        $Sample_1 =~ s/_SE_report.txt$// or die;
        my $TotalSeqPairs_2        = "NA";
        my $UniqueBestHit_3        = "NA";
        my $UniqueBestHit_Ratio_4  = "NA";
        my $NoAlignments_5         = "NA";
        my $NoAlignments_Ratio_6   = "NA";
        my $MultiMapped_7          = "NA";
        my $MultiMapped_Ratio_8    = "NA";
        my $notExtracted_8a        = "NA";
        my $notExtracted_Ratio_8b  = "NA";
        my $TotalFragments_8c      = "NA";
        my $TotalCs_9              = "NA";

        my $mCpG_10                = "NA";
        my $mCHG_11                = "NA";
        my $mCHH_12                = "NA";
        my $mCunknown_13           = "NA";
        my $mCH_14                 = "NA";
        my $mCall_15               = "NA";

        my $unmCpG_16              = "NA";
        my $unmCHG_17              = "NA";
        my $unmCHH_18              = "NA";
        my $unmCunknown_19         = "NA";
        my $unmCH_20               = "NA";
        my $unmCall_21             = "NA";

        my $mCpG_Ratio_22          = "NA";
        my $mCHG_Ratio_23          = "NA";
        my $mCHH_Ratio_24          = "NA";
        my $mCunknown_Ratio_25     = "NA";
        my $mCH_Ratio_26           = "NA";
        my $mCall_Ratio_27         = "NA";

        my $mCpG_Ratio2_28         = "NA";
        my $mCHG_Ratio2_29         = "NA";
        my $mCHH_Ratio2_30         = "NA";
        my $mCunknown_Ratio2_31    = "NA";

        my $Sum_mC_unmC_32         = "NA";
        my $Sum_all_33             = "NA";
        my $Ratio_all_34           = "NA";
        my $Ratio_unmC_35          = "NA";

        for ( my $j=0; $j<=$#lines1; $j++ ) { 
            if($lines1[$j] =~ m/^Sequences\s+analysed\s+in\s+total:\s+(\d+)\n$/i) {$TotalSeqPairs_2 = $1;}
            if($lines1[$j] =~ m/^Number\s+of\s+alignments\s+with\s+a\s+unique\s+best\s+hit\s+from\s+the\s+different\s+alignments:\s+(\d+)\n$/i) {$UniqueBestHit_3 = $1;}
            if($lines1[$j] =~ m/^Sequences\s+with\s+no\s+alignments\s+under\s+any\s+condition:\s+(\d+)\n$/i) {$NoAlignments_5 = $1;}
            if($lines1[$j] =~ m/^Sequences\s+did\s+not\s+map\s+uniquely:\s+(\d+)\n$/i) {$MultiMapped_7 = $1;}
            if($lines1[$j] =~ m/^Sequences\s+which\s+were\s+discarded\s+because\s+genomic\s+sequence\s+could\s+not\s+be\s+extracted:\s+(\d+)\n$/i) {$notExtracted_8a = $1;}
            if($lines1[$j] =~ m/^Total\s+number\s+of\s+C\'s\s+analysed:\s+(\d+)\n$/i) {$TotalCs_9 = $1;}
            if($lines1[$j] =~ m/^Total\s+methylated\s+C\'s\s+in\s+CpG\s+context:\s+(\d+)\n$/i)     {$mCpG_10 = $1;}
            if($lines1[$j] =~ m/^Total\s+methylated\s+C\'s\s+in\s+CHG\s+context:\s+(\d+)\n$/i)     {$mCHG_11 = $1;}
            if($lines1[$j] =~ m/^Total\s+methylated\s+C\'s\s+in\s+CHH\s+context:\s+(\d+)\n$/i)     {$mCHH_12 = $1;}
            if($lines1[$j] =~ m/^Total\s+methylated\s+C\'s\s+in\s+Unknown\s+context:\s+(\d+)\n$/i) {$mCunknown_13 = $1;}
            if($lines1[$j] =~ m/^Total\s+unmethylated\s+C\'s\s+in\s+CpG\s+context:\s+(\d+)\n$/i)       {$unmCpG_16 = $1;}
            if($lines1[$j] =~ m/^Total\s+unmethylated\s+C\'s\s+in\s+CHG\s+context:\s+(\d+)\n$/i)       {$unmCHG_17 = $1;}
            if($lines1[$j] =~ m/^Total\s+unmethylated\s+C\'s\s+in\s+CHH\s+context:\s+(\d+)\n$/i)       {$unmCHH_18 = $1;}
            if($lines1[$j] =~ m/^Total\s+unmethylated\s+C\'s\s+in\s+Unknown\s+context:\s+(\d+)\n$/i)   {$unmCunknown_19 = $1;}
            if($lines1[$j] =~ m/^C\s+methylated\s+in\s+CpG\s+context:\s+([\d\.]+)\%\n$/i)                         {$mCpG_Ratio2_28 = $1;}
            if($lines1[$j] =~ m/^C\s+methylated\s+in\s+CHG\s+context:\s+([\d\.]+)\%\n$/i)                         {$mCHG_Ratio2_29 = $1;}
            if($lines1[$j] =~ m/^C\s+methylated\s+in\s+CHH\s+context:\s+([\d\.]+)\%\n$/i)                         {$mCHH_Ratio2_30 = $1;}
            if($lines1[$j] =~ m/^C\s+methylated\s+in\s+unknown\s+context\s+\(CN\s+or\s+CHN\):\s+([\d\.]+)\%\n$/i) {$mCunknown_Ratio2_31 = $1;}
        }
        
       $UniqueBestHit_Ratio_4 = $UniqueBestHit_3/$TotalSeqPairs_2;   
       $NoAlignments_Ratio_6  = $NoAlignments_5/$TotalSeqPairs_2;
       $MultiMapped_Ratio_8   = $MultiMapped_7/$TotalSeqPairs_2;
       $notExtracted_Ratio_8b = $notExtracted_8a/$TotalSeqPairs_2;
       $TotalFragments_8c     = $UniqueBestHit_3 + $NoAlignments_5 + $MultiMapped_7 + $notExtracted_8a;

       $mCH_14     = $mCHG_11 + $mCHH_12;
       $mCall_15   = $mCH_14 + $mCpG_10;      ## unknown context is excluded.
       $unmCH_20   = $unmCHG_17 + $unmCHH_18;  
       $unmCall_21 = $unmCH_20 + $unmCpG_16;  ## unknown context is excluded.

       $mCpG_Ratio_22 = $mCpG_10/($mCpG_10+$unmCpG_16);
       $mCHG_Ratio_23 = $mCHG_11/($mCHG_11+$unmCHG_17);
       $mCHH_Ratio_24 = $mCHH_12/($mCHH_12+$unmCHH_18);
       $mCunknown_Ratio_25 = $mCunknown_13/($mCunknown_13+$unmCunknown_19);
       $mCH_Ratio_26   = $mCH_14/($mCH_14+$unmCH_20);
       $mCall_Ratio_27 = $mCall_15/($mCall_15+$unmCall_21); ## unknown context is excluded.

       $Sum_mC_unmC_32 = $mCall_15+$unmCall_21;
       $Sum_all_33     = $mCall_15+$unmCall_21+$mCunknown_13+$unmCunknown_19; 
       $Ratio_all_34   = ($mCall_15+$mCunknown_13)/$Sum_all_33;
       $Ratio_unmC_35  = ($unmCall_21+$unmCunknown_19)/$Sum_all_33;

       $mCpG_Ratio2_28      = $mCpG_Ratio2_28/100;
       $mCHG_Ratio2_29      = $mCHG_Ratio2_29/100;
       $mCHH_Ratio2_30      = $mCHH_Ratio2_30/100;
       $mCunknown_Ratio2_31 = $mCunknown_Ratio2_31/100;

       print    FH1  "$Sample_1\t$TotalSeqPairs_2 \t$UniqueBestHit_3\t$UniqueBestHit_Ratio_4 \t$NoAlignments_5\t$NoAlignments_Ratio_6\t" .
                     "$MultiMapped_7\t$MultiMapped_Ratio_8 \t$notExtracted_8a\t$notExtracted_Ratio_8b \t$TotalFragments_8c \t$TotalCs_9\t" .
                     "$mCpG_10\t$mCHG_11\t$mCHH_12\t$mCunknown_13\t$mCH_14\t$mCall_15\t" . 
                     "$unmCpG_16\t$unmCHG_17\t$unmCHH_18\t$unmCunknown_19\t$unmCH_20\t$unmCall_21\t" . 
                     "$mCpG_Ratio_22\t$mCHG_Ratio_23\t$mCHH_Ratio_24\t$mCunknown_Ratio_25\t$mCH_Ratio_26\t$mCall_Ratio_27\t" .
                     "$mCpG_Ratio2_28\t$mCHG_Ratio2_29\t$mCHH_Ratio2_30\t$mCunknown_Ratio2_31\t$Sum_mC_unmC_32\t$Sum_all_33\t$Ratio_all_34\t$Ratio_unmC_35\n";
}

print "\n\nNumber of samples: $numFiles\n\n\n";


















