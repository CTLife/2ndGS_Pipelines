#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.22;
## Perl5 version >= 5.22
## You can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "5_rmDups/1_BWA"
my $output_g = '';  ## such as "7_BW/1_BWA"

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------ 

        Usage:
               perl  Bulk_5hmCSeal_7.pl    [-version]    [-help]    [-in inputDir]    [-out outDir]
        For instance:
               nohup time  perl  Bulk_5hmCSeal_7.pl    -in 5_rmDups/1_BWA   -out 7_BW/1_BWA    > Bulk_5hmCSeal_7.runLog  2>&1   &

        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.

        -help           Show this help message and exit.

        Required arguments:
        -in inputDir        "inputDir" is the name of input path that contains your BAM files.  (no default)

        -out outDir         "outDir" is the name of output path that contains your running results of this step.  (no default)
        -----------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines, please visit https://github.com/CTLife/2ndGS_Pipelines
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';

## Version Infromation
my $version = " version 0.1,  2022-09-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$input_g  = '5_rmDups/1_BWA';  ## This is only an initialization value or suggesting value, not default value.
$output_g = '7_BW/1_BWA';      ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help    -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  Bulk_5hmCSeal_7.pl  -help' \n";
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
 
&myMakeDir($output_g);

opendir(my $DH_input_g, $input_g)  ||  die;
my @inputFiles_g = readdir($DH_input_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Checking all the necessary softwares in this step......" ;
sub printVersion  {
    my $software = $_[0];
    system("echo    '##############################################################################'  >> $output_g/VersionsOfSoftwares.txt   2>&1");
    system("echo    '#########$software'                                                              >> $output_g/VersionsOfSoftwares.txt   2>&1");
    system("$software                                                                                 >> $output_g/VersionsOfSoftwares.txt   2>&1");
    system("echo    '\n\n\n\n\n\n'                                                                    >> $output_g/VersionsOfSoftwares.txt   2>&1");
}
&printVersion("deeptools --version");
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
say   "\n\n\n\n\n\n##################################################################################################";
say   "Using bamCoverage ......";

for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {
        next unless $inputFiles_g[$i] =~ m/\.bam$/;
        next unless $inputFiles_g[$i] !~ m/^[.]/;
        next unless $inputFiles_g[$i] !~ m/[~]$/;
        next unless $inputFiles_g[$i] !~ m/.bai$/;
        say   "\t......$inputFiles_g[$i]" ;
        my $temp = $inputFiles_g[$i];
        $temp =~ s/.bam$//  or  die;
        system(" bamCoverage   --bam $input_g/$temp.bam     --outFileName $output_g/$temp.bw     --outFileFormat bigwig    --extendReads   --binSize 10    -p 16     --verbose   --normalizeUsing RPKM   --ignoreForNormalization chrX chrY chrM         > $output_g/$temp.runLog.txt  2>&1  ");                                                  
}

}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";


 

  
## END







