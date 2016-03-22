#!/usr/bin/perl -w
###################################################################################
#
#  Description: Process FASTQ files.
#
###################################################################################

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use Cwd;
use FindBin;
use lib "$FindBin::Bin";
use Mosaic;
use File::Basename;
use List::Util qw(min max);

my @usage;
push @usage, "\nUsage: process-fastq.pl file1 [ file2 ] [ -c config.txt -o directory ] \n\n";
push @usage, "Options:\n";
push @usage, "  file1 [file2]  Read files (in FASTQ format) <required>.\n";
push @usage, "  -c, --config   Configuration file [will search current/home dir if not given].\n";
push @usage, "  -h, --help     Displays this information.\n";
push @usage, "  -o, --output   An output directory, default is current working directory.\n\n";


my $help;
my $config_file;
my $fastq1;
my $fastq2;
my $output_dir;


GetOptions
(
 'h|help|?' => \$help,
 'config=s' => \$config_file,
 'output=s' => \$output_dir,
);

if ($help) {
   print @usage;
   exit(0);
}

if (defined $output_dir) {
    if (!-e $output_dir || !-d $output_dir){
    	print "\nThe output directory $output_dir does not exist!\n\n";
        print @usage;
        exit;
    }
}else{
    $output_dir = getcwd;
}

($fastq1, $fastq2) = @ARGV;
if (defined $fastq1) {
    if (!-e $fastq1){
        print "\nFastq file $fastq1 does not exist!\n\n";
        print @usage;
        exit;
    }
}else{
    print "Please provide at least a fastq file!\n\n";
    print @usage;
    exit;
}

if (defined $fastq2 && !-e $fastq2) {
    print "\nFastq file $fastq2 does not exist!\n\n";
    print @usage;
    exit;
}

###################################################################################
# The following code processes configuration file

if (!defined $config_file) {
    $config_file = "$FindBin::Bin/config.txt";
}
( -e $config_file ) or die "ERROR: The configuration file $config_file does not exist\n";

my %config;
map{ chomp; /^\s*([^=\s]+)\s*=\s*(.*)$/; $config{$1} = $2 if (defined $1 && defined $2) } `egrep -v \"^#\" $config_file`;

( exists $config{ fastqc_bin } && exists $config{ trimmomatic_jar } ) || die "fastqc or trimmomatic not defined in configuration file!\n\n";


my $thread_n = 1;
$thread_n           = $config{ thread_no } if ( exists $config{ thread_no } );
my $fastqc_bin      = $config{ fastqc_bin };
my $trimmomatic_bin = 'java -jar ' . $config{ trimmomatic_jar } ;
my $trimmomatic_dir = $config{ trimmomatic_dir };


###################################################################################
# Run FASTQC if necessary

die "Can not find file: $fastqc_bin!\n\n" if (!-e $fastqc_bin);

print STDERR "Running Fastqc...\n";

`mkdir -p $output_dir/fastqc`;

my ($name1, $path1) = fileparse($fastq1);
my ($name2, $path2) = fileparse($fastq2) if (defined $fastq2);

my $idx1 = index($name1,'.');
my $idx2 = index($name2,'.') if (defined $fastq2);

my $name1_pre = $name1;
my $name2_pre = $name2 if (defined $fastq2);

$name1_pre = substr($name1, 0, $idx1) if ($idx1 > 0);
$name2_pre = substr($name2, 0, $idx2) if (defined $fastq2 && $idx2 > 0);


if (defined $fastq2){
    `$fastqc_bin --extract -outdir=$output_dir/fastqc $fastq1 $fastq2` if (!-e "$output_dir/fastqc/$name2_pre".'_fastqc');
}else{
    `$fastqc_bin --extract -outdir=$output_dir/fastqc $fastq1` if (!-e "$output_dir/fastqc/$name1_pre".'_fastqc');
}


# Parse FASTQC results
my @qc_report = `find $output_dir/fastqc -name summary.txt | xargs -i% grep FAIL % | cut -f 2`;

my $btrim_base = 0;
my $btrim_adap = 0;

foreach(@qc_report){
    chomp;
    if (/Sequence Duplication Levels/){
        print STDERR "\nWarning: very high sequence duplication level!\n";
    }elsif (/Per base sequence quality/ || /Per base N content/){
        $btrim_base = 1;
    }elsif (/Overrepresented sequences/) {
        $btrim_adap = 1;
    }
}

###################################################################################
# Run trimmomatic if necessary

if ($btrim_base == 0 && $btrim_adap == 0){
    print STDERR "\nNo need to trim fastq file.\n\nDone!\n\n";
    exit(0);
}

print "Detect FASTQ format...\n";
my $phred = detect_qual_format($fastq1);

print STDERR "Running Trimmomatic...\n";

if (!-e "$path1/trim/$name1" || (defined $fastq2 && !-e "$path2/trim/$name2")){

    `mkdir -p $path1/trim`;
    `mkdir -p $path2/trim` if (defined($fastq2) && $path1 ne $path2);
    
    my $pairing = 'SE';
    $pairing = 'PE' if (defined $fastq2);
    
    my $cmd_str;
    if ($btrim_adap == 1){
        print STDERR "Based on FastQC report, we will do adaptor trimming.\n";
        $cmd_str = 'ILLUMINACLIP:'.$trimmomatic_dir.'/adapters/Combine-'.$pairing.'.fa:2:30:10';
    }
        
    if ($btrim_base == 1){
        print STDERR "Based on FastQC report, we will trim low quality bases.\n";
        if ($phred eq 'phred64'){
            $cmd_str .= ' LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:36';
        }else{
            $cmd_str .= ' LEADING:3  TRAILING:3  SLIDINGWINDOW:4:15 MINLEN:36';
        }
    }

    if (defined $fastq2){
        `$trimmomatic_bin PE -threads $thread_n -$phred $fastq1 $fastq2 $path1/trim/$name1 $path1/trim/unpaired_$name1 $path2/trim/$name2 $path2/trim/unpaired_$name2 $cmd_str`;
    }else{
        `$trimmomatic_bin SE -threads $thread_n -$phred $fastq1 $path1/trim/$name1 $cmd_str`;
    }
}

print STDERR "Done!\n\n";
