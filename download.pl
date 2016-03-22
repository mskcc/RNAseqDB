#!/usr/bin/perl -w

use strict;
use warnings FATAL => 'all';
use File::Basename;
use IO::File;
use Getopt::Std;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin";
use Cwd;
use Cwd 'abs_path';

my @usage;
push @usage, "\nUsage:  download.pl <-f file> [options]\n\n";
push @usage, "Options:\n";
push @usage, "  -h | --help       Displays this information.\n";
push @usage, "  -c | --config     A configuration file, default is config.txt in same dir of this script\n";
push @usage, "  -f | --file       A cart file [for dbGaP data] or manifest.xml/urls.txt file [for CGHub]\n";
push @usage, "  -s | --submit     Submit a job for data download if specified a value 'yes'; default: no\n";
push @usage, "  -d | --directory  A directory to store data. Working directory is used if not specified\n";
push @usage, "  -o | --overwrite  Overwrite previous download by default or if specified 'yes'\n";
push @usage, "Example:\n";
push @usage, "  perl download.pl -f manifest.xml -s yes\n";
push @usage, "  perl download.pl -f urls.txt\n\n";


my ( $help, $config_file, $cart_file, $work_dir, $submit, $overwrite);

GetOptions
(
'h|help|?'      => \$help,
'c|config=s'    => \$config_file,
'f|file=s'      => \$cart_file,
's|submit=s'    => \$submit,
'd|directory=s' => \$work_dir,
'o|overwrite=s' => \$overwrite,
);

if ( $help ) {
    print @usage;
    exit(0);
}

if( !defined $cart_file ) {
    print "ERROR: Must provide a file containing information of data to be downloaded\n";
    print @usage;
    exit(0);
}

if (!defined $work_dir) {
    $work_dir = getcwd;
}else{
    ( -e $work_dir ) or die "ERROR: The directory $work_dir does not exist\n";
    $work_dir = abs_path($work_dir);
}

if(defined $submit and (lc($submit) eq 'yes' or lc($submit) eq 'y')){
    $submit = 'yes' ;
}else{
    $submit = 'no' ;
}

if(defined $overwrite and (lc($overwrite) eq 'no' or lc($overwrite) eq 'n')){
    $overwrite = 'no' ;
    ($cart_file =~ /\.xml$/ and -s $cart_file) or die "ERROR: To avoid overwriting previous download, XML format of input file is required\n";
}else{
    $overwrite = 'yes' ;
}
############################# Read configuration file ############################

(defined $config_file) or $config_file = "$FindBin::Bin/config.txt";
( -e $config_file ) or die "ERROR: The configuration file $config_file does not exist\n";
my ( %config, $gtex_path, $tcga_path );
map{ chomp; /^\s*([^=\s]+)\s*=\s*(.*)$/; $config{$1} = $2 if (defined $1 and defined $2) } `egrep -v \"^#\" $config_file`;
$gtex_path  =  $config{ gtex_path };
$tcga_path  =  $config{ tcga_path };

$gtex_path = abs_path($gtex_path);
$tcga_path = abs_path($tcga_path);

############################# Submit job to download data ############################

my $gtex_dir_len = length($gtex_path);
my $tcga_dir_len = length($tcga_path);

chdir $work_dir;

if ($gtex_path eq substr($work_dir, 0, $gtex_dir_len)) {
    
    my $cmd = "prefetch $cart_file";
    if($submit eq 'yes'){
        my $ret = `perl $FindBin::Bin/qsub.pl -s \047$cmd\047 -m 10 -p 1 -t 48`;
        print "$cart_file\t$ret\n";
    }else{
        system($cmd);
    }
    
} elsif ($tcga_path eq substr($work_dir, 0, $tcga_dir_len)) {
    
    (-e $ENV{"HOME"}.'/.ssh/cghub.key') or die "ERROR: Missing file ".$ENV{"HOME"}."/.ssh/cghub.key\n";
    if ($overwrite eq 'no'){
        
        my $cart_fh = IO::File->new( $cart_file );
        while( my $line = $cart_fh->getline ) {
            next if( $line !~ /\<Result id=/);
            my $hit = "<Result id=\"1\">\n";
            my $analysis_id;
            my $download = 0;
            while( $line = $cart_fh->getline ) {
                $hit .= $line;
                if ($line =~ /\<analysis_id\>/){
                    $line =~ /\<analysis_id\>(.+)\<\/analysis_id\>/;
                    $analysis_id = $1 if (defined $1);
                }
                if ($line =~ /\<filename\>(.+)\<\/filename\>/){
                    $download = 1 if(!-s "$analysis_id/$1");
                }
                last if ($line =~ /\<\/Result/);
            }
            if($download){
                print "$analysis_id\n";
                my $hit_fh = IO::File->new( 'temp.xml', ">" );
                $hit_fh->print("<ResultSet>\n <Hits>1</Hits>\n" . $hit . "</ResultSet>");
                $hit_fh->close;
                
                `gtdownload -c $ENV{"HOME"}/.ssh/cghub.key -d temp.xml`;
                `rm temp.xml`;
            }
        }
        $cart_fh->close;
        
    }else{

        my $cmd = "gtdownload -c ".$ENV{"HOME"}."/.ssh/cghub.key -d $cart_file";
        if($submit eq 'yes'){
            my $ret = `perl $FindBin::Bin/qsub.pl -s \047$cmd\047 -m 10 -p 1 -t 48`;
            print "$cart_file\t$ret\n";
        }else{
            system($cmd);
        }
        
    }
    
} else {
    
    die "ERROR: Not sure from where to download: dbGaP or CGHub\n";
    
}

