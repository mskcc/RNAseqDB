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
push @usage, "  -d | --sample-id  Sample ID, i.e. analysis_id for TCGA data or SRA accession for GTEx\n";
push @usage, "                    The sample ID must be included in the cart file\n";
push @usage, "  -s | --submit     Submit a job to the cluster to download data if specified\n";
push @usage, "Example:\n";
push @usage, "  perl download.pl -f manifest.xml -s\n";


my ( $help, $config_file, $cart_file, $sample_id );
my ( $submit ) = 0;

GetOptions
(
'h|help|?'      => \$help,
'c|config=s'    => \$config_file,
'f|file=s'      => \$cart_file,
'd|sample-id=s' => \$sample_id,
's|submit'      => \$submit,
);

if ( $help ) {
    print @usage;
    exit(0);
}

if( !defined $cart_file or !-s $cart_file ) {
    print "ERROR: Must provide a cart file containing the information about the data for download\n";
    print @usage;
    exit(0);
}

if (defined $sample_id){
    my $str = `grep $sample_id $cart_file`;
    chomp $str;
    (defined $str) or die "ERROR: $sample_id is not in file $cart_file\n";

    if ( $cart_file =~ /^cart/ ){
        ($sample_id =~ /^SRR/) or die "ERROR: $sample_id does not look like SRA accession number\n";
    }else{
        ($sample_id =~ /.*-.*-/) or die "ERROR: Illegal TCGA analysis_id $sample_id\n";
    }
}

############################# Read configuration file ############################
#
#(defined $config_file) or $config_file = "$FindBin::Bin/config.txt";
#( -e $config_file ) or die "ERROR: The configuration file $config_file does not exist\n";
#my ( %config, $gtex_path, $tcga_path );
#map{ chomp; /^\s*([^=\s]+)\s*=\s*(.*)$/; $config{$1} = $2 if (defined $1 and defined $2) } `egrep -v \"^#\" $config_file`;
#$gtex_path  =  $config{ gtex_path };
#$tcga_path  =  $config{ tcga_path };
#
#$gtex_path = abs_path($gtex_path);
#$tcga_path = abs_path($tcga_path);
#
#my $gtex_dir_len = length($gtex_path);
#my $tcga_dir_len = length($tcga_path);

#my $work_dir = getcwd;
#chdir $work_dir;

############################# Submit job to download data ############################

if ( $cart_file =~ /^cart/ ) { #or $gtex_path eq substr($work_dir, 0, $gtex_dir_len)) {
    
    if (defined $sample_id){
        my $cmd = "prefetch $sample_id";
        system($cmd);
    }else{
        my $cmd = "prefetch $cart_file";
        if( $submit ){
            my $ret = `perl $FindBin::Bin/qsub.pl -s \047$cmd\047 -m 10 -p 1 -t 48`;
            print "$cart_file\t$ret\n";
        }else{
            system($cmd);
        }
    }
}else{ #or $tcga_path eq substr($work_dir, 0, $tcga_dir_len)) {
    
    (-e $ENV{"HOME"}.'/.ssh/cghub.key') or die "ERROR: Missing file ".$ENV{"HOME"}."/.ssh/cghub.key\n";
    
    if ( defined $sample_id or (!$submit and $cart_file =~ /^manifest/ and $cart_file =~ /\.xml$/ )){
        
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
            next if(!$download or (defined $sample_id and $sample_id ne $analysis_id));
            
            print "$analysis_id\n";
            my $hit_fh = IO::File->new( 'temp.xml', ">" );
            $hit_fh->print("<ResultSet>\n <Hits>1</Hits>\n" . $hit . "</ResultSet>");
            $hit_fh->close;
            
            `gtdownload -c $ENV{"HOME"}/.ssh/cghub.key -d temp.xml`;
            `rm temp.xml`;
        }
        $cart_fh->close;
       
        
    }else{
        
        my $cmd = "gtdownload -c ".$ENV{"HOME"}."/.ssh/cghub.key -d $cart_file";
        if( $submit ){
            my $ret = `perl $FindBin::Bin/qsub.pl -s \047$cmd\047 -m 10 -p 1 -t 48`;
            print "$cart_file\t$ret\n";
        }else{
            system($cmd);
        }
    }
    
}