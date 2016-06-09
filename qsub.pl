#!/usr/bin/perl -w

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin";
use Cwd;


my @usage;
push @usage, "\nUsage:  qsub.pl <-s command> [options]\n\n";
push @usage, "Options:\n";
push @usage, "  -h | --help       Displays this information.\n";
push @usage, "  -c | --config     A configuration file, default is config.txt in same dir of this script\n";
push @usage, "  -s | --command    Input command to be run in the cluster\n";
push @usage, "  -m | --memory     Specifies memory [gb]; Will read configuration file if not specified\n";
push @usage, "  -p | --cpu        Specifies number of CPUs; Will read configuration file if not specified\n";
push @usage, "  -t | --time       Specifies time [hours]; Will read configuration file if not specified\n";
push @usage, "  -e | --error      Appends the standard error output of the job to the specified file\n";
push @usage, "  -o | --output     Appends the standard output of the job to the specified file\n";
push @usage, "Example:\n";
push @usage, "  perl qsub.pl -s 'perl ~/scripts/calc-expression.pl -i SRR1069166.sra'\n";
push @usage, "  perl qsub.pl -s 'gtdownload -c ~/.ssh/cghub.key -d urls.txt' -m 10 -p 1 -t 48\n\n";


my ( $help, $config_file, $input_cmd, $thread_n, $memory, $wall_time, $error_file, $output_file );

GetOptions
(
 'h|help|?'     => \$help,
 'c|config=s'   => \$config_file,
 's|command=s'  => \$input_cmd,
 'p|cpu=s'      => \$thread_n,
 'm|memory=s'   => \$memory,
 't|time=s'     => \$wall_time,
 'e|error=s'    => \$error_file,
 'o|output=s'   => \$output_file,
);

if ( $help ) {
    print @usage;
    exit(0);
}

if( !defined $input_cmd ) {
    print "ERROR: An command must be provided\n";
    print @usage;
    exit(0);
}

######################### Processes configuration file #################################

(defined $config_file) or $config_file = "$FindBin::Bin/config.txt";
( -e $config_file ) or die "ERROR: The configuration file $config_file does not exist\n";

# Read configuration file
my ( %config, $cluster );
map{ chomp; /^\s*([^=\s]+)\s*=\s*(.*)$/; $config{$1} = $2 if (defined $1 and defined $2) } `egrep -v \"^#\" $config_file`;

if ( !defined $thread_n ){
    if ( defined $config{ thread_no } ){
        $thread_n = $config{ thread_no };
    }else{
        $thread_n = 1;
    }
}
$memory    = $config{ memory }    if ( !defined $memory );
$wall_time = $config{ wall_time } if ( !defined $wall_time );

if ( exists $config{ cluster } ) {
    $cluster = $config{ cluster }
}else{
    die "ERROR: variable 'cluster' is not defined in the configuration file\n";
}

######################### Submit jobs to analyze samples #################################


if(!defined $cluster){
    die "ERROR: Variable 'cluster' is not defined in configuration file\n";
}elsif(lc($cluster) eq "luna"){
    RunLunaQsub();
}elsif(lc($cluster) eq "hal"){
    RunHalQsub();
}else{
    die "ERROR: Unknown cluster\n";
}


sub RunLunaQsub {
    my $err_log = '';
    $err_log  = '-e '.$error_file   if( defined $error_file);
    $err_log .= ' -o '.$output_file if( defined $output_file);
    my $ret = `bsub -We $wall_time:00 -n$thread_n -R \047span[hosts=1]\047 -R \047rusage[mem=$memory]\047 $err_log \047$input_cmd\047\n`;
    $ret =~ m/([0-9]+)/;
    if (defined $1){
        print "$1\n";
    }else{
        print $ret;
    }
}

sub RunHalQsub {
    my $err_log = '';
    $err_log  = '-e '.$error_file   if( defined $error_file);
    $err_log .= ' -o '.$output_file if( defined $output_file);
    my $wor_dir = getcwd;
    my $mem = $memory.'G';
    my $ret = `echo \"$input_cmd\" | qsub -V -q batch -d $wor_dir $err_log -l walltime=$wall_time:00:00,nodes=1:ppn=$thread_n,mem=$mem`;
    $ret =~ m/([0-9]+)/;
    if (defined $1){
        print "$1\n";
    }else{
        print $ret;
    }
}
