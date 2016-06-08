#!/usr/bin/perl -w

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin";

my @usage;
push @usage, "\nUsage:  qstat.pl [options]\n\n";
push @usage, "Options:\n";
push @usage, "  -h | --help    Displays this information.\n";
push @usage, "  -c | --config  A configuration file, default is config.txt in script dir\n";
push @usage, "  -f | --full    Print a full status display to standard out\n";


my ($help, $full_output) = (0, 0);
my ( $config_file );

GetOptions
(
 'h|help|?'   => \$help,
 'c|config=s' => \$config_file,
 'f|full'     => \$full_output,
);

if ( $help ) {
    print @usage;
    exit(0);
}

######################### Processes configuration file #################################

(defined $config_file) or $config_file = "$FindBin::Bin/config.txt";
( -e $config_file ) or die "ERROR: The configuration file $config_file does not exist\n";

# Read configuration file
my ( %config, $cluster );
map{ chomp; /^\s*([^=\s]+)\s*=\s*(.*)$/; $config{$1} = $2 if (defined $1 and defined $2) } `egrep -v \"^#\" $config_file`;
$cluster =  $config{ cluster } if ( exists $config{ cluster } );

######################### Print jobs in the cluster #################################

if(!defined $cluster){
    die "ERROR: Variable 'cluster' is not defined in configuration file\n";
}elsif(lc($cluster) eq "luna"){
    RunLunaQStat();
}elsif(lc($cluster) eq "hal"){
    RunHalQStat();
}else{
    die "ERROR: Unknown cluster\n";
}


sub RunLunaQStat {
    my @jobs;
    if ($full_output){
        @jobs = `bjobs 2>/dev/null`;
    }else{
        @jobs = `bjobs 2>/dev/null | grep -v JOBID | awk \047\{print \$1\}\047`;
    }
    map{ print "$_" }@jobs;
}

sub RunHalQStat {
    my $logname = $ENV{ LOGNAME };
    my @jobs;
    if ($full_output){
        @jobs = `qstat -u $logname`;
        map{ print "$_" }@jobs;
    }else{
        @jobs = `qstat -u $logname | grep $logname`;
        foreach(@jobs){
            my @data = split(/[\t ]/, $_);
            my $jobid = $data[0];
            $jobid =~ m/^([0-9]+)/;
            print "$1\n";
        }
    }
}

