#!/usr/bin/perl -w

use strict;
use warnings FATAL => 'all';
use IO::File;
use FindBin;
use lib "$FindBin::Bin";


####################### Read configuration file ############################

my $config_file = "$FindBin::Bin/config.txt";
( -e $config_file ) or die "ERROR: The configuration file $config_file does not exist\n";
my %config;
map{ chomp; /^\s*([^=\s]+)\s*=\s*(.*)$/; $config{$1} = $2 if (defined $1 && defined $2) } `egrep -v \"^#\" $config_file`;
my $gtex_path  =  $config{ gtex_path };
my $tcga_path  =  $config{ tcga_path };

############################ Count samples ################################

my ($total_sum, $filtered_sum, $unused_sum, $degraded_sum, $poorlymapped_sum) = (0,0,0,0,0);

# Process GTEx sampls
print "GTEx\n";
print "\ttotal\tfiltered\tunused\tdegraded\tpoorlymapped\n";

my @dirs = glob("$gtex_path/sra/*");
foreach (@dirs){
    if(!-s "$_/QC/filtering_sum.txt"){
        s/$gtex_path\/sra\///;
        #print "Warning: File $_/QC/filtering_sum.txt does not exist\n";
        next;
    }
    my ($total, $filtered, $unused, $degraded, $poorlymapped) = CountSamples ($_);
    $total_sum += $total;
    $filtered_sum += $filtered;
    $unused_sum += $unused;
    $degraded_sum += $degraded;
    $poorlymapped_sum += $poorlymapped;
    
    my $filtered_n = `wc -l < $_/QC/filtered_samples.txt`;
    chomp $filtered_n;
    if($filtered != $filtered_n){
        print "ERROR\n";
    }
    
    s/$gtex_path\/sra\///;
    print "$_\t$total\t$filtered\t$unused\t$degraded\t$poorlymapped\n";
}

# Process TCGA sampls
print "TCGA\n";

@dirs = glob("$tcga_path/*");
foreach (@dirs){
    if(!-s "$_/QC/filtering_sum.txt") {
        s/$tcga_path\///;
        #print "Warning: File $_/QC/filtering_sum.txt does not exist\n";
        next;
    }
    my ($total, $filtered, $unused, $degraded, $poorlymapped) = CountSamples ($_);
    $total_sum += $total;
    $filtered_sum += $filtered;
    $unused_sum += $unused;
    $degraded_sum += $degraded;
    $poorlymapped_sum += $poorlymapped;
   
    s/$tcga_path\///;
    print "$_\t$total\t$filtered\t$unused\t$degraded\t$poorlymapped\n";
}


print "Summary\n\ttotal\tfiltered\tunused\tdegraded\tpoorlymapped\n";
print "\t$total_sum\t$filtered_sum\t$unused_sum\t$degraded_sum\t$poorlymapped_sum\n";


sub CountSamples {
    my $path = shift;
    my $qc_file = "$path/QC/filtering_sum.txt";
    
    my ($total, $filtered, $unused, $degraded, $poorlymapped) = (0,0,0,0,0);
    foreach my $line (`grep -v ^# $qc_file`){
        chomp $line;
        next if( !defined $line );
        
        my @data = split (/\t/, $line);
        next if(!-e "$path/$data[0]");
        $filtered++     if (defined $data[2] and $data[2] =~ /Yes/i);
        $unused++       if (defined $data[3] and $data[3] =~ /Unused/);
        if ((defined $data[4] and $data[4] =~ /[0-9]/) or (defined $data[5] and $data[5] =~ /[0-9]/)){
            $poorlymapped++;
        }else{
            $degraded++ if (defined $data[6] and $data[6]=~/[0-9]/);
        }
        $total++;
    }
    return ($total, $filtered, $unused, $degraded, $poorlymapped);
}
