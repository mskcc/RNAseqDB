#!/usr/bin/perl -w
# Pipeline for processing RNA-seq from a study.

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use IO::File;
use Cwd;
use FindBin;
use lib "$FindBin::Bin";


my @usage;
push @usage, "\nUsage:  pipeline-single-study.pl -d directory [options]\n\n";
push @usage, "Options:\n";
push @usage, "  -h | --help         Displays this information.\n";
push @usage, "  -c | --config       Configuration file (default: config.txt in the directory of the script\n";
push @usage, "  -i | --input-dir    Directory of input samples. User must provide this parameter\n";
push @usage, "  -s | --submit       Submit jobs for unprocessed samples or incomplete analysis\n\n";


my ( $config_file, $help, $sample_dir );
my $submit = 0;

GetOptions
(
 'h|help|?'       => \$help,
 'c|config=s'     => \$config_file,
 'i|input-dir=s'  => \$sample_dir,
 's|submit'       => \$submit,
);

if ( $help ) {
    print @usage;
    exit(0);
}

if( !defined $sample_dir ){
    print "ERROR: Please provide directory of your samples\n";
    print @usage;
    exit(-1);
}elsif(!-d $sample_dir){
    print "ERROR: The directory is not correct\n";
    print @usage;
    exit(-1);
}

######################### Read configuration file #################################

(defined $config_file) or $config_file = "$FindBin::Bin/config.txt";
( -e $config_file ) or die "ERROR: The configuration file $config_file does not exist\n";

# Read configuration file
my ( %config, $thread_n);
map{ chomp; /^\s*([^=\s]+)\s*=\s*(.*)$/; $config{$1} = $2 if (defined $1 and defined $2) } `egrep -v \"^#\" $config_file`;
$thread_n = 1;
$thread_n = $config{ thread_no }  if ( exists $config{ defined } );

######################### Submit jobs for each sample #################################

my $log_file = 'cluster-jobs.log';

my @files;
@files = glob( "$sample_dir/*/*.tar.gz" );
@files = glob( "$sample_dir/*/UNCID*.bam" ) if (scalar @files == 0);

if (scalar @files > 0){
    warn "Processing samples in $sample_dir\n\n";
    my %sample_job_status = ReadSampleStatus( \@files );
    
    my $log_fh = IO::File->new( "$sample_dir/$log_file", ">" ) or die "ERROR: Couldn't write file $sample_dir/$log_file\n";
    foreach( @files ){
        chomp;
        my ($name, $path) = fileparse($_);
        chop $path;
        my $sample_id = fileparse($path);
        if ( $sample_job_status{$sample_id} eq 'incomplete' and $submit ){
            my $cmd = "$FindBin::Bin/calc-expression.pl -c $config_file -i $_";
            my $ret = `perl $FindBin::Bin/qsub.pl -s \047$cmd\047`;
            print $ret;
            $log_fh->print("$sample_id\t$ret");
        }else{
            $log_fh->print("$sample_id\t$sample_job_status{$sample_id}");
            $log_fh->print("\tPossible RSEM failure") if(-e "$path/Quant.temp" and $sample_job_status{$sample_id} !~ /^[0-9]+$/);
            $log_fh->print("\n");
        }
    }
    $log_fh->close();
}else{
    foreach (glob( "$sample_dir/*" )){
        next if(!-d $_);
        my ($lane, $id);
        ($lane, $id) = ParseSampleSheet("$_/SampleSheet.csv") if (-e "$_/SampleSheet.csv");
        my @fastq_files;
        if (defined $id){
            @fastq_files = glob( "$_/*$id*fastq*" );
            @fastq_files = glob( "$_/*$id*fq*" )if (scalar @fastq_files == 0);
        }else{
            @fastq_files = glob( "$_/*.fastq.gz" );
            @fastq_files = glob( "$_/*.fastq" )  if (scalar @fastq_files == 0);
            @fastq_files = glob( "$_/*.fq.gz" )  if (scalar @fastq_files == 0);
            @fastq_files = glob( "$_/*.fq" )     if (scalar @fastq_files == 0);
        }
        push (@files, $_) if(scalar @fastq_files > 0);
    }
    
    my %sample_job_status;
    if(scalar @files > 0){
        warn "Processing samples in $sample_dir\n\n";
        %sample_job_status = ReadSampleStatus( \@files );
    }
    
    my $log_fh = IO::File->new( "$sample_dir/$log_file", ">" ) or die "ERROR: Couldn't write file $sample_dir/$log_file\n";
    foreach( @files ){
        chomp;
        my ($sample_id, $path) = fileparse($_);
        chop $path;
        if ( $sample_job_status{$sample_id} eq 'incomplete' and $submit ){
            my $cmd = "$FindBin::Bin/calc-expression.pl -c $config_file -i $_";
            my $ret = `perl $FindBin::Bin/qsub.pl -s \047$cmd\047`;
            print $ret;
            $log_fh->print("$sample_id\t$ret");
        }else{
            $log_fh->print("$sample_id\t$sample_job_status{$sample_id}");
            $log_fh->print("\tPossible RSEM failure") if(-e "$path/Quant.temp" and $sample_job_status{$sample_id} !~ /^[0-9]+$/);
            $log_fh->print("\n");
        }
    }
    $log_fh->close();
}


sub ParseSampleSheet {
    my $sample_sheet = shift;
    
    my ($idx, $lane_idx, $sample_idx) = (0, -1, -1);
    map{$lane_idx=$idx if($_ eq "Lane"); $sample_idx=$idx if($_ eq "SampleID"); $idx++}split(/\,/, `head -1 $sample_sheet`);
    if ($lane_idx != -1 and $sample_idx != -1) {
        my $lane = '';
        my $sample_id = '';
        my @data = split(/\,/, `tail -n +2 $sample_sheet`);
        $lane = $data[$lane_idx] if (defined $data[$lane_idx]);
        $sample_id=$data[$sample_idx]  if (defined $data[$sample_idx]);
        return ($lane, $sample_id);
    }else{
        warn "ERROR Unknown file format: $sample_sheet\n";
        return ('','');
    }
}

sub GetClusterJobs {
    my @jobs = `perl $FindBin::Bin/qstat.pl`;

    my %job_list;
    %job_list = map{chomp; ($_, 1)} @jobs;
    return %job_list;
}

sub ReadSampleStatus{
    my $samples = shift;
    my $fist_sample = $samples->[0];
    my $study_path;
    if ($fist_sample =~ /\.sra$/){
        my $name;
        ($name, $study_path) = fileparse($fist_sample);
    }elsif(!-d $fist_sample){
        my $name;
        ($name, $study_path) = fileparse($fist_sample);
        chop $study_path;
        ($name, $study_path) = fileparse($study_path);
    }else{
        my $name;
        ($name, $study_path) = fileparse($fist_sample);
    }
    
    my %job_list = GetClusterJobs();
    my %sample_status = ();
    if (-s "$study_path/$log_file"){
        foreach(`cat $study_path/$log_file`){
            chomp;
            next if(!$_);
            my @data = split(/[\t ]+/, $_);
            if ($data[1] eq 'done' or exists $job_list{$data[1]}){
                $sample_status{ $data[0] } = $data[1];
            }else{
                my $path = "$study_path/$data[0]";
                if (-e $path and -s "$path/Log.final.out" and -s "$path/Quant.genes.results" and -s "$path/ubu-quan/fcounts.fpkm.normalized_results" and -s "$path/ks/sample.ks.txt" and -e "$path/Aligned.sortedByCoord.out_fastqc"){
                    $sample_status{ $data[0] } = 'done';
                }else{
                    $sample_status{ $data[0] } = 'incomplete';
                }
            }
        }
    }
    
    foreach(@$samples){
        chomp;
        my $sample_id;
        if (/\.sra$/){
            $sample_id = fileparse($_);
            $sample_id =~s/\.sra$//;
        }elsif(!-d $_){
            my $temp;
            ($sample_id, $temp) = fileparse($_);
            chop $temp;
            $sample_id = fileparse($temp);
        }else{
            my $temp;
            ($sample_id, $temp) = fileparse($_);
        }
        next if (exists $sample_status{ $sample_id });
        
        my $path = "$study_path/$sample_id";
        if (-e $path and -s "$path/Log.final.out" and -s "$path/Quant.genes.results" and -s "$path/ubu-quan/fcounts.fpkm.normalized_results" and -s "$path/ks/sample.ks.txt" and -e "$path/Aligned.sortedByCoord.out_fastqc"){
            $sample_status{ $sample_id } = 'done';
        }else{
            $sample_status{ $sample_id } = 'incomplete';
        }
    }
    return %sample_status;
}
