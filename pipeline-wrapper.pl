#!/usr/bin/perl -w
# Pipeline wrapper for analyzing RNA-seq of GTEx and TCGA

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
push @usage, "\nUsage:  pipeline-wrapper.pl -t tissue [options]\n\n";
push @usage, "Options:\n";
push @usage, "  -h | --help      Displays this information.\n";
push @usage, "  -t | --tissue    Tissue type, e.g. bladder. User must provide this parameter\n";
push @usage, "  -c | --config    A configuration file, default is config.txt in same dir of this script\n";
push @usage, "  -s | --submit    Submit jobs for incomplete analysis if specified\n\n";


my ( $help, $tissue_type, $config_file );
my $submit = 0;

GetOptions
(
 'h|help|?'    => \$help,
 'i|tissue=s'  => \$tissue_type,
 'c|config=s'  => \$config_file,
 's|submit'    => \$submit,
);

if ( $help ) {
    print @usage;
    exit(0);
}

if( !defined $tissue_type ){
    print "ERROR: Please provide tissue type\n";
    print @usage;
    exit(-1);
}

######################### Read configuration file #################################

(defined $config_file) or $config_file = "$FindBin::Bin/config.txt";
( -e $config_file ) or die "ERROR: The configuration file $config_file does not exist\n";

# Read configuration file
my ( %config, $thread_n, $gtex_path, $tcga_path );
map{ chomp; /^\s*([^=\s]+)\s*=\s*(.*)$/; $config{$1} = $2 if (defined $1 and defined $2) } `egrep -v \"^#\" $config_file`;
$thread_n   =  1;
$thread_n   =  $config{ thread_no }  if ( exists $config{ thread_no } );
$gtex_path  =  $config{ gtex_path };
$tcga_path  =  $config{ tcga_path };


######################### Submit jobs to analyze samples #################################

my $gtex_normal_path = $gtex_path . '/sra/' . $tissue_type;
my $tcga_normal_path = $tcga_path . '/' . $tissue_type;
my $tcga_tumor_path  = $tcga_path . '/' . $tissue_type . '-t';

my $log_file = 'cluster-jobs.log';
my %job_ids;

if (!-d $gtex_normal_path) {
    warn "Warning: There is no directory for GTEx samples\n\n";
}else{
    warn "Processing GTEx samples\n\n";
    my @files = glob( "$gtex_normal_path/*.sra" );
    if (scalar @files == 0){
        warn "Warning: There is no SRA file for GTEx\n";
    }else{
        my %sample_job_status = ReadSampleStatus( \@files );

        my $log_fh = IO::File->new( "$gtex_normal_path/$log_file", ">" ) or die "ERROR: Couldn't write file $gtex_normal_path/$log_file\n";
        foreach( @files ){
            chomp;
            my $path = $_;
            $path =~ s/\.sra$//;
            my $sample_id = fileparse($path);
            if ( $sample_job_status{$sample_id} ne 'incomplete' or !$submit ){
                $log_fh->print("$sample_id\t$sample_job_status{$sample_id}");
                $log_fh->print("\tPossible RSEM failure") if(-e "$path/Quant.temp" and $sample_job_status{$sample_id} !~ /^[0-9]+$/);
                $log_fh->print("\n");
            }else{
                my $cmd = "$FindBin::Bin/calc-expression.pl -i $_";
                my $ret = `perl $FindBin::Bin/qsub.pl -s \047$cmd\047`;
                print $ret;
                $log_fh->print("$sample_id\t$ret");
            }
        }
        $log_fh->close();
        if( $submit ){
            map{chomp; my @f=split(/\t/,$_); $job_ids{$f[1]}="$gtex_normal_path/$f[0]"}`grep -v \047done\|incomplete\047 $gtex_normal_path/$log_file`;
        }
    }
}


if (!-d $tcga_normal_path) {
    warn "Warning: There is no directory for TCGA normal samples\n\n";
}else{
    warn "Processing TCGA normal samples\n\n";
    my @files = glob( "$tcga_normal_path/*/*.tar.gz" );
    if (scalar @files == 0){
        my @files1 = glob( "$tcga_normal_path/*/UNCID*.bam" );
        if (scalar @files1 == 0){
            warn "Warning: There is no FASTQ/BAM file for TCGA normal\n";
        }else{
            @files = @files1;
        }
    }
    
     if (scalar @files > 0){
        my %sample_job_status = ReadSampleStatus( \@files );
        
        my $log_fh = IO::File->new( "$tcga_normal_path/$log_file", ">" ) or die "ERROR: Couldn't write file $tcga_normal_path/$log_file\n";
        foreach( @files ){
            chomp;
            my ($name, $path) = fileparse($_);
            chop $path;
            my $sample_id = fileparse($path);
            if ( $sample_job_status{$sample_id} ne 'incomplete' or !$submit ){
                $log_fh->print("$sample_id\t$sample_job_status{$sample_id}");
                $log_fh->print("\tPossible RSEM failure") if(-e "$path/Quant.temp" and $sample_job_status{$sample_id} !~ /^[0-9]+$/);
                $log_fh->print("\n");
            }else{
                my $cmd = "$FindBin::Bin/calc-expression.pl -i $_";
                my $ret = `perl $FindBin::Bin/qsub.pl -s \047$cmd\047`;
                print $ret;
                $log_fh->print("$sample_id\t$ret");
            }
        }
        $log_fh->close();
        if( $submit ){
            map{chomp; my @f=split(/\t/,$_); $job_ids{$f[1]}="$tcga_normal_path/$f[0]"}`grep -v \047done\|incomplete\047 $tcga_normal_path/$log_file`;
        }
    }
}


if (!-d $tcga_tumor_path) {
    warn "Warning: There is no directory for TCGA tumor samples\n\n";
}else{
    warn "Processing TCGA tumor samples\n\n";
    my @files = glob( "$tcga_tumor_path/*/*.tar.gz" );
    if (scalar @files == 0){
        my @files1 = glob( "$tcga_tumor_path/*/UNCID*.bam" );
        if (scalar @files1 == 0){
            warn "Warning: There is no FASTQ/BAM file for TCGA tumor\n";
        }else{
            @files = @files1;
        }
    }
    
    if (scalar @files > 0){
        my %sample_job_status = ReadSampleStatus( \@files );
        
        my $log_fh = IO::File->new( "$tcga_tumor_path/$log_file", ">" ) or die "ERROR: Couldn't write file $tcga_tumor_path/$log_file\n";
        foreach( @files ){
            chomp;
            my ($name, $path) = fileparse($_);
            chop $path;
            my $sample_id = fileparse($path);
            if ( $sample_job_status{$sample_id} ne 'incomplete' or !$submit ){
                $log_fh->print("$sample_id\t$sample_job_status{$sample_id}");
                $log_fh->print("\tPossible RSEM failure") if(-e "$path/Quant.temp" and $sample_job_status{$sample_id} !~ /^[0-9]+$/);
                $log_fh->print("\n");
            }else{
                my $cmd = "$FindBin::Bin/calc-expression.pl -i $_";
                my $ret = `perl $FindBin::Bin/qsub.pl -s \047$cmd\047`;
                print $ret;
                $log_fh->print("$sample_id\t$ret");
            }
        }
        $log_fh->close();
        if( $submit ){
            map{chomp; my @f=split(/\t/,$_); $job_ids{$f[1]}="$tcga_tumor_path/$f[0]"}`grep -v \047done\|incomplete\047 $tcga_tumor_path/$log_file`;
        }
    }
}


######################### Wait for qsub jobs to finish #################################


if ($submit and  scalar(keys %job_ids) > 0){
    print "Waiting for jobs to terminate...\n";
    foreach my $job_id ( keys %job_ids ){
        my %job_list;
        do{
            sleep 30;
            %job_list = GetClusterJobs();
        }
        while( exists $job_list{$job_id} );
        
        my $path = $job_ids{$job_id};
        if (-e $path and -s "$path/Log.final.out" and -s "$path/Quant.genes.results" and -s "$path/ubu-quan/fcounts.fpkm.normalized_results" and -s "$path/ks/sample.ks.txt" and -e "$path/Aligned.sortedByCoord.out_fastqc"){
            print "$path\tDone\n";
        }else{
            print "$path\tIncomplete\n";
        }
    }
}

######################### Filter samples #################################

my %job_list = GetClusterJobs();

# Process GTEx normals
my ($gtex_qc, $batch_correction) = (1,1);
$gtex_qc = 0 if (-s "$gtex_normal_path/$log_file" and HasJobsRunning("$gtex_normal_path/$log_file") );

if(!-s "$gtex_normal_path/SraRunTable.txt"){
    print "Warning: File $gtex_normal_path/SraRunTable.txt does not exist\n";
    $gtex_qc=0;
}else{
    my @files = glob("$gtex_normal_path/SRR*");
    $gtex_qc=0 if (!@files);
}
if($gtex_qc){
    print "Collecting QC of GTEx normals...\n\n";
    (-d "$gtex_normal_path/QC") or `mkdir -p $gtex_normal_path/QC`;
    `perl $FindBin::Bin/collect-qc.pl -i $gtex_normal_path > $gtex_normal_path/QC/qc_sum.txt`;
    print "\nPlease check $gtex_normal_path/QC/qc_sum.txt for QC summary\n\n";
    my $m = `grep -v ^Assay_Type_s $gtex_normal_path/SraRunTable.txt | wc -l`;
    chomp $m;
    my $n = 0;
    $n = `wc -l < $gtex_normal_path/QC/filtered_samples.txt` if(-s "$gtex_normal_path/QC/filtered_samples.txt");
    chomp $n;
    print "$n (out of $m) samples were kept after filtering\n\n";
}else{
    print "Warning: Skip doing QC for incomplete analysis of GTEx normals\n";
    $batch_correction=0;
}

# Process TCGA normals
my $tcga_qc = 1;
$tcga_qc = 0 if (-s "$tcga_normal_path/$log_file" and HasJobsRunning("$tcga_normal_path/$log_file") );

if(!-s "$tcga_normal_path/summary.tsv"){
    print "Warning: File $tcga_normal_path/summary.tsv does not exist\n";
    $tcga_qc=0;
}else{
    my @files = glob("$tcga_normal_path/*-*-*");
    $tcga_qc=0 if (!@files);
}
if($tcga_qc){
    print "Collecting QC of TCGA normals...\n\n";
    (-d "$tcga_normal_path/QC") or `mkdir -p $tcga_normal_path/QC`;
    `perl $FindBin::Bin/collect-qc.pl -i $tcga_normal_path > $tcga_normal_path/QC/qc_sum.txt`;
    print "\nPlease check $tcga_normal_path/QC/qc_sum.txt for QC summary\n\n";
    my $m = `grep -v ^study $tcga_normal_path/summary.tsv | wc -l`;
    chomp $m;
    my $n = 0;
    $n = `wc -l < $tcga_normal_path/QC/filtered_samples.txt` if(-s "$tcga_normal_path/QC/filtered_samples.txt");
    chomp $n;
    print "$n (out of $m) samples were kept after filtering\n\n";
}else{
    print "Warning: Skip doing QC for incomplete analysis of TCGA normals\n\n";
    $batch_correction=0;
}

# Process TCGA tumors
my $tcga_t_qc = 1;
$tcga_t_qc = 0 if (-s "$tcga_tumor_path/$log_file" and HasJobsRunning("$tcga_tumor_path/$log_file")  );

if(!-s "$tcga_tumor_path/summary.tsv"){
    print "Warning: File $tcga_tumor_path/summary.tsv does not exist\n";
    $tcga_t_qc=0;
}else{
    my @files = glob("$tcga_tumor_path/*-*-*");
    $tcga_t_qc=0 if (!@files);
}
if($tcga_t_qc){
    print "Collecting QC of TCGA tumors...\n\n";
    (-d "$tcga_tumor_path/QC") or `mkdir -p $tcga_tumor_path/QC`;
    `perl $FindBin::Bin/collect-qc.pl -i $tcga_tumor_path > $tcga_tumor_path/QC/qc_sum.txt`;
    print "\nPlease check $tcga_tumor_path/QC/qc_sum.txt for QC summary\n\n";
    my $m = `grep -v ^study $tcga_tumor_path/summary.tsv | wc -l`;
    chomp $m;
    my $n = 0;
    $n = `wc -l < $tcga_tumor_path/QC/filtered_samples.txt` if(-s "$tcga_tumor_path/QC/filtered_samples.txt");
    chomp $n;
    print "$n (out of $m) samples were kept after filtering\n\n";
}else{
    print "Warning: Skip doing QC for incomplete analysis of TCGA tumors\n\n";
    $batch_correction=0;
}


######################### Correct batch bias #################################

if($batch_correction){
    print "Correcting batch bias...\n\n";
    `perl $FindBin::Bin/run-combat.pl -t $tissue_type -u fpkm -r`;
    `perl $FindBin::Bin/run-combat.pl -t $tissue_type -u tpm  -r`;
    `perl $FindBin::Bin/run-combat.pl -t $tissue_type -u count`;
}else{
    print "Skip batch bias correction...\n";
}

print "Done\n";


sub GetClusterJobs {
    my @jobs = `perl $FindBin::Bin/qstat.pl`;

    my %job_list;
    %job_list = map{chomp; ($_, 1)} @jobs;
    return %job_list;
}

sub HasJobsRunning {
    return 0 if (scalar(keys %job_list) == 0);
    
    my $log_file = shift;
    map{ chomp; return 1 if (defined $job_list{$_}) }`cut -f 2 $log_file`;
    return 0;
}

sub ReadSampleStatus{
    my $samples = shift;
    my $fist_sample = $samples->[0];
    my $study_path;
    if ($fist_sample =~ /\.sra$/){
        my $name;
        ($name, $study_path) = fileparse($fist_sample);
    }else{
        my $name;
        ($name, $study_path) = fileparse($fist_sample);
        chop $study_path;
        ($name, $study_path) = fileparse($study_path);
    }
    
    my %job_list = GetClusterJobs();
    my %sample_status = ();
    if (-s "$study_path/$log_file"){
        foreach(`cat $study_path/$log_file`){
            chomp;
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
        }else{
            my $temp;
            ($sample_id, $temp) = fileparse($_);
            chop $temp;
            $sample_id = fileparse($temp);
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
