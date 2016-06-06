#!/usr/bin/perl -w
# Pipeline for processing RNA-seq from multiple studies (currently GTEx and TCGA studies).

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;
use IO::File;
use Cwd;
use FindBin;
use lib "$FindBin::Bin";


my @usage;
push @usage, "\nUsage:  pipeline-wrapper.pl -t tissue [options]\n\n";
push @usage, "Options:\n";
push @usage, "  -h | --help         Displays this information.\n";
push @usage, "  -c | --config       Configuration file (default: config.txt in the directory of the script\n";
push @usage, "  -t | --tissue       Tissue type, e.g. bladder. User must provide this parameter\n";
push @usage, "  -T | --tissue-conf  Input tissue configuration file (default: tissue-conf.txt)\n";
push @usage, "  -s | --submit       Submit jobs for unprocessed samples or incomplete analysis\n\n";


my ( $config_file, $help, $tissue, $tissue_conf );
my $submit = 0;

GetOptions
(
 'h|help|?'           => \$help,
 'c|config=s'         => \$config_file,
 't|tissue=s'         => \$tissue,
 'T|tissue-config=s'  => \$tissue_conf,
 's|submit'           => \$submit,
);

if ( $help ) {
    print @usage;
    exit(0);
}

if( !defined $tissue ){
    print "ERROR: Please provide tissue type\n";
    print @usage;
    exit(-1);
}

(defined $tissue_conf) or $tissue_conf = "$FindBin::Bin/tissue-conf.txt";
if (!-e $tissue_conf){
    print "ERROR: Cannot find file $tissue_conf\n";
    print @usage;
    die;
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


# Read tissue configuration file
my @lines = `grep $tissue $tissue_conf`;
my $flag = 0;

if (scalar @lines == 0) {
    warn "Warning: did not find $tissue in $tissue_conf\n";
}else{
    my $cluster = -1;
    foreach my $line (@lines){
        my @data = split(/\t/, $line);
        next if ($data[1] ne $tissue and $data[3] ne $tissue);
        
        if ($cluster == -1){
            $cluster = $data[0];
        }elsif ($cluster != $data[0]){
            die "ERROR: $tissue does not have a unique cluster #: $cluster,$data[0]\n";
        }
    }
    
    @lines = ();
    my $prefix;
    foreach my $line (`grep -v ^# $tissue_conf`){
        chomp $line;
        next if(!$line);
        
        my @data = split(/\t/, $line);
        next if ($data[0] != $cluster);
        
        (-e "$gtex_path/sra/$data[1]" or -e "$tcga_path/$data[3]" or -e "$tcga_path/$data[3]-t") or die "ERROR: Nonexistent paths for $data[1]/$data[3]\n";
        
        $flag = 1 if (-e "$gtex_path/sra/$data[1]/SraRunTable.txt" and -e "$tcga_path/$data[3]/summary.tsv");
        push @lines, $line;
    }
    
}



######################### Submit jobs to analyze samples #################################

my $log_file = 'cluster-jobs.log';
my %job_ids;

if (scalar @lines == 0){
    SubmitGTExJobs( $gtex_path . '/sra/' . $tissue );
    SubmitTCGAJobs( $tcga_path . '/' . $tissue );
    SubmitTCGAJobs( $tcga_path . '/' . $tissue . '-t' );
}else{
    my %finished_tissues;
    foreach my $line (@lines){
        my @data = split(/\t/, $line);

        if (-e "$gtex_path/sra/$data[1]" and !defined $finished_tissues{ "$gtex_path/sra/$data[1]" }){
            $finished_tissues{ "$gtex_path/sra/$data[1]" } = 1;
            SubmitGTExJobs( "$gtex_path/sra/$data[1]" );
        }
        if (-e "$tcga_path/$data[3]" and !defined $finished_tissues{ "$tcga_path/$data[3]" }){
            $finished_tissues{ "$tcga_path/$data[3]" } = 1;
            SubmitTCGAJobs( "$tcga_path/$data[3]" );
        }
        if (-e "$tcga_path/$data[3]-t" and !defined $finished_tissues{ "$tcga_path/$data[3]-t" }){
            $finished_tissues{ "$tcga_path/$data[3]-t" } = 1;
            SubmitTCGAJobs( "$tcga_path/$data[3]-t" );
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

######################### Quality Control #################################

my %job_list = GetClusterJobs();

# Process GTEx normals
my $batch_correction = 1;

if (scalar @lines == 0){
    $batch_correction = $batch_correction and GetGTExQC($gtex_path . '/sra/' . $tissue);
    $batch_correction = $batch_correction and GetTCGAQC( $tcga_path . '/' . $tissue );
    $batch_correction = $batch_correction and GetTCGAQC( $tcga_path . '/' . $tissue . '-t' );
}else{
    my %finished_tissues = ();
    foreach my $line (@lines){
        my @data = split(/\t/, $line);
        if (-e "$gtex_path/sra/$data[1]" and !defined $finished_tissues{ "$gtex_path/sra/$data[1]" }){
            $finished_tissues{ "$gtex_path/sra/$data[1]" } = 1;
            $batch_correction = $batch_correction and GetGTExQC("$gtex_path/sra/$data[1]");
        }
        if (-e "$tcga_path/$data[3]" and !defined $finished_tissues{ "$tcga_path/$data[3]" }){
            $finished_tissues{ "$tcga_path/$data[3]" } = 1;
            $batch_correction = $batch_correction and GetTCGAQC("$tcga_path/$data[3]");
        }
        
        next if (scalar @data == 6 and defined $data[5] and $data[5] eq 'control');
        
        if (-e "$tcga_path/$data[3]-t" and !defined $finished_tissues{ "$tcga_path/$data[3]-t" }){
            $finished_tissues{ "$tcga_path/$data[3]-t" } = 1;
            $batch_correction = $batch_correction and GetTCGAQC("$tcga_path/$data[3]-t");
        }
    }
}

######################### Correct batch bias #################################

if($flag == 0){
    warn "Waring: Skip batch bias correction for lack of tissue with both GTEx and TCGA normals\n";
}elsif($batch_correction){
    print "Correcting batch bias...\n\n";
    
    my $cmd_fpkm  = "perl $FindBin::Bin/post-process.pl -t $tissue -c $tissue_conf -u fpkm  -p -r";
    my $cmd_tpm   = "perl $FindBin::Bin/post-process.pl -t $tissue -c $tissue_conf -u tpm   -p -r";
    my $cmd_count = "perl $FindBin::Bin/post-process.pl -t $tissue -c $tissue_conf -u count -p -r";
    
    if ($submit){
        my $ret = `perl $FindBin::Bin/qsub.pl -s \047$cmd_fpkm\047 -p 1 -t 12`;
        print $ret;
        $ret = `perl $FindBin::Bin/qsub.pl -s \047$cmd_tpm\047 -p 1 -t 12`;
        print $ret;
        $ret = `perl $FindBin::Bin/qsub.pl -s \047$cmd_count\047 -p 1 -t 12`;
        print $ret;
    }else{
        system($cmd_fpkm);
        system($cmd_tpm);
        system($cmd_count);
    }
}else{
    print "Skip batch bias correction...\n";
}

print "Done\n";


sub GetGTExQC{
    my $sample_path = shift;
    # Process GTEx normals
    return 0 if (-s "$sample_path/$log_file" and HasJobsRunning("$sample_path/$log_file") );
    
    if(!-s "$sample_path/SraRunTable.txt"){
        print "Warning: File $sample_path/SraRunTable.txt does not exist\n";
        return 0;
    }else{
        my @files = glob("$sample_path/SRR*");
        return 0 if (!@files);
    }
    print "Collecting QC of GTEx normals...\n\n";
    (-d "$sample_path/QC") or `mkdir -p $sample_path/QC`;
    `perl $FindBin::Bin/collect-qc.pl -c $config_file -i $sample_path > $sample_path/QC/qc_sum.txt`;
    print "\nPlease check $sample_path/QC/qc_sum.txt for QC summary\n\n";
    my $m = `grep -v ^Assay_Type_s $sample_path/SraRunTable.txt | wc -l`;
    chomp $m;
    my $n = 0;
    $n = `wc -l < $sample_path/QC/filtered_samples.txt` if(-s "$sample_path/QC/filtered_samples.txt");
    chomp $n;
    print "$n (out of $m) samples were kept in $sample_path\n\n";
    return 1;
}


sub GetTCGAQC{
    my $sample_path = shift;
    # Process TCGA samples
    return 0 if (-s "$sample_path/$log_file" and HasJobsRunning("$sample_path/$log_file") );
    
    if(!-s "$sample_path/summary.tsv"){
        print "Warning: File $sample_path/summary.tsv does not exist\n";
        return 0;
    }else{
        my @files = glob("$sample_path/*-*-*");
        return 0 if (!@files);
    }
    print "Collecting QC for $sample_path...\n\n";
    (-d "$sample_path/QC") or `mkdir -p $sample_path/QC`;
    `perl $FindBin::Bin/collect-qc.pl -c $config_file -i $sample_path > $sample_path/QC/qc_sum.txt`;
    print "\nPlease check $sample_path/QC/qc_sum.txt for QC summary\n\n";
    my $m = `grep -v ^study $sample_path/summary.tsv | wc -l`;
    chomp $m;
    my $n = 0;
    $n = `wc -l < $sample_path/QC/filtered_samples.txt` if(-s "$sample_path/QC/filtered_samples.txt");
    chomp $n;
    print "$n (out of $m) samples were kept in $sample_path\n\n";
    return 1;
}


sub SubmitGTExJobs{
    my $sample_path = shift;
    
    if (!-d $sample_path) {
        warn "Warning: $sample_path do not exist\n\n";
    }else{
        my @files = glob( "$sample_path/*.sra" );
        if (scalar @files == 0){
            warn "Warning: There is no SRA file in $sample_path\n";
        }else{
            warn "Processing GTEx samples in $sample_path\n\n";
            my %sample_job_status = ReadSampleStatus( \@files );
            
            my $log_fh = IO::File->new( "$sample_path/$log_file", ">" ) or die "ERROR: Couldn't write file $sample_path/$log_file\n";
            foreach( @files ){
                chomp;
                my $path = $_;
                $path =~ s/\.sra$//;
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
            if( $submit ){
                map{chomp; my @f=split(/\t/,$_); $job_ids{$f[1]}="$sample_path/$f[0]"}`grep -v \047done\|incomplete\047 $sample_path/$log_file`;
            }
        }
    }
    
}

sub SubmitTCGAJobs{
    my $sample_path = shift;
  
    if (!-d $sample_path) {
        warn "Warning: $sample_path do not exist\n\n";
    }else{
        my @files = glob( "$sample_path/*/*.tar.gz" );
        if (scalar @files == 0){
            my @files1 = glob( "$sample_path/*/UNCID*.bam" );
            if (scalar @files1 == 0){
                warn "Warning: There is no FASTQ/BAM file in $sample_path\n";
            }else{
                @files = @files1;
            }
        }
        
        if (scalar @files > 0){
            warn "Processing TCGA samples in $sample_path\n\n";
            my %sample_job_status = ReadSampleStatus( \@files );
            
            my $log_fh = IO::File->new( "$sample_path/$log_file", ">" ) or die "ERROR: Couldn't write file $sample_path/$log_file\n";
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
            if( $submit ){
                map{chomp; my @f=split(/\t/,$_); $job_ids{$f[1]}="$sample_path/$f[0]"}`grep -v \047done\|incomplete\047 $sample_path/$log_file`;
            }
        }
    }
}

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
