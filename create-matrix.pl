#!/usr/bin/perl -w
# Create a sample-gene matrix for all samples in a study

use strict;
use warnings FATAL => 'all';
use File::Basename;
use IO::File;
use Getopt::Std;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin";
use Statistics::Descriptive;
use Cwd;
use File::Temp qw( tempdir );


my @usage;
push @usage, "\n-----------------------------------------------------------\n\n";
push @usage, "Usage:  create-matrix.pl --input directory\n\n";
push @usage, "Options:\n";
push @usage, "  -h | --help             Display this information\n";
push @usage, "  -c | --config           Configuration file\n";
push @usage, "  -i | --input            Input directory\n";
push @usage, "  -o | --output           Output file\n";
push @usage, "  -t | --quan-tool        Expression quantification tool: RSEM (default) | FeatureCounts\n";
push @usage, "  -u | --quan-unit        Unit to measure gene expression: TPM | Count | FPKM (default)\n";
push @usage, "  -g | --transcript-type  Output data types: gene (default) | transcript\n";
push @usage, "  -p | --protein-coding   Output protein-coding genes/transcripts only if specified\n\n";


############################## Input parameters ##############################
# Read input parameters

my ($help, $input, $out_file, $config_file, $quan_tool, $quan_unit, $gene_type);
my $protein_coding = 0;

GetOptions
(
    'h|help|?'       => \$help,
    'i|input=s'      => \$input,
    'o|output=s'     => \$out_file,
    'c|config=s'     => \$config_file,
    't|quan-tool=s'  => \$quan_tool,
    'u|quan-unit=s'  => \$quan_unit,
    'g|gene-type=s'  => \$gene_type,
    'p|protein-coding'  => \$protein_coding,
);

if ($help) {
    print @usage;
    exit(0);
}

if (defined $input) {
    if (!-e $input){
        print "\nError: $input directory does not exist!\n";
        print @usage;
        exit;
    }
}else{
    print "\nError: input directory was not provided!\n";
    print @usage;
    exit(0);
}

if (!defined $quan_tool) {
    $quan_tool = 'rsem';
}else{
    $quan_tool = lc($quan_tool);
    if ($quan_tool ne 'rsem' and $quan_tool ne 'featurecounts'){
        print "ERROR: Unknown tool. Please specify either RSEM or FeatureCounts using -q\n";
        print @usage;
        die;
    }
}


if (!defined $quan_unit) {
    $quan_unit = 'fpkm';
}else{
    $quan_unit = lc($quan_unit);
    if ($quan_unit ne 'tpm' and $quan_unit ne 'count' and $quan_unit ne 'fpkm'){
        print "ERROR: Please specify one of TPM | Count | FPKM using -u\n";
        print @usage;
        die;
    }
}

if (defined $gene_type){
    $gene_type = lc($gene_type);
    ($gene_type eq 'transcript' or $gene_type eq 'gene') or die "ERROR: unknown transcript type\n";
    if($gene_type eq 'transcript' and $quan_tool eq 'featurecounts'){
        print "Warning: transcript result is not available for FeatureCounts\n";
        print @usage;
        die;
    }
}else{
    $gene_type = 'gene';
}


############################# Read configuration file ############################


(defined $config_file) or $config_file = "$FindBin::Bin/config.txt";
( -e $config_file ) or die "ERROR: The configuration file $config_file does not exist\n";
my %config;
map{ chomp; /^\s*([^=\s]+)\s*=\s*(.*)$/; $config{$1} = $2 if (defined $1 && defined $2) } `egrep -v \"^#\" $config_file`;

# Use the configuration file to initialize variables
my $gene_name_file = $config{ ENSG_UCSC_common_genes };
my $ubu_dir = $config{ ubu_dir } if ( defined $config{ ubu_dir } );


########################### Read meta data for sample barcodes ################################


my $work_dir = getcwd;
chdir $input;

# Map sample ID to barcode
my %barcode_hash;
my @sample_list;

if (-d "QC" and -e "QC/filtered_samples.txt"){
    map{chomp; my @d=split(/\t/); push(@sample_list,$d[0]); $barcode_hash{$d[0]}=$d[1]}`cat QC/filtered_samples.txt`;
}else{
    # GTEx barcode
    if (-e "SraRunTable.txt") {
        # Determine the column Sample_Name_s in
        my ($idx, $barcode_idx, $run_idx) = (0, -1, -1);
        map{$idx++; $run_idx = $idx if($_ eq "Run_s"); $barcode_idx = $idx if($_ eq "Sample_Name_s")}split(/\t/, `head SraRunTable.txt | grep ^Assay_Type_s`);
        ($barcode_idx != -1 and $run_idx != -1) or die "Unknown file format: SraRunTable.txt\n";
        
        foreach(`grep -v ^Assay_Type_s SraRunTable.txt | cut -f $run_idx,$barcode_idx`){
            chomp;
            my @f = split(/\t/,$_);
            if (-s "$f[0]/Log.final.out" and -s "$f[0]/Quant.genes.results" and -s "$f[0]/ubu-quan/fcounts.fpkm.normalized_results" and -s "$f[0]/ks/sample.ks.txt" and -e "$f[0]/Aligned.sortedByCoord.out_fastqc"){
                $barcode_hash{$f[0]} = $f[1];
                push(@sample_list,$f[0]);
            }
        }
    # TCGA barcode
    }elsif (-e "summary.tsv"){
        my @lines = `head -2 summary.tsv`;
        if ($lines[0] =~ /^study/ and $lines[1] =~ /^TCGA/){
            my ($idx, $barcode_idx, $analysis_idx) = (0, -1, -1);
            map{$idx++; $barcode_idx = $idx if($_ eq "barcode"); $analysis_idx = $idx if($_ eq "analysis_id")}split(/\t/, `head summary.tsv | grep ^study`);
            ($barcode_idx != -1 and $analysis_idx != -1) or die "Unknown file format: summary.tsv\n";
            
            foreach(`grep -v ^study summary.tsv | cut -f $barcode_idx,$analysis_idx`){
                chomp;
                my @f = split(/\t/,$_);
                $barcode_hash{$f[1]} = $f[0];
                if (-s "$f[1]/Log.final.out" and -s "$f[1]/Quant.genes.results" and -s "$f[1]/ubu-quan/fcounts.fpkm.normalized_results" and -s "$f[1]/ks/sample.ks.txt" and -e "$f[1]/Aligned.sortedByCoord.out_fastqc"){
                    $barcode_hash{$f[0]} = $f[1];
                    push(@sample_list,$f[0]);
                }
            }
        }
    }
    
    # If neither GTEx nor TCGA
    if (0 == scalar @sample_list){
        #my @files = `find .. -name Log.final.out`;
        my @files = glob("./*");
        foreach( @files ){
            s/^\.\///;
            if(-e "$_/Log.final.out"){
                if (-s "$_/Quant.genes.results" and -s "$_/ubu-quan/fcounts.fpkm.normalized_results" and -s "$_/ks/sample.ks.txt" and -e "$_/Aligned.sortedByCoord.out_fastqc"){
                    $barcode_hash{$_} = $_;
                    push(@sample_list,$_);
                }
            }else{
                my @sample_sheet = glob( "$_/*SampleSheet.csv" );
                next if(scalar @sample_sheet == 0);
                
                my ($lane, $id) = ParseSampleSheet( $sample_sheet[0] );
                next if (!defined $id or !defined $lane);
                for (my $i=1; $i<=$lane; $i++){
                    next if(!-e "$_/L$i");
                    if (-s "$_/L$i/Quant.genes.results" and -s "$_/L$i/ubu-quan/fcounts.fpkm.normalized_results" and -s "$_/L$i/ks/sample.ks.txt" and -e "$_/L$i/Aligned.sortedByCoord.out_fastqc"){
                        $barcode_hash{"$_/L$i"} = "$_-L$i";
                        push(@sample_list,"$_/L$i");
                    }
                }
            }
        }
    }
}


( @sample_list ) or die "ERROR: Cannot find quantification file\n";
#( scalar @sample_list > 1 ) or die "ERROR: Skip creating data matrix as sample number is 1\n";


########################### Get genes of interest ################################

if ($gene_type eq 'transcript'){
    $gene_name_file =~ s/_ucsc\.known/\.transcript/;
}else{
    $gene_name_file =~ s/_ucsc\.known/\.gene/;
}

my %ens2hugo;
if ( $protein_coding ){
    map{chomp; my @f=split(/\t/); $ens2hugo{$f[0]}=$f[1]}`grep protein_coding $gene_name_file | cut -f 1-2`;
}else{
    map{chomp; my @f=split(/\t/); $ens2hugo{$f[0]}=$f[1]}`cut -f 1-2 $gene_name_file`;
}


########################### Read data files into memory ################################


my $tmp_dir = tempdir( CLEANUP => 1 );

my (@samples, %data_matrix, %gene_list);
foreach my $line (@sample_list){
    chomp $line;
    
    my $quan_file;
    if($quan_tool eq "rsem"){
        $quan_file = "$line/Quant.genes.results";
    }else{
        $quan_file = "$line/fcounts.".$quan_unit;
    }
    if(!-s $quan_file) {
        print "Warning: Skip $line\n";
        next;
    }

    if($quan_unit eq 'fpkm'){   # Do quantile normalization for fpkm
        if ($gene_type eq 'gene'){ # gene expression
            if( !$protein_coding ){
                if($quan_tool eq "rsem"){
                    $quan_file = "$line/ubu-quan/rsem.genes.normalized_results";
                }else{
                    $quan_file = "$line/ubu-quan/fcounts.fpkm.normalized_results";
                }
            }else{
                if($quan_tool eq "rsem"){
                    my $column = -1;
                    my $idx = 0;
                    map{ $idx++; $column = $idx if($_ eq "FPKM") }split(/\t/, `head $quan_file | grep ^gene_id`);
                    ($column != -1) or die "Unknown file format: $quan_file\n";
                    `grep -v gene_id $quan_file | cut -f 1,$column > $tmp_dir/temp0.txt`;
                    $quan_file = "$tmp_dir/temp0.txt";
                }
                my $tmp_header = IO::File->new( "$tmp_dir/temp1.txt", ">" ) or die "ERROR: Couldn't create file $tmp_dir/temp1.txt\n";
                foreach(`cat $quan_file`){
                    my @e=split(/\t/, $_);
                    next if(! exists $ens2hugo{$e[0]});
                    $tmp_header->print($_);
                }
                $tmp_header->close;
                $quan_file = "$tmp_dir/temp1.txt";
                `perl $ubu_dir/perl/quartile_norm.pl -c 2 -q 75 -t 1000 -o $tmp_dir/temp2.txt $quan_file`;
                $quan_file = "$tmp_dir/temp2.txt";
            }
        }else{ # transcript expression
            if($quan_tool eq "rsem"){
                $quan_file = "$line/Quant.isoforms.results";
                my $column = -1;
                my $idx = 0;
                map{ $idx++; $column = $idx if($_ eq "FPKM") }split(/\t/, `head $quan_file | grep ^transcript_id`);
                ($column != -1) or die "Unknown file format: $quan_file\n";
                `grep -v transcript_id $quan_file | cut -f 1,$column > $tmp_dir/temp0.txt`;
                $quan_file = "$tmp_dir/temp0.txt";
            }
            my $tmp_header = IO::File->new( "$tmp_dir/temp1.txt", ">" ) or die "ERROR: Couldn't create file $tmp_dir/temp1.txt\n";
            foreach(`cat $quan_file`){
                my @e=split(/\t/, $_);
                next if(! exists $ens2hugo{$e[0]});
                $tmp_header->print($_);
            }
            $tmp_header->close;
            $quan_file = "$tmp_dir/temp1.txt";
            `perl $ubu_dir/perl/quartile_norm.pl -c 2 -q 75 -t 300 -o $tmp_dir/temp2.txt $quan_file`;
            $quan_file = "$tmp_dir/temp2.txt";
        }
    }
    
    push @samples, $line;
    my $n=0;
    if($quan_tool eq "rsem" and $quan_unit ne "fpkm"){
        $quan_file = "$line/Quant.isoforms.results" if($gene_type eq "transcript");
        my $column = -1;
        my $idx = 0;
        my $col_name = "TPM";
        $col_name = "expected_count" if ($quan_unit eq "count");
        if($gene_type eq "transcript"){
            map{ $idx++; $column = $idx if($_ eq $col_name) }split(/\t/, `head $quan_file | grep ^transcript_id`);
        }else{
            map{ $idx++; $column = $idx if($_ eq $col_name) }split(/\t/, `head $quan_file | grep ^gene_id`);
        }
        ($column != -1) or die "Unknown file format: $quan_file\n";
        map{chomp; my @e=split(/\t/); if(defined $ens2hugo{$e[0]}){$gene_list{$e[0]}=1; $data_matrix{$line}{$e[0]}=$e[1];$n++}}`cut -f 1,$column $quan_file | tail -n +2`;
    }else{
        map{chomp; my @e=split(/\t/); if(defined $ens2hugo{$e[0]}){$gene_list{$e[0]}=1; $data_matrix{$line}{$e[0]}=$e[1];$n++}}`cat $quan_file`;
    }
    print STDERR "$line\t$n\n";
}


########################### Print data matrix ################################


chdir $work_dir;

my @header_ids;
map{$_ = $barcode_hash{$_} if(defined $barcode_hash{$_}); push (@header_ids, $_)}@samples;

# Store the matrix
if (defined $out_file){
    my $r_fh = IO::File->new( $out_file, ">" ) or die "ERROR: Failed to create file $out_file\n";
    $r_fh->print( "Gene\tDescription\t". join("\t", @header_ids) . "\n" );
    foreach my $gene (sort {$gene_list{$a} <=> $gene_list{$b}} (keys %gene_list)){
        $r_fh->print("$gene\t$ens2hugo{$gene}");
        foreach my $sample ( @samples ){
            my $expr = defined( $data_matrix{$sample}{$gene} ) ? $data_matrix{$sample}{$gene} : 0;
            $r_fh->print( "\t$expr" );
        }
        $r_fh->print( "\n" );
    }
    $r_fh->close;
}else{
    print( "Gene\tDescription\t". join("\t", @header_ids) . "\n" );
    foreach my $gene (sort {$gene_list{$a} <=> $gene_list{$b}} (keys %gene_list)){
        print("$gene\t$ens2hugo{$gene}");
        foreach my $sample ( @samples ){
            my $expr = defined( $data_matrix{$sample}{$gene} ) ? $data_matrix{$sample}{$gene} : 0;
            print( "\t$expr" );
        }
        print( "\n" );
    }
}

sub ParseSampleSheet {
    my $sample_sheet = shift;
    
    my ($idx, $lane_idx, $sample_idx) = (0, -1, -1);
    map{chomp;$lane_idx=$idx if($_ eq "Lane"); $sample_idx=$idx if($_ eq "SampleID"); $idx++}split(/\,/, `head -1 $sample_sheet`);
    if ($lane_idx != -1 and $sample_idx != -1) {
        my ($lane, $sample_id);
        foreach(`tail -n +2 $sample_sheet`){
            chomp;
            my @data = split(/\,/, $_);
            if (defined $data[$lane_idx]){
                $lane = $data[$lane_idx] if (!defined $lane or $lane < $data[$lane_idx]);
            }
            $sample_id=$data[$sample_idx]  if (defined $data[$sample_idx]);
        }
        return ($lane, $sample_id);
    }else{
        warn "ERROR Unknown file format: $sample_sheet\n";
    }
}

