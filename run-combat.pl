#!/usr/bin/perl -w

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use Cwd;
use FindBin;
use lib "$FindBin::Bin";
use IO::File;
use List::Util qw(min max);


### process input parameters
my @usage;
push @usage, "\nUsage:  rm-batch.pl -t tissueType\n\n";
push @usage, "Options:\n";
push @usage, "  -h, --help          Displays this information.\n";
push @usage, "  -t, --tissue        Tissue type\n";
push @usage, "  -a, --tissue2       Another TCGA tissue type to run ComBat together with the 1st tissue\n";
push @usage, "  -b, --tissue3       Another TCGA tissue type to run ComBat together with other tissues\n";
push @usage, "  -n, --norm-method   Method to normalize gene expression: TPM | Count | RSEM | FeatureCounts\n\n";

my ($help, $tissue, $tissue2, $tissue3, $normalization);

GetOptions
(
 'h|help|?'         => \$help,
 't|tissue=s'       => \$tissue,
 'a|tissue2=s'      => \$tissue2,
 'b|tissue3=s'      => \$tissue3,
 'n|norm-method=s'  => \$normalization,
);

if (!defined $normalization) {
    print "ERROR: you did not specify arugment -n\n";
    print @usage;
    die;
}else{
    $normalization = lc($normalization);
}

if ($help) {
   print @usage;
   exit(0);
}

######################### Read configuration file #################################

my $config_file = "$FindBin::Bin/config.txt";
( -e $config_file ) or die "ERROR: The configuration file $config_file does not exist\n";
my %config;
map{ chomp; /^\s*([^=\s]+)\s*=\s*(.*)$/; $config{$1} = $2 if (defined $1 && defined $2) } `egrep -v \"^#\" $config_file`;
# Use configuration file to initialize variables
my ($gtex_path, $tcga_path, $ubu_dir, $gtex_sample_attr, $ensg_2_hugo_file);
$gtex_path        =  $config{ gtex_path };
$tcga_path        =  $config{ tcga_path };
$ubu_dir          =  $config{ ubu_dir } if ( exists $config{ ubu_dir } );
$gtex_sample_attr =  $config{ gtex_sample_attr };
$ensg_2_hugo_file =  $config{ ENSG_UCSC_common_genes };


# Check existence of required directories
( $tissue ) or die "ERROR: please provide tissue type\n";
(-e "$gtex_path/sra/$tissue") or die "ERROR: Non-existent directory $gtex_path/sra/$tissue\n";
(-e "$tcga_path/$tissue")     or die "ERROR: Non-existent directory $tcga_path/$tissue\n";

if ($tissue2){
    (-e "$tcga_path/$tissue2") or die "ERROR: Non-existent directory $tcga_path/$tissue2\n";
}

if ($tissue3){
    (-e "$tcga_path/$tissue3") or die "ERROR: Non-existent directory $tcga_path/$tissue3\n";
}


######################### Create data matrix #################################

my $work_dir = getcwd;
my (@samples, %expr, %gene, $data_matrix_fh);


# declare subroutines
#sub SplitFile($,$,$,$);
#sub ReadSampleExpression($);


my (%ens2hugo, %hugo2entrez);
map{chomp; my @f=split(/\t/); $ens2hugo{$f[0]}=$f[1]; $hugo2entrez{$f[1]}=$f[2]}`cat $ensg_2_hugo_file`;


my ( $prefix, $prefix2, $prefix3, $batch_str );
$prefix  = lc("$tissue-$normalization");
$prefix2 = lc("$tissue2-$normalization") if ($tissue2);
$prefix3 = lc("$tissue3-$normalization") if ($tissue3);


if (!-e "$prefix.txt" or ($normalization ne "count" and !-s "$prefix.adjusted.txt")){
    $data_matrix_fh = IO::File->new( "$prefix.txt", ">" ) or die "ERROR: Couldn't create file $prefix.txt\n";
    # Print header line
    # $data_matrix_fh->print("Gene\tDescription");
    $data_matrix_fh->print("Gene");
    
    # Read GTEx normal samples
    chdir "$gtex_path/sra";
    &ReadSampleExpression($tissue);
    
    my $n = scalar @samples;
    $batch_str = '';
    map{ $batch_str .= "1\n" }@samples;
    
    chdir $tcga_path;
    # Read TCGA normal samples
    &ReadSampleExpression($tissue);
    
    my $m = scalar @samples;
    map{ $batch_str .= "2\n" }@samples[($n+1)..$m];
    $n = $m;
    
    # Read TCGA tumor samples
    &ReadSampleExpression($tissue . '-t');

    $m = scalar @samples;
    map{ $batch_str .= "2\n" }@samples[($n+1)..$m];
    $n = $m;
    
    if ($tissue2){
        # Read TCGA normal samples
        &ReadSampleExpression($tissue2);
        
        $m = scalar @samples;
        map{ $batch_str .= "3\n" }@samples[($n+1)..$m];
        $n = $m;
        
        # Read TCGA tumor samples
        &ReadSampleExpression($tissue2 . '-t');
        
        $m = scalar @samples;
        map{ $batch_str .= "3\n" }@samples[($n+1)..$m];
        $n = $m;
        
        if ($tissue3){
            # Read TCGA normal samples
            &ReadSampleExpression($tissue3);
            
            $m = scalar @samples;
            map{ $batch_str .= "4\n" }@samples[($n+1)..$m];
            $n = $m;
            
            # Read TCGA tumor samples
            &ReadSampleExpression($tissue3 . '-t');
            
            $m = scalar @samples;
            map{ $batch_str .= "4\n" }@samples[($n+1)..$m];
            $n = $m;
        }
    }

    # Finish printing header line
    $data_matrix_fh->print("\n");
    
    # Write data matrix
    foreach my $g (keys %gene){
        next if(! exists $ens2hugo{$g});
        #$data_matrix_fh->print( $g );
        $data_matrix_fh->print( $ens2hugo{$g} );
        map{ $data_matrix_fh->print("\t" . $expr{$_}{$g}) if(defined $_) }@samples;
        $data_matrix_fh->print("\n");
    }
    $data_matrix_fh->close();

    ######################### Correct batch effect #################################
    
    chdir $work_dir;
    
    my $Rscript = "$FindBin::Bin/run-combat.R";
    
    my $batch_fh = IO::File->new( $prefix.'-combat-batch.txt', ">" ) or die "ERROR: Couldn't create file combat-batch.txt\n";
    $batch_fh->print( $batch_str. "\n");
    $batch_fh->close();
    
    if ($normalization ne "count"){
        (-s "$prefix.adjusted.txt") or `Rscript $Rscript $prefix.txt $prefix.adjusted`;
        (-s "$prefix.adjusted.txt") or die "ERROR: Failed to run $Rscript\n";
    }
    
    (!-e $prefix.'-combat-batch.txt') or `rm $prefix-combat-batch.txt`;
    
}

chdir $work_dir;


######################### Split output files #################################
my ($combined_file, $gtex_out_file, $tcga_out_file, $tcga_t_out_file, $tcga_out_file2, $tcga_t_out_file2, $tcga_out_file3, $tcga_t_out_file3);

if ($normalization eq "count"){
    $combined_file    = "$prefix.txt";
    $gtex_out_file    = "$prefix.gtex.txt";
    $tcga_out_file    = "$prefix.tcga.txt";
    $tcga_t_out_file  = "$prefix.tcga-t.txt";
    $tcga_out_file2   = "$prefix2.tcga.txt"    if ($tissue2);
    $tcga_t_out_file2 = "$prefix2.tcga-t.txt"  if ($tissue2);
    $tcga_out_file3   = "$prefix3.tcga.txt"    if ($tissue3);
    $tcga_t_out_file3 = "$prefix3.tcga-t.txt"  if ($tissue3);
}else{
    $combined_file    = "$prefix.adjusted.txt";
    $gtex_out_file    = "$prefix.adjusted.gtex.txt";
    $tcga_out_file    = "$prefix.adjusted.tcga.txt";
    $tcga_t_out_file  = "$prefix.adjusted.tcga-t.txt";
    $tcga_out_file2   = "$prefix2.adjusted.tcga.txt"   if ($tissue2);
    $tcga_t_out_file2 = "$prefix2.adjusted.tcga-t.txt" if ($tissue2);
    $tcga_out_file3   = "$prefix3.adjusted.tcga.txt"   if ($tissue3);
    $tcga_t_out_file3 = "$prefix3.adjusted.tcga-t.txt" if ($tissue3);
}


if (!-s $gtex_out_file or !-s $tcga_out_file or !-s $tcga_t_out_file){
    print "Splitting file\n";
    my $gtex_header_str    = "Hugo_Symbol\tEntrez_Gene_Id";
    my $tcga_header_str    = "Hugo_Symbol\tEntrez_Gene_Id";
    my $tcga_t_header_str  = "Hugo_Symbol\tEntrez_Gene_Id";
    my $tcga_header_str2   = "Hugo_Symbol\tEntrez_Gene_Id";
    my $tcga_t_header_str2 = "Hugo_Symbol\tEntrez_Gene_Id";
    my $tcga_header_str3   = "Hugo_Symbol\tEntrez_Gene_Id";
    my $tcga_t_header_str3 = "Hugo_Symbol\tEntrez_Gene_Id";
    
    my $header = `head -1 $combined_file`;
    chomp $header;

    my ($col, $gtex_col, $tcga_col, $tcga_t_col, $tcga_col2, $tcga_t_col2, $tcga_col3, $tcga_t_col3) = (0, 0, 0, 0, 0, 0, 0, 0);
    foreach (split(/\t/, $header)){
        $col++;
        next if (!$_);
        if (/$tissue/){
            if (/^GTEX/){
                s/.$tissue.//g;
                s/\./-/g;
                $gtex_header_str .= "\t$_";
                $gtex_col = $col;
            }elsif(/^TCGA/){
                if(/$tissue.t./){
                    s/.$tissue.t.//g;
                    s/\./-/g;
                    $tcga_t_header_str .= "\t$_";
                    $tcga_t_col = $col;
                }else{
                    s/.$tissue.//g;
                    s/\./-/g;
                    $tcga_header_str .= "\t$_";
                    $tcga_col = $col;
                }
            }
            next;
        }
        if( $tissue2 and /$tissue2/ ){
            if(/$tissue2.t./){
                s/.$tissue2.t.//g;
                s/\./-/g;
                $tcga_t_header_str2 .= "\t$_";
                $tcga_t_col2 = $col;
            }else{
                s/.$tissue2.//g;
                s/\./-/g;
                $tcga_header_str2 .= "\t$_";
                $tcga_col2 = $col;
            }
            next;
        }
        
        if( $tissue3 and /$tissue3/ ){
            if(/$tissue3.t./){
                s/.$tissue3.t.//g;
                s/\./-/g;
                $tcga_t_header_str3 .= "\t$_";
                $tcga_t_col3 = $col;
            }else{
                s/.$tissue3.//g;
                s/\./-/g;
                $tcga_header_str3 .= "\t$_";
                $tcga_col3 = $col;
            }
        }
    }
    
    &SplitFile( $combined_file, $gtex_out_file,   $gtex_header_str,   "1-$gtex_col" );
    $col = $gtex_col + 1;
    &SplitFile( $combined_file, $tcga_out_file,   $tcga_header_str,   "1,$col-$tcga_col" );
    $col = $tcga_col + 1;
    &SplitFile( $combined_file, $tcga_t_out_file, $tcga_t_header_str, "1,$col-$tcga_t_col" );
    if($tissue2){
        $col = $tcga_t_col + 1;
        &SplitFile( $combined_file, $tcga_out_file2,   $tcga_header_str2,   "1,$col-$tcga_col2" );
        $col = $tcga_col2 + 1;
        &SplitFile( $combined_file, $tcga_t_out_file2, $tcga_t_header_str2, "1,$col-$tcga_t_col2" );

        if($tissue3){
            $col = $tcga_t_col2 + 1;
            &SplitFile( $combined_file, $tcga_out_file3,   $tcga_header_str3,   "1,$col-$tcga_col3" );
            $col = $tcga_col3 + 1;
            &SplitFile( $combined_file, $tcga_t_out_file3, $tcga_t_header_str3, "1,$col-$tcga_t_col3" );
        }
    }
}

sub SplitFile() {
    my $in_file_name    = shift;
    my $out_file_name   = shift;
    my $file_header_str = shift;
    my $cols            = shift;
    
    my $file_handle = IO::File->new( $out_file_name, ">" ) or die "ERROR: Couldn't create file $out_file_name\n";
    $file_handle->print( $file_header_str. "\n");
    
    #`tail -n +2 $in_file_name | cut -f $cols >> $out_file_name`;
    foreach(`tail -n +2 $in_file_name | cut -f $cols`){
        chomp;
        my @F=split(/\t/);
        $F[0] .= "\t".((defined $hugo2entrez{$F[0]}) ? $hugo2entrez{$F[0]} : 0);
        $file_handle->print(join("\t",@F)."\n");
    }
    $file_handle->close();
}

sub CreateDataMatrix() {
    my $tissue = shift;

    my %barcode_hash = ();
    my ($idx, $barcode_idx, $run_idx) = (0, -1, -1);
    if(-e "$tissue/SraRunTable.txt"){
        map{$idx++; $run_idx = $idx if($_ eq "Run_s"); $barcode_idx = $idx if($_ eq "Sample_Name_s")}split(/\t/, `head $tissue/SraRunTable.txt | grep ^Assay_Type_s`);
        ($barcode_idx != -1 and $run_idx != -1) or die "ERROR: Unknown file format: SraRunTable.txt\n";
        map{chomp;my @f=split(/\t/,$_);$barcode_hash{$f[0]}=$f[1]}`grep -v ^Assay_Type_s $tissue/SraRunTable.txt | cut -f $run_idx,$barcode_idx`;
    }elsif(-e "$tissue/summary.tsv") {
        map{$idx++; $barcode_idx = $idx if($_ eq "barcode"); $run_idx = $idx if($_ eq "analysis_id")}split(/\t/, `head $tissue/summary.tsv | grep ^study`);
        ($barcode_idx != -1 and $run_idx != -1) or die "ERROR: Unknown file format: summary.tsv\n";
        map{chomp;my @f=split(/\t/,$_);$barcode_hash{$f[1]}=$f[0]}`grep -v ^study $tissue/summary.tsv | cut -f $barcode_idx,$run_idx`;
    }else{
        $data_matrix_fh->close();
        `rm -f "$prefix.txt"`;
        die "ERROR: File SraRunTable.txt or summary.tsv do not exist\n";
    }

    my $column = -1;
    #foreach(`find $tissue -name Quant.genes.results`){
    foreach(`cut -f 1 $tissue/QC/filtered_samples.txt`){
        chomp;
        if($normalization eq "featurecounts"){
            $_ = "$tissue/$_/ubu-quan/fcounts.fpkm.normalized_results";
        }elsif($normalization eq "rsem"){
            $_ = "$tissue/$_/ubu-quan/rsem.genes.normalized_results";
        }else{
            $_ = "$tissue/$_/Quant.genes.results";
        }
        
        if(!-e $_ or !-s $_) {
            print "Warning: Skip $_\n";
            next;
        }

        my @fields  = split(/\//);
        my $barcode = $barcode_hash{$fields[1]};
        
        my @id = split(/\//);
        push @samples, $id[1];
        $data_matrix_fh->print("\t$barcode($id[0])");
        
        my $n=0;
        if($normalization eq "rsem" or $normalization eq "featurecounts"){
            map{chomp; my @e=split(/\t/); $gene{$e[0]}=1; $expr{$id[1]}{$e[0]}=$e[1];$n++}`cat $_`;
        }else{
            if ($column == -1){
                my $idx = 0;
                my $col_name = "TPM";
                $col_name = "expected_count" if ($normalization eq "count");
                map{ $idx++; $column = $idx if($_ eq $col_name) }split(/\t/, `head $_ | grep ^gene_id`);
                ($column != -1) or die "Unknown file format: $_\n";
            }
            map{chomp; my @e=split(/\t/); $gene{$e[0]}=1; $expr{$id[1]}{$e[0]}=$e[1];$n++}`cut -f 1,$column $_ | tail -n +2`;
        }
        print "$fields[1]\t$barcode($id[0])\t$n\n";
    }
}

sub ReadSampleExpression() {
    my $tissue = shift;
    
    my %barcode_hash = ();
    my ($idx, $barcode_idx, $run_idx) = (0, -1, -1);
    if(-e "$tissue/SraRunTable.txt"){
        map{$idx++; $run_idx = $idx if($_ eq "Run_s"); $barcode_idx = $idx if($_ eq "Sample_Name_s")}split(/\t/, `head $tissue/SraRunTable.txt | grep ^Assay_Type_s`);
        ($barcode_idx != -1 and $run_idx != -1) or die "ERROR: Unknown file format: SraRunTable.txt\n";
        map{chomp;my @f=split(/\t/,$_);$barcode_hash{$f[0]}=$f[1]}`grep -v ^Assay_Type_s $tissue/SraRunTable.txt | cut -f $run_idx,$barcode_idx`;
    }elsif(-e "$tissue/summary.tsv") {
        map{$idx++; $barcode_idx = $idx if($_ eq "barcode"); $run_idx = $idx if($_ eq "analysis_id")}split(/\t/, `head $tissue/summary.tsv | grep ^study`);
        ($barcode_idx != -1 and $run_idx != -1) or die "ERROR: Unknown file format: summary.tsv\n";
        map{chomp;my @f=split(/\t/,$_);$barcode_hash{$f[1]}=$f[0]}`grep -v ^study $tissue/summary.tsv | cut -f $barcode_idx,$run_idx`;
    }else{
        $data_matrix_fh->close();
        `rm -f "$prefix.txt"`;
        die "ERROR: File SraRunTable.txt or summary.tsv do not exist\n";
    }
    
    #foreach(`find $tissue -name Quant.genes.results`){
    foreach my $line (`cut -f 1 $tissue/QC/filtered_samples.txt`){
        chomp $line;
        if(!-e "$tissue/$line/fcounts.fpkm" or !-s "$tissue/$line/Quant.genes.results") {
            print "Warning: Skip $line\n";
            next;
        }

        my $quan_file;
        if($normalization eq "tpm" or $normalization eq "count"){
            $quan_file = "$tissue/$line/Quant.genes.results";
        }else{
            if($normalization eq "featurecounts"){
                $quan_file = "$tissue/$line/fcounts.fpkm";
            }elsif($normalization eq "rsem"){
                $quan_file = "$tissue/$line/Quant.genes.results";
                my $column = -1;
                my $idx = 0;
                map{ $idx++; $column = $idx if($_ eq "FPKM") }split(/\t/, `head $quan_file | grep ^gene_id`);
                ($column != -1) or die "Unknown file format: $_\n";
                `grep -v gene_id $quan_file | cut -f 1,$column > $tissue/$line/ubu-quan/temp0.txt`;
                $quan_file = "$tissue/$line/ubu-quan/temp0.txt";
            }
            my $tmp_header = IO::File->new( "$tissue/$line/ubu-quan/temp1.txt", ">" ) or die "ERROR: Couldn't create file $tissue/$line/ubu-quan/temp1.txt\n";
            foreach(`cat $quan_file`){
                my @e=split(/\t/, $_);
                next if(! exists $ens2hugo{$e[0]});
                $tmp_header->print($_);
            }
            $tmp_header->close;
            $quan_file = "$tissue/$line/ubu-quan/temp1.txt";
            `perl $ubu_dir/perl/quartile_norm.pl -c 2 -q 75 -t 1000 -o $tissue/$line/ubu-quan/temp2.txt $quan_file`;
            $quan_file = "$tissue/$line/ubu-quan/temp2.txt";
        }
        
        my $barcode = $barcode_hash{$line};
        push @samples, $line;
        $data_matrix_fh->print("\t$barcode($tissue)");
        
        my $n=0;
        if($normalization eq "rsem" or $normalization eq "featurecounts"){
            map{chomp; my @e=split(/\t/); $gene{$e[0]}=1; $expr{$line}{$e[0]}=$e[1];$n++}`cat $quan_file`;
            `rm $tissue/$line/ubu-quan/temp*.txt`;
        }else{
            my $column = -1;
            my $idx = 0;
            my $col_name = "TPM";
            $col_name = "expected_count" if ($normalization eq "count");
            map{ $idx++; $column = $idx if($_ eq $col_name) }split(/\t/, `head $quan_file | grep ^gene_id`);
            ($column != -1) or die "Unknown file format: $_\n";
            map{chomp; my @e=split(/\t/); $gene{$e[0]}=1; $expr{$line}{$e[0]}=$e[1];$n++}`cut -f 1,$column $quan_file | tail -n +2`;
        }
        print "$line\t$barcode($tissue)\t$n\n";
    }
}

