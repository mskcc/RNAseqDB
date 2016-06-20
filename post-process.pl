#!/usr/bin/perl -w
# The program creates data matrix and runs combat to remove batch effects

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use Cwd;
use FindBin;
use lib "$FindBin::Bin";
use IO::File;
use File::Temp qw( tempdir );

######################### process input parameters #########################
my @usage;
push @usage, "\nUsage:  post-process.pl -i tissueType [other options]\n\n";
push @usage, "Options:\n";
push @usage, "  -h | --help             Displays this information\n";
push @usage, "  -i | --input-tissue     Input tissue type\n";
push @usage, "  -c | --tissue-conf      Input tissue configuration file (default: path of script file)\n";
push @usage, "  -t | --quan-tool        Expression quantification tool: RSEM (default) | FeatureCounts\n";
push @usage, "  -u | --quan-unit        Unit to measure gene expression: TPM | Count | FPKM (default)\n";
push @usage, "  -g | --transcript-type  Output data types: gene (default) | transcript\n";
push @usage, "  -p | --protein-coding   Output protein-coding genes/transcripts only if specified\n";
push @usage, "  -n | --hugo-gene-name   Output hugo gene names instead of ensembl gene ids\n";
push @usage, "  -r | --run-combat       Run combat to correct batch biases\n\n";


my ($help, $tissue, $tissue_conf, $quan_tool, $quan_unit, $gene_type);
my ($protein_coding, $hugo_gene_name, $run_combat) = (0, 0, 0);

GetOptions
(
 'h|help|?'             => \$help,
 'i|input-tissue=s'     => \$tissue,
 'c|tissue-conf=s'      => \$tissue_conf,
 't|quan-tool=s'        => \$quan_tool,
 'u|quan-unit=s'        => \$quan_unit,
 'g|transcript-type=s'  => \$gene_type,
 'p|protein-coding'     => \$protein_coding,
 'n|hugo-gene-name'     => \$hugo_gene_name,
 'r|run-combat'         => \$run_combat,
);

if (!defined $tissue) {
    print "ERROR: please provide tissue type\n";
    print @usage;
    die;
}

(defined $tissue_conf) or $tissue_conf = "$FindBin::Bin/tissue-conf.txt";
if (!-e $tissue_conf){
    print "ERROR: Cannot find tissue configuration file\n";
    print @usage;
    die;
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
        die;    }
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

if ($help) {
   print @usage;
   exit(0);
}

######################### Read configuration file #################################

# Read configuration file
my $config_file = "$FindBin::Bin/config.txt";
( -e $config_file ) or die "ERROR: The configuration file $config_file does not exist\n";
my %config;
map{ chomp; /^\s*([^=\s]+)\s*=\s*(.*)$/; $config{$1} = $2 if (defined $1 && defined $2) } `egrep -v \"^#\" $config_file`;
# Use configuration file to initialize variables
my ($gtex_path, $tcga_path, $ubu_dir, $gene_name_file);
$gtex_path        =  $config{ gtex_path };
$tcga_path        =  $config{ tcga_path };
$gene_name_file   =  $config{ ENSG_UCSC_common_genes };
$ubu_dir          =  $config{ ubu_dir } if ( exists $config{ ubu_dir } );


# Read tissue configuration file
my @lines = `grep $tissue $tissue_conf`;
(scalar @lines > 0) or die "ERROR: did not find $tissue in $tissue_conf\n";

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
my $flag = 0;
my $prefix;
foreach my $line (`grep -v ^# $tissue_conf`){
    chomp $line;
    next if(!$line);
    
    my @data = split(/\t/, $line);
    next if ($data[0] != $cluster);
    
    (-e "$gtex_path/sra/$data[1]" or -e "$tcga_path/$data[3]" or -e "$tcga_path/$data[3]-t") or die "ERROR: Nonexistent paths for $data[1]/$data[3]\n";
    
    if (-e "$gtex_path/sra/$data[1]/SraRunTable.txt" and -e "$tcga_path/$data[3]/summary.tsv"){
        $flag = 1;
        $prefix = "$data[1]-$data[3]-$quan_tool-$quan_unit";
    }
    push @lines, $line;
}

($flag == 1) or die "ERROR: There should be >=1 tissue with both GTEx and TCGA normals\n";


######################### Get genes of interest #################################

my (%ens2hugo, %hugo2entrez);
if ($hugo_gene_name){
    map{chomp; my @f=split(/\t/); $ens2hugo{$f[0]}=$f[1]; $hugo2entrez{$f[1]}=$f[2]}`cat $gene_name_file`;
}else{
    if ($gene_type eq 'transcript'){
        $gene_name_file =~ s/_ucsc\.known/\.transcript/;
    }else{
        $gene_name_file =~ s/_ucsc\.known/\.gene/;
    }
    
    if ( $protein_coding ){
        map{chomp; my @f=split(/\t/); $ens2hugo{$f[0]}=$f[1]}`grep protein_coding $gene_name_file | cut -f 1-2`;
    }else{
        map{chomp; my @f=split(/\t/); $ens2hugo{$f[0]}=$f[1]}`cut -f 1-2 $gene_name_file`;
    }
}

######################### Create data matrix #################################

my $batch_str = '';
my (%col_ranges, @samples, %expr, %gene);

my $work_dir = getcwd;
#if (!-e "$prefix.txt" or ($quan_unit ne "count" and !-s "$prefix.adjusted.txt")){
my $data_matrix_fh = IO::File->new( "$prefix.txt", ">" ) or die "ERROR: Couldn't create file $prefix.txt\n";
# Print header line
# $data_matrix_fh->print("Gene\tDescription");
$data_matrix_fh->print("Gene");

my %finished_tissues;
foreach my $line (@lines){
    my @data = split(/\t/, $line);
    my $tissue_type = '';
    $tissue_type = "\t$data[5]" if (scalar @data == 6 and defined $data[5]);
    
    if (-e "$gtex_path/sra/$data[1]" and !defined $finished_tissues{ "$gtex_path/sra/$data[1]" } and "$gtex_path/sra/$data[1]/SraRunTable.txt"){
        $finished_tissues{ "$gtex_path/sra/$data[1]" } = 1;
        chdir "$gtex_path/sra";
        # Read GTEx normal samples
        my $n1 = scalar @samples;
        &ReadSampleExpression($data[1]);
        my $n2 = scalar @samples;
        $col_ranges {"$data[1]-$quan_tool-$quan_unit-gtex.txt"} = ($n1+2).'-'.($n2+1);
        
        $n2 -= $n1;
        map{ $batch_str .= "normal\t$data[2]$tissue_type\n" }(1..$n2);
    }
    
    if (-e "$tcga_path/$data[3]" and !defined $finished_tissues{ "$tcga_path/$data[3]" } and -e "$tcga_path/$data[3]/summary.tsv"){
        $finished_tissues{ "$tcga_path/$data[3]" } = 1;

        chdir $tcga_path;
        # Read TCGA normal samples
        my $n1 = scalar @samples;
        &ReadSampleExpression($data[3]);
        my $n2 = scalar @samples;
        $col_ranges {"$data[3]-$quan_tool-$quan_unit-tcga.txt"} = ($n1+2).'-'.($n2+1);
        
        $n2 -= $n1;
        map{ $batch_str .= "normal\t$data[4]$tissue_type\n" }(1..$n2);
    }
    
    next if( scalar @data == 6 and defined $data[5] and $data[5] eq 'control' );
    
    if (-e "$tcga_path/$data[3]-t" and !defined $finished_tissues{ "$tcga_path/$data[3]-t" } and -e "$tcga_path/$data[3]-t/summary.tsv"){
        $finished_tissues{ "$tcga_path/$data[3]-t" } = 1;

        chdir $tcga_path;
        # Read TCGA tumor samples
        my $n1 = scalar @samples;
        &ReadSampleExpression("$data[3]-t");
        my $n2 = scalar @samples;
        $col_ranges {"$data[3]-$quan_tool-$quan_unit-tcga-t.txt"} = ($n1+2).'-'.($n2+1);
        
        $n2 -= $n1;
        map{ $batch_str .= "tumor\t$data[4]$tissue_type\n" }(1..$n2);
    }
}

# Finish printing header line
$data_matrix_fh->print("\n");

# Write data matrix
foreach my $g (sort {$gene{$a} <=> $gene{$b}} (keys %gene)){
    if ( $hugo_gene_name ){
        next if(!exists $ens2hugo{$g});
        $data_matrix_fh->print( $ens2hugo{$g} );
    }else{
        $data_matrix_fh->print( $g );
    }
    map{ $data_matrix_fh->print("\t" . $expr{$_}{$g}) if(defined $_) }@samples;
    $data_matrix_fh->print("\n");
}
$data_matrix_fh->close();


######################### Correct batch effect #################################

chdir $work_dir;

my $batch_fh = IO::File->new( $prefix.'-combat-batch.txt', ">" ) or die "ERROR: Couldn't create file combat-batch.txt\n";
$batch_fh->print( $batch_str. "\n");
$batch_fh->close();

my $Rscript = "$FindBin::Bin/run-combat.R";

if ($run_combat){
    (-s "$prefix.adjusted.txt") or `Rscript $Rscript $prefix.txt $prefix.adjusted`;
    (-s "$prefix.adjusted.txt") or die "ERROR: Failed to run $Rscript\n";
}

#(!-e $prefix.'-combat-batch.txt') or `rm $prefix-combat-batch.txt`;
#}

######################### Split output files #################################

print "Splitting file\n";

foreach my $outfile (keys %col_ranges){
    if ($run_combat){
        &SplitFile( "$prefix.adjusted.txt", $outfile,  $col_ranges{ $outfile } );
    }else{
        &SplitFile( "$prefix.txt", $outfile,  $col_ranges{ $outfile } );
    }
}


sub SplitFile() {
    my $in_file_name    = shift;
    my $out_file_name   = shift;
    my $cols            = shift;
    
    #my @idx = split(/-/, $cols);
    #my $header_str = join( "\t", @samples[ ($idx[0]-2)..($idx[1]-2) ] );
    my $head_str = `head -1 $in_file_name | cut -f $cols`;
    chomp $head_str;
    
    #my @fields = split(/-/, $out_file_name);
    #my $tissue_name = $fields[0];
    #$tissue_name =~ s/-/\./g;
    #$head_str = join(/\t/, map{ s/.$tissue_name.//g; }split(/\t/, $head_str));
    
    my $file_handle = IO::File->new( $out_file_name, ">" ) or die "ERROR: Couldn't create file $out_file_name\n";
    if ( $hugo_gene_name ){
        $file_handle->print( "Hugo_Symbol\tEntrez_Gene_Id\t$head_str\n");
    }else{
        $file_handle->print( "Gene_Id\tHugo_Symbol\t$head_str\n");
    }
    
    #`tail -n +2 $in_file_name | cut -f $cols >> $out_file_name`;
    foreach(`tail -n +2 $in_file_name | cut -f 1,$cols`){
        chomp;
        my @F=split(/\t/);
        if ( $hugo_gene_name ){
            $F[0] .= "\t".((defined $hugo2entrez{$F[0]}) ? $hugo2entrez{$F[0]} : 0);
        }else{
            $F[0] .= "\t".((defined $ens2hugo{$F[0]}) ? $ens2hugo{$F[0]} : '');
        }
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
    foreach my $line (`cut -f 1 $tissue/QC/filtered_samples.txt`){
        chomp;
        if(!-s "$tissue/$line/fcounts.fpkm" or !-s "$tissue/$line/Quant.genes.results") {
            print "Warning: Skip $line\n";
            next;
        }
        
        my $quan_file;
        if($quan_unit eq "fpkm"){
            if($quan_tool eq "rsem"){
                $quan_file = "$tissue/$line/ubu-quan/rsem.genes.normalized_results";
            }else{
                $quan_file = "$tissue/$line/ubu-quan/fcounts.fpkm.normalized_results";
            }
        }else{
            if($quan_tool eq "rsem"){
                $quan_file = "$tissue/$line/Quant.genes.results";
            }else{
                $quan_file = "$tissue/$line/fcounts.".$quan_unit;
            }
        }

        my @fields  = split(/\//, $quan_file);
        my $barcode = $barcode_hash{$fields[1]};
        
        my @id = split(/\//, $quan_file);
        push @samples, $id[1];
        $data_matrix_fh->print("\t$barcode($id[0])");
        
        my $n=0;
        if($quan_tool eq "rsem" and $quan_unit ne "fpkm"){
            if ($column == -1){
                my $idx = 0;
                my $col_name = "TPM";
                $col_name = "expected_count" if ($quan_unit eq "count");
                map{ $idx++; $column = $idx if($quan_file eq $col_name) }split(/\t/, `head $quan_file | grep ^gene_id`);
                ($column != -1) or die "Unknown file format: $quan_file\n";
            }
            map{chomp; my @e=split(/\t/); $gene{$e[0]}=1; $expr{$id[1]}{$e[0]}=$e[1];$n++}`cut -f 1,$column $quan_file | tail -n +2`;
        }else{
            map{chomp; my @e=split(/\t/); $gene{$e[0]}=1; $expr{$id[1]}{$e[0]}=$e[1];$n++}`cat $quan_file`;
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
    
    my $tmp_dir = tempdir( CLEANUP => 1 );
    #foreach(`find $tissue -name Quant.genes.results`){
    foreach my $line (`cut -f 1 $tissue/QC/filtered_samples.txt`){
        chomp $line;

        my $quan_file;
        if($quan_tool eq "rsem"){
            $quan_file = "$tissue/$line/Quant.genes.results";
        }else{
            $quan_file = "$tissue/$line/fcounts.".$quan_unit;
        }

        if(!-s $quan_file) {
            print "Warning: Skip $line\n";
            next;
        }

        if($quan_unit eq 'fpkm'){          # Do quantile normalization for fpkm
            if ($gene_type eq 'gene'){ # gene expression
                if( !$hugo_gene_name and !$protein_coding ){
                    if($quan_tool eq "rsem"){
                        $quan_file = "$tissue/$line/ubu-quan/rsem.genes.normalized_results";
                    }else{
                        $quan_file = "$tissue/$line/ubu-quan/fcounts.fpkm.normalized_results";
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
                    $quan_file = "$tissue/$line/Quant.isoforms.results";
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
        
        my $barcode = $barcode_hash{$line};
        push @samples, $line;
        #$data_matrix_fh->print("\t$barcode($tissue)");
        $data_matrix_fh->print("\t$barcode");
        
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
            map{chomp; my @e=split(/\t/); if(defined $ens2hugo{$e[0]}){$gene{$e[0]}=1; $expr{$line}{$e[0]}=$e[1];$n++}}`cut -f 1,$column $quan_file | tail -n +2`;
        }else{
            map{chomp; my @e=split(/\t/); if(defined $ens2hugo{$e[0]}){$gene{$e[0]}=1; $expr{$line}{$e[0]}=$e[1];$n++}}`cat $quan_file`;
        }
        print "$line\t$barcode($tissue)\t$n\n";
    }
}

