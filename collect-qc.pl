#!/usr/bin/perl -w

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
use Cwd 'abs_path';


my @usage;
push @usage, "\n-----------------------------------------------------------\n\n";
push @usage, "Usage:  collect-qc.pl --input file\n\n";
push @usage, "Options:\n";
push @usage, "  -h | --help         Display this information\n";
push @usage, "  -c | --config       Configuration file\n";
push @usage, "  -i | --input        Input directory\n";
push @usage, "  -m | --mRIN         mRIN cutoff(default: -0.11)\n";
push @usage, "  -d | --dataMatrix   Type of data matrix: TPM, Count, fcount, FPKM\n";
push @usage, "--------------------------------------------------------------\n\n";


my $help;
my $input;
my $config_file;


my $alignment_cutoff       = 40;    # Percentage of aligned reads
my $assigned_reads_cutoff  = 0;     # Ratio of reads assigned by FeatureCount. 0 means keep all samples
my $mRIN_score_cutoff      = -0.11;
my $out_matrix;

############################## Input parameters ##############################
# Read input parameters

GetOptions
(
    'h|help|?'       => \$help,
    'c|config=s'     => \$config_file,
    'i|input=s'      => \$input,
    'm|mRIN=s'       => \$mRIN_score_cutoff,
    'd|dataMatrix=s' => \$out_matrix,
);

if ($help) {
    print @usage;
    exit(0);
}

if ($input) {
    if (!-e $input){
        print "\nError: $input file does not exist!\n";
        print @usage;
        exit;
    }
}else{
    print @usage;
    exit(0);
}

if (defined $out_matrix){
    $out_matrix = uc($out_matrix);
    if ($out_matrix ne "TPM" and $out_matrix ne 'FPKM' and $out_matrix ne 'COUNT' and $out_matrix ne 'FCOUNT'){
        print "\nError: unknown data type for creating dataMatrix!\n";
        print @usage;
        exit;
    }
}

############################# Read configuration file ############################

(defined $config_file) or $config_file = "$FindBin::Bin/config.txt";
( -e $config_file ) or die "ERROR: The configuration file $config_file does not exist\n";
my %config;
map{ chomp; /^\s*([^=\s]+)\s*=\s*(.*)$/; $config{$1} = $2 if (defined $1 && defined $2) } `egrep -v \"^#\" $config_file`;

# Use the configuration file to initialize variables
my ( $thread_n, $rseqc_dir, $sra_dir, $fastqc_bin, $mRIN_dir, $gtex_path, $gtex_sample_attr, $ensg_2_hugo_file);

$thread_n               = 1;
$thread_n               = $config{ thread_no }              if ( exists $config{ thread_no  } );
$rseqc_dir              = $config{ rseqc_dir }              if ( exists $config{ rseqc_dir } );
$sra_dir                = $config{ sra_dir }                if ( exists $config{ sra_dir } );
$fastqc_bin             = $config{ fastqc_bin };
$mRIN_dir               = $config{ mRIN_dir }               if ( exists $config{ mRIN_dir } );
$gtex_path              = $config{ gtex_path }              if ( exists $config{ gtex_path } );
$gtex_sample_attr       = $config{ gtex_sample_attr }       if ( exists $config{ gtex_sample_attr } );
$ensg_2_hugo_file       = $config{ ENSG_UCSC_common_genes };

###################################################################################

`mkdir -p $input/QC`;
chdir "$input/QC";


# Map sample ID to barcode
my (%barcode_hash, %filtered_sample, %filtering_sum);
my $work_dir     = abs_path(getcwd);
$gtex_path       = abs_path($gtex_path);
my $gtex_dir_len = length($gtex_path);

# GTEx barcode
if ($gtex_path eq substr($work_dir, 0, $gtex_dir_len) or lc($work_dir) =~ /dbgap-8936/) {
    # Determine the column Sample_Name_s in
    my ($idx, $barcode_idx, $run_idx) = (0, -1, -1);
    map{$idx++; $run_idx = $idx if($_ eq "Run_s"); $barcode_idx = $idx if($_ eq "Sample_Name_s")}split(/\t/, `head ../SraRunTable.txt | grep ^Assay_Type_s`);
    ($barcode_idx != -1 and $run_idx != -1) or die "Unknown file format: SraRunTable.txt\n";
    
    foreach(`grep -v ^Assay_Type_s ../SraRunTable.txt | cut -f $run_idx,$barcode_idx`){
        chomp;
        my @f = split(/\t/,$_);
        $barcode_hash{$f[0]} = $f[1];
        if (-s "../$f[0]/Log.final.out" and -s "../$f[0]/Quant.genes.results" and -s "../$f[0]/ubu-quan/fcounts.fpkm.normalized_results" and -s "../$f[0]/ks/sample.ks.txt" and -e "../$f[0]/Aligned.sortedByCoord.out_fastqc"){
            $filtered_sample{$f[0]}=1;
        }else{
            $filtered_sample{$f[0]}=0;
        }
    }
}

# TCGA barcode
if (-e "../summary.tsv"){
    my @lines = `head -2 ../summary.tsv`;
    if ($lines[0] =~ /^study/ and $lines[1] =~ /^TCGA/){
        my ($idx, $barcode_idx, $analysis_idx) = (0, -1, -1);
        map{$idx++; $barcode_idx = $idx if($_ eq "barcode"); $analysis_idx = $idx if($_ eq "analysis_id")}split(/\t/, `head ../summary.tsv | grep ^study`);
        ($barcode_idx != -1 and $analysis_idx != -1) or die "Unknown file format: summary.tsv\n";
        
        foreach(`grep -v ^study ../summary.tsv | cut -f $barcode_idx,$analysis_idx`){
            chomp;
            my @f = split(/\t/,$_);
            $barcode_hash{$f[1]} = $f[0];
            if (-s "../$f[1]/Log.final.out" and -s "../$f[1]/Quant.genes.results" and -s "../$f[1]/ubu-quan/fcounts.fpkm.normalized_results" and -s "../$f[1]/ks/sample.ks.txt" and -e "../$f[1]/Aligned.sortedByCoord.out_fastqc"){
                $filtered_sample{$f[1]}=1;
            }else{
                $filtered_sample{$f[1]}=0;
            }
        }
    }
}

# If neither GTEx nor TCGA
if (0 == scalar keys %filtered_sample){
    my @files = `find .. -name Log.final.out`;
    foreach( @files ){
        chomp;
        s/\/Log.final.out//;
        s/\.\.\///;
        $filtered_sample{$_}=1;
    }
}

# Filter GTEx data
if ($gtex_path eq substr($work_dir, 0, $gtex_dir_len) or lc($work_dir) =~ /dbgap-8936/) {
    # Read GTEx sample attribute file
    if (defined $gtex_path and defined $gtex_sample_attr and -e "$gtex_path/files/$gtex_sample_attr") {
        my ($idx, $sample_col_idx) = (0, -1);
        map{$idx++; $sample_col_idx = $idx if($_ eq "SAMPID")}split(/\t/, `head -20 $gtex_path/files/$gtex_sample_attr | grep ^dbGaP_Sample_ID`);
        ($sample_col_idx != -1) or die "Unknown file format: $gtex_path/files/$gtex_sample_attr\n";
        
        my %good_samples = map{ chomp; ($_,1) }`grep USE $gtex_path/files/$gtex_sample_attr|grep TrueSeq|cut -f $sample_col_idx`;
        
        print "\nSamples unused by GTEx:\n";
        foreach my $key (keys %filtered_sample){
            #my $barcode=`grep $key ../SraRunTable.txt | cut -f $barcode_idx`;
            my $barcode = $barcode_hash{$key};
            chomp $barcode;
            if(!exists $good_samples{$barcode}) {
                print "$key\n";
                $filtered_sample{ $key } = 0;
                $filtering_sum{ $key }{ gtex } = 'Unused';
            }
        }
    }
}

# Collect fastqc report
CollectFastQC();

# Collect alignment rate
(-e 'align_sum.txt') or SummarizeAlignment();

print "\nPoorly aligned samples:\n";
foreach(`cat align_sum.txt`){
    next if(/^#/);
    chomp;
    my @data = split(/\t/, $_);
    my $align_rate = $data[1];
    $align_rate =~ s/\%//;

    #next if ( $filtered_sample{ $key } == 0);
    
    if ($align_rate < $alignment_cutoff) {
        s/^\.\.\///;
        print "$_\n";
 
        my $id = $data[0];
        $id=~s/^\.\.\///;
        $id =~ s/^\s+|\s+$|\r|\n//g;
        $filtered_sample { $id } = 0;
        $filtering_sum{ $id }{ align } = $align_rate."%";
    }
}

# Collect FeatureCounts stat
(-e 'fcount_sum.txt') or FeatureCountsStat();

print "\nSamples with low ratio of assigned reads:\n";
foreach(`cat fcount_sum.txt`){
    next if(/^#/);
    chomp;
    my @data = split(/\t/, $_);
    my $assign_rate = $data[1];
    if ($assign_rate < $assigned_reads_cutoff) {
        s/^\.\.\///;
        print "$_\n";

        my $id = $data[0];
        $id =~ s/^\s+|\s+$|\r|\n//g;
        if( exists $filtered_sample{ $id } ){
            $filtered_sample { $id } = 0;
            $filtering_sum{ $id }{ fcount } = $assign_rate;
        }
    }
}

# Collect mRIN report
(-e 'out.mRIN.txt') or SummarizeMRIN();

if (-e 'out.mRIN.txt'){
    #$ret = `awk \047\{if(\$2<-0.11)print\}\047 out.mRIN.txt`;
    #if (defined $ret){
    #    print "\nDegraded samples:\n";
    #    print $ret;
    #}
    print "\nDegraded samples:\n";
    my $mRIN_col_idx = -1;
    foreach (`cat out.mRIN.txt`){
        if (/^Sample/ and /Pvalue$/){
            my $idx = 0;
            map{$mRIN_col_idx = $idx if($_ eq "mRIN"); $idx++}split(/\t/);
            ($mRIN_col_idx != -1) or die "Unknown file format: out.mRIN.txt\n";
            next;
        }
        chomp;
        my @data = split(/\t/, $_);
        next if ( $data[$mRIN_col_idx] >= $mRIN_score_cutoff );
        print "$_\n";

        my $id = $data[0];
        $id =~ s/\./-/g;
        if( exists $filtered_sample{ $id } ){
            $filtered_sample{ $id } = 0;
            $filtering_sum{ $id }{ mRIN } = $data[$mRIN_col_idx];
        }
        # mRIN modifies some sample IDs by adding prefix 'X'.
        elsif( exists $filtered_sample{ substr($id,1) } ){
            $filtered_sample{ substr($id,1) } = 0;
            $filtering_sum{ substr($id,1) }{ mRIN } = $data[$mRIN_col_idx];
        }
    }
}

my $r_fh = IO::File->new( 'filtered_samples.txt', ">" ) or die "ERROR: Cannot overwrite filtered_samples.txt\n";
my $f_fh = IO::File->new( 'filtering_sum.txt', ">" ) or die "ERROR: Cannot overwrite filtering_sum.txt\n";
$f_fh->print( "#Sample\tBarcode\tKept\tGTEx\tAlignment rate\(%\)\tFeatureCounts assigned read rate\tmRIN\n" );
my @sample_list;
foreach my $key (keys %filtered_sample){
    if ( defined $filtered_sample{ $key } ){
        $f_fh->print( $key."\t" );
        $f_fh->print( $barcode_hash{$key} )  if (exists $barcode_hash{$key});
        $f_fh->print( "\t" );
        if ( $filtered_sample{ $key } > 0 ){
            $r_fh->print( $key );
            $r_fh->print( "\t". $barcode_hash{$key} )  if (exists $barcode_hash{$key});
            $r_fh->print( "\n" );
            $f_fh->print( "Yes" );
            push (@sample_list, $key);
        }
        foreach ( qw(gtex  align  fcount  mRIN) ){
            if ( exists $filtering_sum{ $key }{ $_ } ){
                $f_fh->print( "\t$filtering_sum{ $key }{ $_ }" );
            }else{
                $f_fh->print( "\t" );
            }
        }
        $f_fh->print( "\n" );
    }
}
$r_fh->close;
$f_fh->close;


# Calculate gene coverage
# (-e 'rseqc.geneBodyCoverage.curves2.pdf') or GetGeneBodyCov(0, 'rseqc.geneBodyCoverage.curves2.pdf');
# (-e 'rseqc.geneBodyCoverage.curves.pdf')  or GetGeneBodyCov(1, 'rseqc.geneBodyCoverage.curves.pdf');

if (defined $out_matrix){
    ( @sample_list ) or die "ERROR: No file left after filtering\n";

    if ($out_matrix eq "FCOUNT"){
        if (scalar (@sample_list) <= 1) {
            print "Warning: " . scalar (@sample_list) . " data file; skip creating data matrix\n";
        }elsif( !-e 'data-matrix.txt.gz' ){
            CreateDataMatrix(\@sample_list, 'fcounts.count', 2, 'data-matrix.txt');
            `gzip 'data-matrix.txt'`;
        }
    }else{
        $out_matrix = 'expected_count' if ($out_matrix ne "TPM" and $out_matrix ne 'FPKM');
        my $file_name = "../$sample_list[0]/Quant.genes.results";
        ( -s $file_name ) or die "ERROR: file $file_name do not exist\n";
        
        my ($column, $idx) = (-1, 0);
        map{ $idx++; $column = $idx if($_ eq $out_matrix) }split(/\t/, `head $file_name | grep ^gene_id`);
        ($column != -1) or die "Unknown file format: $file_name\n";
        
        
        if (scalar (@sample_list) <= 1) {
            print "Warning: " . scalar (@sample_list) . " data file; skip creating data matrix\n";
        }elsif( !-e 'data-matrix.txt.gz' ){
            CreateDataMatrix(\@sample_list, 'Quant.genes.results', $column, 'data-matrix.txt');
            `gzip 'data-matrix.txt'`;
        }
    }
}

sub SummarizeMRIN {
    
    print "\nSummarize mRIN scores...\n\n";
    
    my $r_fh = IO::File->new( 'mRIN_samples.txt', ">" );
    my $n = 0;
    foreach my $file (`find ../ -name sample.ks.txt`){
        chomp $file;
        my @f  = split(/\//,$file);
        my $id = $f[1];
        $r_fh->print( "$file\t$id\n" );
        $n++;
    }
    $r_fh->close;
    
    return if($n==0);
    
    `perl $mRIN_dir/gen_ks_matrix.pl -v -base . --min-avg-cov 2 -v mRIN_samples.txt mRIN_samples.mat.txt`;
    `Rscript $mRIN_dir/cal_mrin.R -k mRIN_samples.mat.txt -m out.mRIN.txt -G out.GIS.txt -v` if(!-e 'out.mRIN.txt');
    #`Rscript $mRIN_dir/cal_mrin.R -k mRIN_samples.mat.txt -x mRIN_samples.RPKM.mat.txt -r 2 -s 0.05 -e 0.5 -b -m out.mRIN.txt -G out.GIS.txt -v`;
}

sub GetGeneBodyCov {
    my $filtering = shift;
    my $out_pdf   = shift;

    print "\nPlot gene body coverage...\n\n";

    my @r_files = `find .. -name rseqc.geneBodyCoverage.r`;
    
    my $d_matrix = 'data_matrix <- matrix(c(';
    my $r_label  = 'rowLabel <- c(';
    my @ids;
    
    my $r_fh = IO::File->new( 'rseqc.geneBodyCoverage2.r', ">" );
    foreach my $r_file ( @r_files ){
        next if ($r_file =~ /..\/QC/);
        my $id = $r_file;
        chomp $id;
        $id =~ s/\/rseqc.geneBodyCoverage.r//;
        $id = substr($id, 3);
        
        if($filtering and (!defined $filtered_sample{ $id } or $filtered_sample{ $id } == 0)){
            print "Skip $id\n";
            next;
        }
        
        $id =~ s/-/_/g;
        $id =  'S_' . $id;
        push @ids, $id;
        
        my $r_str = `grep \'Aligned.sortedByCoord.out <-\' $r_file`;
        $r_str =~ s/Aligned.sortedByCoord.out/$id/;
        $r_fh->print( $r_str );
        
        my $idx1 = index($r_str, '(');
        my $idx2 = index($r_str, ')');
        $r_str = substr($r_str, $idx1+1, $idx2-$idx1-1);
        
        my ($sum_50, $sum_100, $idx) = (0, 0, 0);
        map{ $sum_100+=$_; $sum_50+=$_ if($idx<50);$idx++ }split(/,/, $r_str);
        
        $d_matrix .= "$id,";
        $r_label  .= "\"$id\",";
    }
    
    $d_matrix =~ s/,$/\)/;
    $r_fh->print( $d_matrix, ", byrow=T, ncol=100)\n" );
    
    $r_label  =~ s/,$/\)/;
    $r_fh->print( $r_label, "\n" );
    
    $r_fh->print( "pdf(\"$out_pdf\")\n" );
    $r_fh->print( "x=1:100\n" );
    $r_fh->print( "icolor = colorRampPalette(c(\"#7fc97f\",\"#beaed4\",\"#fdc086\",\"#ffff99\",\"#386cb0\",\"#f0027f\"))(".scalar(@ids).")\n" );
    $r_fh->print( "layout(matrix(c(1,1,1,2,1,1,1,2,1,1,1,2), 4, 4, byrow = TRUE))\n" );
    $r_fh->print( "plot(x,". $ids[0] .",type='l',xlab=\"Gene body percentile (5'->3')\", ylab=\"Coverage\",lwd=0.8,col=icolor[1])\n" );
    foreach(2..scalar(@ids)){
        $r_fh->print( "lines(x,". $ids[$_-1] .",type='l',col=icolor[$_])\n" );
    }
    $r_fh->print( "par(mar=c(1,0,2,1))\n" );
    $r_fh->print( "plot.new()\n" );
    
    my $legend = "legend(0,1,fill=icolor[1:". scalar(@ids) ."],legend=c(";
    map {$legend .= "'$_',";} @ids;
    $legend  =~ s/,$/\)\)/;
    
    $r_fh->print( $legend, "\n" );
    $r_fh->print( "dev.off()\n" );
    $r_fh->close;
    
    `Rscript rseqc.geneBodyCoverage2.r`;
    
}

sub CollectFastQC {
    print "\nSummarizing FastQC results...\n\n";
    
    `echo -e '\tBasic Statistics\tPer base sequence quality\tPer tile sequence quality\tPer sequence quality scores\tPer base sequence content\tPer sequence GC content\tPer base N content\tSequence Length Distribution\tSequence Duplication Levels\tOverrepresented sequences\tAdapter Content\tKmer Content' > fastqc_sum.txt`;
    
    my @files = `find ../*/Aligned.sortedByCoord.out_fastqc -name summary.txt`;
    print "Sample that failed Per base sequence quality:\n";
    foreach( @files ){
        chomp;
        my $id = $_;
        $id =~ s/\/Aligned.sortedByCoord.out_fastqc\/summary.txt//;
        $id =  substr($id, 3);
        
        my $str = $id;
        my @lines = `cat $_`;
        foreach my $line (@lines) {
            chomp $line;
            next if ( length( $line )==0 );
            
            my @temp = split(/\t/, $line);
            print "$id\n" if ($temp[1] eq 'Per base sequence quality' and $temp[0] eq 'FAIL');
            $str .= "\t".$temp[0];
        }
        `echo \'$str\' >> fastqc_sum.txt`;
    }
}

sub SummarizeAlignment {
    print "\nSummarize rate of uniquely aligned reads...\n\n";
    
    `echo \'#Sample\tAlignment rate\tTotal number of reads\' > align_sum.txt`;
    
    foreach( keys %filtered_sample ){
        next if(!defined $_);
        if (!-s "../$_/Log.final.out"){
            if (exists $filtering_sum{ $_ }{ gtex }){
                warn "Error: $_ (unused by GTEx) was not aligned correctly\n";
            }else{
                warn "Error: $_ was not aligned correctly\n";
            }
            next;
        }
        if (-e "../$_/Quant.temp"){
            warn "Error: $_ was not quantified by RSEM correctly\n";
            next;
        }
        my @temp_fcount = glob("../$_/temp-sort-*.bin");
        if (@temp_fcount){
            warn "Error: $_ was not quantified by FeatureCounts correctly\n";
            next;
        }
        #my $str  = 'Uniquely mapped reads %';
        #my $line = `grep \'$str\' ../$_/Log.final.out`;
        #chomp $line;
        #$line =~ s/\s+Uniquely mapped reads . \|//;
        
        #$str  = 'Number of input reads';
        #my $line2 = `grep \'$str\' ../$_/Log.final.out`;
        #chomp $line2;
        #$line2 =~ s/\s+Number of input reads \|//;
        
        #$line = "$_$line$line2";
        #`echo \'$line\' >> align_sum.txt`;

        my $echo_str = $_;
        my ($map_rate, $total_read) = (0, 0);
        foreach my $line (`cat ../$_/Log.final.out`){
            chomp $line;
            if( $line =~ /\s+Uniquely mapped reads . \|/ ){
                $line =~ s/\s+Uniquely mapped reads . \|//;
                $line =~ s/\%//;
                $map_rate += $line;
            }
            if( $line =~  /\s+\% of reads mapped to multiple loci \|/ ){
                $line =~ s/\s+\% of reads mapped to multiple loci \|//;
                $line =~ s/\%//;
                $map_rate += $line;
            }
            if( $line =~  /\s+\% of reads mapped to too many loci \|/ ){
                $line =~ s/\s+\% of reads mapped to too many loci \|//;
                $line =~ s/\%//;
                $map_rate += $line;
            }
            if( $line =~ /\s+Number of input reads \|/ ){
                $line =~ s/\s+Number of input reads \|//;
                $total_read = $line;
            }
        }
        $map_rate = sprintf("%.2f", $map_rate);
        `echo \'$_\t$map_rate%$total_read\' >> align_sum.txt`;
    }
}

sub FeatureCountsStat {
    print "\nSummarize FeatureCounts stat...\n";
    
    `echo \'#Sample\tAssigned\tAmbiguity\tMultiMapping\tNoFeatures\' > fcount_sum.txt`;
    
    my @files = `find .. -name fcounts.stat`;
    foreach( @files ){
        chomp;
        my ($assigned, $ambiguity, $multimap, $noFeature, $total_reads) = (0,0,0,0,0);
        foreach(`cat $_`)
        {
            chomp;
            my @f=split(/\t/);
            $total_reads += $f[1];
            $assigned  = $f[1] if($f[0] eq "Assigned");
            $ambiguity = $f[1] if($f[0] eq "Unassigned_Ambiguity");
            $multimap  = $f[1] if($f[0] eq "Unassigned_MultiMapping");
            $noFeature = $f[1] if($f[0] eq "Unassigned_NoFeatures");
        }
        $assigned  = sprintf("%.2f", $assigned  / $total_reads);
        $ambiguity = sprintf("%.2f", $ambiguity / $total_reads);
        $multimap  = sprintf("%.2f", $multimap  / $total_reads);
        $noFeature = sprintf("%.2f", $noFeature / $total_reads);
        
        my $line = $_;
        $line =~ s/\/fcounts.stat/ /;
        $line =~ s/\.\.\///;
        `echo \'$line\t$assigned\t$ambiguity\t$multimap\t$noFeature\' >> fcount_sum.txt`;
    }
}


sub CreateDataMatrix {
    my $file_list  = shift;
    my $quant_file = shift;
    my $column     = shift;
    my $out_file   = shift;
    
    print "\nCreate data matrix $out_file...\n\n";
    
    $ensg_2_hugo_file =~ s/_ucsc\.known/\.gene/;
    my %ens2hugo;
    map{chomp; my @f=split(/\t/); $ens2hugo{$f[0]}=$f[1]}`cat $ensg_2_hugo_file`;

    # Create data matrix
    my (@samples, %data_matrix, %gene_list);
    foreach my $file (@{$file_list}){
        next if (!-s "../$file/$quant_file");
        push @samples, $file;
        if ($quant_file eq 'fcounts.count'){
            map{ chomp; my @e=split(/\t/); $data_matrix{$file}{$e[0]}=$e[1]; $gene_list{$e[0]}++}`cat ../$file/$quant_file`;
        }else{
            map{ chomp; my @e=split(/\t/); $data_matrix{$file}{$e[0]}=$e[1]; $gene_list{$e[0]}++}`cut -f 1,$column ../$file/$quant_file | tail -n +2`;
        }
    }
    
    # Store the matrix
    my $r_fh = IO::File->new( $out_file, ">" ) or die "ERROR: Failed to creast file $out_file\n";
    $r_fh->print( "Gene\tDescription\t". join("\t", map{exists $barcode_hash{$_} ? $barcode_hash{$_} : $_}@samples) . "\n" );
    foreach my $gene (keys %gene_list){
        $r_fh->print("$gene\t$ens2hugo{$gene}");
        foreach my $sample ( @samples ){
            my $expr = defined( $data_matrix{$sample}{$gene} ) ? $data_matrix{$sample}{$gene} : 0;
            $r_fh->print( "\t$expr" );
        }
        $r_fh->print( "\n" );
    }
    $r_fh->close;
}

