#!/usr/bin/perl -w

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use IO::File;

(-e 'unc_hg19.gtf')               or die "ERROR: could not find file unc_hg19.gtf\n";
(-e 'gencode.v19.annotation.gtf') or die "ERROR: could not find file gencode.v19.annotation.gtf\n";
(-e 'unc_knownToLocus.txt')       or die "ERROR: could not find file unc_knownToLocus.txt\n";


# Read known genes used by UNC
my %unc_gene = ();
my %hugo2entrez = ();
print "Reading file unc_knownToLocus.txt\n";
map{chomp; my @cols=split(/\t/); my @f=split(/\|/, $cols[0]); $unc_gene{$f[0]}=$cols[1]; $hugo2entrez{$f[0]}=$f[1]}`cat unc_knownToLocus.txt`;

#`awk '{if($3=="gene")print $10"\t"$18"\t"$14}' gencode.v19.annotation.gtf | sed 's/["|;]//g' > gencode.v19.gene.txt`
#`awk '{if($3=="transcript")print $12"\t"$18"\t"$14}' gencode.v19.annotation.gtf | sed 's/["|;]//g' > gencode.v19.transcript.txt`

# Find uniq mapping between ENSEMBL and UNC genes
my %unc_2_ensg = ();
print "Processing file gencode.v19.gene.txt\n";
foreach(`cat gencode.v19.gene.txt`){
    chomp;
    my @fields=split(/\t/);
    next if( !exists $unc_gene{$fields[1]} );
    if ( !exists $unc_2_ensg{ $fields[1] } ){
        $unc_2_ensg{ $fields[1] } = $fields[0];
    }else{
        my $ret = `grep '$unc_gene{$fields[1]}' unc_hg19.gtf`;
        my @data = split(/\t/, $ret);
        my ($chr, $pos) = ($data[0], $data[3]);
        
        $ret = `awk '\$3==\"gene\"' gencode.v19.annotation.gtf | grep '$fields[0]'`;
        @data = split(/\t/, $ret);
        if ($chr eq $data[0] and $pos eq $data[3]){
            $unc_2_ensg{ $fields[1] } = $fields[0];
        }
    }
}

# uniquely map ENSEMBL genes to UNC genes
my %ensg_2_unc = ();
foreach my $key (keys %unc_2_ensg){
    $ensg_2_unc{ $unc_2_ensg{$key} } = $key;
}

# Write (ENSEMBL gene, hugo gene) pairs to file gencode.v19_ucsc.known.txt
my $gene_fh = IO::File->new( "gencode.v19_ucsc.known.txt", ">" ) or die "ERROR: Couldn't create file gencode.v19_ucsc.known.txt\n";
print "Writing file gencode.v19_ucsc.known.txt\n";
foreach(`cat gencode.v19.gene.txt`){
    chomp;
    my @f=split(/\t/);
    if(exists $ensg_2_unc{$f[0]}){
        if(!exists $hugo2entrez{$f[1]}){
            $gene_fh->print( "$_\t0\n" );
        }else{
            $gene_fh->print( "$_\t$hugo2entrez{$f[1]}\n" );
        }
    }
}
$gene_fh->close();
