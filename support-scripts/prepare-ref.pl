#!/usr/bin/perl -w

use strict;
use warnings FATAL => 'all';
use IO::File;

my $gencode_file = $ARGV[0];

die "ERROR: Gencode file (GTF) not provided\n" if (!defined $gencode_file);

# The following code was adapted from http://www.gencodegenes.org/gencodeformat.html
my $gene_fh = IO::File->new( "gencode.v19.gene.txt", ">" ) or die "ERROR: Couldn't overwrite file gencode.v19.gene.txt\n";
my $transcript_fh = IO::File->new( "gencode.v19.transcript.txt", ">" ) or die "ERROR: Couldn't overwrite file gencode.v19.transcript.txt\n";

my %hugo_symbols;
my @gene_names;
open(IN, "<$gencode_file") or die "ERROR: Can't open $gencode_file.\n";
while(<IN>){
    next if(/^##/); #ignore header
    chomp;
    my %attribs = ();
    my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split("\t");
    #store nine columns in hash
    my %fields = (
        chr        => $chr,
        source     => $source,
        type       => $type,
        start      => $start,
        end        => $end,
        score      => $score,
        strand     => $strand,
        phase      => $phase,
        attributes => $attributes,
    );
    
    my @add_attributes = split(";", $attributes);
    # store ids and additional information in second hash
    foreach my $attr ( @add_attributes ) {
        next unless $attr =~ /^\s*(.+)\s(.+)$/;
        my $c_type  = $1;
        my $c_value = $2;
        $c_value =~ s/\"//g;
        if($c_type && $c_value){
            if(!exists($attribs{$c_type})){
                $attribs{$c_type} = [];
            }
            $attribs{$c_type} = $c_value;
        }
        $attribs{start} = $fields{start};
        $attribs{end} = $fields{end};
    }
    
    if ($type eq "gene"){
        $gene_fh->print( $attribs{gene_id}."\t".$attribs{gene_name}."\t".$attribs{gene_type}."\n" );
        next if($attribs{transcript_type} ne "protein_coding");
        
        if(!exists($hugo_symbols{$attribs{gene_name}})){
            $hugo_symbols{$attribs{gene_name}} = \%attribs;
            push @gene_names,$attribs{gene_name};
        }elsif ($hugo_symbols{$attribs{gene_name}}->{gene_status} ne $attribs{gene_status}){
            ($hugo_symbols{$attribs{gene_name}}->{gene_status} eq 'KNOWN') or $hugo_symbols{$attribs{gene_name}} = \%attribs;
        }elsif ($hugo_symbols{$attribs{gene_name}}->{level} != $attribs{level}){
            ($hugo_symbols{$attribs{gene_name}}->{level} > $attribs{level}) or $hugo_symbols{$attribs{gene_name}} = \%attribs;
        }else{
            my $length1 = $hugo_symbols{$attribs{gene_name}}->{end} - $hugo_symbols{$attribs{gene_name}}->{start};
            my $length2 = $attribs{end} - $attribs{start};
            $hugo_symbols{$attribs{gene_name}} = \%attribs if( $length1<$length2 );
        }
    }elsif ($type eq "transcript"){
        $transcript_fh->print( $attribs{transcript_id}."\t".$attribs{gene_name}."\t".$attribs{gene_type}."\n" );
    }
}

$gene_fh->close();
$transcript_fh->close();

my %hugo_2_entrez = map{chomp; split(/\t/)}`cat hugo-entrez.txt`;
my $hugo_fh = IO::File->new( "gencode.v19_ucsc.known.txt", ">" ) or die "ERROR: Couldn't overwrite file gencode.v19_ucsc.known.txt\n";
foreach(@gene_names){
    my $entrez = (defined $hugo_2_entrez{$hugo_symbols{$_}->{gene_name}}) ? $hugo_2_entrez{$hugo_symbols{$_}->{gene_name}} : 0;
    $hugo_fh->print ($hugo_symbols{$_}->{gene_id}."\t".$hugo_symbols{$_}->{gene_name}."\t".$entrez."\n");
}
$hugo_fh->close();
