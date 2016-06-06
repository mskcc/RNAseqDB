#!/usr/bin/perl -w

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use Cwd;
use FindBin;
use lib "$FindBin::Bin";
use Switch;

my $rsem_rnd_seed = 12345;

my @usage;
push @usage, "\nUsage:  calc-expression.pl [-i input]\n\n";
push @usage, "Options:\n";
push @usage, "  -h, --help     Displays this information.\n";
push @usage, "  -c | --config  Configuration file\n";
push @usage, "  -i | --input   An input sra/fastq/bam file or directory. For sra/fastq/bam,\n";
push @usage, "                 FASTQ will be extracted and processed(realignment and RPKM\n";
push @usage, "                 quantification). For directory, unzip step will be skipped.\n\n";
push @usage, "  --fq1          Forward end of paired-end reads. Will override -i if provided\n";
push @usage, "  --fq2          Reverse end of paired-end reads. Will override -i if provided\n";


my ( $help, $config_file, $input, $fq1, $fq2);

GetOptions
(
 'h|help|?'   => \$help,
 'c|config=s' => \$config_file,
 'i|input=s'  => \$input,
 'fq1=s'      => \$fq1,
 'fq2=s'      => \$fq2,
);

if ($help or (!defined $input and !defined $fq1 and !defined $fq2)) {
   print @usage;
   exit(0);
}

######################### Processes configuration file #################################
(defined $config_file) or $config_file = "$FindBin::Bin/config.txt";
( -e $config_file ) or die "ERROR: The configuration file $config_file does not exist\n";

my ( %config, $thread_n, $rsem_dir,  $index_rsem, $index_star, $picard_dir, $gencode, $house_keeping_genes, $ubu_dir );
my ( $bedtools_dir, $star_bin, $rseqc_dir, $fastqc_bin, $sratool_dir, $bowtie_dir, $kallisto_bin, $kallisto_index, $mRIN_dir );

# Initialize variables
map{ chomp; /^\s*([^=\s]+)\s*=\s*(.*)$/; $config{$1} = $2 if (defined $1 and defined $2) } `egrep -v \"^#\" $config_file`;
$thread_n           = 1;
$thread_n           = $config{ thread_no }           if ( exists $config{ thread_no  } );
$rsem_dir           = $config{ rsem_dir }            if ( exists $config{ rsem_dir   } );
$gencode            = $config{ gencode }             if ( exists $config{ gencode } );
$index_rsem         = $config{ index_rsem }          if ( exists $config{ index_rsem } );
$index_star         = $config{ index_star }          if ( exists $config{ index_star } );
$picard_dir         = $config{ picard_dir }          if ( exists $config{ picard_dir } );
$star_bin           = $config{ star_bin }            if ( exists $config{ star_bin } );
$rseqc_dir          = $config{ rseqc_dir }           if ( exists $config{ rseqc_dir } );
$ubu_dir            = $config{ ubu_dir }             if ( exists $config{ ubu_dir } );
$bedtools_dir       = $config{ bedtools_dir }        if ( exists $config{ bedtools_dir } );
$house_keeping_genes= $config{ house_keeping_genes } if ( exists $config{ house_keeping_genes } );
$fastqc_bin         = $config{ fastqc_bin }          if ( exists $config{ fastqc_bin } );
$sratool_dir        = $config{ sratool_dir }         if ( exists $config{ sratool_dir } );
$bowtie_dir         = $config{ bowtie_dir }          if ( exists $config{ bowtie_dir } );
$kallisto_bin       = $config{ kallisto }            if ( exists $config{ kallisto } );
$kallisto_index     = $config{ kallisto_index }      if ( exists $config{ kallisto_index } );
$mRIN_dir           = $config{ mRIN_dir }            if ( exists $config{ mRIN_dir } );

################################## Process input parameter #################################

if (defined $fq1 and defined $fq2) {
    (-e $fq1 and -e $fq2) or die "ERROR: $fq1 or $fq2 does not exist\n";
    
    # print "Running pipeline STAR+RSEM...\n";
    Run_Star($fq1, $fq2);
    Run_RSEM() if (defined $index_rsem);

    #print "Running pipeline Bowtie2+RSEM...\n";
    #Run_Bowtie_RSEM($file1, $file2) if(!-e "Aligned.toTranscriptome.out.bam");

    #print "Running kallisto...\n";
    #Run_kallisto($fq1, $fq2);
    
}elsif (defined $input) {
    (-e $input) or die "ERROR: $input does not exist\n";
    print "Processing $input...\n";

    my $new_fastq = 0;
    if (-d $input){
        chdir $input;
    }else{
        # Extract FASTQ
        print "Extract FASTQ...\n";
        ExtractFastQ($input);
        $new_fastq = 1;
    }

    my @fqs;
    @fqs = glob( "*.fastq.gz" );
    @fqs = glob( "*.fastq" ) if (scalar @fqs == 0);
    @fqs = glob( "*.fq.gz" ) if (scalar @fqs == 0);
    @fqs = glob( "*.fq" )    if (scalar @fqs == 0);
    
    # Perform alignment
    if ( scalar( @fqs ) > 1 and (!-e "Log.final.out" or !-e "Aligned.sortedByCoord.out.bam")){
        
        # print "Running pipeline STAR+RSEM...\n";
        if (!-e 'SampleSheet.csv'){
            Run_Star($fqs[0], $fqs[1]);
        }else{
            my ($lane, $id);
            ($lane, $id) = ParseSampleSheet("SampleSheet.csv");
            my ($fq1, $fq2);
            map{ $fq1=$_ if(/L0+1/ and /R1/); $fq2=$_ if(/L0+1/ and /R2/)}@fqs;
            for (my $i=2; $i<=$lane; $i++){
                map{ $fq1 .= ",$_" if(/L0+$i/ and /R1/); $fq2 .= ",$_" if(/L0+$i/ and /R2/) }@fqs;
            }
            Run_Star($fq1, $fq2);
        }
        `rm *_?.fastq.gz` if ($new_fastq and -s "Aligned.sortedByCoord.out.bam");
        
        Run_RSEM() if (defined $index_rsem);
        
        #print "Running pipeline Bowtie2+RSEM...\n";
        #Run_Bowtie_RSEM($fqs[0], $fqs[1]);
        
        #print "Running kallisto...\n";
        #Run_kallisto($fqs[0], $fqs[1]);
    }
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

################################ Quantify Expression #######################################

if(-e 'Log.final.out' and -e 'Aligned.out.bam' and !-e 'Aligned.sortedByCoord.out.bam'){
    `samtools sort -o Aligned.sortedByCoord.out.bam -O bam -T Aligned.sortedByCoord.out.bam.0 Aligned.out.bam`;
    `samtools index Aligned.sortedByCoord.out.bam`;
}

(-e 'Log.final.out' and -e 'Aligned.sortedByCoord.out.bam') or die "ERROR: STAR aligner failed to align reads\n";

if (defined $index_rsem){
    (-e 'Quant.genes.results' and -e 'Quant.isoforms.results') or Run_RSEM(); #die "ERROR: Did not find RSEM output\n";
}

(-e 'fcounts.tpm' and -e 'fcounts.fpkm') or Run_FeatureCounts(); #die "ERROR: Did not find FeatureCounts output\n";


if( !-e 'Aligned.sortedByCoord.out.bam.bai' ){
    print "Create index for Aligned.sortedByCoord.out.bam\n";
    `samtools index Aligned.sortedByCoord.out.bam`;
}


(-e 'ks/sample.ks.txt') or Run_mRIN('Aligned.sortedByCoord.out.bam');

# Do not follow TCGA pipeline to normalize rsem output
QuantileNorm('Quant.genes.results', -1, 'rsem.genes.normalized_results') if (defined $index_rsem and !-s 'ubu-quan/rsem.genes.normalized_results');

QuantileNorm('fcounts.fpkm', 2, 'fcounts.fpkm.normalized_results') if (!-s 'ubu-quan/fcounts.fpkm.normalized_results');


##################################### Calculate QC #####################################

# Collect RNA-seq metrics
#(-e 'PicardRNASeqMetrics.txt') or `java -Xmx2g -jar $picard_dir/picard.jar CollectRnaSeqMetrics REF_FLAT=$gencode.genePred INPUT=Aligned.sortedByCoord.out.bam OUTPUT=PicardRNASeqMetrics.txt CHART=PicardRNASeqMetrics.pdf STRAND_SPECIFICITY=NONE VALIDATION_STRINGENCY=LENIENT`;

# Calculate gene coverage
if (defined $house_keeping_genes){
    (-e 'rseqc.geneBodyCoverage.txt') or `geneBody_coverage.py -r $house_keeping_genes -i Aligned.sortedByCoord.out.bam -o rseqc`;
}

# Run FastQC
(-e 'Aligned.sortedByCoord.out_fastqc') or `$fastqc_bin --extract Aligned.sortedByCoord.out.bam`;


############################### Extract FASTQ #######################################

# Extract FASTQ files
sub ExtractFastQ {
    my $inFile = shift;
    if ($inFile =~ /.sra$/){
        my $path = $inFile;
        $path =~ s/.sra$//;
        `$sratool_dir/fastq-dump --split-3 $inFile -O $path/` if (!-e "$path/Aligned.sortedByCoord.out.bam");
        chdir $path;
    }elsif ($inFile =~ /.tar.gz$/){
        my ($name, $path) = fileparse($inFile);
        chdir $path;
        `tar xzf $name` if (!-e "Aligned.sortedByCoord.out.bam");
    }elsif ($inFile =~ /.bam$/){
        my ($name, $path) = fileparse($inFile);
        chdir $path;
        if (!-e "Aligned.sortedByCoord.out.bam"){
            #`java -jar -Xmx40g $picard_dir/picard.jar SamToFastq INPUT=$name FASTQ=read_1.fastq SECOND_END_FASTQ=read_2.fastq INCLUDE_NON_PF_READS=True VALIDATION_STRINGENCY=SILENT`;
            `samtools view -b -o paired.bam -\@ $thread_n -f 1 $name`;
            `samtools sort -n -o namesort.bam -T namesort_pre -\@ $thread_n -m 3G -O bam paired.bam`;
            `java -Xmx512M -jar $ubu_dir/ubu-1.2-jar-with-dependencies.jar sam2fastq --in namesort.bam --fastq1 read_1.fastq --fastq2 read_2.fastq --end1 /1 --end2 /2`;
            #`$bedtools_dir/bamToFastq -i namesort.bam -fq read_1.fastq -fq2 read_2.fastq`;
            `rm paired.bam namesort.bam` if(-s "read_1.fastq" and "read_2.fastq");
        }
    }else{
        die "ERROR: unknow input data format!\n";
    }
    my @fastq_files = glob("*.fastq");
    `gzip *.fastq` if(scalar @fastq_files > 0);
}

############################### Determine Strandness #######################################

# Get Strandness of RNA-seq experiment
sub GetStrandness {
    my $bedFile = shift;
    my $bamFile = shift;
    my @ret = `$rseqc_dir/infer_experiment.py -r $bedFile -i $bamFile`;
    
    my $paired = 0;
    my ($frac1, $frac2) = (0, 0);
    
    foreach( @ret ){
        if( /This is PairEnd Data/ ){
            $paired = 1;
        }elsif( /Fraction of reads explained by/ ){
            if ($frac1 == 0) {
                $frac1 = m/\": (\S+)$/;
            } else {
                $frac2 = m/\": (\S+)$/;
            }
        }
    }
    my $ret = '_SE';
    $ret = '_PE' if( $paired );
    
    if ($frac1/$frac2 > 2 or $frac1/$frac2 < 0.5){
        $ret = 'str' . $ret;
    }else{
        $ret = 'unstr' . $ret;
    }
    return $ret;
}


############################### STAR + RSEM #######################################

sub Run_kallisto {
    my $fastq1 = shift;
    my $fastq2 = shift;
    (-e 'out.bam') or `$kallisto_bin quant -i $kallisto_index -o . --pseudobam --plaintext $fastq1 $fastq2 | samtools view -Sb - -o out.bam`;
}

sub Run_mRIN {
    my $bam = shift;
    (-e $bam) or die "ERROR: $bam does not exist\n";

    #`samtools view -b -F 1548 -q 30 $bam | $bedtools_dir/bamToBed -i stdin > sample.bed` if (!-e 'sample.bed');
    if (!-e 'sample.bedGraph') {
        `samtools view -b -F 1548 -q 30 $bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM | $bedtools_dir/bamToBed -i stdin > sample.bed` if (!-e 'sample.bed');
    
        `perl $mRIN_dir/tag2profile.pl -big -exact -of bedgraph -v sample.bed sample.bedGraph`;
        `rm -f sample.bed` if (-e 'sample.bed' and -e 'sample.bedGraph');
    }
    `mkdir -p ks`;
    #`perl ~/bin/mRIN/gen_transcript_cdf.pl -v ~/bin/mRIN/genes/refGene.rep.uniq.hg38.bed sample.bedGraph cdf/sample.cdf.bedGraph`;
    #`perl ~/bin/mRIN/ks_test_uniform.pl    -v cdf/sample.cum.bedGraph ks/sample.ks.txt}`;
    #`perl  /cbio/ski/schultz/home/wangq/bin/mRIN/gen_transcript_cdf.pl -v /cbio/ski/schultz/home/wangq/bin/mRIN/genes/refGene.rep.uniq.hg38.bed sample.bedGraph - | perl /cbio/ski/schultz/home/wangq/bin/mRIN/ks_test_uniform.pl -v - ks/sample.ks.txt`;
    #`perl  /cbio/ski/schultz/home/wangq/bin/mRIN/gen_transcript_cdf.pl -v /cbio/ski/schultz/home/wangq/ref/transcript/gencode.v21/gencode.v21.annotation.bed sample.bedGraph - | perl /cbio/ski/schultz/home/wangq/bin/mRIN/ks_test_uniform.pl -v - ks/sample.ks.txt`;
    `perl  $mRIN_dir/gen_transcript_cdf.pl -v $mRIN_dir/genes/refGene.rep.uniq.hg19.bed sample.bedGraph - | perl $mRIN_dir/ks_test_uniform.pl -v - ks/sample.ks.txt`;
    `rm -f sample.bedGraph` if (-e 'sample.bedGraph' and -e 'ks/sample.ks.txt');
}

############################### STAR + RSEM #######################################

# Run STAR to align reads
sub Run_Star {
    my $fastq1 = shift;
    my $fastq2 = shift;

    return if (-e "Aligned.sortedByCoord.out.bam");

    # STAR parameters: common (options used here are similar as https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/dnanexus/align-star-pe/resources/usr/bin/lrna_align_star_pe.sh)
    my $star_para_comm = "--genomeDir $index_star --readFilesIn $fastq1 $fastq2  --outFilterType BySJout --outSAMattributes NH HI AS NM MD  --outFilterMultimapNmax 20   --outFilterMismatchNmax 999  --outFilterMismatchNoverLmax 0.04   --alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000 --alignSJoverhangMin 8  --alignSJDBoverhangMin 1";
    
    my $read_fq_cmd = '';
    if ($fastq1 =~ /.bz2$/) {
        $read_fq_cmd = '--readFilesCommand bzip2 -c';
    } elsif ($fastq1 =~ /.gz$/) {
        $read_fq_cmd = '--readFilesCommand zcat';
    } elsif ($fastq1 =~ /.fq|.fastq|.txt$/) {
    } else {
        die ("ERROR: Input file is not fastq format\n");
    }
    
    # STAR parameters: run-time, controlled by DCC
    my $star_run = "--runThreadN $thread_n";
    
    # STAR parameters: type of BAM output: quantification or sorted BAM or both
    #     OPTION: sorted BAM output
    ## STARparBAM="--outSAMtype BAM SortedByCoordinate"
    #     OPTION: transcritomic BAM for quantification
    ## STARparBAM="--outSAMtype None --quantMode TranscriptomeSAM"
    #     OPTION: both
    ## STARparBAM="--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM"
    
    # my $star_bam = "--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM";
    ## STAR parameters: metadata
    # my $star_meta = '--outSAMheaderHD @HD VN:1.4 SO:coordinate';
    # print "$star_bin $star_para_comm $read_fq_cmd $star_run $star_bam $star_meta\n";
    # `$star_bin $star_para_comm $read_fq_cmd $star_run $star_bam $star_meta`;
    
    my $star_bam = "--outSAMtype BAM Unsorted --quantMode TranscriptomeSAM";
    
    # STAR parameters: metadata
    my $star_meta = '--outSAMheaderHD @HD VN:1.4 SO:coordinate';
    
    print "$star_bin $star_para_comm $read_fq_cmd $star_run $star_bam $star_meta\n";
    `$star_bin $star_para_comm $read_fq_cmd $star_run $star_bam $star_meta`;
    
    if(-e "Aligned.out.bam"){
        `samtools sort -o Aligned.sortedByCoord.out.bam -T Aligned.sortedByCoord.out.bam.0 -O bam Aligned.out.bam`;
        `samtools index Aligned.sortedByCoord.out.bam`;
        `rm -f Aligned.out.bam`;
    }
}

sub Run_RSEM {
    
    (-e 'Aligned.toTranscriptome.out.bam' and -e 'Log.final.out') or die "ERROR: STAR aligner failed to align reads\n";
    
    # Run RSEM on BAM file
    my $rsem_para_comm = "--bam --estimate-rspd  --calc-ci --no-bam-output --seed $rsem_rnd_seed";
    
    # RSEM parameters: run-time, number of threads and RAM in MB
    my $rsem_run = "-p $thread_n --ci-memory 30000";
    
    # RSEM parameters: data type dependent
    my $gencode_bed = $gencode;
    $gencode_bed =~ s/gtf/bed/;
    
    my $data_type = GetStrandness( $gencode_bed, 'Aligned.sortedByCoord.out.bam' );
    my $rsem_para_type;
    
    switch ($data_type) {
        case 'str_SE'   { $rsem_para_type = "--forward-prob 0" }
        case 'str_PE'   { $rsem_para_type = "--paired-end --forward-prob 0" }
        case 'unstr_SE' { $rsem_para_type = "" }
        case 'unstr_PE' { $rsem_para_type = "--paired-end" }
        else		    { die "ERROR: unknown data type\n" }
    }
    
    # RSEM command
    print "Calculate expression using RSEM...\n";
    print "\n$rsem_dir/rsem-calculate-expression $rsem_para_comm $rsem_run $rsem_para_type Aligned.toTranscriptome.out.bam $index_rsem Quant >& Log.rsem\n";
    (-e 'Quant.genes.results') or `$rsem_dir/rsem-calculate-expression $rsem_para_comm $rsem_run $rsem_para_type Aligned.toTranscriptome.out.bam $index_rsem Quant >& Log.rsem`;
    
    # Filtering using ubu.jar does not work
    #`java -Xmx20g -jar $ubu_dir/ubu.jar sam-filter --in Aligned.toTranscriptome.out.bam --out Aligned.toTranscriptome.filtered.bam --strip- indels --max-insert 10000 --mapq 1 > sam_filter.log`;
    #print "\n$rsem_dir/rsem-calculate-expression $rsem_para_comm $rsem_run $rsem_para_type Aligned.toTranscriptome.filtered.bam $index_rsem Quant.filtered >& Log.filtered.rsem\n";
    #`$rsem_dir/rsem-calculate-expression $rsem_para_comm $rsem_run $rsem_para_type Aligned.toTranscriptome.filtered.bam $index_rsem Quant.filtered >& Log.filtered.rsem`;
 
}

############################### FeatureCounts #######################################

# Run STAR to align reads
sub Run_FeatureCounts {
    my $bamFile = 'Aligned.sortedByCoord.out.bam';
    $bamFile = 'Aligned.out.bam' if (-e 'Aligned.out.bam');
    
    (-e $bamFile) or die "ERROR: Fail to find alignment file $bamFile\n";

    print "Quantifying expression using FeatureCounts...\n";
    
    `Rscript $FindBin::Bin/run-featurecounts.R $bamFile $gencode $thread_n fcounts`;
}

############################### Bowtie + RSEM #######################################

# Run STAR to align reads
sub Run_Bowtie_RSEM {
    my $fastq1 = shift;
    my $fastq2 = shift;
    `$rsem_dir/rsem-calculate-expression -p $thread_n --bowtie2 --output-genome-bam  --calc-ci --ci-memory 30000 --seed $rsem_rnd_seed --bowtie2-path $bowtie_dir --paired-end $fastq1 $fastq2 $index_rsem Quant >& Log.rsem`;
    `ln -s Quant.transcript.bam Aligned.toTranscriptome.out.bam`;
    `samtools sort -o Aligned.sortedByCoord.out.bam -T Aligned.sortedByCoord.out.bam.0 -O bam Quant.transcript.bam`;
    `rm -rf Quant.temp`;
}

############################### Quantile Normalization #######################################

# Follow TCGA pipeline to normalize rsem output
sub QuantileNorm  {
    my $quant_file = shift;
    my $col_idx    = shift;
    my $out_file   = shift;
    
    print "Perform quantile normalization...\n";
    `mkdir ubu-quan` if (!-e 'ubu-quan');

    # Read file header and find out the FPKM column
    if ($col_idx == -1) {
        my $idx = 0;
        map{ $idx++; $col_idx = $idx if($_ eq "FPKM") }split(/\t/, `head $quant_file | grep ^gene_id`);
        ($col_idx != -1) or die "Unknown file format: $quant_file\n";
    }
    
    #`perl $ubu_dir/perl/strip_trailing_tabs.pl --input Quant.isoforms.results --temp ubu-quan/orig.isoforms.results > ubu-quan/trim_isoform_tabs.log`;
    #`perl $ubu_dir/perl/quartile_norm.pl -c $col_idx -q 75 -t 300  -o ubu-quan/rsem.isoforms.normalized_results ubu-quan/orig.isoforms.results`;
    
    `grep -v gene_id $quant_file > ubu-quan/rsem.genes.results`;
    `perl $ubu_dir/perl/quartile_norm.pl -c $col_idx -q 75 -t 1000 -o ubu-quan/$out_file ubu-quan/rsem.genes.results`;
    `rm -f ubu-quan/rsem.genes.results`;
}
