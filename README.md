# RNAseqDB

Background
----------
GTEx is the largest database of RNA-seq normal tissues. TCGA provides a large pool of tumor RNA-seq in addition to RNA-seq of matched normals. The level 3 data of the two studies are not directly comparable. To provide a reference for studies of abnormal gene expressions in diseases, here we build a database of gene expression in healthy human tissues.

Study design
----------
We collect RNA-seq of healthy human tissues from GTEx study, and both tumors and normals from TCGA. We uniformly reanalyze the samples: calculate expression and remove study-specific bias using the same pipeline.

Related software
----------
To run our pipeline, user needs to install the following third-party tools (versions used by us were provided) and related modules that they depend on:

 STAR aligner v2.4.2a 
 
 rsem v1.2.20 (http://deweylab.biostat.wisc.edu/rsem/src/rsem-1.2.20.tar.gz)
 
 SubRead v1.5.0 (http://sourceforge.net/projects/subread/files/subread-1.5.0-p1/subread-1.5.0-p1-source.tar.gz/download)

 SVAseq, (R package downloaded from http://bioconductor.org/biocLite.R)
 
 samtools (v1.2)

 RSeQC v1.1.8 (http://www.broadinstitute.org/cancer/cga/tools/rnaseqc/RNA-SeQC_v1.1.8.jar)

 ubu v1.2 (https://github.com/mozack/ubu/archive/v1.2.tar.gz)

 mRIN

 picard-tools v1.126

 FastQC v0.11.3

 bedtools v2.23

Finally, to download FASTQ files from CGHub and dbGaP, user also needs tools [gtdownload]() and [sratoolkit]().

Summary of pipeline
----------
We used STAR to align RNA-seq reads, RSEM to quantify gene expression, mRIN to discard degraded samples, RSeQC
to control quality, and SVAseq to correct batch bias.  Our pipeline also run FeatureCounts to quantify read counts. 

Although this pipeline is not as fast as the latest RNA-seq tools such as Kallisto, it provides the functionality we need:
 1. It generates BAM files that are required by QC tools.
 2. We can use it to exclude low quality samples.
 3. It provides satisfactory accuracy and gene expression comparable accross studies.


Contact
----------

    Qingguo Wang
    Nikolaus Schultz Lab
    Memorial Sloan Kettering Cancer Center
    New York, NY 10065
