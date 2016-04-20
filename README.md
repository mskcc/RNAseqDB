# RNAseqDB

Background
----------
GTEx is the largest study of RNA-seq healthy tissues. TCGA is the most comprehensive study of human cancer, providing a great number of normal and tumor RNA-seq samples. The level 3 data of the two studies are not directly comparable. To provide a reference for investigation of abnormal gene expressions in diseases, here we propose a pipeline to quantify RNA-seq gene expression and remove study-specific batch effects.

Study design
----------
We collect RNA-seq of healthy human tissues from GTEx study, both tumors and normals from TCGA. Starting with raw reads downloaded from CGHub and dbGaP, we uniformly reanalyze the samples: calculating expression and removing study-specific biases.

Related software
----------
To run the pipeline, user needs to install the following third-party tools and related modules that they depend on:

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

Finally, to be able to download FASTQ files from CGHub and dbGaP, user also needs tools [gtdownload]() and [sratoolkit]().

Summary of pipeline
----------
We used STAR to align RNA-seq reads, RSEM and FeatureCounts to quantify gene expression, mRIN to evaluate sample degradation, RSeQC
to measure sample strandness and quality, and SVAseq to correct batch biases.  

Comparing with other methods, e.g. the latest Kallisto, our pipeline provides the following functionality that are important for RNA-seq analysis:
 1. It creates BAM files that are required by QC tools.
 2. It measures sample degradation and excludes degraded samples.
 3. It performs batch effect correction to make samples comparable accross studies.


Contact
----------

    Qingguo Wang
    Nikolaus Schultz Lab
    Memorial Sloan Kettering Cancer Center
    New York, NY 10065
