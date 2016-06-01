# RNAseqDB

Background
----------
A multitude of large-scale studies, e.g. TCGA and GTEx, have recently generated an unprecedented volume of RNA-seq data. The RNA-seq expression data from different studies typically are not directly comparable, due to differences in sample and data processing and other batch effects. Here, we developped a pipeline that processes and unifies RNA-seq data from different studies. Using the pipeline, we have processed data from the GTEx and TCGA and have successfully corrected for study-specific biases, allowing comparative analysis across studies. 

Methods
----------
The input of the pipeline is raw sequencing reads (in FASTQ format). The raw reads of the RNA-seq samples for the TCGA and GTEx projects were retrieved from the Cancer Genomics Hub (CGHub, https://cghub.ucsc.edu) and the Database of Genotypes and Phenotypes (dbGaP, http://www.ncbi.nlm.nih.gov/gap), respectively.

We used STAR to align sequencing reads, RSEM and FeatureCounts to quantify gene expression, mRIN to evaluate sample degradation, RSeQC
to measure sample strandness and quality, and SVAseq to correct batch biases.  

Comparing with other methods, e.g. Kallisto, our pipeline provides the following important functionality:
 1. It creates BAM files that are required by QC tools.
 2. It measures RNA-seq degradation and excludes degraded samples.
 3. It performs batch effect correction to make samples comparable accross studies.

Related software
----------
The pipeline requires the following third-party tools, which were already installed in the HAL and LUNA cluster (see the last section). So if you are a MSKCC user, probably you don't need install anything.

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

Finally, to download FASTQ files from CGHub and dbGaP, user also needs to install tools [gtdownload]() and [sratoolkit]().

Summary of pipeline
----------
We used STAR to align RNA-seq reads, RSEM and FeatureCounts to quantify gene expression, mRIN to evaluate sample degradation, RSeQC
to measure sample strandness and quality, and SVAseq to correct batch biases.  

Comparing with other methods, e.g. the latest Kallisto, our pipeline provides the following functionality that are important for RNA-seq analysis:
 1. It creates BAM files that are required by QC tools.
 2. It measures sample degradation and excludes degraded samples.
 3. It performs batch effect correction to make samples comparable accross studies.

Installation
----------
The pipeline was installed in the LUNA cluster (luna.cbio.mskcc.org) under the following directory:

/ifs/e63data/schultzlab/wangq/bin/RNAseqDB

The following paths are needed to run [RSeQC]() and [gtdownload]() in LUNA. So it is better to add to your file ~/.bashrc. 
 1. export PATH=/opt/common/CentOS_6/python/python-2.7.8/bin/:$PATH
 2. export PYTHONPATH=/ifs/e63data/schultzlab/bin/RSeQC-2.6.1/opt/common/CentOS_6/python/python-2.7.8/lib/python2.7/site-packages:$PYTHONPATH
 3. export PATH=/ifs/e63data/schultzlab/bin/RSeQC-2.6.1/opt/common/CentOS_6/python/python-2.7.8/bin:$PATH
 4. export PATH=/ifs/e63data/schultzlab/wangq/bin/GeneTorrent-download-3.8.7-207/bin:$PATH


In the HAL cluster (hal.cbio.mskcc.org), it was installed under this directory: 

/cbio/ski/schultz/home/wangq/scripts

Contact
----------

    Qingguo Wang
    Nikolaus Schultz Lab
    Memorial Sloan Kettering Cancer Center
    New York, NY 10065
