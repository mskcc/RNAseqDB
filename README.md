# RNAseqDB

Background
----------
GTEx is the largest database of RNA-seq normal tissues. TCGA provides a large pool of tumor RNA-seq in addition to RNA-seq of matched normals. The level 3 data of the two studies are not directly comparable. To provide a reference for studies of abnormal gene expressions in diseases, here we want to build a database of gene expression in healthy human tissues.

Study design
----------
We collect RNA-seq of healthy human tissues from GTEx study, and both tumors and normals from TCGA. We uniformly reanalyze the samples: calculate expression and remove study-specific bias using the same pipeline.

Pipeline
----------
We used STAR to align RNA-seq reads, RSEM to quantify gene expression, mRIN to discard degraded samples, and SVAseq to correct batch bias.  Our pipeline also run FeatureCounts to quantify read counts. 

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
