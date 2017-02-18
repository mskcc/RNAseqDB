# RNAseqDB

Background
----------
A multitude of large-scale studies, e.g. TCGA and GTEx, have recently generated an unprecedented volume of RNA-seq data. The RNA-seq expression data from different studies typically are not directly comparable, due to differences in sample and data processing and other batch effects. Here, we developped a pipeline that processes and unifies RNA-seq data from different studies. Using the pipeline, we have processed data from the GTEx and TCGA and have successfully corrected for study-specific biases, allowing comparative analysis across studies. 

Methods
----------
The input of the pipeline is paired-end raw sequencing reads (in FASTQ format). The raw reads of the RNA-seq samples for the TCGA and GTEx projects were retrieved from the Cancer Genomics Hub (CGHub, https://cghub.ucsc.edu) and the Database of Genotypes and Phenotypes (dbGaP, http://www.ncbi.nlm.nih.gov/gap), respectively.

We used STAR to align sequencing reads, RSEM and FeatureCounts to quantify gene expression, mRIN to evaluate sample degradation, RSeQC to measure sample strandness and quality, and SVAseq to correct batch biases.  

Related software
----------
The pipeline requires the following third-party tools, whose full directories needs to be specified in a configuration file config.txt. Several example configuration files are provided in folder [configuration](https://github.com/mskcc/RNAseqDB/tree/master/configuration).

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

These software have already been installed in two computer clusters at MSKCC: HAL and LUNA. So users of these two clusters don't need install them (except Rsubread) to run the pipeline.

Finally, the pipeline wraps [gtdownload]() and [sratoolkit]() to allow user to conveniently download GTEx and TCGA samples from CGHub and dbGaP. To download the data through the pipeline, user needs install [gtdownload]() and [sratoolkit](). 

Availability
----------
The pipeline was designed to run on various computer clusters, e.g. the HAL cluster [hal.cbio.mskcc.org]() and LUNA cluster [luna.cbio.mskcc.org]().  

For LUNA user, the complete pipeline was under directory /ifs/e63data/schultzlab/wangq/bin/RNAseqDB. 

To run the pipeline in LUNA, user needs to specify some paths using the environment variables PATH, PYTHONPATH, and PERL5LIB. The following is a summary of the paths user can add into personal '.bashrc' file. These paths were also summarized in file [configuration/sourceme](https://github.com/mskcc/RNAseqDB/blob/master/configuration/sourceme). 
 1. export PERL5LIB=/ifs/e63data/schultzlab/wangq/perl5:/opt/common/CentOS_6/perl/perl-5.22.0/lib/5.22.0:/ifs/e63data/schultzlab/opt/perl5/lib/perl5:/ifs/e63data/schultzlab/opt/perl5/lib/perl5/czplib
 2. export PATH=/opt/common/CentOS_6-dev/perl/perl-5.22.0/bin:$PATH
 3. export PATH=/opt/common/CentOS_6/python/python-2.7.8/bin/:$PATH
 4. export PYTHONPATH=/ifs/e63data/schultzlab/bin/RSeQC-2.6.1/opt/common/CentOS_6/python/python-2.7.8/lib/python2.7/site-packages:$PYTHONPATH
 5. export PATH=/ifs/e63data/schultzlab/bin/RSeQC-2.6.1/opt/common/CentOS_6/python/python-2.7.8/bin:$PATH
 6. export PATH=/ifs/e63data/schultzlab/wangq/bin/GeneTorrent-download-3.8.7-207/bin:$PATH

In the HAL cluster [hal.cbio.mskcc.org](), a copy of the pipeline was under this directory: /cbio/ski/schultz/home/wangq/scripts

For users outside MSKCC, the source code of the pipeline is freely accessible through GitHub (at https://github.com/mskcc/RNAseqDB/).

Quick start
----------
To demonstrate how to run the pipeline, suppose we have a set of RNA-seq samples under a directory ~/data/RNA-seq/ and we want to quantify expression levels of genes in them. 

If samples are in SRA format, you should put them directly in ~/data/RNA-seq/. For FASTQ or BAM files, you should save each sample in a sub-directory under ~/data/RNA-seq/. The script will automatically find and analyze all the samples under the given directory ~/data/RNA-seq/. For samples downloaded from GTEx or TCGA, user should also put meta data file, SraRunTable.txt for GTEx and manifest.xml and summary.tsv for TCGA, under ~/data/RNA-seq/.

To run the pipeline, firstly initialize the environment variables using the following command:

    source /ifs/e63data/schultzlab/wangq/bin/RNAseqDB/sourceme

Then, use the command below to analyze all the samples under the directory ~/data/RNA-seq/:

    perl /ifs/e63data/schultzlab/wangq/bin/RNAseqDB/pipeline.pl -i ~/data/RNA-seq/ -s

The argument '-s' means to submit a job for each sample. 

The script pipeline.pl requires a configuration file to find and execute other software it needs. If we do not specify it in the command line (like above), the script will look for a default one, config.txt, under its own directory /ifs/e63data/schultzlab/wangq/bin/RNAseqDB/.

Another script pipeline-wrapper.pl provides more functionality, e.g. batch bias correction, than pipeline.pl. When executed with no argument or with the argument '-h', pipeline-wrapper.pl, as well as other script files, will print detailed instructions on how to use it.

After all the jobs terminate, you can (optionally) use another script file collect-qc.pl to create a report on the quality of the samples and to filter out low quality one. Downstream scripts read QC report created by collect-qc.pl and filter out low-quality samples by default. 

    perl /ifs/e63data/schultzlab/wangq/bin/RNAseqDB/collect-qc.pl -i ~/data/RNA-seq/

For expression of the genes in all samples, if you want to create a sample-gene matrix, you can run create-matrix.pl:

    perl /ifs/e63data/schultzlab/wangq/bin/RNAseqDB/create-matrix.pl -i ~/data/RNA-seq/ -o data-matrix-file.txt -p

The script create-matrix.pl provides both raw or normalized outputs: read count, TPM, and FPKM, for both gene or transcript expression. 

Batch bias correction
----------

Our pipeline can correct biases specific to the TCGA and GTEx projects so that TCGA samples are directly comparable with the GTEx samples. 

To do it, you need the configuration file, e.g. [config-luna.txt](https://github.com/mskcc/RNAseqDB/blob/master/configuration/config-luna.txt), to specify paths of the data. In [config-luna.txt](https://github.com/mskcc/RNAseqDB/blob/master/configuration/config-luna.txt), we utilize two variables, 'gtex_path' and 'tcga_path', to point to two directories for storing GTEx and TCGA data, respectively.

The script file, post-process.pl and run-combat.R, do the actual work of batch effect correction. Another script, pipeline-wrapper.pl, wraps all necessary steps together to make the analysis easy. The following is a command I used to analyze bladder tissue:

    perl /ifs/e63data/schultzlab/wangq/bin/RNAseqDB/pipeline-wrapper.pl -t bladder -s

This is what happen when running the command above: the script pipeline-wrapper.pl firstly reads a configuration file [tissue-conf.txt](https://github.com/mskcc/RNAseqDB/blob/master/configuration/tissue-conf.txt) for lines containing the word 'bladder'. Directories storing GTEx and TCGA samples are then obtained by concatenating 'gtex_path' and 'tcga_path' with the keywords found in [tissue-conf.txt](https://github.com/mskcc/RNAseqDB/blob/master/configuration/tissue-conf.txt). Then, the script submits jobs for the samples and does batch bias correction after all job terminates. Finally, the script creates sample-gene matrices from the samples.

Handling replicates
----------

If a sample has more than 2 FASTQ files or has multiple replicates, user should provide a file named SampleSheet.csv together with the FASTQ files. The file SampleSheet.csv should be CSV format and include at least two columns: 'SampleID' and 'Lane'. 

Below is the content of an example SampleSheet.csv file. It is the simplest sample sheet file allowed, as it only contains two required columns. Given this file, the program will look for all FASTQ files with file names matching prefix '130723_7001407_0116_AC2AEVACXX' and lane number <= 8 (each lane is treated as a replicate).

    SampleID,Lane
    130723_7001407_0116_AC2AEVACXX,8

[pipeline.pl](https://github.com/mskcc/RNAseqDB/blob/master/pipeline.pl) provides an arugment '-m | --merge-replicates' to allow user either to merge all replicates of a sample (if they are technical replicates) or analyze each replicate separately (if they are biological replicates).

The following is an example command to merge all replicates of each sample under the directory ~/data/RNA-seq/, 

    perl /ifs/e63data/schultzlab/wangq/bin/RNAseqDB/pipeline.pl -i ~/data/RNA-seq/ -s -m

Handling species other than human
----------
The pipeline can be applied to species other than human. An example configuration file, [config-luna-mouse.txt](https://github.com/mskcc/RNAseqDB/blob/master/configuration/config-luna-mouse.txt), for mouse is provided under folder [configuration](https://github.com/mskcc/RNAseqDB/tree/master/configuration). To quantify gene/transcript expression for mouse samples, user need to specify mouse configuration file (not necessarily in command line):

    perl /ifs/e63data/schultzlab/wangq/bin/RNAseqDB/pipeline.pl -i ~/data/RNA-seq/ -s -m -c /ifs/e63data/schultzlab/wangq/bin/RNAseqDB/config-luna-mouse.txt

Contact
----------

    Qingguo Wang, Jianjiong Gao
    Nikolaus Schultz Lab
    Memorial Sloan Kettering Cancer Center
    New York, NY 10065
