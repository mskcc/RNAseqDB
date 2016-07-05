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

The pipeline also needs [gtdownload]() and [sratoolkit](). The pipeline wraps [gtdownload]() and [sratoolkit]() to allow user to conveniently download GTEx and TCGA samples from CGHub and dbGaP.

These software have already been installed in the HAL and LUNA cluster (see section below). If you would like to run the pipeline in LUNA, you don't need install any of these software except Rsubread.

Availability
----------
The pipeline is designed to run on various computer clusters, e.g. the HAL cluster [hal.cbio.mskcc.org]() and LUNA cluster [luna.cbio.mskcc.org]().  

For LUNA user, the complete pipeline was installed under directory /ifs/e63data/schultzlab/wangq/bin/RNAseqDB. 

To run the pipeline in LUNA, user needs to specify some paths using the environment variables PATH, PYTHONPATH, and PERL5LIB. The following is a summary of the paths user can add into personal '.bashrc' file. These paths were also summarized in file [configuration/sourceme](https://github.com/mskcc/RNAseqDB/blob/master/configuration/sourceme). 
 1. export PERL5LIB=/ifs/e63data/schultzlab/wangq/perl5:/opt/common/CentOS_6/perl/perl-5.22.0/lib/5.22.0:/ifs/e63data/schultzlab/opt/perl5/lib/perl5:/ifs/e63data/schultzlab/opt/perl5/lib/perl5/czplib
 2. export PATH=/opt/common/CentOS_6-dev/perl/perl-5.22.0/bin:$PATH
 3. export PATH=/opt/common/CentOS_6/python/python-2.7.8/bin/:$PATH
 4. export PYTHONPATH=/ifs/e63data/schultzlab/bin/RSeQC-2.6.1/opt/common/CentOS_6/python/python-2.7.8/lib/python2.7/site-packages:$PYTHONPATH
 5. export PATH=/ifs/e63data/schultzlab/bin/RSeQC-2.6.1/opt/common/CentOS_6/python/python-2.7.8/bin:$PATH
 6. export PATH=/ifs/e63data/schultzlab/wangq/bin/GeneTorrent-download-3.8.7-207/bin:$PATH

In the HAL cluster [hal.cbio.mskcc.org](), the program was installed under this directory: 

/cbio/ski/schultz/home/wangq/scripts

Quick start
----------
Suppose you want to analyze a set of RNA-seq samples under a directory ~/data/RNA-seq/. 

If samples are in SRA format, you should put them directly in ~/data/RNA-seq/. For FASTQ or BAM files, you should save each sample in one sub-directory under ~/data/RNA-seq/. The script will automatically find and analyze all the samples under the given directory ~/data/RNA-seq/. For samples downloaded from GTEx or TCGA, user should also put meta data file, SraRunTable.txt for GTEx and manifest.xml for TCGA, under ~/data/RNA-seq/.

To run the pipeline, firstly initialize the environment variables using the following command:

    source /ifs/e63data/schultzlab/wangq/bin/RNAseqDB/sourceme

Then, use the command below to analyze all the samples under the directory ~/data/RNA-seq/:

    perl /ifs/e63data/schultzlab/wangq/bin/RNAseqDB/pipeline.pl -i ~/data/RNA-seq/ -s

The argument '-s' means to submit a job for each sample. 

The script pipeline.pl requires a configuration file to find and execute other software it needs. But we did not specified it in the command line above. In this case, The script utilizes a default configuration file, config.txt, under its own directory /ifs/e63data/schultzlab/wangq/bin/RNAseqDB/.

Another script pipeline-wrapper.pl provides more functionality, e.g. batch bias correction, than pipeline.pl. When executed with no argument or with the argument “-h”, pipeline-wrapper.pl, as well as other script file, will print detailed instructions on how to use it.

After all the jobs terminate, you can (optionally) use another script file collect-qc.pl to create a report on the quality of the samples and to filter out low quality one. 

    perl /ifs/e63data/schultzlab/wangq/bin/RNAseqDB/collect-qc.pl -i ~/data/RNA-seq/

For expression of the genes in all samples, if you want to create a sample-gene matrix, you can run create-matrix.pl:

    perl /ifs/e63data/schultzlab/wangq/bin/RNAseqDB/create-matrix.pl -i ~/data/RNA-seq/ -o data-matrix-file.txt -p

The script create-matrix.pl provides both raw or normalized outputs: read count, TPM, and FPKM, for both gene or transcript expression. 

Batch bias correction
----------

If you want to compare TCGA samples with GTEx samples, our pipeline can help you correct biases specific to the TCGA and GTEx projects. 

To do it, you need specify paths of your data in the configuration file. In the example configuration file [config-luna.txt](https://github.com/mskcc/RNAseqDB/blob/master/configuration/config-luna.txt), variables 'gtex_path' and 'tcga_path' point to two directories that are used to store GTEx and TCGA data, respectively.

The script file, post-process.pl and run-combat.R, do the actual work of batch effect correction. Another script, pipeline-wrapper.pl, wraps all necessary steps together to make analysis easy. The following is the command I used to analyze bladder tissue:

    perl /ifs/e63data/schultzlab/wangq/bin/RNAseqDB/pipeline-wrapper.pl -t bladder -s

This is what happen when running the command above: the script pipeline-wrapper.pl firstly reads a configuration file [tissue-conf.txt](https://github.com/mskcc/RNAseqDB/blob/master/configuration/tissue-conf.txt) for lines containing the word 'bladder'. Directories storing GTEx and TCGA samples are then obtained by concatenating 'gtex_path' and 'tcga_path' with the keywords found in [tissue-conf.txt](https://github.com/mskcc/RNAseqDB/blob/master/configuration/tissue-conf.txt). Then, the script submits jobs for the samples and does batch bias correction after all jobs terminates. Finally, sample-gene matrices are created.

Handling replicates
----------

If a sample has more than 2 FASTQ files or has multiple replicates, user should provide a file named SampleSheet.csv under the same directory of the sample FASTQ files. The file SampleSheet.csv should be CSV format and include at least two columns: 'SampleID' and 'Lane'. 

Below is the content of an example SampleSheet.csv file. It is the simplest sample sheet file allowed, as it contains only two necessary columns. With it, the program will look for all FASTQ files matching name prefix '130723_7001407_0116_AC2AEVACXX' and with lane number<=8 (each lane is regarded a replicate).      
    SampleID,Lane
    130723_7001407_0116_AC2AEVACXX,8

In the script file [pipeline.pl](https://github.com/mskcc/RNAseqDB/blob/master/pipeline.pl), an arugment '-m | --merge-replicates' is provided to allow user either to merge all replicates of each sample in the analysis (if '-m' is specified) or analyze each replicate separately.

The following is an example command to merge all replicates of each sample under directory ~/data/RNA-seq/, 

    perl /ifs/e63data/schultzlab/wangq/bin/RNAseqDB/pipeline.pl -i ~/data/RNA-seq/ -s -m


Contact
----------

    Qingguo Wang
    Nikolaus Schultz Lab
    Memorial Sloan Kettering Cancer Center
    New York, NY 10065
