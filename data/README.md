expected_count
--------------
The maximum likelihood gene expression levels estimated using RSEM, i.e. the expected_count in RSEM’s output. There are 52 data files in this subdirectory, each being a sample-gene matrix of a certain tissue type. These files can be provided to programs such as EBSeq, DESeq, or edgeR for identifying differentially expressed genes. The expression of only protein coding genes was provided here.


unnormalized
--------------
The gene expression levels calculated from fpkm of RSEM’s output. The data matrices in this subdirectory were not the direct output of RSEM. They underwent quantile normalization, but were not corrected for batch effects.


normalized
--------------
The normalized gene expression levels (FPKM) calculated using RSEM. This set of data files not only was quantile normalized, but was corrected for batch effects (using tool ComBat).


file naming
--------------
For each tissue type, TCGA tumor and matched normal samples are stored in two separate files: tumor in file `xxxx-xxxx-xxxx-`**tcga-t**`.txt.gz` while normal in `xxxx-xxxx-xxxx-`**tcga**`.txt.gz`. For example, TCGA bladder urothelial carcinoma are in  file `blca-rsem-fpkm-`**tcga-t**`.txt.gz` and the corresponding bladder normal tissue in file `blca-rsem-fpkm-`**tcga**`.txt.gz`. 


GTEx/TCGA tissue mapping
--------------
The mapping of GTEx tissues to TCGA tissues is provided in [this table](https://www.nature.com/articles/sdata201861/tables/1).
