expected_count
--------------
The maximum likelihood gene expression levels estimated using RSEM, i.e. the expected_count in RSEM’s output. It can be provided to programs such as EBSeq, DESeq, or edgeR to identify differentially expressed genes. Here, expression of only protein coding genes was provided.


unnormalized
--------------
Gene expression levels (FPKM) calculated using RSEM. The data in this folder underwent quantile normalization, but was not corrected batch effects.


normalized
--------------
Normalized gene expression levels (FPKM) calculated using RSEM. The data in this folder not only underwent quantile normalization, but also was batch effect corrected (using ComBat).
