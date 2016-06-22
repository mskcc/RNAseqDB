#!/opt/common/CentOS_6/bin/current/Rscript
# Perform differential expression analysis

library(limma)
library(edgeR)

# read in target file
options(digits=2)
options(max.print=10000000)


ReadDataTable <- function(filename){
    tempdata = read.table(filename, row.names=1, header=T, sep="\t", strip.white=T)
    s_names = names(tempdata)
    n = length(s_names)
    x = rownames(tempdata)
    if(s_names[1] == 'Entrez_Gene_Id' | s_names[1] == 'Description')
    {
        tempdata = tempdata[order(x), 2:n]
    }else{
        tempdata = tempdata[order(x), ]
    }
    return(tempdata)
}

# Input parameters
args <- commandArgs( trailingOnly = TRUE )
# inFile1 should be TCGA normal samples
# inFile2 should be TCGA tumor samples
inFile1 <- args[1]
inFile2 <- args[2]

# Load two data tables
mydata1 <- ReadDataTable(inFile1)
#dim(mydata1)
mydata2 <- ReadDataTable(inFile2)
#dim(mydata2)

# Filter out genes appearing in only one file
list_of_data = list(mydata1, mydata2)
common_names = Reduce(intersect, lapply(list_of_data, row.names))
list_of_data = lapply(list_of_data, function(x) { x[row.names(x) %in% common_names,] })

mydata1 = list_of_data[[1]]
mydata2 = list_of_data[[2]]

mydata  = cbind(mydata1, mydata2)
#dim(mydata)

# Create variable 'design'
n1 = length(names(mydata1))
n2 = length(names(mydata2))

design <- cbind( Grp1=1, Grp2vs1=c(rep(1,n1), rep(2,n2)) )

v_max = max(mydata)
v_min = min(mydata)

if(v_min >= 0 & v_max >= 1000){  # If the expression data was not normalized before
    # filter out low-count genes
    isexpr <- rowSums(cpm(mydata) > 5) >= 10

    DGE=DGEList(mydata[isexpr,])
    DGE=calcNormFactors(DGE,method =c("TMM"))

    # perform voom normalization
    mydata <- voom(DGE,design,plot=TRUE)
}

# fit linear model and assess differential expression
fit <- eBayes(lmFit(mydata,design))

# Print results
# topTable(fit,coef=2)
topTable(fit,coef=2,sort="p",n=Inf)

