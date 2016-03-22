library(sva)
library(genefilter)
library(limma)

Iteration = 15


# Find cutoff that make ComBat work
GetFilterCutoff <- function(mydata, batch, modcombat, turnedOffGene){
    steps = (1:20)*5;
    for (i in steps){
        print(i)
        filter = apply(mydata, 1, function(x) length(x[x>5])>=i)
        filter = filter & turnedOffGene
        filtered <- as.matrix( log2(mydata[filter,]+1) )
        result <- tryCatch({
            combat_edata1 = ComBat(dat=filtered, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
        }, warning = function(war) {
        }, error = function(err) {
        }, finally = {
            if(exists("combat_edata1")){
                print("ComBat Sucess")
            }else{
                print("ComBat failure")
            }
        })
    
        if(exists("combat_edata1")){
            for (j in (i-4):(i-1)){
                print(j)
                filter = apply(mydata, 1, function(x) length(x[x>5])>=j)
                filter = filter & turnedOffGene
                filtered <- as.matrix( log2(mydata[filter,]+1) )
                
                result2 <- tryCatch({
                    combat_edata2 = ComBat(dat=filtered, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
                }, warning = function(war) {
                }, error = function(err) {
                }, finally = {
                    if(exists("combat_edata2")){
                        print("ComBat Sucess")
                    }else{
                        print("ComBat failure")
                    }
                })
                
                if(exists("combat_edata2")){
                    return(j)
                }
            }
            return(i)
        }
    }
    return(0)
}

# Find genes that fails ComBat
GetTroubleGene <- function (mydata, batch, modcombat, filterCutoff, turnedOffGene){
    filterPass = apply(mydata, 1, function(x) length(x[x>5])>=filterCutoff)
    filterPass = filterPass & turnedOffGene
    filterFail = apply(mydata, 1, function(x) length(x[x>5])>=(filterCutoff-1))
    filterFail = filterFail & turnedOffGene
    
    diff = xor(filterPass,filterFail)
    dub_genes = which(diff)
    
    good_genes = c()
    for (i in dub_genes){
        print(i)
        filter = filterPass
        filter[i] = TRUE
        filtered <- as.matrix( log2(mydata[filter,]+1) )
        
        result2 <- tryCatch({
            ComBat(dat=filtered, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
            good_genes = c(good_genes, i)
        }, warning = function(war) {
        }, error = function(err) {
        }, finally = {
        })
    }
    bad_genes = setdiff(dub_genes, good_genes)
    return(bad_genes)
    
}

# Filter out genes that collapse ComBat
ComBatWrapper <- function (mydata, batch, modcombat, iter_n){
    # Find cutoff that makes ComBat work
    turnedOffGene <- rep(TRUE, length(rownames(mydata)))
    filterCutoff = GetFilterCutoff(mydata, batch, modcombat, turnedOffGene)
    print(paste('Cutoff to filter genes = ', filterCutoff, sep = ''))
    
    for (i in 1:iter_n){
        print(paste('Iteration', i, sep = ''))
        if (filterCutoff <= 2){
            filter = apply(mydata, 1, function(x) length(x[x>5])>=filterCutoff)
            filter[which(turnedOffGene==FALSE)] = FALSE
            filtered <- as.matrix( log2(mydata[filter,]+1) )
            edata = ComBat(dat=filtered, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
            return(edata)
        }else{
            filteredGene = GetTroubleGene(mydata, batch, modcombat, filterCutoff, turnedOffGene)
            turnedOffGene [filteredGene] = FALSE
            
            # Find new cutoff
            newfilterCutoff = GetFilterCutoff(mydata, batch, modcombat, turnedOffGene)
            if(newfilterCutoff != filterCutoff){
                filterCutoff = newfilterCutoff
                print(paste('New cutoff to filter genes = ', filterCutoff, sep = ''))
            }else{
                filter = apply(mydata, 1, function(x) length(x[x>5])>=filterCutoff)
                filter[which(turnedOffGene==FALSE)] = FALSE
                filtered <- as.matrix( log2(mydata[filter,]+1) )
                edata = ComBat(dat=filtered, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
                return(edata)
            }
        }
    }
    filter = apply(mydata, 1, function(x) length(x[x>5])>=filterCutoff)
    filter[which(turnedOffGene==FALSE)] = FALSE
    filtered <- as.matrix( log2(mydata[filter,]+1) )
    edata = ComBat(dat=filtered, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
    return(edata)
}

# Input parameters
args <- commandArgs( trailingOnly = TRUE )
inFile <- args[1]
outPrefix <- args[2]

mydata <- read.table(inFile,row.names=1, header=TRUE, sep="\t", strip.white=TRUE)

# Create model to run ComBat
idx.gtex   = grep("^GTEX",names(mydata))
idx.tcga   = grep("^TCGA",names(mydata))
idx.tumor  = grep(".01[A|B|C].", names(mydata))
idx.tcga_t = intersect(idx.tumor, idx.tcga)

if( file.exists('combat-batch.txt')){
    batch <- read.table('combat-batch.txt',row.names=0, header=F, sep="\t", strip.white=TRUE)
}else{
    batch <- rep(2,length(names(mydata)))
    batch[idx.gtex] = 1
}

if (length (idx.tcga_t) > 0 ){
    cancer = rep('normal',length(names(mydata)))
    cancer[idx.tcga_t] = 'tumor'
    model = data.frame(batch, row.names=names(mydata))
    modcombat = model.matrix(~cancer, data=model)
}else{
    model = data.frame(batch, row.names=names(mydata))
    modcombat = model.matrix(~1, data=model)
}

#filtered <- ComBatWrapper(mydata, batch, modcombat, Iteration)
#filename <- paste(outPrefix, '.png', sep = '');
#png(filename)
#combat_edata = ComBat(dat=filtered, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)

combat_edata = ComBatWrapper(mydata, batch, modcombat, Iteration)

#dev.off()

filename <- paste(outPrefix, '.txt', sep = '');
write.table(combat_edata, file=filename, row.names=TRUE, col.names=NA, quote=F, sep="\t")


