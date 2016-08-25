#!/usr/bin/Rscript

args <- commandArgs(TRUE)
argsfile <- args
if(argsfile[1]=='help'){
  cat('\n')
  cat('Parameters/Argument needed \n')
  cat('1). Number of markers \n')
  cat('2). Number of animals \n')
  cat('3). % variance eigen value captures \n')

  cat('\n')
  cat('Example run as ... \n')
  cat('./dReduction_pvaluesArgs.R 1000 100 0.999 \n')
} else {

##########################################################
#                 The code below                         #
##########################################################

nsnp=as.numeric(argsfile[1])   # number of SNPs
nanim=as.numeric(argsfile[2])   # number of animals

#### generate genotype data
geno <- matrix(rbinom(n=nsnp*nanim,size=2,prob=c(0.3,0.3,0.4)),nrow=nanim,ncol=nsnp)

##compute the correlation between markers (similarity index)
ZZsnp <- cor(geno)

## compute eigen value decomposition
eigenDsnp <- eigen(ZZsnp)

# store the eigen values 
# Note that the number of eihen values equals to number of markers
eigenval <- data.frame(eigenval=eigenDsnp$values)

#eigenvec <- eigenDsnp$vectors
# calculate the proportion of varaince captured by each eigen value
eigenval$inc <- eigenval$eigenval/sum(eigenval$eigenval)
## calculate the cummulative sum of the proportion
eigenval$incSum <- cumsum(eigenval$inc)

#head(eigenval)
#determine how many markers (eigen values) explain e.g. 99.999% of the variation
howmanySNPs <- length(which(eigenval$incSum<=as.numeric(argsfile[3])))

origPthresh <- round(-log10(0.05/ncol(geno)),3)     ## what is the original P-value threshold
origDreducthresh <- round(-log10(0.05/howmanySNPs),3)    ## what is the new P-value threshold

cat('\n')
cat('... Original P-value Threshold based on \n 1). norminal P-value of 0.05 \n 2). Number of markers ...\n')
cat(paste('-loq10 P-value threshold of ',origPthresh),'\n')

cat('\n')
cat(paste('... New P-value Threshold based on \n 1). norminal P-value of 0.05 \n 2). Number of eigen values explaining',args[3],'% of the variation'),'...\n')
cat(paste('Number of eigen values (implicitely number of independent markers) ==',howmanySNPs),'\n')
cat(paste('-loq10 P-value threshold of ',origDreducthresh),'\n')
}
