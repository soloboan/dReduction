#!/usr/bin/Rscript

nsnp=1000   # number of SNPs
nanim=100   # number of animals

#### generate genotype data
geno <- matrix(rbinom(n=nsnp*nanim,size=2,prob=c(0.3,0.3,0.4)),
               nrow=nanim,ncol=nsnp)

#compute the correlation between markers (similarity index)
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

head(eigenval)
#determine how many markers (eigen values) explain 99.999% of the variation
howmanySNPs <- length(which(eigenval$incSum<=0.99999))

-log10(0.05/ncol(geno))     ## what is the original P-value threshold
-log10(0.05/howmanySNPs)    ## what is the new P-value threshold


