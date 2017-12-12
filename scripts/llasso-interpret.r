## Note that at the end, we will identify the original SNPs by
## X0[common[kept[i]]], the ith column in X is actually the kept[i]
## column in X prior to elimination of repeated columns, and
## the kept[i] column is then the common[kept[i]] in the original
## bedfile before eliminating rare variants
## Later use http://genome.ucsc.edu to identify which genes have hits

setwd("/Users/Clauberry/Dropbox/Documents/gwas/projects/22q_new/llasso/results/test2")
indkept = read.table("22q-chr22.indices-kept", header=TRUE)
str(indkept)
indcommon = read.table("22q-chr22.common", header=TRUE)
str(indcommon)
bim = read.table("/Users/Clauberry/Documents/gwas/data/22q/22q_files_NEW/bedfiles/22q-chr22.bim", header=FALSE)
str(bim) ## V2 = chr22_position, V4= position
names(bim) = c("chr", "chr_position", "x", "position", "a1", "a2")
beta = read.table("22q-chr22.beta", header=TRUE, sep=",")
str(beta)

## significant with Lasso.jl
sigL = which(beta$betaL!=0)
dfL = data.frame(position=bim$position[indcommon$ind[indkept$ind[sigL]]], beta=beta$betaL[sigL])
head(dfL)
library(ggplot2)
ggplot(dfL,aes(x=position, y=beta))+geom_point()

## Now we want to do the same plot but will all positions that were significant with any of the methods
## maybe truncate IHT at 50 because of phase transition plot
sig = beta[,1] != 0
for(i in 2:52){ ## L, SR, 50 IHT
    sig = sig | beta[,i] != 0
}
length(which(sig))
df = beta[sig,1:52]
str(df)
df$position = bim$position[indcommon$ind[indkept$ind[sig]]]
ggplot(df,aes(x=position[position<20000000], y=betaL[position<20000000]))+
    geom_point()+
    geom_point(aes(x=position[position<20000000],y=betaSR[position<20000000]), col="red")



## Later when we do the interpretation:
## We might need the R package for Manhattan plot: `qqman`, and the code that I did for Jai. The function needs BP to be in the order position in each chromosome. Info [here](http://www.gettinggeneticsdone.com/2014/05/qqman-r-package-for-qq-and-manhattan-plots-for-gwas-results.html)
## MANHATTAN PLOT
dat = read.table("results/OptimalSortSkatWChrPos.txt",sep="\t")
dat2 = subset(dat,select=c("V2","V1","V4"))
names(dat2) = c("SNP","CHR","P")
str(dat2)

dat.ord = data.frame()
for(i in 1:max(dat2$CHR)){
  dat.chr = subset(dat2,CHR==i)
  dat.chr = within(dat.chr,BP<-order(dat.chr$SNP))
  str(dat.chr)
  dat.chr = dat.chr[order(dat.chr$BP),]
  dat.ord = rbind(dat.ord, dat.chr)
}
