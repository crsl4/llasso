## Note that at the end, we will identify the original SNPs by
## X0[common[kept[i]]], the ith column in X is actually the kept[i]
## column in X prior to elimination of repeated columns, and
## the kept[i] column is then the common[kept[i]] in the original
## bedfile before eliminating rare variants
## Later use http://genome.ucsc.edu to identify which genes have hits
## BARPLOT: https://stackoverflow.com/questions/11938293/how-to-label-a-barplot-bar-with-positive-and-negative-bars-with-ggplot2
## COLORBREW: http://colorbrewer2.org/#type=diverging&scheme=RdBu&n=5


library(ggplot2)
library(viridis)
setwd("/Users/Clauberry/Dropbox/Documents/gwas/projects/22q_new/llasso/results/hgcc")
chr = "23"
rootname = paste0("22q_caucasian-chr",chr)
datafolder = "/Users/Clauberry/Documents/gwas/data/22q/22q_files_NEW/bedfiles/"

indkept = read.table(paste0(rootname,".indices-kept"), header=TRUE)
str(indkept)
indcommon = read.table(paste0(rootname,".common"), header=TRUE)
str(indcommon)
bim = read.table(paste0(datafolder,rootname,".bim"), header=FALSE)
str(bim) ## V2 = chr22_position, V4= position
names(bim) = c("chr", "chr_position", "x", "position", "a1", "a2")
beta = read.table(paste0(rootname,".beta"), header=TRUE, sep=",")
str(beta) ##36610 obs for chr22
## need to remove the intercept:
beta = beta[2:nrow(beta),]

tol = 1e-8 ##chr10
tol = 1e-6 ##chr22,17,23
tol = (1e-4)*6 ##chr15
tol = 1e-3 ##chr1,3,4,8,9,11,12,13,16,20,21
tol = 1e-2 ##chr6
sigL = which(abs(beta$betaL)>tol | abs(beta$betaSR)>tol)
length(sigL) ##163 for chr22, need to check if it is too big

newpos = bim$position[indcommon$ind[indkept$ind[sigL]]]
newbetaL = beta$betaL[sigL]
newbetaSR = beta$betaSR[sigL]

## we want to include the first and last position too
r = range(bim$position)
if(newpos[1] != r[1]){
    newpos = c(r[1],newpos)
    newbetaL = c(0,newbetaL)
    newbetaSR = c(0,newbetaSR)}
if(newpos[length(newpos)] != r[2]){
    newpos = c(newpos,r[2])
    newbetaL = c(newbetaL,0)
    newbetaSR = c(newbetaSR,0)}

df1 = data.frame(position=newpos, beta=newbetaL)
df1$who = "Lasso"
df2 = data.frame(position=newpos, beta=newbetaSR)
df2$who = "SR"
df = rbind(df1,df2)
head(df)


pdf(paste0(rootname,".pdf"),height=4,width=2)
ggplot(df, aes(x=position, y=beta, colour = who)) +
    geom_linerange(aes(ymin = 0, ymax = beta)) +coord_flip()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    geom_hline(yintercept=0, color = "black") +
    ylim(c(-0.35,0.35))+guides(col=FALSE)+
    xlab(paste("Chromosome",chr)) + ylab("SNP effects") +
    scale_color_manual(values=c("#0571b0","#ca0020"))
dev.off()





#####################################################################################################
## Old plots
#####################################################################################################

## the problem with the following graph is that position is treated
## as categorical, so there is no proper spacing unless we manually
## add a bunch of zeroes: untractable
ggplot(df, aes(x=as.factor(position), y=beta)) +
  geom_bar(stat = "identity",aes(fill=who)) + coord_flip() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    geom_hline(yintercept=0, color = "black") +
    ylim(c(-0.35,0.35))+guides(fill=FALSE)+
    xlab("Chromosome") + ylab("SNP effects") +
    scale_fill_manual(values=c("#0571b0","#ca0020"))

## another approach trying to keep all zeroes too
## but does not work because too many points to plot
## sigL = which(beta$betaL!=0 | beta$betaSR!=0)
## length(sigL) ##163 for chr22, need to check if it is too big
## bbetaL = rep(0,length(bim$position))
## bbetaL[indcommon$ind[indkept$ind[sigL]]] = beta$betaL[sigL]
## bbetaSR = rep(0,length(bim$position))
## bbetaSR[indcommon$ind[indkept$ind[sigL]]] = beta$betaSR[sigL]
## df1 = data.frame(position=bim$position, beta=bbetaL)
## df1$who = "Lasso"
## df2 = data.frame(position=bim$position, beta=bbetaSR)
## df2$who = "SR"
## df = rbind(df1,df2)


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
