# ----------------------------------------------- #
# post selection inference of lasso coefficients
# ----------------------------------------------- #
## our julia functions are still under development, so we
## will use the R package selectiveInference
bedfile = "22q-chr22"
outfile = paste0(bedfile,".beta")

@rput X
@rput y
if(resL)
    beta_hatL = convert(Array,betaL)[:,1]
    unshift!(beta_hatL,interceptL)
    @rput beta_hatL
    @rput lambdaL
end
if(resSR)
    beta_hatSR = convert(Array,betaSR)[:,1]
    unshift!(beta_hatSR,interceptSR)
    @rput beta_hatSR
    @rput lambdaSR
end
@rput bedfile
@rput resL
@rput resSR

R"""
if(resL && !resSR){
    save(X,y,beta_hatL,lambdaL,file=paste0(bedfile,".Rda"))
}else if(resSR && !resL){
    save(X,y, beta_hatSR,lambdaSR,file=paste0(bedfile,".Rda"))
}else if(resL && resSR){
    save(X,y,beta_hatL,lambdaL, beta_hatSR,lambdaSR,file=paste0(bedfile,".Rda"))
}
"""

load("22q-chr22.Rda")
library(selectiveInference)
outL = fixedLassoInf(X,y,beta_hatL,lambdaL,family="binomial")
outSR = fixedLassoInf(X,y,beta_hatSR,lambdaSR,family="binomial")
#out=fixedLogitLassoInf(X,y,beta,lambda,alpha=alpha, type="partial", tol.beta=tol.beta, tol.kkt=tol.kkt,
#                           gridrange=gridrange, bits=bits, verbose=verbose,this.call=this.call)

