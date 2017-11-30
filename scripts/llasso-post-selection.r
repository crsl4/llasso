bedfile = "22q-chr22"
load(paste0(bedfile, ".Rda"))

library(selectiveInference)
outL = fixedLassoInf(X, y, beta_hatL, lambdaL, family = "binomial")
outSR = fixedLassoInf(X, y, beta_hatSR, lambdaSR, family = "binomial")
#out=fixedLogitLassoInf(X,y,beta,lambda,alpha=alpha, type="partial", tol.beta=tol.beta, tol.kkt=tol.kkt,
#                           gridrange=gridrange, bits=bits, verbose=verbose,this.call=this.call)

