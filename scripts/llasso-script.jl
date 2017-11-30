# Julia script to run multiple steps on post selection
# for a given bedfile
# missing for efficiency:
# - change the post selection julia functions to the ones in R CRAN
#   now in llasso-post-selection.r
# - can we modify SnpArray so that we do not have to convert to matrix?
# Claudia October 2017

include("llasso-functions.jl")
using JLD
datafolder = "/Users/Clauberry/Documents/gwas/data/22q/22q_files_NEW/bedfiles/"
bedfile = "22q-chr22"
## to save output
logfile = string(bedfile,".log")
outfile = string(bedfile,".beta")

## Here we read the BED file, and remove rare variants
chr = readBedfile(datafolder, bedfile);

# -------------------------- #
# choosing the best lambda
# -------------------------- #
info("Choosing the best lambda with cross validation")
# Lasso.jl will fit 100 models, so we do CV
lassoconv = true
## initialize outside of try-catch
dfL = nothing
bestlL = nothing
f = nothing
try
    dfL, bestlL, f = chooseLambdaLasso(chr, datafolder, bedfile)
catch
    warn("Lasso fit did not converge in cross validation, so we do not have best lambda")
    lassoconv = false
end

## if lasso did not converge to give us a vector of lambdas to test,
## we use one based on experience
lambdas = lassoconv ? f.λ : collect(0.000854897:0.00001:0.002)
dfSR, bestlSR = chooseLambdaSparseRegression(chr, datafolder, bedfile, lambdas)

@save string(bedfile,"-crossval.jld")

# ----------------------------- #
# penalized logistic lasso fit
# ----------------------------- #
info("Cleaning up data for penalized likelihood")
## Until this point, we have run CV on Lasso.jl and SparseRegression.jl,
## now we will run Lasso.jl and SparseRegression.jl on the full data
## with the best lambda.
## For chromosome 22, we have 45k common variants, so we can convert the
## whole SnpArray
X,y = convertBedfile(chr, datafolder,bedfile)


df = DataFrame()
## now we run Lasso.jl, which will fit by default 100 models
## we will choose at the end the corresponding the best lambda we found
info("Running penalized likelihood with Lasso.jl")
resL,fall = fitRepeatedLlasso(X, y)
betaL = nothing
lambdaL = nothing
interceptL = nothing
if(resL)
    @show fall
    ind = find(x->isapprox(x,bestlL,atol=0.001),fall.λ)
    isempty(ind) && warn("We could not find the Lasso model with the best lambda")
    if(length(ind)>1)
        ind = ind[end]
    end
    betaL = fall.coefs[:,ind]
    lambdaL = fall.λ[ind][1]
    interceptL = fall.b0[ind][1]
    ## writing output
    yy = convert(Array,betaL)
    df[:betaL] = yy[:,1]
    writetable(outfile,df)
    f = open(logfile,"w")
    write(f,string("LASSO implementation of best lambda with Lasso.jl\nlambda = ",lambdaL[1],"\nintercept = ",interceptL[1],"\nEstimated betas in file ",outfile))
    close(f)
else
    warn("Lasso fit never converged for best lambda")
end


## Now we will run SparseRegression.jl for the best lambda
info("Running penalized likelihood with SparseRegression.jl")
X2 = hcat(fill(1.0,size(X,1)),X) ## we need to add the intercept, because SparseReg does not do this
y[y.==0] = -1 ##SparseRegression needs response -1,1
sall = SModel(X2,y, LogitDistLoss(), L1Penalty(),fill(bestlSR[1],size(X2,2)))
resSR = true
try
    learn!(sall, strategy(ProxGrad(sall,0.001), MaxIter(10000), Converged(coef)));
catch
    warn("SparseRegression fit did not converge for best lambda")
    resSR = false
end
## we need here beta, lambda, intercept from SR
betaSR = sall.β[2:end]
lambdaSR = bestlSR[1]
interceptSR = sall.β[1]
## writing output
yy = convert(Array,betaSR)
df[:betaSR] = yy[:,1]
writetable(outfile,df)
f = open(logfile,"a")
write(f,string("\nLASSO implementation of best lambda with SparseRegression.jl\nlambda = ",lambdaSR[1],"\nintercept = ",interceptSR[1],"\nEstimated betas in file ",outfile))
close(f)

@save string(bedfile,"-llasso.jld")

# ----------------------------------------------- #
# penalized logistic model with IHT penalty
# ----------------------------------------------- #
## This penalty does not force beta coefficients to be small
ll = sqrt(log(size(X,2))/size(X,1))
info("Default lambda used in IHT: $(ll)")
success = true
cvout = nothing
try
    cvout = @time cv_log(X2,y,tol=0.01, tolG=0.01, tolrefit=0.01)
catch
    success = false
end
m = min(sum(beta_hatL .!= 0)-1, sum(beta_hatSR .!= 0)-1)
df0 = DataFrame()
for i in 1:m
    output = @time L0_log(X2, y, i)
    df0[Symbol(eval(string("IHT",i)))] = output.beta[2:end]
end

## writing output to file --------------------------------------------
df = [df df0]
writetable(outfile,df)


## Removed part for RCall because it does not work in HGCC
# using RCall
# @rput X
# @rput y
# if(resL)
#     beta_hatL = convert(Array,betaL)[:,1]
#     unshift!(beta_hatL,interceptL)
#     @rput beta_hatL
#     @rput lambdaL
# end
# if(resSR)
#     beta_hatSR = convert(Array,betaSR)[:,1]
#     unshift!(beta_hatSR,interceptSR)
#     @rput beta_hatSR
#     @rput lambdaSR
# end
# @rput bedfile
# @rput resL
# @rput resSR

# R"""
# if(resL && !resSR){
#     save(X,y,beta_hatL,lambdaL,file=paste0(bedfile,".Rda"))
# }else if(resSR && !resL){
#     save(X,y, beta_hatSR,lambdaSR,file=paste0(bedfile,".Rda"))
# }else if(resL && resSR){
#     save(X,y,beta_hatL,lambdaL, beta_hatSR,lambdaSR,file=paste0(bedfile,".Rda"))
# }
# """