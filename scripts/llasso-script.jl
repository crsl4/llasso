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
X,y = convertBedfile(chr, datafolder,bedfile);

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
    betaL = convert(Array,fall.coefs[:,ind])
    betaL = betaL[:,]
    lambdaL = fall.λ[ind][1]
    interceptL = fall.b0[ind][1]
    unshift!(betaL,interceptL)
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
## we keep intercept in beta to save this way
betaSR = sall.β
lambdaSR = bestlSR[1]
## writing output
yy = convert(Array,betaSR)
df[:betaSR] = yy[:,1]
writetable(outfile,df)
f = open(logfile,"a")
write(f,string("\nLASSO implementation of best lambda with SparseRegression.jl\nlambda = ",lambdaSR[1],"\nintercept = ",betaSR[1],"\nEstimated betas in file ",outfile))
close(f)



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
if(success)
    f = open(logfile,"a")
    write(f,"\n")
    write(f,string(cvout))
    close(f)
end
m = min(sum(betaL .!= 0)-1, sum(betaSR .!= 0)-1,50)
## 50 from phase transition plot
df0 = DataFrame()
for i in 1:m
    output = @time L0_log(X2, y, i)
    df0[Symbol(eval(string("IHT",i)))] = output.beta
end

## writing output to file --------------------------------------------
df = [df df0]
writetable(outfile,df)
@save string(bedfile,"-llasso.jld")
#@load string(bedfile,"-llasso.jld")

## Now that RCall works in HGCC, we can call:
#include("llasso-post-sel-preparation.jl")
## no, it does not really work

## ------------------------------------------------------------------
## Overfitting the data: fitting the logistic model with the candidate
## SNPs (not great, but we have so little data that we cannot do 
## data splitting)

## We include sex as the only covariate that was significant
## in the exploration
cov = readtable("/Users/Clauberry/Documents/gwas/data/22q/22q_files_NEW/justfinalset/22q_all.covariates.txt", separator='\t')
## sizes don't really match (526 in cov, and 519 in chr), but we
## ignore this for now as we will use other data
sex = cov[:Sex][1:519]

sumbeta = sum(Array(df),2)
df[:sumbeta] = sumbeta[:,1]
df[:rownum] = collect(1:1:length(df[:betaL])) ## to keep track of indices
sort!(df, cols=[:sumbeta], rev=true)

## candidate SNPs, in "kept" indices
## using 50 because of phase transition plot
candidate = Array(df[:rownum][1:50])
dfind = DataFrame(ind=candidate)
writetable(string(bedfile,".candidate"),dfind)
Xsub = X[:,candidate]
y[y.==-1] = 0
dat = DataFrame(Xsub)
dat[:y] = y
dat[:sex] = sex

model1 = glm(allvarsformula(names(dat),:y), dat, Normal(), IdentityLink())
f = open(string(bedfile,".glm-output"),"w")
write(f,string(model1))
close(f)

cc = coef(model1)
se = stderr(model1)
zz = cc ./ se
dat2 = DataFrame(beta=cc,se=se,z=zz,pvalue=2.0 * ccdf.(Normal(), abs.(zz)))
dat2[:candidate] = push!(unshift!(candidate, 0),0)
writetable(string(bedfile,".glm"),dat2)
    