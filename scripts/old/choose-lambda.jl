## julia script to do pseudo cross validation on the 22-q deletion data
## (subset) to choose the best lambda for LASSO
## Note: we cannot open this in emacs because of β.
## This is not really CV because we are not splitting the data,
## just reusing the models fit by Lasso.jl
## Claudia September 2017
## NOTE: not used anymore. True cross validation is implemented in
##       file llasso-functions.jl
## Claudia October 2017

using SnpArrays, Lasso, DataFrames, Distributions

## Read in the data
datafolder = "/Users/Clauberry/Documents/gwas/data/22q/22q_files_NEW/bedfiles/"
bedfile = "22q-chr22"
chr = SnpArray(string(datafolder,bedfile))
chr = chr[:, maf .≥ 0.05]
r = 12345
srand(r)
ind = sample(1:size(chr,2),10000,replace=false)
chr22sub = chr[:,ind]
@time X=convert(Matrix{Float64},chr22sub,impute=true)
dat = readtable(string(datafolder,"bedfiles/22q-chr22.fam"), separator=' ', header=false)
y = convert(Array{Float64,1},dat[:,6])
y = y-1
# get rid of constant columns, indispensable for convergence
#ctecol = [std(X[:,i])==0 for i in 1:size(X,2)]
#sum(ctecol)
#X2 = X[:,.!ctecol];
# fit 100 models by default
f2=@time fit(LassoPath,X,y,Bernoulli(),LogitLink())

# now compute score per lambda value:
scores = fill(0.0,size(f2.coefs,2))

for i in 1:size(f2.coefs,2)
    a = exp.(X*f2.coefs[:,i])
    p1 = a./(1+a)
    yest = rand.(Bernoulli.(p1))
    scores[i] = sum(abs.(y-yest))
end

bestmod = find(x->x==minimum(scores),scores)
## best lambda:
f2.λ[bestmod]
f2.coefs[:,bestmod] ## 427 covariates




