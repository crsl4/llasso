# Understanding Lasso.jl
We will follow the steps in the `test/lasso.jl` file in [here](https://github.com/simonster/Lasso.jl/blob/master/test/lasso.jl). Now in  `test-lasso.jl`
```julia
using Lasso, GLM, Distributions, GLMNet
function makeX(ρ, nsamples, nfeatures, sparse)
    Σ = fill(ρ, nfeatures, nfeatures)
    Σ[diagind(Σ)] = 1
    X = rand(MvNormal(Σ), nsamples)'
    sparse && (X[randperm(length(X))[1:round(Int, length(X)*0.95)]] = 0)
    β = [(-1)^j*exp(-2*(j-1)/20) for j = 1:nfeatures]
    (X, β)
end

randdist(::Normal, x) = rand(Normal(x))
randdist(::Binomial, x) = rand(Bernoulli(x))
randdist(::Poisson, x) = rand(Poisson(x))
function genrand(T::DataType, d::Distribution, l::GLM.Link, nsamples::Int, nfeatures::Int, sparse::Bool)
    X, coef = makeX(0.0, nsamples, nfeatures, sparse)
    y = X*coef
    for i = 1:length(y)
        y[i] = randdist(d, linkinv(l, y[i]))
    end
    (X, y)
end
function gen_penalty_factors(X,nonone_penalty_factors;sparcity=0.7)
    if nonone_penalty_factors
        penalty_factor = ones(size(X,2))
        nonone = Int(floor(size(X,2)*(1-sparcity)))
        srand(7337)
        penalty_factor[1:nonone] = rand(Float64,nonone)
        penalty_factor_glmnet = penalty_factor
    else
        penalty_factor = nothing
        penalty_factor_glmnet = ones(size(X,2))
    end
    penalty_factor, penalty_factor_glmnet
end
#dist = Normal()
#link = IdentityLink()
dist = Binomial()
link = LogitLink()
#dist = Poisson()
#link = LogLink()
sp = false #sparse
(X, y) = genrand(Float64, dist, link, 1000, 10, sp)
typeof(X)
#Array{Float64,2}
typeof(y)
#Array{Float64,1}
```
So, we have the right format of X and y for Lasso.
```julia
l = fit(LassoPath, X, y, dist, link)
```
This works!! But with my dataset, I get an error (this was fixed, see below):
```julia
using SnpArrays
using DataFrames
using Lasso, GLM, Distributions

chr22 = SnpArray("bedfiles/22q-chr22")
chr22mat=convert(Matrix{Float64},chr22)
dat = readtable("bedfiles/22q-chr22.fam", separator=' ', header=false)
y = convert(Vector,dat[:,6])
y = y-1 ## 0=controls, 1=cases
y = convert(Array{Float64,1},y)
fit(LassoPath,chr22mat,y,d=Bernoulli(),link=LogitLink())
```
So, I will fork the `Lasso.jl` package, and delete the `Lasso` folder in `.julia` with `Pkg.rm("Lasso")` and `rm -rf Lasso`, to git clone my fork: 
```julia
Pkg.clone("https://github.com/crsl4/Lasso.jl.git")
```

**This was actually undone, because I found the source of my mistake.** So, now I do not have a fork of `Lasso.jl`.

We will send a subset of the data to professor Bates for testing Lasso:
```julia
#519×450438 Array{Float64,2}
srand(1234)
ind = rand(1:519,100)
xx = chr22mat[ind,rand(1:end,10000)]
yy = y[ind]
df = DataFrame(hcat(yy,xx))
writetable("test.txt",df)
```
Now, if we want to run lasso:
```julia
using DataFrames, Lasso
dat = readtable("test.txt", header=true)
y = convert(Array{Float64,1},dat[:,1])
X = convert(Array{Float64,2},dat[:,2:end])
l=fit(LassoPath,X,y,d=Bernoulli(),link=LogitLink())
```
We get an error: this has to do with missing data. So, thanks to prof Bates (`HumanGenetics.ipynb`), we can get rid of the columns with missing data. Apparently we can do this directly also with `SnpArray`.
```julia
Pkg.status("DataFrames")
 #- DataFrames                    0.10.1
using DataFrames, Lasso
dat = readtable("test.txt", header=true, nastrings=["NaN"]);
nacol = [any(isna.(v)) for v in dat.columns];
dat = dat[:,.!nacol];
size(dat) #(100, 9126)
any(completecases(dat)) #true
y = convert(Array{Float64,1},dat[:,1])
X = convert(Array{Float64,2},dat[:,2:end])
l=fit(LassoPath,X,y,d=Bernoulli(),link=LogitLink())
```
Still get the same error:
```
WARNING: One of the predicators (columns of X) is a constant, so it can not be standardized. 
To include a constant predicator set standardize = false and intercept = false
WARNING: Array{T}(::Type{T}, m::Int) is deprecated, use Array{T}(m) instead.
Stacktrace:
 [1] depwarn(::String, ::Symbol) at ./deprecated.jl:70
 [2] Array(::Type{Float64}, ::Int64) at ./deprecated.jl:57
 [3] Lasso.NaiveCoordinateDescent{Float64,true,Array{Float64,2},Lasso.RandomCoefficientIterator,Void}(::Array{Float64,2}, ::Float64, ::Int64, ::Float64, ::Lasso.RandomCoefficientIterator, ::Void) at /Users/Clauberry/.julia/v0.6/Lasso/src/coordinate_descent.jl:40
 [4] #fit#1(::Array{Float64,1}, ::Array{Float64,1}, ::Float64, ::Int64, ::Float64, ::Void, ::Bool, ::Bool, ::Type{T} where T, ::Bool, ::Float64, ::Bool, ::Int64, ::Void, ::Array{Any,1}, ::StatsBase.#fit, ::Type{Lasso.LassoPath}, ::Array{Float64,2}, ::Array{Float64,1}, ::Distributions.Normal{Float64}, ::GLM.IdentityLink) at /Users/Clauberry/.julia/v0.6/Lasso/src/Lasso.jl:328
 [5] (::StatsBase.#kw##fit)(::Array{Any,1}, ::StatsBase.#fit, ::Type{Lasso.LassoPath}, ::Array{Float64,2}, ::Array{Float64,1}, ::Distributions.Normal{Float64}, ::GLM.IdentityLink) at ./<missing>:0 (repeats 2 times)
 [6] eval(::Module, ::Any) at ./boot.jl:235
 [7] eval_user_input(::Any, ::Base.REPL.REPLBackend) at ./REPL.jl:66
 [8] macro expansion at ./REPL.jl:97 [inlined]
 [9] (::Base.REPL.##1#2{Base.REPL.REPLBackend})() at ./event.jl:73
while loading no file, in expression starting on line 0
WARNING: Array{T}(::Type{T}, m::Int) is deprecated, use Array{T}(m) instead.
Stacktrace:
 [1] depwarn(::String, ::Symbol) at ./deprecated.jl:70
 [2] Array(::Type{Float64}, ::Int64) at ./deprecated.jl:57
 [3] build_model(::Array{Float64,2}, ::Array{Float64,1}, ::Distributions.Normal{Float64}, ::GLM.IdentityLink, ::Lasso.NaiveCoordinateDescent{Float64,true,Array{Float64,2},Lasso.RandomCoefficientIterator,Void}, ::Float64, ::Void, ::Array{Float64,1}, ::Array{Float64,1}, ::Float64, ::Int64, ::Void, ::Bool, ::Float64) at /Users/Clauberry/.julia/v0.6/Lasso/src/Lasso.jl:234
 [4] #fit#1(::Array{Float64,1}, ::Array{Float64,1}, ::Float64, ::Int64, ::Float64, ::Void, ::Bool, ::Bool, ::Type{T} where T, ::Bool, ::Float64, ::Bool, ::Int64, ::Void, ::Array{Any,1}, ::StatsBase.#fit, ::Type{Lasso.LassoPath}, ::Array{Float64,2}, ::Array{Float64,1}, ::Distributions.Normal{Float64}, ::GLM.IdentityLink) at /Users/Clauberry/.julia/v0.6/Lasso/src/Lasso.jl:332
 [5] (::StatsBase.#kw##fit)(::Array{Any,1}, ::StatsBase.#fit, ::Type{Lasso.LassoPath}, ::Array{Float64,2}, ::Array{Float64,1}, ::Distributions.Normal{Float64}, ::GLM.IdentityLink) at ./<missing>:0 (repeats 2 times)
 [6] eval(::Module, ::Any) at ./boot.jl:235
 [7] eval_user_input(::Any, ::Base.REPL.REPLBackend) at ./REPL.jl:66
 [8] macro expansion at ./REPL.jl:97 [inlined]
 [9] (::Base.REPL.##1#2{Base.REPL.REPLBackend})() at ./event.jl:73
while loading no file, in expression starting on line 0
ERROR: ArgumentError: start and stop must be finite, got NaN and NaN
Stacktrace:
 [1] _linspace(::Float64, ::Float64, ::Int64) at ./twiceprecision.jl:342
 [2] linspace(::Float64, ::Float64, ::Int64) at ./twiceprecision.jl:338
 [3] computeλ(::Array{Float64,1}, ::Float64, ::Float64, ::Int64, ::Void) at /Users/Clauberry/.julia/v0.6/Lasso/src/Lasso.jl:213
 [4] build_model(::Array{Float64,2}, ::Array{Float64,1}, ::Distributions.Normal{Float64}, ::GLM.IdentityLink, ::Lasso.NaiveCoordinateDescent{Float64,true,Array{Float64,2},Lasso.RandomCoefficientIterator,Void}, ::Float64, ::Void, ::Array{Float64,1}, ::Array{Float64,1}, ::Float64, ::Int64, ::Void, ::Bool, ::Float64) at /Users/Clauberry/.julia/v0.6/Lasso/src/Lasso.jl:242
 [5] #fit#1(::Array{Float64,1}, ::Array{Float64,1}, ::Float64, ::Int64, ::Float64, ::Void, ::Bool, ::Bool, ::Type{T} where T, ::Bool, ::Float64, ::Bool, ::Int64, ::Void, ::Array{Any,1}, ::StatsBase.#fit, ::Type{Lasso.LassoPath}, ::Array{Float64,2}, ::Array{Float64,1}, ::Distributions.Normal{Float64}, ::GLM.IdentityLink) at /Users/Clauberry/.julia/v0.6/Lasso/src/Lasso.jl:332
 [6] (::StatsBase.#kw##fit)(::Array{Any,1}, ::StatsBase.#fit, ::Type{Lasso.LassoPath}, ::Array{Float64,2}, ::Array{Float64,1}, ::Distributions.Normal{Float64}, ::GLM.IdentityLink) at ./<missing>:0 (repeats 2 times)
```
So, maybe we need to get rid of constant columns:
```julia
ctecol = [std(v)==0 for v in dat.columns]
dat = dat[:,.!ctecol];  
size(dat) #(100, 2164)
y2 = convert(Array{Float64,1},dat[:,1])
X2 = convert(Array{Float64,2},dat[:,2:end])
ll=fit(LassoPath,X2,y2,d=Bernoulli(),link=LogitLink())
```
I get a different error now:
```
ERROR: MethodError: no method matching fit!(::Lasso.LassoPath{GLM.LinearModel{GLM.LmResp{Array{Float64,1}},Lasso.NaiveCoordinateDescent{Float64,true,Array{Float64,2},Lasso.RandomCoefficientIterator,Void}},Float64}; irls_tol=1.0e-7, d=Distributions.Bernoulli{Float64}(p=0.5), link=GLM.LogitLink())
```
This does not make sense!

**Error fix** This was my mistake. I should call without the `d=` and `link=`, like this: `fit(LassoPath,X,y,Bernoulli(),LogitLink())`.

# Using SparseRegression.jl
Documentation [here](http://joshday.github.io/SparseRegression.jl/latest/smodel.html)

First a toy example:
```julia
x, y = randn(100, 10000), randn(100)
obs = Obs(x, y)
s = SModel(obs)
# Learn the model using Proximal Gradient Method
# - step size: 0.01 (default=1)
# - maximum of 500 iterations
# - convergence criteria: norm(β - βold) < 1e-6
learn!(s, ProxGrad(obs, .01), MaxIter(500), Converged(coef))
```

```julia
Pkg.clone("https://github.com/joshday/SparseRegression.jl.git")
using SparseRegression
using DataFrames
dat = readtable("test.txt", header=true, nastrings=["NaN"]);
nacol = [any(isna.(v)) for v in dat.columns];
dat = dat[:,.!nacol];
size(dat) #(100, 9126)
any(completecases(dat)) #true
y = convert(Array{Float64,1},dat[:,1])
X = convert(Array{Float64,2},dat[:,2:end])
obs = Obs(X, y)
s = SModel(obs, LogisticRegression(), L1Penalty())
learn!(s, ProxGrad(obs,0.001), MaxIter(10000), Converged(coef))
```
Converges! But now I need to identify the SNPs that have beta>0:
```julia
find(x -> x >0,s.β)
#12-element Array{Int64,1}:
#  564
#  819
#  987
# 1047
# 1150
# 5131
# 5607
# 5739
# 7228
# 7880
# 8059
# 8505
```
 Also, need to look at the code to see what it is doing and how to add dependence in any way. Probably we also want to check the code to see what it is doing.
 
# Using GLMNet
This is a wrapper for fortran code. It is strange that you cannot decide if you want Lasso or Elastic Net.
```julia
using GLMNet
y = collect(1:100) + randn(100)*10;
X = [1:100 (1:100)+randn(100)*5 (1:100)+randn(100)*10 (1:100)+randn(100)*20];
path = glmnet(X, y)
cv = glmnetcv(X, y)
```
We don't like this package because it is not being updated, and it is kept for julia 0.5. Also, no documentation. And we need a different $y$ for the response of logistic regression (it has to be a matrix!).