# Bug in Lasso.jl

I am not able to set one given value of lambda for the Lasso fit.
I am not sure if the lambda parameter is one lambda per covariate, or just a vector of lambdas to try (if we want only one lambda, we would still need to put in an array, because an array is expected). Either way, I get an error.
```julia
using DataFrames, Lasso
dat = readtable("test.txt", header=true)
y = convert(Array{Float64,1},dat[:,1]) ##1st column: response
X = convert(Array{Float64,2},dat[:,2:end]) ##10,000 SNPs
f2=@time fit(LassoPath,X,y,Bernoulli(),LogitLink())
## 11.347873 seconds (6.19 M allocations: 363.455 MiB, 1.41% gc time)
##Bernoulli LassoPath (100 solutions for 10000 predictors in 1387 iterations):
##                 λ    pct_dev ncoefs
##  [1]    0.0893361        0.0      0
##  [2]    0.0852756 0.00204586      1
##  [3]    0.0813997 0.00391051      1
##  [4]       0.0777 0.00655943      3
##  [5]    0.0741684  0.0114346      5
##  [6]    0.0707973  0.0185944      8
##  [7]    0.0675795  0.0293905     12
##  [8]    0.0645079  0.0414231     16
##  [9]    0.0615759  0.0549327     19
## [10]    0.0587772  0.0695958     24
## ...
## [92]   0.00129611   0.976892    416
## [93]    0.0012372   0.977948    419
## [94]   0.00118097   0.978963    420
## [95]   0.00112729   0.979923    423
## [96]   0.00107606   0.980836    418
## [97]   0.00102715   0.981715    422
## [98]  0.000980462    0.98255    421
## [99]  0.000935899   0.983345    420
##[100]  0.000893361   0.984107    427

## using the smalles lambda:
f3=@time fit(LassoPath,X,y,Bernoulli(),LogitLink(),λ=[0.000893361],nλ=1)
#ERROR: maximum number of coefficients 1038 exceeded at λ = 0.000893361 (λωj=0.000893361)
#Stacktrace:
 #[1] cycle!(::Lasso.SparseCoefficients{Float64}, ::Lasso.NaiveCoordinateDescent{Float64,true,Array{Float64,2},Lasso.RandomCoefficientIterator,Void}, ::Float64, ::Bool) at /Users/Clauberry/.julia/v0.6/Lasso/src/coordinate_descent.jl:237
 #[2] cdfit!(::Lasso.SparseCoefficients{Float64}, ::Lasso.NaiveCoordinateDescent{Float64,true,Array{Float64,2},Lasso.RandomCoefficientIterator,Void}, ::Float64, ::Symbol) at /Users/Clauberry/.julia/v0.6/Lasso/src/coordinate_descent.jl:596
 #[3] #fit!#23(::Bool, ::Int64, ::Int64, ::Float64, ::Float64, ::Symbol, ::Float64, ::Function, ::Lasso.LassoPath{GLM.GeneralizedLinearModel{GLM.GlmResp{Array{Float64,1},Distributions.Bernoulli{Float64},GLM.LogitLink},Lasso.NaiveCoordinateDescent{Float64,true,Array{Float64,2},Lasso.RandomCoefficientIterator,Void}},Float64}) at /Users/Clauberry/.julia/v0.6/Lasso/src/coordinate_descent.jl:701
 #[4] (::StatsBase.#kw##fit!)(::Array{Any,1}, ::StatsBase.#fit!, ::Lasso.LassoPath{GLM.GeneralizedLinearModel{GLM.GlmResp{Array{Float64,1},Distributions.Bernoulli{Float64},GLM.LogitLink},Lasso.NaiveCoordinateDescent{Float64,true,Array{Float64,2},Lasso.RandomCoefficientIterator,Void}},Float64}) at ./<missing>:0
 #[5] #fit#1(::Array{Float64,1}, ::Array{Float64,1}, ::Float64, ::Int64, ::Float64, ::Array{Float64,1}, ::Bool, ::Bool, ::Type{T} where T, ::Bool, ::Float64, ::Bool, ::Int64, ::Void, ::Array{Any,1}, ::StatsBase.#fit, ::Type{Lasso.LassoPath}, ::Array{Float64,2}, ::Array{Float64,1}, ::Distributions.Bernoulli{Float64}, ::GLM.LogitLink) at /Users/Clauberry/.julia/v0.6/Lasso/src/Lasso.jl:338
 #[6] (::StatsBase.#kw##fit)(::Array{Any,1}, ::StatsBase.#fit, ::Type{Lasso.LassoPath}, ::Array{Float64,2}, ::Array{Float64,1}, ::Distributions.Bernoulli{Float64}, ::GLM.LogitLink) at ./<missing>:0

## one lambda per covariate:
f3=@time fit(LassoPath,X,y,Bernoulli(),LogitLink(),λ=fill(0.000893361,size(X,2)),nλ=1)
#ERROR: maximum number of coefficients 1038 exceeded at λ = 0.000893361 (λωj=0.000893361)
#Stacktrace:
 #[1] cycle!(::Lasso.SparseCoefficients{Float64}, ::Lasso.NaiveCoordinateDescent{Float64,true,Array{Float64,2},Lasso.RandomCoefficientIterator,Void}, ::Float64, ::Bool) at /Users/Clauberry/.julia/v0.6/Lasso/src/coordinate_descent.jl:237
 #[2] cdfit!(::Lasso.SparseCoefficients{Float64}, ::Lasso.NaiveCoordinateDescent{Float64,true,Array{Float64,2},Lasso.RandomCoefficientIterator,Void}, ::Float64, ::Symbol) at /Users/Clauberry/.julia/v0.6/Lasso/src/coordinate_descent.jl:596
 #[3] #fit!#23(::Bool, ::Int64, ::Int64, ::Float64, ::Float64, ::Symbol, ::Float64, ::Function, ::Lasso.LassoPath{GLM.GeneralizedLinearModel{GLM.GlmResp{Array{Float64,1},Distributions.Bernoulli{Float64},GLM.LogitLink},Lasso.NaiveCoordinateDescent{Float64,true,Array{Float64,2},Lasso.RandomCoefficientIterator,Void}},Float64}) at /Users/Clauberry/.julia/v0.6/Lasso/src/coordinate_descent.jl:701
 #[4] (::StatsBase.#kw##fit!)(::Array{Any,1}, ::StatsBase.#fit!, ::Lasso.LassoPath{GLM.GeneralizedLinearModel{GLM.GlmResp{Array{Float64,1},Distributions.Bernoulli{Float64},GLM.LogitLink},Lasso.NaiveCoordinateDescent{Float64,true,Array{Float64,2},Lasso.RandomCoefficientIterator,Void}},Float64}) at ./<missing>:0
 #[5] #fit#1(::Array{Float64,1}, ::Array{Float64,1}, ::Float64, ::Int64, ::Float64, ::Array{Float64,1}, ::Bool, ::Bool, ::Type{T} where T, ::Bool, ::Float64, ::Bool, ::Int64, ::Void, ::Array{Any,1}, ::StatsBase.#fit, ::Type{Lasso.LassoPath}, ::Array{Float64,2}, ::Array{Float64,1}, ::Distributions.Bernoulli{Float64}, ::GLM.LogitLink) at /Users/Clauberry/.julia/v0.6/Lasso/src/Lasso.jl:338
 #[6] (::StatsBase.#kw##fit)(::Array{Any,1}, ::StatsBase.#fit, ::Type{Lasso.LassoPath}, ::Array{Float64,2}, ::Array{Float64,1}, ::Distributions.Bernoulli{Float64}, ::GLM.LogitLink) at ./<missing>:0
```

The error is strange because for that exact value of lambda, the initial fit would include 427 covariates in the model (less than the limit of 1038)