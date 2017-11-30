# Bug 1 in SparseRegression

I was trying to code cross validation, when I realized that no matter which subset of the data I selected, I would always get all $\beta$ positive.
This could be a bug in the optimization?
I don't get the same results when I run with `GLMNet`.

```julia
using DataFrames, SparseRegression
dat = readtable("test.txt", header=true)
y = convert(Array{Float64,1},dat[:,1]) ##1st column: response
X = convert(Array{Float64,2},dat[:,2:end]) ##10,000 SNPs
obs = Obs(X, y)
s = SModel(obs, LogisticRegression(), L1Penalty())
@time learn!(s, ProxGrad(obs,0.001), MaxIter(10000), Converged(coef))
#23.880103 seconds (106.87 k allocations: 765.675 MiB, 0.69% gc time)
#SparseRegression.SModel{LossFunctions.LogitMarginLoss,PenaltyFunctions.L1Penalty}
 #> β        : [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
 #> λ factor : [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1  …  0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
 #> Loss     : LogitMarginLoss
 #> Penalty  : L1Penalty
find(x -> x != 0,s.β) ##103 SNPs for λ=0.1
mean(s.β[find(x -> x != 0,s.β)])
find(x -> x < 0,s.β) ## 0
```

Using  `GLMNet`:

```julia
using GLMNet, Distributions
y1 = Bool.(y)
y2 = .!y1
yy = [y1 y2]
path = @time glmnet(X, convert(Matrix{Float64}, yy), Binomial())
```
In `GLMNet`, the biggest $\lambda$ used is 0.089, and with this lambda, all betas are equal to 0 (no covariates in the model). `SparseRegression` is using $\lambda=0.1$ (even bigger $\lambda$), so I would expect all betas to be close to zero (but the mean is equal to 0.013, not that close to zero).
Could there by any rounding error somewhere, or something, making these estimated coefficients further from zero?

When trying with a different $\lambda$ (the best $\lambda$ by CV in `GLMNet` which yields 19 covariates in the model), we still get all positive betas for 170 covariates:

```julia
cv = @time glmnetcv(X, yy, Binomial())
#11.709583 seconds (299.73 k allocations: 623.351 MiB, 1.03% gc time)
#Logistic GLMNet Cross Validation
#100 models for 10000 predictors in 10 folds
#Best λ 0.062 (mean loss 1.380, std 0.008)
obs = Obs(X, y)
s = SModel(obs, LogisticRegression(), L1Penalty(), fill(0.062,size(X,2)))
@time learn!(s, ProxGrad(obs,0.001), MaxIter(10000), Converged(coef))
find(x -> x != 0,s.β) ##170 SNPs for λ=0.1
mean(s.β[find(x -> x != 0,s.β)])
find(x -> x < 0,s.β) ## 0
```

This bug was fixed when using a response -1,1 instead of 0,1.

# Bug 2 in SparseRegression
It would be great if we could use a `SnpArray` object as input for an `SModel`, instead of having to convert to a matrix which can be very heavy: @klkeys @Hua-Zhou
```julia
using DataFrames, SparseRegression, SnpArrays, Distributions
chr = SnpArray("SNP_data29a_missing"; people = 212, snps = 253141)
typeof(chr) #SnpArrays.SnpArray{2}
SnpArray <: AbstractArray #true
y = convert(Vector{Float64},randn(212).>0)
y[y.==0] = -1 # SparseRegression needs outcome -1,1
isa(chr,AbstractArray) #true
s = SModel(chr,y, LogitDistLoss(), L1Penalty())
@time learn!(s)
## ERROR: MethodError: no method matching *(::Tuple{Bool,Bool}, ::Float64)
```
So, the multiplication operations does not work. It would be nice it they did, as we would avoid converting to float matrix a big dataset.

Josh says that the following would fix the problem:
```julia
Base.A_mul_B!(x::Matrix{Float64}, a::SnpArray, b::Vector{Float64}) = (x[:] = a * b)
Base.At_mul_B!(x::Matrix{Float64}, a::SnpArray, b::Vector{Float64}) = (x[:] = x'b)
```
But this does not work! We get the same error.

## Old:
This cannot be done even when` Obs` is defined for `AbstractArray`: `function Obs(x::AbstractArray, y::AbstractArray)`
```julia
using DataFrames, SparseRegression, SnpArrays, Distributions
chr = SnpArray("SNP_data29a_missing"; people = 212, snps = 253141)
typeof(chr) #SnpArrays.SnpArray{2}
SnpArray <: AbstractArray #true
y = convert(Vector{Float64},randn(212).>0)
y[y.==0] = -1 # SparseRegression needs outcome -1,1
isa(chr,AbstractArray) #true
obs = Obs(chr, y)
#ERROR: TypeError: Type: in T, expected T<:Number, got Type{Tuple{Bool,Bool}}
#Stacktrace:
# [1] SparseRegression.Obs(::SnpArrays.SnpArray{2}, ::Array{Float64,1}) at #/Users/Clauberry/.julia/v0.6/SparseRegression/src/obs.jl:18
```
Later, I discovered that the package had been updated, so I did `Pkg.update`. Now, we do not need the `Obs` definition.

