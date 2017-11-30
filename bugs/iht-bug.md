# Bug in logistic IHT

```julia
using DataFrames, IHT
dat = readtable("test.txt", header=true)
y = convert(Array{Float64,1},dat[:,1]) ##1st column: response
X = convert(Array{Float64,2},dat[:,2:end]) ##10,000 SNPs
output1 = @time L0_log(X, y, 1)
#0.004151 seconds (389 allocations: 1.165 MiB)
#IHT results:

#Compute time (sec):   0.0041
#Final loss:           0.6860874
#Iterations:           4
#IHT estimated 1 nonzero coefficients.
#1×2 DataFrames.DataFrame
#│ Row │ Predictor │ Estimated_β │
#├─────┼───────────┼─────────────┤
#│ 1   │ 8550      │ 0.331812    │
output2 = @time L0_log(X, y, 2)
#0.006286 seconds (445 allocations: 1.176 MiB)
#IHT results:

#Compute time (sec):   0.0062
#Final loss:           0.6840915
#Iterations:           4
#IHT estimated 2 nonzero coefficients.
#2×2 DataFrames.DataFrame
#│ Row │ Predictor │ Estimated_β │
#├─────┼───────────┼─────────────┤
#│ 1   │ 4307      │ 0.193256    │
#│ 2   │ 8550      │ 0.194642    │
output3 = @time L0_log(X, y, 3)
#ERROR: BoundsError: attempt to access 519×3 Array{Float64,2} at index [Base.OneTo(519), Base.OneTo(4)]
#Stacktrace:
# [1] throw_boundserror(::Array{Float64,2}, ::Tuple{Base.OneTo{Int64},Base.OneTo{Int64}}) at ./abstractarray.jl:433
# [2] checkbounds at ./abstractarray.jl:362 [inlined]
# [3] copy!(::Array{Float64,2}, ::SubArray{Float64,2,Array{Float64,2},Tuple{Base.Slice{Base.OneTo{Int64}},Array{Int64,1}},false}) at ./multidimensional.jl:803
# [4] #L0_log#90(::IHT.IHTLogVariables{Float64,Array{Float64,1}}, ::Float64, ::Float64, ::Float64, ::Float64, ::Float64, ::Int64, ::Int64, ::Bool, ::Bool, ::IHT.#L0_log, ::Array{Float64,2}, ::Array{Float64,1}, ::Int64) at /Users/Clauberry/.julia/v0.6/IHT/src/log.jl:343
# [5] L0_log(::Array{Float64,2}, ::Array{Float64,1}, ::Int64) at /Users/Clauberry/.julia/v0.6/IHT/src/log.jl:207
```

This error was fixed with commit `238c44c912c886c574ebf0f10fe39d85c2fec0d1`.