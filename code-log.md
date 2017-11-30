# To do
  - leave llasso-script running in HGCC
  - create new bed files (with plink, or julia? maybe plink to do it just once with mike's list of IDs) keeping only caucasian
  - run julia script for all chromosomes: parallelize
  - check with rich/jay how to identify in which genes the significant SNPs are
  - we need to add confounding covariates at some point: sex, age? not sure, as they are not significant

Later:
  - meet with Rich to learn about amazon web
  - make functions more efficient: see `llasso-script.jl`, put functions in PostSelectionInference

# Scripts
- `llasso-script.jl`: for a given bedfile name, it runs:
  - cross validation with Lasso.jl to find the best lambda
  - cross validation with SparseRegression.jl to find the best lambda
  - penalized likelihodd with Lasso.jl with best lambda
  - penalized likelihood with SparseRegression.jl with best lambda
  - penalized likelihood with IHT.jl
  - saves all betas in *.beta file to be read later and interpreted
- `llasso-post-sel-preparation.jl`: will read the output files from  `llasso-script.jl`, create an Rda file
- `llasso-post-selection.r`: will run post selection methods in R (still under development)
- `llasso-interpret.r`: will read the *.beta file created by `llasso-script.jl` and will identify the significant variants (still under development)

# Results
- `test1` has the results when running the whole script in julia (before removing the RCall part that does not work in HGCC). It also uses the whole dataset (not just europeans)
- `test2`: run the whole `llasso-script.jl`, but after removing the RCall part that cannot be run in HGCC.

# HGCC

We need to install the necessary julia packages in HGCC in an interactive job:
```
ssh csolislemus@hgcc.genetics.emory.edu
qlogin -q i.q
module load julia
module load R
which julia
##/sw/hgcc/Pkgs/julia/0.6.0/bin/julia
julia
Pkg.add("DataFrames")
## INFO: Initializing package repository /home/csolislemus/.julia/v0.6
Pkg.add("Distributions")
Pkg.add("Lasso") 
Pkg.clone("https://github.com/OpenMendel/SnpArrays.jl.git")
Pkg.clone("https://github.com/joshday/SparseRegression.jl")
Pkg.clone("https://github.com/klkeys/PLINK.jl.git")
Pkg.clone("https://github.com/klkeys/RegressionTools.jl.git")
Pkg.clone("https://github.com/klkeys/IHT.jl.git")
Pkg.add("RCall")
```
I get an error with `RCall`: "LoadError: Could not find R installation. Either set the "R_HOME" environmental variable, or ensure the R executable is available in "PATH".", even when I can call R:
```
$ which R
/sw/hgcc/Pkgs/R/3.3.2/bin/R
```
But inside R, we get:
```R
R.home()
##/sw/hgcc/Pkgs/R/3.3.2/lib64/R
```

From [here](https://github.com/JuliaInterop/RCall.jl/blob/master/docs/src/installation.md), we have the following instructions:
```
Firstly, try setting the R_HOME environmental variable to the location of your R installation, which can be found by running R.home() from within R. This can be set in your ~/.juliarc.jl file via the ENV global variable, e.g.

ENV["R_HOME"] = "/sw/hgcc/Pkgs/R/3.3.2/lib64/R"
```
I have done this, but I keep getting the same error. In [this issue](https://github.com/JuliaInterop/RCall.jl/issues/198), people mention that this is a problem when running in a cluster, and they suggest the solution of:
```
the solution was to rebuild the R installation with --enable-R-shlib and to add $R_HOME/lib to the LD_LIBRARY_PATH.
```
So, I need to open another HGCC ticket. We are waiting for this.

### Running script without RCall
We modified `llasso-script.jl` to eliminate the part with `RCall`.