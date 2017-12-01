# To do
  - test llasso-submit.sh
  - create new bed files (with plink, or julia? maybe plink to do it just once with mike's list of IDs) keeping only caucasian
  - modify the scripts so that they take a number (the chromosome) as input
  - run julia script for all chromosomes: parallelize
  - check with rich/jai how to identify in which genes the significant SNPs are
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
So, I need to open another HGCC ticket. We are waiting for this. It was fixed by changing the version of R (needed 3.4+). It is ok that the two scripts are separated, but I will run `llasso-post-sel-preparation.jl` inside `llasso-script.jl`.

So, we need to copy all the scripts to HGCC and the data:
```shell
ssh csolislemus@hgcc.genetics.emory.edu
mkdir 22q

cd Dropbox/Documents/gwas/projects/22q_new/llasso/scripts/
scp *.jl csolislemus@hgcc.genetics.emory.edu:/home/csolislemus/22q
```
We also need to copy the data:
```shell
cd Documents/gwas/data/22q/22q_files_NEW/bedfiles
ls
##22q-chr22.bed 22q-chr22.bim 22q-chr22.fam
scp * csolislemus@hgcc.genetics.emory.edu:/home/csolislemus/22q
```
Finally, we need to change the paths in the script to assume that everything is in the same folder.

Now we need to create the bash script for the job submission SGE: `llasso-submit.sh`. We need to copy it:
```
cd Dropbox/Documents/gwas/projects/22q_new/llasso/scripts/
scp *.jl csolislemus@hgcc.genetics.emory.edu:/home/csolislemus/22q
```




# Later when we do the interpretation
We might need the R package for Manhattan plot: `qqnan`, and the code that I did for Jai. The function needs BP to be in the order position in each chromosome. Info [here](http://www.gettinggeneticsdone.com/2014/05/qqman-r-package-for-qq-and-manhattan-plots-for-gwas-results.html )
```
dat = read.table("results/OptimalSortSkatWChrPos.txt",sep="\t")
dat2 = subset(dat,select=c("V2","V1","V4"))
names(dat2) = c("SNP","CHR","P")
str(dat2)

dat.ord = data.frame()

for(i in 1:max(dat2$CHR)){
  dat.chr = subset(dat2,CHR==i)
  dat.chr = within(dat.chr,BP<-order(dat.chr$SNP))
  str(dat.chr)
  dat.chr = dat.chr[order(dat.chr$BP),]
  dat.ord = rbind(dat.ord, dat.chr)
}
``