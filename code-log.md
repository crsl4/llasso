# To do
  - copy the llasso script to hgcc again, and modify path; copy covariates data too
  - run a submit script test again (make sure to delete all previous output files)
  - modify the scripts so that they take a number (the chromosome) as input, also copy the caucasian data (modify paths in scripts to read this new data); run julia script for all chromosomes: parallelize
  - after this, copy back the scripts used and results

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
This was fixed by updating R to 3.4.

Now, we need to copy all the scripts to HGCC and the data:
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
scp llasso-submit.sh csolislemus@hgcc.genetics.emory.edu:/home/csolislemus/22q
```

Now we can give it a try
```shell
ssh csolislemus@hgcc.genetics.emory.edu
qsub -q b.q -cwd -j y 22q/llasso-submit.sh
##Your job 240967 ("llasso-submit.sh") has been submitted
```
Note that we cannot copy and paste into the terminal, we need to type it, otherwise we get the error:
`Unable to read script file because of error: error opening –q: No such file or directory`

However, the job does not appear when doing `qstat`, and there are a bunch of errors:
```
/bin/mktemp: too many templates
Try '/bin/mktemp --help' for more information.
cp: target ‘/home/csolislemus/22q/llasso-submit.sh’ is not a directory
ERROR: could not open file /mnt/icebreaker/data2/home/csolislemus/llasso-script.jl
Stacktrace:
 [1] include_from_node1(::String) at ./loading.jl:569
 [2] include(::String) at ./sysimg.jl:14
 [3] process_options(::Base.JLOptions) at ./client.jl:305
 [4] _start() at ./client.jl:371
/bin/rm: cannot remove ‘/*.fam’: No such file or directory
/bin/rm: cannot remove ‘/*.bim’: No such file or directory
/bin/rm: cannot remove ‘/*.bed’: No such file or directory
rsync: link_stat "/mnt/icebreaker/data2/home/csolislemus/–av" failed: No such file or directory (2)
skipping directory .
rsync error: some files/attrs were not transferred (see previous errors) (code 23) at main.c(1052) [sender=3.0.9]
/bin/rm: cannot remove ‘–fr’: No such file or directory
```

The problem is that it cannot create the tmp folder.
It turns out that it was a hidden character in the submit file from copying from the slides. I had to rewrite the whole `llasso-submit.sh` and will copy again to HGCC:
```
cd Dropbox/Documents/gwas/projects/22q_new/llasso/scripts/
scp llasso-submit.sh csolislemus@hgcc.genetics.emory.edu:/home/csolislemus/22q
```

So, we can give it another try:
```shell
ssh csolislemus@hgcc.genetics.emory.edu
cd 22q
qsub -q b.q -cwd -j y llasso-submit.sh
##Your job 240976 ("llasso-submit.sh") has been submitted
```
It appears to be running now. Need to check back with `qstat`.
It ran without error, but no output files were produced. I will rerun without deleting the tmp folder, and then ask if it is possible to check if the files were produced in there.
```shell
hgcc:node00:[22q] % qsub -q b.q -cwd -j y llasso-submit.sh
Your job 240979 ("llasso-submit.sh") has been submitted
hgcc:node00:[22q] % qstat
job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID 
-----------------------------------------------------------------------------------------------------------------
 240979 0.09375 llasso-sub csolislemus  r     12/04/2017 13:15:13 b.q@node09.local                   1 
```
Again, we do not have output files, and there are none in the scratch tmp file:
```shell
qlogin -q i.q@node09
cd /scratch/
node09:[scratch] % ls -l
drwx------ 2 csolislemus    epsteinlab          4096 Dec  4 13:15 qfXh6c
cd qfXh6c/
 ```
So, let's do a small test `test.jl` that only creates a data frame"
```shell
qsub -q b.q -cwd -j y test.sh
```
Everything works, and the test.csv file is created!

So, I will run line by line the submit llasso-submit.sh, and see if there is any error. 
`qlogin -q i.q`: node08. I think that the problem is that I do not have the package JLD. I am installing it now, but I am unable to because I need sudo permissions to install it. Viren installed HDF5, so I will test again if I can install JLD.
```shell
ssh csolislemus@hgcc.genetics.emory.edu
qlogin -q i.q
module load julia
julia
```
Perfect! I was able to install and use `JLD` package. So now we will check the `llasso-submit.sh` script line by line. This will be run with the old 22q chromosome 22 data.

We got an error:
```
ERROR: LoadError: LoadError: Could not find R installation. Either set the "R_HOME" environmental variable, or ensure the R executable is available in "PATH".

ERROR: LoadError: LoadError: LoadError: Failed to precompile RCall to /home/csolislemus/.julia/lib/v0.6/RCall.ji.
Stacktrace:
 [1] compilecache(::String) at ./loading.jl:703
 [2] _require(::Symbol) at ./loading.jl:490
 [3] require(::Symbol) at ./loading.jl:398
 [4] include_from_node1(::String) at ./loading.jl:569
 [5] include(::String) at ./sysimg.jl:14
 [6] include_from_node1(::String) at ./loading.jl:569
 [7] include(::String) at ./sysimg.jl:14
 [8] include_from_node1(::String) at ./loading.jl:569
 [9] include(::String) at ./sysimg.jl:14
 [10] process_options(::Base.JLOptions) at ./client.jl:305
 [11] _start() at ./client.jl:371
while loading /scratch/PoR0DJ/llasso-post-sel-functions.jl, in expression starting on line 4
while loading /scratch/PoR0DJ/llasso-post-sel-preparation.jl, in expression starting on line 10
while loading /scratch/PoR0DJ/llasso-script.jl, in expression starting on line 144
```
Could this be a problem with running things in a screen? Maybe. Let's start the shell submit script and see what we get. First, we save whatever output from this first trial in folder `test1` in `output/22q`.
```shell
ssh csolislemus@hgcc.genetics.emory.edu
qsub -q b.q llasso-submit.sh
## Your job 241199 ("llasso-submit.sh") has been submitted
```
and I have the same RCall error!

So, I will modify the script to avoid the R part. Also, I will do the overfitting part (not great) of refitting a model only with the candidate SNPs (in lieu of the post selection). This is not great because we are using the same data, but we have so little data (~500) to begin with.