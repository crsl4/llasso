## Separate function to llasso-functions.jl because we do not
## want to have RCall in there, because of problems with HGCC

using RCall

## Function to read output files from llasso-script.jl
## and save an Rda file with the necessary objects for post selection
## the name of the Rda file will be bedfile.Rda
## WARNING: we cannot convert the bedfile inside because the convert function
## of SnpArrays.jl gives a different X every run
function saveRdaFile(X::AbstractArray, y::AbstractArray,datafolder::AbstractString, bedfile::AbstractString)
    @rput bedfile

    #chr = readBedfile(datafolder, bedfile);
    #X,y = convertBedfile(chr, datafolder,bedfile);
    @show size(X)
    @rput X
    @rput y


    ## Reading lambdas:--------------------------------------------------
    logfile = string(bedfile,".log")
    outfolder = ""
    lines = readlines(string(outfolder,logfile))
    who = String[]
    lambdas = Float64[]
    for l in lines
        if(contains(l,"Lasso.jl"))
            push!(who,"Lasso")
        elseif(contains(l,"SparseRegression.jl"))
            push!(who,"SR")
        end
        if(contains(l,"lambda ="))
            push!(lambdas,parse(Float64,string(split(l,"lambda = ")[2])))
        end
    end
    length(lambdas) <= 2 || error("Found more than 2 lambdas")

    obsl = [findfirst(who, "Lasso"),findfirst(who, "SR")]
    ## Lasso.jl lambda:
    resL = false
    if(obsl[1] != 0) ## found lambda for Lasso
        lambdaL = lambdas[obsl[1]]
        @rput lambdaL
        resL = true
    end

    ## SparseRegression.jl lambda:
    resSR = false
    if(obsl[2] != 0) ## found lambda for SR
        lambdaSR = lambdas[obsl[2]]
        @rput lambdaSR
        resSR = true
    end

    ## Beta coefficients:------------------------------------------------
    ##outfolder = "/Users/Clauberry/Dropbox/Documents/gwas/projects/22q_new/llasso/results/test1/"
    outfolder = ""
    outfile = string(bedfile,".beta")
    df = readtable(string(outfolder,outfile));
    obscol = [findfirst(names(df), :betaL),findfirst(names(df), :betaSR)]

    ## Lasso.jl lasso coefficients:
    if(obscol[1] != 0) ## found :betaL
        resL || warn("We did not find Lasso.jl lambda, but we found Lasso betas")
        beta_hatL = convert(Array,df[:betaL])
        @rput beta_hatL
        resL &= true
    end

    ## SparseRegression.jl lasso coefficients:
    if(obscol[2] != 0) ## found :betaSR
        resSR || warn("We did not find SR.jl lambda, but we found SR betas")
        beta_hatSR = convert(Array,df[:betaSR])
        @rput beta_hatSR
        resSR &= true
    end

    R"""
    save.image(file=paste0(bedfile,".Rda"))
    """
end
