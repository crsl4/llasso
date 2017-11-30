using SnpArrays, DataFrames, SparseRegression, Distributions, Lasso, IHT

## function that reads a bedfile, removes rare variants (by default) and check the missingness
## We keep a list of the SNPs that are still in the model
## (to be able to give interpretations at the end)
function readBedfile(datafolder::AbstractString, bedfile::AbstractString; 
    eliminateRare=true::Bool)
    chr = SnpArray(string(datafolder,bedfile))
    maf, minor_allele, missings_by_snp, missings_by_person = summarize(chr)
    if(eliminateRare)
        info("Eliminating rare variants (MAF<5%)")
        ind = find(x->x>=0.05,maf)
        writetable(string(bedfile,".common"),DataFrame(ind=ind))
        chr = chr[:, maf .>= 0.05]
    else
        warn("We are not eliminating rare variants (MAF<5%), and this can cause problems in the post selection functions")
    end 
    miss = sum(missings_by_snp) / length(chr)
    info("Percentage of missingness in matrix: $(round(miss*100,2))%")
    miss > 0.15 && warn("you have missingness of more than 15%")
    return chr
end 

## function that reads the fam file, and extracts the 6th column (case-control)
## and returns it as an Array{Float64}
## assumes that the name of the file is: datafile/bedfile.fam
function readResponse(datafolder::AbstractString, bedfile::AbstractString)
    dat = readtable(string(datafolder,bedfile,".fam"), separator=' ', header=false)
    warn("Reading fam file: assuming that the 6th column is the case-control status")
    y = convert(Array{Float64,1},dat[:,6]) ## assumes 6th column is case-control status
    y = y-1 ##0=controls, 1=cases
    return y
end 

## function that fits Lasso, but instead of returning an error, it returns a false if failed
function fitLLasso(X::AbstractArray, y::AbstractArray;itol = 0.0000001::Float64)
    try
        f = fit(LassoPath,X,y,Bernoulli(),LogitLink(),irls_tol=itol)
        return true,f
    catch
        return false,nothing
    end 
end

## fit Lasso function for different tolerance values until it converges
## default numtry=20 for largest tolerance = 0.1
function fitRepeatedLlasso(X::AbstractArray, y::AbstractArray;
    itol = 0.0000001::Float64, numtry = 20::Int)
    res,f = fitLLasso(X,y)
    i = 1
    while(!res && i<numtry)
        info("Lasso.fit did not converge for irls_tol of $itol, will try now for $(itol*2)")
        itol *= 2
        res,f = fitLLasso(X,y, itol=itol)
        i += 1
    end
    !res && warn("Lasso fit did not converge after $numtry attempts, with relaxed tolerance of $itol")
    return res,f
end

## function to get the data split for cross validation
## we always take a sample of snps because we need to fit multiple times the model
function getSubsetXY4CV(chr::SnpArray, datafolder::AbstractString, bedfile::AbstractString,
    prc::Number, rng::Int, nsnps::Int)

    nsub = Int(floor(size(chr,1)*prc))
    srand(rng)
    ## Sample SNPs if there are more than 10,000
    indsnp = size(chr,2) > nsnps ? sample(1:size(chr,2),nsnps,replace=false) : collect(1:1:size(chr,2))
    ## Sample the individuals to fit the model
    ind = sample(1:size(chr,1),nsub,replace=false)
    chrsub = chr[ind,indsnp]
    Xsub = convert(Matrix{Float64},chrsub,impute=true) ## need to impute for SparseRegression   
    ## Sample the individuals to predict the response
    nall = collect(1:1:size(chr,1))
    other = setdiff(nall,sort(ind))
    chrother = chr[other,indsnp]
    Xother = convert(Matrix{Float64},chrother,impute=true)  
    ## Read in the response
    y = readResponse(datafolder, bedfile)   
    ysub = y[ind]
    yother = y[other]
    return Xsub, Xother, ysub, yother
end

## function to estimate y from an estimated beta
## instead of cutoff, we are sampling bernoulli with the estimated p
## returns score = ||y-yest||_1
function predictY(X::AbstractArray,y::AbstractArray,beta::AbstractArray)
    size(X,1) == size(y,1) || error("Dimension mismatch: X,y")
    size(X,2) == size(beta,1) || error("Dimension mismatch: X, beta")
    a = exp.(X*beta)
    p1 = a./(1+a)
    yest = rand.(Bernoulli.(p1))
    score = sum(abs.(y-yest))
    return score
end

## function to choose lambda with Lasso.jl
## the fit function already fits for multiple lambdas
## - chr: SnpArray read from bedfile, and removed rare variants already
## - bedfile (assumes that the fam file is of the form bedfile.fam)
## - prc = percentage of the data split to validate a given lambda (0.5, default)
## - rng = seed for random numbers (12345, default)
## - nsnps = number of SNPs to sample for speed (10000, default)
## returns the df with scores, best lambda, lasso fit
function chooseLambdaLasso(chr::SnpArray, datafolder::AbstractString, 
    bedfile::AbstractString; 
    prc=0.5::Number, rng=12345::Int, nsnps=10000::Int)
    
    Xsub, Xother, ysub, yother = getSubsetXY4CV(chr, datafolder, bedfile, prc, rng, nsnps)
    res,f = fitRepeatedLlasso(Xsub, ysub)
    @show f
    res || error("Lasso fit never converged, cannot get best lambda")

    df = DataFrame(lambda=Number[], score=Number[], numcov=Number[])
    for i in 1:length(f.λ)
        l = f.λ[i]
        @show l
        beta = f.coefs[:,i]
        score = predictY(Xother,yother,beta)
        ncov = length(beta[beta .!= 0])
        push!(df,[l,score,ncov])
    end
    ## find the lambda with smaller score:
    bestmod = find(x->x==minimum(df[:score]),df[:score])
    if(length(bestmod)>1)
        bestmod = bestmod[end]
    end
    bestl = df[:lambda][bestmod]
    println("Best model: $(df[bestmod,:])")
    return df,bestl,f
end

#datafolder = "/Users/Clauberry/Documents/gwas/data/22q/22q_files_NEW/bedfiles/"
#bedfile = "22q-chr22"
#lambdas = [0.0007,0.0008,0.0009]
#prc = 0.5
#rng = 12345
#nsnps = 10000
#chooseLambda(datafolder,bedfile,lambdas)

## function to do pseudo cross validation on SparseRegression for lambda
## - lambdas: a given vector of lambda values
## - chr: SnpArray read from bedfile, and removed rare variants already
## - bedfile (assumes that the fam file is of the form bedfile.fam)
## - prc = percentage of the data splitted to validate a given lambda (0.5, default)
## - rng = seed for random numbers (12345, default)
## - nsnps = number of SNPs to sample for speed (10000, default)
## returns the best lambda
function chooseLambdaSparseRegression(chr::SnpArray, datafolder::AbstractString, 
    bedfile::AbstractString, lambdas::AbstractArray; 
    prc=0.5::Number, rng=12345::Int, nsnps=10000::Int)

    Xsub, Xother, ysub, yother = getSubsetXY4CV(chr, datafolder, bedfile, prc, rng, nsnps)
    Xsub = hcat(fill(1.0,size(Xsub,1)),Xsub) ## we need to add the intercept, because SparseReg does not do this
    Xother = hcat(fill(1.0,size(Xother,1)),Xother) ## we need to add the intercept, because SparseReg does not do this
    ysub[ysub.==0] = -1 ##SparseRegression needs response -1,1
    yother[yother.==0] = -1 ##SparseRegression needs response -1,1

    df = DataFrame(lambda=Number[], score=Number[], numcov=Number[])
    for l in lambdas
        ## fitting lasso on prc% of data:
        info("Fitting llasso for lambda = $l")
        s = SModel(Xsub,ysub, LogitDistLoss(), L1Penalty(),fill(l,size(Xsub,2)))
        try
            learn!(s, strategy(ProxGrad(s,0.001), MaxIter(10000), Converged(coef)));
        catch
            warn("SparseRegression fit did not converge for lambda = $l, will skip it")
            continue
        end       
        ## prediting the (1-prc)% of data
        beta = s.β
        score = predictY(Xother,yother,beta)
        cov = length(beta[beta .!= 0])
        push!(df,[l,score,cov])
    end
    ## find the lambda with smaller score:
    bestmod = find(x->x==minimum(df[:score]),df[:score])
    if(length(bestmod)>1)
        bestmod = bestmod[end]
    end
    bestl = df[:lambda][bestmod]
    println("Best model: $(df[bestmod,:])")
    return df,bestl
end

## This function identifies the repeated columns in X
## It writes a file with column kept index, 
## followed by the identical column indices that should be excluded
## It also returns a vector with the indices to exclude
function identifyRepeatedColumns(X::AbstractArray, filename::String)
    y = sum(X,1)
    excluded = Int[]
    text = "col-kept col-excluded\n"
    for i in 1:length(y)
        ind = find(x->x==y[i], y[i+1:end])
        isempty(ind) && continue
        firsttime = true
        written = false
        for j in 1:length(ind)
            if( all(X[:,i] .== X[:,ind[j]+i]) )
                if(firsttime)
                    text *= string(i," ")
                    firsttime = false
                end
                push!(excluded,ind[j]+i)
                text *= string(ind[j]+i," ")
                written = true
            end
        end
        text *= written ? "\n" : ""
    end
    f = open(filename,"w")
    write(f,text)
    close(f)
    return unique(excluded)
end

function removeRepeatedColumns(X::AbstractArray; filename="indices.txt"::String)
    excluded = identifyRepeatedColumns(X,filename)
    kept = setdiff(1:size(X,2),excluded)
    writetable(string(filename,"-kept"),DataFrame(ind=kept))
    return X[:,kept]
end

## We need a function that will convert the SnpArray into a Float matrix
## fixit: do we need Float or can Int work?
## next fixit: we want to change the functions to avoid this step
function convertBedfile(chr::SnpArray, datafolder::AbstractString,bedfile::AbstractString)
    size(chr,2) > 100000 && warn("SnpArray has more than 100k SNPs (after removing rare variants). Conversion to float matrix will be heavy")
    if(size(chr,2) > 200000)
        warn("SnpArray has more than 200k SNPs (after removing rare variants) so we will take only the first 100k for this run")
        chr = chr[:,1:100000]
    end
    y = readResponse(datafolder, bedfile)
    X = @time convert(Matrix{Float64},chr,impute=true)
    ## Note that R post selection functions do not allow repeated columns,
    ## we could get rid of the repeated columns with DataFrame, but we want to
    ## keep track of which columns we are eliminating
    X = removeRepeatedColumns(X, filename=string(bedfile,".indices"))
    return X,y
end




