# ----------------------------------------------- #
# post selection inference of lasso coefficients
# ----------------------------------------------- #
## our julia functions are still under development, so we
## will use the R package selectiveInference
## this script has to be in julia because we have in julia
## the functions to read bedfile and convert.

## called inside llasso-script.jl, so no need to define variables
##include("llasso-functions.jl")
include("llasso-post-sel-functions.jl")
##datafolder = "/Users/Clauberry/Documents/gwas/data/22q/22q_files_NEW/bedfiles/"
##bedfile = string("22q_caucasian-chr",chrnum)

saveRdaFile(X,y,datafolder,bedfile)

