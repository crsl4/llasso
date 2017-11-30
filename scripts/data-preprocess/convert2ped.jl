## julia script that will convert chromosomes BED files into PED files.
## we need to be in folder: Documents/gwas/projects/22q/22q_files_NEW/llasso
## Claudia September 2017

## chromosome 1 already done manually
for i in 2:23
    println("chromosome $i")
    try
        out = readstring(`plink --bfile ../justfinalset/22q_all_dropped_519_pruned --chr $i --make-bed --out 22q-chr$i`);
    catch
        warn("Error in plink when extracting chromosome $i")
    end
    try
        out = readstring(`plink --bfile 22q-chr$i --recode tab --out 22q-chr$(i)-ped`);
    catch
        warn("Error in plink when converting to PED for chromosome $i")
    end
    try
        out = readstring(`plink --bfile 22q-chr$i --recodeAD --out 22q-chr$(i)-add`);
    catch
        warn("Error in plink when converting to additive format for chromosome $i")
    end
end


