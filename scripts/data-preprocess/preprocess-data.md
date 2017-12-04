# Convert to more intuitive genotype data
The `bed` file has the genotype data, but it is not intuitive to read. I think that the ped file is a bit more intuitive. We need a more intuitive file if we want to use julia to make plots and stuff.
## MAP files: The fields in a MAP file are:
  - Chromosome
  - Marker ID
  - Genetic distance
  - Physical position
Example of a MAP file of the standard PLINK format:
```
21	rs11511647	0	26765
X	rs3883674	0	32380
X	rs12218882	0	48172
9	rs10904045	0	48426
9	rs10751931	0	49949
8	rs11252127	0	52087
10	rs12775203	0	52277
8	rs12255619	0	52481
```

## PED files: the fields in a PED file are
  - Family ID
  - Sample ID
  - Paternal ID
  - Maternal ID
  - Sex (1=male; 2=female; other=unknown)
  - Affection (0=unknown; 1=unaffected; 2=affected)
  - Genotypes (space or tab separated, 2 for each marker. 0=missing)
Example of a PED file of the standard PLINK format:
```
FAM1	NA06985	0	0	1	1	A	T	T	T	G	G	C	C	A	T	T	T	G	G	C	C
FAM1	NA06991	0	0	1	1	C	T	T	T	G	G	C	C	C	T	T	T	G	G	C	C
0	NA06993	0	0	1	1	C	T	T	T	G	G	C	T	C	T	T	T	G	G	C	T
0	NA06994	0	0	1	1	C	T	T	T	G	G	C	C	C	T	T	T	G	G	C	C
0	NA07000	0	0	2	1	C	T	T	T	G	G	C	T	C	T	T	T	G	G	C	T
0	NA07019	0	0	1	1	C	T	T	T	G	G	C	C	C	T	T	T	G	G	C	C
0	NA07022	0	0	2	1	C	T	T	T	G	G	0	0	C	T	T	T	G	G	0	0
0	NA07029	0	0	1	1	C	T	T	T	G	G	C	C	C	T	T	T	G	G	C	C
FAM2	NA07056	0	0	0	2	C	T	T	T	A	G	C	T	C	T	T	T	A	G	C	T
FAM2	NA07345	0	0	1	1	C	T	T	T	G	G	C	C	C	T	T	T	G	G	C	C
```
## Additive format
Using the --recodeAD option generates the file plink-recode.raw:
```
     FID IID PAT MAT SEX PHENOTYPE snp1_2 snp1_HET snp2_G snp2_HET
     1 1 0 0 1 1  0  0   2 0
     1 2 0 0 2 1  NA NA  1 1
     1 3 0 0 1 1  0  0   1 1
     1 4 0 0 2 1  1  1   0 0
```
The column labels reflect the snp name (e.g. snp1) with the name of the minor allele appended (i.e. snp1_2 in the first instance, as 2 is the minor allele) for the additive component. The dominant component ( a dummy variable reflecting heterozygote state) is coded with the _HET suffix.
This file can be easily loaded into R: for example: `d <- read.table("plink.raw",header=T)`

## Converting data:
1. First we extract that chromosome:
```shell
cd Documents/gwas/data/22q/22q_files_NEW/bedfiles
plink --bfile ../justfinalset/22q_all_dropped_519_pruned --chr 1 --make-bed --out 22q-chr1
PLINK v1.90b4.4 64-bit (21 May 2017)           www.cog-genomics.org/plink/1.9/
(C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to 22q-chr1.log.
Options in effect:
  --bfile ../justfinalset/22q_all_dropped_519_pruned
  --chr 1
  --make-bed
  --out 22q-chr1

8192 MB RAM detected; reserving 4096 MB for main workspace.
2394339 out of 30747425 variants loaded from .bim file.
519 people (229 males, 290 females) loaded from .fam.
519 phenotype values loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 519 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.996703.
2394339 variants and 519 people pass filters and QC.
Among remaining phenotypes, 259 are cases and 260 are controls.
--make-bed to 22q-chr1.bed + 22q-chr1.bim + 22q-chr1.fam ... done.
```
2. We convert to PED:
```shell
plink --bfile 22q-chr1 --recode tab --out 22q-chr1-ped
```
3. We convert to additive format that could be read in R:
```shell
plink --file 22q-chr1-ped --recodeAD --out 22q-chr1-add
PLINK v1.90b4.4 64-bit (21 May 2017)           www.cog-genomics.org/plink/1.9/
(C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
Note: --recodeAD flag deprecated.  Use 'recode AD ...'.
Logging to 22q-chr1-add.log.
Options in effect:
  --file 22q-chr1-ped
  --out 22q-chr1-add
  --recode AD

8192 MB RAM detected; reserving 4096 MB for main workspace.
Possibly irregular .ped line.  Restarting scan, assuming multichar alleles.
.ped scan complete (for binary autoconversion).
Performing single-pass .bed write (2394339 variants, 519 people).
--file: 22q-chr1-add-temporary.bed + 22q-chr1-add-temporary.bim +
22q-chr1-add-temporary.fam written.
2394339 variants loaded from .bim file.
519 people (229 males, 290 females) loaded from .fam.
519 phenotype values loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 519 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.996703.
2394339 variants and 519 people pass filters and QC.
Among remaining phenotypes, 259 are cases and 260 are controls.
--recode AD to 22q-chr1-add.raw ... done.
```
We did these steps manually for chromosome 1 only, but we use the julia script `convert2ped.jl` to convert the data for all chromosomes (except 24-26).
It is actually more efficient to do `plink --bfile 22q-chr1 --recodeAD --out 22q-chr1-add`, so we don't really need the PED files.
Got these warnings:
```
Warning: 53716 het. haploid genotypes present (see 22q-chr23.hh ); many
commands treat these as missing.
Warning: 53716 het. haploid genotypes present (see 22q-chr23-ped.hh ); many
commands treat these as missing.
Warning: 53716 het. haploid genotypes present (see 22q-chr23-add.hh ); many
commands treat these as missing.
```

I created `logfiles` folder for the log files (and hh files), and `bedfiles` folder for the `bed,fam,bim` files. Also, folder `pedfiles` for the `ped,map` files;  `rawfiles` folder for `raw` files. The PED files are not needed and very heavy (~70Gb), so I will delete them.

Now, we want to read these files into julia. It turns out that the [OpenMendel](https://github.com/OpenMendel) package, in particular [SnpArrays.jl](https://github.com/OpenMendel/SnpArrays.jl) has a function `SnpArray` that read the BED file directly into julia:
```julia
using SnpArrays
chr1 = SnpArray("bedfiles/22q-chr1");
```
We get this output:
```
519×2394339 SnpArrays.SnpArray{2}:
 (true, true)  (true, true)  (true, true)   (true, true)   …  (true, true)  (true, true)  (true, true)   (true, true) 
 (true, true)  (true, true)  (true, true)   (true, true)      (true, true)  (true, true)  (true, true)   (true, true) 
 (true, true)  (true, true)  (true, true)   (true, true)      (true, true)  (true, true)  (true, true)   (true, true) 
 (true, true)  (true, true)  (true, true)   (true, false)     (true, true)  (true, true)  (true, true)   (true, true) 
 (true, true)  (true, true)  (true, true)   (true, true)      (true, true)  (true, true)  (true, true)   (true, true) 
 (true, true)  (true, true)  (true, true)   (true, true)   …  (true, true)  (true, true)  (true, true)   (true, true) 
 (true, true)  (true, true)  (true, true)   (true, true)      (true, true)  (true, true)  (true, true)   (true, true) 
 (true, true)  (true, true)  (true, true)   (true, true)      (true, true)  (true, true)  (true, true)   (true, true) 
 (true, true)  (true, true)  (true, true)   (true, true)      (true, true)  (true, true)  (true, true)   (true, true) 
 (true, true)  (true, true)  (true, true)   (true, true)      (true, true)  (true, true)  (true, true)   (true, true) 
 (true, true)  (true, true)  (true, true)   (true, true)   …  (true, true)  (true, true)  (true, true)   (true, true) 
 (true, true)  (true, true)  (true, true)   (true, true)      (true, true)  (true, true)  (true, false)  (true, false)
 (true, true)  (true, true)  (true, true)   (true, true)      (true, true)  (true, true)  (true, true)   (true, true) 
 ```
 which has the right dimensions, and represents (A1,A1)=(false,false), (A1,A2)=(false,true), (A2,A2)=(true,true), missing value=(true,false). 
 Also, we want to convert to numerical matrix with minor allele count: 0,1,2. First, it finds which one is the minor allele (because at this point we only denote A1,A2):
 ```julia
 chr1mat=convert(Matrix{Float64},chr1);
 ```
 We want to compare this matrix to the one we get with `plink`:
 ```julia
 using DataFrames
 chr1plink = readtable("rawfiles/22q-chr1-add.raw", separator=' ', header=true)
 ```
We cannot compare the two matrices, because we cannot convert the DataFrame into a matrix because of the NAs.

Since we can import data into Julia directly from the bed files, we don't need the raw files anymore (and they are heavy ~70Gb). I will delete the rawfiles as well.
Now, we have the input dataset. We need the model. See `llasso-model.jmd`.

# Using sample of europeans only
We will test llasso method with a sample of europeans only (because we do not want to include the PCs as population structure, actually we cannot include any covariates because of the penalization).

We need the first two columns of `22q_caucasians_outliers_removed.cov` which has the IDs of the caucasian individuals.
```shell
cd Dropbox/Documents/gwas/projects/22q_new/caucasian-ind
cut -d' ' -f1-2 22q_caucasian_outliers_removed.cov > 22q_caucasian_id.txt
cp 22q_caucasian_id.txt ~/Documents/gwas/data/22q/22q_files_NEW/bedfiles
```

Now that we have the IDs to keep, we run:
```shell
cd Documents/gwas/data/22q/22q_files_NEW/
plink --bfile justfinalset/22q_all_dropped_519_pruned --keep 22q_caucasian_id.txt --make-bed --out 22q_caucasian

PLINK v1.90b4.4 64-bit (21 May 2017)           www.cog-genomics.org/plink/1.9/
(C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to 22q_caucasian.log.
Options in effect:
  --bfile justfinalset/22q_all_dropped_519_pruned
  --keep 22q_caucasian_id.txt
  --make-bed
  --out 22q_caucasian

8192 MB RAM detected; reserving 4096 MB for main workspace.
30747425 variants loaded from .bim file.
519 people (229 males, 290 females) loaded from .fam.
519 phenotype values loaded from .fam.
--keep: 435 people remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 435 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: 48792 het. haploid genotypes present (see 22q_caucasian.hh ); many
commands treat these as missing.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate in remaining samples is 0.996878.
30747425 variants and 435 people pass filters and QC.
Among remaining phenotypes, 214 are cases and 221 are controls.
--make-bed to 22q_caucasian.bed + 22q_caucasian.bim + 22q_caucasian.fam ...
done.

mv 22q_caucasian* bedfiles/
```

Now, we want to extract one chromosome at a time (again):
```shell
cd Documents/gwas/data/22q/22q_files_NEW/bedfiles
plink --bfile 22q_caucasian --chr 1 --make-bed --out 22q_caucasian-chr1

PLINK v1.90b4.4 64-bit (21 May 2017)           www.cog-genomics.org/plink/1.9/
(C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to 22q_caucasian-chr1.log.
Options in effect:
  --bfile 22q_caucasian
  --chr 1
  --make-bed
  --out 22q_caucasian-chr1

8192 MB RAM detected; reserving 4096 MB for main workspace.
2394339 out of 30747425 variants loaded from .bim file.
435 people (196 males, 239 females) loaded from .fam.
435 phenotype values loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 435 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.996494.
2394339 variants and 435 people pass filters and QC.
Among remaining phenotypes, 214 are cases and 221 are controls.
--make-bed to 22q_caucasian-chr1.bed + 22q_caucasian-chr1.bim +
22q_caucasian-chr1.fam ... done.
```
We do this for every chromosome, and then move the logfiles to log folder.

We have the data ready now for the llasso script!