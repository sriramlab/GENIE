# GENIE
**G**ene-**EN**viroment **I**nteraction **E**stimator



## Prerequisites
The following packages are required on a Linux machine to compile and use the software package.
```
g++
cmake
make
```

## How to install :

```
git clone https://github.com/sriramlab/GENIE.git
cd GENIE
mkdir build
cd build/
cmake ..
make
```

# Documentation for GENIE
An executable file named GENIE will be in the build folder after the installation steps. Run GENIE as follows:
 ```
 ./GENIE <command_line arguments>
```
If there are many annotations, then we recommend to run the following : 
```
 ./GENIE_mem <command_line arguments>
```
GENIE_mem is slower than GENIE, but it uses less memory than GENIE.


## Parameters

```
Genotype (-g): The path of  genotype file
phenotype (-p): The path of phenotype file
covariate (-c): The path of covariate file
environment (-e): The path of environment file
annotation (-annot): The path of genotype annotation file.
num_vec (-k): The number of random vectors (10 is recommended). 
num_block (-jn): The number of jackknife blocks (100 is recommended). The higher the number of jackknife blocks, the higher the memory usage.
output (-o): The path of the output file.
model (-m): Specification of the model: it reduces to [RHE-mc](https://www.nature.com/articles/s41467-020-17576-9) if the model only fits the additive genetic (G) component; users can also consider to estimate GxE heritability with or without the noise heterogeneous component (NxE) (G/G+GxE/G+GxE+NxE). 
num_threads (-t): The number of threads.
seed (-s): The random seed.
verbose (-v): Whether to output extra information or not (0/1).

By default, GENIE fits a single GxE variance component. To partition the GxE component w.r.t the annotation file, add "-eXannot" flag. The phenotype vector is standardized after regressing covariates. To turn this off, add "-norm_proj_pheno 0". In addition, a one's vector is appended to the covariates (the intercept term). To remove this intercept term, add "-cov_add_intercept 0".

```
## File formats
```
Genotype: It must be in bed format.
Phenotype: It must have a header in the following format: FID IID name_of_phenotype
Covariate: It must have a header in the following format: FID IID name_of_cov_1 name_of_cov_2  . . .   name_of_cov_n
Environment: It must have a header in the following format: FID IID name_of_env_1 name_of_env_2  . . .   name_of_env_n
Annotation: It has M rows (M=number  of SNPs) and K columns (K=number of annotations). If SNP i belongs to annotation j, then there is  "1" in row i and column j. Otherwise, there is "0". (delimiter is " ")

1) The number and order of individuals must be the same in phenotype, genotype, environment, and covariate files.
2) The number and order of SNPs must be the same in bim file and annotation file.
3) The annotation file does not have a header. 
4) SNPs with MAF=0 must be excluded from the genotype file.
5) GENIE excludes individuals with NA or -9 values in the phenotype or environment files from the analysis.
```

## GxE linear mixed models
```
GENIE is able to fit single/multiple additive and GxE variance components. 
1) SNPs are partitioned with respect to the annotation file.
2) By default, GENIE fits a single GxE variance component. To partition the GxE component w.r.t the annotation file, please add "-eXannot" flag.
3) GENIE fits noise by environment (NxE) variance component (heterogeneous noise).
4) By default, GENIE adds environmental variables to covariates as fixed effects.
```

## Toy example 
Sample files are provided in the example directory to ensure that everything works well. Look at test.sh file and run it  :
```
chmod +x test.sh
./test.sh
```

## Citation
```
```


