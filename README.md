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

Alternatively, you may run either ```GENIE``` with a newline-separated config file:
```
./GENIE --config <config file>
```
When using a config file, the full keys (e.g., genotype="") must be used instead of shortcut flags (e.g., -g).

If there are many annotations, then we recommend to run a memory-efficient version of GENIE (that is slower). This can be enabled by setting
the flags -mem, --memeff on the command-line or memeff = 1 in the configuration file. This memory-efficient version reads in a specified block of SNPs that can be set 
using the flags -mem_Nsnp, --mem_Nsnp (mem_Nsnp in the configuration file). The default value for mem_Nsnp = 10. 
Alternatively, this memory efficient version can read an entire block of SNPs specified by the size of Jackknife blocks and this can be specified
using the flags -opt2 0, --opt2 0 (opt2 = 0 in the configuration file).



## Parameters

```
  -h, --help                                    Show this message and exit
  -g, --genotype                                The path of PLINK BED genotype file
  -p, --phenotype                               The path of phenotype file
  -c, --covariate                               The path of covariate file
  -a, --annot                                   The path of input annotation for partitioned heritability
  -o, --output                                  Output file path (prefix)
  -e, --environment                             The path of environment file
  -m, --model                                   Specification of the model. Currently there are 3 options:
                                                1. additive genetic (G) effects only (arg: G)
                                                        The model reduces to RHE-mc (Pazokitoroudi et al. Nat Commun (2020). https://doi.org/10.1038/s41467-020-17576-9).
                                                2. additive genetic (G) and gene-environment (GxE) effects (arg: G+GxE)
                                                        The model treats noise/environment effects as homogeneous.
                                                3. additive genetic (G), gene-environment (GxE) and heterogeneous noise (NxE) effects (arg: G+GxE+NxE)
                                                        The model treats noise/environment effects as heterogeneous.
  -k, --num-vec                                 The number of random vectors (10 is recommended).
  -jn, --num-jack                               The number of jackknife blocks (100 is recommended).
  -t, --nthreads                                The number of threads for multithreading
  -tr, --trace                                  Flag for printing trace summary files (.trace).
  -tr_input                                     Read in the trace estimates.
  -v, --verbose                                 Verbose modes; Output extra information (Normal equation, number of samples, etc.). Default = 0. 
                                                    Setting '--verbose V' where V is 1,2,3,4, or 5 prints out increasingly more information during the run.
  -eXa, --eXannot                               By default, GENIE fits a single GxE variance component. To partition the GxE component w.r.t the annotation file, add '-eXannot' flag.
  -np, --norm-proj-pheno                        By default, the phenotype vector is standardized after regressing covariates. Turn this off by setting '--norm-proj-pheno 0'.
  -i, --cov-add-intercept                       By default, a vector of ones is appended to the covariates (the intercept term). Turn this off by setting '--cov-add-intercept 0'.

```
## File formats
```
Genotype: must be in PLINK BED format.
Phenotype: must have a header in the following format 
(multiple phenotypes only supported by GENIE_multi_pheno): 
    FID IID name_of_pheno_1 name_of_pheno_2  . . .   name_of_pheno_n
Covariate: must have a header in the following format: 
    FID IID name_of_cov_1 name_of_cov_2  . . .   name_of_cov_n
Environment: must have a header in the following format: 
    FID IID name_of_env_1 name_of_env_2  . . .   name_of_env_n
Annotation: must have M rows (M=number  of SNPs) and K columns (K=number of annotations).
    If SNP i belongs to annotation j, then there is  "1" in row i and column j.
    Otherwise, there is "0". (delimiter is " ")

1) The number and order of individuals must be the same in phenotype, 
    genotype, environment, and covariate files.
2) The number and order of SNPs must be the same in bim file and annotation file.
3) The annotation file does not have a header. 
4) SNPs with MAF=0 must be excluded from the genotype file.
5) GENIE excludes individuals with NA or -9 values in the phenotype or 
    environment files from the analysis.
```

## GxE linear mixed models
```
GENIE is able to fit single/multiple additive and GxE variance components. 
1) SNPs are partitioned with respect to the annotation file.
2) By default, GENIE fits a single GxE variance component. 
    To partition the GxE component w.r.t the annotation file, please add "-eXannot" flag.
3) GENIE fits noise by environment (NxE) variance component (heterogeneous noise).
4) By default, GENIE adds environmental variables to covariates as fixed effects.
```

## Toy example 
Sample files are provided in the example directory to ensure that everything works well. 
Look at test.sh file and run it  :
```
chmod +x test.sh
./test.sh
```

## Simulator
To simulate phenotypes with GxE effects, see https://github.com/sriramlab/Simulator.

## Citation
```
1. Ali Pazokitoroudi, Zhengtong Liu, Andrew Dahl, Noah Zaitlen, Saharon Rosset, Sriram Sankararaman.
   AJHG (2024); doi: [10.1016/j.ajhg.2024.05.015](https://doi.org/10.1016/j.ajhg.2024.05.015)

2. Ali Pazokitoroudi, Yue Wu, Kathryn S. Burch, Kangcheng Hou, Aaron Zhou, Bogdan Pasaniuc, Sriram Sankararaman.
   Nature Communications (2020); doi: [10.1038/s41467-020-17576-9](https://doi.org/10.1038/s41467-020-17576-9)
```

## Version
```
v1.2.0
```

## Zenodo
https://zenodo.org/records/14618542
