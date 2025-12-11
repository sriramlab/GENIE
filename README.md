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

## CLI Quick Start

### Model G: Additive Heritability Only
```bash
./GENIE -g data/mydata -p data/mydata.pheno -e data/mydata.env -a data/mydata.annot -o results/out -m G -k 10 -jn 100
```

### Model G+GxE: Gene-Environment Interaction
```bash
./GENIE -g data/mydata -p data/mydata.pheno -e data/mydata.env -a data/mydata.annot -o results/out -m G+GxE -k 10 -jn 100
```

### Model G+GxE+NxE: Full Model
```bash
./GENIE -g data/mydata -p data/mydata.pheno -e data/mydata.env -a data/mydata.annot -o results/out -m "G+GxE+NxE" -k 10 -jn 100
```

### With Covariates
```bash
./GENIE -g data/mydata -p data/mydata.pheno -e data/mydata.env -c data/mydata.cov -a data/mydata.annot -o results/out -m "G+GxE+NxE" -k 10 -jn 100
```

See below for usage with config files.

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

## Config File Examples

### Basic Model G (Additive genetic effects only)
```
genotype=data/mydata
phenotype=data/mydata.pheno
covariate=data/mydata.cov
annotation=data/mydata.annot
output=results/model_g
model=G
num_vec=10
num_jack=100
nthreads=4
seed=42
```

### Full GxE Model (G+GxE+NxE) - Recommended for GxE analysis
```
genotype=data/mydata
phenotype=data/mydata.pheno
covariate=data/mydata.cov
environment=data/mydata.env
annotation=data/mydata.annot
output=results/model_gxe
model=G+GxE+NxE
num_vec=10
num_jack=100
nthreads=4
seed=42
```

### Memory-Efficient Mode (for large datasets, Model G only)
```
genotype=data/mydata
phenotype=data/mydata.pheno
annotation=data/mydata.annot
output=results/model_g_memeff
model=G
num_vec=10
num_jack=100
nthreads=4
memory_mode=1
mem_Nsnp=10
seed=42
```

### Config File Parameters Reference
| Config Key | CLI Flag | Description |
|------------|----------|-------------|
| genotype | -g | Path to PLINK BED file (without .bed extension) |
| phenotype | -p | Path to phenotype file |
| covariate | -c | Path to covariate file |
| environment | -e | Path to environment file |
| annotation | -a | Path to annotation file |
| output | -o | Output file prefix |
| model | -m | Model type: G, G+GxE, or G+GxE+NxE |
| num_vec | -k | Number of random vectors (recommended: 10) |
| num_jack | -jn | Number of jackknife blocks (recommended: 100) |
| nthreads | -t | Number of threads |
| seed | -s | Random seed for reproducibility |
| memory_mode | -mm | Memory mode: 0 (standard), 1 (memory-efficient), 2 (most memory-efficient) |
| memeff | -- | Legacy: Memory mode 1 (use memory_mode=1 instead) |
| opt1 | -- | Legacy: Memory mode 2 (use memory_mode=2 instead) |
| opt2 | -- | Legacy: Used with opt1 for mode 2 (use memory_mode=2 instead) |
| mem_Nsnp | -- | SNPs per block in memory mode 1 (default: 10) |
| trace | -tr | Print trace summary files (0 or 1) |
| verbose | -v | Verbosity level (0-5) |
| eXannot | -eXa | Partition GxE by annotation (0 or 1) |
| norm_proj_pheno | -np | Normalize projected phenotype (default: 1) |
| cov_add_intercept | -i | Add intercept to covariates (default: 1) |
| no_match_ids | --no-match-ids | Skip sample ID matching (faster, requires pre-aligned files) |

## Memory Modes

GENIE supports three memory modes for handling large datasets. Use `-mm`/`--memory-mode` on the command line or `memory_mode` in config files:

| Mode | CLI | Config Setting | Description |
|------|-----|----------------|-------------|
| 0 (default) | `-mm 0` | `memory_mode=0` | Standard mode. Fastest but uses most memory. |
| 1 | `-mm 1` | `memory_mode=1` | Memory-efficient. Reads SNPs in blocks of size `mem_Nsnp` (default: 10). |
| 2 | `-mm 2` | `memory_mode=2` | Most memory-efficient. Reads SNPs in jackknife-sized blocks. For very large datasets. |

**Note**: Memory modes 1 and 2 are only supported for Model G. GxE models (`G+GxE`, `G+GxE+NxE`) require the default memory mode (`memory_mode=0`).

## Sample Matching

By default, GENIE matches samples across genotype, phenotype, covariate, and environment files by FID/IID and analyzes only the intersection.

For pre-aligned files where row order is guaranteed to match across all input files, you can use `--no-match-ids` (CLI) or `no_match_ids=1` (config) to skip ID matching. This is faster but will produce incorrect results if files are not perfectly aligned.

## Parameters

```
  -h, --help                                    Show this message and exit
  -V, --version                                 Show version and exit
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
  -s, --seed                                    Set the random seed for reproducibility. Using the same seed produces identical results.
  -mm, --memory-mode                            Memory mode: 0 (standard), 1 (memory-efficient), 2 (most memory-efficient).
                                                Modes 1 and 2 only supported for Model G.
  --no-match-ids                                Disable sample ID matching (use faster legacy position-based matching).
                                                WARNING: May produce incorrect results if files are not pre-aligned.

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

1) Samples are matched by FID/IID across all input files; only the intersection is analyzed.
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
v1.2.2
```

## Zenodo
https://zenodo.org/records/14618542
