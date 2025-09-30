# FLEX
**F**ast **L**D-aware **E**stimation of e**X**ome-wide and gene-level heritability

FLEX provides fast, LD-aware estimation of gene-level and exome-wide heritability from either **individual-level genotypes** or **GWAS summary statistics**. It also includes a conditional testing module to perform a calibrated and sensitive test for gene-level heritability.

---

## Prerequisites

Install these packages on a Linux machine to build and use FLEX:
```bash
g++
cmake
make
python3
```

You will also need **PLINK 2** in your `PATH` to generate GWAS summary statistics in the example pipeline.

---

## Installation

FLEX consists of two parts:
1. **FLEX_h2 (C++ core binary)** for individual-level heritability estimation.  
2. **Python scripts** (`flex_cond_test.py` and `flex_summ_h2/main.py`) for conditional testing and summary-statistics mode.  

### 1. FLEX_h2 (C++ core binary)

Clone and build:

```bash
git clone https://github.com/sriramlab/GENIE.git
cd GENIE/FLEX
mkdir build && cd build
cmake ..
make
```

This produces the main binary (named `FLEX_h2`) under `build/`.

---

### 2. Python Dependencies (for FLEX-cond-test and FLEX-summ-h2)

If you plan to run the Python scripts, you need the following Python packages:

```text
numpy>=1.20
pandas>=2.0
scipy>=1.10
tqdm>=4.66
bed-reader<=1.0
matplotlib>=3.7
scikit-learn>=1.3
statsmodels>=0.14
chi2comb>=0.1
```

#### Option A: Using Conda with environment.yml

We provide an `environment.yml` for convenience:

```bash
conda env create -f environment.yml
conda activate flex
```

#### Option B: Using pip and requirements.txt

If you prefer a virtual environment:

```bash
python -m venv flex
# On Linux/macOS
source flex/bin/activate
# On Windows
flex\Scripts\activate

pip install -r requirements.txt
```
---

## Command-line Interfaces

### 1) FLEX (individual-level heritability; a.k.a. `FLEX_h2`)

```text
  -h, --help                 Show this help message and exit
  -b, --binary               Treat phenotype as binary (default: off)
  -g, --genotype             Prefix to PLINK BED/BIM/FAM genotype files
  -p, --phenotype            Path to phenotype file
  -c, --covariate            Path to covariate file
  -a, --annot                Path to SNP annotation file for partitioned heritability
  -o, --output               Output file prefix
  -k, --num-vec              Number of random vectors for stochastic trace (recommend: 10)
  -jn, --num-jack            Number of jackknife blocks (recommend: 100)
  -t, --nthreads             Number of threads
```

### 2) FLEX-cond-test (conditional testing; Python)

```text
  -g, --geno_prefix          Genotype file prefix (PLINK)
  -a, --annot_file           Annotation file
  -p, --pheno_file           Phenotype file
  -o, --output_dir           Output directory
  -c, --covar                Covariate file
  -b, --binary               Treat phenotype as binary (flag)
```

### 3) FLEX-summ-h2 (summary-statistics mode; Python)

FLEX-summ-h2 operates in **two stages**.

**(i) LD scores & stochastic LD sketches from genotypes**  
(Use this to prepare LD inputs once per reference panel / annotation.)

```text
  --geno                     Prefix to PLINK BED/BIM files (reference panel)
  --annot                    SNP annotation file
  --out                      Output prefix for LD scores and sketches
  --nworkers                 Number of workers for multiprocessing (default: 4)
  --nvecs                    Number of random vectors for LD sketches (default: 10)
  --step_size                SNPs per block during LD computation (default: 1000)
  --seed                     Random seed (optional)
```

**(ii) Heritability from GWAS summary statistics**  
(Consumes (i) outputs along with the GWAS z-scores.)

```text
  --pheno                    GWAS summary statistics path (or a directory containing *.sumstat)
  --annot                    SNP annotation file (must match the LD inputs)
  --ldscores                 Path to LD scores (.ldscore.gz)
  --ld-sketch                Path to stochastic LD sketches (.npy)
  --out                      Output prefix (writes <out>.log with all console prints)
  --filter-both-sides        Apply SNP filtering to both T and q if needed (requires LD scores)
  --verbose                  Verbose logging
```

---

## File Formats

**Genotype (PLINK)**  
- Input genotypes must be in PLINK **BED/BIM/FAM** format.  
- SNPs with **MAF = 0** should be removed.

**Phenotype**  
- Header format:
  ```text
  FID IID phenotype_1 phenotype_2 ... phenotype_n
  ```
- `NA` or `-9` values are excluded.  

**Covariates**  
- Header format:
  ```text
  FID IID cov_1 cov_2 ... cov_k
  ```

**Annotation**  
- Matrix with `M` rows (SNPs) and `K` columns (annotations).  
- Entry is `1` if SNP *i* belongs to annotation *j*, else `0`.  
- Space-delimited.
- The **first annotation** is used as flanking/LD region to condition on (see manuscript).

> **Ordering:** The number and order of **individuals** must be consistent across phenotype, genotype, and covariate files.  
> The number and order of **SNPs** must match between `.bim` and the annotation file.

---

## End-to-end Example Pipeline (`example/test.sh`)

The example script demonstrates the usage of three algorithms:

**1) FLEX on individual-level data**  
Estimates exome-wide and gene-level heritability (optionally partitioned) from genotypes.

**2) FLEX-cond-test**  
Performs a conditional test to quantify additional variance explained by a target annotation after conditioning on others.

**3) FLEX-summ-h2 in two stages**  
(i) **Generate LD scores and stochastic LD sketches** from genotype and annotation.  
(ii) **Estimate gene-level heritability** (and SEs) from **GWAS summary statistics** and the outputs of (i).

> Ultra-rare variants have already been collapsed upstream. Ensure your input genotype/annotation reflect this preprocessing.

Below is a minimal version of the commands used in `example/test.sh` (adjust paths as needed).

```bash
#!/usr/bin/env bash

# --- Paths & names (toy example) ---
geno_prefix="gene"            # PLINK prefix: gene.bed/.bim/.fam
annot_file="gene.annot"       # SNP x annotation matrix
pheno_cont="cont.pheno"       # continuous phenotype
pheno_bin="bin.pheno"         # binary phenotype (optional)
pheno_null="null.pheno"       # null phenotype (example)
out_prefix="output"           # shared output prefix for FLEX
threads=6

# 1) Individual-level heritability with FLEX
build/FLEX_h2   -g "${geno_prefix}"   -p "${pheno_null}"   -a "${annot_file}"   -o "${out_prefix}"   -k 10   -jn 5   -t "${threads}"

# 2) Conditional test (Python)
python3 src/flex_cond_test.py   -g "${geno_prefix}"   -p "${pheno_null}"   -a "${annot_file}"   -o "cond_test_output"   -c ""   -b    # drop -b if trait is continuous

# 3) FLEX-summ-h2 in two stages
outdir="summ_h2_output"
mkdir -p "${outdir}"

plink2 --bfile "${geno_prefix}"        --pheno "${pheno_null}"        --glm 'allow-no-covars'        --threads "${threads}"        --out "${outdir}/gene"

mv  "${outdir}/gene.pheno.glm.linear"  "${outdir}/gene.${pheno_null}.glm.linear"
sed -i '1s/\tOBS_CT\t/\tN\t/; 1s/\tT_STAT\t/\tZ\t/' "${outdir}/gene.${pheno_null}.glm.linear"

python3 src/flex_summ_h2/main.py   --geno "${geno_prefix}"   --annot "${annot_file}"   --out "${outdir}/gene"   --nvecs 100   --nworkers 8   --step_size 1000

python3 src/flex_summ_h2/main.py   --pheno     "${outdir}/gene.${pheno_null}.glm.linear"   --annot     "${annot_file}"   --ldscores  "${outdir}/gene.ldscore.gz"   --ld-sketch "${outdir}/gene.sketches.npy"   --out       "${outdir}/output.summ"
```

---

## Notes & Tips

- **Annotation**: The first column often denotes flanking LD regions to be conditioned on; match your study design.  
- **Ultra-rare variant collapsing**: Ensure collapsing is performed **before** running FLEX, consistent with the manuscript’s specification.

---

## Citation

If you use FLEX, please cite the accompanying manuscript.

---

## Version

```
v1.0.0
```
