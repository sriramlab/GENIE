# FLEX
**F**ast **L**D-aware **E**stimation of e**X**ome-wide and gene-level heritability

## Prerequisites
The following packages are required on a Linux machine to compile and use FLEX:
```
g++
cmake
make
python3
```


## Installation

```
git clone https://github.com/sriramlab/GENIE.git
cd GENIE/FLEX
mkdir build
cd build/
cmake ..
make
```

# Usage
Running FLEX-h2
 ```
 bash build/FLEX_h2 <command-line arguments>
 ```

Running FLEX-cond-test
```
 python3 src/flex_cond_test.py <command-line arguments>
```

## Parameters

For FLEX_h2

```
  -h, --help                Show this help message and exit
  -b, --binary              Flag indicating whether the phenotype is binary
  -g, --genotype            Path to PLINK BED genotype file
  -p, --phenotype           Path to phenotype file
  -c, --covariate           Path to covariate file
  -a, --annot               Path to annotation file for partitioned heritability
  -o, --output              Output file path (prefix)
  -k, --num-vec             Number of random vectors (recommended: 10)
  -jn, --num-jack           Number of jackknife blocks (recommended: 100)
  -t, --nthreads            Number of threads for multithreading
```

## File formats

Genotype
    * Must be in PLINK BED format.
    * SNPs with MAF = 0 must be excluded.
Phenotype
    * Must include a header in the format: 
    ```
        FID IID name_of_pheno_1 name_of_pheno_2  . . .   name_of_pheno_n
    ```
    * Individuals with ``NA`` or ``-9`` values are excluded.
Covariate
    * Must include a header in the format:
    ``` 
        FID IID name_of_cov_1 name_of_cov_2  . . .   name_of_cov_n
    ```
Annotation
    * Matrix with `M` rows (number of SNPs) and `K` columns (number of annotations).
    * Entry = ``1`` if SNP i belongs to annotation j, otherwise ``0``.
    * Delimiter: space ``" "``.
    * No header row.
    * The **first annotation** is treated as the surrounding LD region to be conditioned on.


Important notes
1) The number and order of individuals must be the same in phenotype, 
    genotype, environment, and covariate files.
2) The number and order of SNPs must be the same in ``.bim`` file and annotation file.



## Example 
A toy example is provided in the ``example/`` directory. To run it: 
```
chmod +x example/test.sh
bash example/test.sh
```

## Simulator
To simulate phenotypes with GxE effects, see https://github.com/sriramlab/Simulator.

## Citation
If you use FLEX in your work, please cite: 
```

```

## Version
```
v1.0.0
```