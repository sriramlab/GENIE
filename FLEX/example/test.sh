#!/bin/bash
# test.sh — End-to-end example pipeline for running FLEX, FLEX-cond-test, and FLEX-summ-h2

geno=gene                     # PLINK BED/BIM/FAM prefix for genotype data
cont_pheno=cont.pheno         # Continuous phenotype file
bin_pheno=bin.pheno           # Binary phenotype file
null_pheno=null.pheno         # Null phenotype file for testing

annot=${geno}.annot           # Annotation file for partitioned heritability
output=output.txt             # Output prefix

# --- 1. Run FLEX directly on genotype + phenotype ---
pheno=${null_pheno}
../build/FLEX_h2 -g ${geno} -p ${pheno} -k 10 -jn 5 -o ${output} -annot ${annot} -t 6

# --- 2. Run FLEX conditional test (FLEX-cond-test) ---
# The conditional test evaluates additional variance explained by a target set after conditioning.
output_dir=cond_test_output
python3 ../src/flex_cond_test.py -g ${geno} -p ${pheno} -o ${output_dir} -a ${annot}

# --- 3. Run FLEX-summ-h2 on GWAS summary statistics ---
output_dir=summ_h2_output
mkdir -p ${output_dir}

# 3a. Run PLINK2 to compute GWAS summary statistics for the phenotype.
plink2=/u/project/sriram/zhengton/softwares/plink2 ## update to your PLINK path
${plink2} --pheno ${pheno} --bfile ${geno} --glm 'allow-no-covars' --out ${output_dir}/gene --threads 6
mv  ${output_dir}/gene.pheno.glm.linear  ${output_dir}/gene.${pheno}.glm.linear
# Adjust header fields to expected format for FLEX-summ-h2
sed -i '1s/\tOBS_CT\t/\tN\t/; 1s/\tT_STAT\t/\tZ\t/' ${output_dir}/gene.${pheno}.glm.linear

# 3b. Generate LD scores and stochastic LD sketches from genotype + annotation.
python3 ../src/flex_summ_h2/main.py --geno ${geno} \
    --annot ${annot} \
    --out ${output_dir}/gene \
    --nvecs 100 \
    --nworkers 8

# 3c. Estimate gene-level heritability from GWAS summary statistics and sketches.
python3 ../src/flex_summ_h2/main.py \
        --pheno ${output_dir}/gene.${pheno}.glm.linear \
        --ldscores ${output_dir}/gene.ldscore.gz \
        --ld-sketch ${output_dir}/gene.sketches.npy \
        --annot ${annot} \
        --out ${output_dir}/output.summ.txt
