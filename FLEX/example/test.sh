#!/bin/bash

data_folder=/u/scratch/z/zhengton/gene_h2/

geno=gene
cont_pheno=cont.pheno
bin_pheno=bin.pheno
null_pheno=null.pheno

annot=${geno}.annot
output=output.txt

# ../build/FLEX_h2 -g ${geno} -p ${null_pheno} -k 10 -jn 5 -o ${output} -annot ${annot} -t 6
# ../build/FLEX_h2 -g ${geno} -p ${cont_pheno} -k 10 -jn 5 -o ${output} -annot ${annot} -t 6
# ../build/FLEX_h2 -b -g ${geno} -p ${bin_pheno} -k 10 -jn 5 -o ${output} -annot ${annot} -t 6

output_dir=cond_test_output
python3 ../src/flex_cond_test.py -g ${geno} -p ${cont_pheno} -o ${output_dir} -a ${annot}