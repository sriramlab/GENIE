#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os, sys
import numpy as np
import pandas as pd
from bed_reader import open_bed

# sys.path.append('/u/home/z/zhengton/project-sriram/wes/scripts/gene_h2/')
# sys.path.append('/u/project/sriram/zhengton/wes/scripts/gene_h2/')
from regional_her import RegionalHer

from chi2comb import chi2comb_cdf, ChiSquared
import time
import scipy

import argparse

def score_test(sigma, S, dofs):
    k = len(S)
    ncents = np.zeros(k)
    chi2s = [ChiSquared(S[i], ncents[i], dofs[i]) for i in range(k)]
    t0 = time.time()
    p, error, info = chi2comb_cdf(sigma, chi2s, 0, lim= 10000000, atol=1e-30)
    # p = qf.qf(0, Phi, acc = 1e-7)[0] 
    t1 = time.time()
    return (1-p,error)


parser = argparse.ArgumentParser()
parser.add_argument("-g", "--geno_prefix", help="Genotype file prefix.")
parser.add_argument("-a", "--annot_file", help="Annotation file.")
parser.add_argument("-p", "--pheno_file", help="Phenotype file.")
parser.add_argument("-o", "--output_dir", help = "output directory.")
parser.add_argument("-c", "--covar", help='Covariance file.')
parser.add_argument("-b", "--binary", action='store_true', help='If binary trait or not.')
args = parser.parse_args()

# chrom = sys.argv[1]
# gene = sys.argv[2]
# data_folder = "/u/scratch/z/zhengton/gene_h2/"
# covar_file="/u/project/sriram/ukbiobank/data/wes/pheno/stratified_PCs.bin1_5.covar"
# geno_file=f"{data_folder}/geno/chr_10kb_mac_1/chr{chrom}/{gene}.new"
# pheno_dir=f"{data_folder}/pheno/continuous/"
# outdir=f"{data_folder}/results/gene_test/chr{chrom}/{gene}/"

geno_file = args.geno_prefix
annot_file = args.annot_file
pheno_file = args.pheno_file
covar_file = args.covar
outdir = args.output_dir
binary_trait = args.binary


# In[4]:
print(f"current geno: {geno_file}")

annot_mat = np.loadtxt(f"{annot_file}")


regional_her = RegionalHer(f"{geno_file}.bed", pheno_file, covar_file=covar_file, multi_pheno=False, binary_trait=binary_trait)


h2_dict, norm_A, norm_b = regional_her.estimate_partitioned_regional_her2(annot_mat, num_vecs=10, debug_flag=True)
print(h2_dict)

results = pd.DataFrame(h2_dict).T
results.columns = ['common', 'low_freq', 'rare', 'ultra-rare', 'noise']
results['all'] = results.iloc[:, :4].sum(axis=1)


P = np.linalg.inv(norm_A)
v = np.array([1, 1, 1, 1, 0])
N = regional_her.G.shape[0]


annot_mat = annot_mat[:, 1:]
num_annot = annot_mat.shape[1]
eigenvalues = []
dofs = []
L = 200
for k in np.arange(num_annot+1):
    print(f'start annotation {k}')
    if k == num_annot:
#         eigens = np.ones(L+1)
#         dof = np.ones(L)
#         dof = np.append(dof, N-L)
        eigens = np.array([1])
        dof = np.array([N])
    else:
        indices = np.where(annot_mat[:, k] == 1)[0]
        geno = regional_her.G.read(index=np.s_[:, indices])
        means = np.nanmean(geno, axis=0)
        stds = np.nanstd(geno, axis=0)
        ## standardization
        geno = (geno - means) / stds
        ## mean imputation
        geno[np.isnan(geno)] = 0
        print('start eigen')
        K = (geno.T @ geno) / geno.shape[1]
        K2 = K @ K
        tr = np.trace(K)
        tr2 = np.trace(K @ K.T)
        eigens = scipy.linalg.eig(K, left=False, right=False).real
        eigens[eigens <= 1e-6] = 0
        eigens = eigens[np.nonzero(eigens)]
        eigens = eigens[~np.isnan(eigens)]
        eigens = eigens[:L]
        if (tr - np.sum(eigens) > 0) and (tr2 - np.sum(np.square(eigens)) > 0):
            a = (tr2 - np.sum(np.square(eigens)))/(tr - np.sum(eigens))
            d = np.square((tr - np.sum(eigens)))/(tr2 - np.sum(np.square(eigens)))
            dof = np.append(np.ones(len(eigens)), d)
            eigens = np.append(eigens, a)
        else:
            dof = np.ones(len(eigens))
    dofs.append(dof)
    eigenvalues.append(eigens)


# eigenvalues = np.concatenate(eigenvalues).flatten()
# print(eigenvalues)
weights = v.T @ P
# print(weights)
weights = weights.flatten()
weighted_eigens = ([weights[k]*eigenvalues[k] for k in range(num_annot+1)])
coef_list = np.concatenate(weighted_eigens).flatten()
dof_list = np.concatenate(dofs).flatten()
print(len(coef_list), len(dof_list))

pvalues = []
for h2gene in results['all'].values:
    pvalue = score_test(h2gene, coef_list, dof_list)[0]
    pvalues.append(pvalue)
results['pvalue'] = pvalues

os.makedirs(outdir, exist_ok=True)
np.savetxt(f"{outdir}/coefs.txt", coef_list)
np.savetxt(f"{outdir}/dofs.txt", dof_list)
results.index.name = 'pheno'
results = results.reset_index()
results.to_csv(f"{outdir}/all_continuous_2.tsv", sep='\t', index=False)




