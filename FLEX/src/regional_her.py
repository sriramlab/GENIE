import numpy as np
import pandas as pd
from bed_reader import open_bed
from sklearn.preprocessing import StandardScaler
import multiprocessing as mp
from multiprocessing import Process, Manager
from tqdm import tqdm
from functools import partial
import scipy.stats as stats 
from scipy.stats import norm
import matplotlib.pyplot as plt
from time import time
import scipy
from glob import glob
from functools import reduce
import statsmodels.api as sm
import os
import warnings
from scipy.linalg import pinvh
from sklearn.decomposition import PCA

class RegionalHer:
    def __init__(self, bed_file, pheno_file, covar_file, pre_compute_Xty=False, save_folder=None, Xty_file=None,
                multi_pheno=False, pheno_file_list=None, binary_trait=False):
        self.G = open_bed(bed_file)
        self.N, self.M = self.G.shape
        self.num_vecs, self.num_workers, self.step_size = 10, 8, 1000
        self.region_indices = np.arange(self.M)
        self.partitioned_indices = None
        self.multi_pheno = multi_pheno
        self.binary_trait = binary_trait
        print(f"The number of individuals: {self.N}")
        print(f"The number of SNPs: {self.M}\n")
        
        if covar_file is not None:
            covar_df = pd.read_csv(covar_file, sep=' ')
            if 'ethnic' in covar_df.columns:
                covar_df = covar_df.drop("ethnic", axis=1)
            
        scaler = StandardScaler()
        if multi_pheno:
            if self.binary_trait:
                self.prevalence_dict = {}
            if pheno_file_list is not None:
                full_pheno_path_list = pheno_file_list
            else:
                pheno_dir_path = pheno_file
                full_pheno_path_list = glob(f"{pheno_dir_path}/*.pheno")
            pheno_list = ['.'.join(file.split('/')[-1].split('.')[:-1]) for file in full_pheno_path_list]
            pheno_name_to_file = {'.'.join(file.split('/')[-1].split('.')[:-1]) : file for file in full_pheno_path_list}
            self.pheno_names = pheno_list
            pheno_list = self.pheno_names
            self.num_pheno = len(pheno_list)
            pheno_df_list = []
            if covar_file is not None:
                all_covar = scaler.fit_transform(covar_df.iloc[:, 2:].values)
            for pheno_name, pheno_file in tqdm(pheno_name_to_file.items(), total=self.num_pheno):
                if covar_file is not None:
                    covar = all_covar
                pheno_df = pd.read_csv(f"{pheno_file}", sep='\s+')
                pheno_df.columns = ['FID', 'IID', pheno_name]
                pheno = pheno_df[pheno_name].values.flatten()
                if self.binary_trait:
                    notnan_pheno = pheno[pheno != -9]
                    # print(f"prevalence for {pheno_name} is {np.nanmean(notnan_pheno)}")
                    self.prevalence_dict[pheno_name] = np.nanmean(notnan_pheno)
                pheno[np.isnan(pheno)] = -9
                if covar_file is not None:
                    covar = covar[pheno != -9, :]
                    not_na_pheno = pheno[pheno != -9]
                    proj = covar @ (np.linalg.inv(covar.T @ covar) @ (covar.T @ not_na_pheno))
                    res = pheno[pheno != -9] - proj
                    pheno_df.loc[pheno_df[pheno_name] != -9, pheno_name] = res
                ind_idxs = pheno_df[pheno_df[pheno_name] != -9].index.values
                pheno = pheno_df.loc[ind_idxs, pheno_name].values.flatten()
                pheno = (pheno - np.nanmean(pheno))/np.nanstd(pheno)
                pheno_df.loc[ind_idxs, pheno_name] = pheno
                pheno_df.loc[pheno_df[pheno_name] == -9, pheno_name] = 0

                pheno_df_list.append(pheno_df)

            warnings.filterwarnings('ignore')
            # pheno_df_merged = reduce(lambda left, right: pd.merge(left, right, on=['FID', 'IID']), pheno_df_list)
            pheno_df_merged = pheno_df[['FID', 'IID']]
            for k, pheno_name in enumerate(pheno_list):
                pheno_df = pheno_df_list[k]
                pheno_df_merged[pheno_name] = pheno_df[pheno_name]

            self.pheno = pheno_df_merged.iloc[:, 2:].values
            self.ind_idxs = np.arange(len(pheno_df_merged))
        else:
            if binary_trait:
                self.prevalence_dict = {}
            self.num_pheno = 1
            pheno_df = pd.read_csv(pheno_file, sep='\s+')
            pheno_name = pheno_df.columns[-1]
            if covar_file is not None:
                covar = scaler.fit_transform(covar_df.iloc[:, 2:].values)
                pheno_name = pheno_df.columns[-1]
                pheno = pheno_df[pheno_name].values.flatten()
                if binary_trait:
                    notnan_pheno = pheno[pheno != -9]
                    print(notnan_pheno[:10])
                    self.prevalence_dict[pheno_name] = np.nanmean(notnan_pheno)
                pheno[np.isnan(pheno)] = -9
                covar = covar[pheno != -9, :]
                not_na_pheno = pheno[pheno != -9]
                proj = covar @ (np.linalg.inv(covar.T @ covar) @ (covar.T @ not_na_pheno))
                res = pheno[pheno != -9] - proj
                pheno_df.loc[pheno_df[pheno_name] != -9, pheno_name] = res
            else:
                pheno_name = pheno_df.columns[-1]
                pheno = pheno_df[pheno_name].values.flatten()
                if binary_trait:
                    notnan_pheno = pheno[pheno != -9]
                    print(notnan_pheno[:10])
                    self.prevalence_dict[pheno_name] = np.nanmean(notnan_pheno)
            trait_pre_df = pd.read_csv(f"/u/scratch/z/zhengton/gene_h2/annotation/bin_disease_phenos_prevalence.txt",
                sep='\t')
            self.prevalence_dict = trait_pre_df[['trait', 'k']].set_index('trait').to_dict()['k']
            ind_idxs = pheno_df[pheno_df[pheno_name] != -9].index.values
            pheno = pheno_df.loc[ind_idxs, pheno_name].values.flatten()
            pheno = (pheno - np.nanmean(pheno))/np.nanstd(pheno)
            pheno_df.loc[ind_idxs, pheno_name] = pheno
            pheno_df.loc[pheno_df[pheno_name] == -9, pheno_name] = 0

            self.pheno_name = pheno_name
            self.pheno = pheno_df.iloc[:, 2:].values
            # self.ind_idxs = ind_idxs
            self.ind_idxs = np.arange(len(pheno_df))
        
        self.total_Xty_results = None
        if pre_compute_Xty:
            with open(f"{save_folder}/{Xty_file}", 'rb') as f:
                self.total_Xty_results = np.load(f)
            print(f"Loaded pre-computed Xty...")
            
        
    def _compute_Xz(self, current_input):
        j, current_indices = current_input
        num_snps = len(current_indices)

        Zs = np.random.normal(size=(num_snps, self.num_vecs))
        
        geno = self.G.read(index=np.s_[:, current_indices])
        

        
        means = np.nanmean(geno, axis=0)
        stds = np.nanstd(geno, axis=0)

        ## standardization
        geno = (geno - means) / stds
        ## mean imputation
        geno[np.isnan(geno)] = 0


        geno = np.array(geno, order='F')
        Zs = np.array(Zs, order='F')
        cur_Xz = scipy.linalg.blas.sgemm(1.0, geno, Zs)

        
        return cur_Xz

    def _compute_XtXz(self, current_input):
        j, current_indices = current_input
        num_snps = len(current_indices)

        geno = self.G.read(index=np.s_[:, current_indices])
        means = np.nanmean(geno, axis=0)
        stds = np.nanstd(geno, axis=0)

        ## standardization
        geno = (geno - means) / stds
        ## mean imputation
        geno[np.isnan(geno)] = 0

        geno_t = np.array(geno.T, order='F')
        cur_XtXz = scipy.linalg.blas.sgemm(1.0, geno_t, self.temp_results)
        return (j, cur_XtXz)

    def _compute_Xty(self, current_input):
        j, current_indices = current_input
        num_snps = len(current_indices)
        geno = self.G.read(index=np.s_[self.ind_idxs, current_indices])

        means = np.nanmean(geno, axis=0)
        stds = np.nanstd(geno, axis=0)

        ## standardization
        geno = (geno - means) / stds
        ## mean imputation
        geno[np.isnan(geno)] = 0
        geno_t = np.array(geno.T, order='F')
        if self.num_pheno > 100:
            cur_Xty = np.zeros((geno.shape[1], self.num_pheno))
            all_indices = np.arange(self.num_pheno)[::100]
            for k in tqdm(range(len(all_indices)), total=len(all_indices), desc='Partitioned calculations of Xty'):
                if k == len(all_indices)-1:
                    temp_pheno = np.array(self.pheno[:, all_indices[k]:], order='F')
                    cur_Xty[:, all_indices[k]:] = \
                        scipy.linalg.blas.sgemm(1.0, geno_t, temp_pheno)
                else:
                    temp_pheno = np.array(self.pheno[:, all_indices[k]:all_indices[k+1]], order='F')
                    cur_Xty[:, all_indices[k]:all_indices[k+1]] = \
                        scipy.linalg.blas.sgemm(1.0, geno_t, temp_pheno)
        else:
            pheno = np.array(self.pheno, order='F')
            cur_Xty = scipy.linalg.blas.sgemm(1.0, geno_t, pheno)
        return (j, cur_Xty)
    
    def _compute_X1(self, current_input):
        j, current_indices = current_input
        num_snps = len(current_indices)

        Zs = np.ones((num_snps, 1))
        
        geno = self.G.read(index=np.s_[:, current_indices])
        

        
        means = np.nanmean(geno, axis=0)
        stds = np.nanstd(geno, axis=0)

        ## standardization
        geno = (geno - means) / stds
        ## mean imputation
        geno[np.isnan(geno)] = 0


        geno = np.array(geno, order='F')
        Zs = np.array(Zs, order='F')
        cur_Xz = scipy.linalg.blas.sgemm(1.0, geno, Zs)

        
        return cur_Xz
    

    def compute_ld_proj(self, region_indices=None, num_vecs=None, num_workers=None, step_size=None):
        if region_indices is not None: self.region_indices = region_indices
        if num_vecs: self.num_vecs = num_vecs
        if num_workers: self.num_workers = num_workers
        if step_size: self.step_size = step_size
            
        M_region = len(self.region_indices)
        
        print(f"...Start computing LD projections, with the set of parameters: ")
        print(f"\t num_vecs={self.num_vecs}, num_workers={self.num_workers}, step_size={self.step_size}")
        print(f"\t number of SNPs in this region: {M_region}")
        
        
        all_indices = np.arange(len(self.region_indices))[::self.step_size]
        num_blocks = len(all_indices)
        self.partitioned_indices = []
        for j in range(num_blocks):
            if j == (num_blocks - 1):
                self.partitioned_indices.append((j, self.region_indices[all_indices[j]:]))
            else:
                self.partitioned_indices.append((j, self.region_indices[all_indices[j]:all_indices[j+1]]))


        self.temp_results = np.zeros((self.N, self.num_vecs))
        results = np.zeros((M_region, self.num_vecs))

        for cur_indices in tqdm(self.partitioned_indices, desc='Calculating Xz', total=len(self.partitioned_indices)):
            cur_res = self._compute_Xz(cur_indices)
            self.temp_results += cur_res

        results_dict = {}

        for cur_indices in tqdm(self.partitioned_indices, desc='Calculating XtXz', total=len(self.partitioned_indices)):
            label, cur_res = self._compute_XtXz(cur_indices)
            results_dict[label] = cur_res

        
        for j in range(num_blocks):
            if j == (num_blocks-1):
                start_index = all_indices[j]
                results[start_index:, :] = results_dict[j]
            else:
                start_index, end_index = all_indices[j], all_indices[j+1]
                results[start_index:end_index, :] = results_dict[j]
        return results
    
    def compute_partitioned_ld_proj(self, partitioned_region_indices, region_indices=None, num_vecs=None, 
                                    num_workers=None, step_size=None):
        if region_indices is not None: self.region_indices = region_indices
        if num_vecs: self.num_vecs = num_vecs
        if num_workers: self.num_workers = num_workers
        if step_size: self.step_size = step_size
            
        M_region = len(self.region_indices)
        
        print(f"...Start computing LD projections, with the set of parameters: ")
        print(f"\t num_vecs={self.num_vecs}, num_workers={self.num_workers}, step_size={self.step_size}")
        print(f"\t number of SNPs in this region: {M_region}")
        
        
        all_indices = np.arange(len(self.region_indices))[::self.step_size]
        num_blocks = len(all_indices)
        self.partitioned_indices = []
        for j in range(num_blocks):
            if j == (num_blocks - 1):
                self.partitioned_indices.append((j, self.region_indices[all_indices[j]:]))
            else:
                self.partitioned_indices.append((j, self.region_indices[all_indices[j]:all_indices[j+1]]))

        all_results = []
        for cur_regional_indices in partitioned_region_indices:
            self.temp_results = np.zeros((self.N, self.num_vecs))
            results = np.zeros((M_region, self.num_vecs))
        
            all_partitioned_indices = np.arange(len(cur_regional_indices))[::self.step_size]
            num_blocks = len(all_partitioned_indices)
            cur_partitioned_indices = []
            for j in range(num_blocks):
                if j == (num_blocks - 1):
                    cur_partitioned_indices.append((j, cur_regional_indices[all_partitioned_indices[j]:]))
                else:
                    cur_partitioned_indices.append((j, 
                                cur_regional_indices[all_partitioned_indices[j]:all_partitioned_indices[j+1]]))
        
            for cur_indices in tqdm(cur_partitioned_indices, desc='Calculating Xz', 
                                    total=len(cur_partitioned_indices)):
                cur_res = self._compute_Xz(cur_indices)
                self.temp_results += cur_res

            results_dict = {}

            for cur_indices in tqdm(self.partitioned_indices, desc='Calculating XtXz', total=len(self.partitioned_indices)):
                label, cur_res = self._compute_XtXz(cur_indices)
                results_dict[label] = cur_res

            num_blocks = len(all_indices)
            for j in range(num_blocks):
                if j == (num_blocks-1):
                    start_index = all_indices[j]
                    results[start_index:, :] = results_dict[j]
                else:
                    start_index, end_index = all_indices[j], all_indices[j+1]
                    results[start_index:end_index, :] = results_dict[j]
            all_results.append(results)
        return all_results
    
    
    
    def compute_Xty(self, region_indices=None, num_workers=None, step_size=None):
        if region_indices is not None: self.region_indices = region_indices
        if num_workers: self.num_workers = num_workers
        if step_size: self.step_size = step_size
        
        M_region = len(self.region_indices)
        
        Xty_results = np.zeros((M_region, self.num_pheno))
        
        # if self.total_Xty_results is not None:
        #     Xty_results = self.total_Xty_results[self.region_indices]
        #     return Xty_results
        
        print(f"...Start computing Xty, with the set of parameters: ")
        print(f"\t num_workers={self.num_workers}, step_size={self.step_size}")
        print(f"\t number of SNPs in this region: {M_region}")
        

        all_indices = np.arange(len(self.region_indices))[::self.step_size]
        num_blocks = len(all_indices)
        self.partitioned_indices = []
        for j in range(num_blocks):
            if j == (num_blocks - 1):
                self.partitioned_indices.append((j, self.region_indices[all_indices[j]:]))
            else:
                self.partitioned_indices.append((j, self.region_indices[all_indices[j]:all_indices[j+1]]))
        
        Xty_results_dict = {}
        # for cur_indices in tqdm(self.partitioned_indices, desc='Calculating Xty', total=len(self.partitioned_indices)):
        for cur_indices in self.partitioned_indices:
            label, cur_res = self._compute_Xty(cur_indices)
            Xty_results_dict[label] = cur_res

        for j in range(num_blocks):
            if j == (num_blocks-1):
                start_index = all_indices[j]
                Xty_results[start_index:] = Xty_results_dict[j]
            else:
                start_index, end_index = all_indices[j], all_indices[j+1]
                Xty_results[start_index:end_index] = Xty_results_dict[j]
        return Xty_results
    
    def estimate_regional_her(self, region_indices=None, num_vecs=None, num_workers=None, step_size=None,
                             jackknife_blocks=None):
        if region_indices is not None: self.region_indices = region_indices
        if num_vecs: self.num_vecs = num_vecs
        if num_workers: self.num_workers = num_workers
        if step_size: self.step_size = step_size
        
        B = self.num_vecs
        M_region = len(self.region_indices)
        
        region_ld_proj = self.compute_ld_proj()
    
        if self.multi_pheno:
            all_region_Xty = self.compute_Xty()
            tr_KK = (region_ld_proj.T @ region_ld_proj).sum()/(B*M_region*M_region)

            b_trK = self.N
            Nc = self.N
            
            norm_eq_A = np.array([[tr_KK, b_trK],
                      [b_trK,   Nc]])
            
            h2_vec_dict = {}
            for k in tqdm(range(self.num_pheno), total=self.num_pheno, desc='Estimating h2 for multiple phenotypes'):
                cur_pheno = self.pheno[:, k]
                region_Xty = all_region_Xty[:, k]
                cur_pheno_name = self.pheno_names[k]

                c_yKy = np.square(region_Xty).sum()/M_region
                yy = np.inner(cur_pheno, cur_pheno)
                norm_eq_b = np.array([c_yKy, yy])


                sigma_vec = np.linalg.solve(norm_eq_A, norm_eq_b)

                h2_vec = sigma_vec/sum(sigma_vec)
                h2_vec_dict[cur_pheno_name] = sigma_vec


            if jackknife_blocks != None:
                print(f"Start calculating jackknife subsample estimates, num of jack: {jackknife_blocks}")
                sigma_vec_jack_list = []
                num_snps = M_region // jackknife_blocks
                h2_vec_std_dict = {}
                for j in tqdm(range(jackknife_blocks), total=jackknife_blocks):
                    if j == (jackknife_blocks - 1):
                        num_snps = (M_region // jackknife_blocks) + (M_region % jackknife_blocks)
                    start_index, end_index = j*(M_region // jackknife_blocks), j*(M_region // jackknife_blocks)+num_snps

                    region_ld_proj_jack = np.delete(region_ld_proj, np.arange(start_index, end_index), 0)
                    all_region_Xty_jack = np.delete(all_region_Xty, np.arange(start_index, end_index), 0)

                    M_region_jack = M_region - num_snps
                    tr_KK_jack = (region_ld_proj_jack.T @ region_ld_proj_jack).sum()/(B*M_region_jack*M_region_jack)

                    norm_eq_A = np.array([[tr_KK_jack, b_trK],
                                        [b_trK,   Nc]])
                    
                    for k in tqdm(range(self.num_pheno), total=self.num_pheno, desc='Estimating h2 for multiple phenotypes'):
                        cur_pheno = self.pheno[:, k]
                        region_Xty_jack = all_region_Xty_jack[:, k]
                        cur_pheno_name = self.pheno_names[k]

                        c_yKy = np.square(region_Xty_jack).sum()/M_region
                        yy = np.inner(cur_pheno, cur_pheno)
                        norm_eq_b = np.array([c_yKy, yy])


                        sigma_vec = np.linalg.solve(norm_eq_A, norm_eq_b)

                        h2_vec = sigma_vec/sum(sigma_vec)
                        if cur_pheno_name not in h2_vec_std_dict:
                            h2_vec_std_dict[cur_pheno_name] = [h2_vec]
                        else:
                            h2_vec_std_dict[cur_pheno_name].append(h2_vec)
                
                for k in range(self.num_pheno):
                    cur_pheno_name = self.pheno_names[k]
                    h2_vec_std_list = h2_vec_std_dict[cur_pheno_name]
                    sigma_vec_std = np.std(np.vstack(h2_vec_std_list), axis=0) * np.sqrt(jackknife_blocks-1)
                    h2_vec_std_dict[cur_pheno_name] = sigma_vec_std
                
                if self.binary_trait:
                    for k in range(self.num_pheno):
                        cur_pheno_name = self.pheno_names[k]
                        h2_vec = h2_vec_dict[cur_pheno_name]
                        h2_vec_std = h2_vec_std_dict[cur_pheno_name]
                        K = self.prevalence_dict[cur_pheno_name]
                        z = norm.pdf(norm.ppf(K))
                        h2_vec = h2_vec*K*(1-K)/np.square(z)
                        h2_vec[-1] = 1 - np.sum(h2_vec[:-1])
                        h2_vec_std = h2_vec_std*K*(1-K)/np.square(z)
                        h2_vec_std[-1] = 1 - np.sum(h2_vec_std[:-1])
                        h2_vec_dict[cur_pheno_name] = h2_vec
                        h2_vec_std_dict[cur_pheno_name] = h2_vec_std

                return h2_vec_dict, h2_vec_std_dict


            return h2_vec_dict
        else:
            region_Xty = self.compute_Xty()



            tr_KK = (region_ld_proj.T @ region_ld_proj).sum()/(B*M_region*M_region)

            b_trK = self.N
            Nc = self.N
            norm_eq_A = np.array([[tr_KK, b_trK],
                                  [b_trK,   Nc]])

            c_yKy = np.square(region_Xty).sum()/M_region
            yy = np.inner(self.pheno, self.pheno)
            norm_eq_b = np.array([c_yKy, yy])

            print(f"norm_eq_A: \n{norm_eq_A}")
            print(f"norm_eq_b: \n{norm_eq_b}")
            sigma_vec = np.linalg.solve(norm_eq_A, norm_eq_b)
            h2_vec = sigma_vec/sum(sigma_vec)

            if self.binary_trait:
                cur_pheno_name = self.pheno_name
                K = self.prevalence_dict[cur_pheno_name]
                z = norm.pdf(norm.ppf(K))
                h2_vec = h2_vec*K*(1-K)/np.square(z)
                h2_vec[-1] = 1 - np.sum(h2_vec[:-1])

            return h2_vec
    
    def estimate_joint_her(self, region_indices=None, num_vecs=None, num_workers=None, step_size=None,
                             jackknife_blocks=None):
        if region_indices is not None: self.region_indices = region_indices
        if num_vecs: self.num_vecs = num_vecs
        if num_workers: self.num_workers = num_workers
        if step_size: self.step_size = step_size
        
        B = self.num_vecs
        M_region = len(self.region_indices)
        
        region_ld_proj = self.compute_ld_proj()
    
        if self.multi_pheno:
            C = self.num_pheno
            all_region_Xty = self.compute_Xty()
            tr_KK = (region_ld_proj.T @ region_ld_proj).sum()/(B*M_region*M_region)

            N = self.N

            norm_eq_A = np.zeros((2*C+2, 2*C+2))
            for i in range(2*C+2):
                if i == 0:
                    norm_eq_A[0, i] = norm_eq_A[i, 0] = (C*C)*tr_KK
                elif i < (C+1):
                    norm_eq_A[0, i] = norm_eq_A[i, 0] =  tr_KK
                else:
                    norm_eq_A[0, i] = norm_eq_A[i, 0] = N
            for i in range(2*C+2):
                if i == 0:
                    norm_eq_A[(C+1), i] = norm_eq_A[i, (C+1)] = C*C*N
                else:
                    norm_eq_A[(C+1), i] = norm_eq_A[i, (C+1)] = N
            norm_eq_A[(C+1),(C+1)] = C*C*N
            
            norm_eq_A[1:(C+1), 1:(C+1)] = np.eye(C)*tr_KK
            norm_eq_A[1:(C+1), (C+2):(2*C+2)] = norm_eq_A[(C+2):(2*C+2), 1:(C+1)] = np.eye(C)*N
            norm_eq_A[(C+2):(2*C+2), (C+2):(2*C+2)] = np.eye(C)*N
            
            norm_eq_b = np.zeros((2*C+2, 1))
            y_col_sum = self.pheno.sum(axis=1)

            norm_eq_b[C+1] = np.inner(y_col_sum, y_col_sum)
            for i in range(C):
                norm_eq_b[C+2+i] = np.inner(self.pheno[:, i], self.pheno[:, i])
                norm_eq_b[i+1] = np.square(all_region_Xty[:, i]).sum()/M_region

            all_pheno = self.pheno
            all_num_pheno = self.num_pheno
            self.pheno = y_col_sum
            self.num_pheno = 1
            Xty_sum = self.compute_Xty()
            self.pheno = all_pheno
            self.num_pheno = all_num_pheno
            norm_eq_b[0] = np.square(Xty_sum).sum()/M_region

            # print(f"norm_eq_A: \n{norm_eq_A}")
            # print(f"norm_eq_b: \n{norm_eq_b}")
            sigma_vec = np.linalg.solve(norm_eq_A, norm_eq_b)
            
            if jackknife_blocks != None:
                print(f"Start calculating jackknife subsample estimates, num of jack: {jackknife_blocks}")
                sigma_vec_jack_list = []
                num_snps = M_region // jackknife_blocks
                for j in tqdm(range(jackknife_blocks), total=jackknife_blocks):
                    if j == (jackknife_blocks - 1):
                        num_snps = (M_region // jackknife_blocks) + (M_region % jackknife_blocks)
                    start_index, end_index = j*(M_region // jackknife_blocks), j*(M_region // jackknife_blocks)+num_snps

                    region_ld_proj_jack = np.delete(region_ld_proj, np.arange(start_index, end_index), 0)
                    all_region_Xty_jack = np.delete(all_region_Xty, np.arange(start_index, end_index), 0)
                    Xty_sum_jack = np.delete(Xty_sum, np.arange(start_index, end_index), None)

                    M_region_jack = M_region - num_snps
                    tr_KK_jack = (region_ld_proj_jack.T @ region_ld_proj_jack).sum()/(B*M_region_jack*M_region_jack)

                    norm_eq_A = np.zeros((2*C+2, 2*C+2))
                    for i in range(2*C+2):
                        if i == 0:
                            norm_eq_A[0, i] = norm_eq_A[i, 0] = (C*C)*tr_KK_jack
                        elif i < (C+1):
                            norm_eq_A[0, i] = norm_eq_A[i, 0] =  tr_KK_jack
                        else:
                            norm_eq_A[0, i] = norm_eq_A[i, 0] = N
                    for i in range(2*C+2):
                        if i == 0:
                            norm_eq_A[(C+1), i] = norm_eq_A[i, (C+1)] = C*C*N
                        else:
                            norm_eq_A[(C+1), i] = norm_eq_A[i, (C+1)] = N
                    norm_eq_A[(C+1),(C+1)] = C*C*N
                    
                    norm_eq_A[1:(C+1), 1:(C+1)] = np.eye(C)*tr_KK_jack
                    norm_eq_A[1:(C+1), (C+2):(2*C+2)] = norm_eq_A[(C+2):(2*C+2), 1:(C+1)] = np.eye(C)*N
                    norm_eq_A[(C+2):(2*C+2), (C+2):(2*C+2)] = np.eye(C)*N
                    
                    norm_eq_b = np.zeros((2*C+2, 1))

                    norm_eq_b[C+1] = np.inner(y_col_sum, y_col_sum)
                    for i in range(C):
                        norm_eq_b[C+2+i] = np.inner(self.pheno[:, i], self.pheno[:, i])
                        norm_eq_b[i+1] = np.square(all_region_Xty_jack[:, i]).sum()/M_region_jack

                    norm_eq_b[0] = np.square(Xty_sum_jack).sum()/M_region_jack

                    sigma_vec_jack = np.linalg.solve(norm_eq_A, norm_eq_b)
                    sigma_vec_jack_list.append(sigma_vec_jack)

                sigma_vec_std = np.std(np.hstack(sigma_vec_jack_list), axis=1) * np.sqrt(jackknife_blocks-1)
                return sigma_vec, sigma_vec_std

            return sigma_vec
            
        else:
            return self.estimate_regional_her(region_indices, num_vecs, num_workers, step_size,
                             jackknife_blocks)


    def estimate_joint_her_full(self, region_indices=None, num_vecs=None, num_workers=None, step_size=None,
                             jackknife_blocks=None):
        if region_indices is not None: self.region_indices = region_indices
        if num_vecs: self.num_vecs = num_vecs
        if num_workers: self.num_workers = num_workers
        if step_size: self.step_size = step_size
        
        B = self.num_vecs
        M_region = len(self.region_indices)
        
        region_ld_proj = self.compute_ld_proj()
    
        if self.multi_pheno:
            C = self.num_pheno
            all_region_Xty = self.compute_Xty()
            tr_KK = (region_ld_proj.T @ region_ld_proj).sum()/(B*M_region*M_region)

            N = self.N

            dim_norm_eq = C*C+C+2
            norm_eq_A = np.zeros((dim_norm_eq, dim_norm_eq))
            norm_eq_A[0, 0] = (C*C)*tr_KK
            norm_eq_A[int(dim_norm_eq/2),int(dim_norm_eq/2)] = (C*C)*N

            tmp_idx = 0
            for i in range(C):
                for j in range(i, C):
                    if i == j:
                        norm_eq_A[0, int(tmp_idx+1)] = norm_eq_A[int(tmp_idx+1), 0] =  tr_KK
                        norm_eq_A[0, int(dim_norm_eq/2+tmp_idx+1)] = norm_eq_A[int(dim_norm_eq/2+tmp_idx+1), 0] =  N
                    else:
                        norm_eq_A[0, int(tmp_idx+1)] = norm_eq_A[int(tmp_idx+1), 0] =  2*tr_KK
                        norm_eq_A[0, int(dim_norm_eq/2+tmp_idx+1)] = norm_eq_A[int(dim_norm_eq/2+tmp_idx+1), 0] = 2*N
                    tmp_idx += 1

            norm_eq_A[0, int(dim_norm_eq/2)] = norm_eq_A[int(dim_norm_eq/2), 0] = (C*C)*N

            tmp_idx = 0
            for i in range(C):
                for j in range(i, C):
                    if i == j:
                        norm_eq_A[int(dim_norm_eq/2), int(tmp_idx+1)] = norm_eq_A[int(tmp_idx+1), int(dim_norm_eq/2)] =  N
                        norm_eq_A[int(dim_norm_eq/2), int(dim_norm_eq/2+tmp_idx+1)] = norm_eq_A[int(dim_norm_eq/2+tmp_idx+1), int(dim_norm_eq/2)] =  N
                    else:
                        norm_eq_A[int(dim_norm_eq/2), int(tmp_idx+1)] = norm_eq_A[int(tmp_idx+1), int(dim_norm_eq/2)] =  2*N
                        norm_eq_A[int(dim_norm_eq/2), int(dim_norm_eq/2+tmp_idx+1)] = norm_eq_A[int(dim_norm_eq/2+tmp_idx+1), int(dim_norm_eq/2)] = 2*N
                    tmp_idx += 1
            

            norm_eq_A[1:int(dim_norm_eq/2), 1:int(dim_norm_eq/2)] = np.diag(norm_eq_A[0, 1:int(dim_norm_eq/2)])
            norm_eq_A[1:int(dim_norm_eq/2), (int(dim_norm_eq/2)+1):] = np.diag(norm_eq_A[0, (int(dim_norm_eq/2)+1):])
            norm_eq_A[(int(dim_norm_eq/2)+1):, 1:int(dim_norm_eq/2)] = np.diag(norm_eq_A[0, (int(dim_norm_eq/2)+1):])
            norm_eq_A[(int(dim_norm_eq/2)+1):, (int(dim_norm_eq/2)+1):] = np.diag(norm_eq_A[0, (int(dim_norm_eq/2)+1):])
            

            norm_eq_b = np.zeros((dim_norm_eq, 1))
            y_col_sum = self.pheno.sum(axis=1)

            norm_eq_b[int((dim_norm_eq)/2)] = np.inner(y_col_sum, y_col_sum)

            tmp_idx = 0
            for i in range(C):
                pheno_i = self.pheno[:, i]
                pheno_K_i = all_region_Xty[:, i]
                for j in range(i, C):
                    pheno_j = self.pheno[:, j]
                    pheno_K_j = all_region_Xty[:, j]
                    if i == j:
                        norm_eq_b[int(dim_norm_eq/2 + tmp_idx + 1)] = np.inner(pheno_i, pheno_j)
                        norm_eq_b[int(tmp_idx + 1)] = np.inner(pheno_K_i, pheno_K_j)/M_region
                    else:
                        norm_eq_b[int(dim_norm_eq/2 + tmp_idx + 1)] = 2*np.inner(pheno_i, pheno_j)
                        norm_eq_b[int(tmp_idx + 1)] = 2*np.inner(pheno_K_i, pheno_K_j)/M_region
                    tmp_idx += 1

            all_pheno = self.pheno
            all_num_pheno = self.num_pheno
            self.pheno = y_col_sum
            self.num_pheno = 1
            Xty_sum = self.compute_Xty()
            self.pheno = all_pheno
            self.num_pheno = all_num_pheno
            norm_eq_b[0] = np.square(Xty_sum).sum()/M_region

            # print(f"norm_eq_A: \n{norm_eq_A}")
            # print(f"norm_eq_b: \n{norm_eq_b}")
            norm_eq_A = np.delete(norm_eq_A, [0, int(dim_norm_eq/2)], axis=0)
            norm_eq_A = np.delete(norm_eq_A, [0, int(dim_norm_eq/2)], axis=1)
            norm_eq_b = np.delete(norm_eq_b, [0, int(dim_norm_eq/2)], axis=0)

            sigma_vec = np.linalg.solve(norm_eq_A, norm_eq_b)
            all_norm_eq_A = norm_eq_A
            all_norm_eq_b = norm_eq_b
            
            if jackknife_blocks != None:
                print(f"Start calculating jackknife subsample estimates, num of jack: {jackknife_blocks}")
                sigma_vec_jack_list = []
                num_snps = M_region // jackknife_blocks
                for j in tqdm(range(jackknife_blocks), total=jackknife_blocks):
                    if j == (jackknife_blocks - 1):
                        num_snps = (M_region // jackknife_blocks) + (M_region % jackknife_blocks)
                    start_index, end_index = j*(M_region // jackknife_blocks), j*(M_region // jackknife_blocks)+num_snps

                    region_ld_proj_jack = np.delete(region_ld_proj, np.arange(start_index, end_index), 0)
                    all_region_Xty_jack = np.delete(all_region_Xty, np.arange(start_index, end_index), 0)
                    Xty_sum_jack = np.delete(Xty_sum, np.arange(start_index, end_index), None)

                    M_region_jack = M_region - num_snps
                    tr_KK_jack = (region_ld_proj_jack.T @ region_ld_proj_jack).sum()/(B*M_region_jack*M_region_jack)

                    norm_eq_A = np.zeros((dim_norm_eq, dim_norm_eq))
                    norm_eq_A[0, 0] = (C*C)*tr_KK
                    norm_eq_A[int(dim_norm_eq/2),int(dim_norm_eq/2)] = (C*C)*N

                    tmp_idx = 0
                    for i in range(C):
                        for j in range(i, C):
                            if i == j:
                                norm_eq_A[0, int(tmp_idx+1)] = norm_eq_A[int(tmp_idx+1), 0] =  tr_KK_jack
                                norm_eq_A[0, int(dim_norm_eq/2+tmp_idx+1)] = norm_eq_A[int(dim_norm_eq/2+tmp_idx+1), 0] =  N
                            else:
                                norm_eq_A[0, int(tmp_idx+1)] = norm_eq_A[int(tmp_idx+1), 0] =  2*tr_KK_jack
                                norm_eq_A[0, int(dim_norm_eq/2+tmp_idx+1)] = norm_eq_A[int(dim_norm_eq/2+tmp_idx+1), 0] = 2*N
                            tmp_idx += 1


                    norm_eq_A[0, int(dim_norm_eq/2)] = norm_eq_A[int(dim_norm_eq/2), 0] = (C*C)*N

                    tmp_idx = 0
                    for i in range(C):
                        for j in range(i, C):
                            if i == j:
                                norm_eq_A[int(dim_norm_eq/2), int(tmp_idx+1)] = norm_eq_A[int(tmp_idx+1), int(dim_norm_eq/2)] =  N
                                norm_eq_A[int(dim_norm_eq/2), int(dim_norm_eq/2+tmp_idx+1)] = norm_eq_A[int(dim_norm_eq/2+tmp_idx+1), int(dim_norm_eq/2)] =  N
                            else:
                                norm_eq_A[int(dim_norm_eq/2), int(tmp_idx+1)] = norm_eq_A[int(tmp_idx+1), int(dim_norm_eq/2)] =  2*N
                                norm_eq_A[int(dim_norm_eq/2), int(dim_norm_eq/2+tmp_idx+1)] = norm_eq_A[int(dim_norm_eq/2+tmp_idx+1), int(dim_norm_eq/2)] = 2*N
                            tmp_idx += 1


                    norm_eq_A[1:int(dim_norm_eq/2), 1:int(dim_norm_eq/2)] = np.diag(norm_eq_A[0, 1:int(dim_norm_eq/2)])
                    norm_eq_A[1:int(dim_norm_eq/2), (int(dim_norm_eq/2)+1):] = np.diag(norm_eq_A[0, (int(dim_norm_eq/2)+1):])
                    norm_eq_A[(int(dim_norm_eq/2)+1):, 1:int(dim_norm_eq/2)] = np.diag(norm_eq_A[0, (int(dim_norm_eq/2)+1):])
                    norm_eq_A[(int(dim_norm_eq/2)+1):, (int(dim_norm_eq/2)+1):] = np.diag(norm_eq_A[0, (int(dim_norm_eq/2)+1):])
                    

                    norm_eq_b = np.zeros((dim_norm_eq, 1))

                    norm_eq_b[0] = np.square(Xty_sum_jack).sum()/M_region
                    norm_eq_b[int((dim_norm_eq)/2)] = np.inner(y_col_sum, y_col_sum)

                    tmp_idx = 0
                    for i in range(C):
                        pheno_i = self.pheno[:, i]
                        pheno_K_i = all_region_Xty_jack[:, i]
                        for j in range(i, C):
                            pheno_j = self.pheno[:, j]
                            pheno_K_j = all_region_Xty_jack[:, j]
                            if i == j:
                                norm_eq_b[int(dim_norm_eq/2 + tmp_idx + 1)] = np.inner(pheno_i, pheno_j)
                                norm_eq_b[int(tmp_idx + 1)] = np.inner(pheno_K_i, pheno_K_j)/M_region
                            else:
                                norm_eq_b[int(dim_norm_eq/2 + tmp_idx + 1)] = 2*np.inner(pheno_i, pheno_j)
                                norm_eq_b[int(tmp_idx + 1)] = 2*np.inner(pheno_K_i, pheno_K_j)/M_region
                            tmp_idx += 1

                    norm_eq_A = np.delete(norm_eq_A, [0, int(dim_norm_eq/2)], axis=0)
                    norm_eq_A = np.delete(norm_eq_A, [0, int(dim_norm_eq/2)], axis=1)
                    norm_eq_b = np.delete(norm_eq_b, [0, int(dim_norm_eq/2)], axis=0)

                    sigma_vec_jack = np.linalg.solve(norm_eq_A, norm_eq_b)
                    sigma_vec_jack_list.append(sigma_vec_jack)

                sigma_vec_std = np.std(np.hstack(sigma_vec_jack_list), axis=1) * np.sqrt(jackknife_blocks-1)
                return sigma_vec, sigma_vec_std, all_norm_eq_A, all_norm_eq_b

            return sigma_vec, all_norm_eq_A, all_norm_eq_b
            
        else:
            return self.estimate_regional_her(region_indices, num_vecs, num_workers, step_size,
                             jackknife_blocks)



    def estimate_partitioned_regional_her(self, annot_mat, region_indices=None, num_vecs=None, num_workers=None, 
                                          step_size=None, jackknife_blocks=None, debug_flag=False):
        if region_indices is not None: self.region_indices = region_indices
        if num_vecs: self.num_vecs = num_vecs
        if num_workers: self.num_workers = num_workers
        if step_size: self.step_size = step_size
        
        B = self.num_vecs
        
        num_partition = annot_mat.shape[1]
        partitioned_region_indices = []
        actual_region_indices = []
        for k in range(num_partition):
            indices = np.where(annot_mat[:, k] == 1)[0]
            actual_region_indices.append(indices)
            partitioned_region_indices.append([
                np.where(self.region_indices == idx)[0][0] for idx in indices
            ])
#         print(partitioned_region_indices)
        
        M_region = len(self.region_indices)
        
        all_region_ld_proj = self.compute_partitioned_ld_proj(actual_region_indices)
    
        if self.multi_pheno:
            all_region_Xty = self.compute_Xty()
            
            norm_eq_A = np.zeros((num_partition+1, num_partition+1))
            
            norm_eq_A[:, num_partition] = [self.N]*(num_partition+1)
            norm_eq_A[num_partition, :] = [self.N]*(num_partition+1)
            
            for k in range(num_partition):
                cur_region_ld_proj = all_region_ld_proj[k]
                cur_partitioned_region_indices = partitioned_region_indices[k]
                cur_M_region = len(cur_partitioned_region_indices)
                cur_region_ld_proj = cur_region_ld_proj[cur_partitioned_region_indices]
                norm_eq_A[k, k] = (cur_region_ld_proj.T @ 
                                   cur_region_ld_proj).sum()/(B*cur_M_region*cur_M_region)

            for k in range(num_partition):
                M_region_k = len(partitioned_region_indices[k])
                cur_region_ld_proj = all_region_ld_proj[k]
                for j in range(k+1, num_partition):
                    M_region_j = len(partitioned_region_indices[j])
                    cross_region_ld_proj = cur_region_ld_proj[partitioned_region_indices[j]]
                    norm_eq_A[k, j] = norm_eq_A[j, k] = (cross_region_ld_proj.T @
                                                        cross_region_ld_proj).sum()/(B*M_region_k*M_region_j)
            
            h2_vec_dict = {}
            for k in tqdm(range(self.num_pheno), total=self.num_pheno):
                cur_pheno = self.pheno[:, k]
                region_Xty = all_region_Xty[:, k]
                cur_pheno_name = self.pheno_names[k]
                
                yy = np.inner(cur_pheno, cur_pheno)
                norm_eq_b = np.zeros((num_partition+1,))
                
                norm_eq_b[-1] = yy
                for j in range(num_partition):
                    M_region_j = len(partitioned_region_indices[j])
                    region_Xty_j = region_Xty[partitioned_region_indices[j]]
                    norm_eq_b[j] = np.square(region_Xty_j).sum()/M_region_j


                sigma_vec = np.linalg.solve(norm_eq_A, norm_eq_b)

                h2_vec = sigma_vec/sum(sigma_vec)
                h2_vec_dict[cur_pheno_name] = h2_vec

            if self.binary_trait:
                for k in range(self.num_pheno):
                    cur_pheno_name = self.pheno_names[k]
                    h2_vec = h2_vec_dict[cur_pheno_name]
                    K = self.prevalence_dict[cur_pheno_name]
                    z = norm.pdf(norm.ppf(K))
                    h2_vec = h2_vec*K*(1-K)/np.square(z)
                    h2_vec[-1] = 1 - np.sum(h2_vec[:-1])
                    h2_vec_dict[cur_pheno_name] = h2_vec


            if debug_flag:
                return h2_vec_dict, norm_eq_A, norm_eq_b
            else:
                return h2_vec_dict
        else:
            all_region_Xty = self.compute_Xty()
            
            norm_eq_A = np.zeros((num_partition+1, num_partition+1))
            
            norm_eq_A[:, num_partition] = [self.N]*(num_partition+1)
            norm_eq_A[num_partition, :] = [self.N]*(num_partition+1)
            
            for k in range(num_partition):
                cur_region_ld_proj = all_region_ld_proj[k]
                cur_partitioned_region_indices = partitioned_region_indices[k]
                cur_M_region = len(cur_partitioned_region_indices)
                cur_region_ld_proj = cur_region_ld_proj[cur_partitioned_region_indices]
                norm_eq_A[k, k] = (cur_region_ld_proj.T @ 
                                   cur_region_ld_proj).sum()/(B*cur_M_region*cur_M_region)
#                 print(cur_M_region)

            for k in range(num_partition):
                M_region_k = len(partitioned_region_indices[k])
                cur_region_ld_proj = all_region_ld_proj[k]
                for j in range(k+1, num_partition):
                    M_region_j = len(partitioned_region_indices[j])
                    cross_region_ld_proj = cur_region_ld_proj[partitioned_region_indices[j]]
                    norm_eq_A[k, j] = norm_eq_A[j, k] = (cross_region_ld_proj.T @
                                                        cross_region_ld_proj).sum()/(B*M_region_k*M_region_j)
#                     print(M_region_k, M_region_j, k, j)
                
            
            

            pheno = self.pheno
            region_Xty = all_region_Xty

            yy = np.inner(pheno, pheno)
            norm_eq_b = np.zeros((num_partition+1,))

            norm_eq_b[-1] = yy
            for j in range(num_partition):
                M_region_j = len(partitioned_region_indices[j])
                region_Xty_j = region_Xty[partitioned_region_indices[j]]
                norm_eq_b[j] = np.square(region_Xty_j).sum()/M_region_j
#                 print(M_region_j)

            print(f"norm_eq_A: \n{norm_eq_A}")
            print(f"norm_eq_b: \n{norm_eq_b}")
            sigma_vec = np.linalg.solve(norm_eq_A, norm_eq_b)

            h2_vec = sigma_vec/sum(sigma_vec)

            if self.binary_trait:
                cur_pheno_name = self.pheno_name
                K = self.prevalence_dict[cur_pheno_name]
                z = norm.pdf(norm.ppf(K))
                h2_vec = h2_vec*K*(1-K)/np.square(z)
                h2_vec[-1] = 1 - np.sum(h2_vec[:-1])

            if debug_flag:
                return h2_vec, norm_eq_A, norm_eq_b
            else:
                return h2_vec
    

    def inverse(self, X):
        return pinvh(X.T@X)

    def projection(self, Z,X,P1):
        print(f"X shape: {X.shape}, Z shape: {Z.shape}")
    # Perform (I-X(X^TX)^-1 X^T)Z
        Z = np.array(Z,order='F')
        X = np.array(X,order='F')
        P1 = np.array(P1,order='F')
        t1 = scipy.linalg.blas.sgemm(1.,X,Z,trans_a=True)
        t2 = scipy.linalg.blas.sgemm(1.,X,P1)
        t3 = scipy.linalg.blas.sgemm(1.,t2,t1)
        Z  = Z - t3
        return Z

    ## remove flanking effects
    def estimate_partitioned_regional_her2(self, annot_mat, region_indices=None, num_vecs=None, num_workers=None, 
                                          step_size=None, jackknife_blocks=None, debug_flag=False):
        if region_indices is not None: self.region_indices = region_indices
        if num_vecs: self.num_vecs = num_vecs
        if num_workers: self.num_workers = num_workers
        if step_size: self.step_size = step_size
        
        B = self.num_vecs
        
        indices = np.argwhere(annot_mat[:, 0] == 1).flatten()
        cur_geno = self.G.read(index=np.s_[:, indices])
        means = np.nanmean(cur_geno, axis=0)
        stds = np.nanstd(cur_geno, axis=0)

        ## standardization
        cur_geno = (cur_geno - means) / stds
        ## mean imputation
        cur_geno[np.isnan(cur_geno)] = 0

        pca = PCA(n_components=min(200, cur_geno.shape[1]))
        cur_geno_pc = pca.fit_transform(cur_geno)
        self.surround_geno_pc = cur_geno_pc
        P1 = self.inverse(cur_geno_pc)
        self.surround_geno_pc_inv = P1
        if self.multi_pheno:
            pheno_vec = self.pheno
        else:
            pheno_vec = self.pheno.reshape(-1, 1)
        res_pheno = self.projection(pheno_vec, cur_geno_pc, P1)
        if self.multi_pheno:
            self.pheno = res_pheno
        else:
            self.pheno = res_pheno.flatten()
        annot_mat = annot_mat[:, 1:]
        print(f"removed flanking effects...")

        num_partition = annot_mat.shape[1]
        partitioned_region_indices = []
        actual_region_indices = []
        for k in range(num_partition):
            indices = np.where(annot_mat[:, k] == 1)[0]
            actual_region_indices.append(indices)
            partitioned_region_indices.append([
                np.where(self.region_indices == idx)[0][0] for idx in indices
            ])
#         print(partitioned_region_indices)
        
        M_region = len(self.region_indices)
        
        all_region_ld_proj = self.compute_partitioned_ld_proj(actual_region_indices)
    
        if self.multi_pheno:
            all_region_Xty = self.compute_Xty()
            
            norm_eq_A = np.zeros((num_partition+1, num_partition+1))
            
            norm_eq_A[:, num_partition] = [self.N]*(num_partition+1)
            norm_eq_A[num_partition, :] = [self.N]*(num_partition+1)
            
            for k in range(num_partition):
                cur_region_ld_proj = all_region_ld_proj[k]
                cur_partitioned_region_indices = partitioned_region_indices[k]
                cur_M_region = len(cur_partitioned_region_indices)
                cur_region_ld_proj = cur_region_ld_proj[cur_partitioned_region_indices]
                norm_eq_A[k, k] = (cur_region_ld_proj.T @ 
                                   cur_region_ld_proj).sum()/(B*cur_M_region*cur_M_region)

            for k in range(num_partition):
                M_region_k = len(partitioned_region_indices[k])
                cur_region_ld_proj = all_region_ld_proj[k]
                for j in range(k+1, num_partition):
                    M_region_j = len(partitioned_region_indices[j])
                    cross_region_ld_proj = cur_region_ld_proj[partitioned_region_indices[j]]
                    norm_eq_A[k, j] = norm_eq_A[j, k] = (cross_region_ld_proj.T @
                                                        cross_region_ld_proj).sum()/(B*M_region_k*M_region_j)
            
            h2_vec_dict = {}
            for k in tqdm(range(self.num_pheno), total=self.num_pheno):
                cur_pheno = self.pheno[:, k]
                region_Xty = all_region_Xty[:, k]
                cur_pheno_name = self.pheno_names[k]
                
                yy = np.inner(cur_pheno, cur_pheno)
                norm_eq_b = np.zeros((num_partition+1,))
                
                norm_eq_b[-1] = yy
                for j in range(num_partition):
                    M_region_j = len(partitioned_region_indices[j])
                    region_Xty_j = region_Xty[partitioned_region_indices[j]]
                    norm_eq_b[j] = np.square(region_Xty_j).sum()/M_region_j


                sigma_vec = np.linalg.solve(norm_eq_A, norm_eq_b)

                h2_vec = sigma_vec/sum(sigma_vec)
                h2_vec_dict[cur_pheno_name] = h2_vec

            if self.binary_trait:
                for k in range(self.num_pheno):
                    cur_pheno_name = self.pheno_names[k]
                    h2_vec = h2_vec_dict[cur_pheno_name]
                    K = self.prevalence_dict[cur_pheno_name]
                    z = norm.pdf(norm.ppf(K))
                    h2_vec = h2_vec*K*(1-K)/np.square(z)
                    h2_vec[-1] = 1 - np.sum(h2_vec[:-1])
                    h2_vec_dict[cur_pheno_name] = h2_vec


            if debug_flag:
                return h2_vec_dict, norm_eq_A, norm_eq_b
            else:
                return h2_vec_dict
        else:
            all_region_Xty = self.compute_Xty()
            
            norm_eq_A = np.zeros((num_partition+1, num_partition+1))
            
            norm_eq_A[:, num_partition] = [self.N]*(num_partition+1)
            norm_eq_A[num_partition, :] = [self.N]*(num_partition+1)
            
            for k in range(num_partition):
                cur_region_ld_proj = all_region_ld_proj[k]
                cur_partitioned_region_indices = partitioned_region_indices[k]
                cur_M_region = len(cur_partitioned_region_indices)
                cur_region_ld_proj = cur_region_ld_proj[cur_partitioned_region_indices]
                norm_eq_A[k, k] = (cur_region_ld_proj.T @ 
                                   cur_region_ld_proj).sum()/(B*cur_M_region*cur_M_region)
#                 print(cur_M_region)

            for k in range(num_partition):
                M_region_k = len(partitioned_region_indices[k])
                cur_region_ld_proj = all_region_ld_proj[k]
                for j in range(k+1, num_partition):
                    M_region_j = len(partitioned_region_indices[j])
                    cross_region_ld_proj = cur_region_ld_proj[partitioned_region_indices[j]]
                    norm_eq_A[k, j] = norm_eq_A[j, k] = (cross_region_ld_proj.T @
                                                        cross_region_ld_proj).sum()/(B*M_region_k*M_region_j)
#                     print(M_region_k, M_region_j, k, j)

            pheno = self.pheno
            region_Xty = all_region_Xty

            yy = np.inner(pheno, pheno)
            norm_eq_b = np.zeros((num_partition+1,))

            norm_eq_b[-1] = yy
            for j in range(num_partition):
                M_region_j = len(partitioned_region_indices[j])
                region_Xty_j = region_Xty[partitioned_region_indices[j]]
                norm_eq_b[j] = np.square(region_Xty_j).sum()/M_region_j
#                 print(M_region_j)

            print(f"norm_eq_A: \n{norm_eq_A}")
            print(f"norm_eq_b: \n{norm_eq_b}")
            sigma_vec = np.linalg.solve(norm_eq_A, norm_eq_b)
            h2_vec = sigma_vec/sum(sigma_vec)

            if self.binary_trait:
                cur_pheno_name = self.pheno_name
                K = self.prevalence_dict[cur_pheno_name]
                z = norm.pdf(norm.ppf(K))
                h2_vec = h2_vec*K*(1-K)/np.square(z)
                h2_vec[-1] = 1 - np.sum(h2_vec[:-1])

            if debug_flag:
                return h2_vec, norm_eq_A, norm_eq_b
            else:
                return h2_vec
       
    
    def compute_ld_score(self, gene_indices=None, region_indices=None, num_vecs=None, 
                                    num_workers=None, step_size=None):
        if gene_indices is None: gene_indices = self.region_indices
        if region_indices is not None: self.region_indices = region_indices
        if num_vecs: self.num_vecs = num_vecs
        if num_workers: self.num_workers = num_workers
        if step_size: self.step_size = step_size
            
        M_region = len(self.region_indices)
        
        print(f"...Start computing LD projections, with the set of parameters: ")
        print(f"\t num_vecs={self.num_vecs}, num_workers={self.num_workers}, step_size={self.step_size}")
        print(f"\t number of SNPs in this region: {M_region}")
        
        
        all_indices = np.arange(len(self.region_indices))[::self.step_size]
        num_blocks = len(all_indices)
        self.partitioned_indices = []
        for j in range(num_blocks):
            if j == (num_blocks - 1):
                self.partitioned_indices.append((j, self.region_indices[all_indices[j]:]))
            else:
                self.partitioned_indices.append((j, self.region_indices[all_indices[j]:all_indices[j+1]]))


        cur_regional_indices = gene_indices
        self.temp_results = np.zeros((self.N, 1))
        results = np.zeros((M_region, 1))

        all_partitioned_indices = np.arange(len(cur_regional_indices))[::self.step_size]
        num_blocks = len(all_partitioned_indices)
        cur_partitioned_indices = []
        for j in range(num_blocks):
            if j == (num_blocks - 1):
                cur_partitioned_indices.append((j, cur_regional_indices[all_partitioned_indices[j]:]))
            else:
                cur_partitioned_indices.append((j, 
                            cur_regional_indices[all_partitioned_indices[j]:all_partitioned_indices[j+1]]))

        for cur_indices in tqdm(cur_partitioned_indices, desc='Calculating X1', 
                                total=len(cur_partitioned_indices)):
            cur_res = self._compute_X1(cur_indices)
            self.temp_results += cur_res

        results_dict = {}

        for cur_indices in tqdm(self.partitioned_indices, desc='Calculating XtX1', total=len(self.partitioned_indices)):
            label, cur_res = self._compute_XtXz(cur_indices)
            results_dict[label] = cur_res

        num_blocks = len(all_indices)
        for j in range(num_blocks):
            if j == (num_blocks-1):
                start_index = all_indices[j]
                results[start_index:, :] = results_dict[j]
            else:
                start_index, end_index = all_indices[j], all_indices[j+1]
                results[start_index:end_index, :] = results_dict[j]
                
        results = results/(self.N*len(gene_indices))
        results = np.abs(results.flatten())
        return results

   
    def permute_regional_her(self, num_perm, region_indices=None, num_vecs=None, num_workers=None, step_size=None,
                             jackknife_blocks=None):
        if self.multi_pheno:
            h2_vec_dict = self.estimate_regional_her(region_indices=region_indices, num_vecs=num_vecs, 
                            num_workers=num_workers, step_size=step_size,
                            jackknife_blocks=jackknife_blocks)
            original_pheno = self.pheno
            original_pheno_names = self.pheno_names
            num_pheno = self.num_pheno
            # num_pheno = original_pheno.shape[1]
            permuted_pheno = np.zeros((self.N, num_pheno*num_perm))
            permuted_pheno_names = []
            for perm_idx in tqdm(range(num_perm), total=num_perm, desc='Permuting phenotypes'):
                np.random.seed(perm_idx)
                indices = np.random.permutation(self.N)
                cur_perm_pheno = original_pheno[indices]
                permuted_pheno[:, (num_pheno*perm_idx):(num_pheno*(perm_idx+1))] = cur_perm_pheno
                permuted_pheno_names.extend([f"{pheno_name}_perm{perm_idx}" for pheno_name in original_pheno_names])

            self.pheno = permuted_pheno
            self.num_pheno = num_pheno*num_perm
            self.pheno_names = np.array(permuted_pheno_names)
            perm_h2_vec_dict = self.estimate_regional_her(region_indices=region_indices, num_vecs=num_vecs, 
                            num_workers=num_workers, step_size=step_size,
                            jackknife_blocks=jackknife_blocks)
            self.pheno = original_pheno
            self.num_pheno = num_pheno
            self.pheno_names = original_pheno_names

            h2_pval_dict = {}
            perm_h2_mat = np.zeros((num_perm+1, len(self.pheno_names)))
            for k, pheno_name in enumerate(self.pheno_names):
                h2_value = h2_vec_dict[pheno_name][0]
                perm_h2_values = np.array([perm_h2_vec_dict[f"{pheno_name}_perm{k}"][0] for k in range(num_perm)])
                perm_pval = (np.sum(h2_value < perm_h2_values)+1)/(num_perm+1)
                h2_pval_dict[pheno_name] = np.array([h2_value, perm_pval])
                perm_h2_mat[:, k] = np.concatenate([[h2_value], perm_h2_values])
        
            return h2_pval_dict, perm_h2_mat
        else:
            h2_vec = self.estimate_regional_her(region_indices=region_indices, num_vecs=num_vecs, 
                            num_workers=num_workers, step_size=step_size,
                            jackknife_blocks=jackknife_blocks)
            original_pheno = self.pheno
            permuted_pheno = np.zeros((self.N, num_perm))
            permuted_pheno_names = []
            for perm_idx in range(num_perm):
                np.random.seed(perm_idx)
                indices = np.random.permutation(self.N)
                cur_perm_pheno = original_pheno[indices]
                permuted_pheno[:, perm_idx] = cur_perm_pheno
                permuted_pheno_names.append(f"perm{perm_idx}")

            self.pheno = permuted_pheno
            self.pheno_names = np.array(permuted_pheno_names)
            self.num_pheno = num_perm
            self.multi_pheno = True
            perm_h2_vec_dict = self.estimate_regional_her(region_indices=region_indices, num_vecs=num_vecs, 
                            num_workers=num_workers, step_size=step_size,
                            jackknife_blocks=jackknife_blocks)
            self.pheno = original_pheno
            self.pheno_names = None
            self.multi_pheno = False

            h2_value = h2_vec[0]
            perm_h2_values = np.array([perm_h2_vec_dict[f"perm{k}"][0] for k in range(num_perm)])
            perm_pval = (np.sum(h2_value < perm_h2_values)+1)/(num_perm+1)
            # perm_pval = np.mean(h2_value < perm_h2_values)
            h2_pval = np.array([h2_value, perm_pval])

            return h2_pval, perm_h2_values

    