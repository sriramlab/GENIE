This folder summarized some of the results presented in GENIE

```bash
results
├── info    # mapping tables
│   ├── tissue_info.csv # UKB trait to category
│   └── trait_info.csv  # GTEx tissue to category
├── README.md
├── real_data   # UKB real data results 
│   ├── array_snps  # N = 291,273 white British unrelated individuals, M = 454,207 array SNPs
│   │   ├── all.age.norm.txt            # E = age, standardized binary encoding
│   │   ├── all.sex.bin.def.txt         # E = sex, binary (0/1) coding
│   │   ├── all.sex.bin.flip.txt        # flipped binary (1/0) coding
│   │   ├── all.sex.norm.def.txt        # standardized binary coding
│   │   ├── all.sex.norm.flip.txt       # standardized flipped binary encoding
│   │   ├── all.smok.bin.def.txt        # E = smoking status
│   │   ├── all.smok.bin.flip.txt       
│   │   ├── all.smok.norm.def.txt       
│   │   ├── all.smok.norm.flip.txt      
│   │   ├── all.smok.norm.lev3.txt      # standardized ternary (0/1/2) encoding
│   │   ├── all.statin.bin.def.txt      # E = statin
│   │   ├── all.statin.bin.flip.txt     
│   │   ├── all.statin.norm.def.txt     
│   │   ├── all.statin.norm.flip.txt    
│   │   └── sex_stratified_analysis     # sex specific anlysis
│   │       ├── stratified_by_sex_age.tsv
│   │       ├── stratified_by_sex_smoking.tsv
│   │       └── stratified_by_sex_statin.tsv
│   ├── imp_snps    # M = 7,774,235 SNPs (MAF >= 0.1%)
│   │   ├── mafld   # paritioned by MAF/LD (20 annotations)
│   │   │   ├── imp.mafld.gxsex.def.txt
│   │   │   ├── imp.mafld.gxsex.flip.txt
│   │   │   ├── imp.mafld.gxsmok.txt
│   │   │   ├── imp.mafld.gxstatin.def.txt
│   │   │   └── imp.mafld.gxstatin.flip.txt
│   │   └── single  # single annotation
│   │       ├── imp.gxsex.def.txt
│   │       ├── imp.gxsex.flip.txt
│   │       ├── imp.gxsmok.singlevc.txt
│   │       ├── imp.gxstatin.def.txt
│   │       └── imp.gxstatin.flip.txt
│   ├── med_adjust  # phenotype adjusted for medication usage
│   │   ├── bp_diastolic_adj_results.tsv    # DBP, adjusted for blood pressure medication
│   │   ├── bp_systolic_adj_results.tsv     # SBP, adjusted for blood pressure medication
│   │   ├── cholesterol_adj_results.tsv     # TC, adjusted for statin usage
│   │   ├── ldl_adj_results.tsv             # LDL-C, adjusted for statin usage
│   │   ├── multi_env_results.tsv           # Es = (statin, TC, LDL-C), Y = HA1c
│   │   └── tissue.cholesterol.adj.age.tsv  # GxAge results after TC adjusted for statin
│   ├── robustness  # robustness of GENIE in simulation and real data analyses
│   │   ├── ldsc_50_traits_results.tsv  # compare h2g estimates with LDSC
│   │   ├── nxe_effect      # NxE effect (basic: G+GxE model, full: G+GxE+NxE model)
│   │   │   ├── age
│   │   │   │   ├── age.basic.model.txt
│   │   │   │   └── age.full.model.txt
│   │   │   ├── sex
│   │   │   │   ├── sex.basic.model.def.txt
│   │   │   │   ├── sex.basic.model.flip.txt
│   │   │   │   ├── sex.full.model.def.txt
│   │   │   │   └── sex.full.model.flip.txt
│   │   │   ├── smok
│   │   │   │   ├── smok.basic.model.txt
│   │   │   │   └── smok.full.model.txt
│   │   │   └── statin
│   │   │       ├── gxe.basic.model.txt
│   │   │       └── gxe.full.model.txt
│   │   └── permute_test    # results after permutation of G (genotype), E (environment) or P (phenotype)
│   │       ├── sex
│   │       │   ├── basic.permute.E.txt
│   │       │   ├── basic.permute.G.txt
│   │       │   ├── basic.permute.P.txt
│   │       │   ├── full.permute.E.txt
│   │       │   ├── full.permute.G.txt
│   │       │   └── full.permute.P.txt
│   │       ├── smok
│   │       │   ├── basic.permute.E.txt
│   │       │   ├── basic.permute.G.txt
│   │       │   ├── basic.permute.P.txt
│   │       │   ├── full.permute.E.txt
│   │       │   ├── full.permute.G.txt
│   │       │   └── full.permute.P.txt
│   │       └── statin
│   │           ├── basic.permute.E.txt
│   │           ├── basic.permute.G.txt
│   │           └── basic.permute.P.txt
│   └── tissue      # enrichment in tissue-specific gene annotations
│       ├── age         # Over: based on overlapping annotations, noOver: non-overlapping annotations
│       │   ├── all.g.noOver.txt
│       │   ├── all.g.Over.txt
│       │   ├── all.gxe.noOver.txt
│       │   └── all.gxe.Over.txt
│       ├── sex
│       │   ├── all.g.noOver.txt
│       │   ├── all.g.Over.txt
│       │   ├── all.gxe.noOver.txt
│       │   ├── all.gxe.Over.txt
│       │   └── analysed.phen.txt
│       ├── smok
│       │   ├── all.g.noOver.txt
│       │   ├── all.g.Over.txt
│       │   ├── all.gxe.noOver.txt
│       │   └── all.gxe.Over.txt
│       └── statin
│           ├── all.g.noOver.txt
│           ├── all.g.Over.txt
│           ├── all.gxe.noOver.txt
│           └── all.gxe.Over.txt
└── simulation
    ├── benchmark               # benchmark with existing GxE methods on individual-level genotype
    │   ├── benchmark_results_filtered_geno_discrete.tsv    # independent SNPs, discrete E
    │   ├── benchmark_results_filtered_geno.tsv             # continuous E, fpr: P(p<0.05), fpr.2: P(p<0.05/200)
    │   ├── benchmark_results_fpr_adj.tsv                   # fpr: P(p<0.05/200)
    │   ├── benchmark_results.tsv                           # fpr: P(p<0.05)
    │   └── benchmark_runtime.tsv  
    ├── jackknife_SE_vs_true_SE # estimated vs true SE
    │   ├── estimated_SE_vs_true_SE_summary_gxe.tsv         # varying h2gxe
    │   ├── estimated_SE_vs_true_SE_summary_sample_size.tsv # varying sample size
    │   ├── estimated_SE_vs_true_SE_summary.tsv             # varying sample size/h2gxe
    │   ├── estimated_vs_true_SE_no_nxe.tsv                 # estimated SE over 100 replicates when no NxE in simulation
    │   └── estimated_vs_true_SE.tsv                        # with NxE in simulation (1st row: h2gxe, 2nd row: sample size, last row: true SE)
    ├── mean_impute_effect      # impact of mean imputation
    │   ├── simul_missing_env_results.tsv                   # if impute missing E
    │   └── simul_missing_pheno_results.tsv                 # if impute missing P
    ├── nxe_effect              # NxE effect (G+GxE vs G+GxE+NxE)
    │   └── fig2.txt
    ├── power                   # Power in simulation
    │   ├── power_vary_h2gxe.tsv                            # varying h2gxe
    │   ├── power_vary_h2gxe_values.tsv                     # estimated h2gxe over 100 replicates                     
    │   ├── sample_size_power_results_no_nxe.tsv            # varying sample size (no NxE)
    │   └── sample_size_power_results.tsv                   # with NxE
    ├── rand_impact             # impact of randomization      
    │   ├── compare_w_exactMoM.rand_SE.tsv                  # comparison with exact MoM
    │   ├── rand_SE_sample_size.tsv                         # fraction of SE due to randomization
    │   └── rand_vec_10_100.tsv                             # B=10 vs B=100
    └── robustness              # robustness in simulation (various scenarios)
        └── various_conditions_results.tsv
```
