# RHE

Randomized Haseman–Elston regression for Multi-variance Components: Python wrappers for **G**ene-**EN**viroment **I**nteraction **E**stimator


## Installation

```
pip install rhe
```

This package provides basic Python wrappers for the GENIE executables.
Full documentation for the GENIE project, as well as sample data necessary to run the toy examples below, is available here:
https://github.com/sriramlab/GENIE


## Interface

This package provides the following functions:

- `rhe.run_genie`


### Run GENIE with a list of command line arguments

```py
from rhe import run_genie

gen = 'path/to/test'
phen = 'path/to/test.pheno'
covar = 'path/to/test.cov'
annot = 'path/to/single.annot'
env = 'path/to/test.env'

run_genie(
    ['-g', gen, '-p', phen, '-c', covar, '-e', env, '-m', 'G+GxE+NxE', '-k', '10', '-jn', '10', '-o',
     'path/to/analysis.out', '-annot', annot, '-t', '6']
)
```


### Run GENIE with a single string of command line arguments

```py
from rhe import run_genie

gen = 'path/to/test'
phen = 'path/to/test.pheno'
covar = 'path/to/test.cov'
annot = 'path/to/single.annot'
env = 'path/to/test.env'

run_genie(
    f'-g {gen} -p {phen} -c {covar} -e {env} -m G+GxE+NxE -k 10 -jn 10 -o path/to/analysis.out -annot {annot} -t 6'
)
```


### Run GENIE with using a config file

```py
from rhe import run_genie_mem

with open('config.txt', 'w') as f:
    f.write('genotype=path/to/test\n')
    f.write('phenotype=path/to/test.pheno\n')
    f.write('covariate=path/to/test.cov\n')
    f.write('environment=path/to/test.env\n')
    f.write('annotation=path/to/single.annot\n')
    f.write('output=path/to/test.py.list.3.out\n')
    f.write('nthreads=6\n')
    f.write('num_vec=20\n')
    f.write('num_jack=10\n')
    f.write('trace=1\n')
    f.write('model=G\n')
    f.write('verbose=1\n')

run_genie_mem(
    ['--config', 'config.list.txt']
)
```


## Citation
```
1. Ali Pazokitoroudi, Andrew Dahl, Noah Zaitlen, Saharon Rosset, Sriram Sankararaman.
bioRxiv 2023.12.12.571316; doi: https://doi.org/10.1101/2023.12.12.571316
2. Ali Pazokitoroudi, Yue Wu, Kathryn S. Burch, Kangcheng Hou, Aaron Zhou, 
Bogdan Pasaniuc, Sriram Sankararaman. Nature Communications (2020); doi: https://doi.org/10.1101/522003
```


## Version
```
v1.2.0
```
