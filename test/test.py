from pathlib import Path
from rhe import run_genie

REPO_BASE_DIR = Path(__file__).parent.parent
EXAMPLE_DIR = REPO_BASE_DIR / 'example'

gen = EXAMPLE_DIR / 'test'
phen = EXAMPLE_DIR / 'test.pheno'
covar = EXAMPLE_DIR / 'test.cov'
annot = EXAMPLE_DIR / 'single.annot'
env = EXAMPLE_DIR / 'test.env'

run_genie(
    ['-g', gen, '-p', phen, '-c', covar, '-e', env, '-m', 'G+GxE+NxE', '-k', '10', '-jn', '10', '-o',
     EXAMPLE_DIR / 'test.py.list.1.out', '-annot', annot, '-t', '6']
)

run_genie(
    f'-g {gen} -p {phen} -c {covar} -e {env} -m G+GxE+NxE -k 10 -jn 10 -o {EXAMPLE_DIR / "test.py.str.1.out"} --annot {annot} -t 6'
)
