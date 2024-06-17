from pathlib import Path
from rhe import run_genie, run_genie_mem, run_genie_multi_pheno

REPO_BASE_DIR = Path(__file__).parent.parent
EXAMPLE_DIR = REPO_BASE_DIR / 'example'


def test_genie():
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
        f'-g {gen} -p {phen} -c {covar} -e {env} -m G+GxE+NxE -k 10 -jn 10 -o {EXAMPLE_DIR / "test.py.str.1.out"} -annot {annot} -t 6'
    )


def test_genie_multi_pheno():
    gen = EXAMPLE_DIR / 'test'
    phen = EXAMPLE_DIR / 'test.2.pheno'
    covar = EXAMPLE_DIR / 'test.cov'
    annot = EXAMPLE_DIR / 'single.annot'
    env = EXAMPLE_DIR / 'test.env'

    run_genie_multi_pheno(
        ['-g', gen, '-p', phen, '-c', covar, '-e', env, '-m', 'G+GxE+NxE', '-k', '10', '-jn', '10', '-o',
         EXAMPLE_DIR / 'test.py.list.2.out', '-annot', annot, '-t', '6']
    )

    run_genie_multi_pheno(
        f'-g {gen} -p {phen} -c {covar} -e {env} -m G+GxE+NxE -k 10 -jn 10 -o {EXAMPLE_DIR / "test.py.str.2.out"} -annot {annot} -t 6'
    )


def test_genie_mem():

    with open('config.list.txt', 'w') as f:
        f.write(f'genotype={EXAMPLE_DIR}/test\n')
        f.write(f'phenotype={EXAMPLE_DIR}/test.pheno\n')
        f.write(f'covariate={EXAMPLE_DIR}/test.cov\n')
        f.write(f'environment={EXAMPLE_DIR}/test.env\n')
        f.write(f'annotation={EXAMPLE_DIR}/single.annot\n')
        f.write(f'output={EXAMPLE_DIR}/test.py.list.3.out\n')
        f.write('nthreads=6\n')
        f.write('num_vec=20\n')
        f.write('num_jack=10\n')
        f.write('trace=1\n')
        f.write('model=G\n')
        f.write('verbose=1\n')

    with open('config.str.txt', 'w') as f:
        f.write(f'genotype={EXAMPLE_DIR}/test\n')
        f.write(f'phenotype={EXAMPLE_DIR}/test.pheno\n')
        f.write(f'covariate={EXAMPLE_DIR}/test.cov\n')
        f.write(f'environment={EXAMPLE_DIR}/test.env\n')
        f.write(f'annotation={EXAMPLE_DIR}/single.annot\n')
        f.write(f'output={EXAMPLE_DIR}/test.py.str.3.out\n')
        f.write('nthreads=6\n')
        f.write('num_vec=20\n')
        f.write('num_jack=10\n')
        f.write('trace=1\n')
        f.write('model=G\n')
        f.write('verbose=1\n')

    run_genie_mem(
        ['--config', 'config.list.txt']
    )

    run_genie_mem('--config config.str.txt')

    Path('config.list.txt').unlink()
    Path('config.str.txt').unlink()


if __name__ == "__main__":
    test_genie()
    test_genie_multi_pheno()
    test_genie_mem()
