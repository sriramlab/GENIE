from ldscore import PartitionedLDScore
from flex_summ_h2 import FLEX
import argparse
import sys

parser = argparse.ArgumentParser(description='FLEX-summ-h2 (stochastic SE only)')

# --- FLEX run inputs ---
parser.add_argument("--pheno", type=str, default=None,
                    help="Path to phenotype-specific summary statistics. If a directory, all *.sumstat files are used.")
parser.add_argument("--filter-both-sides", action='store_true', default=False,
                    help="When filtering SNPs, remove effects on both trace and yKy. Requires LD scores of the SNPs used in trace.")
parser.add_argument("--ldscores", type=str, default=None,
                    help="Path to LD scores of the reference SNPs (traditional .l2.ldscore.gz or stochastic .ldscore.gz).")
parser.add_argument("--ld-sketch", type=str, default=None,
                    help="Path to stochastic LD matrix sketches (.npy) used for stochastic SE.")
parser.add_argument("--out", type=str, default=None,
                    help="Output prefix to save logs/results or LD scores.")
parser.add_argument("--verbose", action="store_true", default=False,
                    help="Verbose logging.")
parser.add_argument("--annot", type=str, default=None,
                    help="Path to SNP annotation file for partitioned analyses. Optional.")

# --- LD-score/sketch generation mode (if --geno is given) ---
parser.add_argument("--geno", type=str, default=None,
                    help="Prefix to PLINK BED/BIM files to calculate stochastic LD sketches and LD scores.")
parser.add_argument("--nworkers", type=int, default=4,
                    help="Workers for multiprocessing during LD-score/sketch computation. Default 4.")
parser.add_argument("--nvecs", type=int, default=10,
                    help="Number of random vectors for stochastic LD sketches. Default 10.")
parser.add_argument("--step_size", type=int, default=1000,
                    help="Number of SNPs per block for LD computation. Default 1000.")
parser.add_argument("--seed", type=int, default=None,
                    help="Random seed for stochastic LD sketches. Optional.")

if __name__ == '__main__':
    args = parser.parse_args()

    log_msgs = []  # collect log lines to optionally save later

    def logprint(msg):
        """Print and collect log messages."""
        print(msg)
        log_msgs.append(msg + "\n")

    logprint(">>> FLEX-summ-h2 arguments")
    logprint("python3 main.py " + " ".join(sys.argv[1:]))

    # --- Mode 1: LD-score/sketch generation ---
    if args.geno is not None:
        if args.out is None:
            logprint("!!! An output prefix (--out) is required to save LD scores/sketches !!!")
            sys.exit(1)

        gwld = PartitionedLDScore(
            bed_path=args.geno,
            annot_path=args.annot,
            out_path=args.out,
            num_vecs=args.nvecs,
            num_workers=args.nworkers,
            step_size=args.step_size,
            seed=args.seed
        )
        gwld._compute_ldscore()
        # Save log if desired
        if args.out is not None:
            with open(args.out + ".log", "w") as f:
                f.writelines(log_msgs)
        sys.exit(0)

    # --- Mode 2: FLEX run with stochastic SE ---
    if args.pheno is None or args.ld_sketch is None or args.ldscores is None:
        logprint("!!! For a FLEX run, please provide --pheno, --ld-sketch, and --ldscores !!!")
        sys.exit(1)

    if args.out is None:
        logprint(">>> Warning: --out not provided; logs/results may not be saved to disk.")

    # Save log file at the end if requested
    if args.out is not None:
        with open(args.out + ".log", "w") as f:
            f.writelines(log_msgs)

    flex = FLEX(
        pheno_path=args.pheno,
        out=args.out,
        verbose=args.verbose,
        filter_both=args.filter_both_sides,
        ldscores=args.ldscores,
        XtXz=args.ld_sketch,
        annot=args.annot
    )

    flex._run()

