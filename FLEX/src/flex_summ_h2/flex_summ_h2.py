from sumstats import Sumstats
from trace import Trace
import time
import datetime

import numpy as np
import os

def _get_time():
    """Return current epoch time (float)."""
    return time.time()


def _get_timestr(current_time):
    """Format a timestamp with local timezone."""
    timezone = datetime.datetime.now(datetime.timezone.utc).astimezone().tzinfo
    return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(current_time)) + " " + str(timezone)

class FLEX:
    """
    FLEX: Method-of-moments variance component estimation with stochastic-LD SEs only.

    Required inputs:
      - pheno_path: file or directory of *.sumstat files
      - (optional) out: output prefix; if provided, logs are written to <out>.log (and also printed)
      - (optional) verbose: bool for extra debug prints
      - (optional) filter_both: bool, filter SNPs on both sides
      - (optional) ldscores: path for LD scores (if needed by Trace)
      - (optional) XtXz: path to stochastic LD sketches .npy
      - (optional) annot: SNP annotation file path

    Internally fixed:
      - chisq_threshold = 10000
      - njack = 100
      - mem = False
      - allsnp = False
      - bim_path = None, sum_path = None, save_path = None
    """

    # ---- Internal defaults ----
    _INTERNAL_CHISQ_THRESHOLD = 10000
    _INTERNAL_NJACK = 100
    _INTERNAL_MEM = False
    _INTERNAL_ALLSNP = False
    _INTERNAL_BIM = None
    _INTERNAL_SUM = None
    _INTERNAL_SAVE = None

    def __init__(self,
                 pheno_path,
                 out=None,
                 verbose=False,
                 filter_both=False,
                 ldscores=None,
                 XtXz=None,
                 annot=None):

        # print helper: stdout always; to file if out provided
        self.out = out
        self._logfh = open(f"{out}.log", "a") if out else None

        def _print(msg):
            print(msg)
            if self._logfh:
                self._logfh.write(msg + "\n")
                self._logfh.flush()
        self._print = _print

        # Basic state
        self.mem = self._INTERNAL_MEM
        self.verbose = verbose
        self.filter_both = filter_both

        # Start time
        self.start_time = _get_time()
        self._print("Analysis started at: " + _get_timestr(self.start_time))

        # Trace (provides traces T) and metadata
        self.tr = Trace(
            bimpath=self._INTERNAL_BIM,
            sumpath=self._INTERNAL_SUM,
            savepath=self._INTERNAL_SAVE,
            ldscores=ldscores,
            nblks=self._INTERNAL_NJACK,
            annot=annot,
            verbose=verbose
        )
        self.snplist = self.tr.snplist
        self.nblks = self.tr.nblks
        self.annot = self.tr.annot
        self.annot_header = self.tr.annot_header
        self.nbins = self.tr.nbins

        # Sumstats
        if self._INTERNAL_ALLSNP:
            self.sums = Sumstats(
                nblks=self.nblks,
                chisq_threshold=self._INTERNAL_CHISQ_THRESHOLD,
                annot=self.annot,
                nbins=self.nbins
            )
        else:
            self.sums = Sumstats(
                nblks=self.nblks,
                snplist=self.snplist,
                chisq_threshold=self._INTERNAL_CHISQ_THRESHOLD,
                annot=self.annot,
                nbins=self.nbins
            )

        # Phenotype file(s)
        if os.path.exists(pheno_path):
            if os.path.isdir(pheno_path):
                self._print("Reading phenotype sumstat files from a directory...")
                phen_files = sorted([f for f in os.listdir(pheno_path) if f.endswith(".sumstat")])
                if len(phen_files) == 0:
                    raise ValueError(f"--pheno path {pheno_path} has no valid phenotype summary files (.sumstat)")
                self.phen_dir = [pheno_path.rstrip("/") + "/" + name for name in phen_files]
            elif os.path.isfile(pheno_path):
                self._print("Reading a single phenotype sumstat file")
                self.phen_dir = [pheno_path]
            else:
                raise ValueError(f"--pheno path {pheno_path} is invalid")
        else:
            raise ValueError(f"--pheno path {pheno_path} is invalid")

        self.npheno = len(self.phen_dir)
        self.nsamp = []  # one per phenotype, filled in _run()
        self.phen_names = [os.path.basename(name)[:-8] for name in self.phen_dir]

        # Output containers:
        # hsums[phenotype, bin or total, (estimate, SE)]
        # bins: 0..nbins-1 are structural bins, index -1 is total (over structural bins)
        self.hsums = np.zeros((self.npheno, self.nbins + 2, 2))

        # LD sketches XtXz (nsnps, nbins, B)
        self.lds = np.load(XtXz) if XtXz is not None else None


    def _run(self):
        """
        Run analysis for all phenotypes using stochastic-LD SE.

        Key points for SE:
        • self.lds has shape (nsnps, nbins_ext, B) where nbins_ext = nbins_genetic + 1
        and the LAST column (index -1) is the **noise sketch** (X^T z_b)/N.
        • Noise row in S uses  z^T X^T β (replicated across columns):
            s_noise = sum_s β_s * ((X^T z_b)/N)_s
        • M_noise = 1.
        • Cov[q] uses  Σ_l (σ_l^2 / M_l) * S[:, l] S[:, l]^T, scaled by (2/B) * N^2 / (M_i M_j).
        """
        if self.lds is None:
            raise ValueError("Cannot run: stochastic LD sketches (XtXz/lds) are not provided.")
        self.start_time = _get_time()

        def _info(msg):  self._print(msg)
        def _debug(msg):
            if self.verbose: self._print(f"[DEBUG] {msg}")
        def _warn(msg):
            if self.verbose: self._print(f"[WARN]  {msg}")

        XtXz_all = self.lds
        nsnps_all, nbins_ext_all, B_all = XtXz_all.shape
        _debug(f"XtXz global shape: nsnps={nsnps_all}, nbins_ext={nbins_ext_all}, B={B_all}")

        for idx in range(self.npheno):
            t0 = _get_time()
            pheno_path = self.phen_dir[idx]
            _info(f"=== Phenotype {idx+1}/{self.npheno}: {self.phen_names[idx]} ===")
            removesnps = self.sums._process(pheno_path, self.phen_names[idx])
            self.tr._reset()
            if removesnps is not None and self.filter_both:
                self.tr._filter_snps(removesnps)

            N = float(self.sums.nsamp)
            self.nsamp.append(N)

            nbins_gen = int(self.nbins)
            XtXz = self.lds  # (nsnps, nbins_ext, B)
            nsnps, nbins_ext, B = XtXz.shape
            if nbins_ext != nbins_gen + 1:
                raise ValueError(f"XtXz second dim {nbins_ext} must equal nbins_gen+1 {nbins_gen+1} (last col = noise).")

            betas = self.sums.zscores
            annot = self.sums.annot
            M = np.asarray(self.sums.nsnps_bin, float)
            _info(f"idx={idx} | nbins_gen={nbins_gen} | nsnps={nsnps} | B={B} | N={N:.0f}")

            # trace and RHS
            T_full = self.tr._calc_trace(N)[self.nblks]
            q_full = self.sums.rhs[self.nblks]

            # theta_hat
            theta_hat, *_ = np.linalg.lstsq(T_full, q_full, rcond=None)
            Ktheta = T_full.shape[1]

            # M_ext
            M_ext = np.zeros(nbins_gen+1)
            M_ext[:nbins_gen] = M
            M_ext[-1] = 1.0  # M_noise=1

            # sigma_l^2/M_l
            if Ktheta == nbins_gen+1:
                sig2_l_ext = np.asarray(theta_hat).ravel()
            elif Ktheta == nbins_gen:
                sig2_l_ext = np.concatenate([np.asarray(theta_hat).ravel(), [0.0]])
            else:
                raise ValueError(f"Unexpected #theta cols: {Ktheta}")

            sig2_l_ext = np.clip(sig2_l_ext, 0.0, None)
            sig2_divM = sig2_l_ext / M_ext


            # precompute masks
            bin_masks = [annot[:, i] == 1 for i in range(nbins_gen)]
            betas_by_bin = [betas[mask] for mask in bin_masks]

            # build Cov[q]
            cov_q = np.zeros((nbins_gen+1, nbins_gen+1))
            for b in range(B):
                XtXz_b = XtXz[:, :, b]
                noise_vec = XtXz_b[:, -1]  # (X^T z_b)/N

                S = np.zeros((nbins_gen+1, nbins_gen+1))
                for i in range(nbins_gen):
                    Xi = XtXz_b[bin_masks[i], :]
                    if Xi.size:
                        S[i, :] = (Xi * betas_by_bin[i][:, None]).sum(axis=0)
                S[-1, :] = (XtXz_b * betas[:, None]).sum(axis=0)

                cov_q += S @ np.diag(sig2_divM) @ S.T

            # scale
            scale = (2.0 / B) * N
            for i in range(nbins_gen+1):
                for j in range(nbins_gen+1):
                    cov_q[i, j] *= scale / (M_ext[i] * M_ext[j])
            cov_q = 0.5*(cov_q+cov_q.T)

            cov_q_used = cov_q if (Ktheta == nbins_gen+1) else cov_q[:nbins_gen,:nbins_gen]

            # Var(theta)
            TtT_inv = np.linalg.pinv(T_full.T @ T_full)
            var_theta = TtT_inv @ (T_full.T @ cov_q_used @ T_full) @ TtT_inv
            var_theta = 0.5*(var_theta+var_theta.T)
            se_theta = np.sqrt(np.clip(np.diag(var_theta), 0.0, None))

            # store per-bin
            self.hsums[idx,:nbins_gen,0] = np.asarray(theta_hat).ravel()[:nbins_gen]
            self.hsums[idx,:nbins_gen,1] = se_theta[:nbins_gen]

            # combined (gene) bins
            w = np.zeros((nbins_gen,1))
            w[1:] = 1.0  # for example, sum bins 1..end (skip flanking)
            h2_gene = float(np.asarray(theta_hat).ravel()[:nbins_gen].dot(w.ravel()))
            var_gene = float(w.T @ var_theta[:nbins_gen,:nbins_gen] @ w)
            self.hsums[idx,-1,0] = h2_gene
            self.hsums[idx,-1,1] = float(np.sqrt(max(var_gene,0.0)))

            self._print("Heritabilities:")
            for i in range(nbins_gen):
                self._print(f"h2_g[{i}] : {self.hsums[idx,i,0]:.6g} SE: {self.hsums[idx,i,1]:.6g}")
            self._print(f"h2_g[gene] : {self.hsums[idx,-1,0]:.6g} SE: {self.hsums[idx,-1,1]:.6g}")
            _debug(f"Phenotype runtime: {_get_time()-t0:.3f}s")

        self.end_time = _get_time()
        self._print("Analysis ended at: " + _get_timestr(self.end_time))
        self._print("run time: " + format(self.end_time - self.start_time, '.3f') + " s")
        if self._logfh: self._logfh.close()
        return self.hsums
