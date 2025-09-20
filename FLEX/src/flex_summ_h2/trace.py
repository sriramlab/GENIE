"""
Trace: read trace summaries (.tr/.MN) or LD scores (.ldscore.gz) and
compute blockwise traces for FLEX.

Reference: https://github.com/sriramlab/SUMRHE/blob/main/src/trace.py
"""

import numpy as np
import pandas as pd
from os import listdir, path


def _calc_trace_from_ld(ldsum, n, m1, m2):
    """
    Calculate the trace from the sum of LD scores.
    ldsum: sum of LD scores for (k,l) block
    n: sample size
    m1, m2: SNP counts in bins k and l for the (sub)sample
    """
    if (m1 == 0) or (m2 == 0):
        return n
    return ldsum * (n ** 2) / (m1 * m2) + n


def _read_multiple_lines(file_path, num_lines, sep=','):
    """
    Yield numpy arrays by reading a file in chunks of `num_lines` rows.
    """
    for val in pd.read_csv(file_path, chunksize=num_lines, sep=sep):
        yield val.to_numpy()


def _read_with_optional_header(file_path):
    """
    Read a whitespace-delimited numeric matrix, optionally with a header line.
    If the first line is non-numeric, it is treated as a header (list of column names).
    Returns (header_or_None, data_ndarray).
    """
    with open(file_path, 'r') as fd:
        first = fd.readline().strip()
        try:
            _ = [float(x) for x in first.split()]
            has_header = False
        except ValueError:
            has_header = True

    if has_header:
        header = first.split()
        data = np.loadtxt(file_path, skiprows=1)
        return header, data
    else:
        data = np.loadtxt(file_path)
        return None, data

class Trace:
    def __init__(self, bimpath=None, sumpath=None, savepath=None,
                 ldscores=None, nblks=100, annot=None, verbose=False):
        self.sumpath = sumpath
        self.savepath = savepath
        self.ldscorespath = ldscores
        self.ldscores = None
        self.sums = []
        self.nblks = nblks
        self.ntrace = 0
        self.K = []
        self.snplist = []
        self.nsamp = []           # per-trace sample sizes
        self.nsnps = 0
        self.nsnps_blk = None
        self.nsnps_bin = None
        self.nsnps_blk_filt = None
        self.verbose = verbose

        if bimpath and bimpath.endswith(".bim"):
            with open(bimpath, "r") as fd:
                for line in fd:
                    self.snplist.append(line.split()[1])
        elif bimpath:  # provided but not a .bim file
            raise ValueError(f"{bimpath} is not a .bim file")

        if self.sumpath is not None:
            self._read_all_trace()
        elif self.ldscorespath is not None:
            self._read_ldscores()

        self._read_annot(annot)

    def _read_annot(self, annot_path):
        if annot_path is None:
            self.annot_header = np.array(['L2'])
            self.annot = np.ones((self.nsnps, 1), dtype=int)
        else:
            self.annot_header, self.annot = _read_with_optional_header(annot_path)
            if self.annot.ndim == 1:
                self.annot = self.annot.reshape(-1, 1)
            if self.annot_header is None:
                self.annot_header = np.array([f"L2_{i}" for i in range(self.annot.shape[1])])

        if hasattr(self, "nbins"):
            expected_bins = self.nbins
        else:
            expected_bins = self.annot.shape[1]
            self.nbins = expected_bins

        if self.annot.shape[1] != expected_bins or self.annot.shape[0] != self.nsnps:
            raise ValueError("Annotation dimensions do not match the loaded trace/LD scores.")

        self.blk_size = self.nsnps // self.nblks if self.nblks else self.nsnps
        self.nsnps_bin = self.annot.sum(axis=0)
        self.nsnps_blk = np.full((self.nblks + 1, self.nbins), self.nsnps_bin, dtype=float)
        for i in range(self.nblks):
            idx_start = self.blk_size * i
            idx_end = self.nsnps if i == self.nblks - 1 else self.blk_size * (i + 1)
            self.nsnps_blk[i] -= self.annot[idx_start:idx_end].sum(axis=0)

    def _read_trace(self, filename, idx):
        with open(filename + ".MN", "r") as fd:
            next(fd)
            nsamp, nsnps, nblks, nbins, K = map(int, fd.readline().split(","))
        self.nsamp.append(nsamp)
        self.K.append(K)

        if idx == 0:
            self.nsnps = nsnps
            self.nblks = nblks
            self.nbins = nbins
        else:
            if nblks != self.nblks:
                raise ValueError(f"Trace {filename} has inconsistent jackknife blocks")
            if nbins != self.nbins:
                raise ValueError(f"Trace {filename} has inconsistent bins")

        sums = np.zeros((self.nblks + 1, self.nbins, self.nbins))
        nsnps_blk = np.zeros((self.nblks + 1, self.nbins))
        for cnt, vals in enumerate(_read_multiple_lines(filename + ".tr", self.nbins)):
            sums[cnt] = vals[:, :-1]
            nsnps_blk[cnt] = vals[:, -1].T

        if idx == 0:
            self.nsnps_blk = nsnps_blk
        else:
            if not np.array_equal(self.nsnps_blk, nsnps_blk):
                raise ValueError(f"Trace {filename} has different annotation SNP counts")

        self.sums.append(sums)
        self.ntrace += 1
        return self.sums

    def _read_all_trace(self):
        trace_files = [self.sumpath]
        if path.isdir(self.sumpath):
            prefix = self.sumpath.rstrip('/') + "/"
            trace_files = sorted([prefix + f.rstrip('.tr') for f in listdir(self.sumpath) if f.endswith('.tr')])

        for i, f in enumerate(trace_files):
            self._read_trace(f, i)

        if self.ntrace > 1:
            if np.std(self.sums, axis=0).sum() == 0.0:
                self.effective_K = int(np.min(self.K))
            else:
                self.effective_K = int(np.sum(self.K))
            self.sums = np.average(self.sums, axis=0, weights=self.nsamp)
            self.nsamp = float(np.mean(self.nsamp))
        else:
            self.effective_K = int(self.K[0])
            self.sums = self.sums[0]
            self.nsamp = float(self.nsamp[0])

        return self.sums

    def _read_ldscores(self):
        df = pd.read_csv(self.ldscorespath, compression='gzip', sep='\s+', index_col=False)
        self.ldscores = df.iloc[:, 3:].to_numpy()
        self.snplist = df['SNP'].to_numpy() if 'SNP' in df.columns else []
        self.nsnps = self.ldscores.shape[0]
        self.nbins = self.ldscores.shape[1]

    def _calc_trace(self, nsample):
        if self.ldscores is not None:
            return self._calc_trace_from_ldscores(nsample)
        elif isinstance(self.sums, np.ndarray) and self.sums.ndim == 3:
            return self._calc_trace_from_sums(nsample)
        else:
            raise ValueError("Trace not initialized: provide .tr/.MN or .ldscore.gz")

    def _calc_trace_from_sums(self, N):
        trace = np.full((self.nblks + 1, self.nbins + 1, self.nbins + 1), N, dtype=float)
        for k in range(self.nbins):
            for l in range(self.nbins):
                for j in range(self.nblks + 1):
                    trace[j, k, l] = _calc_trace_from_ld(
                        self.sums[j, k, l], N,
                        self.nsnps_blk[j, k], self.nsnps_blk[j, l])
        return trace

    def _calc_trace_from_ldscores(self, N):
        trace = np.full((self.nblks + 1, self.nbins + 1, self.nbins + 1), N, dtype=float)
        for k in range(self.nbins):
            for l in range(self.nbins):
                ld_sum = self.ldscores[self.annot[:, k] == 1][:, l].sum()
                for j in range(self.nblks + 1):
                    idx_start = self.blk_size * j
                    idx_end = self.nsnps if j == self.nblks - 1 else self.blk_size * (j + 1)
                    annot_jn = self.annot[idx_start:idx_end]
                    ldscores_jn = self.ldscores[idx_start:idx_end]
                    ld_sum_jn = ld_sum - ldscores_jn[annot_jn[:, k] == 1][:, l].sum()
                    if j == self.nblks:
                        ld_sum_jn = ld_sum
                    trace[j, k, l] = _calc_trace_from_ld(
                        ld_sum_jn, N,
                        self.nsnps_blk[j, k], self.nsnps_blk[j, l])
        return trace

    def _reset(self):
        self.nsnps_blk_filt = self.nsnps_blk
        self.sums_filt = self.sums

    def _filter_snps(self, removelist):
        """
        Placeholder for removing SNPs from trace calculation. No real filtering logic implemented.
        """
        # Add filtering logic if needed later
        return
