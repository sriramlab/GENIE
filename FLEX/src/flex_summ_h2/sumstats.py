"""
Read the phenotype-specific summary statistics (generally PLINK).
Standard LDSC-like format (.sumstat) works.

Required columns:
    - SNP ID: one of [ID, id, snp, SNP]
    - Sample size: one of [N, n]
    - Z-score: one of [Z, z]

Reference: https://github.com/sriramlab/SUMRHE/blob/main/src/sumstats.py
"""

import numpy as np
import pandas as pd


def _replace_None(li: list):
    """
    Replace None elements in the list with 0.0 (in place) and return the list.
    """
    for i, val in enumerate(li):
        if val is None:
            li[i] = 0.0
    return li


def _partition_bin_non_overlapping(jn_values: np.ndarray, jn_annot: np.ndarray, nbins: int):
    """
    Partition a 1D array of values by a one-hot annotation matrix (rows = SNPs, cols = bins).
    Assumes each SNP belongs to exactly one bin.
    Returns (list_of_bins, snp_counts_per_bin).
    """
    partitions = {i: [] for i in range(nbins)}
    for z, row in zip(jn_values, jn_annot):
        bin_idx = int(np.argmax(row))
        partitions[bin_idx].append(z)
    snp_cnts = [len(partitions[i]) - sum(1 for s in partitions[i] if s is None) for i in range(nbins)]
    return [_replace_None(partitions[i]) for i in range(nbins)], snp_cnts

def _parse_column(df, letters, min_index=3):
    """
    Select a single column from `df` whose name (from column index `min_index` onward)
    contains any of the substrings in `letters`. Returns its values.

    Raises if none or multiple matches are found.
    """
    matching = [col for col in df.columns[min_index:] if any(letter in col for letter in letters)]
    if len(matching) == 0:
        raise ValueError(f"No column containing any of {letters} found starting from column {min_index}.")
    if len(matching) > 1:
        raise ValueError(f"Multiple matching columns found: {matching}. Expected exactly one.")
    return df[matching[0]].values

class Sumstats:
    """
    Holds a single phenotype's summary stats at a time (to save memory).

    Exposes:
      - zscores, zscores_bin, zscores_blk
      - nsamp, nsnps, nsnps_blk, nsnps_bin
      - rhs (for normal equation)
      - removesnps (indices to remove, if chi-sq filtering applied)
    """

    def __init__(self, nblks=100, snplist=None, chisq_threshold=0, both_side=False, annot=None, nbins=1):
        self.nblks = nblks
        self.nbins = nbins
        self.snplist = snplist
        if snplist is not None:
            self.nsnps_trace = len(snplist)

        self.annot = annot
        self.zscores = []
        self.zscores_bin = []   # list per bin
        self.zscores_blk = []   # [blk][bin] -> array
        self.nsamp = 0
        self.nsnps = 0
        self.nsnps_blk = None   # (#blks, #bins)
        self.nsnps_bin = None   # (#bins,)
        self.chisq_threshold = chisq_threshold
        self.both_side = both_side
        self.matched_snps = None
        if self.both_side and snplist is None:
            raise ValueError("SNP list (.bim/.snplist) must be provided to use both-side SNP filtering.")
        self.name = None
        self.removesnps = None
        self.snpids = []

    def _filter_snps(self):
        """
        Remove SNPs with chi-sq above threshold by zeroing their z-scores and
        decrementing nsnps. Return the list of matched snp indices (if available)
        that were removed.
        """
        chisq = np.square(self.zscores)
        removesnps = []
        for i in range(self.nsnps):
            if chisq[i] > self.chisq_threshold:
                self.zscores[i] = 0.0
                self.nsnps -= 1
                if self.snplist is not None and self.matched_snps is not None:
                    if self.matched_snps[i] is not None:
                        removesnps.append(self.matched_snps[i])
        return removesnps

    def _read_sumstats(self, path, name):
        """Read one phenotype's summary statistics."""
        self.name = name
        sumdf = pd.read_csv(path, sep=r"\s+")
        Nmiss = _parse_column(sumdf, ['N', 'n'], 3)
        zscores = _parse_column(sumdf, ['Z', 'z'], 3)
        self.snpids = _parse_column(sumdf, ['ID', 'id', 'snp', 'SNP'], 0)

        nsamp = float(np.max(Nmiss))
        zscores = np.asarray(zscores, dtype=float) * np.sqrt(np.asarray(Nmiss, dtype=float) / nsamp)

        self.zscores = zscores
        self.nsamp = nsamp
        self.nsnps = len(zscores)
        return zscores

    def _match_snps(self):
        """
        Match SNPs between trace SNP list and summary stats, then:
          - split into blocks
          - partition by annotation within each block
        Produces: zscores_blk (#blks x #bins), nsnps_blk (#blks x #bins),
                  zscores_bin (list per bin), nsnps_bin (per bin)
        """
        matched_zscores = []
        nsnps_blk = np.zeros((self.nblks, self.nbins), dtype=float)

        if self.snplist is None:
            # Use all SNPs in sumstats, contiguous block split
            for i in range(self.nblks):
                blk_size = self.nsnps // self.nblks
                if i < self.nblks - 1:
                    sl = slice(blk_size * i, blk_size * (i + 1))
                else:
                    sl = slice(blk_size * i, None)
                blk_zscores = np.array(self.zscores[sl])
                blk_annot = self.annot[sl]
                part, ns_part = _partition_bin_non_overlapping(blk_zscores, blk_annot, self.nbins)
                matched_zscores.append(part)
                nsnps_blk[i] = ns_part
        else:
            # Map SNP IDs to z-scores; then take blocks over snplist order
            zscore_dict = dict(zip(self.snpids, self.zscores))
            # record matched indices of sumstats to snplist order (for potential filtering bookkeeping)
            self.matched_snps = np.arange(len(self.snplist))

            for i in range(self.nblks):
                blk_size = len(self.snplist) // self.nblks
                if i < self.nblks - 1:
                    blk_ids = self.snplist[blk_size * i: blk_size * (i + 1)]
                    blk_annot = self.annot[blk_size * i: blk_size * (i + 1)]
                else:
                    blk_ids = self.snplist[blk_size * i:]
                    blk_annot = self.annot[blk_size * i:]
                blk_zscores = np.array([zscore_dict.get(pid) for pid in blk_ids], dtype=float)
                part, ns_part = _partition_bin_non_overlapping(blk_zscores, blk_annot, self.nbins)
                matched_zscores.append(part)
                nsnps_blk[i] = ns_part

        self.nsnps_blk = nsnps_blk
        self.zscores_blk = matched_zscores
        self.zscores_bin, self.nsnps_bin = _partition_bin_non_overlapping(self.zscores, self.annot, self.nbins)

        # Optional chi-sq filtering
        if self.chisq_threshold is not None:
            self.removesnps = self._filter_snps()

    def _calc_RHS(self):
        """Compute RHS of the normal equation for all blocks and bins."""
        self.rhs = np.full((self.nblks + 1, self.nbins + 1), self.nsamp, dtype=float)

        for i in range(self.nbins):
            total_zTz = float(np.dot(self.zscores_bin[i], self.zscores_bin[i]))
            for j in range(self.nblks + 1):
                if j < self.nblks:
                    blk_zTz = float(np.dot(self.zscores_blk[j][i], self.zscores_blk[j][i]))
                    # if the block covers the whole bin, RHS bin entry is 0
                    if self.nsnps_bin[i] == self.nsnps_blk[j][i]:
                        self.rhs[j, i] = 0.0
                    else:
                        self.rhs[j, i] = (total_zTz - blk_zTz) * self.nsamp / (self.nsnps_bin[i] - self.nsnps_blk[j][i])
                else:
                    self.rhs[j, i] = total_zTz * self.nsamp / self.nsnps_bin[i]

    def _process(self, path, name):
        """High-level: read → match/partition → compute RHS. Returns indices to remove (if any)."""
        self._read_sumstats(path, name)
        self._match_snps()
        self._calc_RHS()
        return self.removesnps
