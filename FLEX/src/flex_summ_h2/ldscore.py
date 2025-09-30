"""
Stochastically estimate (partitioned) LD scores.

Reference: https://github.com/sriramlab/SUMRHE/blob/main/src/gw_ldscore.py
"""

import numpy as np
import pandas as pd
from bed_reader import open_bed
import multiprocessing as mp
from tqdm import tqdm
import scipy
import gc

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


class PartitionedLDScore:
    def __init__(self, bed_path, annot_path, out_path, num_vecs=10, num_workers=4, step_size=1000, seed=None):
        self.G = open_bed(bed_path + ".bed")
        self.nsamp, self.nsnps = self.G.shape
        self.nvecs = num_vecs
        self.nworkers = num_workers
        self.step_size = step_size
        self.snp_idx = np.arange(self.nsnps)
        self.outpath = out_path
        self._read_bim(bed_path + ".bim")
        self._read_annot(annot_path)
        self.root_seed = seed
        # will be filled in _compute_ldscore
        self.Z_indiv = None  # (N, nvecs) per-individual random vectors for the noise sketch

    def _compute_Xz_blk(self, blk_idxs):
        """
        Compute X_k z for one block: returns (nbins, N, nvecs).
        Reads the genotype block once, then slices per-annotation bin.
        """
        j, blk_start, blk_end, idxs = blk_idxs  # j is block index; idxs: list of SNP index arrays per bin
        nsnps = sum(len(binidx) for binidx in idxs)
        Xz = np.zeros((self.nbins, self.nsamp, self.nvecs), dtype=float)

        # different RNG per block (deterministic if root_seed is set)
        rng = np.random.default_rng() if self.root_seed is None else np.random.default_rng([j, self.root_seed])
        Zs = rng.standard_normal(size=(nsnps, self.nvecs))  # SNP-space random vectors

        geno = self.G.read(index=np.s_[:, blk_start:blk_end])
        means = np.nanmean(geno, axis=0)
        stds = np.nanstd(geno, axis=0)
        geno = (geno - means) / stds
        geno[np.isnan(geno)] = 0

        geno = np.array(geno, order='F')
        Zs = np.array(Zs, order='F')

        # Multiply per bin in this block
        offset = 0
        for k, binidx in enumerate(idxs):
            m = len(binidx)
            if m == 0:
                continue
            Zs_k = Zs[offset:offset + m, :]
            offset += m
            # (N x m) @ (m x nvecs) -> (N x nvecs)
            Xz[k, :, :] = scipy.linalg.blas.sgemm(1.0, geno[:, binidx], Zs_k)
        return Xz

    def _compute_XtXz_blk(self, blk_idx):
        """
        For a genotype block, multiply with precomputed Xz to get XtX_k z (first nbins columns),
        and compute the noise column as X^T z using the shared per-individual Z_indiv.

        Returns (blk_start, blk_end, XtXz_blk) with shape (block_len, nbins+1, nvecs).
        """
        blk_start, blk_end = blk_idx

        geno = self.G.read(index=np.s_[:, blk_start:blk_end])
        means = np.nanmean(geno, axis=0)
        stds = np.nanstd(geno, axis=0)
        geno = (geno - means) / stds
        geno[np.isnan(geno)] = 0

        geno_t = np.array(geno.T, order='F')  # (block_len, N)

        # allocate with +1 column for the noise sketch
        XtXz_blk = np.zeros((blk_end - blk_start, self.nbins + 1, self.nvecs), dtype=float)

        # structural bins: (block_len, nvecs) = (block_len, N) @ (N, nvecs)
        for k in range(self.nbins):
            XtXz_blk[:, k, :] = scipy.linalg.blas.sgemm(1.0, geno_t, self.Xz[k])

        # noise column: X^T z  (z is per-individual random; same across genome/blocks)
        # (block_len, nvecs) = (block_len, N) @ (N, nvecs)
        XtXz_blk[:, self.nbins, :] = scipy.linalg.blas.sgemm(1.0, geno_t, self.Z_indiv)

        return (blk_start, blk_end, XtXz_blk)

    def _read_annot(self, annot_path):
        """
        Read annotation (optional header). If None, use a single-bin (all ones).
        """
        if annot_path is None:
            self.l2cols = None
            self.annot = np.ones((self.nsnps, 1), dtype=int)
        else:
            self.l2cols, self.annot = _read_with_optional_header(annot_path)
            if self.annot.ndim == 1:
                self.annot = self.annot.reshape(-1, 1)

        if self.annot.shape[0] != self.nsnps:
            raise ValueError(
                f"Number of SNPs in annotation ({self.annot.shape[0]}) "
                f"does not match genotype ({self.nsnps})."
            )

        self.nbins = self.annot.shape[1]
        if self.l2cols is None:
            self.l2cols = [f"L2_{i}" for i in range(self.nbins)]
        else:
            self.l2cols = [f"{c}L2" for c in self.l2cols]
        self.nsnps_bin = self.annot.sum(axis=0)

    def _read_bim(self, bim_path):
        """
        Read .bim for SNP metadata if available; otherwise keep snplist=None.
        """
        try:
            self.snplist = pd.read_csv(bim_path, header=None, sep='\t')
        except FileNotFoundError:
            self.snplist = None
            return

        if len(self.snplist) != self.nsnps:
            raise ValueError(
                f"Number of SNPs in .bed ({self.nsnps}) does not match .bim ({len(self.snplist)})."
            )

    def _partition_index(self, snpidx, annot):
        """Partition SNP indices by annotation for a given block."""
        return [snpidx[annot[:, c] == 1] for c in range(self.nbins)]

    def _compute_ldscore(self):
        """
        Compute stochastic LD sketches and LD scores.

        Sketch arrays saved to disk:
          XtXz (RAW) shape = (nsnps, nbins + 1, nvecs)
            - columns 0..nbins-1: X^T X_k z_b  (structural)
            - column nbins:       X^T z_b      (noise)

        LD scores are computed (and saved) only for the structural bins, from the
        normalized sketches XtXz / N.
        """
        # Prepare a shared per-individual random matrix for the noise sketch
        rng = np.random.default_rng(self.root_seed) if self.root_seed is not None else np.random.default_rng()
        self.Z_indiv = rng.standard_normal(size=(self.nsamp, self.nvecs))  # (N, nvecs)

        # Block definitions
        self.nblks = len(np.arange(self.nsnps)[::self.step_size])
        Xz_input = []
        XtXz_input = []
        for j in range(self.nblks):
            idx_start = self.step_size * j
            idx_end = self.nsnps if j == self.nblks - 1 else self.step_size * (j + 1)
            annot_blk = self.annot[idx_start:idx_end]
            Xz_input.append((j, idx_start, idx_end,
                             self._partition_index(np.arange(len(annot_blk)), annot_blk)))
            XtXz_input.append((idx_start, idx_end))

        # Aggregate Xz over all blocks (structural bins only)
        self.Xz = np.zeros((self.nbins, self.nsamp, self.nvecs), dtype=float)
        with mp.Pool(self.nworkers) as pool, tqdm(total=self.nblks, desc="Calculating Xz") as pbar:
            for result in pool.imap_unordered(self._compute_Xz_blk, Xz_input):
                self.Xz += result
                gc.collect()
                pbar.update()
            pool.close()
            pool.join()

        # Compute XtXz by re-reading genotype blocks (now with noise column as well)
        self.XtXz = np.zeros((self.nsnps, self.nbins + 1, self.nvecs), dtype=float)  # RAW sketches
        with mp.Pool(self.nworkers) as pool, tqdm(total=self.nblks, desc="Calculating XtXz (+ noise)") as pbar:
            for result in pool.imap_unordered(self._compute_XtXz_blk, XtXz_input):
                idx_start, idx_end, XtXz_blk = result
                # XtXz_blk is (block_len, nbins+1, nvecs) -> fits the destination
                self.XtXz[idx_start:idx_end, :, :] = XtXz_blk
                gc.collect()
                pbar.update()
            pool.close()
            pool.join()

        # Normalized sketches (per LDS score convention)
        self.XtXz_norm = self.XtXz / self.nsamp

        # Convert sketches to LD scores (per structural bin only)
        ldscore_struct = self.nsamp / (self.nsamp + 1.0) * (
            np.square(self.XtXz_norm[:, :self.nbins, :]).mean(axis=2) - self.nsnps_bin / self.nsamp
        )

        # Prepare output table and files
        snpcols = ['CHR', 'SNP', 'BP']
        if self.snplist is None:
            snpdf = pd.DataFrame(np.nan * np.ones((self.nsnps, 3)), columns=snpcols)
        else:
            snpdf = self.snplist.iloc[:, :3].copy()
            snpdf.columns = snpcols

        ld_df = pd.DataFrame(ldscore_struct, columns=self.l2cols)
        ld_df = pd.concat([snpdf, ld_df], axis=1)

        # Save RAW sketches (with noise column) and LD scores
        np.save(f"{self.outpath}.sketches.npy", self.XtXz_norm)
        ld_df.to_csv(f"{self.outpath}.ldscore.gz", index=False, compression="gzip", sep="\t", float_format="%.3f")
