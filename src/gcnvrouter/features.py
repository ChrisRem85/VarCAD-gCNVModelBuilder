"""Preflight checks, feature selection and normalization for clustering.

Features are autosomal intervals only, so clusters reflect technical batches rather than sex.
Normalization: per-sample median scaling -> log2 -> residual against the reference median
-> robust (floored MAD) scaling -> winsorizing. See docs/clustering_router_plan.md §3.1.
"""

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from .io_counts import COORDS, CountMatrix
from .profiles import Profile

PSEUDOCOUNT = 0.01
AUTOSOMES = {"chr%d" % i for i in range(1, 23)} | {str(i) for i in range(1, 23)}
CHRX = {"chrX", "X"}
CHRY = {"chrY", "Y"}
# GRCh38 PARs, excluded from the sex check
PAR = {"chrX": [(10001, 2781479), (155701383, 156030895)],
       "chrY": [(10001, 2781479), (56887903, 57217415)]}


def autosomal_mask(intervals: pd.DataFrame) -> np.ndarray:
    return intervals["CONTIG"].isin(AUTOSOMES).to_numpy()


def _outside_par(intervals: pd.DataFrame) -> np.ndarray:
    keep = np.ones(len(intervals), bool)
    for contig, spans in PAR.items():
        on = intervals["CONTIG"].to_numpy() == contig
        for s, e in spans:
            keep &= ~(on & (intervals["START"].to_numpy() <= e) & (intervals["END"].to_numpy() >= s))
    return keep


# --------------------------------------------------------------------------------------------
# preflight
# --------------------------------------------------------------------------------------------

def preflight(cm: CountMatrix, profile: Profile) -> pd.DataFrame:
    """Per-sample sanity checks done before any clustering.

    ASSAY_MISMATCH: too many targets with ~zero depth relative to the cohort, which usually
    means a different enrichment design or probe pool, not a batch effect.
    """
    auto = autosomal_mask(cm.intervals)
    covered = auto & (np.median(cm.counts, axis=0) >= profile.min_median_count)
    if covered.sum() == 0:
        raise ValueError("no autosomal interval reaches min_median_count=%s" % profile.min_median_count)
    med = np.median(cm.counts[:, covered], axis=1)
    med[med == 0] = np.nan
    norm = cm.counts / med[:, None]
    zero_frac = np.mean(norm[:, covered] < profile.zero_depth, axis=1)

    non_par = _outside_par(cm.intervals)
    contigs = cm.intervals["CONTIG"].to_numpy()
    x_mask = np.isin(contigs, list(CHRX)) & non_par
    y_mask = np.isin(contigs, list(CHRY)) & non_par
    x_ratio = np.median(norm[:, x_mask], axis=1) if x_mask.any() else np.full(len(cm.samples), np.nan)
    y_ratio = np.median(norm[:, y_mask], axis=1) if y_mask.any() else np.full(len(cm.samples), np.nan)

    sex = np.full(len(cm.samples), "NA", dtype=object)
    sex[(x_ratio < 0.8) & (y_ratio > 0.1)] = "male"
    sex[(x_ratio >= 0.8) & (y_ratio < 0.05)] = "female"

    status = np.where(np.isnan(med), "NO_COVERAGE",
                      np.where(zero_frac > profile.zero_frac_max, "ASSAY_MISMATCH", "OK"))
    return pd.DataFrame({
        "sample": cm.samples,
        "total_counts": cm.counts.sum(axis=1),
        "median_count": med,
        "dropout_frac": zero_frac,
        "chrX_ratio": x_ratio,
        "chrY_ratio": y_ratio,
        "inferred_sex": sex,
        "preflight": status,
    })


# --------------------------------------------------------------------------------------------
# feature selection
# --------------------------------------------------------------------------------------------

def read_annotated_intervals(path: str, intervals: pd.DataFrame) -> pd.DataFrame:
    """GATK AnnotateIntervals TSV, re-indexed onto the count intervals (NaN where absent)."""
    n_header = 0
    with open(path) as fh:
        for line in fh:
            if not line.startswith("@"):
                break
            n_header += 1
    ann = pd.read_csv(path, sep="\t", skiprows=n_header, dtype={"CONTIG": str})
    return intervals.merge(ann, on=COORDS, how="left")


def read_bed(path: str) -> pd.DataFrame:
    bed = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[0, 1, 2],
                      names=["contig", "start0", "end"], dtype={"contig": str})
    return bed[~bed["contig"].str.startswith(("track", "browser"))]


def bed_overlap_mask(intervals: pd.DataFrame, bed: pd.DataFrame) -> np.ndarray:
    """True for intervals overlapping any BED record."""
    hit = np.zeros(len(intervals), bool)
    contigs = intervals["CONTIG"].to_numpy()
    starts = intervals["START"].to_numpy()
    ends = intervals["END"].to_numpy()
    for contig, grp in bed.groupby("contig"):
        on = np.flatnonzero(contigs == contig)
        if on.size == 0:
            continue
        g = grp.sort_values("start0")
        bs = g["start0"].to_numpy() + 1
        be = g["end"].to_numpy()
        # merge overlapping records so a binary search is exact
        ms, me = [bs[0]], [be[0]]
        for s, e in zip(bs[1:], be[1:]):
            if s <= me[-1] + 1:
                me[-1] = max(me[-1], e)
            else:
                ms.append(s)
                me.append(e)
        ms, me = np.array(ms), np.array(me)
        j = np.searchsorted(ms, ends[on], side="right") - 1
        ok = j >= 0
        hit[on[ok]] = me[j[ok]] >= starts[on[ok]]
    return hit


def select_features(cm: CountMatrix, fit_mask: np.ndarray, profile: Profile,
                    annotated: Optional[pd.DataFrame] = None,
                    mask_beds: Sequence[pd.DataFrame] = (), seed: int = 1729) -> np.ndarray:
    """Indices (into cm.intervals) of clustering features, in genomic file order."""
    elig = autosomal_mask(cm.intervals)
    elig &= np.median(cm.counts[fit_mask], axis=0) >= profile.min_median_count
    gc = None
    if annotated is not None:
        if "MAPPABILITY" in annotated:
            elig &= annotated["MAPPABILITY"].fillna(0).to_numpy() >= 0.9
        if "SEGMENTAL_DUPLICATION_CONTENT" in annotated:
            elig &= annotated["SEGMENTAL_DUPLICATION_CONTENT"].fillna(1).to_numpy() == 0
        if "GC_CONTENT" in annotated:
            gc = annotated["GC_CONTENT"].to_numpy()
            elig &= (gc >= 0.25) & (gc <= 0.75)
    for bed in mask_beds:
        elig &= ~bed_overlap_mask(cm.intervals, bed)
    idx = np.flatnonzero(elig)
    if idx.size == 0:
        raise ValueError("no eligible features left after filtering")
    if profile.max_features is not None and idx.size > profile.max_features:
        idx = _stratified_subsample(idx, cm.intervals, gc, profile.max_features, seed)
    return idx


def _stratified_subsample(idx, intervals, gc, n_target, seed):
    rng = np.random.default_rng(seed)
    strata = intervals["CONTIG"].to_numpy()[idx].astype(object)
    if gc is not None:
        dec = np.clip((gc[idx] * 10).astype(int), 0, 9).astype(str)
        strata = strata + ":" + dec
    keep = []
    frac = n_target / idx.size
    for s in pd.unique(strata):
        members = idx[strata == s]
        n = int(round(members.size * frac))
        if n:
            keep.append(rng.choice(members, size=n, replace=False))
    return np.sort(np.concatenate(keep))


# --------------------------------------------------------------------------------------------
# normalization
# --------------------------------------------------------------------------------------------

def log_normalize(counts: np.ndarray, feature_idx: np.ndarray) -> np.ndarray:
    """log2(count / per-sample median over features + pseudocount). Works on 1 or many samples."""
    c = np.atleast_2d(counts)[:, feature_idx]
    med = np.median(c, axis=1, keepdims=True)
    med[med == 0] = np.nan
    return np.log2(c / med + PSEUDOCOUNT)


@dataclass
class FeatureModel:
    feature_idx: np.ndarray
    ref_median: np.ndarray
    scale: np.ndarray
    winsor_z: float

    @classmethod
    def fit(cls, X: np.ndarray, feature_idx: np.ndarray, winsor_z: float) -> "FeatureModel":
        ref = np.median(X, axis=0)
        mad = 1.4826 * np.median(np.abs(X - ref), axis=0)
        floor = max(np.quantile(mad, 0.25), 1e-6)
        return cls(feature_idx, ref, np.maximum(mad, floor), winsor_z)

    def residual(self, X: np.ndarray) -> np.ndarray:
        return np.atleast_2d(X) - self.ref_median

    def transform(self, X: np.ndarray) -> np.ndarray:
        return np.clip(self.residual(X) / self.scale, -self.winsor_z, self.winsor_z)


def gc_fingerprint(residual: np.ndarray, gc: np.ndarray, n_bins: int = 20):
    """Median residual log-depth per GC bin for each sample (rows). Returns (bin centres, matrix)."""
    edges = np.linspace(0.25, 0.75, n_bins + 1)
    centres = (edges[:-1] + edges[1:]) / 2
    which = np.clip(np.digitize(gc, edges) - 1, 0, n_bins - 1)
    out = np.full((residual.shape[0], n_bins), np.nan)
    for b in range(n_bins):
        sel = which == b
        if sel.sum() >= 5:
            out[:, b] = np.median(residual[:, sel], axis=1)
    return centres, out
