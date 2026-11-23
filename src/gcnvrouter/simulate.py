"""Synthetic CollectReadCounts-like data with batch effects, for tests and demos."""

from typing import Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .io_counts import CountMatrix


def simulate(batch_sizes: Sequence[int] = (40, 40, 40), n_intervals: int = 1500, seed: int = 0,
             batch_effect_sd: float = 0.15, gc_effect: float = 1.0, dispersion: float = 0.01,
             depth_range: Tuple[float, float] = (1.5e6, 3e6), dup_chromosome: Optional[str] = None,
             dup_fraction: float = 0.0, intervals: Optional[pd.DataFrame] = None):
    """Return (CountMatrix, batch labels, sex labels, GC content per interval).

    Each batch gets its own GC-bias curve and its own per-interval capture deviations, which is
    how wet-lab batches look in real panels/exomes. Sex only affects chrX/chrY.
    `dup_chromosome`: for `dup_fraction` of samples, all intervals on that chromosome get 1.5x
    depth (a whole-chromosome event shared by a subgroup -> biology that should be flagged).
    """
    rng = np.random.default_rng(seed)
    if intervals is None:
        contigs = np.array(["chr%d" % (1 + i % 22) for i in range(n_intervals - 24)] + ["chrX"] * 20 + ["chrY"] * 4)
        contigs = np.array(sorted(contigs, key=lambda c: (c[3:].isdigit() is False, int(c[3:]) if c[3:].isdigit() else c)))
        starts = np.zeros(len(contigs), dtype=np.int64)
        for c in pd.unique(contigs):
            sel = contigs == c
            starts[sel] = 10_000_000 + 10_000 * np.arange(sel.sum())
        intervals = pd.DataFrame({"CONTIG": contigs, "START": starts, "END": starts + 299})
    p = len(intervals)
    gc = rng.uniform(0.3, 0.7, p)
    efficiency = rng.lognormal(0, 0.4, p)
    contigs = intervals["CONTIG"].to_numpy()
    rows, batches, sexes = [], [], []
    for b, n in enumerate(batch_sizes):
        a, q = rng.normal(0, gc_effect), rng.normal(0, 3 * gc_effect)
        capture = rng.normal(0, batch_effect_sd, p)
        log_batch = a * (gc - 0.5) + q * (gc - 0.5) ** 2 + capture
        for _ in range(n):
            sex = rng.choice(["female", "male"])
            dose = np.ones(p)
            if sex == "male":
                dose[contigs == "chrX"] = 0.5
                dose[contigs == "chrY"] = 0.5
            else:
                dose[contigs == "chrY"] = 0.0
            mu = efficiency * np.exp(log_batch + rng.normal(0, 0.03, p)) * dose
            mu = mu / mu.sum() * rng.uniform(*depth_range)
            lam = mu * rng.gamma(1 / dispersion, dispersion, p)
            rows.append(rng.poisson(lam))
            batches.append(b)
            sexes.append(sex)
    counts = np.vstack(rows).astype(float)
    if dup_chromosome:
        n = counts.shape[0]
        pick = rng.choice(n, size=int(round(dup_fraction * n)), replace=False)
        counts[np.ix_(pick, contigs == dup_chromosome)] *= 1.5
    names = ["S%03d" % i for i in range(counts.shape[0])]
    return CountMatrix(names, intervals.reset_index(drop=True), counts), np.array(batches), np.array(sexes), gc


def add_samples(cm: CountMatrix, counts: np.ndarray, prefix: str) -> CountMatrix:
    names = cm.samples + ["%s%03d" % (prefix, i) for i in range(counts.shape[0])]
    return CountMatrix(names, cm.intervals, np.vstack([cm.counts, counts]))
