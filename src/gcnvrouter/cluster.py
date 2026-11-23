"""Gaussian-mixture sweep with a minimum-cluster-size constraint and subsampling stability."""

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture

COV_TYPES = ("tied", "diag")


def _gmm(K: int, cov: str, P: np.ndarray, seed: int, n_init: int = 10) -> GaussianMixture:
    reg = 1e-3 * float(np.median(P.var(axis=0))) + 1e-9
    return GaussianMixture(n_components=K, covariance_type=cov, n_init=n_init, reg_covar=reg,
                           random_state=seed, max_iter=500).fit(P)


@dataclass
class SweepResult:
    table: pd.DataFrame                              # K, cov, bic, sizes, min_size, valid
    models: Dict[Tuple[int, str], GaussianMixture]


def sweep(P: np.ndarray, k_max: int, n_min: int, seed: int = 1729,
          cov_types=COV_TYPES) -> SweepResult:
    n = P.shape[0]
    rows, models = [], {}
    for K in range(1, max(1, min(k_max, n // max(n_min, 1))) + 1):
        for cov in cov_types:
            g = _gmm(K, cov, P, seed)
            sizes = np.bincount(g.predict(P), minlength=K)
            rows.append({"K": K, "cov": cov, "bic": g.bic(P),
                         "sizes": ",".join(str(int(s)) for s in sorted(sizes, reverse=True)),
                         "min_size": int(sizes.min()), "valid": bool(sizes.min() >= n_min)})
            models[(K, cov)] = g
    return SweepResult(pd.DataFrame(rows), models)


SUBSAMPLE_FRACTION = 0.8


def bootstrap_stability(P: np.ndarray, K: int, cov: str, ref_labels: np.ndarray,
                        n_boot: int = 100, seed: int = 1729):
    """Consensus-clustering stability: refit on 80% subsamples (without replacement), predict all
    samples, and compare with the full fit (ARI). Also returns a co-clustering frequency matrix.

    Subsampling rather than the classic bootstrap: with small cohorts a bootstrap resample holds
    only ~63% unique samples, so small but real clusters often cannot be refitted at all.
    """
    n = P.shape[0]
    co = np.zeros((n, n))
    if K == 1:
        return np.ones(n_boot), np.ones((n, n))
    rng = np.random.default_rng(seed)
    aris = []
    m = max(K + 1, int(round(SUBSAMPLE_FRACTION * n)))
    for b in range(n_boot):
        idx = rng.choice(n, size=m, replace=False)
        try:
            g = _gmm(K, cov, P[idx], seed + b + 1, n_init=3)
        except ValueError:            # degenerate resample
            continue
        lab = g.predict(P)
        aris.append(adjusted_rand_score(ref_labels, lab))
        co += lab[:, None] == lab[None, :]
    return np.array(aris), co / max(len(aris), 1)


@dataclass
class ClusterResult:
    K: int
    cov: str
    labels: np.ndarray           # 0..K-1, ordered by size (desc) then mean PC1
    posteriors: np.ndarray       # n x K, columns in label order
    model: GaussianMixture
    order: np.ndarray            # model component index for each label
    aris: np.ndarray
    coassign: np.ndarray
    sweep: pd.DataFrame
    notes: List[str]


def choose_clusters(P: np.ndarray, k_max: int, n_min: int, ari_min: float,
                    n_boot: int = 100, seed: int = 1729) -> ClusterResult:
    sw = sweep(P, k_max, n_min, seed)
    table = sw.table.copy()
    table["chosen"] = False
    table["median_ari"] = np.nan
    notes = []
    # every valid (K, covariance) solution in BIC order; the first stable one wins
    candidates = table[table["valid"]].sort_values(["bic", "K"])
    chosen = None
    for _, row in candidates.iterrows():
        K, cov = int(row["K"]), row["cov"]
        g = sw.models[(K, cov)]
        labels = g.predict(P)
        aris, co = bootstrap_stability(P, K, cov, labels, n_boot, seed)
        med = float(np.median(aris)) if aris.size else 0.0
        table.loc[(table["K"] == K) & (table["cov"] == cov), "median_ari"] = med
        if K == 1 or med >= ari_min:
            chosen = (K, cov, g, aris, co)
            break
        notes.append("K=%d (%s) rejected: subsample median ARI %.2f < %.2f" % (K, cov, med, ari_min))
    if chosen is None:   # nothing valid and stable -> one cluster
        g = sw.models.get((1, "tied")) or _gmm(1, "tied", P, seed)
        chosen = (1, "tied", g, np.ones(n_boot), np.ones((P.shape[0], P.shape[0])))
        notes.append("no valid stable solution; falling back to K=1")
    K, cov, g, aris, co = chosen
    table.loc[(table["K"] == K) & (table["cov"] == cov), "chosen"] = True

    raw = g.predict(P)
    post = g.predict_proba(P)
    sizes = np.bincount(raw, minlength=K)
    pc1 = np.array([P[raw == j, 0].mean() if sizes[j] else np.inf for j in range(K)])
    order = np.lexsort((pc1, -sizes))            # label 0 = largest cluster
    remap = np.empty(K, int)
    remap[order] = np.arange(K)
    return ClusterResult(K, cov, remap[raw], post[:, order], g, order, aris, co, table, notes)


def label_name(i: int) -> str:
    return "C%d" % (i + 1)
