"""Out-of-distribution gates: Hotelling T² inside the PC subspace, Q (SPE) outside it."""

from dataclasses import dataclass
from typing import List

import numpy as np
from scipy import stats
from sklearn.covariance import LedoitWolf, MinCovDet


def hotelling_threshold(n: int, p: int, alpha: float = 0.01) -> float:
    """T² limit for a NEW observation given a reference of n samples in p dimensions."""
    if n <= p:
        return float("inf")
    return p * (n + 1) * (n - 1) / (n * (n - p)) * stats.f.ppf(1 - alpha, p, n - p)


@dataclass
class ClusterStats:
    location: np.ndarray
    precision: np.ndarray
    n: int
    method: str
    t2_threshold: float


def cluster_stats(P: np.ndarray, labels: np.ndarray, K: int, alpha: float, seed: int = 1729) -> List[ClusterStats]:
    out = []
    p = P.shape[1]
    for j in range(K):
        X = P[labels == j]
        n = X.shape[0]
        if n >= 5 * p:
            est = MinCovDet(random_state=seed).fit(X)
            method = "MCD"
        else:
            est = LedoitWolf().fit(X)
            method = "LedoitWolf"
        out.append(ClusterStats(est.location_, est.precision_, n, method, hotelling_threshold(n, p, alpha)))
    return out


def t2(P: np.ndarray, cs: ClusterStats) -> np.ndarray:
    d = np.atleast_2d(P) - cs.location
    return np.einsum("ij,jk,ik->i", d, cs.precision, d)


def status(posterior_max: float, margin: float, t2_ratio: float, q_ratio: float,
           strong_downweight: float, profile) -> str:
    if strong_downweight > profile.strong_downweight_max:
        return "QC_FAIL"
    if t2_ratio > 1 or q_ratio > 1:
        return "OOD"
    if posterior_max < profile.posterior_min or margin < profile.margin_min:
        return "AMBIGUOUS"
    return "ROUTED"
