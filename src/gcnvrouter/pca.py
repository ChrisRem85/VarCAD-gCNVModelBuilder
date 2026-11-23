"""PCA with parallel analysis, a loading-locality check and robust (IRLS) projection."""

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd


@dataclass
class PCAModel:
    mean: np.ndarray            # p
    components: np.ndarray      # k x p, orthonormal rows, sign-canonicalized
    eigvals: np.ndarray         # all eigenvalues (variance per PC)
    pa_threshold: np.ndarray    # parallel-analysis 95th percentile per PC (NaN if not run)
    k: int
    k_pa: int                   # PCs above the parallel-analysis threshold (0 = no structure)

    @property
    def explained(self) -> np.ndarray:
        return self.eigvals / self.eigvals.sum()

    def scores(self, Z: np.ndarray) -> np.ndarray:
        return (np.atleast_2d(Z) - self.mean) @ self.components.T

    def q_plain(self, Z: np.ndarray) -> np.ndarray:
        zc = np.atleast_2d(Z) - self.mean
        e = zc - (zc @ self.components.T) @ self.components
        return np.sum(e ** 2, axis=1)


def _eig(Zc: np.ndarray):
    _, s, vt = np.linalg.svd(Zc, full_matrices=False)
    return s ** 2 / max(Zc.shape[0] - 1, 1), vt


def fit_pca(Z: np.ndarray, k_max: int, n_perm: int = 50, seed: int = 1729,
            k_fixed: Optional[int] = None) -> PCAModel:
    mean = Z.mean(axis=0)
    Zc = Z - mean
    eigvals, vt = _eig(Zc)
    n_pc = len(eigvals)
    pa = np.full(n_pc, np.nan)
    k_pa = 0
    if n_perm and k_fixed is None:
        rng = np.random.default_rng(seed)
        perm_eigs = np.empty((n_perm, n_pc))
        for i in range(n_perm):
            Zp = np.column_stack([rng.permutation(Zc[:, j]) for j in range(Zc.shape[1])])
            perm_eigs[i] = _eig(Zp - Zp.mean(axis=0))[0]
        pa = np.quantile(perm_eigs, 0.95, axis=0)
        above = eigvals > pa
        k_pa = int(np.argmin(above)) if not above.all() else n_pc
    if k_fixed is not None:
        k = int(k_fixed)
    else:
        k = int(np.clip(k_pa, 2, k_max))   # >=2 so there is always something to look at
    k = min(k, n_pc)
    comps = vt[:k].copy()
    for i in range(k):                       # sign convention: largest |loading| positive
        if comps[i, np.argmax(np.abs(comps[i]))] < 0:
            comps[i] = -comps[i]
    return PCAModel(mean, comps, eigvals, pa, k, k_pa)


def locality(model: PCAModel, feature_contigs: np.ndarray, ratio_max: float = 3.0,
             share_min: float = 0.15) -> pd.DataFrame:
    """Flag PCs whose loading mass concentrates on one chromosome (likely biology, not batch)."""
    contigs = pd.unique(feature_contigs)
    expected = np.array([np.mean(feature_contigs == c) for c in contigs])
    rows = []
    for i in range(model.k):
        mass = model.components[i] ** 2
        share = np.array([mass[feature_contigs == c].sum() for c in contigs])
        ratio = share / expected
        j = int(np.argmax(ratio))
        rows.append({"pc": "PC%d" % (i + 1), "top_contig": contigs[j], "share": share[j],
                     "expected_share": expected[j], "ratio": ratio[j],
                     "flagged": bool(ratio[j] > ratio_max and share[j] >= share_min)})
    return pd.DataFrame(rows)


@dataclass
class RobustScore:
    scores: np.ndarray
    weights: np.ndarray
    q: float                    # residual energy over features with weight >= 0.5, rescaled to p
    strong_downweight: float    # fraction of features with Huber weight < 0.5 (|u| > 2c)


def robust_project(z: np.ndarray, model: PCAModel, c: float = 1.345, max_iter: int = 30,
                   tol: float = 1e-6) -> RobustScore:
    """Huber-IRLS projection of one sample onto the PCs.

    A contiguous block of deviating features (a large CNV) gets down-weighted instead of
    dragging the PC scores. Under Gaussian noise ~0.7% of features have weight < 0.5.
    """
    zc = np.asarray(z, float).ravel() - model.mean
    V = model.components
    s = V @ zc
    w = np.ones_like(zc)
    for _ in range(max_iter):
        e = zc - s @ V
        sigma = 1.4826 * np.median(np.abs(e)) + 1e-12
        u = np.abs(e) / sigma
        w = np.where(u <= c, 1.0, c / u)
        Vw = V * w
        s_new = np.linalg.solve(Vw @ V.T, Vw @ zc)
        if np.max(np.abs(s_new - s)) < tol:
            s = s_new
            break
        s = s_new
    e = zc - s @ V
    good = w >= 0.5
    q = float(np.sum(e[good] ** 2) * zc.size / max(good.sum(), 1))
    return RobustScore(s, w, q, float(np.mean(~good)))
