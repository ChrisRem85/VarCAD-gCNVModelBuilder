"""Phase-1 exploration: preflight -> features -> PCA -> clustering -> gates -> outputs.

`explore()` returns an `Exploration` that holds every intermediate result, so the report,
tests and (later) the router serialization all read the same objects.
"""

import json
import os
import platform
from dataclasses import dataclass, field
from typing import List, Optional, Sequence

import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hier

from . import __version__
from .cluster import ClusterResult, choose_clusters, label_name
from .features import (FeatureModel, gc_fingerprint, log_normalize, preflight, read_annotated_intervals,
                       read_bed, select_features)
from .gates import ClusterStats, cluster_stats, status, t2
from .io_counts import CountMatrix, read_counts
from .pca import PCAModel, fit_pca, locality, robust_project
from .profiles import Profile


@dataclass
class Exploration:
    profile: Profile
    cm: CountMatrix
    pre: pd.DataFrame
    fit_mask: np.ndarray
    feature_idx: np.ndarray
    feature_model: FeatureModel
    X: np.ndarray                      # log-normalized features, fit samples
    Z: np.ndarray                      # scaled/winsorized features, fit samples
    pca: PCAModel
    loc: pd.DataFrame
    pcs_used: List[int]                # PCs used for clustering (non-flagged)
    clusters: ClusterResult
    cstats: List[ClusterStats]
    cv_q: np.ndarray
    cv_strong: np.ndarray
    q_threshold: float
    samples: pd.DataFrame
    seed: int
    inputs: List[str]
    metadata_columns: List[str] = field(default_factory=list)
    gc_centres: Optional[np.ndarray] = None
    gc_profile: Optional[np.ndarray] = None
    warnings: List[str] = field(default_factory=list)

    @property
    def fit_samples(self) -> List[str]:
        return [s for s, m in zip(self.cm.samples, self.fit_mask) if m]

    def heatmap_order(self) -> np.ndarray:
        """Fit-sample order: by cluster, then hierarchical (correlation) within cluster."""
        if self.Z.shape[0] > 2:
            leaves = hier.leaves_list(hier.linkage(self.Z, "average", metric="correlation"))
        else:
            leaves = np.arange(self.Z.shape[0])
        rank = np.empty_like(leaves)
        rank[leaves] = np.arange(leaves.size)
        return np.lexsort((rank, self.clusters.labels))

    def score(self, counts_row: np.ndarray) -> dict:
        """Gate one sample (counts over the full interval set) against this exploration's model."""
        x = log_normalize(counts_row, self.feature_idx)
        z = self.feature_model.transform(x)[0]
        rs = robust_project(z, self.pca, self.profile.huber_c)
        p_use = rs.scores[self.pcs_used][None, :]
        post = self.clusters.model.predict_proba(p_use)[0][self.clusters.order]
        j = int(np.argmax(post))
        srt = np.sort(post)[::-1]
        margin = float(srt[0] - srt[1]) if post.size > 1 else 1.0
        cs = self.cstats[j]
        t2v = float(t2(p_use, cs)[0])
        return {"scores": rs.scores, "cluster": j, "posterior": float(post[j]), "margin": margin,
                "t2": t2v, "t2_threshold": cs.t2_threshold, "q": rs.q, "q_threshold": self.q_threshold,
                "strong_downweight": rs.strong_downweight,
                "status": status(post[j], margin, t2v / cs.t2_threshold, rs.q / self.q_threshold,
                                 rs.strong_downweight, self.profile)}


def cluster_spread(P: np.ndarray, labels: np.ndarray, K: int) -> np.ndarray:
    """Median distance of members to their cluster median (PC units)."""
    out = np.full(K, np.nan)
    for j in range(K):
        X = P[labels == j]
        if len(X) >= 2:
            out[j] = float(np.median(np.linalg.norm(X - np.median(X, axis=0), axis=1)))
    return out


def cross_fitted_q(X: np.ndarray, feature_idx: np.ndarray, profile: Profile, k: int, seed: int):
    """Q (robust SPE) of each fit sample computed from a model that did not see it."""
    n = X.shape[0]
    folds = min(10, n)
    rng = np.random.default_rng(seed)
    assign = rng.permutation(np.arange(n) % folds)
    q = np.empty(n)
    strong = np.empty(n)
    for f in range(folds):
        test = assign == f
        fm = FeatureModel.fit(X[~test], feature_idx, profile.winsor_z)
        pm = fit_pca(fm.transform(X[~test]), profile.k_pc_max, n_perm=0, k_fixed=min(k, int((~test).sum()) - 1))
        for i in np.flatnonzero(test):
            rs = robust_project(fm.transform(X[i])[0], pm, profile.huber_c)
            q[i], strong[i] = rs.q, rs.strong_downweight
    return q, strong


def _read_metadata(path: str) -> pd.DataFrame:
    md = pd.read_csv(path, sep="\t", dtype=str)
    key = "sample_id" if "sample_id" in md.columns else ("sample" if "sample" in md.columns else None)
    if key is None:
        raise ValueError("metadata %s needs a 'sample_id' or 'sample' column" % path)
    md = md.rename(columns={key: "sample"})
    for c in md.columns:
        if c != "sample":
            num = pd.to_numeric(md[c], errors="coerce")
            if num.notna().all():
                md[c] = num
    return md


def explore(paths: Sequence[str], profile: Profile, metadata_path: Optional[str] = None,
            annotated_path: Optional[str] = None, mask_bed_paths: Sequence[str] = (),
            exclude: Sequence[str] = (), seed: int = 1729, n_boot: int = 100,
            cm: Optional[CountMatrix] = None) -> Exploration:
    cm = cm if cm is not None else read_counts(paths)
    warnings = []
    pre = preflight(cm, profile)
    if exclude:
        sel = pre["sample"].isin(list(exclude)) & (pre["preflight"] == "OK")
        pre.loc[sel, "preflight"] = "EXCLUDED"
    fit_mask = (pre["preflight"] == "OK").to_numpy()
    n_fit = int(fit_mask.sum())
    if n_fit < 4:
        raise ValueError("only %d samples passed preflight; need at least 4" % n_fit)
    for st, n in pre["preflight"].value_counts().items():
        if st != "OK":
            warnings.append("%d sample(s) excluded from the fit by preflight: %s" % (n, st))

    annotated = read_annotated_intervals(annotated_path, cm.intervals) if annotated_path else None
    beds = [read_bed(p) for p in mask_bed_paths]
    fidx = select_features(cm, fit_mask, profile, annotated, beds, seed)
    X = log_normalize(cm.counts[fit_mask], fidx)
    fm = FeatureModel.fit(X, fidx, profile.winsor_z)
    Z = fm.transform(X)

    pca = fit_pca(Z, profile.k_pc_max, profile.n_perm, seed)
    if pca.k_pa == 0:
        warnings.append("no PC exceeds the parallel-analysis threshold: no detectable batch structure")
    contigs = cm.intervals["CONTIG"].to_numpy()[fidx]
    loc = locality(pca, contigs)
    pcs_used = [i for i in range(pca.k) if not loc["flagged"].iloc[i]]
    for _, r in loc[loc["flagged"]].iterrows():
        warnings.append("%s concentrates %.0f%% of its loading mass on %s (expected %.0f%%): likely biology, "
                        "excluded from clustering" % (r["pc"], 100 * r["share"], r["top_contig"],
                                                      100 * r["expected_share"]))
    if not pcs_used:
        pcs_used = list(range(pca.k))
        warnings.append("all PCs flagged by the locality check; using them anyway")
    scores = pca.scores(Z)
    P = scores[:, pcs_used]

    cl = choose_clusters(P, profile.k_clusters_max, profile.explore_n_min, profile.stability_ari_min,
                         n_boot, seed)
    warnings.extend(cl.notes)
    cs = cluster_stats(P, cl.labels, cl.K, profile.alpha, seed)
    spread = cluster_spread(P, cl.labels, cl.K)
    if cl.K > 1 and np.nanmax(spread) > 3 * np.nanmin(spread):
        j = int(np.nanargmax(spread))
        warnings.append("%s is %.1fx more spread out than the tightest cluster: probably several small batches "
                        "or outliers merged into one catch-all cluster" % (label_name(j), spread[j] / np.nanmin(spread)))
    small =[label_name(j) for j in range(cl.K) if cs[j].n < profile.n_min]
    if small:
        warnings.append("cluster(s) %s are below the production minimum of %d samples: not yet trainable as "
                        "separate cohort models" % (", ".join(small), profile.n_min))

    cvq, cvs = cross_fitted_q(X, fidx, profile, pca.k, seed)
    q_thr = float(np.quantile(cvq, 1 - profile.alpha, method="higher"))

    # ---------------- per-sample table ----------------
    rows = {s: {} for s in cm.samples}
    fit_names = [s for s, m in zip(cm.samples, fit_mask) if m]
    for i, s in enumerate(fit_names):
        j = int(cl.labels[i])
        post = cl.posteriors[i]
        srt = np.sort(post)[::-1]
        margin = float(srt[0] - srt[1]) if post.size > 1 else 1.0
        t2v = float(t2(P[i:i + 1], cs[j])[0])
        r = rows[s]
        r.update({"used_in_fit": True, "cluster": label_name(j), "posterior": float(post[j]), "margin": margin,
                  "t2": t2v, "t2_threshold": cs[j].t2_threshold, "q": cvq[i], "q_threshold": q_thr,
                  "strong_downweight": cvs[i]})
        for k in range(pca.k):
            r["PC%d" % (k + 1)] = scores[i, k]
        r["status"] = status(post[j], margin, t2v / cs[j].t2_threshold, cvq[i] / q_thr, cvs[i], profile)
    exp_tmp = Exploration(profile, cm, pre, fit_mask, fidx, fm, X, Z, pca, loc, pcs_used, cl, cs, cvq, cvs,
                          q_thr, pd.DataFrame(), seed, list(paths))
    for idx in np.flatnonzero(~fit_mask):
        s = cm.samples[idx]
        r = rows[s]
        r["used_in_fit"] = False
        try:
            sc = exp_tmp.score(cm.counts[idx])
        except Exception as exc:      # e.g. zero coverage
            r["status"] = pre.loc[idx, "preflight"]
            warnings.append("could not score %s: %s" % (s, exc))
            continue
        r.update({"cluster": "(" + label_name(sc["cluster"]) + ")", "posterior": sc["posterior"],
                  "margin": sc["margin"], "t2": sc["t2"], "t2_threshold": sc["t2_threshold"], "q": sc["q"],
                  "q_threshold": sc["q_threshold"], "strong_downweight": sc["strong_downweight"]})
        for k in range(pca.k):
            r["PC%d" % (k + 1)] = sc["scores"][k]
        # preflight verdict wins; the gate result is kept for the plots
        r["gate_status"] = sc["status"]
        r["status"] = pre.loc[idx, "preflight"]
    table = pre.merge(pd.DataFrame.from_dict(rows, orient="index").rename_axis("sample").reset_index(),
                      on="sample", how="left")
    table["t2_ratio"] = table["t2"] / table["t2_threshold"]      # 0 when the limit is inf (cluster too small)
    table["q_ratio"] = table["q"] / table["q_threshold"]
    md_cols = []
    if metadata_path:
        md = _read_metadata(metadata_path)
        md_cols = [c for c in md.columns if c != "sample"]
        table = table.merge(md, on="sample", how="left")
        missing = sorted(set(cm.samples) - set(md["sample"]))
        if missing:
            warnings.append("%d sample(s) have no metadata row" % len(missing))

    gc_c = gc_p = None
    if annotated is not None and "GC_CONTENT" in annotated:
        gc_c, gc_p = gc_fingerprint(fm.residual(X), annotated["GC_CONTENT"].to_numpy()[fidx])

    exp_tmp.samples = table
    exp_tmp.metadata_columns = md_cols
    exp_tmp.gc_centres, exp_tmp.gc_profile = gc_c, gc_p
    exp_tmp.warnings = warnings
    return exp_tmp


# --------------------------------------------------------------------------------------------
# outputs
# --------------------------------------------------------------------------------------------

def summary(exp: Exploration) -> dict:
    cl = exp.clusters
    sizes = [int(c.n) for c in exp.cstats]
    return {
        "gcnvrouter_version": __version__,
        "python": platform.python_version(),
        "profile": exp.profile.name,
        "seed": exp.seed,
        "n_samples": len(exp.cm.samples),
        "n_fit": int(exp.fit_mask.sum()),
        "preflight": exp.pre["preflight"].value_counts().to_dict(),
        "n_intervals": int(len(exp.cm.intervals)),
        "n_features": int(exp.feature_idx.size),
        "k_parallel_analysis": int(exp.pca.k_pa),
        "k_pcs": int(exp.pca.k),
        "pcs_used_for_clustering": ["PC%d" % (i + 1) for i in exp.pcs_used],
        "K": int(cl.K),
        "covariance": cl.cov,
        "cluster_sizes": dict(zip([label_name(j) for j in range(cl.K)], sizes)),
        "bootstrap_median_ari": float(np.median(cl.aris)) if cl.aris.size else None,
        "explore_n_min": exp.profile.explore_n_min,
        "production_n_min": exp.profile.n_min,
        "all_clusters_trainable": bool(min(sizes) >= exp.profile.n_min),
        "q_threshold": exp.q_threshold,
        "t2_thresholds": {label_name(j): c.t2_threshold for j, c in enumerate(exp.cstats)},
        "status_counts": exp.samples["status"].value_counts().to_dict(),
        "warnings": exp.warnings,
        "inputs": exp.inputs,
    }


def write_outputs(exp: Exploration, out_dir: str, html: bool = True) -> dict:
    os.makedirs(out_dir, exist_ok=True)
    exp.samples.to_csv(os.path.join(out_dir, "samples.tsv"), sep="\t", index=False, float_format="%.6g")
    exp.clusters.sweep.to_csv(os.path.join(out_dir, "sweep.tsv"), sep="\t", index=False, float_format="%.6g")
    n_pc = len(exp.pca.eigvals)
    pd.DataFrame({"pc": ["PC%d" % (i + 1) for i in range(n_pc)], "eigenvalue": exp.pca.eigvals,
                  "explained": exp.pca.explained, "pa_threshold_95": exp.pca.pa_threshold}).to_csv(
        os.path.join(out_dir, "pca.tsv"), sep="\t", index=False, float_format="%.6g")
    exp.loc.to_csv(os.path.join(out_dir, "locality.tsv"), sep="\t", index=False, float_format="%.6g")
    feats = exp.cm.intervals.iloc[exp.feature_idx].reset_index(drop=True)
    for i in range(exp.pca.k):
        feats["loading_PC%d" % (i + 1)] = exp.pca.components[i]
    feats.to_csv(os.path.join(out_dir, "features.tsv"), sep="\t", index=False, float_format="%.6g")
    summ = summary(exp)
    with open(os.path.join(out_dir, "summary.json"), "w") as fh:
        json.dump(summ, fh, indent=2, default=str)
    if html:
        from .report import build_report
        build_report(exp, os.path.join(out_dir, "report.html"))
    return summ
