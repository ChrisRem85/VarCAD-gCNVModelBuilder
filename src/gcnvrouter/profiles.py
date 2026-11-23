"""Assay profiles: the only place where WGS / WES / panel behaviour differs."""

from dataclasses import dataclass, replace
from typing import Optional


@dataclass(frozen=True)
class Profile:
    name: str
    # feature eligibility
    min_median_count: float          # cohort median count an interval needs to be a feature
    max_features: Optional[int]      # stratified subsample cap (None = use all eligible)
    # preflight
    zero_depth: float = 0.05         # normalized depth below which a target counts as "dropped out"
    zero_frac_max: float = 0.05      # more dropped-out targets than this -> ASSAY_MISMATCH
    # PCA
    k_pc_max: int = 10
    n_perm: int = 50                 # permutations for parallel analysis
    # clustering
    k_clusters_max: int = 8
    n_min: int = 30                  # minimum cluster size to train a production cohort model
    explore_n_min: int = 10          # minimum cluster size accepted during exploration
    stability_ari_min: float = 0.8
    # gates
    posterior_min: float = 0.90
    margin_min: float = 0.20
    alpha: float = 0.01
    winsor_z: float = 3.0
    huber_c: float = 1.345
    strong_downweight_max: float = 0.15

    def with_overrides(self, **kw) -> "Profile":
        kw = {k: v for k, v in kw.items() if v is not None}
        return replace(self, **kw)


PROFILES = {
    "panel": Profile(name="panel", min_median_count=50, max_features=None,
                     k_pc_max=5, k_clusters_max=6, n_min=30, explore_n_min=3),
    "wes": Profile(name="wes", min_median_count=30, max_features=20000,
                   k_pc_max=10, n_min=50, explore_n_min=10),
    "wgs": Profile(name="wgs", min_median_count=30, max_features=20000,
                   k_pc_max=10, n_min=30, explore_n_min=10),
}


def get_profile(name: str) -> Profile:
    try:
        return PROFILES[name]
    except KeyError:
        raise ValueError("unknown profile '%s' (choose from %s)" % (name, ", ".join(PROFILES)))
