"""Synthetic tests for gcnvrouter phase 1 (exploration)."""

import os

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import adjusted_rand_score

from gcnvrouter.cluster import label_name
from gcnvrouter.explore import explore, write_outputs
from gcnvrouter.features import bed_overlap_mask
from gcnvrouter.io_counts import read_counts
from gcnvrouter.pca import robust_project
from gcnvrouter.profiles import get_profile
from gcnvrouter.simulate import simulate

PANEL = get_profile("panel").with_overrides(explore_n_min=10)


def _labels(exp):
    return np.array([label_name(j) for j in exp.clusters.labels])


def _write_tsv(path, sample, intervals, counts):
    with open(path, "w") as fh:
        fh.write("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n@RG\tID:GATKCopyNumber\tSM:%s\n" % sample)
        fh.write("CONTIG\tSTART\tEND\tCOUNT\n")
        for (c, s, e), n in zip(intervals[["CONTIG", "START", "END"]].itertuples(index=False), counts):
            fh.write("%s\t%d\t%d\t%d\n" % (c, s, e, n))


# ---------------------------------------------------------------- io

def test_read_counts_and_alignment(tmp_path):
    cm, *_ = simulate((3,), n_intervals=60, seed=1)
    paths = []
    for i, s in enumerate(cm.samples):
        p = tmp_path / ("%s.tsv" % s)
        _write_tsv(p, "SM_" + s, cm.intervals, cm.counts[i])
        paths.append(str(p))
    m = read_counts(paths)
    assert m.samples == ["SM_" + s for s in cm.samples]          # name comes from @RG SM
    assert np.array_equal(m.counts, cm.counts)
    bad = cm.intervals.copy()
    bad.loc[5, "START"] += 1
    _write_tsv(tmp_path / "bad.tsv", "bad", bad, cm.counts[0])
    with pytest.raises(ValueError, match="interval mismatch"):
        read_counts(paths + [str(tmp_path / "bad.tsv")])


def test_bed_overlap_mask():
    iv = pd.DataFrame({"CONTIG": ["chr1"] * 4, "START": [100, 200, 300, 400], "END": [150, 250, 350, 450]})
    bed = pd.DataFrame({"contig": ["chr1", "chr1"], "start0": [149, 340], "end": [160, 345]})
    assert bed_overlap_mask(iv, bed).tolist() == [True, False, True, False]


# ---------------------------------------------------------------- clustering

def test_three_batches_recovered():
    cm, batch, *_ = simulate((40, 40, 40), seed=2)
    exp = explore([], PANEL, cm=cm, n_boot=20)
    assert exp.clusters.K == 3
    assert adjusted_rand_score(batch, exp.clusters.labels) >= 0.9


def test_no_batch_gives_one_cluster():
    cm, *_ = simulate((60,), seed=3)
    exp = explore([], PANEL, cm=cm, n_boot=20)
    assert exp.clusters.K == 1


def test_clusters_independent_of_sex():
    cm, batch, sex, _ = simulate((40, 40), seed=4)
    exp = explore([], PANEL, cm=cm, n_boot=20)
    assert adjusted_rand_score(batch, exp.clusters.labels) >= 0.9
    assert abs(adjusted_rand_score(sex, exp.clusters.labels)) < 0.1


def test_sex_inference():
    cm, _, sex, _ = simulate((30,), seed=5)
    exp = explore([], PANEL, cm=cm, n_boot=5)
    assert (exp.pre["inferred_sex"].to_numpy() == sex).mean() == 1.0


# ---------------------------------------------------------------- preflight and gates

def test_assay_mismatch_flagged_and_gated():
    cm, *_ = simulate((40, 40), seed=6)
    drop = np.arange(0, cm.counts.shape[1], 3)                  # a third of the targets missing
    cm.counts[:2, drop] = 0
    exp = explore([], PANEL, cm=cm, n_boot=10)
    assert exp.pre["preflight"].iloc[:2].tolist() == ["ASSAY_MISMATCH"] * 2
    assert not exp.fit_mask[:2].any()
    assert (exp.samples["q_ratio"].iloc[:2] > 1).all()           # Q also catches them independently


def test_q_detects_novel_batch():
    cm, *_ = simulate((40, 40), seed=7)
    exp = explore([], PANEL, cm=cm, n_boot=10)
    novel, *_ = simulate((20,), seed=99, batch_effect_sd=0.3, intervals=cm.intervals)
    results = [exp.score(row) for row in novel.counts]
    q_flag = np.mean([r["q"] > r["q_threshold"] for r in results])
    assert q_flag >= 0.9


def test_in_distribution_samples_mostly_pass():
    cm, batch, *_ = simulate((50, 50), seed=8)
    held = np.zeros(len(batch), bool)
    held[40:50] = held[90:100] = True                            # 10 per batch never seen by the fit
    exp = explore([], PANEL, cm=cm.subset(~held), n_boot=10)
    results = [exp.score(row) for row in cm.counts[held]]
    assert np.mean([r["status"] == "OOD" for r in results]) <= 0.1
    assigned = np.array([r["cluster"] for r in results])
    assert adjusted_rand_score(batch[held], assigned) == 1.0


def test_large_cnv_keeps_cluster():
    cm, batch, *_ = simulate((40, 40), seed=9)
    exp = explore([], PANEL, cm=cm, n_boot=10)
    i = 3
    base = exp.score(cm.counts[i])
    row = cm.counts[i].copy()
    auto = np.flatnonzero(cm.intervals["CONTIG"].str.match(r"chr\d+$").to_numpy())
    block = auto[200:200 + int(0.05 * auto.size)]               # contiguous 5% of features, het deletion
    row[block] *= 0.5
    cnv = exp.score(row)
    assert cnv["cluster"] == base["cluster"]
    assert cnv["strong_downweight"] > base["strong_downweight"]
    assert cnv["status"] != "OOD"


def test_robust_projection_matches_plain_without_outliers():
    cm, *_ = simulate((40, 40), seed=10)
    exp = explore([], PANEL, cm=cm, n_boot=5)
    rs = robust_project(exp.Z[0], exp.pca)
    plain = exp.pca.scores(exp.Z[0:1])[0]
    assert np.allclose(rs.scores, plain, atol=0.15 * np.abs(plain).max())


def test_locality_flags_chromosome_wide_event():
    cm, *_ = simulate((80,), seed=11, dup_chromosome="chr5", dup_fraction=0.5)
    exp = explore([], PANEL, cm=cm, n_boot=5)
    flagged = exp.loc[exp.loc["flagged"]]
    assert "chr5" in flagged["top_contig"].tolist()
    assert exp.clusters.K == 1                                   # biology must not become a "batch"


# ---------------------------------------------------------------- outputs

def test_outputs_and_determinism(tmp_path):
    cm, *_ = simulate((30, 30), seed=12)
    a = write_outputs(explore([], PANEL, cm=cm, n_boot=10), str(tmp_path / "a"), html=True)
    b = write_outputs(explore([], PANEL, cm=cm, n_boot=10), str(tmp_path / "b"), html=False)
    for f in ("samples.tsv", "sweep.tsv", "pca.tsv", "features.tsv", "locality.tsv"):
        assert (tmp_path / "a" / f).read_bytes() == (tmp_path / "b" / f).read_bytes(), f
    assert a["K"] == b["K"]
    html = (tmp_path / "a" / "report.html").read_text(encoding="utf-8")
    assert "Batch exploration" in html and "plotly" in html.lower()
    assert os.path.getsize(tmp_path / "a" / "report.html") > 100_000   # plotly.js is inlined
