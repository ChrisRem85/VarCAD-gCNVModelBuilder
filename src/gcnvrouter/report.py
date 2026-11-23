"""Self-contained HTML report (plotly.js inlined, works offline).

Colours follow the dataviz reference palette: categorical slots in fixed order with marker
symbols as secondary encoding, a single-hue blue ramp for magnitude, blue<->red with a grey
midpoint for correlation, and the reserved status palette (always with a text label).
"""

import datetime
import html
from typing import List

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import get_plotlyjs
from plotly.subplots import make_subplots
from scipy import stats

from .cluster import label_name

SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT2 = "#52514e"
MUTED = "#9e9d97"
GRID = "#ebeae6"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SYMBOLS = ["circle", "square", "diamond", "triangle-up", "cross", "x", "star", "hexagon"]
STATUS = {"ROUTED": "#0ca30c", "AMBIGUOUS": "#fab219", "OOD": "#ec835a", "QC_FAIL": "#d03b3b",
          "ASSAY_MISMATCH": "#d03b3b", "NO_COVERAGE": "#d03b3b", "EXCLUDED": MUTED}
STATUS_SYMBOL = {"ROUTED": "circle", "AMBIGUOUS": "diamond", "OOD": "square", "QC_FAIL": "x",
                 "ASSAY_MISMATCH": "x", "NO_COVERAGE": "x", "EXCLUDED": "circle-open"}
SEQ = [[0.0, "#cde2fb"], [0.25, "#86b6ef"], [0.5, "#3987e5"], [0.75, "#1c5cab"], [1.0, "#0d366b"]]
DIV = [[0.0, "#1c5cab"], [0.25, "#86b6ef"], [0.5, "#f0efec"], [0.75, "#f08a89"], [1.0, "#b8302f"]]
CONFIG = {"displaylogo": False, "responsive": True}


def _style(fig: go.Figure, height: int = 420, **kw) -> go.Figure:
    fig.update_layout(paper_bgcolor=SURFACE, plot_bgcolor=SURFACE, height=height,
                      font=dict(family="system-ui, -apple-system, 'Segoe UI', sans-serif", size=12, color=TEXT),
                      margin=dict(l=64, r=24, t=48, b=56),
                      legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0, font=dict(color=TEXT2)),
                      hoverlabel=dict(bgcolor="white", font=dict(color=TEXT)), **kw)
    fig.update_xaxes(gridcolor=GRID, linecolor=GRID, zeroline=False, showline=True, tickfont=dict(color=TEXT2))
    fig.update_yaxes(gridcolor=GRID, linecolor=GRID, zeroline=False, showline=True, tickfont=dict(color=TEXT2))
    return fig


def _div(fig: go.Figure) -> str:
    return fig.to_html(full_html=False, include_plotlyjs=False, config=CONFIG)


def _marker(i: int, size: int = 10) -> dict:
    return dict(color=SERIES[i % len(SERIES)], symbol=SYMBOLS[i % len(SYMBOLS)], size=size,
                line=dict(color=SURFACE, width=1.5))


def _fmt(v, digits=3):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return ""
    if isinstance(v, (float, np.floating)):
        if np.isinf(v):
            return "∞"
        return ("%." + str(digits) + "g") % v
    return str(v)


def _table(df: pd.DataFrame, cls: str = "") -> str:
    head = "".join("<th>%s</th>" % html.escape(str(c)) for c in df.columns)
    body = "".join("<tr>%s</tr>" % "".join("<td>%s</td>" % html.escape(_fmt(v)) for v in row)
                   for row in df.itertuples(index=False))
    return '<table class="%s"><thead><tr>%s</tr></thead><tbody>%s</tbody></table>' % (cls, head, body)


# --------------------------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------------------------

def fig_dropout(exp) -> go.Figure:
    pre = exp.pre.sort_values("dropout_frac", ascending=False)
    fig = go.Figure()
    for ok, name, color in [(True, "passed preflight", SERIES[0]), (False, "failed preflight", STATUS["QC_FAIL"])]:
        sel = (pre["preflight"].isin(["OK", "EXCLUDED"])) == ok
        if sel.any():
            fig.add_bar(x=pre.loc[sel, "sample"], y=100 * pre.loc[sel, "dropout_frac"], name=name,
                        marker=dict(color=color), customdata=pre.loc[sel, ["preflight"]],
                        hovertemplate="%{x}<br>%{y:.1f}% targets dropped out<br>%{customdata[0]}<extra></extra>")
    fig.add_hline(y=100 * exp.profile.zero_frac_max, line=dict(color=TEXT2, width=1),
                  annotation_text="assay-mismatch limit %.0f%%" % (100 * exp.profile.zero_frac_max),
                  annotation_font_color=TEXT2)
    fig.update_layout(barcornerradius=3, xaxis=dict(categoryorder="array", categoryarray=list(pre["sample"])))
    fig.update_yaxes(title="targets with depth < %.0f%% of median (%%)" % (100 * exp.profile.zero_depth))
    return _style(fig, 380)


def fig_sex(exp) -> go.Figure:
    pre = exp.pre
    fig = go.Figure()
    for i, sx in enumerate(["female", "male", "NA"]):
        sel = pre["inferred_sex"] == sx
        if not sel.any():
            continue
        m = _marker(i) if sx != "NA" else dict(color=MUTED, symbol="diamond", size=10)
        fig.add_scatter(x=pre.loc[sel, "chrX_ratio"], y=pre.loc[sel, "chrY_ratio"], mode="markers", name=sx,
                        marker=m, text=pre.loc[sel, "sample"],
                        hovertemplate="%{text}<br>chrX %{x:.2f} · chrY %{y:.3f}<extra>" + sx + "</extra>")
    fig.update_xaxes(title="chrX depth / autosomal median")
    fig.update_yaxes(title="chrY depth / autosomal median")
    return _style(fig, 380)


def fig_scree(exp) -> go.Figure:
    n = min(len(exp.pca.eigvals), 15)
    x = ["PC%d" % (i + 1) for i in range(n)]
    fig = go.Figure()
    fig.add_bar(x=x, y=exp.pca.eigvals[:n], name="observed variance", marker=dict(color=SERIES[0]),
                customdata=100 * exp.pca.explained[:n],
                hovertemplate="%{x}: %{y:.1f} (%{customdata:.1f}% of total)<extra></extra>")
    if not np.all(np.isnan(exp.pca.pa_threshold)):
        fig.add_scatter(x=x, y=exp.pca.pa_threshold[:n], mode="lines+markers", name="noise level (parallel analysis, 95th pct)",
                        line=dict(color=TEXT2, width=2), marker=dict(size=8, color=TEXT2),
                        hovertemplate="%{x}: noise %{y:.1f}<extra></extra>")
    fig.update_layout(barcornerradius=3)
    fig.update_yaxes(title="variance (eigenvalue)")
    return _style(fig, 380)


def _color_variables(exp) -> List[str]:
    vars_ = ["cluster", "status", "inferred_sex", "total_counts", "median_count", "dropout_frac"]
    return vars_ + [c for c in exp.metadata_columns if c not in vars_]


def fig_splom(exp, var: str) -> go.Figure:
    df = exp.samples[exp.samples["used_in_fit"] == True].copy()  # noqa: E712
    pcs = ["PC%d" % (i + 1) for i in range(max(2, min(exp.pca.k, 4)))]
    dims = [dict(label="%s (%.0f%%)" % (p, 100 * exp.pca.explained[int(p[2:]) - 1]), values=df[p]) for p in pcs]
    hover = ("%{text}<br>cluster %{customdata[0]} · %{customdata[1]}<br>" + var + ": %{customdata[2]}<extra></extra>")
    custom = np.column_stack([df["cluster"], df["status"], df[var].map(_fmt)])
    fig = go.Figure()
    values = df[var]
    numeric = pd.api.types.is_numeric_dtype(values) and values.nunique() > 8
    common = dict(dimensions=dims, diagonal_visible=False, showupperhalf=False, hovertemplate=hover)
    if numeric:
        fig.add_trace(go.Splom(**common, text=df["sample"], customdata=custom, name=var,
                               marker=dict(color=values, colorscale=SEQ, size=10, showscale=True,
                                           colorbar=dict(title=var, thickness=12),
                                           line=dict(color=SURFACE, width=1.5))))
    else:
        cats = list(pd.unique(values.astype(str)))
        if var == "cluster":
            cats = sorted(cats)
        if len(cats) > 8:
            top = values.astype(str).value_counts().index[:7]
            values = values.astype(str).where(values.astype(str).isin(top), "Other")
            cats = list(top) + ["Other"]
        for i, c in enumerate(cats):
            sel = (values.astype(str) == c).to_numpy()
            if var == "status":
                m = dict(color=STATUS.get(c, MUTED), symbol=STATUS_SYMBOL.get(c, "circle"), size=10,
                         line=dict(color=SURFACE, width=1.5))
            else:
                m = _marker(i)
            d = [dict(label=x["label"], values=np.asarray(x["values"])[sel]) for x in dims]
            fig.add_trace(go.Splom(dimensions=d, diagonal_visible=False, showupperhalf=False, name=c,
                                   text=df["sample"][sel], customdata=custom[sel], marker=m,
                                   hovertemplate=hover))
    side = 220 * (len(pcs) - 1) + 80
    return _style(fig, max(side, 420))


def fig_bic(exp) -> go.Figure:
    sw = exp.clusters.sweep
    fig = go.Figure()
    for i, cov in enumerate(sorted(sw["cov"].unique(), key=lambda c: ["tied", "diag"].index(c) if c in ("tied", "diag") else 9)):
        s = sw[sw["cov"] == cov]
        symbols = [SYMBOLS[i] if v else SYMBOLS[i] + "-open" for v in s["valid"]]
        fig.add_scatter(x=s["K"], y=s["bic"], mode="lines+markers", name="%s covariance" % cov,
                        line=dict(color=SERIES[i], width=2),
                        marker=dict(color=SERIES[i], symbol=symbols, size=11, line=dict(color=SERIES[i], width=1.5)),
                        customdata=np.column_stack([s["sizes"], s["valid"], s["median_ari"].map(_fmt)]),
                        hovertemplate="K=%{x} · BIC %{y:.1f}<br>sizes %{customdata[0]}<br>valid: %{customdata[1]}"
                                      "<br>stability ARI %{customdata[2]}<extra></extra>")
    ch = sw[sw["chosen"]]
    if len(ch):
        r = ch.iloc[0]
        fig.add_annotation(x=r["K"], y=r["bic"], text="chosen: K=%d (%s)" % (r["K"], r["cov"]), showarrow=True,
                           arrowcolor=TEXT2, font=dict(color=TEXT), ay=-40)
    fig.update_xaxes(title="number of clusters K", dtick=1)
    fig.update_yaxes(title="BIC (lower is better)")
    return _style(fig, 380)


def fig_ari(exp) -> go.Figure:
    fig = go.Figure()
    fig.add_histogram(x=exp.clusters.aris, xbins=dict(start=-0.05, end=1.05, size=0.05), name="subsample refits",
                      marker=dict(color=SERIES[0], line=dict(color=SURFACE, width=2)),
                      hovertemplate="ARI %{x}<br>%{y} refits<extra></extra>")
    fig.add_vline(x=exp.profile.stability_ari_min, line=dict(color=TEXT2, width=1),
                  annotation_text="stability limit", annotation_font_color=TEXT2)
    fig.update_xaxes(title="adjusted Rand index vs. full fit", range=[-0.05, 1.05])
    fig.update_yaxes(title="subsample refits")
    return _style(fig, 340, showlegend=False)


def _ordered_labels(exp):
    order = exp.heatmap_order()
    names = np.array(exp.fit_samples)[order]
    labs = np.array([label_name(j) for j in exp.clusters.labels])[order]
    return order, ["%s · %s" % (l, n) for l, n in zip(labs, names)]


def fig_coassign(exp) -> go.Figure:
    order, ticks = _ordered_labels(exp)
    m = exp.clusters.coassign[np.ix_(order, order)]
    fig = go.Figure(go.Heatmap(z=m, x=ticks, y=ticks, colorscale=SEQ, zmin=0, zmax=1, xgap=1, ygap=1,
                               colorbar=dict(title="same cluster", thickness=12),
                               hovertemplate="%{y}<br>%{x}<br>together in %{z:.0%} of refits<extra></extra>"))
    fig.update_yaxes(autorange="reversed", showgrid=False)
    fig.update_xaxes(showgrid=False, tickangle=-60)
    return _style(fig, 300 + 16 * len(ticks))


def fig_corr(exp) -> go.Figure:
    order, ticks = _ordered_labels(exp)
    c = np.corrcoef(exp.Z)[np.ix_(order, order)]
    fig = go.Figure(go.Heatmap(z=c, x=ticks, y=ticks, colorscale=DIV, zmin=-1, zmax=1, xgap=1, ygap=1,
                               colorbar=dict(title="Pearson r", thickness=12),
                               hovertemplate="%{y}<br>%{x}<br>r = %{z:.2f}<extra></extra>"))
    fig.update_yaxes(autorange="reversed", showgrid=False)
    fig.update_xaxes(showgrid=False, tickangle=-60)
    return _style(fig, 300 + 16 * len(ticks))


def fig_loadings(exp) -> go.Figure:
    k = exp.pca.k
    contigs = exp.cm.intervals["CONTIG"].to_numpy()[exp.feature_idx]
    uniq = list(pd.unique(contigs))
    parity = np.array([uniq.index(c) % 2 for c in contigs])
    x = np.arange(contigs.size)
    titles = []
    for i in range(k):
        r = exp.loc.iloc[i]
        t = "PC%d · largest share: %s holds %.0f%% of loading mass (expected %.0f%%)" % (
            i + 1, r["top_contig"], 100 * r["share"], 100 * r["expected_share"])
        if r["flagged"]:
            t += " · flagged as likely biology"
        titles.append(t)
    fig = make_subplots(rows=k, cols=1, shared_xaxes=True, vertical_spacing=0.06, subplot_titles=titles)
    for i in range(k):
        for par, color in [(0, SERIES[0]), (1, MUTED)]:
            sel = parity == par
            fig.add_scatter(x=x[sel], y=exp.pca.components[i][sel], mode="markers", row=i + 1, col=1,
                            marker=dict(color=color, size=4), showlegend=False,
                            text=contigs[sel], hovertemplate="%{text}<br>loading %{y:.3f}<extra></extra>")
    # label only chromosomes wide enough to hold a label (tiny ones would collide)
    shown = [c for c in uniq if np.mean(contigs == c) >= 0.012]
    mids = [float(np.mean(x[contigs == c])) for c in shown]
    fig.update_xaxes(tickvals=mids, ticktext=[c.replace("chr", "") for c in shown], row=k, col=1,
                     title="feature (genomic order, chromosome)")
    fig.update_annotations(font=dict(size=12, color=TEXT2), x=0, xanchor="left")
    return _style(fig, 170 * k + 90)


def fig_gates(exp) -> go.Figure:
    df = exp.samples[exp.samples["q_ratio"].notna()].copy()
    # T² is not testable when the cluster has <= p members (limit = inf): pin those to the left edge
    finite = df["t2_ratio"].replace(0, np.nan).dropna()
    floor = float(finite.min() / 3) if len(finite) else 1e-3
    df["t2_plot"] = df["t2_ratio"].where(df["t2_ratio"] > 0, floor)
    df["t2_text"] = [("%.3g" % v) if v > 0 else "not testable (cluster too small)" for v in df["t2_ratio"].fillna(0)]
    fig = go.Figure()
    fit = df[df["used_in_fit"] == True]  # noqa: E712
    for j in range(exp.clusters.K):
        s = fit[fit["cluster"] == label_name(j)]
        fig.add_scatter(x=s["t2_plot"], y=s["q_ratio"], mode="markers", name="%s (in fit)" % label_name(j),
                        marker=_marker(j, 11), text=s["sample"], customdata=s[["status", "t2_text"]],
                        hovertemplate="%{text}<br>T² / limit %{customdata[1]} · Q / limit %{y:.2f}"
                                      "<br>%{customdata[0]}<extra></extra>")
    other = df[df["used_in_fit"] != True]  # noqa: E712
    for st in pd.unique(other["status"]):
        s = other[other["status"] == st]
        fig.add_scatter(x=s["t2_plot"], y=s["q_ratio"], mode="markers", name="%s (not in fit)" % st,
                        marker=dict(color=STATUS.get(st, MUTED), symbol=STATUS_SYMBOL.get(st, "x"), size=12,
                                    line=dict(color=STATUS.get(st, MUTED), width=2)),
                        text=s["sample"], customdata=s[["cluster", "t2_text"]],
                        hovertemplate="%{text}<br>T² / limit %{customdata[1]} · Q / limit %{y:.2f}"
                                      "<br>nearest %{customdata[0]}<extra>" + st + "</extra>")
    if (df["t2_ratio"].fillna(0) <= 0).any():
        fig.add_vline(x=floor, line=dict(color=GRID, width=1))
        fig.add_annotation(x=np.log10(floor), y=1, yref="paper", xanchor="left", yanchor="top", showarrow=False,
                           text="T² not testable", font=dict(color=TEXT2, size=11))
    fig.add_vline(x=1, line=dict(color=TEXT2, width=1))
    fig.add_hline(y=1, line=dict(color=TEXT2, width=1))
    fig.add_annotation(x=0, y=0, xref="paper", yref="paper", xanchor="left", yanchor="bottom", showarrow=False,
                       text="inside both limits = in distribution", font=dict(color=TEXT2))
    fig.update_xaxes(type="log", title="Hotelling T² / limit (inside the PC subspace)")
    fig.update_yaxes(type="log", title="Q / limit (outside the PC subspace)")
    return _style(fig, 460)


def fig_posterior(exp) -> go.Figure:
    df = exp.samples[exp.samples["used_in_fit"] == True].sort_values(["cluster", "posterior"])  # noqa: E712
    fig = go.Figure()
    for j in range(exp.clusters.K):
        s = df[df["cluster"] == label_name(j)]
        fig.add_bar(x=s["sample"], y=s["posterior"], name=label_name(j), marker=dict(color=SERIES[j % 8]),
                    customdata=s[["margin"]],
                    hovertemplate="%{x}<br>posterior %{y:.3f} · margin %{customdata[0]:.3f}<extra></extra>")
    fig.add_hline(y=exp.profile.posterior_min, line=dict(color=TEXT2, width=1),
                  annotation_text="ambiguity limit", annotation_font_color=TEXT2)
    fig.update_layout(barcornerradius=3, xaxis=dict(categoryorder="array", categoryarray=list(df["sample"])))
    fig.update_yaxes(title="posterior of assigned cluster", range=[0, 1.02])
    return _style(fig, 360)


def fig_gc(exp) -> go.Figure:
    fig = go.Figure()
    for j in range(exp.clusters.K):
        rows = np.flatnonzero(exp.clusters.labels == j)
        for n, i in enumerate(rows):
            fig.add_scatter(x=exp.gc_centres, y=exp.gc_profile[i], mode="lines", legendgroup=label_name(j),
                            showlegend=False, line=dict(color=SERIES[j % 8], width=1), opacity=0.45,
                            name=exp.fit_samples[i], hovertemplate="%{fullData.name}<br>GC %{x:.2f}: %{y:.3f}<extra></extra>")
        fig.add_scatter(x=exp.gc_centres, y=np.nanmedian(exp.gc_profile[rows], axis=0), mode="lines",
                        name="%s median" % label_name(j), legendgroup=label_name(j),
                        line=dict(color=SERIES[j % 8], width=3))
    fig.update_xaxes(title="GC content of interval")
    fig.update_yaxes(title="median residual log2 depth")
    return _style(fig, 380)


# --------------------------------------------------------------------------------------------
# tables
# --------------------------------------------------------------------------------------------

def _contingency(exp, col: str) -> str:
    df = exp.samples[exp.samples["used_in_fit"] == True]  # noqa: E712
    if col not in df or df[col].isna().all():
        return ""
    tab = pd.crosstab(df["cluster"], df[col].astype(str))
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        return "<h4>cluster × %s</h4>%s<p class='caption'>Only one level: no test.</p>" % (
            html.escape(col), _table(tab.reset_index()))
    chi2, p, _, _ = stats.chi2_contingency(tab, correction=False)
    v = np.sqrt(chi2 / (tab.values.sum() * (min(tab.shape) - 1)))
    test = "χ²"
    if tab.shape == (2, 2):
        p = stats.fisher_exact(tab.values)[1]
        test = "Fisher exact"
    return "<h4>cluster × %s</h4>%s<p class='caption'>%s p = %s · Cramér's V = %.2f</p>" % (
        html.escape(col), _table(tab.reset_index()), test, _fmt(p, 2), v)


def _cluster_table(exp) -> str:
    rows = []
    df = exp.samples
    for j, c in enumerate(exp.cstats):
        s = df[(df["cluster"] == label_name(j)) & (df["used_in_fit"] == True)]  # noqa: E712
        rows.append({"cluster": label_name(j), "samples": c.n,
                     "female": int((s["inferred_sex"] == "female").sum()),
                     "male": int((s["inferred_sex"] == "male").sum()),
                     "trainable (n ≥ %d)" % exp.profile.n_min: "yes" if c.n >= exp.profile.n_min else "no",
                     "covariance estimate": c.method, "T² limit": c.t2_threshold})
    return _table(pd.DataFrame(rows))


def _samples_table(exp) -> str:
    cols = ["sample", "status", "cluster", "posterior", "margin", "t2_ratio", "q_ratio", "strong_downweight",
            "dropout_frac", "inferred_sex", "total_counts"] + exp.metadata_columns
    cols = [c for c in cols if c in exp.samples]
    return _table(exp.samples[cols].sort_values(["status", "cluster", "sample"]), "samples")


# --------------------------------------------------------------------------------------------
# page
# --------------------------------------------------------------------------------------------

CSS = """
.viz-root{--surface-1:%s;--text-primary:%s;--text-secondary:%s;--grid:%s;color-scheme:light;}
body{margin:0;background:var(--surface-1);color:var(--text-primary);
     font-family:system-ui,-apple-system,'Segoe UI',sans-serif;font-size:14px;line-height:1.45}
main{max-width:1180px;margin:0 auto;padding:24px 28px 64px}
h1{font-size:22px;margin:0 0 4px} h2{font-size:17px;margin:40px 0 4px;padding-top:12px;border-top:1px solid var(--grid)}
h4{margin:18px 0 6px;font-size:14px}
.caption{color:var(--text-secondary);margin:2px 0 10px;max-width:900px}
.tiles{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));gap:12px;margin:18px 0}
.tile{border:1px solid var(--grid);border-radius:8px;padding:10px 12px}
.tile .v{font-size:24px;font-weight:600} .tile .l{color:var(--text-secondary);font-size:12px}
.warn{border-left:3px solid #fab219;padding:6px 10px;margin:6px 0;background:#fffaf0}
table{border-collapse:collapse;margin:6px 0 4px;font-size:12.5px}
th,td{border-bottom:1px solid var(--grid);padding:4px 10px;text-align:left;white-space:nowrap}
th{color:var(--text-secondary);font-weight:600}
.scroll{overflow-x:auto} select{font-size:13px;padding:3px 6px;margin:6px 0}
nav a{color:var(--text-secondary);margin-right:12px;font-size:12.5px}
""" % (SURFACE, TEXT, TEXT2, GRID)

JS = """
function showColor(v){document.querySelectorAll('.splom').forEach(function(d){
  var on=d.dataset.var===v; d.style.display=on?'block':'none';
  if(on){d.querySelectorAll('.plotly-graph-div').forEach(function(g){Plotly.Plots.resize(g);});}});}
"""


def _section(anchor: str, title: str, caption: str, body: str) -> str:
    return '<h2 id="%s">%s</h2><p class="caption">%s</p>%s' % (anchor, title, caption, body)


def build_report(exp, path: str) -> None:
    from .explore import summary
    s = summary(exp)
    k = exp.clusters.K
    tiles = [("samples", s["n_samples"]), ("in fit", s["n_fit"]), ("features", s["n_features"]),
             ("PCs above noise", s["k_parallel_analysis"]), ("clusters K", k),
             ("stability ARI", _fmt(s["bootstrap_median_ari"], 2) if k > 1 else "n/a")]
    tiles_html = "".join('<div class="tile"><div class="v">%s</div><div class="l">%s</div></div>' % (v, l)
                         for l, v in tiles)
    warn_html = "".join('<div class="warn">%s</div>' % html.escape(w) for w in exp.warnings)

    color_vars = _color_variables(exp)
    options = "".join('<option value="%s">%s</option>' % (v, v) for v in color_vars)
    sploms = "".join('<div class="splom" data-var="%s" style="display:%s">%s</div>' % (
        v, "block" if v == "cluster" else "none", _div(fig_splom(exp, v))) for v in color_vars)

    parts = [
        _section("preflight", "1 · Preflight",
                 "Samples with many targets at ~zero depth were probably captured with a different design or "
                 "probe pool. They are excluded from the fit, then scored against it for comparison.",
                 _div(fig_dropout(exp))),
        _section("sex", "2 · Sex check",
                 "Sex inferred from chrX/chrY depth (PARs excluded). Clustering uses autosomes only, so sex "
                 "should not drive the clusters (see section 9).", _div(fig_sex(exp))),
        _section("scree", "3 · How many technical axes?",
                 "Bars above the grey noise line are PCs with more structure than random permutations of the "
                 "same data. %d PC(s) exceed it; %d are kept (min 2 for display)." % (exp.pca.k_pa, exp.pca.k),
                 _div(fig_scree(exp))),
        _section("pcs", "4 · Samples in PC space",
                 "Lower-triangle scatter matrix of the kept PCs, fit samples only. Recolour to check whether "
                 "clusters coincide with run, kit, sex or depth. Hover for sample details.",
                 'Colour by <select onchange="showColor(this.value)">%s</select>%s' % (options, sploms)),
        _section("bic", "5 · Choosing K",
                 "BIC per K and covariance type. Open markers violate the exploration minimum cluster size of "
                 "%d. The chosen K is the lowest valid BIC that also passes the subsampling stability check."
                 % exp.profile.explore_n_min, _div(fig_bic(exp))),
        _section("stability", "6 · Are the clusters stable?",
                 "Refits on 80% subsamples (without replacement). ARI near 1 and a block-diagonal co-clustering matrix mean the "
                 "clusters are reproducible rather than arbitrary splits.",
                 (_div(fig_ari(exp)) + _div(fig_coassign(exp))) if k > 1 else "<p>K = 1: not applicable.</p>"),
        _section("corr", "7 · Sample–sample correlation",
                 "Correlation of the normalized feature profiles, ordered by cluster. Independent of the PCA, so "
                 "it is a cross-check of the batch structure.", _div(fig_corr(exp))),
        _section("loadings", "8 · What drives each PC?",
                 "Loadings along the genome. A technical axis spreads over all chromosomes; a PC concentrated on "
                 "one chromosome is likely biology (a common CNV or an aneuploid sample) and is excluded from "
                 "clustering.", _div(fig_loadings(exp))),
        _section("bio", "9 · Is it biology?",
                 "Cluster composition against sex and any metadata columns. A strong association with sex or a "
                 "genotype is a red flag; association with run, kit or date supports a technical batch.",
                 _cluster_table(exp) + _contingency(exp, "inferred_sex")
                 + "".join(_contingency(exp, c) for c in exp.metadata_columns)),
        _section("gates", "10 · Out-of-distribution gates",
                 "T² measures distance inside the PC subspace (relative to the assigned cluster), Q the residual "
                 "outside it. Both are shown relative to their limits. T² for fit samples is in-sample and "
                 "therefore optimistic; Q is cross-fitted. Small clusters have very wide T² limits.",
                 _div(fig_gates(exp))),
        _section("posterior", "11 · Assignment confidence",
                 "Posterior probability of each fit sample's cluster. Samples below the line are AMBIGUOUS.",
                 _div(fig_posterior(exp)) if k > 1 else "<p>K = 1: all posteriors are 1.</p>"),
    ]
    if exp.gc_profile is not None:
        parts.append(_section("gc", "12 · GC-bias fingerprint",
                              "Median residual depth per GC bin. Clusters with different curve shapes differ in "
                              "GC bias, a typical library-prep or enrichment signature.", _div(fig_gc(exp))))
    parts.append(_section("table", "Samples", "Full per-sample results (also in samples.tsv).",
                          '<div class="scroll">%s</div>' % _samples_table(exp)))

    nav = "".join('<a href="#%s">%s</a>' % (a, t) for a, t in [
        ("preflight", "preflight"), ("pcs", "PCs"), ("bic", "K"), ("stability", "stability"), ("corr", "correlation"),
        ("loadings", "loadings"), ("bio", "biology"), ("gates", "gates"), ("table", "samples")])
    doc = """<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<title>gcnvrouter exploration · %s</title><style>%s</style><script>%s</script><script>%s</script></head>
<body class="viz-root"><main>
<h1>Batch exploration · %s profile</h1>
<p class="caption">gcnvrouter %s · generated %s · seed %d · %d count files</p>
<nav>%s</nav>
<div class="tiles">%s</div>%s
%s
</main></body></html>""" % (exp.profile.name, CSS, get_plotlyjs(), JS, exp.profile.name, s["gcnvrouter_version"],
                           datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), exp.seed, s["n_samples"], nav,
                           tiles_html, warn_html, "\n".join(parts))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(doc)
