#!/usr/bin/env python3
# 08_biotic_time_series_pub.py
# Pollen-only biotic time-series with period-aware sub-bins, rolling window, LOESS smoothing.
# Publication-ready figures: legends outside (inside reserved margin), subtle in-panel shading,
# and distinct colors for each geological period.

import os, math, json, warnings, random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from collections import deque

warnings.filterwarnings("ignore", category=UserWarning)

# ---------------- PATHS ----------------
BASE   = os.path.dirname(os.path.abspath(__file__))
XLSX   = os.path.join(BASE, "03_nordic_analysis", "nordic_pollen_data.xlsx")
OUTDIR = os.path.join(BASE, "05_biotic_time_series")

# ---------------- CONFIG ----------------
USECOLS = ['site_name','latitude','longitude','dataset_type','sample_id','age_BP','taxon','count']

EXCLUDE_TAXA = {
    'Cerealia','Secale','Hordeum','Triticum','Avena','Zea','Fagopyrum','Cannabis','Humulus',
    'Plantago','Plantago lanceolata','Plantago major',
    'Rumex','Rumex acetosa','Rumex acetosella','Centaurea cyanus',
    'Nymphaea','Potamogeton','Myriophyllum','Sparganium','Typha',
    'Pediastrum','Botryococcus','Concentricystes'
}

CAP_PCT = {'Pinus':40.0,'Betula':40.0,'Alnus':30.0,'Corylus':30.0,'Poaceae':30.0}

ARBOR_TAXA = ['Pinus','Picea','Betula','Alnus','Corylus','Quercus','Tilia','Ulmus','Fraxinus',
              'Fagus','Carpinus','Acer','Salix','Populus','Larix','Juniperus','Abies']
HERB_TAXA  = ['Poaceae','Artemisia','Chenopodiaceae','Amaranthaceae','Asteraceae','Cyperaceae',
              'Calluna','Ericaceae','Brassicaceae','Caryophyllaceae','Ranunculaceae','Rosaceae']

# Named periods with distinct colors - NON-OVERLAPPING, realistic data range
UPDATED_PERIODS = [
    ("Late Pleistocene",     11_700,    27_500,    "#d4b9da"),  # light purple
    ("Early Holocene",        8_200,    11_700,    "#fdd49e"),  # light orange
    ("Mid Holocene",          4_200,     8_200,    "#c7e9c0"),  # light green
    ("Late Holocene",         1_000,     4_200,    "#fff4e6"),  # very light orange
    ("Historical",              100,     1_000,    "#7ec8e3"),  # medium cyan
    ("Recent",                    0,       100,    "#b3e0f2"),  # light cyan
]

# Binning & smoothing
SUBBIN_WIDTH_YR     = 500
MIN_SAMPLES_PER_BIN = 30
ROLLING_WINDOW_YR   = 2000
ROLLING_STEP_YR     = 200
LOESS_FRAC          = 0.25

# Uncertainty ribbon (bootstrap across sites per sub-bin)
BOOTSTRAP_REPS = 1000
RANDOM_STATE   = 42

# Figure style
DPI = 300
FIGSIZE = (8.2, 4.8)
RIGHT_MARGIN = 0.72
PERIOD_SHADE_ALPHA = 0.35
SHOW_PERIOD_STRIP  = False
ADD_PERIOD_LEGEND  = True

# ---------------- UTILS ----------------
def seed_everything(s=RANDOM_STATE):
    np.random.seed(s); random.seed(s)

def safe_mkdir(d): os.makedirs(d, exist_ok=True)

def to_percent_matrix(df_counts):
    wide = df_counts.pivot_table(index='sample_id', columns='taxon', values='count',
                                 aggfunc='sum', fill_value=0.0)
    total = wide.sum(axis=1).replace(0, np.nan)
    pct = (wide.div(total, axis=0) * 100.0).fillna(0.0)
    return pct, list(pct.columns)

def apply_caps_then_renorm(wide_pct):
    wp = wide_pct.copy()
    for t, cap in CAP_PCT.items():
        if t in wp.columns:
            wp[t] = wp[t].clip(upper=cap)
    rs = wp.sum(axis=1).replace(0, np.nan)
    return (wp.div(rs, axis=0)*100.0).fillna(0.0)

def richness_from_pct(row, thresh=0.5): return float((row >= thresh).sum())

def shannon_H_from_pct(row):
    p = np.clip(row.values/100.0,0,1); p = p[p>0]
    return float(-(p*np.log(p)).sum()) if p.size else np.nan

def pielou_evenness(H, S):
    return float(H/np.log(S)) if S and S>1 and not np.isnan(H) else np.nan

def openness_from_pct(row):
    ap = sum(row.get(t,0.0) for t in ARBOR_TAXA if t in row.index)
    if any(t in row.index for t in HERB_TAXA):
        nap = sum(row.get(t,0.0) for t in HERB_TAXA if t in row.index)
        denom = ap + nap
        if denom <= 0:
            nap = max(0.0, 100.0-ap); denom = ap+nap
    else:
        nap = max(0.0, 100.0-ap); denom = ap+nap
    return float(nap/denom) if denom>0 else np.nan

def chord_distance(a,b):
    va = np.asarray(a,float); vb = np.asarray(b,float)
    if va.size==0 or vb.size==0: return np.nan
    va = va/(np.linalg.norm(va)+1e-12)
    vb = vb/(np.linalg.norm(vb)+1e-12)
    return float(np.linalg.norm(va-vb))

def assign_period(age_bp):
    for item in UPDATED_PERIODS:
        name, lo, hi = item[0], item[1], item[2]
        if lo <= age_bp < hi: return name
    return "Outside range"

def get_period_color(period_name):
    for item in UPDATED_PERIODS:
        if item[0] == period_name:
            return item[3] if len(item) > 3 else "#cccccc"
    return "#cccccc"

def build_period_subbins(periods, subw, max_age):
    edges = list(range(0, int(max_age)+subw, subw))
    out = []
    for i in range(len(edges)-1):
        lo, hi = edges[i], edges[i+1]
        label = assign_period((lo+hi)/2.0)
        out.append((lo,hi,label))
    merged = []
    for lo,hi,label in out:
        if not merged or merged[-1][2]!=label:
            merged.append([lo,hi,label])
        else:
            merged[-1][1]=hi
    return [(a,b,l) for a,b,l in merged]

def count_bin_samples(df_samples, lo, hi):
    sub = df_samples[(df_samples['age_BP']>=lo) & (df_samples['age_BP']<hi)]
    return int(sub['sample_id'].nunique())

def merge_small_bins(bins, df_samples, min_n):
    dq = deque(bins)
    merged = []
    while dq:
        lo,hi,label = dq.popleft()
        n = count_bin_samples(df_samples, lo, hi)
        if n >= min_n or not dq:
            merged.append([lo,hi,label])
        else:
            cur_lo, cur_hi, cur_label = lo, hi, label
            n_acc = n
            while dq and n_acc < min_n:
                lo2,hi2,label2 = dq.popleft()
                cur_hi = hi2
                cur_label = cur_label if cur_label==label2 else f"{cur_label} / {label2}"
                n_acc = count_bin_samples(df_samples, cur_lo, cur_hi)
            merged.append([cur_lo,cur_hi,cur_label])
    return [(a,b,l) for a,b,l in merged]

def loess_1d(x, y, frac=0.25, iters=1):
    x = np.asarray(x, float); y = np.asarray(y, float)
    n = len(x)
    if n==0: return np.array([])
    r = max(2, int(math.ceil(frac * n)))
    yfit = np.empty(n); yfit[:] = np.nan
    robust = np.ones(n)
    for _ in range(max(1,iters)):
        for i in range(n):
            dist = np.abs(x - x[i])
            idx = np.argsort(dist)[:r]
            xw, yw, rw = x[idx], y[idx], robust[idx]
            dmax = dist[idx[-1]] if dist[idx[-1]]>0 else 1.0
            w = (1 - (dist[idx]/dmax)**3)**3
            w *= rw
            X = np.vstack([np.ones_like(xw), xw]).T
            W = np.diag(w)
            try:
                beta = np.linalg.pinv(X.T @ W @ X) @ (X.T @ W @ yw)
                yfit[i] = beta[0] + beta[1]*x[i]
            except Exception:
                yfit[i] = np.nan
        resid = y - yfit
        s = np.nanmedian(np.abs(resid)) + 1e-12
        robust = 1.0 / (1.0 + (resid/(6.0*s))**2)
    return yfit

def bootstrap_iqr(values, reps=BOOTSTRAP_REPS):
    vals = np.asarray(values, float)
    vals = vals[~np.isnan(vals)]
    if len(vals) == 0 or reps <= 0:
        return np.nan, np.nan
    rng = np.random.default_rng(RANDOM_STATE)
    qs25, qs75 = [], []
    for _ in range(reps):
        samp = vals[rng.integers(0, len(vals), len(vals))]
        qs25.append(np.percentile(samp, 25))
        qs75.append(np.percentile(samp, 75))
    return float(np.nanmedian(qs25)), float(np.nanmedian(qs75))

def get_visible_periods(xlim, periods=UPDATED_PERIODS):
    """Return only periods that are visible in the current x-axis range"""
    xmin, xmax = sorted(xlim)
    visible = []
    for item in periods:
        name, lo, hi = item[0], item[1], item[2]
        color = item[3] if len(item) > 3 else "#cccccc"
        if lo < hi and not (hi < xmin or lo > xmax):
            visible.append((name, lo, hi, color))
    return visible

def draw_period_strip(fig, ax_main, periods=UPDATED_PERIODS):
    pos = ax_main.get_position()
    h = pos.height * 0.12
    gap = h * 0.25
    ax_strip = fig.add_axes([pos.x0, pos.y1 + gap, pos.width, h])
    ax_strip.set_xlim(ax_main.get_xlim())
    ax_strip.set_ylim(0, 1)
    ax_strip.axis("off")
    for i, item in enumerate(periods):
        name, lo, hi = item[0], item[1], item[2]
        color = item[3] if len(item) > 3 else "#cccccc"
        if lo >= hi: continue
        ax_strip.add_patch(Rectangle((lo, 0), hi-lo, 1,
                                     facecolor=color,
                                     alpha=PERIOD_SHADE_ALPHA,
                                     edgecolor='none'))
        xm = (lo + hi) / 2.0
        ax_strip.text(xm, 0.5, name, ha='center', va='center', fontsize=8)
    if ax_main.xaxis_inverted():
        ax_strip.invert_xaxis()

# ---------------- MAIN ----------------
def main():
    seed_everything()
    safe_mkdir(OUTDIR)
    if not os.path.exists(XLSX):
        raise FileNotFoundError(f"Cannot find Excel at: {XLSX}")

    df = pd.read_excel(XLSX, usecols=USECOLS)
    df = df.dropna(subset=['sample_id','taxon','count'])
    df['count'] = pd.to_numeric(df['count'], errors='coerce').fillna(0)
    df = df[df['count']>0].copy()
    df = df[~df['taxon'].isin(EXCLUDE_TAXA)].copy()

    ages = (df.groupby('sample_id', as_index=False)['age_BP']
              .median().rename(columns={'age_BP':'age_BP'}))
    meta = (df[['sample_id','site_name','latitude','longitude']]
            .drop_duplicates('sample_id')
            .merge(ages, on='sample_id', how='left'))

    wide_pct, taxa_list = to_percent_matrix(df)
    wide_pct = apply_caps_then_renorm(wide_pct)
    wide_pct = wide_pct.loc[meta['sample_id'].values]
    wide_pct.index.name = 'sample_id'

    recs = []
    for sid, row in wide_pct.iterrows():
        S  = richness_from_pct(row, 0.5)
        H  = shannon_H_from_pct(row)
        J  = pielou_evenness(H, S)
        OP = openness_from_pct(row)
        recs.append((sid,S,H,J,OP,float(row.sum())))
    per_sample = (meta.merge(pd.DataFrame(recs, columns=['sample_id','richness','shannon_H','evenness','openness','pollen_sum']),
                             on='sample_id', how='left')
                       .sort_values(['site_name','age_BP']))

    per_sample['roc_chord'] = np.nan
    for site, sub in per_sample.groupby('site_name', observed=False):
        sids = sub['sample_id'].values
        Xs = wide_pct.loc[sids].values
        d = [np.nan]
        for i in range(1, Xs.shape[0]):
            d.append(chord_distance(Xs[i-1], Xs[i]))
        per_sample.loc[sub.index, 'roc_chord'] = d
    per_sample.to_csv(os.path.join(OUTDIR, "per_sample_metrics.csv"), index=False)

    max_age = int(np.nanmax(per_sample['age_BP'])) if per_sample['age_BP'].notna().any() else 500
    raw_bins = build_period_subbins(UPDATED_PERIODS, SUBBIN_WIDTH_YR, max_age)
    bins_ok  = merge_small_bins(raw_bins, per_sample, MIN_SAMPLES_PER_BIN)

    mids = [ (lo+hi)/2.0 for lo,hi,_ in bins_ok ]
    iv   = pd.IntervalIndex.from_tuples([(lo,hi) for lo,hi,_ in bins_ok], closed='left')
    per_sample['bin_idx'] = iv.get_indexer(per_sample['age_BP'])
    per_sample = per_sample[per_sample['bin_idx']>=0].copy()
    per_sample['age_bin']     = per_sample['bin_idx'].apply(lambda i: iv[i])
    per_sample['age_mid']     = per_sample['bin_idx'].apply(lambda i: mids[i])
    per_sample['period_name'] = per_sample['bin_idx'].apply(lambda i: bins_ok[i][2])

    site_bin = (per_sample
        .groupby(['site_name','latitude','longitude','age_bin','age_mid','period_name'], observed=False)
        [['richness','shannon_H','evenness','openness','roc_chord']]
        .mean().reset_index())
    site_bin.to_csv(os.path.join(OUTDIR, "site_bin_metrics_period_subbins.csv"), index=False)

    agg = (site_bin
        .groupby(['age_bin','age_mid','period_name'], observed=False)
        .agg(
            n_sites=('site_name','nunique'),
            richness_med=('richness','median'),
            shannon_med=('shannon_H','median'),
            evenness_med=('evenness','median'),
            openness_med=('openness','median'),
            roc_med=('roc_chord','median')
        )
        .reset_index().sort_values('age_mid')
    )
    if BOOTSTRAP_REPS > 0:
        q25, q75 = [], []
        for (a,b,p) in agg[['age_bin','age_mid','period_name']].itertuples(index=False, name=None):
            vals = site_bin.loc[site_bin['age_bin']==a, 'openness'].values
            lo, hi = bootstrap_iqr(vals, reps=BOOTSTRAP_REPS)
            q25.append(lo); q75.append(hi)
        agg['openness_q25'] = q25
        agg['openness_q75'] = q75
    agg.to_csv(os.path.join(OUTDIR, "regional_bins_period_subbins.csv"), index=False)

    centers = list(range(0, max_age+1, ROLLING_STEP_YR))
    half = ROLLING_WINDOW_YR/2.0
    rows = []
    for c in centers:
        lo,hi = c-half, c+half
        sub = per_sample[(per_sample['age_BP']>=lo) & (per_sample['age_BP']<hi)]
        if len(sub)==0:
            rows.append((c,0,np.nan,np.nan,np.nan,np.nan)); continue
        bysite = (sub.groupby('site_name', observed=False)
                    [['richness','shannon_H','evenness','openness','roc_chord']].median())
        rows.append((c, bysite.shape[0],
                     bysite['richness'].median(),
                     bysite['shannon_H'].median(),
                     bysite['evenness'].median(),
                     bysite['openness'].median()))
    rolling = pd.DataFrame(rows, columns=['center_BP','n_sites','richness_med','shannon_med','evenness_med','openness_med'])
    rolling = rolling.sort_values('center_BP')
    rolling.to_csv(os.path.join(OUTDIR, "rolling_window_medians.csv"), index=False)

    x = agg['age_mid'].values
    y = agg['openness_med'].values
    mask = ~np.isnan(y)
    loess_y = np.full_like(y, np.nan, dtype=float)
    if mask.sum() >= 5:
        loess_y[mask] = loess_1d(x[mask], y[mask], frac=LOESS_FRAC, iters=2)
    agg['openness_loess'] = loess_y
    agg.to_csv(os.path.join(OUTDIR, "regional_with_loess.csv"), index=False)

    # ---- plotting helpers ----
    def shade_periods(ax):
        for item in UPDATED_PERIODS:
            name, lo, hi = item[0], item[1], item[2]
            color = item[3] if len(item) > 3 else "#cccccc"
            if lo < hi:
                ax.axvspan(lo, hi, alpha=PERIOD_SHADE_ALPHA, color=color, linewidth=0)

    def finalize(fig, ax, xlabel, ylabel, add_strip=True):
        ax.invert_xaxis()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        
        xlim = ax.get_xlim()
        visible_periods = get_visible_periods(xlim)
        
        fig.subplots_adjust(right=RIGHT_MARGIN)
        
        handles, labels = ax.get_legend_handles_labels()
        
        if ADD_PERIOD_LEGEND and visible_periods:
            handles.append(Patch(facecolor='none', edgecolor='none'))
            labels.append('')
            
            for name, lo, hi, color in visible_periods:
                period_patch = Patch(facecolor=color, alpha=PERIOD_SHADE_ALPHA, 
                                    edgecolor='gray', linewidth=0.5)
                handles.append(period_patch)
                labels.append(f'{name}: {lo:,}-{hi:,} BP')
        
        if handles:
            ax.legend(handles, labels, loc='upper left',
                      bbox_to_anchor=(1.01, 1.0), borderaxespad=0., frameon=False,
                      fontsize=8.5)
        
        if add_strip and SHOW_PERIOD_STRIP:
            draw_period_strip(fig, ax)

    # ---- Plot 1: Openness ----
    fig, ax = plt.subplots(figsize=FIGSIZE)
    shade_periods(ax)
    if 'openness_q25' in agg.columns:
        ax.fill_between(agg['age_mid'], agg['openness_q25'], agg['openness_q75'],
                        alpha=0.25, label='IQR (bootstrap)')
    ax.plot(agg['age_mid'], agg['openness_med'], lw=1.8, label='Median across sites')
    if np.isfinite(agg['openness_loess']).any():
        ax.plot(agg['age_mid'], agg['openness_loess'], lw=2.0, label=f'LOESS (frac={LOESS_FRAC})')
    finalize(fig, ax, "Age (years BP)", "Openness (NAP / (AP+NAP))")
    fig.savefig(os.path.join(OUTDIR, "plot_openness_period_loess.png"), dpi=DPI, bbox_inches='tight')
    plt.close(fig)

    # ---- Plot 2: Shannon H′ ----
    fig, ax = plt.subplots(figsize=FIGSIZE)
    shade_periods(ax)
    ax.plot(agg['age_mid'], agg['shannon_med'], lw=1.8, label='Median across sites')
    ax.plot(rolling['center_BP'], rolling['shannon_med'], lw=1.2, alpha=0.9,
            label=f'Rolling median ({ROLLING_WINDOW_YR} yr)')
    finalize(fig, ax, "Age (years BP)", "Shannon H′")
    fig.savefig(os.path.join(OUTDIR, "plot_shannon_period_rolling.png"), dpi=DPI, bbox_inches='tight')
    plt.close(fig)

    # ---- Plot 3: Rate of change ----
    fig, ax = plt.subplots(figsize=FIGSIZE)
    shade_periods(ax)
    ax.plot(agg['age_mid'], agg['roc_med'], lw=1.8, label='Median across sites')
    finalize(fig, ax, "Age (years BP)", "Rate of change (chord distance)")
    fig.savefig(os.path.join(OUTDIR, "plot_rate_of_change_period.png"), dpi=DPI, bbox_inches='tight')
    plt.close(fig)

    meta_out = {
        "input_xlsx": XLSX,
        "n_samples": int(per_sample['sample_id'].nunique()),
        "n_sites":   int(per_sample['site_name'].nunique()),
        "n_taxa":    int(len(taxa_list)),
        "subbin_width_years": SUBBIN_WIDTH_YR,
        "min_samples_per_bin": MIN_SAMPLES_PER_BIN,
        "rolling_window_years": ROLLING_WINDOW_YR,
        "rolling_step_years":   ROLLING_STEP_YR,
        "loess_frac": LOESS_FRAC,
        "bootstrap_reps": BOOTSTRAP_REPS,
        "periods_used": UPDATED_PERIODS,
    }
    with open(os.path.join(OUTDIR, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(meta_out, f, indent=2)

    print("Done. Results written to:", OUTDIR)

if __name__ == "__main__":
    main()