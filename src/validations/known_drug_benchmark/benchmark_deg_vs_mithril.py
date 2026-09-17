#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Known-drug positioning: where the known ALS/IPF drugs rank among all drugs,
under each disease-signature representation.

The analysis is descriptive: it reports where a small set of known drugs falls
in the drug ranking under each representation, and tests no hypothesis.

DESIGN
------
Held fixed across every configuration:
  * CS_METHOD  = 'bin_chen'
  * rank_on    = 'magnitude'
  * drug input = the same LMM-collapsed LINCS L1000 drug set (3,221 drugs
                 x 3 timepoints)
  * disease    = one at a time (als_NYGC = GEO GSE153960; ipf)
The only thing that varies is the disease-signature representation:
  DEG      gene-level differential expression, all genes      (8,409 / 8,669)
  DEG_LM   gene-level differential expression, landmark only  (826 / 862)
  mith     MITHrIL perturbation, all nodes                    (15,825)
  mith_LM  MITHrIL perturbation, landmark only                (763)
  mith_pw  MITHrIL perturbation aggregated to pathways        (2,159)
DEG / DEG_LM are the differential-expression family; mith / mith_LM / mith_pw
are the network-enhanced (MITHrIL) family.

SIGN CONVENTION
---------------
More NEGATIVE connectivity_score = stronger reversal = better candidate. The
runs listed in RUN_FILES all use this convention; check_sign_convention()
asserts corr(connectivity_score, gene-level pearson) > 0 on every run before
anything is reported.

DRUG ROSTER TIERING (labels only -- no statistic uses the tiers)
----------------------------------------------------------------
  Tier A  approved and disease-modifying
  Tier B  approved but purely symptomatic
  Tier C  failed phase 3 / withdrawn
Drugs relevant to the diseases but absent from this drug set are listed as
comments next to ROSTER so the coverage gap is explicit.

INPUT
-----
The connectivity-score TSVs named in RUN_FILES, read from
connectivity_score/output/ under the repo root that resolve_base() locates.

OUTPUT
------
Percentile/rank tables to the results output directory, and two figures (a
percentile heatmap and per-disease score distributions) to the image output
directory. Both directories come from conf.py; nothing is written under src/.

USAGE
-----
    python -m src.validations.known_drug_benchmark.benchmark_deg_vs_mithril
    python benchmark_deg_vs_mithril.py --base <repo root>
From Spyder: run the file. With no arguments the repo root comes from
conf.BASE_DIR (src/local.py); see resolve_base() for the full fallback chain.
"""
import argparse, re, sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
# Interactive sessions (Spyder / IPython) keep their own backend so the
# figures appear inline as well as being written to disk; a plain
# command-line run gets the headless one.
_INTERACTIVE = hasattr(sys, 'ps1') or ('IPython' in sys.modules)
if not _INTERACTIVE:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# --- repo bootstrap ---------------------------------------------------------
# Same convention as src/validations/binchen2017/test_sRGES_replication.py:
# put .../src on sys.path so `import conf` works however the file is invoked
# (as a script, with -m, or from Spyder).
try:
    THIS_DIR = Path(__file__).resolve().parent
except NameError:                       # pasted into a console: no __file__
    THIS_DIR = Path.cwd()
SRC_DIR = THIS_DIR.parent.parent        # .../src
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

# ---------------------------------------------------------------- roster ----
# Tier A: approved + disease-modifying. Tier B: approved, symptomatic only.
# Tier C: failed phase 3 / withdrawn. Labels only -- no statistic uses them.
ROSTER = {
    'als_NYGC': {
        'riluzole':         ('A', 'approved 1995, disease-modifying (survival benefit)'),
        'dextromethorphan': ('B', 'Nuedexta component; treats pseudobulbar affect only'),
        'quinidine':        ('B', 'Nuedexta component; CYP2D6 inhibitor, not ALS-directed'),
        'phenylbutyrate':   ('C', 'Relyvrio = phenylbutyrate + taurursodiol; only THIS '
                                   'component is scored (taurursodiol absent from the drug '
                                   'set); PHOENIX ph3 failed, withdrawn 2024. Caveat: "failed '
                                   'clinically" != "does not reverse the signature", and only '
                                   'half the failed combination is represented here.'),
    },
    'ipf': {
        'pirfenidone':      ('A', 'approved 2014, antifibrotic, pleiotropic (TGF-beta etc.)'),
        'nintedanib':       ('A', 'approved 2014; oncology angiokinase inhibitor repurposed'),
    },
}
# Known-but-unavailable, recorded so the coverage gap is explicit and auditable:
#   ALS  edaravone (approved, disease-modifying) -- absent from this LINCS set
#        tofersen (approved, ASO)                -- not a small molecule
#        taurursodiol (Relyvrio component)       -- absent
#   IPF  nerandomilast (Jascayd, approved 2025)  -- postdates LINCS L1000

CONFIGS = ['DEG', 'DEG_LM', 'mith', 'mith_LM', 'mith_pw']
DEG_FAMILY  = ['DEG', 'DEG_LM']
MITH_FAMILY = ['mith', 'mith_LM', 'mith_pw']

RUN_FILES = {
    'als_NYGC': {
        'DEG':     'als_NYGC/30_06_2026_11_55_DEG_bin_chen_connectivity_score.tsv',
        'DEG_LM':  'als_NYGC/30_06_2026_12_02_DEG_LM_bin_chen_connectivity_score.tsv',
        'mith':    'als_NYGC/30_06_2026_12_03_mith_bin_chen_connectivity_score.tsv',
        'mith_pw': 'als_NYGC/30_06_2026_12_08_mith_pw_bin_chen_connectivity_score.tsv',
        'mith_LM': 'als_NYGC/30_06_2026_12_09_mith_LM_bin_chen_connectivity_score.tsv',
    },
    'ipf': {
        'DEG':     'ipf/30_06_2026_10_49_DEG_bin_chen_connectivity_score.tsv',
        'DEG_LM':  'ipf/30_06_2026_13_00_DEG_LM_bin_chen_connectivity_score.tsv',
        'mith':    'ipf/30_06_2026_13_00_mith_bin_chen_connectivity_score.tsv',
        'mith_pw': 'ipf/21_06_2026_17_28_mith_pw_bin_chen_connectivity_score.tsv',
        'mith_LM': 'ipf/30_06_2026_13_05_mith_LM_bin_chen_connectivity_score.tsv',
    },
}
TIMEPOINTS = ('6h_24h', '24h', '6h')


# ------------------------------------------------------- base-dir resolution -
def _looks_like_repo(p):
    """A NetCoS root, for this script's purposes, is one that holds the CS
    outputs the benchmark reads."""
    try:
        return (Path(p) / 'connectivity_score' / 'output').is_dir()
    except OSError:
        return False


def base_dir_from_conf(verbose=True):
    """conf.BASE_DIR (which conf.py takes from src/local.py).

    Returns None rather than raising, so a missing or machine-specific
    local.py degrades to the next candidate instead of killing the run.
    """
    try:
        import conf
    except Exception as exc:
        if verbose:
            print(f'[base] conf.py not importable '
                  f'({exc.__class__.__name__}: {exc})')
        return None
    b = getattr(conf, 'BASE_DIR', None)
    if b is None:
        if verbose:
            print('[base] conf.py imported but defines no BASE_DIR')
        return None
    return Path(b)


def resolve_base(cli_base=None, verbose=True):
    """Resolve the repo root, in priority order:

      1. ``--base`` on the command line -- explicit, and validated strictly:
         if it is given and wrong, that is an error, not a reason to guess.
      2. ``conf.BASE_DIR`` (from src/local.py) -- the machine-local setting.
         This is the one that makes a bare "Run file" in Spyder work.
      3. This file's own location (<root>/src/validations/known_drug_benchmark).
      4. The current working directory.

    Candidates 2-4 are each checked for connectivity_score/output before being
    accepted, so a local.py pointing somewhere stale falls through with a note
    instead of failing later inside read_csv with a confusing path. That is not
    hypothetical: local.py may say 'G:/My Drive/...' while the same Drive
    mounts as 'G:/Il mio Drive/...' on an Italian-locale machine.
    """
    if cli_base:
        p = Path(cli_base)
        if not _looks_like_repo(p):
            raise RuntimeError(
                f'--base was given as {p!s} but there is no '
                f'connectivity_score/output under it. Fix the path, or drop '
                f'--base to fall back to conf.BASE_DIR.')
        if verbose:
            print(f'[base] using --base: {p}')
        return p

    candidates = []
    c = base_dir_from_conf(verbose=verbose)
    if c is not None:
        candidates.append(('conf.BASE_DIR', c))
    if len(THIS_DIR.parents) >= 3:
        candidates.append(("this file's location", THIS_DIR.parents[2]))
    candidates.append(('cwd', Path.cwd()))

    for name, p in candidates:
        if _looks_like_repo(p):
            if verbose:
                print(f'[base] using {name}: {p}')
            return Path(p)
        if verbose:
            print(f'[base] skipped {name} (no connectivity_score/output): {p}')

    tried = '\n  '.join(f'{n}: {p}' for n, p in candidates)
    raise RuntimeError(
        'Could not locate the NetCoS repo root. Tried:\n  ' + tried +
        '\nEither fix BASE_DIR in src/local.py or pass --base explicitly.')


# ---------------------------------------------------- output-dir resolution -
# src/ holds code only. Images go to imgs/, every other result artifact to
# other_outputs/, under one canonical slug mirrored across the two output trees.
# Paths come from conf.py, never from a literal in this file.
ANALYSIS_SLUG = 'known_drug_benchmark'


def output_dirs(base, verbose=True):
    """Return (tables_dir, images_dir) for this analysis.

    Taken from conf.OOUT_DIR / conf.IMG_DIR when conf is importable; otherwise
    derived from the resolved repo root using the same layout, so the script
    still writes to the right place on a machine where conf cannot load.
    Nothing is ever written under src/.
    """
    img = oout = None
    try:
        import conf
        img = getattr(conf, 'IMG_DIR', None)
        oout = getattr(conf, 'OOUT_DIR', None)
    except Exception as exc:
        if verbose:
            print(f'[out] conf.py not importable for paths '
                  f'({exc.__class__.__name__}); deriving from the repo root')
    img = Path(img) if img else Path(base) / 'imgs'
    oout = Path(oout) if oout else Path(base) / 'other_outputs'
    return (oout / 'validations' / ANALYSIS_SLUG,
            img / 'validations' / ANALYSIS_SLUG)


def normalize(name):
    """Same normalization as exploration_utils / load_cs drug_key."""
    return re.sub(r'[^a-z0-9]', '', str(name).lower())


def load_run(base, disease, config):
    df = pd.read_csv(Path(base) / 'connectivity_score' / 'output' / RUN_FILES[disease][config],
                     sep='\t')
    suffix = '|'.join(TIMEPOINTS)
    df['drug'] = df['pert_id'].str.replace(rf'_({suffix})$', '', regex=True)
    df['timepoint'] = df['pert_id'].str.extract(rf'_({suffix})$')[0]
    df['drug_key'] = df['drug'].map(normalize)
    return df


def to_drug_level(df, collapse):
    """Collapse the 3 per-drug timepoint signatures to one score per drug.

    collapse='best' -> min (most negative) ; 'mean' -> mean.
    Percentile is ascending in CS, so LOW percentile = strong reverser.
    """
    g = df.groupby('drug_key')['connectivity_score']
    s = (g.min() if collapse == 'best' else g.mean()).rename('cs')
    out = s.reset_index()
    out['percentile'] = out['cs'].rank(pct=True) * 100
    out['rank'] = out['cs'].rank(method='min').astype(int)
    out['n_drugs'] = len(out)
    return out.set_index('drug_key')


def check_sign_convention(base):
    """A drug that reverses the disease signature must have BOTH a negative CS
    and a negative gene-level correlation, so corr(CS, pearson) must be > 0.

    This is a cheap correctness guard, not an inferential result: it protects
    the whole positioning read-out from a silent sign flip. Kept for that
    reason even though the enrichment analysis that used to follow it is gone.
    """
    rows = []
    for disease in RUN_FILES:
        for config in CONFIGS:
            d = load_run(base, disease, config)
            rows.append({'disease': disease, 'config': config,
                         'corr_CS_pearson': d['connectivity_score'].corr(d['pearson']),
                         'cs_min': d['connectivity_score'].min(),
                         'cs_max': d['connectivity_score'].max(),
                         'n_drugs': d['drug_key'].nunique()})
    res = pd.DataFrame(rows)
    bad = res[res['corr_CS_pearson'] <= 0]
    if len(bad):
        raise AssertionError(f"SIGN CONVENTION VIOLATED in:\n{bad}")
    # Also guard the cross-config comparison: every config must rank the SAME
    # drug set, or a percentile in one is not comparable to a percentile in
    # another. Warn rather than raise, so an odd run is visible but not fatal.
    per_disease = res.groupby('disease')['n_drugs'].nunique()
    if (per_disease > 1).any():
        print('[warn] drug-set size differs across configs within a disease; '
              'cross-config percentiles are only loosely comparable:\n'
              + res[['disease', 'config', 'n_drugs']].to_string(index=False))
    return res


def build_table(base, collapse):
    rows = []
    for disease, drugs in ROSTER.items():
        levels = {c: to_drug_level(load_run(base, disease, c), collapse) for c in CONFIGS}
        for drug, (tier, note) in drugs.items():
            key = normalize(drug)
            rec = {'disease': disease, 'drug': drug, 'tier': tier, 'note': note}
            for c in CONFIGS:
                dl = levels[c]
                rec[f'{c}_pct'] = round(dl.loc[key, 'percentile'], 1) if key in dl.index else np.nan
                rec[f'{c}_rank'] = int(dl.loc[key, 'rank']) if key in dl.index else np.nan
                rec[f'{c}_n'] = int(dl['n_drugs'].iloc[0])
            rows.append(rec)
    return pd.DataFrame(rows)


# ============================================================== plotting ====
# Palette: validated categorical/sequential steps (blue slot 1, orange slot 2).
# Sequential blue is used light->dark for magnitude; the heatmap deliberately
# maps LOW percentile (= strong reverser = good) to the DARK end, and says so on
# the colourbar. Figures render on the light surface only (static PNG).
SURFACE   = '#fcfcfb'
INK       = '#0b0b0b'
INK_2     = '#52514e'
MUTED     = '#898781'
GRID      = '#e1e0d9'
AXIS      = '#c3c2b7'
BLUE      = '#2a78d6'   # categorical slot 1
BLUE_LT   = '#86b6ef'   # blue 250
BLUE_DK   = '#1c5cab'   # blue 550
ORANGE    = '#eb6834'   # categorical slot 2
BAR_MUTED = '#c3c2b7'
BLUE_RAMP = ['#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7',
             '#3987e5', '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281',
             '#0d366b']
DISEASE_LABEL = {'als_NYGC': 'ALS', 'ipf': 'IPF'}


def _style():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['DejaVu Sans', 'Segoe UI', 'Helvetica', 'Arial'],
        'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE,
        'savefig.facecolor': SURFACE,
        'text.color': INK, 'axes.labelcolor': INK_2,
        'xtick.color': MUTED, 'ytick.color': MUTED,
        'axes.edgecolor': AXIS, 'axes.linewidth': 0.8,
        'grid.color': GRID, 'grid.linewidth': 0.8, 'grid.linestyle': '-',
        'axes.grid': False, 'axes.spines.top': False, 'axes.spines.right': False,
        'font.size': 9, 'axes.titlesize': 10, 'legend.frameon': False,
        'figure.dpi': 150,
    })


def _save(fig, path):
    """Write the figure, and close it only when nothing is watching.

    Under Spyder's inline backend a closed figure never gets rendered, so in
    an interactive session the figure is left open and shows up in the Plots
    pane as well as on disk.
    """
    fig.savefig(path, bbox_inches='tight')
    if not _INTERACTIVE:
        plt.close(fig)


def _row_order(table):
    """Tier A first, then B, then C - the reading order of the argument."""
    order = {'A': 0, 'B': 1, 'C': 2}
    t = table.copy()
    t['_o'] = t['tier'].map(order)
    return t.sort_values(['_o', 'disease', 'drug']).drop(columns='_o')


def plot_percentile_heatmap(ax, table):
    """The whole positioning result in one grid: every drug x every config."""
    t = _row_order(table)
    M = t[[f'{c}_pct' for c in CONFIGS]].values.astype(float)
    cmap = LinearSegmentedColormap.from_list('bluerev', list(reversed(BLUE_RAMP)))
    im = ax.imshow(M, cmap=cmap, vmin=0, vmax=100, aspect='auto')

    ax.set_xticks(range(len(CONFIGS)))
    ax.set_xticklabels(CONFIGS, color=INK_2)
    ax.set_yticks(range(len(t)))
    ax.set_yticklabels([f"{r.drug}\n({r.tier}) {DISEASE_LABEL[r.disease]}"
                        for r in t.itertuples()], color=INK_2)
    ax.set_xticks(np.arange(-.5, len(CONFIGS), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(t), 1), minor=True)
    # 2px surface gap between cells, not borders around marks
    ax.grid(which='minor', color=SURFACE, linewidth=2)
    ax.tick_params(which='minor', length=0)
    ax.tick_params(which='major', length=0)
    for s in ax.spines.values():
        s.set_visible(False)

    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            ax.text(j, i, f'{v:.1f}', ha='center', va='center', fontsize=8.5,
                    color='#ffffff' if v < 45 else INK,
                    fontweight='bold' if v < 20 else 'normal')

    # family brackets
    ax.text(0.5, -1.05, 'differential expression', ha='center', va='bottom',
            fontsize=8.5, color=MUTED)
    ax.text(3.0, -1.05, 'NetCoS (MITHrIL)', ha='center', va='bottom',
            fontsize=8.5, color=MUTED)
    ax.plot([-0.4, 1.4], [-0.75, -0.75], color=AXIS, lw=0.8, clip_on=False)
    ax.plot([1.6, 4.4], [-0.75, -0.75], color=AXIS, lw=0.8, clip_on=False)

    cb = ax.figure.colorbar(im, ax=ax, fraction=0.032, pad=0.02)
    cb.set_label('percentile among ~3,221 drugs\n(dark = low = stronger reverser)',
                 fontsize=8, color=INK_2)
    cb.ax.tick_params(labelsize=8, color=MUTED)
    cb.outline.set_visible(False)
    ax.set_title('Where each known drug ranks, by disease-signature representation',
                 loc='left', color=INK, pad=26)


def plot_distributions(base, table, collapse, out_path):
    """Where the known drugs sit inside the score distribution.

    Faceted disease x config. The two diseases are NEVER pooled: each drug's
    percentile is computed within its own disease, so the histogram it is drawn
    against must be that same disease's.
    """
    _style()
    diseases = list(ROSTER.keys())
    fig, axes = plt.subplots(len(CONFIGS), len(diseases),
                             figsize=(12.4, 11.0), squeeze=False)
    for row, cfg in enumerate(CONFIGS):
        for col, disease in enumerate(diseases):
            ax = axes[row][col]
            dl = to_drug_level(load_run(base, disease, cfg), collapse)
            ax.hist(dl['cs'].values, bins=70, color=BLUE_LT,
                    histtype='stepfilled', edgecolor=BLUE, linewidth=0.8)
            top = ax.get_ylim()[1]
            ax.set_ylim(0, top * 1.30)     # headroom for the labels

            marked = table[(table['disease'] == disease) &
                           (table['tier'].isin(['A', 'C']))]
            marked = marked.assign(cs_val=[dl.loc[normalize(d), 'cs']
                                        for d in marked['drug']]).sort_values('cs_val')
            for k, r in enumerate(marked.itertuples()):
                colour = ORANGE if r.tier == 'A' else MUTED
                ax.axvline(r.cs_val, color=colour, lw=2, zorder=4)
                # stagger heights so neighbouring markers never overprint
                yl = top * (1.22 if k % 2 == 0 else 1.05)
                lo, hi = ax.get_xlim()
                right_half = r.cs_val > (lo + hi) / 2
                ax.annotate(f'{r.drug} ({getattr(r, f"{cfg}_pct"):.0f}%)',
                            xy=(r.cs_val, yl),
                            xytext=(-5 if right_half else 5, 0),
                            textcoords='offset points', fontsize=8,
                            color=colour, va='center',
                            ha='right' if right_half else 'left')

            if row == 0:
                ax.set_title(DISEASE_LABEL[disease], color=INK, fontsize=10.5,
                             loc='center', pad=8)
            if col == 0:
                ax.set_ylabel(cfg, color=INK, fontsize=9.5)
            ax.set_yticks([])
            ax.xaxis.grid(True); ax.set_axisbelow(True)
            for s in ('top', 'right', 'left'):
                ax.spines[s].set_visible(False)
    for col in range(len(diseases)):
        axes[-1][col].set_xlabel('drug-level connectivity score\n(more negative = stronger reversal)',
                                 color=INK_2)
    handles = [plt.Line2D([], [], color=ORANGE, lw=2),
               plt.Line2D([], [], color=MUTED, lw=2)]
    fig.legend(handles, ['Tier A - approved & disease-modifying',
                         'Tier C - failed phase 3 / withdrawn'],
               loc='upper right', bbox_to_anchor=(0.99, 0.975), fontsize=8.5,
               labelcolor=INK_2)
    fig.suptitle('Where the known drugs sit in each score distribution\n'
                 'all ~3,221 drugs, per disease - qualitative positioning only',
                 x=0.02, ha='left', color=INK, fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    _save(fig, out_path)


def make_plots(base, table, collapse, plot_dir):
    """Two descriptive figures: the percentile heatmap, and the per-disease
    score distributions with the known drugs marked. No inferential panels."""
    _style()
    fig, ax = plt.subplots(figsize=(9.0, 5.4))
    plot_percentile_heatmap(ax, table)
    fig.suptitle('Where do the known ALS/IPF drugs rank, per signature representation?\n'
                 'Descriptive only - too few positives (n=3) for any method claim.',
                 x=0.02, ha='left', color=INK, fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    p1 = Path(plot_dir) / f'known_drug_benchmark_heatmap_{collapse}.png'
    _save(fig, p1)

    p2 = Path(plot_dir) / f'known_drug_benchmark_distributions_{collapse}.png'
    plot_distributions(base, table, collapse, p2)
    return [p1, p2]


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='Known-drug positioning: where the known ALS/IPF drugs '
                    'rank among all drugs, per disease-signature representation.')
    ap.add_argument('--base', default=None,
                    help='NetCoS repo root. Omit to take it from conf.BASE_DIR '
                         '(src/local.py), then from this file\'s location.')
    # Defaults are resolved from conf.py AFTER the repo root is known: result
    # tables -> other_outputs/validations/<slug>, figures -> imgs/validations/<slug>.
    # Never beside this file, which lives under src/.
    ap.add_argument('--out', default=None,
                    help='override directory for result TSVs '
                         '(default: other_outputs/validations/<slug> via conf.py)')
    ap.add_argument('--plots', default=None,
                    help='override directory for PNG figures '
                         '(default: imgs/validations/<slug> via conf.py)')
    # parse_known_args, not parse_args: an IDE may inject argv entries of its
    # own, and an unrecognised one should not abort the run.
    args, unknown = ap.parse_known_args(argv)
    if unknown:
        print(f'[args] ignoring unrecognised argument(s): {unknown}')

    base = resolve_base(args.base)
    default_out, default_plots = output_dirs(base)
    out_dir = Path(args.out) if args.out else default_out
    plot_dir = Path(args.plots) if args.plots else default_plots
    print(f'[out] tables  -> {out_dir}')
    print(f'[out] figures -> {plot_dir}')

    print('=' * 78); print('SIGN-CONVENTION CHECK (corr(CS, pearson) must be > 0)'); print('=' * 78)
    print(check_sign_convention(base).to_string(index=False))

    tables = {}
    for collapse in ('best', 'mean'):
        print('\n' + '=' * 78)
        print(f'KNOWN-DRUG PERCENTILES  --  timepoint collapse = {collapse}')
        print('lower percentile = ranked as a stronger reverser. DESCRIPTIVE ONLY.')
        print('=' * 78)
        t = build_table(base, collapse); tables[collapse] = t
        print(t[['disease', 'drug', 'tier'] + [f'{c}_pct' for c in CONFIGS]].to_string(index=False))

    if out_dir:
        o = Path(out_dir); o.mkdir(parents=True, exist_ok=True)
        for collapse, t in tables.items():
            t.to_csv(o / f'known_drug_benchmark_percentiles_{collapse}.tsv', sep='\t', index=False)
        print(f'\nWrote percentile TSVs to {o}')

    if plot_dir:
        p = Path(plot_dir); p.mkdir(parents=True, exist_ok=True)
        # 'mean' is the headline collapse; 'best' is regenerated for the record.
        for collapse in ('mean', 'best'):
            paths = make_plots(base, tables[collapse], collapse, p)
            for q in paths:
                print(f'Wrote figure {q}')

    # Returned so an interactive session can poke at the results:
    #     res = main()
    #     res['tables']['mean']
    return {'base': base, 'tables': tables,
            'sign_check': check_sign_convention(base)}

#%%
if __name__ == '__main__':
    result = main()
