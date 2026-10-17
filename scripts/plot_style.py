# Vendored, unchanged below this block, from the figure style of the author's CTEPH
# analyses (plot_style.py, copied 2026-10-17) so that this repository
# stands alone; draw_figures.py imports it as a sibling module.
#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Single source of truth for the assoc_plink2 figure aesthetic.
#             - one canonical publication style (Arial, 600 dpi, despined)
#             - a SEMANTIC palette: colour is assigned by what a mark MEANS, so
#               the same meaning is never two colours across figures
#             - the measured, collision-free layout engine (caption block, left
#               margin, value labels, label separation, legend placement)
#             - the GWAS-specific helpers (LD r^2 bins, thresholds, lambda_GC)
#           Imported as a sibling module (ScriptDir is sys.path[0] when Nextflow
#           runs `python3 ${ScriptDir}/x.py`), so no packaging is needed.
# Component: assoc_plink2 — common-variant association scan
# ---------------------------------------------------------------------------
import sys
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

# ── Semantic palette ─────────────────────────────────────────────────────────
# FOUR LAYERS, and a mark is coloured by which layer it belongs to — never by
# what looked good in one panel. Wong (2011) colourblind-safe hues, taken down
# ~15% in saturation from the previous set: at 7.2 in journal width the old
# fully-saturated blue/vermilion vibrated against the grid and pulled the eye to
# whichever panel happened to use the most of it.
#
#   DATA       the measurement itself — points, bars, the estimate
#   ACCENT     what the reader must not miss — the lead variant, genome-wide
#   REFERENCE  what the measurement is judged against — null lines, expectation
#   NEUTRAL    structure that carries no information — grid, bands, fills
#
# Any figure needing more than one DATA colour takes them from COHORT_RAMP (an
# ordered variable) or SERIES (an unordered one), never by inventing a hex.
DATA      = '#3C6E9F'
DATA_DARK = '#1B3B5F'
ACCENT    = '#B0413E'
ACCENT_LT = '#D98C7F'
REFERENCE = '#8A9199'
NEUTRAL   = '#D4D8DC'
NEUTRAL_D = '#9AA1A8'
INK       = '#222222'          # text that is not a label
INK_SOFT  = '#454A50'          # secondary text: subordinate to INK, still fully
                               # legible at small sizes. NOT the same as
                               # REFERENCE, which colours reference LINES and
                               # is too light to read as body text.

# Ordered 3-level ramp: narrow -> intermediate -> full. Light-to-dark is the
# nesting direction, so "more samples" always reads as "lighter" in every figure.
COHORT_RAMP = ['#1B3B5F', '#3C6E9F', '#8FB4D4']

# Unordered discrete series (models, sources) — hue-distinct and Wong-safe.
SERIES = ['#3C6E9F', '#B0413E', '#3F8A6E', '#B08A3C', '#7A5C9E']

# Alternating chromosome bands for whole-genome scatters.
CHROM_BANDS = ['#31577F', '#89A9C6']

# ── Academic math symbols (single source; identical in figures + docs) ───────
NEGLOG10P = r'$-\log_{10}\mathit{P}$'
LAMBDA_GC = r'$\lambda_{\mathrm{GC}}$'
R2_LD     = r'$\mathit{r}^{2}$'
OR_SYM    = r'$\mathrm{OR}$'
MAF_SYM   = r'$\mathrm{MAF}$'
EAF_SYM   = r'$\mathrm{EAF}$'
GW_ALPHA  = 5e-8                       # genome-wide significance
SUGGEST_ALPHA = 1e-5                   # suggestive threshold

# Panel titles are DEFAULTS for the three LD roles; the configured resource for
# each is named by --panel-label so no particular panel is written into the code.
LD_SOURCE_TITLE = {
    'cohort':    'in-sample cohort LD',
    'tommo':     'population co-occurrence panel',
    '1000g_eas': 'out-of-sample reference panel',
}


def sci(x):
    """Render a threshold the way the parameter is written: 5e-8, 1e-6.

    Captions quote the thresholds they were drawn with, so they have to be
    formatted from the value rather than spelled out in the prose.
    """
    mant, exp = f'{float(x):.0e}'.split('e')
    return f'{mant}e{int(exp)}'


def mathsci(x):
    """The same threshold as mathtext, for drawing beside its line."""
    mant, exp = f'{float(x):.0e}'.split('e')
    return rf'$\mathdefault{{{mant}\times10^{{{int(exp)}}}}}$'


def tier_lines(ax, alpha, suggestive, *, decision=None, label=False):
    """Draw both tier thresholds on a -log10 P axis.

    A locus figure shows one peak out of a genome-wide scan, and the reader
    arrives from that scan. Drawing only the genome-wide line leaves a
    suggestive peak with no visible reference at all — the peak floats below
    the only line in the panel and looks like noise. Both lines keep the
    scan figure's vocabulary: genome-wide solid accent, suggestive dashed
    grey.

    decision : 'genome_wide' | 'suggestive' | None. Names the line a stepwise
               procedure actually compared against, drawn heavier so the
               figure states which rule produced its result rather than
               leaving the reader to infer it from the peak's tier.
    label    : write each threshold's value at the right edge. Off by
               default because stacked panels sharing an axis need it once.
    """
    tiers = (('genome_wide', alpha, ACCENT, '-', 1.1),
             ('suggestive', suggestive, REFERENCE, '--', 0.9))
    for name, thr, colour, style, lw in tiers:
        if thr is None:
            continue
        chosen = (name == decision)
        y = -np.log10(thr)
        ax.axhline(y, color=colour, lw=lw + (0.6 if chosen else 0), ls=style, zorder=2)
        if label:
            txt = mathsci(thr) + (' decision' if chosen else '')
            ax.text(0.995, y, txt, transform=ax.get_yaxis_transform(),
                    ha='right', va='bottom', color=colour,
                    fontsize=plt.rcParams['legend.fontsize'] - 1.5,
                    fontweight='bold' if chosen else 'normal', zorder=6,
                    clip_on=False)
POS_MB    = 'position (Mb)'

# ── Peak tiers — defined ONCE, used by every figure that draws a peak ────────
TIER_STYLE = {
    'genome_wide': dict(color=ACCENT, marker='D', size=26, label='genome-wide peak'),
    'suggestive':  dict(color=DATA,   marker='o', size=20, label='suggestive peak'),
}

# What ONE COHORT made of a variant, for the cross-cohort forest. Shape and fill
# carry the state; colour is free for the cohort, so the two encodings never
# compete. Three states, not two: collapsing `suggestive` and `not_a_peak` into a
# single "open" marker throws away the distinction the reader most wants —
# whether a cohort saw anything there at all.
CALLED_STYLE = {
    'genome_wide':     dict(marker='D', size=36, filled=True,
                            label='genome-wide in this cohort'),
    'suggestive':      dict(marker='o', size=30, filled=True,
                            label='suggestive in this cohort'),
    'not_a_peak':      dict(marker='o', size=22, filled=False,
                            label='not significant in this cohort'),
    'not_in_call_set': dict(marker='x', size=26, filled=False,
                            label='not in this cohort\'s call set'),
}


def called_legend_handles(states=None):
    """Legend keys for the three (or four) cross-cohort states, in canonical order."""
    keys = [k for k in CALLED_STYLE if states is None or k in states]
    out = []
    for k in keys:
        st = CALLED_STYLE[k]
        out.append(Line2D([], [], ls='none', marker=st['marker'], markersize=6.5,
                          markerfacecolor=NEUTRAL_D if st['filled'] else 'white',
                          markeredgecolor=NEUTRAL_D, markeredgewidth=1.3, label=st['label']))
    return out

# ── LD r^2 bins — ONE scale for every LD source in every figure ──────────────
# Binning is what makes a truncated LD panel honest: co-occurrence tables
# typically publish only pairs with
# r^2 >= 0.2, so a variant that is IN the panel but returns no pair with the lead
# is not missing — it is bounded below 0.2 and therefore falls entirely inside the
# lowest bin. Only variants absent from the panel altogether are truly unknown,
# and those (and only those) get LD_UNKNOWN_COLOR.
LD_BIN_EDGES  = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
LD_BIN_LABELS = ['< 0.2', '0.2 – 0.4', '0.4 – 0.6', '0.6 – 0.8', '0.8 – 1.0']
LD_BIN_COLORS = ['#4A5FA8', '#6FA0C4', '#5FA87A', '#D9A03C', '#B0413E']
LD_UNKNOWN_COLOR = '#BFC4C9'
LEAD_COLOR = '#7D4C9E'


def ld_bin_colors(r2, known=None):
    """Map r^2 to the shared bin colours; `known=False` entries render grey.

    `known` distinguishes "measured below 0.2" from "no information at all". Pass
    it for a truncated panel (membership) and leave it None elsewhere, where a missing
    r^2 really is missing.
    """
    r2 = pd.to_numeric(pd.Series(r2), errors='coerce')
    idx = np.digitize(r2.fillna(-1.0).to_numpy(), LD_BIN_EDGES[1:-1], right=False)
    out = np.array([LD_BIN_COLORS[min(i, len(LD_BIN_COLORS) - 1)] for i in idx], dtype=object)
    unknown = r2.isna().to_numpy() if known is None else ~np.asarray(known, dtype=bool)
    out[unknown] = LD_UNKNOWN_COLOR
    return out


def ld_legend_handles(include_unknown=True, lead=True):
    """Legend keys for the shared r^2 scale."""
    h = [Line2D([], [], marker='o', linestyle='', markersize=7, markeredgecolor='white',
                markeredgewidth=0.5, color=c, label=l)
         for c, l in zip(LD_BIN_COLORS, LD_BIN_LABELS)]
    if include_unknown:
        h.append(Line2D([], [], marker='o', linestyle='', markersize=7, markeredgecolor='white',
                        markeredgewidth=0.5, color=LD_UNKNOWN_COLOR, label='not in panel'))
    if lead:
        h.append(Line2D([], [], marker='D', linestyle='', markersize=8, color=LEAD_COLOR,
                        markeredgecolor='black', markeredgewidth=0.6, label='lead'))
    return h


def tier_legend_handles(tiers=None):
    """Legend keys for the peak tiers, in the canonical order."""
    keys = tiers or list(TIER_STYLE)
    return [Line2D([], [], ls='none', marker=TIER_STYLE[k]['marker'],
                   color=TIER_STYLE[k]['color'], markersize=7,
                   markeredgecolor='white', markeredgewidth=0.6,
                   label=TIER_STYLE[k]['label']) for k in keys if k in TIER_STYLE]


# ── What each symbol MEANS IN THIS STUDY (not a generic glossary) ────────────
# key -> (display label, definition). Rendered one-per-line into figure captions
# by caption_block(defs=[...]) and into every sidecar by figure_doc. Every entry
# states the quantity, its units where applicable, and what it does and does not
# license in THIS design. Keep in lock-step with docs/METHODS.md §2.
SYMBOL_DEFS = {
    'model': ('genetic model',
              r'the encoding of the genotype in the GLM: ADD counts alt alleles (0/1/2), DOM contrasts '
              r'carriers against non-carriers, REC contrasts alt-homozygotes against the rest. Three '
              r'separate genome-wide scans, not one joint test; ADD is the primary.'),
    'or': (r'$\mathrm{OR}$',
           r'odds ratio per copy of the A1 (effect) allele from logistic regression, with its 95% CI. '
           r'Compared on the log scale, so a protective and a risk allele of equal strength are equally '
           r'far from 1.'),
    'lambda_gc': (r'$\lambda_{\mathrm{GC}}$',
                  r'genomic-control inflation factor = median $\chi^2$ / 0.4549. $\lambda>1$ is inflation, '
                  r'$\lambda<1$ deflation; neither is corrected for here — it is reported as a calibration '
                  r'read-out, and its interpretation is in METHODS §7.'),
    'gw_sig': ('genome-wide significance',
               r'$P<5\times10^{-8}$, applied identically to every cohort and model. The three cohorts are '
               r'nested and the three models correlated, so these are not independent tests and no further '
               r'multiplicity adjustment is made — stated, not silently assumed.'),
    'ld_r2': (r'$r^{2}$',
              r'linkage disequilibrium with the lead variant, binned on one scale for all sources. A '
              r'co-occurrence panel typically '
              r'publishes only pairs with $r^{2}\geq0.2$, so a variant present in that panel but absent from '
              r'the query is bounded below 0.2 and sits in the lowest bin; only variants absent from the '
              r'panel entirely are unknown (grey).'),
    'neff': (r'$N_{\mathrm{eff}}$',
             r'$4/(1/N_{\mathrm{case}}+1/N_{\mathrm{ctrl}})$, the balanced-design equivalent sample size. '
             r'Reported as a cohort descriptor — it is the size any later meta-analysis or mixed-model '
             r'run quotes — and it collapses toward 4x the smaller arm as the design becomes unbalanced.'),
    'eaf': (r'$\mathrm{EAF}$',
            r'frequency of the EFFECT allele (plink2 A1) among called samples of that group. Case and '
            r'control EAF are the two numbers the odds ratio is computed from, so plotting them directly '
            r'shows what drives an effect estimate and whether one group carries the whole difference.'),
    'pip': ('PIP',
            r'SuSiE posterior inclusion probability — the probability a variant is causal given the locus '
            r'summary statistics and an in-sample LD matrix. Credible sets are the smallest variant groups '
            r'covering 95% posterior mass.'),
    'conditional': ('conditional analysis',
                    r'the locus re-fitted with the lead variant (then each further signal) added as a '
                    r'covariate. A signal that survives conditioning is independent of the lead; one that '
                    r'vanishes was tagging it.'),
    'errcode': ('ERRCODE',
                r"plink2's per-variant fit diagnostic. Variants with a non-'.' code (e.g. VIF_INFINITE, "
                r'SEPARATION) did not fit cleanly and are excluded from $\lambda_{\mathrm{GC}}$ and from the '
                r'hit list rather than silently carried.'),
    'called_peak': ('called_peak',
                    r'whether the variant was a peak *in that cohort* (genome_wide / suggestive) or not '
                    r'(not_a_peak). An estimate is reported either way — a cohort that called no peak still '
                    r'has an odds ratio there, and a blank cell would be indistinguishable from a missing one.'),
}

# ── Shared formula bank (academic notation; single source for figures + docs) ─
def glm_formula(pc_label=None, n_pcs=None):
    """The fitted model, in mathtext, with the PC space as a RUN PARAMETER.

    The number of PCs and the space they were computed in are configuration, not
    properties of the method, so they arrive from the caller (ultimately from the
    pipeline's PcLabel / NPcs) rather than being written into this module. With
    neither given the formula degrades to a generic K and names no space, which is
    what a reader of an unconfigured figure should see rather than another study's.

    `pc_label` is the SPACE the PCs were computed in, not the whole covariate set:
    the formula already writes SEX and the PC sum, so passing the full CovarLabel
    here renders the redundant 'PCs = SEX + N <space> PCs'.
    """
    k = str(int(n_pcs)) if n_pcs else 'K'
    tail = f'; PCs = {pc_label}' if pc_label else ''
    return (r'$\mathrm{logit}\,\Pr(\mathrm{case}_i) = \beta_0 + \beta\,g_i'
            r' + \gamma_{\mathrm{sex}}\,\mathrm{SEX}_i'
            r' + \sum_{k=1}^{' + k + r'}\gamma_k\,\mathrm{PC}_{k,i}$      '
            r'($g_i$ = genotype under the stated model' + tail + ')')


FORMULAS = {
    'glm': glm_formula(),
    'lambda': r'$\lambda_{\mathrm{GC}} = \mathrm{median}(\chi^{2})\,/\,\chi^{2}_{1,0.5}$,   '
              r'$\chi^{2}_{1,0.5} = 0.4549$',
    'neff': r'$N_{\mathrm{eff}} = 4\,/\,(1/N_{\mathrm{case}} + 1/N_{\mathrm{ctrl}})$',
    'susie': r'$\mathrm{susie\_rss}(\hat{\beta}, \mathrm{se}(\hat{\beta}), R, n)$,   '
             r'$R$ = in-sample LD (signed $r$), $L = 10$, credible-set coverage 0.95',
}


# ── Formatters ────────────────────────────────────────────────────────────────
def _sci_fmt(v, _):
    """Scientific axis labels in superscript form (e.g. 1.6×10⁷), publication-style.

    Small magnitudes (|v| < 1000) are shown as plain integers; 0 stays '0'.
    """
    import math
    if v == 0:
        return '0'
    if abs(v) < 1000:
        return f'{v:g}'
    e = int(math.floor(math.log10(abs(v))))
    m = v / 10 ** e
    ms = f'{m:.1f}'.rstrip('0').rstrip('.')      # 1-decimal mantissa, drop trailing .0
    return rf'${ms}\times10^{{{e}}}$'


COUNT_FMT = FuncFormatter(_sci_fmt)
# Plain thousands separators, for count axes small enough to read as integers.
INT_FMT = FuncFormatter(lambda v, _: format(int(v), ','))


def p_tex(p):
    """A P-value as publication mathtext: 1.2×10⁻⁸, or '< 1×10⁻³⁰⁰' on underflow."""
    try:
        p = float(p)
    except (TypeError, ValueError):
        return '—'
    if not np.isfinite(p) or p <= 0:
        return r'$<1\times10^{-300}$'
    if p >= 0.01:
        return f'{p:.3g}'
    e = int(np.floor(np.log10(p)))
    m = p / 10 ** e
    return rf'${m:.2f}\times10^{{{e}}}$'


# ── Type scales: one definition, two sizes ───────────────────────────────────
#   slide : larger type so a figure stays legible when projected (PPT / talk).
#   paper : tighter type for a dense multi-panel manuscript figure.
_SCALES = {
    'slide': dict(base=13, title=15, label=14, tick=12, legend=11, panel=16),
    # Journal double-column: 183 mm / 7.2 in. Sizes follow the usual 7-9 pt range
    # for printed figures.
    'paper': dict(base=8, title=9, label=8, tick=7, legend=7, panel=9),
}

# The scale in force, so panel_tag()/panel_letters() size themselves from the
# figure they are actually drawing on rather than from a hard-coded 'slide'.
_ACTIVE = _SCALES['slide']

# Journal column widths, in inches — the width every figure is authored against.
COL_SINGLE, COL_ONEHALF, COL_DOUBLE = 3.50, 4.72, 7.20


def setup_style(scale='slide'):
    """Apply the canonical rcParams at the requested type scale; return the scale dict."""
    global _ACTIVE
    s = _SCALES[scale]
    _ACTIVE = s
    plt.rcParams.update({
        'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
        'savefig.dpi': 600, 'savefig.bbox': None, 'figure.dpi': 150,
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'mathtext.fontset': 'dejavusans',
        'font.size': s['base'], 'axes.titlesize': s['title'], 'axes.labelsize': s['label'],
        'xtick.labelsize': s['tick'], 'ytick.labelsize': s['tick'], 'legend.fontsize': s['legend'],
        'axes.linewidth': 0.9, 'axes.edgecolor': INK,
        'text.color': INK, 'axes.labelcolor': INK,
        'xtick.color': INK, 'ytick.color': INK,
        'xtick.direction': 'out', 'ytick.direction': 'out',
        'axes.grid': True, 'grid.color': NEUTRAL, 'grid.linewidth': 0.7, 'grid.alpha': 1.0,
        'axes.spines.top': False, 'axes.spines.right': False,
    })
    return s


def shorten(names):
    """{name: short name} by dropping the suffix EVERY name shares.

    'a_suffix'/'b_suffix' -> 'a'/'b'. Computed from the names given rather
    than from a hard-coded suffix, so any project's cohort naming collapses the
    same way and none has to be known here. Falls back to the full names when
    there is no shared suffix, when trimming would empty a name, or when two
    names would collide.
    """
    names = list(names)
    if len(names) < 2:
        return {n: n for n in names}
    rev = [n[::-1] for n in names]
    k = 0
    while k < min(map(len, rev)) and len({r[k] for r in rev}) == 1:
        k += 1
    # Cut back to a separator so a shared letter run is not sliced mid-word.
    while k > 0 and names[0][-k] not in '_-.':
        k -= 1
    out = {n: (n[:-k] or n) for n in names} if k else {n: n for n in names}
    return out if len(set(out.values())) == len(names) else {n: n for n in names}


# ── Axis / annotation helpers ─────────────────────────────────────────────────
def despine(ax, grid_axis='both'):
    """Recessive grid, top/right spines off, outward ticks — the shared axis look."""
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if grid_axis in ('both', 'x', 'y'):
        ax.grid(True, which='major', axis=grid_axis)
    else:
        ax.grid(False)
    ax.tick_params(direction='out', length=4.0, width=0.9)


def panel_letters(axes, letters='abcdefgh', size=None, dx=-0.07, dy=1.04, ha='right', inside=False):
    """Bold (a)/(b)/… panel tags in one consistent position/size across all figures.

    inside=True places the tag in the panel's empty upper-left corner (use when an
    outside tag would collide with a wide/2-line y-axis label or the suptitle).
    """
    size = size or _ACTIVE['panel']
    for ax, letter in zip(axes, letters):
        if inside:
            ax.text(0.02, 0.97, f'({letter})', transform=ax.transAxes, fontsize=size,
                    fontweight='bold', ha='left', va='top', zorder=25)
        else:
            ax.text(dx, dy, f'({letter})', transform=ax.transAxes,
                    fontsize=size, fontweight='bold', ha=ha, va='bottom')


def panel_tag(ax, letter, title=None, *, size=None, tsize=None, pad=7):
    """Nature/Science panel label: a **bold letter** flush-left in the title slot,
    with an optional short centred title. Both live ABOVE the axes box (matplotlib
    reserves the space), so neither can ever overlap the data, ticks, or y-label.

    Keep the centred title short so its left edge stays clear of the letter on
    narrow panels; below ~2.2 in, pass no title at all and describe the panel in
    the caption, as journals do.
    """
    size = size or (_ACTIVE['panel'] - 1)
    ax.set_title(f'{letter}', loc='left', fontweight='bold', fontsize=size, pad=pad)
    if title:
        ax.set_title(title, loc='center', fontsize=tsize or plt.rcParams['axes.titlesize'], pad=pad)


def thin_tick_labels(ax, axis='x', pad_px=2.0):
    """Blank tick labels that would overlap their kept neighbour.

    Chromosome axes are the motivating case: chr19-22 are short enough that at
    any readable font their labels collide, and matplotlib will happily draw
    them on top of one another. Measure the rendered boxes and drop whichever
    labels cannot fit, keeping the leftmost of each colliding run.
    """
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    labels = ax.get_xticklabels() if axis == 'x' else ax.get_yticklabels()
    texts, kept_end = [], None
    for lab in labels:
        txt = lab.get_text()
        if not txt:
            texts.append('')
            continue
        bb = lab.get_window_extent(renderer=r)
        lo, hi = (bb.x0, bb.x1) if axis == 'x' else (bb.y0, bb.y1)
        if kept_end is not None and lo < kept_end + pad_px:
            texts.append('')
        else:
            texts.append(txt)
            kept_end = hi
    # Reinstall through the formatter, not by set_text on the artists: with a
    # FixedFormatter in place every later draw regenerates the labels from the
    # formatter and would undo an artist-level edit.
    setter = ax.set_xticklabels if axis == 'x' else ax.set_yticklabels
    size = labels[0].get_fontsize() if labels else None
    setter(texts, fontsize=size)
    return ax


def or_log_axis(ax, axis='x'):
    """Label a log OR axis with odds ratios, not powers of ten.

    A log axis over a range like 1.2-3.5 gets one decade tick ($10^0$) from the
    default locator and nothing else, and ScalarFormatter alone does not add the
    missing ticks back. Place readable OR values explicitly.
    """
    from matplotlib.ticker import FixedLocator, NullFormatter, ScalarFormatter
    a = ax.xaxis if axis == 'x' else ax.yaxis
    lo, hi = ax.get_xlim() if axis == 'x' else ax.get_ylim()
    cand = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.35, 1.5, 1.75,
            2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0]
    ticks = [c for c in cand if lo <= c <= hi]
    # A wide range (say 0.5-5) keeps a dozen candidates and the labels collide.
    while len(ticks) > 7:
        ticks = ticks[::2]
    if len(ticks) < 3:
        ticks = list(np.geomspace(lo, hi, 5))
    a.set_major_locator(FixedLocator(ticks))
    a.set_major_formatter(ScalarFormatter())
    a.set_minor_formatter(NullFormatter())
    a.set_minor_locator(FixedLocator([]))


# ── Collision control: four helpers, one per defect class ────────────────────
# Everything here MEASURES with the renderer instead of guessing offsets. The
# hand-rolled alternatives these replace (a fixed left=0.088, a `k % 2` label
# stagger, a corner legend chosen by eye) each failed silently the moment the
# data changed shape — a clipped tick label or a legend over the data reads as
# a rendering bug, not as a layout parameter that needs tuning.

def _renderer(fig):
    fig.canvas.draw()
    return fig.canvas.get_renderer()


def fit_left_margin(fig, axes=None, *, pad_in=0.07, min_frac=0.05, max_frac=0.34):
    """Left margin (figure fraction) that actually fits the y tick labels + y label.

    A hardcoded fraction is a guess about the widest label, and the guess breaks
    the moment a label gets longer: `cohort_compare` panel (c) rendered
    'intermediate' clipped off the canvas edge under left=0.088. Measure instead.

    Only left-labelled axes in the leftmost column are considered, so a twinx()
    (whose ticks sit on the right) cannot inflate the left margin.
    """
    try:
        r = _renderer(fig)
    except Exception:
        return min_frac
    axes = list(axes) if axes is not None else list(fig.axes)
    axes = [ax for ax in axes if ax.yaxis.get_label_position() == 'left']
    if not axes:
        return min_frac
    x0min = min(ax.get_position().x0 for ax in axes)
    tick_in = lab_in = 0.0
    for ax in axes:
        if ax.get_position().x0 > x0min + 0.01:
            continue
        for t in ax.get_yticklabels():
            if t.get_text():
                tick_in = max(tick_in, t.get_window_extent(renderer=r).width / fig.dpi)
        yl = ax.yaxis.get_label()
        if yl.get_text():
            lab_in = max(lab_in, yl.get_window_extent(renderer=r).width / fig.dpi)
        # Annotations parked left of the axes — grouped_row_labels' group names.
        # They are not tick labels, so without this they measure as nothing and
        # get clipped off the canvas edge.
        ax_x0 = ax.get_window_extent(renderer=r).x0
        for t in ax.texts:
            if not t.get_text():
                continue
            try:
                bb = t.get_window_extent(renderer=r)
            except Exception:
                continue
            if bb.x1 <= ax_x0 + 1:
                tick_in = max(tick_in, (ax_x0 - bb.x0) / fig.dpi)
    # 3x the pad, not 2.2: at 2.2 the y-axis label of the tallest-label panel
    # ended up touching the canvas edge (measured at -0.9 px).
    need = tick_in + lab_in + 3.0 * pad_in
    return float(min(max_frac, max(min_frac, need / fig.get_size_inches()[0])))


def value_labels(ax, xs, ys, texts, *, axis='y', offset=9.0, dx=0.0, expand=True, **kw):
    """Annotate points with their own values, offset off the marker, WITH headroom.

    Two defects this closes, both present in the previous cohort_compare: a value
    label drawn past the top of the axes (lambda = 1.132), and a label landing on
    the stem it annotates (the N_eff series). The axis limit is grown after
    measuring, so the label is inside the panel whatever the data range is.

    Returns the Annotation artists, so they can be handed to spread_labels().
    """
    style = dict(ha='center', va='bottom' if axis == 'y' else 'center',
                 fontsize=plt.rcParams['legend.fontsize'] - 0.5, zorder=12)
    style.update(kw)
    # `dx` nudges the label off a vertical element it would otherwise sit on — a
    # stem, an error bar, the marker's own connecting line.
    dx, dy = (dx, offset) if axis == 'y' else (offset, dx)
    out = [ax.annotate(str(t), xy=(x, y), xytext=(dx, dy), textcoords='offset points', **style)
           for x, y, t in zip(xs, ys, texts)]
    if expand and out:
        _grow_axis(ax, out, axis)
    return out


def _grow_axis(ax, artists, axis='y'):
    """Extend the axis limit so every artist's rendered box fits inside the axes."""
    fig = ax.figure
    try:
        r = _renderer(fig)
    except Exception:
        return
    inv = ax.transData.inverted()
    lo, hi = ax.get_ylim() if axis == 'y' else ax.get_xlim()
    box = ax.get_window_extent()
    need = None
    for a in artists:
        try:
            bb = a.get_window_extent(renderer=r)
        except Exception:
            continue
        px = bb.y1 if axis == 'y' else bb.x1
        if px > (box.y1 if axis == 'y' else box.x1):
            pt = inv.transform((bb.x1, bb.y1))
            v = pt[1] if axis == 'y' else pt[0]
            need = v if need is None else max(need, v)
    if need is not None and np.isfinite(need):
        span = hi - lo
        if axis == 'y':
            ax.set_ylim(lo, need + 0.04 * span)
        else:
            ax.set_xlim(lo, need + 0.04 * span)


def repel_from_centre(pos, min_gap, lo, hi, max_iter=300):
    """Separate 1-D positions by pushing overlapping pairs APART FROM THEIR MIDPOINT.

    gwaslab's `adjust_text_position` rule: find the first pair closer than
    `min_gap`, move both outward by half the shortfall, repeat. Symmetric, so the
    displacement is shared between neighbours instead of accumulating.

    A greedy left-to-right sweep (what spread_labels does) only ever pushes right,
    so a dense cluster drags every later label with it — on a dense scan the
    chr2-6 run ended up shifted past chr6, pointing at peaks it did not belong to.

    Finally clamped to [lo, hi]: shifted back if it overhangs one edge, and
    uniformly compressed if the whole run no longer fits, so a label can never be
    pushed off the axis.
    """
    p = np.asarray(pos, dtype=float)
    if not len(p):
        return p
    order = np.argsort(p)
    q = p[order].astype(float).copy()
    for _ in range(max_iter):
        gaps = np.diff(q)
        bad = np.where(gaps < min_gap)[0]
        if not len(bad):
            break
        i = bad[0]
        need = (min_gap - gaps[i]) / 2.0
        q[:i + 1] -= need
        q[i + 1:] += need
    # Fitting the result into [lo, hi]. Symmetric repelling pushes the OUTERMOST
    # labels past both edges — a dense middle cluster propagates outward — and
    # scaling the whole run back to fit shrinks every gap below min_gap again,
    # which is how ten dominant-scan labels ended up 10 % too close and two of
    # them printed on top of each other. So enforce the gap directly instead.
    need = min_gap * (len(q) - 1)
    if need <= (hi - lo):
        q[0] = max(q[0], lo)
        for i in range(1, len(q)):                 # forward: honour min_gap
            q[i] = max(q[i], q[i - 1] + min_gap)
        if q[-1] > hi:                             # overran: pull back from the right
            q[-1] = hi
            for i in range(len(q) - 2, -1, -1):
                q[i] = min(q[i], q[i + 1] - min_gap)
    else:
        # Genuinely more labels than the axis can hold at a readable size; spread
        # them evenly, which is the least-bad option, and let the caller cap the
        # count instead.
        q = np.linspace(lo, hi, len(q))
    out = np.empty_like(p, dtype=float)
    out[order] = q
    return out


def _text_size_in(fig, s, fontsize, weight='normal', style='italic'):
    """(width, height) of a string in inches, without leaving it on the figure."""
    try:
        r = _renderer(fig)
        t = fig.text(0, 0, s, fontsize=fontsize, fontweight=weight, fontstyle=style)
        bb = t.get_window_extent(renderer=r)
        t.remove()
        return bb.width / fig.dpi, bb.height / fig.dpi
    except Exception:
        return 0.6 * len(s) * fontsize / 72.0, 1.2 * fontsize / 72.0


# gwaslab's annotation styles, by name. Read from the user's install
# (gwaslab/viz/viz_aux_annotate_plot.py, and the plot_mqq signature): `right` is
# gwaslab's own default and sweeps greedily left to right; `expand` runs the
# symmetric repel and forces the text vertical; `tight` is `expand` with a short
# fixed arm. THIS COMPONENT DEFAULTS TO `expand`.
ANNO_STYLES = {
    'right':  dict(rotation=40, sweep='greedy'),
    'expand': dict(rotation=90, sweep='repel'),
    'tight':  dict(rotation=90, sweep='repel'),
    # `auto` is this component's default: measure, then use the least intrusive
    # rotation the labels actually fit in. Always-vertical is right for eleven
    # names on a peak panel and needlessly ugly for the two on a Manhattan.
    'auto':   dict(rotation=None, sweep='repel'),
}
REPEL_FORCE = 0.03          # gwaslab's default; a fraction of the plotted x-span


def _greedy_forward(pos, min_gap, lo, hi):
    """gwaslab's `right` sweep: keep a position if it clears the last kept one by
    min_gap, otherwise push it to exactly that."""
    p = np.asarray(pos, dtype=float)
    order = np.argsort(p)
    q = p[order].astype(float).copy()
    last = -np.inf
    for i in range(len(q)):
        q[i] = q[i] if q[i] > last + min_gap else last + min_gap
        last = q[i]
    if q[-1] > hi:                      # ran off the right edge: pull the run back
        q -= q[-1] - hi
        q[0] = max(q[0], lo)
        for i in range(1, len(q)):
            q[i] = max(q[i], q[i - 1] + min_gap)
    out = np.empty_like(p, dtype=float)
    out[order] = q
    return out


def gene_labels(ax, xs, ys, texts, *, anno_style='auto', repel_force=REPEL_FORCE,
                color=INK, arm_color=None, fontsize=None, weight='normal',
                rotation=None, gap_in=0.055):
    """Locus labels in the gwaslab idiom, in a strip ABOVE the data area.

    Matches `gwaslab/viz/viz_aux_annotate_plot.py`:
      - labels sit in one band, not offset from each point;
      - text is italic and rotated — `anno_style` picks the angle;
      - the leader is `arc,angleA=0,armA=0,angleB=90,armB=<px>` — VERTICAL up from
        the peak, then horizontal to the text. The vertical segment is what marks
        the peak's true position, so the horizontal segment can absorb as much
        displacement as the repel step needs without the reader losing the anchor;
      - `expand` repels symmetrically, `right` sweeps greedily forward, `auto`
        (the default) measures and takes the least intrusive rotation that fits.

    THE LABELS DO NOT LIVE IN THE DATA AREA. The axes box is shrunk by the strip
    height and the text drawn above it, so the data limits are never inflated to
    make room for text. They used to be: the peak panel's ceiling was pushed to
    ~2x the largest -log10 P and the peaks ended up in the bottom third of the
    panel, so two cohorts' panels differed by however much whitespace their label
    counts happened to need.

    THE MINIMUM GAP IS `max(measured footprint, repel_force * span)`. gwaslab
    spaces on `repel_force * span` alone, which overlaps as soon as the labels are
    wider than that fraction of the axis; the measured footprint is the floor that
    makes a collision impossible. Keeping repel_force on top lets a caller spread
    labels FURTHER than they strictly need, which is what the parameter is for.

    MUST be called after the final layout (caption_block resizes the figure) — the
    band height and the arm length are computed from the rendered axes box.
    """
    fig = ax.figure
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    if not len(xs):
        return []
    style = ANNO_STYLES.get(anno_style, ANNO_STYLES['auto'])
    fs = fontsize if fontsize is not None else plt.rcParams['legend.fontsize']
    # colour, weight and size may be given PER LABEL, so one band can carry two
    # tiers — the distinction then lives in the ink rather than in a second strip
    # at a second rotation.
    n = len(xs)
    cols = list(color) if isinstance(color, (list, tuple, np.ndarray)) else [color] * n
    wts = list(weight) if isinstance(weight, (list, tuple, np.ndarray)) else [weight] * n
    fss = list(fontsize) if isinstance(fontsize, (list, tuple, np.ndarray)) else [fs] * n
    arm_color = arm_color or NEUTRAL_D
    try:
        _renderer(fig)
    except Exception:
        pass
    box = ax.get_window_extent()
    ax_w_in = box.width / fig.dpi

    # Horizontal footprint at a rotation: the width of the AXIS-ALIGNED box around
    # rotated text is w*cos(t) + h*sin(t). Dropping the h*sin(t) term understates a
    # 40-degree label by its own line height, which was enough for two names in a
    # dense band to print on top of each other.
    sizes = [_text_size_in(fig, s, f, w) for s, f, w in zip(texts, fss, wts)]

    def _gap(deg):
        rad = np.radians(deg)
        c, s = abs(np.cos(rad)), abs(np.sin(rad))
        return max(w * c + h * s for w, h in sizes) + gap_in

    rot = rotation if rotation is not None else style['rotation']
    if rot is None:                      # anno_style='auto'
        # The least intrusive rotation the labels actually fit in, measured
        # against the axis width — so the choice can never drift from the spacing.
        rot = next((d for d in (0, 40) if _gap(d) * len(texts) <= ax_w_in), 90)
    min_gap_in = _gap(rot)

    x0, x1 = ax.get_xlim()
    per_in = (x1 - x0) / ax_w_in if ax_w_in else 1.0
    min_gap = max(min_gap_in * per_in, repel_force * (x1 - x0))
    lo, hi = x0 + 0.5 * min_gap_in * per_in, x1 - 0.5 * min_gap_in * per_in
    lab_x = (repel_from_centre(xs, min_gap, lo, hi) if style['sweep'] == 'repel'
             else _greedy_forward(xs, min_gap, lo, hi))

    # ── the strip: shrink the AXES, do not inflate the DATA ────────────────
    # Labels used to live inside the data limits, so the peak panel's ceiling was
    # pushed to ~2x the largest -log10 P and the data ended up occupying a third
    # of the panel. Reserving the space by shrinking the axes box instead leaves
    # the data range untouched: two labels cost two labels' worth of height, and
    # twenty cost twenty's, on any panel size.
    strip_in = max(h * abs(np.cos(np.radians(rot))) + w * abs(np.sin(np.radians(rot)))
                   for w, h in sizes) + 2.2 * gap_in
    pos = ax.get_position()
    fig_h_in = fig.get_size_inches()[1]
    strip_frac = min(strip_in / fig_h_in, 0.55 * pos.height)
    prev_frac = getattr(ax, '_label_strip_frac', 0.0)
    ax.set_position([pos.x0, pos.y0, pos.width, pos.height - strip_frac])
    # Record what this axes gave up, so a sibling in the same row can be brought
    # down to match and the panel letter can be pushed above the strip.
    ax._label_strip_frac = prev_frac + strip_frac
    # A SECOND band on the same axes stacks above the first. Each call shrinks the
    # axes and draws just above it, so without this offset the genome-wide names
    # would be drawn straight into the strip the suggestive names already occupy.
    new_h = ax.get_position().height
    base_y = 1.015 + (prev_frac / new_h if new_h > 0 else 0.0)

    out = []
    trans = ax.get_xaxis_transform()          # x in data, y in axes fraction
    for x, y, lx, t, col, wt, f in zip(xs, ys, lab_x, texts, cols, wts, fss):
        # The arm rises from the point to the top of the (shrunken) axes, then
        # runs horizontally to the label. armB is that rise, in pixels.
        arm_px = max(2.0, ax.transAxes.transform((0, base_y))[1]
                     - ax.transData.transform((0, y))[1])
        out.append(ax.annotate(
            t, xy=(x, y), xytext=(lx, base_y), textcoords=trans,
            rotation=rot, rotation_mode='anchor', ha='left', va='bottom',
            fontstyle='italic', fontweight=wt, fontsize=f, color=col,
            zorder=9, annotation_clip=False,
            arrowprops=dict(arrowstyle='-', color=arm_color, lw=0.6, shrinkA=0, shrinkB=2,
                            relpos=(0, 0),
                            connectionstyle=f'arc,angleA=0,armA=0,angleB=90,'
                                            f'armB={arm_px:.1f},rad=0')))
        out[-1].set_clip_on(False)
    _report_label_overlaps(ax, out)
    return out


def grouped_row_labels(ax, groups, *, italic=True, band=True, band_alpha=0.13,
                       name_pad_in=0.26, fontsize=None):
    """Row ticks carry the MEMBER; the group name is written ONCE beside its block.

    A cross-cohort forest repeats every locus once per cohort, and labelling each
    row `GENE - cohort` says the gene name three times while leaving the reader to
    work out which three rows belong together. Here the ticks carry only the
    cohort, the gene is written once at its block's vertical centre, and
    alternate blocks get a faint band so the grouping is visible without a rule.

    groups : [(name, [member_label, ...]), ...] in TOP-TO-BOTTOM draw order.
             Rows are assigned y = n-1 down to 0, matching the forests' own
             convention of drawing the strongest locus at the top.

             `name` is either a plain string, or a list of (text, kwargs) lines
             stacked at the block's centre — which is how a locus carries its
             gene symbol, the variant actually plotted and its rsID without
             repeating any of them on all three rows.

    Returns the y positions used, so the caller can place per-row artists on the
    same coordinates.
    """
    n = sum(len(m) for _g, m in groups)
    fs = fontsize or plt.rcParams['ytick.labelsize']
    ticks, labels, y = [], [], n - 1
    blocks = []
    for gi, (name, members) in enumerate(groups):
        top = y
        for m in members:
            ticks.append(y)
            labels.append(m)
            y -= 1
        bottom = y + 1
        blocks.append((name, top, bottom, gi))
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=fs)
    ax.set_ylim(-0.7, n - 0.3)
    offset_pt = _tick_offset_pt(ax, name_pad_in)

    line_h = fs * 1.35 / 72.0            # inches per stacked line
    for name, top, bottom, gi in blocks:
        if band and gi % 2 == 1:
            # Behind everything, and clipped, so it reads as a ground rather than
            # a datum. axhspan spans the full axes width by construction.
            ax.axhspan(bottom - 0.5, top + 0.5, color=NEUTRAL, alpha=band_alpha,
                       linewidth=0, zorder=0)
        # The name sits left of the tick column, in FIGURE-independent offset
        # points so it keeps its distance whatever the tick labels turn out to be.
        lines = ([(name, dict(fontstyle='italic' if italic else 'normal',
                              fontweight='bold', fontsize=fs))]
                 if isinstance(name, str) else list(name))
        # Stacked about the block's centre so the group reads as one label.
        y0 = (top + bottom) / 2
        span = (len(lines) - 1) * line_h * 72.0 / 2.0
        for k, (txt, kw) in enumerate(lines):
            if not txt:
                continue
            ax.annotate(txt, xy=(0, y0), xycoords=('axes fraction', 'data'),
                        xytext=(-offset_pt, span - k * line_h * 72.0),
                        textcoords='offset points', ha='right', va='center',
                        annotation_clip=False,
                        **{'fontsize': fs, 'color': INK, **kw})
    return ticks


def _tick_offset_pt(ax, pad_in):
    """Points to clear the widest y tick label, MEASURED once it is rendered.

    Estimating the width from the string is what failed here: the fallback was
    narrower than 'intermediate' and every group name landed on top of a tick.
    The labels exist by this point, so measure them.
    """
    fig = ax.get_figure()
    w = 0.0
    try:
        r = _renderer(fig)
        for t in ax.get_yticklabels():
            if t.get_text():
                w = max(w, t.get_window_extent(renderer=r).width / fig.dpi)
    except Exception:
        pass
    if w <= 0:
        w = 0.9                      # a wide default: overlap is worse than a gap
    return (w + pad_in) * 72.0


def align_panel_tops(*axes):
    """Bring a row of axes to the same top edge.

    gene_labels() reserves its strip by SHRINKING the axes, which is what keeps the
    data limits honest — but it shrinks only the panel being labelled, so its
    unlabelled neighbour in the same gridspec row is left standing taller and the
    row stops reading as a row.
    """
    axes = [a for a in axes if a is not None]
    if len(axes) < 2:
        return
    top = min(a.get_position().y1 for a in axes)
    for a in axes:
        p = a.get_position()
        if p.y1 > top:
            a.set_position([p.x0, p.y0, p.width, top - p.y0])


def equalise_row_heights(*axes):
    """Shrink a STACK of axes to a common data height, keeping each one's top.

    The companion to align_panel_tops(), which matches panels ACROSS a row. This
    one matches them DOWN a column, and it exists because gene_labels() takes its
    strip out of only the panel it labels: on a stacked comparison figure one row
    carries the label band and the others do not, so their data areas end up
    different heights.

    That is not cosmetic. Stacked panels of the same quantity are drawn on ONE
    shared y-limit precisely so heights are comparable between them — and unequal
    axes heights re-introduce exactly the distortion the shared limit removes, by
    giving the same -log10 P a different number of pixels in each row.

    Each axes keeps its BOTTOM edge and gives up height at the top, so every row
    loses the same amount and the gaps between them stay uniform. Keeping the top
    instead would raise the lower rows' bottoms and open a gap that grows down the
    figure.
    """
    axes = [a for a in axes if a is not None]
    if len(axes) < 2:
        return
    h = min(a.get_position().height for a in axes)
    for a in axes:
        p = a.get_position()
        if p.height > h:
            a.set_position([p.x0, p.y0, p.width, h])


def strip_pad(ax, base=7.0):
    """Title pad, in points, that clears whatever label strip `ax` gave up.

    The panel letter lives in the title slot, which is attached to the axes, so
    after the strip is reserved the letter would sit BELOW the gene names.
    """
    frac = getattr(ax, '_label_strip_frac', 0.0)
    return base + frac * ax.figure.get_size_inches()[1] * 72.0


def _report_label_overlaps(ax, annots):
    """Measure the placed labels and say so on stderr if any two touch.

    Self-verifying, because the spacing depends on rendered text metrics and so
    cannot be checked by reading the code. Measured on the TEXT ONLY:
    Annotation.get_window_extent() unions in the arrow patch, so using it reports
    the arms as collisions and every figure looks broken.

    Reports rather than raises — with more labels than the axis can hold at a
    readable size there is no non-overlapping answer, and the caller should cap
    the count instead of losing the figure.
    """
    if len(annots) < 2:
        return 0
    from matplotlib.text import Text as _Text
    try:
        r = _renderer(ax.figure)
        boxes = sorted((_Text.get_window_extent(a, renderer=r) for a in annots),
                       key=lambda b: b.x0)
    except Exception:
        return 0
    n = sum(1 for b1, b2 in zip(boxes, boxes[1:]) if b2.x0 < b1.x1)
    if n:
        print(f'[gene_labels] WARNING: {n} of {len(annots) - 1} adjacent label pairs overlap',
              file=sys.stderr)
    return n


def spread_labels(ax, annots, *, axis='y', pad_px=2.5):
    """Push apart text annotations whose rendered boxes overlap, along one axis.

    Greedy 1-D separation on measured boxes. `annots` must be Annotations created
    with textcoords='offset points' — the shift is applied to the offset, so the
    anchor (and any leader line supplied via arrowprops) still points at the data.

    This replaces the per-figure `k % 2` alternating-row hacks, which staggered
    labels on a fixed rhythm whether or not they actually collided and therefore
    both wasted space and still overlapped when three labels landed together.
    """
    fig = ax.figure
    if not annots:
        return annots
    try:
        r = _renderer(fig)
    except Exception:
        return annots
    px_to_pt = 72.0 / fig.dpi
    boxes = []
    for a in annots:
        try:
            boxes.append(a.get_window_extent(renderer=r))
        except Exception:
            boxes.append(None)
    order = sorted((i for i, b in enumerate(boxes) if b is not None),
                   key=lambda i: boxes[i].y0 if axis == 'y' else boxes[i].x0)
    last = None
    for i in order:
        bb = boxes[i]
        lo, hi = (bb.y0, bb.y1) if axis == 'y' else (bb.x0, bb.x1)
        if last is not None and lo < last + pad_px:
            shift = last + pad_px - lo
            dx, dy = annots[i].xyann
            annots[i].xyann = (dx, dy + shift * px_to_pt) if axis == 'y' \
                else (dx + shift * px_to_pt, dy)
            hi += shift
        last = hi
    _grow_axis(ax, annots, axis)
    return annots


def _emptiest_corner(ax, data=None):
    """Which axes corner holds the fewest data points (default: upper right)."""
    corners = {'upper right': (0.62, 1.0, 0.66, 1.0), 'upper left': (0.0, 0.38, 0.66, 1.0),
               'lower right': (0.62, 1.0, 0.0, 0.34), 'lower left': (0.0, 0.38, 0.0, 0.34)}
    if data is None:
        return 'upper right'
    x, y = (np.asarray(a, dtype=float) for a in data)
    ok = np.isfinite(x) & np.isfinite(y)
    if not ok.any():
        return 'upper right'
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    fx = (x[ok] - x0) / (x1 - x0) if x1 != x0 else np.zeros(ok.sum())
    fy = (y[ok] - y0) / (y1 - y0) if y1 != y0 else np.zeros(ok.sum())
    best, best_n = 'upper right', None
    for name, (a, b, c, d) in corners.items():
        n = int(((fx >= a) & (fx <= b) & (fy >= c) & (fy <= d)).sum())
        if best_n is None or n < best_n:
            best, best_n = name, n
    return best


def place_legend(ax, handles, *, title=None, ncol=None, data=None, y=1.02):
    """Put the legend where it cannot cover data: above the panel if a single row
    fits its width, otherwise inside the corner holding the fewest data points.

    The measurement decides between two PLACEMENTS, not between two column counts.
    Wrapping an over-wide strip to a second row above the axes only moves the
    overflow into the vertical, where it runs past the space caption_block
    reserved and is clipped — four cohort keys above a half-width panel did
    exactly that. The corner is chosen by counting the points that actually fall
    in each one, never by assuming 'upper right'.
    """
    n = ncol or len(handles)
    kw = dict(loc='lower center', bbox_to_anchor=(0.5, y), frameon=False,
              fontsize=plt.rcParams['legend.fontsize'], handletextpad=0.5,
              columnspacing=1.4, borderaxespad=0.0, title=title)
    leg = ax.legend(handles=handles, ncol=n, **kw)
    try:
        fits = (leg.get_window_extent(renderer=_renderer(ax.figure)).width
                <= ax.get_window_extent().width)
    except Exception:
        fits = True
    if fits:
        return leg
    leg.remove()
    return legend_inside(ax, handles, loc=_emptiest_corner(ax, data), title=title,
                         ncol=2 if len(handles) > 3 else 1)


def legend_above(ax, handles, *, title=None, ncol=None, y=1.02, fontsize=None,
                 full_width=False, pad_in=None):
    """Legend ABOVE the axes (off the data region), as a single horizontal strip.

    Deliberately does NOT wrap to a second row. `caption_block` reserves a fixed
    slice of the plot block above the axes (`top_pad`), so extra legend rows grow
    past the canvas and are silently clipped — a 7-key LD legend wrapped to two
    columns rendered as two visible keys. Callers that cannot guarantee the strip
    fits should use place_legend(), which switches PLACEMENT rather than shape.

    full_width : centre the strip on the CANVAS instead of on `ax`. Needed when
        `ax` is narrower than the labels — a strip centred on a quarter-width
        panel runs off the left edge and the first label is clipped mid-word,
        which is how this surfaced. The y position stays in `ax` coordinates so
        the strip still sits directly above that panel.
    """
    if pad_in is not None:
        # A FIXED distance above the axes, not a fraction of its height. With a
        # proportional y a tall panel throws its legend far above itself: in
        # cohort_compare panel (c) is several times taller than (b) and the same
        # y=1.11 put (c)'s legend under (b)'s x-axis label, where it read as (b)'s.
        try:
            h_in = ax.get_window_extent(renderer=_renderer(ax.get_figure())).height \
                   / ax.get_figure().dpi
        except Exception:
            h_in = 1.0
        y = 1.0 + pad_in / max(h_in, 0.2)
    kw = {}
    if full_width:
        kw['bbox_transform'] = blended_transform_factory(
            ax.get_figure().transFigure, ax.transAxes)
    return ax.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, y),
                     ncol=ncol or len(handles), frameon=False,
                     fontsize=fontsize or plt.rcParams['legend.fontsize'],
                     handletextpad=0.5, columnspacing=1.4, borderaxespad=0.0,
                     title=title, **kw)


def legend_inside(ax, handles, *, loc='upper right', title=None, ncol=1, fontsize=None):
    """Compact, lightly-framed legend INSIDE the axes, over an empty corner."""
    fs = fontsize or (plt.rcParams['legend.fontsize'] - 0.5)
    leg = ax.legend(handles=handles, loc=loc, ncol=ncol, fontsize=fs, title=title,
                    frameon=True, framealpha=0.92, borderpad=0.5, labelspacing=0.32,
                    handletextpad=0.5, borderaxespad=0.5)
    fr = leg.get_frame()
    fr.set_edgecolor('#C8CDD2'); fr.set_linewidth(0.7); fr.set_facecolor('white')
    if leg.get_title() is not None and leg.get_title().get_text():
        leg.get_title().set_fontweight('bold')
        leg.get_title().set_fontsize(fs + 0.5)
    leg.set_zorder(20)
    return leg


def _wrap_chars(fig, fs, left, right):
    """Measure how many caption characters actually fit on one figure-wide line.

    A fixed chars-per-inch guess is wrong whenever the type scale or font changes,
    and guessing low silently inflates the caption by extra wrapped lines. Measure
    the real advance width with the renderer and keep a margin for the mathtext
    runs (β, subscripts) that are wider than plain text.
    """
    probe = ('OLS beta and the depth-adjusted Case-Control burden difference across '
             'the full minAC range, then a steep climb.')
    try:
        r = fig.canvas.get_renderer()
        t = fig.text(0, 0, probe, fontsize=fs)
        per_char = (t.get_window_extent(renderer=r).width / fig.dpi) / len(probe)
        t.remove()
        usable = fig.get_size_inches()[0] * (right - 0.012)   # canvas x=0.012 → right margin
        return max(80, int(0.93 * usable / per_char))         # 0.93: slack for mathtext runs
    except Exception:
        return max(80, int(fig.get_size_inches()[0] * 12.5))


SEP_DOT  = '  \u00b7  '   # separator between parts of one identity string
HEADER_H = 0.24        # inches reserved for a caption_block(header=…) row


def rsid_text(row):
    """How an rsID is DISPLAYED, given an annotation row. The only such rule.

    `variant_annot.lookup_rsids` attributes an accession to a variant only when
    the source VCF resolved it to that ALLELE. When it did not, `rsID` is '.' and
    the candidates sit in `rsID_unresolved` — saying so is the honest line,
    because the variant id is the identifier that does hold. Several accessions
    on one allele is a genuine dbSNP merge: both name the variant, so both are
    shown.
    """
    rs = str(row.get('rsID', '') or '')
    if rs not in ('.', 'nan', ''):
        return ' / '.join(rs.split(';'))
    un = str(row.get('rsID_unresolved', '') or '')
    return 'rsID unresolved' if un not in ('.', 'nan', '') else ''


def lead_identity(annotation, peak_id, lead_id):
    """(gene, variant_id, rsid) for one locus, for `caption_block(header=…)`.

    Keyed on `peak_id`, not on the variant: a variant that leads two peaks would
    make a variant key ambiguous, while `peak_id` is one row per peak by
    construction. Degrades to the variant id alone rather than failing — a
    missing annotation must cost a figure its gene name, not its existence.
    """
    if not annotation or not Path(annotation).exists():
        return (None, lead_id, '')
    try:
        d = pd.read_csv(annotation, sep='\t', dtype=str)
        r = d[d['peak_id'].astype(str) == str(peak_id)]
        if not len(r):
            return (None, lead_id, '')
        row = r.iloc[0]
        gene = str(row.get('Gene', '') or '')
        return (gene if gene not in ('.', 'nan', '') else None,
                lead_id, rsid_text(row))
    except Exception:
        return (None, lead_id, '')


def caption_block(fig, *, plot_h, title=None, panels=None,
                  notes=None, letters='abcdefgh', wrap=None,
                  top_pad=0.30, left='auto', right=0.985,
                  wspace=0.26, hspace=0.42, extra_bottom=0.0, cap_fontsize=None,
                  margin_axes=None, header=None):
    """Journal-style caption rendered BELOW the axes. THREE lines of content only:

        <title>          bold declarative sentence — the figure's finding
        a, <panels[0]>   bold run-in panel letter + ONE sentence on what it shows
        Data: <notes>    one line of provenance

    Deliberately nothing else. The interpretation paragraph, the symbol
    definitions and the estimator used to live here too, and the caption grew to
    63-70 % of the canvas — a scan figure was 16.2 in tall of which 10.2 in
    was text. These figures are laid out as slides, so the plot has to be the
    figure. All three blocks moved to the sidecar `.md`
    (`figure_doc.write_doc(interpretation=…, defs=…, model=…)`), which already
    carried the value table, the reading order and the limits; nothing was deleted.

    THE PLOT BLOCK IS A FIXED PHYSICAL SIZE. `plot_h` is the height in inches of
    the whole croppable block, measured from the top of the canvas: panel titles,
    the axes themselves, and the x tick labels and x-axis label beneath them. Crop
    `plot_h` inches off the top of the PNG and you have the publication figure.

        figure height = plot_h + rule gap + caption height + pad

    `top_pad` is the slice of plot_h reserved above the axes for panel titles; pass
    a larger value when the figure also carries a suptitle.

    `header=(gene, variant_id, rsid)` names the variant the figure is ABOUT, on one
    line at the very top. It belongs here rather than in the caption because the
    caption is below the crop line — a reader who crops the publication figure
    would otherwise be holding a locus plot that never says which locus it is. The
    row adds HEADER_H inches to the croppable block and shifts the axes down by the
    same amount, so panels keep the size they had before and callers do not retune
    `plot_h` or `top_pad`.

    `left='auto'` measures the y tick labels and y-axis label and fits the margin to
    them (see fit_left_margin); pass a float only to override a measurement that is
    known to be wrong. Caption type is deliberately NOT tied to the plot scale — the
    caption is read on screen, not printed.
    """
    import re
    import textwrap
    auto_left = (left == 'auto')
    if auto_left:
        left = fit_left_margin(fig, margin_axes)
    fw, _fh0 = fig.get_size_inches()
    fs = cap_fontsize if cap_fontsize is not None else max(8.5, plt.rcParams['legend.fontsize'])
    if wrap is None:
        wrap = _wrap_chars(fig, fs, left, right)

    def _wrap(txt, width=None):
        """Wrap WITHOUT ever breaking inside a $…$ mathtext span.

        textwrap breaks on whitespace, and mathtext such as $p \\geq \\alpha$ contains
        spaces — a break there renders as literal '$p' / '\\geq \\alpha$' garbage. Mask
        the spaces inside math spans, wrap, then restore them.
        """
        masked = re.sub(r'\$[^$]*\$', lambda m: m.group(0).replace(' ', '\x00'), txt)
        return [ln.replace('\x00', ' ') for ln in textwrap.wrap(masked, width or wrap)]

    def _runin(label, txt, sep, width=None, indent=''):
        """One bold run-in entry ('a, …' or 'Data: …'), wrapped to figure width."""
        lines = _wrap(f'{label}{sep}{txt}', width) or [f'{label}{sep}'.rstrip()]
        lines[0] = f'$\\mathbf{{{label}}}${sep}' + lines[0][len(label) + len(sep):]
        return [lines[0]] + [indent + ln for ln in lines[1:]]

    def _mathbf(s):
        """Bold a plain-text label without letting mathtext mangle it.

        Inside a $...$ span mathtext DROPS literal spaces and renders '-' as a
        minus, so $\\mathbf{genome-wide significance}$ comes out as
        'genome-widesignificance'. Bold each alphanumeric run separately and leave
        the separators outside the math spans, where they render as themselves.
        """
        return ''.join(f'$\\mathbf{{{p}}}$' if p and p.isalnum() else p
                       for p in re.split(r'([^A-Za-z0-9])', s))


    def _build_body(width):
        """Caption lines, top → bottom, at a given wrap width."""
        b = []
        for i, p in enumerate(panels or []):
            b.extend(_runin(letters[i], p, ', ', width))   # Nature panel run-in: "a, …"
        if notes:
            b.append('')
            b.extend(_runin('Data', notes, ': ', width))
        return b

    def _fs(key, dflt):
        v = plt.rcParams.get(key, dflt)
        return v if isinstance(v, (int, float)) else dflt

    def _height_in(artist):
        """Rendered height of a text artist, in inches — a font metric, so it does
        not change when the figure is resized. That is what makes the two-pass
        layout below stable."""
        try:
            bb = artist.get_window_extent(renderer=fig.canvas.get_renderer())
            return bb.height / fig.dpi
        except Exception:
            return 0.0

    try:
        fig.set_layout_engine('none')           # deterministic margins from here
    except Exception:
        pass

    # ── Pass 1: draw the caption, measure what it actually needs (in inches) ──
    max_w_in = fw * (right - 0.012)          # usable text width, canvas x=0.012 → right

    def _fit(render, fontsize, weight, color, seed):
        """Draw text, then shrink the wrap width until it really fits the canvas.

        A chars-per-inch estimate cannot be trusted across weights and mathtext —
        bold runs ~15 % wider, which is what pushed the bold caption title off the
        right edge. `render(width)` returns the lines at a candidate width; measure,
        re-wrap, repeat. Converges in one or two passes.
        """
        w, t = seed, None
        for _ in range(6):
            if t is not None:
                t.remove()
            t = fig.text(0.012, 0.0, '\n'.join(render(w)), ha='left', va='bottom',
                         fontsize=fontsize, fontweight=weight, color=color, linespacing=1.5)
            try:
                wid = t.get_window_extent(renderer=fig.canvas.get_renderer()).width / fig.dpi
            except Exception:
                return t
            if wid <= max_w_in or w <= 40:
                return t
            w = max(40, int(w * max_w_in / wid) - 1)
        return t

    t_body = t_title = None
    if _build_body(wrap):
        t_body = _fit(_build_body, fs, 'normal', '#333333', wrap)
    if title:
        # The title is set larger AND bold, so it needs its own narrower wrap width —
        # left unwrapped it silently ran off the right edge of the canvas.
        t_title = _fit(lambda w: _wrap(title, w), fs + 1.5, 'bold', '#111111',
                       max(40, int(wrap * fs / (fs + 1.5))))

    pad_in = 0.09                                        # canvas edge → caption text
    gap_in = 0.55 * fs / 72.0                            # body ↔ title
    sep_in = 0.10                                        # caption ↔ rule ↔ plot block
    # x tick labels and the x-axis label live BELOW the axes rectangle and are part
    # of the figure, not of the caption — they must be inside plot_h.
    xlab_in = (_fs('axes.labelsize', 8) + _fs('xtick.labelsize', 7) + 8) * 1.5 / 72.0
    h_body = _height_in(t_body) if t_body is not None else 0.0
    h_title = _height_in(t_title) if t_title is not None else 0.0

    cap_in = pad_in + h_body + ((gap_in + h_title) if t_title is not None else 0.0)

    # ── Size the canvas around a FIXED plot block ────────────────────────────
    # The header row GROWS the block instead of eating into it: subtracting it
    # from top_pad would silently shrink every panel that gained a header.
    head_in = HEADER_H if header else 0.0
    below_in = cap_in + 2 * sep_in + extra_bottom
    fh = plot_h + head_in + below_in
    fig.set_size_inches(fw, fh, forward=True)

    if t_body is not None:
        t_body.set_position((0.012, pad_in / fh))
    if t_title is not None:
        t_title.set_position((0.012, (pad_in + h_body + gap_in) / fh))

    sep_y = (cap_in + sep_in) / fh                      # thin rule under the panels
    bottom = (below_in + xlab_in) / fh
    top = 1.0 - (top_pad + head_in) / fh
    try:
        fig.subplots_adjust(left=left, right=right, top=top, bottom=bottom,
                            wspace=wspace, hspace=hspace)
    except Exception:
        pass
    # SECOND margin pass. The first fit measures tick labels as they stand before
    # subplots_adjust; moving the axes changes the tick locator, so a panel whose
    # labels widen (5.0 -> 12.5) can end up needing more room than was reserved
    # and its y-axis label is clipped off the canvas. Re-measure and grow only.
    if auto_left:
        try:
            need = fit_left_margin(fig, margin_axes)
            if need > left + 0.002:
                fig.subplots_adjust(left=need, right=right, top=top, bottom=bottom,
                                    wspace=wspace, hspace=hspace)
                left = need
        except Exception:
            pass
    # Drawn last: it is anchored to `left`, which the second margin pass may have
    # moved. Segments are placed left to right by MEASURING each one, because a
    # single text artist cannot carry two colours and two weights.
    if header:
        gene, vid, rsid = (list(header) + ['', '', ''])[:3]
        segs = []
        if gene and str(gene) not in ('.', 'nan'):
            segs.append((str(gene), dict(fontsize=fs, fontstyle='italic',
                                         fontweight='bold', color=INK)))
        if vid:
            segs.append((f'{SEP_DOT}{vid}' if segs else str(vid),
                         dict(fontsize=fs - 0.5, color=INK_SOFT)))
        if rsid:
            segs.append((f'{SEP_DOT}{rsid}' if segs else str(rsid),
                         dict(fontsize=fs - 0.5, color=INK_SOFT)))
        x_in = left * fw
        y = 1.0 - (head_in / 2.0) / fh
        for txt, kw in segs:
            a = fig.text(x_in / fw, y, txt, ha='left', va='center', **kw)
            try:
                x_in += a.get_window_extent(
                    renderer=fig.canvas.get_renderer()).width / fig.dpi
            except Exception:
                break

    # thin separator rule delineating the figure from its caption (print-figure style)
    fig.add_artist(Line2D([left, right], [sep_y, sep_y], transform=fig.transFigure,
                          color=NEUTRAL, linewidth=0.8, zorder=0))
