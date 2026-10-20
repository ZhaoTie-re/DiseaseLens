#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The figures of one snapshot of one disease (PNG + PDF), in the
#           shared figure style (scripts/plot_style.py, a vendored copy), and
#           one catalog README for the family:
#             evidence_grid      the genes with genetic evidence for the disease:
#                                evidence level (rows, strongest first) x
#                                function (gene theme, columns)
#             condition_map      the disease at the center; the related
#                                conditions in sectors by category; line style =
#                                kind of evidence, node size = amount of evidence
#             shared_regions     conditions x the GWAS regions of the disease:
#                                colocalization (filled) or lead-variant overlap
#                                (open)
#             gene_map           every gene of the evidence ladder, by evidence
#                                level x the disease and the related conditions,
#                                each cell the kinds of evidence (G, G*, C, L),
#                                with how each kind is derived printed beside it
#             literature_genes   the genetic studies of the disease: studies per
#                                gene and result, for genes reported by two or
#                                more studies, grouped by gene theme
#             genetic_studies    every genetic study of the disease by year and
#                                design class, colored by place (country of the
#                                first author's institution), sized by cases
#           A figure with nothing to show (no GWAS region, no gene reported by
#           two studies, no genetic study) is not drawn; the README says so.
#           Names come from the disease profile; text specific to a disease from
#           its text/figures.yaml (a generic sentence when it has none).
#           Presentation only: every mark is computed from the tables.
# Inputs  : the disease profile, the tidy tables (with the ladder and the gene
#           themes of the snapshot), plot_style.py
# Outputs : <name>.png, <name>.pdf, README.md, figures.json
# ---------------------------------------------------------------------------
import argparse
import json
import math
import re
import sys
import textwrap
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import (NO_GENE, a_, load_glosses, load_status, OTHER_GENES, SHARED_SIGNAL, assign_themes, gene_tallies, load_ladder, relates,  # noqa: E402
                       load_profile, load_themes, log, not_in_gene_map, num, read_tsv)

# colors of the categories of related conditions: the profile's clinical categories in their order, then any other
# category of the table (both from one list, extended in turn); the ontology category and the shared signal keep
# their own
CATEGORY_COLORS = ['#C2627A', '#B08A3C', '#7A5C9E', '#2E8B8B', '#6E8B3D', '#C46B2E', '#A0527F', '#5B7DB1', '#8C6D46',
                   '#3F8A6E']
ONTOLOGY_COLOR = '#4D5B6A'
NUMBER_WORDS = ['no', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine']
# short labels of the design classes (config/curation.yaml, design_classes) in the genetic_studies figure
DESIGN_SHORT = {'genome_wide': 'Genome-wide (GWAS, genome\nsequencing, MR)', 'clonal': 'Clonal hematopoiesis',
                'sequencing': 'Exome or gene panel', 'hla': 'HLA typing', 'candidate': 'Candidate gene',
                'other': 'Other designs', 'case_report': 'Case reports'}


def number(n):
    """Numbers below ten in words (the nine regions), larger ones in digits."""
    return NUMBER_WORDS[n] if 0 <= n < len(NUMBER_WORDS) else f'{n:,}'


def label(name, labels):
    """A short, sentence-case label: the disease's labels (figures text, then profile), else Open Targets' trait
    name tidied."""
    if name in labels:
        return labels[name]
    s = re.sub(r'^amount of |^level of | in blood$| measurement$|\(human\) ', '', name)
    s = re.sub(r'\s+', ' ', s).strip() or name
    return s[0].upper() + s[1:]


def wrap(text, width=30):
    """A label of the condition map, broken into lines of at most `width` characters (a long name without a short
    label would otherwise run into the legends)."""
    return textwrap.fill(text, width, break_long_words=False)


def save(fig, out, name, S):
    for ext in ('png', 'pdf'):
        fig.savefig(out / f'{name}.{ext}', dpi=300 if ext == 'png' else None, bbox_inches='tight', pad_inches=0.08)
    S.plt.close(fig)


def by_condition(d, cond, col):
    """GWAS signals of each condition at each region of the disease (Open Targets IDs mapped to the condition):
    number of studies and the best value of `col`."""
    alias = {i: c.node_id for c in cond.itertuples() for i in c.ot_ids.split()}
    d = d.assign(node_id=d.node_id.map(lambda i: alias.get(i, i)), v=pd.to_numeric(d[col]))
    return d.groupby(['node_id', 'disease_gene']).agg(n=('study_id', 'nunique'), best=('v', 'max'))


def regions_kind(rg):
    """(all published, any published, any from Open Targets) of the regions table."""
    pub = rg.source.isin(['', 'published'])
    return bool(len(rg)) and bool(pub.all()), bool(pub.any()), bool((~pub).any())


def gwas_levels(C, rg, L, curated_min, readme=False):
    """What levels 1-2 of the ladder hold, for the grid's footnote and the README: the disease's text (figures.yaml,
    grid_levels / readme_levels) for its GWAS regions, else a sentence built from the regions; and the curated gene
    evidence when a gene has it."""
    s = C.short
    allpub, anypub, anyot = regions_kind(rg)
    gs = rg.code.eq('G*').any()
    txt = C.text('readme_levels' if readme else 'grid_levels', '') if len(rg) else ''
    if not txt and anypub:
        txt = (f'Levels 1–2: the regions of the {s} GWAS as published (PMID {C.gwas_pmid}); G, genome-wide '
               'significant' + ('; G*, genome-wide significant only in a secondary analysis' if gs else '') + '.')
    elif not txt and anyot:
        txt = (f'Level 1: G, a genome-wide significant region of the {s} GWAS in Open Targets, named by its top '
               'locus-to-gene gene.')
    cur = (f'C, curated gene–disease evidence for {s} (Open Targets source score ≥ {curated_min})'
           if L.code.eq('C').any() or L.curated_sources.ne('').any() else '')
    if not txt:
        return (f'Level 1: {cur}; ' if cur else '') + f'{s} has no GWAS region in this snapshot (no G or G*).'
    return txt + (f' Level 1 also holds {cur}.' if cur else '')


# ── evidence_grid ────────────────────────────────────────────────────────────
def evidence_grid(T, out, S, pal, C, themes, curated_min):
    """Evidence level (rows) x function (columns): each gene with genetic evidence for the disease, at its
    strongest level, in the column of its gene theme."""
    plt = S.plt
    L = T['ladder'].copy()
    if not len(L):
        return None
    s = C.short
    L['level'] = L.level.astype(int)
    short = {t['name']: t['short'] for t in themes}
    short[OTHER_GENES] = OTHER_GENES
    cols = [t['name'] for t in themes if (L.theme == t['name']).any()] + \
        ([OTHER_GENES] if (L.theme == OTHER_GENES).any() else []) + \
        sorted(set(L.theme) - {t['name'] for t in themes} - {OTHER_GENES})
    # every level, for every disease; an empty level says why it is empty
    levels = sorted(C.LADDER)
    LH, CW, LW, HH, NW, FS = 0.27, 1.75, 4.1, 0.78, 0.55, 10.5      # sizes in inches; gene font size
    cell = {(lv, c): list(L[(L.level == lv) & (L.theme == c)].itertuples()) for lv in levels for c in cols}
    tag = lambda g: g.symbol  # noqa: E731
    probe = plt.figure()                                               # gene labels flow left to right, wrapping
    rend = probe.canvas.get_renderer()                                 # at the cell edge: measure them once
    wid = {tag(g): probe.text(0, 0, tag(g), fontsize=FS, style='italic').get_window_extent(rend).width / probe.dpi
           for g in L.itertuples()}
    plt.close(probe)

    def flow(genes):
        out, x, line = [], 0.0, 0
        for g in genes:
            w = wid[tag(g)]
            if x and x + w > CW - 0.24:
                x, line = 0.0, line + 1
            out.append((g, x, line))
            x += w + 0.16
        return out
    place = {k: flow(v) for k, v in cell.items()}
    lines = {lv: max([1] + [p[-1][2] + 1 for (lv_, _c), p in place.items() if lv_ == lv and p]) for lv in levels}
    labs = {lv: textwrap.wrap((f'{C.LADDER[lv][0]} · ' if C.LADDER[lv][0] else '') + C.LADDER[lv][1], 31)
            + (textwrap.wrap(f'({C.GLOSS[lv]})', 31) if C.GLOSS.get(lv) else [])
            + (textwrap.wrap(f'({C.STATUS[lv]["brief"]})', 31) if (C.STATUS.get(lv) or {}).get('brief') else [])
            for lv in levels}
    # a row is as tall as its gene lines, or as its label when that runs to three lines or more
    rh = {lv: max(lines[lv] * LH, 0.19 * len(labs[lv]) - 0.11) + 0.2 for lv in levels}
    n_st = int(T['studies'].pmid.nunique())
    automatic = all(t.get('automatic') for t in themes)       # every theme from Reactome (no curated ones)
    mixed = not automatic and any(t.get('automatic') for t in themes)
    foot = ('Each gene sits at its strongest level. ' + gwas_levels(C, T['regions'], L, curated_min) +
            f' Levels 3–7: what the {n_st} genetic studies of {s} in PubMed report '
            f'for the gene, each finding checked against its abstract; an association in two or more studies ranks '
            f'above one, and levels 3–6 count only findings about {s} itself (a gene associated only with an outcome '
            f'within {s} stays at level 7). '
            + ('Columns: the function of the gene, its Reactome top-level pathway (workbook sheet 11). '
               if automatic else 'Columns: the function of the gene: the curated themes, each with its source, and for a '
               'gene they do not list its Reactome top-level pathway (workbook sheet 11). ' if mixed else
               'Columns: the function of the gene; each theme has its source (workbook sheet 11). ')
            + 'East Asian evidence: workbook sheet 9.')
    fl = textwrap.wrap(foot, 215)
    W = LW + CW * len(cols) + NW
    H = HH + sum(rh.values()) + 0.25 + 0.2 * len(fl)
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis('off')
    fill = pal['level']
    for j, c in enumerate(cols):                                       # column headers
        x = LW + j * CW + CW / 2
        n = int((L.theme == c).sum())
        ax.text(x, HH - 0.3, '\n'.join(textwrap.wrap(short.get(c, c), 17, break_long_words=False)), ha='center',
                va='bottom', fontsize=11, fontweight='bold', color=S.INK, linespacing=1.1)
        ax.text(x, HH - 0.08, f'{n} gene{"s" if n != 1 else ""}', ha='center', va='bottom', fontsize=9.5,
                color=S.INK_SOFT)
    ax.text(LW + len(cols) * CW + NW / 2, HH - 0.08, 'Genes', ha='center', va='bottom', fontsize=9.5,
            color=S.INK_SOFT)
    ax.text(LW - 0.1, HH - 0.3, 'Function (gene theme) →', ha='right', va='bottom', fontsize=10.5, color=S.INK_SOFT,
            style='italic')
    ax.text(0.95, HH - 0.08, 'Evidence level ↓ (strongest first)', ha='left', va='bottom', fontsize=10.5,
            color=S.INK_SOFT, style='italic')
    ax.plot([0.05, W - 0.05], [HH, HH], color=S.INK, lw=1.2)
    y, ys = HH, {}
    for lv in levels:
        h = rh[lv]
        ys[lv] = (y, y + h)
        col = fill.get(lv, S.NEUTRAL)
        ax.add_patch(S.matplotlib.patches.FancyBboxPatch((0.95, y + h / 2 - 0.14), 0.36, 0.28,
                                                         boxstyle='round,pad=0,rounding_size=0.05', facecolor=col,
                                                         edgecolor=S.NEUTRAL_D, lw=0.6))
        ax.text(1.13, y + h / 2, str(lv), ha='center', va='center', fontsize=11, fontweight='bold',
                color='white' if lv in pal['dark'] else S.INK)
        ax.text(1.42, y + h / 2, '\n'.join(labs[lv]), ha='left', va='center', fontsize=10.5, color=S.INK,
                linespacing=1.1)
        for j, c in enumerate(cols):
            genes = cell[(lv, c)]
            x0 = LW + j * CW
            if genes:
                ax.add_patch(S.matplotlib.patches.Rectangle((x0 + 0.04, y + 0.05), CW - 0.08, h - 0.1,
                                                            color=col, lw=0, zorder=1))
            for g, gx, line in place[(lv, c)]:
                ax.text(x0 + 0.14 + gx, y + 0.1 + line * LH + LH / 2, tag(g), ha='left', va='center', fontsize=FS,
                        style='italic', color='white' if lv in pal['dark'] else S.INK, zorder=2)
        n = int((L.level == lv).sum())
        ax.text(LW + len(cols) * CW + NW / 2, y + h / 2, str(n), ha='center', va='center', fontsize=11,
                color=S.INK)
        y += h
        ax.plot([LW, LW + len(cols) * CW], [y, y], color=S.NEUTRAL, lw=0.6, zorder=0)
    for j in range(1, len(cols)):
        ax.plot([LW + j * CW] * 2, [HH, y], color=S.NEUTRAL, lw=0.6, zorder=0)
    ax.plot([0.05, W - 0.05], [y, y], color=S.INK, lw=1.2)
    # the two sources of the ladder, the same for every disease
    for (a, b), text in (((1, 2), 'GWAS,\ncurated'), ((3, 7), f'Genetic studies of {s}')):
        if a not in ys or b not in ys:
            continue
        y0, y1 = ys[a][0] + 0.08, ys[b][1] - 0.08
        ax.plot([0.62, 0.55, 0.55, 0.62], [y0, y0, y1, y1], color=S.INK_SOFT, lw=1.0)
        ax.text(0.38, (y0 + y1) / 2, text, rotation=90, ha='center', va='center', fontsize=10.5, color=S.INK_SOFT,
                fontweight='bold')
    for k, t in enumerate(fl):
        ax.text(0.05, y + 0.28 + 0.2 * k, t, ha='left', va='top', fontsize=9.5, color=S.INK_SOFT)
    save(fig, out, 'evidence_grid', S)
    return {'genes': len(L), 'themes': len(cols), 'by_level': {str(k): int((L.level == k).sum()) for k in levels},
            'studies': n_st}


# ── condition_map ────────────────────────────────────────────────────────────
def condition_map(T, out, S, pal, C):
    plt = S.plt
    s = C.short
    cond, rel = T['conditions'], T['relations']
    rel = rel[[relates(t, r) for t, r in zip(rel.evidence_type, rel.relation)]]   # contrary findings do not relate
    nl = num(cond.set_index('node_id').n_shared_regions).fillna(0).astype(int)
    keep = cond[(cond.category != SHARED_SIGNAL) | (cond.node_id.map(nl) >= 2)].copy()
    single = cond[(cond.category == SHARED_SIGNAL) & (cond.node_id.map(nl) == 1)]
    no_region = cond[(cond.category == SHARED_SIGNAL) & (cond.node_id.map(nl) < 1)]
    # the region most single-region traits share (ties: alphabetical); the profile's checks may name the one expected
    vc = sorted(single.shared_regions.value_counts().items(), key=lambda x: (-x[1], x[0]))
    top, n_top = vc[0] if vc else ('', 0)
    want = C.P.checks.get('single_region_top')
    if want and vc and want != top:
        log(f'condition_map: the most frequent single shared region is {top}, not {want} (profile checks)')
    ev = rel.groupby('node_id').evidence_type.apply(set).to_dict()
    ncl = rel[rel.evidence_type == 'clinical'].groupby('node_id').size().to_dict()
    eas = eas_kinds(cond, rel)
    keep['_o'] = keep.category.map(C.cats.index)
    keep['label'] = keep.name.map(lambda n: label(n, C.labels))
    keep = keep.sort_values(['_o', 'label'], key=lambda x: x.str.lower() if x.dtype == object else x)
    groups = [(c, keep[keep.category == c]) for c in C.cats if (keep.category == c).any()]
    if len(single):
        groups.append(('one region', None))                            # the single-region traits, as one node
    n = len(keep) + (1 if len(single) else 0)
    gap = 2.0                                                          # empty slots between sectors
    step = 2 * math.pi / max(n + gap * len(groups), 1)
    # labels fit the room each condition has on the circle (as in the viewer): the slot's arc at the label radius,
    # in points (12.5 in across 3.7 units); the font shrinks to it (6.5 to 10.5 pt), two lines only where two fit,
    # else one line cut with an ellipsis
    def room_at(lim):                                                  # the slot's arc in points, axes ±lim
        return step * 1.17 * (12.5 * 72 * 0.96 / (2 * lim))
    lim = 1.85
    lfs = max(6.5, min(10.5, room_at(lim) / 1.25))
    two = room_at(lim) >= 2.3 * lfs
    if not two:                     # one line per label: a wider frame, so that the longest label stays inside it
        lim = 2.1
        lfs = max(6.5, min(10.5, room_at(lim) / 1.25))

    def fit(text):
        return wrap(text) if two else (text if len(text) <= 32 else text[:31].rstrip() + '…')
    fig, ax = plt.subplots(figsize=(12.5, 12.5))
    ax.set_aspect('equal')
    ax.axis('off')
    R, slot, pos = 1.0, -math.pi / 2 - 0.5 * step, {}
    sectors = []
    for cat, g in groups:
        slot += gap * step
        start = slot
        rows = list(g.itertuples()) if g is not None else [None]
        for x in rows:
            slot += step
            pos[x.node_id if x is not None else '__single'] = slot
        sectors.append((cat, start + 0.5 * step, slot + 0.5 * step, len(rows)))
    for cat, a0, a1, k in sectors:                                     # a faint wedge behind each sector
        col = pal['cat'].get(cat, S.NEUTRAL_D)
        ax.add_patch(S.matplotlib.patches.Wedge((0, 0), 1.07, math.degrees(a0), math.degrees(a1), width=0.8,
                                                color=col, alpha=0.09, lw=0, zorder=0))
    for x in keep.itertuples():                                       # edges, then nodes
        a = pos[x.node_id]
        e = ev.get(x.node_id, set())
        ls = '-' if 'clinical' in e else ('--' if 'ontology' in e else ':')
        lw = 0.8 + 0.45 * ncl.get(x.node_id, 0) if 'clinical' in e else 1.0
        ax.plot([0, R * math.cos(a)], [0, R * math.sin(a)], ls=ls, lw=lw, color=S.INK_SOFT, alpha=0.55, zorder=1)
    for x in keep.itertuples():
        a = pos[x.node_id]
        k = ncl.get(x.node_id, 0) + int(nl.get(x.node_id, 0)) + (1 if 'ontology' in ev.get(x.node_id, ()) else 0)
        size = 45 + 32 * k
        multi = nl.get(x.node_id, 0) >= 2
        ax.scatter(R * math.cos(a), R * math.sin(a), s=size, color=pal['cat'].get(x.category, S.NEUTRAL_D), zorder=3,
                   edgecolors=S.ACCENT if multi else 'white', linewidths=2.6 if multi else 1.0)
        if x.node_id in eas:
            ax.scatter(1.095 * math.cos(a), 1.095 * math.sin(a), s=30, marker='D', color=S.INK, zorder=4,
                       linewidths=0)
        deg = math.degrees(a)
        flip = 90 < (deg % 360) < 270
        ax.text(1.17 * math.cos(a), 1.17 * math.sin(a), fit(x.label), rotation=deg + 180 if flip else deg,
                rotation_mode='anchor', ha='right' if flip else 'left', va='center', fontsize=lfs, color=S.INK)
    if len(single):
        a = pos['__single']
        ax.plot([0, R * math.cos(a)], [0, R * math.sin(a)], ls=':', lw=1.0, color=S.INK_SOFT, alpha=0.55, zorder=1)
        ax.scatter(R * math.cos(a), R * math.sin(a), s=260, color=S.NEUTRAL, zorder=3, edgecolors='white')
        ax.text(R * math.cos(a), R * math.sin(a), str(len(single)), ha='center', va='center', fontsize=9.5,
                color=S.INK, zorder=4)
        deg = math.degrees(a)
        flip = 90 < (deg % 360) < 270
        ax.text(1.17 * math.cos(a), 1.17 * math.sin(a), f'One {s} region only', rotation=deg + 180 if flip else deg,
                rotation_mode='anchor', ha='right' if flip else 'left', va='center', fontsize=10.5,
                color=S.INK_SOFT, style='italic')
    ax.scatter(0, 0, s=5200, color=S.DATA_DARK, zorder=5)
    ax.text(0, 0, s, ha='center', va='center', color='white', fontsize=17 if len(s) <= 6 else max(9, 102 / len(s)),
            fontweight='bold', zorder=6)
    L = S.Line2D
    shared = bool(len(T['regions'])) or bool((nl >= 1).any())         # a shared signal is possible at all
    # the legend names only what the map draws
    lines = {('clinical' if 'clinical' in ev.get(n, ()) else 'ontology' if 'ontology' in ev.get(n, ()) else 'shared')
             for n in keep.node_id}
    handles = ([L([0], [0], color=S.INK_SOFT, lw=2.2, ls='-', label='Clinical evidence (width: items)')]
               if 'clinical' in lines else []) + \
        ([L([0], [0], color=S.INK_SOFT, lw=1.2, ls='--', label='Disease ontology only')] if 'ontology' in lines else []) + \
        ([L([0], [0], color=S.INK_SOFT, lw=1.2, ls=':', label='Shared GWAS signal only')]
         if 'shared' in lines or len(single) else []) + \
        ([L([0], [0], marker='o', ls='', markersize=11, markerfacecolor='white', markeredgecolor=S.ACCENT,
            markeredgewidth=2.4, label=f'Shares ≥ 2 {s} GWAS regions')] if (keep.node_id.map(nl) >= 2).any() else []) + \
        ([L([0], [0], marker='D', ls='', markersize=6, color=S.INK, label='East Asian evidence of the relation')]
         if set(eas) & set(keep.node_id) else []) + \
        [L([0], [0], marker='o', ls='', markersize=8, color=S.NEUTRAL_D, label='Node size: amount of evidence')]
    if len(single):
        handles.append(S.matplotlib.patches.Patch(color=S.NEUTRAL, label=(
            f'{len(single)} trait{"s" if len(single) != 1 else ""} sharing one {s} region'
            + (f' ({n_top} of them {top})' if len(single) > 1 else f' ({top})') + ': weak')))
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.075), ncol=3, frameon=False,
               fontsize=10.5, handlelength=2.6, columnspacing=1.4)
    names = {SHARED_SIGNAL: f'Shared GWAS signal only (≥ 2 {s} regions)'}
    cats = [S.matplotlib.patches.Patch(color=pal['cat'].get(c, S.NEUTRAL_D), label=f'{names.get(c, c)} ({k})')
            for c, _a0, _a1, k in sectors if c != 'one region']
    fig.legend(handles=cats, loc='lower center', bbox_to_anchor=(0.5, 0.94), ncol=3, frameon=False, fontsize=11,
               columnspacing=1.6)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.93, bottom=0.08)
    save(fig, out, 'condition_map', S)
    # the same map as data, for the interactive map of the viewer: what is drawn, where, and how
    items = []
    for x in keep.itertuples():
        e = ev.get(x.node_id, set())
        items.append({'id': x.node_id, 'label': x.label, 'name': x.name, 'category': x.category,
                      'angle': round(pos[x.node_id], 5),
                      'line': 'clinical' if 'clinical' in e else ('ontology' if 'ontology' in e else 'shared'),
                      'n_clinical': int(ncl.get(x.node_id, 0)), 'n_shared_regions': int(nl.get(x.node_id, 0)),
                      'ontology': 'ontology' in e,
                      'evidence': int(ncl.get(x.node_id, 0) + int(nl.get(x.node_id, 0)) + ('ontology' in e)),
                      'east_asian': sorted(eas.get(x.node_id, ()))})
    names = {SHARED_SIGNAL: f'Shared GWAS signal only (≥ 2 {s} regions)'}
    data = {'short': s, 'core': C.P.id,
            'colors': {'core': S.DATA_DARK, 'multi': S.ACCENT, 'line': S.INK_SOFT, 'east_asian': S.INK,
                       'single': S.NEUTRAL},
            'sectors': [{'category': c, 'label': names.get(c, c), 'color': pal['cat'].get(c, S.NEUTRAL_D),
                         'a0': round(a0, 5), 'a1': round(a1, 5), 'n': k} for c, a0, a1, k in sectors
                        if c != 'one region'],
            'items': items,
            'single': ({'n': len(single), 'angle': round(pos['__single'], 5), 'top': top, 'n_top': int(n_top),
                        'ids': sorted(single.node_id)} if len(single) else None),
            'not_shown': {'shared_signal_without_region': len(no_region)},
            'shared_possible': shared,
            'east_asian_kinds': EAS_KINDS}
    (out / 'condition_map.json').write_text(json.dumps(data, ensure_ascii=False, indent=1))
    return {'conditions_shown': len(keep), 'single_region_traits': len(single), 'single_region_top': int(n_top),
            'single_region_top_name': top, 'shared_signal_without_region': len(no_region),
            'east_asian': len(set(eas) & set(keep.node_id))}


# the kinds of East Asian evidence of a relation (the diamond of the condition map): the evidence that relates the
# condition to the disease, never the condition's other studies
EAS_KINDS = {'clinical': 'clinical evidence from East Asia (a study done there, or an East Asian population it '
                            'names)',
             'shared signal': 'a shared GWAS signal whose two GWAS (of the disease and of the condition) have mostly '
                              'East Asian samples'}


def eas_kinds(cond, rel):
    """{condition: {kind of East Asian evidence of its relation}}, from the relations' own flag (tidy_graph.py)."""
    out = {}
    for x in rel[rel.east_asian == 'True'].itertuples():
        out.setdefault(x.node_id, set()).add('clinical' if x.evidence_type == 'clinical' else 'shared signal')
    return out


# ── shared_regions ───────────────────────────────────────────────────────────
def shared_rows(cond, cats):
    """The conditions of the regions figure: sharing two or more regions, or clinically related and sharing one."""
    cond = cond.copy()
    cond['nreg'] = num(cond.n_shared_regions).fillna(0).astype(int)
    rows = cond[(cond.nreg >= 2) | ((cond.category != SHARED_SIGNAL) & (cond.nreg >= 1))].copy()
    rows['_o'] = rows.category.map(cats.index)
    return cond, rows


def shared_regions(T, out, S, pal, C, clpp_min):
    plt = S.plt
    s = C.short
    rg = T['regions']
    if not len(rg):
        return None
    cond, rows = shared_rows(T['conditions'], C.cats)
    rows['label'] = rows.name.map(lambda n: label(n, C.labels))
    rows = rows.sort_values(['_o', 'nreg', 'label'], ascending=[True, False, True])
    ca, oa = by_condition(T['coloc'], cond, 'clpp'), by_condition(T['overlap'], cond, 'pip')
    regs = list(rg.assign(_c=rg.code.eq('G*')).sort_values('_c', kind='stable').itertuples())
    xs = [j + (0.5 if r.code == 'G*' else 0) for j, r in enumerate(regs)]  # a gap between G and G* regions
    nr = max(len(rows), 1)                                             # rows drawn; the axes hold at least 20, so
    na = max(nr, 20)                                                   # that the legends at the right fit
    h = 0.34 * na + 1.6
    w = max(16.5, 1.1 * (xs[-1] + 1.1) / 0.65)                        # at least 1.1 inch per region column
    fig, ax = plt.subplots(figsize=(w, h))
    fig.subplots_adjust(left=0.15, right=0.8, top=1 - 1.25 / h, bottom=0.2 / h)
    S.despine(ax, grid_axis='none')
    ys = {nid: i for i, nid in enumerate(rows.node_id)}
    cmap = S.matplotlib.colors.LinearSegmentedColormap.from_list('clpp', ['#C9DCEC', S.DATA, S.DATA_DARK])
    norm = S.matplotlib.colors.LogNorm(vmin=0.01, vmax=1.0)
    size = lambda n: 30 + 38 * math.sqrt(n)  # noqa: E731
    for x, r in zip(xs, regs):
        if r.disease_study_locus_id:                                   # column lines down to the last row
            ax.axvline(x, ymin=(na - nr) / (na + 0.1), color=S.NEUTRAL, lw=0.8, zorder=0)
        elif r.variant_id:
            ax.axvline(x, ymin=(na - nr) / (na + 0.1), color=S.NEUTRAL, lw=0.8, ls=(0, (2, 2)), zorder=0)
        else:                                                          # nothing to compare with in Open Targets
            ax.add_patch(S.matplotlib.patches.Rectangle((x - 0.42, -0.5), 0.84, nr, color=S.NEUTRAL,
                                                        alpha=0.35, lw=0, zorder=0))
            ax.text(x, (nr - 1) / 2, 'lead variant not in Open Targets', rotation=90, ha='center',
                    va='center', fontsize=10.5, color=S.INK_SOFT, style='italic')
        for nid in rows.node_id:
            if (nid, r.gene) in ca.index:
                v = ca.loc[(nid, r.gene)]
                ax.scatter(x, ys[nid], s=size(v.n), color=cmap(norm(max(v.best, 0.01))), zorder=3,
                           edgecolors='white', linewidths=0.6)
            elif (nid, r.gene) in oa.index:
                v = oa.loc[(nid, r.gene)]
                ax.scatter(x, ys[nid], s=size(v.n), facecolors='white', edgecolors=S.DATA_DARK, linewidths=1.4,
                           zorder=3)
    n_single = int(((cond.nreg == 1) & (cond.category == SHARED_SIGNAL)).sum())   # not rows, but counted
    if not len(rows):
        ax.text((xs[0] + xs[-1]) / 2, 0, (f'No related condition shares two or more {s} regions, and no clinically '
                                          f'related condition shares one; {n_single} trait'
                                          f'{"s" if n_single != 1 else ""} related only by a shared GWAS signal '
                                          f'share{"s" if n_single == 1 else ""} one region.') if n_single else
                f'No related condition shares {a_(s)} region in this snapshot.', ha='center',
                va='center', fontsize=11, color=S.INK_SOFT, style='italic')
    status = lambda r: ('credible set' if r.disease_study_locus_id else 'lead only' if r.variant_id  # noqa: E731
                        else 'not in OT')
    blend = S.matplotlib.transforms.blended_transform_factory(ax.transData, ax.transAxes)
    up = lambda pt: S.matplotlib.transforms.offset_copy(blend, fig=fig, y=pt, units='points')  # noqa: E731
    for x, r in zip(xs, regs):                                         # column headers, drawn so that the region
        head = [r.rsid or f'chr{r.chrom}:{r.pos_grch38 or r.pos_grch37}'] + ([f'tier {r.tier}'] if r.tier else []) \
            + [status(r)]                                              # name can be bold italic
        ax.text(x, 1.0, '\n'.join(head), transform=up(5), ha='center', va='bottom', fontsize=9.5, color=S.INK_SOFT,
                linespacing=1.2)
        ax.text(x, 1.0, r.gene, transform=up(44), ha='center', va='bottom', fontsize=11.5, fontweight='bold',
                style='normal' if r.gene.startswith('chr') else 'italic', color=S.INK)   # a position: upright
    ax.set_xticks([])
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows.label, fontsize=11)
    ax.set_ylim(na - 0.5, -0.6)
    if na > nr:                                                        # few rows: the frame ends with the last row
        ax.spines['left'].set_bounds(nr - 0.5, -0.6)
        ax.spines['bottom'].set_visible(False)
    ax.set_xlim(-0.6, xs[-1] + 0.5)
    ax.tick_params(length=0)
    for code, text in (('G', C.text('regions_g', '') or f'G · genome-wide significant in {a_(s)} GWAS'),
                       ('G*', C.text('regions_g_star', '') or 'G* · significant only in a secondary analysis')):
        sel = [x for x, r in zip(xs, regs) if r.code == code]
        if not sel:
            continue
        ax.plot([min(sel) - 0.4, max(sel) + 0.4], [1.0, 1.0], transform=up(64), color=S.INK, lw=1.0, clip_on=False)
        ax.text((min(sel) + max(sel)) / 2, 1.0, text, transform=up(67), ha='center', va='bottom', fontsize=11,
                fontweight='bold', color=S.INK)
    prev = None                                                        # category bands
    for i, x in enumerate(rows.itertuples()):
        if x.category != prev:
            k = (rows.category == x.category).sum()
            ax.add_patch(S.matplotlib.patches.Rectangle((-0.6, i - 0.5), 0.0125 * (xs[-1] + 1.1), k,
                                                        color=pal['cat'].get(x.category, S.NEUTRAL_D),
                                                        clip_on=False, transform=ax.transData, lw=0))
            if prev is not None:
                ax.axhline(i - 0.5, color=S.NEUTRAL_D, lw=0.7)
            prev = x.category
    sm = S.matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.03, shrink=0.3, anchor=(0, 0.98))
    cb.set_label('Strongest CLPP (●)', fontsize=11)
    cb.ax.tick_params(labelsize=10)
    L = S.Line2D
    kinds = [L([0], [0], marker='o', ls='', markersize=10, color=S.DATA,
               label=f'Colocalizes with the\n{s} credible set')]
    if any(r.variant_id and not r.disease_study_locus_id for r in regs):
        kinds.append(L([0], [0], marker='o', ls='', markersize=10, markerfacecolor='white',
                       markeredgecolor=S.DATA_DARK, markeredgewidth=1.4,
                       label=f'Its credible set holds the\n{s} lead variant\n(PIP ≥ {clpp_min}; weaker)'))
    b0, t0 = 0.2 / h, 1 - 1.25 / h                                     # legends at figure level (a tight bounding
    fy = lambda ya: b0 + ya * (t0 - b0)  # noqa: E731                    # box drops axes-level extra legends)
    fig.legend(handles=kinds, loc='upper left', bbox_to_anchor=(0.815, fy(0.68)), frameon=False, fontsize=10.5,
               labelspacing=1.1, borderpad=0.2, title='Shared region', title_fontsize=11, alignment='left')
    sizes = [L([0], [0], marker='o', ls='', markersize=math.sqrt(size(k)), color=S.NEUTRAL_D, label=f'{k}')
             for k in (1, 5, 20)]
    fig.legend(handles=sizes, title='GWAS signals', loc='upper left', bbox_to_anchor=(0.815, fy(0.43)),
               frameon=False, fontsize=10.5, title_fontsize=11, labelspacing=1.3, borderpad=0.2, alignment='left')
    handles = [S.matplotlib.patches.Patch(color=pal['cat'].get(c, S.NEUTRAL_D), label=c) for c in C.cats
               if (rows.category == c).any()]
    if handles:
        fig.legend(handles=handles, title='Category (bar at left)', loc='upper left',
                   bbox_to_anchor=(0.815, fy(0.17)), frameon=False, fontsize=10.5, title_fontsize=11,
                   alignment='left', borderpad=0.2)
    save(fig, out, 'shared_regions', S)
    return {'conditions': len(rows), 'regions': [r.gene for r in regs], 'single_region_traits': n_single,
            'colocalisation_cells': int(sum((n, r.gene) in ca.index for n in rows.node_id for r in regs)),
            'overlap_cells': int(sum((n, r.gene) not in ca.index and (n, r.gene) in oa.index
                                     for n in rows.node_id for r in regs))}


# ── gene_map ─────────────────────────────────────────────────────────────────
def gene_map(T, out, S, pal, C, l2g_min, curated_min):
    """Every gene of the evidence ladder, by level, x the disease and the related conditions with GWAS or curated
    evidence for any of them; each cell lists the kinds of evidence (G, G*, C, L), with how each kind is derived
    beside the map. Genes shared with no related condition stay in the map (an empty row)."""
    plt = S.plt
    s, core = C.short, C.P.id
    gc, cond, L = T['gene_condition'], T['conditions'], T['ladder'].copy()
    if not len(L):
        return None
    rank = lambda c: C.CODES.index(c) if c in C.CODES else len(C.CODES)  # noqa: E731
    L['level'] = L.level.astype(int)
    gc = gc.assign(code=gc.evidence.map(C.KIND).fillna(''))
    gc = gc[gc.code.ne('')]
    rel = gc[(gc.node_id != core) & gc.symbol.isin(set(L.symbol))]
    L['nrel'] = L.symbol.map(rel.groupby('symbol').node_id.nunique()).fillna(0).astype(int)
    L = L.sort_values(['level', 'nrel', 'symbol'], ascending=[True, False, True])
    genes, level, nrel = list(L.symbol), dict(zip(L.symbol, L.level)), dict(zip(L.symbol, L.nrel))
    cat = dict(zip(cond.node_id, cond.category))
    names = dict(zip(cond.node_id, cond.name))
    order = {c: i for i, c in enumerate(C.cats)}
    cols = (rel.groupby('node_id').symbol.nunique().rename('n').reset_index()
            .assign(_o=lambda d: d.node_id.map(lambda n: order.get(cat.get(n), len(order)))))
    cols = cols.sort_values(['_o', 'n'], ascending=[True, False])
    nodes = [core] + list(cols.node_id)
    cell = gc[gc.symbol.isin(genes) & gc.node_id.isin(nodes)].groupby(['symbol', 'node_id']).code.apply(
        lambda x: ''.join(sorted(set(x), key=rank))).to_dict()
    # the related conditions that are no column, by reason (the same counts in the workbook and the knowledge base)
    hidden, hidden_text = not_in_gene_map(cond, cols.node_id, 'the regions figure' if len(T['regions']) else '')
    ng, nn = len(genes), len(nodes)
    H = 0.3 * ng + 3.6
    fig = plt.figure(figsize=(17.5, H))
    ax = fig.add_axes([0.11, 0.04, 0.52, 0.78])
    ax.set_xlim(-1.5, nn + 0.6)
    ax.set_ylim(ng - 0.5, -0.5)
    for i, g in enumerate(genes):
        lv = level[g]
        ax.add_patch(S.matplotlib.patches.Rectangle((-1.45, i - 0.42), 0.9, 0.84,
                                                    facecolor=pal['level'].get(lv, S.NEUTRAL),
                                                    edgecolor=S.NEUTRAL_D, lw=0.4))
        ax.text(-1.0, i, str(lv), ha='center', va='center', fontsize=8.5, color=S.INK)
        for j, n in enumerate(nodes):
            k = cell.get((g, n), '')
            if not k:
                if j == 0:                                             # levels 6-7: studied, no positive result
                    ax.text(j, i, '–', ha='center', va='center', fontsize=9, color=S.INK_SOFT)
                continue
            fc = pal['kind'].get(k, pal['multi'])
            ax.add_patch(S.matplotlib.patches.Rectangle((j - 0.45, i - 0.42), 0.9, 0.84, color=fc, lw=0))
            ax.text(j, i, k, ha='center', va='center', fontsize=8.5 if len(k) < 3 else 7.5,
                    color=S.INK if k in pal['kind'] else 'white', fontweight='bold')
        ax.text(nn + 0.15, i, str(nrel[g]), ha='center', va='center', fontsize=9,
                color=S.INK if nrel[g] else S.INK_SOFT)
    ax.set_xticks([-1] + list(range(nn)) + [nn + 0.15])
    ax.set_xticklabels(['Evidence level'] + [s] + [label(names.get(n, n), C.labels) for n in cols.node_id] +
                       ['Related conditions'], rotation=90, fontsize=10.5)
    ax.xaxis.tick_top()
    ax.set_yticks(range(ng))
    ax.set_yticklabels(genes, fontsize=10.5, fontstyle='italic')
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.axvline(0.5, color=S.INK, lw=1.0)
    ax.axvline(-0.5, color=S.INK, lw=1.0)
    ax.axvline(nn - 0.4, color=S.INK, lw=0.6)
    prev, j0 = None, 1                                                 # category bars below the columns
    for j, n in enumerate(list(cols.node_id) + [None], start=1):
        c = cat.get(n, '') if n else None
        if c != prev and prev is not None:
            ax.add_patch(S.matplotlib.patches.Rectangle((j0 - 0.5, ng - 0.5 + 0.25), j - j0, 0.35,
                                                        color=pal['cat'].get(prev, S.NEUTRAL_D), clip_on=False,
                                                        lw=0))
            if c is not None:
                ax.axvline(j - 0.5, color=S.NEUTRAL_D, lw=0.7)
            j0 = j
        prev = c
    for i in range(ng):
        same = i + 1 < ng and level[genes[i + 1]] == level[genes[i]]
        ax.axhline(i + 0.5, color=S.NEUTRAL if same else S.INK_SOFT, lw=0.4 if same else 1.0, zorder=0)
    tx = fig.add_axes([0.66, 0.04, 0.33, 0.78])
    tx.axis('off')
    inch = 1 / (H * 0.78)                                              # one inch of the panel, in axes fraction
    rg = T['regions']
    _allpub, anypub, anyot = regions_kind(rg)
    has_gs = 'G*' in set(L.code) or any('G*' in k for k in cell.values())
    note = C.text('gene_map_g_note', '') or (f'(PMID {C.gwas_pmid})' if C.gwas_pmid else '')
    g_core = (f'the gene the publication names for a genome-wide significant region {note}'.strip() if anypub else
              f'the top locus-to-gene gene (score ≥ {l2g_min}) of a genome-wide significant region of its GWAS in '
              'Open Targets' if anyot else 'none in this snapshot (no GWAS region)')
    lines = [('How each cell is derived', 'h'),
             (f'G   GWAS. {s}: {g_core}. Related condition: the top locus-to-gene gene (score ≥ {l2g_min}) of at '
              'least one of its GWAS credible sets in Open Targets.', 'G')]
    if has_gs:
        lines.append((f'G*  GWAS, {s} only: ' + (C.text('gene_map_g_star', '') or 'a region genome-wide significant '
                                                 'only in a secondary analysis of the published GWAS; marked apart.'),
                      'G*'))
    lines += [('C   Curated: curated gene–disease evidence for the condition or its descendants, from a source whose '
               f'Open Targets score for the gene is ≥ {curated_min} (ClinGen, PanelApp, Orphanet, ClinVar, UniProt, '
               'Gene2Phenotype, gene burden).', 'C'),
              (f'L   Literature, {s} only: a genetic study of {s} in PubMed reports an association, an enrichment '
               'or a variant in the gene; the sentence is checked against the abstract.', 'L'),
              (f'Gray: two or more kinds (letters shown). –: studied in {s} without a positive result '
               '(levels 6–7).', ''),
              ('', ''),
              ('Rows and columns', 'h'),
              (f'Rows: all {ng} genes studied genetically in {s}, by evidence level (square at the left: '
               + ', '.join(f'{k} {(c + " ") if c else ""}{n[:1].lower() + n[1:]}' for k, (c, n) in C.LADDER.items())
               + '), then by the number of related '
               f'conditions in which the gene has G or C evidence (right). {sum(1 for g in genes if not nrel[g])} '
               f'genes have no such condition: they are specific to {s} in this map.', ''),
              (f'Columns: {s}, then the {nn - 1} related conditions with G or C evidence for at least one of these '
               f'genes, by category (bar below) and number of genes. ' + hidden_text, ''),
              ('The evidence behind each cell is in sheet 6 of the workbook.', '')]
    y = 1.0
    for t, k in lines:
        if not t:
            y -= 0.3 * inch
            continue
        if k == 'h':
            tx.text(0, y, t, fontsize=12.5, fontweight='bold', va='top', color=S.INK)
            y -= 0.38 * inch
            continue
        w = textwrap.fill(t, 58)
        if k in pal['kind']:
            w = textwrap.fill(t[4:], 54)
            tx.add_patch(S.matplotlib.patches.Rectangle((0, y - 0.26 * inch), 0.045, 0.26 * inch,
                                                        color=pal['kind'][k], transform=tx.transAxes,
                                                        clip_on=False))
            tx.text(0.0225, y - 0.13 * inch, k, fontsize=10, color=S.INK, fontweight='bold', ha='center',
                    va='center')
            tx.text(0.065, y, w, fontsize=11, va='top', color=S.INK)
        else:
            tx.text(0, y, w, fontsize=11, va='top', color=S.INK_SOFT)
        y -= 0.2 * inch * (w.count('\n') + 1) + 0.18 * inch
    shown = {cat.get(n, '') for n in cols.node_id}
    handles = [S.matplotlib.patches.Patch(color=pal['cat'].get(c, S.NEUTRAL_D), label=c) for c in C.cats
               if c in shown]
    if handles:
        tx.legend(handles=handles, loc='upper left', bbox_to_anchor=(0, y - 0.2 * inch), frameon=False,
                  fontsize=10.5, title='Category of the related condition', title_fontsize=11, alignment='left')
    save(fig, out, 'gene_map', S)
    return {'genes': ng, 'genes_shared': sum(1 for g in genes if nrel[g]), 'conditions': nn - 1,
            'not_shown': {k: n for k, n, _x in hidden}, 'not_shown_text': hidden_text,
            'cells': {k: int(v) for k, v in pd.Series(list(cell.values()), dtype=object).value_counts().items()}}


# ── literature_genes ─────────────────────────────────────────────────────────
def literature_genes(T, out, S, pal, C, themes):
    plt = S.plt
    s = C.short
    fx = assign_themes(T['findings'], themes)
    gt = gene_tallies(fx, C.P.about, C.P.within)
    gt = gt[(gt.n_studies >= 2) & (gt.theme != NO_GENE)].copy()
    if not len(gt):
        return None
    order = {t['name']: i for i, t in enumerate(themes)}
    short = {t['name']: t['short'] for t in themes}
    gt['_t'] = gt.theme.map(lambda t: order.get(t, len(order)))
    gt = gt.sort_values(['_t', 'n_studies', 'label'], ascending=[True, False, True]).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(10.5, 0.42 * len(gt) + 1.8))
    S.despine(ax, grid_axis='x')
    results = [('associated', 'Associated or enriched'), ('within_disease', f'Associated with an outcome within {s}'),
               ('not_associated', 'Not associated'), ('identified', 'Variant identified'), ('other', 'Other result')]
    cols = {'associated': S.ACCENT, 'within_disease': S.ACCENT_LT, 'not_associated': S.NEUTRAL_D,
            'identified': S.DATA, 'other': S.NEUTRAL}
    left = np.zeros(len(gt))
    for key, lab in results:
        v = gt[key].to_numpy(dtype=float)
        ax.barh(range(len(gt)), v, left=left, color=cols[key], label=lab, height=0.68, edgecolor='white', lw=0.6)
        left += v
    ax.set_yticks(range(len(gt)))
    ax.set_yticklabels(gt.label, fontsize=11.5, fontstyle='italic')
    ax.invert_yaxis()
    ax.set_xlabel(f'Genetic studies of {s} (a study can report several results)', fontsize=11.5)
    ax.xaxis.set_major_locator(S.matplotlib.ticker.MaxNLocator(integer=True))
    prev = None
    for i, x in enumerate(gt.itertuples()):
        if x.theme != prev:
            if prev is not None:
                ax.axhline(i - 0.5, color=S.NEUTRAL_D, lw=0.8)
            ax.text(left.max() + 0.3, i - 0.35, short.get(x.theme, x.theme), ha='left', va='top', fontsize=10.5,
                    color=S.INK_SOFT, fontweight='bold')
            prev = x.theme
    ax.set_xlim(0, left.max() + 4.6)
    ax.tick_params(axis='x', labelsize=11)
    ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.01), ncol=3, frameon=False, fontsize=10.5)
    save(fig, out, 'literature_genes', S)
    return {'genes': len(gt)}


# the place groups of a study (country of the first author's institution)
PLACE_GROUPS = [('Japan', ['Japan']), ('China', ['China']), ('Korea, Taiwan', ['Korea', 'Taiwan', 'Hong Kong']),
                ('Europe', ['United Kingdom', 'France', 'Germany', 'Austria', 'Switzerland', 'Czech Republic',
                            'Hungary', 'Italy', 'Spain', 'Netherlands', 'Belgium', 'Denmark', 'Sweden', 'Norway',
                            'Finland', 'Poland', 'Ireland', 'Portugal', 'Greece']),
                ('North America', ['United States', 'Canada']), ('Other', None), ('Not known', None)]


def design_rows(P, labels):
    """The design classes, top to bottom: those of the curation rules (config/curation.yaml, design_classes) in
    their order, the class set by study type (case reports) last, then any other class the table holds; each with its
    short label and whether it is the genome-wide class (its studies are labeled)."""
    dc = (P.curation or {}).get('design_classes') or []
    dc = [c for c in dc if not c.get('study_type')] + [c for c in dc if c.get('study_type')]
    rows = [(c['label'], DESIGN_SHORT.get(c['key'], textwrap.fill(c['label'], 28)), c['key'] == 'genome_wide')
            for c in dc if c['label'] in labels]
    known = {r[0] for r in rows}
    return rows + [(x, textwrap.fill(x, 28), False) for x in sorted(labels - known)]


def genetic_studies(T, out, S, pal, C):
    """Every genetic study of the disease: year (x) x design class (rows); color, the place; size, the cases."""
    plt = S.plt
    st = T['studies'].copy()
    st['year_n'] = pd.to_numeric(st.year, errors='coerce')
    st = st[st.year_n.notna()]
    if not len(st):
        return None
    st['cases'] = pd.to_numeric(st.n_cases, errors='coerce')

    def group(c):
        if not c:
            return 'Not known'
        return next((g for g, cs in PLACE_GROUPS if cs and c in cs), 'Other')
    st['group'] = st.place.map(group)
    color = {'Japan': S.ACCENT, 'China': '#B08A3C', 'Korea, Taiwan': '#D9A03C', 'Europe': S.DATA,
             'North America': '#3F8A6E', 'Other': '#7A5C9E', 'Not known': S.NEUTRAL_D}
    rows = design_rows(C.P, set(st.design_class_label))
    short = {r: sh for r, sh, _g in rows}

    def size(n):
        return 30 + 9 * math.sqrt(n) if n == n else 30
    fig, ax = plt.subplots(figsize=(12.5, 0.85 * len(rows) + 2.6))
    S.despine(ax, grid_axis='x')
    lo, hi = int(st.year_n.min()), int(st.year_n.max())
    pts = []
    for i, (r, _sh, genome_wide) in enumerate(rows):
        sub = st[st.design_class_label == r].sort_values(['year_n', 'pmid'])
        for _y, g in sub.groupby('year_n'):
            offs = (np.arange(len(g)) - (len(g) - 1) / 2) * 0.22
            for x, o in zip(g.itertuples(), offs):
                ax.scatter(x.year_n, i + o, s=size(x.cases), color=color[x.group], alpha=0.85, edgecolor='white',
                           lw=0.6, zorder=3)
                if genome_wide:
                    pts.append((x, i + o))
    # the genome-wide studies cluster in a few recent years: their labels stand in a column to the left of the
    # cluster, in the order of the points' heights, so that the leader lines do not cross
    pts.sort(key=lambda q: (q[1], q[0].year_n))
    x0 = min(x.year_n for x, _y in pts) - 2.5 if pts else hi
    for k, (x, y) in enumerate(pts):
        name = (str(x.first_author).split() or [''])[0]       # a record without an author list: its PMID
        lab = (f'{name} {x.year}' if name else f'PMID {x.pmid}') + (f', {int(x.cases):,} cases' if x.cases == x.cases
                                                                    else '')
        ty = -0.62 + k * (1.05 / max(len(pts) - 1, 1))
        ax.annotate(lab, (x.year_n, y), xytext=(x0, ty), textcoords='data', ha='right', va='center', fontsize=9.5,
                    color=S.INK_SOFT, arrowprops=dict(arrowstyle='-', color=S.NEUTRAL_D, lw=0.7, shrinkA=3,
                                                      shrinkB=4))
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f'{short[r]} ({int((st.design_class_label == r).sum())})' for r, _s, _g in rows], fontsize=11)
    ax.set_ylim(len(rows) - 0.4, -1.0)
    ax.set_xlim(lo - 1.5, hi + 1.5)
    ax.set_xlabel('Year of publication', fontsize=11.5)
    ax.tick_params(axis='x', labelsize=11)
    n_by = st.group.value_counts()
    h1 = [S.Line2D([0], [0], marker='o', lw=0, markersize=9, color=color[g], label=f'{g} ({int(n_by[g])})')
          for g, _c in PLACE_GROUPS if g in n_by]
    h2 = [plt.scatter([], [], s=size(n), color=S.NEUTRAL_D, label=f'{n:,} cases') for n in (10, 100, 1000)]
    # figure-level legends: a tight bounding box keeps them (it can drop an axes-level legend added as an artist)
    fig.legend(handles=h1, loc='upper left', bbox_to_anchor=(1.0, 0.97), frameon=False, fontsize=10.5,
               title="Place (first author's institution)", title_fontsize=10.5, alignment='left')
    fig.legend(handles=h2, loc='lower left', bbox_to_anchor=(1.0, 0.1), frameon=False, fontsize=10.5,
               title='Cases genotyped', title_fontsize=10.5, alignment='left', labelspacing=1.3)
    save(fig, out, 'genetic_studies', S)
    return {'studies': len(st), 'places': {g: int(n_by[g]) for g, _c in PLACE_GROUPS if g in n_by},
            'designs': {short[r].replace(chr(10), ' '): int((st.design_class_label == r).sum()) for r, _s, _g in rows}}


# ── the catalog README ───────────────────────────────────────────────────────
TITLES = {'evidence_grid': 'how strong is the genetic evidence for each {short} gene, and what does the gene do',
          'condition_map': 'conditions related to {short}',
          'shared_regions': 'which conditions share which {short} GWAS region',
          'gene_map': 'genes shared by {short} and its related conditions',
          'literature_genes': 'what the genetic studies of {short} report, gene by gene',
          'genetic_studies': 'every genetic study of {short}, by year, design and place'}


def readme(st, D, C, tidy, a, themes):
    """The catalog: one paragraph per figure, built from its numbers; a figure not drawn says why."""
    s = C.short
    rg, L = D['regions'], D['ladder']
    allpub, anypub, _anyot = regions_kind(rg)
    lead_only = rg.variant_id.ne('') & rg.disease_study_locus_id.eq('')
    absent = rg.variant_id.eq('')
    gs = rg.code.eq('G*').any() or L.code.eq('G*').any()
    out = {}
    eg = st.get('evidence_grid')
    if eg:
        lv = num(L.level)
        within = sorted(L[(lv == 7) & (num(L.within_disease) > 0)].symbol)
        auto = all(t.get('automatic') for t in themes)
        mixed = not auto and any(t.get('automatic') for t in themes)
        out['evidence_grid'] = (
            'Rows: the evidence ladder, strongest first. ' + gwas_levels(C, rg, L, a.curated_min, readme=True) +
            f' Levels 3–7 come from the {eg["studies"]} genetic studies of {s} in PubMed: '
            + '; '.join(f'{k} {C.LADDER[k][1][:1].lower() + C.LADDER[k][1][1:]}' for k in range(3, 8))
            + ' (the reason for each gene at level 7 is in the workbook). Every disease shows all seven levels; an '
            'empty level says why. Levels 3–6 count only findings about {s} itself: a gene associated only with an outcome within '
            f'{s}' + (f' (e.g. {", ".join(within[:2])})' if within else '') + ' stays at level 7. ' +
            ('Columns: the function of the gene, its Reactome top-level pathway (gene theme). ' if auto else
             ('Columns: the function of the gene (gene theme): the curated themes, each with its source, and for a '
              'gene they do not list its Reactome top-level pathway. ' if mixed else
              'Columns: the function of the gene (gene theme), each theme with its source. ')) +
            'A gene sits at its strongest level, in one cell; East Asian evidence is in sheet 9 of the workbook. '
            f'Shown: {eg["genes"]} genes in {eg["themes"]} columns ('
            + ', '.join(f'level {k}: {v}' for k, v in eg['by_level'].items()) + ').')
    cm = st['condition_map']
    if cm['single_region_traits']:
        adj = (C.text('region_notes', {}) or {}).get(cm['single_region_top_name'], '')
        single = (f'Traits sharing a single {s} region ({cm["single_region_traits"]} traits, '
                  f'{cm["single_region_top"]} of them the {adj + " " if adj else ""}{cm["single_region_top_name"]} '
                  'region) are shown as one gray node: they are weak evidence of a relation. ')
    else:
        single = (f'No trait shares a single {s} region. ' if len(rg) else
                  f'{s} has no GWAS region in this snapshot, so no condition shares one. ')
    if cm['shared_signal_without_region']:
        single += (f'Not shown: {cm["shared_signal_without_region"]} traits of the shared-signal category without a '
                   'shared region recorded. ')
    out['condition_map'] = (
        f'{s} at the center; each related condition sits in the sector of its category. The line to {s} shows the '
        'kind of evidence: solid, clinical evidence quoted from PubMed (line width grows with the number of items); '
        'dashed, the disease ontology only; dotted, a shared GWAS signal only. Node size grows with the amount of '
        f'evidence (clinical items, {s} GWAS regions shared, ontology). A red outline marks a condition sharing two '
        f'or more {s} GWAS regions; a small diamond, a relation whose own evidence is East Asian (clinical evidence '
        f'from a study done in East Asia or of an East Asian population, or a shared signal whose two GWAS, of {s} and '
        'of the condition, have mostly East Asian samples). ' +
        single +
        f'Shown: {cm["conditions_shown"]} conditions.')
    sr = st.get('shared_regions')
    if sr:
        out['shared_regions'] = (
            f'Rows: conditions sharing two or more of the {number(len(rg))} {"published " if allpub else ""}{s} '
            'regions, and the clinically related conditions sharing one. Columns: the regions, ' +
            ('named as published, with their lead variant, tier and status in Open Targets. ' if allpub else
             ('named as published or ' if anypub else 'named ') + 'by their top locus-to-gene gene (by the position '
             f'of the lead variant when its score is below {a.l2g_min}), with their lead variant and status in Open '
             'Targets. ') +
            f'A filled dot: a GWAS signal of the condition colocalizes with the {s} credible set. The colocalization '
            'is computed by Open Targets (coloc H4 and eCAVIAR CLPP for each pair of credible sets) and kept here at '
            f'CLPP ≥ {a.clpp_min} or H4 ≥ {a.h4_min}; color, the strongest CLPP. It is possible only where Open '
            f'Targets holds {a_(s)} credible set. ' +
            (f'An open dot, checked by this pipeline: a GWAS credible set of the condition contains the {s} lead '
             f'variant (PIP ≥ {a.clpp_min}), for regions whose lead variant Open Targets holds without a credible '
             f'set; weaker, since the {s} signal is not fine-mapped. ' if lead_only.any() else '') +
            'Dot size: number of GWAS signals, the GWAS studies of the condition with such a credible set at the '
            'region (one per study). ' +
            ('Gray columns: lead variant not in Open Targets, nothing to compare. ' if absent.any() else '') +
            f'Rows: {sr["conditions"]}; {sr["colocalisation_cells"]} filled and {sr["overlap_cells"]} open cells.' +
            ((' Not rows: 1 trait related only by a shared GWAS signal that shares one region (in the gray node of '
              'the condition map).' if sr['single_region_traits'] == 1 else
              f' Not rows: {sr["single_region_traits"]} traits related only by a shared GWAS signal that share one '
              'region each (the gray node of the condition map).') if sr.get('single_region_traits') else ''))
    else:
        out['shared_regions'] = (f'Not drawn: {s} has no GWAS region in this snapshot, so no related condition '
                                 'can share one.')
    gm = st.get('gene_map')
    if gm:
        g_core = ('none in this snapshot (no GWAS region)' if not len(rg) else
                  'the gene the publication names for a genome-wide significant region' if allpub else
                  f'the top locus-to-gene gene (score ≥ {a.l2g_min}) of a genome-wide significant region of its GWAS'
                  if not anypub else 'the gene that names a genome-wide significant region (as published, or its top '
                  'locus-to-gene gene)')
        gstar = (f'G*, {s} only, ' + (C.text('g_star', '') or 'a region significant only in a secondary analysis of '
                                       'the published GWAS') + '; ') if gs else ''
        out['gene_map'] = (
            f'Rows: all {gm["genes"]} genes studied genetically in {s}, by evidence level (the square at the left), '
            'then by the number of related conditions in which the gene has GWAS or curated evidence (the column at '
            f'the right; {gm["genes_shared"]} genes have at least one, the others are specific to {s} in this map). '
            f'Columns: {s}, then the {gm["conditions"]} related conditions with such evidence for at least one of '
            f'these genes, by category. {gm["not_shown_text"]} A cell lists the kinds of evidence: '
            f'G, for {s} {g_core}, for a related condition the top locus-to-gene gene (score ≥ {a.l2g_min}) of a '
            f'GWAS credible set; {gstar}C, curated gene–disease evidence from a source whose score is '
            f'≥ {a.curated_min}; L ({s} only), a genetic study of {s} reports an association, an enrichment or a '
            f'variant in the gene; – ({s} only), studied without a positive result (levels 6–7). Gray cells hold two '
            'or more kinds. The derivation is printed beside the map; the evidence behind each cell is in sheet 6 of '
            'the workbook.')
    lg = st.get('literature_genes')
    out['literature_genes'] = (
        f'The {lg["genes"]} genes reported by two or more genetic studies of {s}, grouped by gene theme. Each bar '
        f'counts studies by result: associated (or enriched), associated with an outcome within {s} (not counted as '
        f'evidence for {s} itself), not associated, a variant identified without a test (case reports, sequencing), '
        'other. A study that reports two results for one gene counts in both. Studies are not independent samples: '
        'some share cases.' if lg else
        f'Not drawn: no gene is reported by two or more genetic studies of {s} in this snapshot.')
    gsd = st.get('genetic_studies')
    out['genetic_studies'] = (
        f'The {gsd["studies"]} genetic studies of {s} that the curation keeps (retracted articles and conference '
        'material are left out), by year of publication (x) and design class (rows: the first class that applies, '
        'case reports last). Color: the place of the study, the country of its first author\'s institution ('
        + ', '.join(f'{k} {v}' for k, v in gsd['places'].items()) + '). Size: the number of cases genotyped, where '
        'the abstract states it. The genome-wide studies are labeled with their first author, year and cases.'
        if gsd else f'Not drawn: no genetic study of {s} is kept in this snapshot.')
    if not eg:
        out['evidence_grid'] = f'Not drawn: no gene has genetic evidence for {s} in this snapshot.'
    if not gm:
        out['gene_map'] = f'Not drawn: no gene has genetic evidence for {s} in this snapshot.'
    md = [f'# Figures of the snapshot ({tidy["fetch_date"]})', '',
          textwrap.fill(f'Generated by `scripts/draw_figures.py` from the tables of this snapshot (Open Targets '
                        f'{tidy["release"]}); do not edit. Each figure is written as PNG (300 dpi) and PDF. Every '
                        'threshold named here has its basis in sheet 11 of the workbook (Choices and basis).', 115)]
    for name in TITLES:
        md += ['', f'## {name} — {C.P.fill(TITLES[name])}', textwrap.fill(out[name], 115, break_long_words=False,
                                                                           break_on_hyphens=False)]
    return '\n'.join(md) + '\n'


def palette(S, P, cats):
    """The fills of the workbook, so that one color means one thing in every output: kinds of evidence G, G*, C, L
    (gray: several kinds) and the levels of the evidence ladder (GWAS blue, literature association orange, the weaker
    levels gray); the categories of related conditions take hues of their own (CATEGORY_COLORS)."""
    clinical = [c for c in cats if c not in (P.ontology_category, SHARED_SIGNAL)]
    cat = {c: CATEGORY_COLORS[i % len(CATEGORY_COLORS)] for i, c in enumerate(clinical)}
    cat.update({P.ontology_category: ONTOLOGY_COLOR, SHARED_SIGNAL: S.REFERENCE, 'one region': S.NEUTRAL_D})
    return {'cat': cat,
            'kind': {'G': '#C6DBEF', 'G*': '#E3EEF8', 'C': '#C7E9C0', 'L': '#FCBBA1'}, 'multi': '#8C8C8C',
            'level': {1: '#C6DBEF', 2: '#E3EEF8', 3: '#FCE4CC', 4: '#FDF2E6', 5: '#F2F2F2', 6: '#F2F2F2',
                      7: '#F2F2F2'},
            'dark': set()}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--profile', required=True, help='the disease profile (profile.py)')
    ap.add_argument('--tables', required=True)
    ap.add_argument('--plot-style', required=True, help='scripts/plot_style.py')
    ap.add_argument('--l2g-min', type=float, default=0.5)
    ap.add_argument('--curated-min', type=float, default=0.5)
    ap.add_argument('--clpp-min', type=float, default=0.01)
    ap.add_argument('--h4-min', type=float, default=0.8)
    ap.add_argument('--out-dir', default='.')
    a = ap.parse_args()
    sys.path.insert(0, str(Path(a.plot_style).resolve().parent))
    import plot_style as S
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    S.plt, S.matplotlib, S.Line2D = plt, matplotlib, Line2D
    import matplotlib.colors  # noqa: F401
    import matplotlib.patches  # noqa: F401
    import matplotlib.ticker  # noqa: F401
    import matplotlib.transforms  # noqa: F401
    S.setup_style('paper')
    plt.rcParams.update({'axes.grid': False, 'pdf.fonttype': 42})
    P = load_profile(a.profile)
    T = Path(a.tables)
    # the columns each table must have, so that a table without rows (a disease without them) reads as empty
    D = {'conditions': read_tsv(T / 'conditions.tsv', ['node_id', 'name', 'category', 'ot_ids', 'genetics',
                                                       'genetics_note', 'shared_regions', 'n_shared_regions',
                                                       'n_eas_studies']),
         'relations': read_tsv(T / 'relations.tsv', ['node_id', 'evidence_type', 'relation', 'east_asian']),
         'coloc': read_tsv(T / 'colocalisation.tsv', ['node_id', 'disease_gene', 'study_id', 'clpp']),
         'overlap': read_tsv(T / 'overlap.tsv', ['node_id', 'disease_gene', 'study_id', 'pip']),
         'regions': read_tsv(T / 'disease_regions.tsv', ['gene', 'rsid', 'chrom', 'pos_grch37', 'pos_grch38', 'tier',
                                                         'code', 'source', 'variant_id', 'disease_study_locus_id']),
         'ladder': read_tsv(T / 'disease_genes.tsv', ['symbol', 'level', 'code', 'theme', 'within_disease',
                                                      'curated_sources']),
         'gene_condition': read_tsv(T / 'gene_condition.tsv', ['symbol', 'node_id', 'evidence']),
         'findings': read_tsv(T / 'disease_findings.tsv', ['pmid', 'gene', 'symbol', 'result', 'association_with']),
         'studies': read_tsv(T / 'disease_studies.tsv', ['pmid', 'year', 'first_author', 'n_cases', 'place',
                                                         'design_class_label'])}
    themes = load_themes(T)
    LADDER, codes, evidence_codes = load_ladder(T)
    cats = list(P.category_order)
    cats += sorted(set(D['conditions'].category) - set(cats))       # categories the profile does not name: last
    C = SimpleNamespace(P=P, short=P.short, LADDER=LADDER, STATUS=load_status(T), GLOSS=load_glosses(T),
                        CODES=list(codes),
                        KIND=evidence_codes, cats=cats,
                        labels={**(P.d.get('labels') or {}), **(P.text('figures', 'labels', {}) or {})},
                        text=lambda k, d='': P.text('figures', k, d) or d,
                        gwas_pmid=str(((P.gwas or {}).get('source') or {}).get('pmid', '') or ''))
    pal = palette(S, P, cats)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    st = {'evidence_grid': evidence_grid(D, out, S, pal, C, themes, a.curated_min),
          'condition_map': condition_map(D, out, S, pal, C),
          'shared_regions': shared_regions(D, out, S, pal, C, a.clpp_min),
          'gene_map': gene_map(D, out, S, pal, C, a.l2g_min, a.curated_min),
          'literature_genes': literature_genes(D, out, S, pal, C, themes),
          'genetic_studies': genetic_studies(D, out, S, pal, C)}
    tidy = json.loads((T / 'tidy.json').read_text())
    (out / 'README.md').write_text(readme(st, D, C, tidy, a, themes))
    drawn = {k: v for k, v in st.items() if v is not None}
    stats = {**drawn, 'figures': {k: P.fill(TITLES[k]) for k in drawn},
             'not_drawn': [k for k, v in st.items() if v is None]}
    (out / 'figures.json').write_text(json.dumps(stats, indent=1, ensure_ascii=False))
    log(f'figures: {", ".join(drawn)}' + (f'; not drawn: {", ".join(stats["not_drawn"])}' if stats['not_drawn']
                                          else ''))


if __name__ == '__main__':
    main()
