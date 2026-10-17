#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The figures of one snapshot (PNG + PDF), in the shared figure style
#           of the CTEPH project (scripts/plot_style.py, a vendored copy), and one
#           catalog README for the family:
#             evidence_grid      the genes studied genetically in CTEPH:
#                                evidence level (rows, strongest first) x
#                                function (gene theme, columns)
#             condition_map      CTEPH at the center; the related conditions in
#                                sectors by category; line style = kind of
#                                evidence, node size = amount of evidence
#             shared_regions     conditions x the nine published CTEPH GWAS
#                                regions: colocalization (filled) or lead-variant
#                                overlap (open)
#             gene_map           every gene studied genetically in CTEPH, by
#                                evidence level x CTEPH and the related
#                                conditions, each cell the kinds of evidence
#                                (G, G*, C, L), with how each kind is derived
#                                printed beside it
#             literature_genes   the genetic studies of CTEPH: studies per gene
#                                and result, for genes reported by two or more
#                                studies, grouped by gene theme
#             genetic_studies    every genetic study of CTEPH by year and design
#                                class, colored by place (country of the first
#                                author's institution), sized by cases
#           Presentation only: every mark is computed from the tables.
# Inputs  : the tidy tables, the gene themes, plot_style.py
# Outputs : <name>.png, <name>.pdf, README.md, figures.json
# ---------------------------------------------------------------------------
import argparse
import json
import math
import re
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import (LADDER, NO_GENE, OTHER_GENES, assign_themes, gene_tallies, load_themes, log,  # noqa: E402
                       read_tsv)

CORE = 'EFO_0022954'
CATEGORIES = ['Precursor / spectrum', 'Risk factor', 'Comorbidity', 'Intermediate trait', 'PH group',
              'Shared genetic signal']
SHORT = {  # labels of the conditions whose names are too long for a figure
    'Chronic thromboembolic pulmonary disease without PH (CTEPD)': 'CTEPD',
    'Antiphospholipid syndrome / antiphospholipid antibodies': 'Antiphospholipid syndrome',
    'Ventriculoatrial shunt, infected pacemaker and intravascular devices': 'Shunts, pacemakers, devices',
    'Hypothyroidism / thyroid replacement therapy': 'Hypothyroidism',
    'Malignancy (history of cancer)': 'Malignancy',
    'Sickle cell disease and hemoglobinopathies': 'Sickle cell disease',
    'Fibrinogen (level and variants)': 'Fibrinogen',
    'Factor XI level (F11 locus)': 'Factor XI level',
    'Platelet count and activation': 'Platelets',
    'Pulmonary arterial hypertension associated with another disease': 'Associated PAH',
    'drug- or toxin-induced pulmonary arterial hypertension': 'Drug-induced PAH',
    'heritable pulmonary arterial hypertension': 'Heritable PAH',
    'idiopathic pulmonary arterial hypertension': 'Idiopathic PAH',
    'Pulmonary arterial hypertension': 'PAH',
    'hyperuricemia-pulmonary hypertension-renal failure-alkalosis syndrome': 'HUPRA syndrome',
    'pulmonary veno-occlusive disease and/or pulmonary capillary hemangiomatosis': 'PVOD / PCH',
    'pulmonary hypertension, neonatal': 'Neonatal PH',
    'pulmonary hypertension': 'Pulmonary hypertension',
    'plasminogen activator inhibitor 1 measurement': 'PAI-1 level',
    'tissue plasminogen activator amount': 'tPA level',
    'amount of von Willebrand factor (human) in blood': 'VWF in blood',
}
RESULTS = [('associated', 'Associated or enriched'), ('within_cteph', 'Associated with an outcome within CTEPH'),
           ('not_associated', 'Not associated'),
           ('identified', 'Variant identified'), ('other', 'Other result')]
CODES = ['G', 'G*', 'C', 'L']                                          # order of the codes in a cell
KIND = {'GWAS': 'G', 'GWAS (PE-adjusted)': 'G*', 'Mendelian': 'C', 'Literature': 'L'}
SUP = str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')


def label(name):
    """A short, sentence-case label: the SHORT map, else Open Targets' trait name tidied."""
    if name in SHORT:
        return SHORT[name]
    s = re.sub(r'^amount of |^level of | in blood$| measurement$|\(human\) ', '', name)
    s = re.sub(r'\s+', ' ', s).strip()
    return s[0].upper() + s[1:]


def p_unicode(p):
    m, e = f'{float(p):.1e}'.split('e')
    return f'{m} × 10{str(int(e)).translate(SUP)}'


def save(fig, out, name, S):
    for ext in ('png', 'pdf'):
        fig.savefig(out / f'{name}.{ext}', dpi=300 if ext == 'png' else None, bbox_inches='tight', pad_inches=0.08)
    S.plt.close(fig)


def by_condition(d, cond, col):
    """GWAS signals of each condition at each CTEPH region (Open Targets IDs mapped to the condition):
    number of studies and the best value of `col`."""
    alias = {i: c.node_id for c in cond.itertuples() for i in c.ot_ids.split()}
    d = d.assign(node_id=d.node_id.map(lambda i: alias.get(i, i)), v=pd.to_numeric(d[col]))
    return d.groupby(['node_id', 'cteph_gene']).agg(n=('study_id', 'nunique'), best=('v', 'max'))


# ── evidence_grid ────────────────────────────────────────────────────────────
def evidence_grid(T, out, S, P, themes):
    """Evidence level (rows) x function (columns): each gene studied genetically in CTEPH, at its strongest
    level, in the column of its gene theme."""
    plt = S.plt
    L = T['ladder'].copy()
    L['level'] = L.level.astype(int)
    short = {t['name']: t['short'] for t in themes}
    short[OTHER_GENES] = 'Other genes'
    cols = [t['name'] for t in themes if (L.theme == t['name']).any()] + \
        ([OTHER_GENES] if (L.theme == OTHER_GENES).any() else [])
    levels = sorted(LADDER)
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
    rh = {lv: lines[lv] * LH + 0.2 for lv in levels}
    n_st = int(T['studies'].pmid.nunique())
    foot = (f'Each gene sits at its strongest level. Levels 1–2: the regions of the CTEPH GWAS as published '
            f'(Liley 2024, Table 1); G, meta-analysis P < 5 × 10⁻⁸; G*, genome-wide significant only after leveraging '
            f'pulmonary embolism (V < 5 × 10⁻⁸). Levels 3–7: what the {n_st} genetic studies of CTEPH in PubMed report '
            f'for the gene, each finding checked against its abstract; an association in two or more studies ranks '
            f'above one, and levels 3–6 count only findings about CTEPH itself (a gene associated only with an outcome '
            f'within CTEPH stays at level 7). Columns: the function of the gene; each theme has its source (workbook '
            f'sheet 11). '
            f'East Asian evidence: workbook sheet 9.')
    fl = textwrap.wrap(foot, 215)
    W = LW + CW * len(cols) + NW
    H = HH + sum(rh.values()) + 0.25 + 0.2 * len(fl)
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis('off')
    fill = P['level']
    for j, c in enumerate(cols):                                       # column headers
        x = LW + j * CW + CW / 2
        n = int((L.theme == c).sum())
        ax.text(x, HH - 0.3, '\n'.join(textwrap.wrap(short[c], 17, break_long_words=False)), ha='center',
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
        code, text = LADDER[lv]
        ax.add_patch(S.matplotlib.patches.FancyBboxPatch((0.95, y + h / 2 - 0.14), 0.36, 0.28,
                                                         boxstyle='round,pad=0,rounding_size=0.05', facecolor=fill[lv],
                                                         edgecolor=S.NEUTRAL_D, lw=0.6))
        ax.text(1.13, y + h / 2, str(lv), ha='center', va='center', fontsize=11, fontweight='bold',
                color='white' if lv in P['dark'] else S.INK)
        lab = (f'{code} · ' if code else '') + text
        ax.text(1.42, y + h / 2, '\n'.join(textwrap.wrap(lab, 31)), ha='left', va='center', fontsize=10.5,
                color=S.INK, linespacing=1.1)
        for j, c in enumerate(cols):
            genes = cell[(lv, c)]
            x0 = LW + j * CW
            if genes:
                ax.add_patch(S.matplotlib.patches.Rectangle((x0 + 0.04, y + 0.05), CW - 0.08, h - 0.1,
                                                            color=fill[lv], lw=0, zorder=1))
            for g, gx, line in place[(lv, c)]:
                ax.text(x0 + 0.14 + gx, y + 0.1 + line * LH + LH / 2, tag(g), ha='left', va='center', fontsize=FS,
                        style='italic', color='white' if lv in P['dark'] else S.INK, zorder=2)
        n = int((L.level == lv).sum())
        ax.text(LW + len(cols) * CW + NW / 2, y + h / 2, str(n), ha='center', va='center', fontsize=11,
                color=S.INK)
        y += h
        ax.plot([LW, LW + len(cols) * CW], [y, y], color=S.NEUTRAL, lw=0.6, zorder=0)
    for j in range(1, len(cols)):
        ax.plot([LW + j * CW] * 2, [HH, y], color=S.NEUTRAL, lw=0.6, zorder=0)
    ax.plot([0.05, W - 0.05], [y, y], color=S.INK, lw=1.2)
    for (a, b), text in (((1, 2), 'Published GWAS'), ((3, 7), 'Genetic studies of CTEPH')):   # the two sources
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
def condition_map(T, out, S, P):
    plt = S.plt
    cond, rel = T['conditions'], T['relations']
    nl = cond.set_index('node_id').n_shared_regions.astype(int)
    keep = cond[(cond.category != 'Shared genetic signal') | (cond.node_id.map(nl) >= 2)].copy()
    single = cond[(cond.category == 'Shared genetic signal') & (cond.node_id.map(nl) < 2)]
    single_abo = int((single.shared_regions == 'ABO').sum())
    ev = rel.groupby('node_id').evidence_type.apply(set).to_dict()
    ncl = rel[rel.evidence_type == 'clinical'].groupby('node_id').size().to_dict()
    eas = set(rel[rel.east_asian == 'True'].node_id)
    keep['_o'] = keep.category.map(CATEGORIES.index)
    keep['label'] = keep.name.map(label)
    keep = keep.sort_values(['_o', 'label'], key=lambda s: s.str.lower() if s.dtype == object else s)
    groups = [(c, keep[keep.category == c]) for c in CATEGORIES if (keep.category == c).any()]
    n = len(keep) + 1                                                  # + the summary node of single-region traits
    gap = 2.0                                                          # empty slots between sectors
    step = 2 * math.pi / (n + gap * (len(groups) + 1))
    fig, ax = plt.subplots(figsize=(12.5, 12.5))
    ax.set_aspect('equal')
    ax.axis('off')
    R, slot, pos = 1.0, -math.pi / 2 - 0.5 * step, {}
    sectors = []
    for cat, g in groups + [('one region', None)]:
        slot += gap * step
        start = slot
        rows = list(g.itertuples()) if g is not None else [None]
        for x in rows:
            slot += step
            pos[x.node_id if x is not None else '__single'] = slot
        sectors.append((cat, start + 0.5 * step, slot + 0.5 * step, len(rows)))
    for cat, a0, a1, k in sectors:                                     # a faint wedge behind each sector
        col = P['cat'].get(cat, S.NEUTRAL_D)
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
        ax.scatter(R * math.cos(a), R * math.sin(a), s=size, color=P['cat'][x.category], zorder=3,
                   edgecolors=S.ACCENT if multi else 'white', linewidths=2.6 if multi else 1.0)
        if x.node_id in eas:
            ax.scatter(1.095 * math.cos(a), 1.095 * math.sin(a), s=30, marker='D', color=S.INK, zorder=4,
                       linewidths=0)
        deg = math.degrees(a)
        flip = 90 < (deg % 360) < 270
        ax.text(1.17 * math.cos(a), 1.17 * math.sin(a), x.label, rotation=deg + 180 if flip else deg,
                rotation_mode='anchor', ha='right' if flip else 'left', va='center', fontsize=10.5, color=S.INK)
    a = pos['__single']                                                # the single-region traits, as one node
    ax.plot([0, R * math.cos(a)], [0, R * math.sin(a)], ls=':', lw=1.0, color=S.INK_SOFT, alpha=0.55, zorder=1)
    ax.scatter(R * math.cos(a), R * math.sin(a), s=260, color=S.NEUTRAL, zorder=3, edgecolors='white')
    ax.text(R * math.cos(a), R * math.sin(a), str(len(single)), ha='center', va='center', fontsize=9.5,
            color=S.INK, zorder=4)
    deg = math.degrees(a)
    ax.text(1.17 * math.cos(a), 1.17 * math.sin(a), 'One CTEPH region only', rotation=deg + 180,
            rotation_mode='anchor', ha='right', va='center', fontsize=10.5, color=S.INK_SOFT, style='italic')
    ax.scatter(0, 0, s=5200, color=S.DATA_DARK, zorder=5)
    ax.text(0, 0, 'CTEPH', ha='center', va='center', color='white', fontsize=17, fontweight='bold', zorder=6)
    L = S.Line2D
    handles = [L([0], [0], color=S.INK_SOFT, lw=2.2, ls='-', label='Clinical evidence (width: items)'),
               L([0], [0], color=S.INK_SOFT, lw=1.2, ls='--', label='Disease ontology only'),
               L([0], [0], color=S.INK_SOFT, lw=1.2, ls=':', label='Shared GWAS signal only'),
               L([0], [0], marker='o', ls='', markersize=11, markerfacecolor='white', markeredgecolor=S.ACCENT,
                 markeredgewidth=2.4, label='Shares ≥ 2 CTEPH GWAS regions'),
               L([0], [0], marker='D', ls='', markersize=6, color=S.INK, label='East Asian clinical evidence'),
               L([0], [0], marker='o', ls='', markersize=8, color=S.NEUTRAL_D, label='Node size: amount of evidence')]
    handles.append(S.matplotlib.patches.Patch(color=S.NEUTRAL, label=f'{len(single)} traits sharing one CTEPH '
                                              f'region ({single_abo} of them ABO): weak'))
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.075), ncol=3, frameon=False,
               fontsize=10.5, handlelength=2.6, columnspacing=1.4)
    names = {'Shared genetic signal': 'Shared GWAS signal only (≥ 2 CTEPH regions)'}
    cats = [S.matplotlib.patches.Patch(color=P['cat'][c], label=f'{names.get(c, c)} ({k})')
            for c, _a0, _a1, k in sectors if c in P['cat'] and c != 'one region']
    fig.legend(handles=cats, loc='lower center', bbox_to_anchor=(0.5, 0.94), ncol=3, frameon=False, fontsize=11,
               columnspacing=1.6)
    lim = 1.85
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.93, bottom=0.08)
    save(fig, out, 'condition_map', S)
    return {'conditions_shown': len(keep), 'single_region_traits': len(single), 'single_region_abo': single_abo,
            'east_asian': len(eas & set(keep.node_id))}


# ── shared_regions ───────────────────────────────────────────────────────────
def shared_regions(T, out, S, P, clpp_min):
    plt = S.plt
    cond, rg = T['conditions'].copy(), T['regions']
    cond['nreg'] = cond.n_shared_regions.astype(int)
    rows = cond[(cond.nreg >= 2) | ((cond.category != 'Shared genetic signal') & (cond.nreg >= 1))].copy()
    rows['_o'] = rows.category.map(CATEGORIES.index)
    rows['label'] = rows.name.map(label)
    rows = rows.sort_values(['_o', 'nreg', 'label'], ascending=[True, False, True])
    ca, oa = by_condition(T['coloc'], cond, 'clpp'), by_condition(T['overlap'], cond, 'pip')
    regs = list(rg.itertuples())
    xs = [j + (0.5 if r.evidence == 'pe' else 0) for j, r in enumerate(regs)]  # a gap between G and G* regions
    h = 0.34 * len(rows) + 1.6
    fig, ax = plt.subplots(figsize=(16.5, h))
    fig.subplots_adjust(left=0.15, right=0.8, top=1 - 1.25 / h, bottom=0.2 / h)
    S.despine(ax, grid_axis='none')
    ys = {nid: i for i, nid in enumerate(rows.node_id)}
    cmap = S.matplotlib.colors.LinearSegmentedColormap.from_list('clpp', ['#C9DCEC', S.DATA, S.DATA_DARK])
    norm = S.matplotlib.colors.LogNorm(vmin=0.01, vmax=1.0)
    size = lambda n: 30 + 38 * math.sqrt(n)  # noqa: E731
    for x, r in zip(xs, regs):
        if r.cteph_study_locus_id:
            ax.axvline(x, color=S.NEUTRAL, lw=0.8, zorder=0)
        elif r.variant_id:
            ax.axvline(x, color=S.NEUTRAL, lw=0.8, ls=(0, (2, 2)), zorder=0)
        else:                                                          # nothing to compare with in Open Targets
            ax.add_patch(S.matplotlib.patches.Rectangle((x - 0.42, -0.5), 0.84, len(rows), color=S.NEUTRAL,
                                                        alpha=0.35, lw=0, zorder=0))
            ax.text(x, (len(rows) - 1) / 2, 'lead variant not in Open Targets', rotation=90, ha='center',
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
    status = lambda r: ('credible set' if r.cteph_study_locus_id else 'lead only' if r.variant_id  # noqa: E731
                        else 'not in OT')
    blend = S.matplotlib.transforms.blended_transform_factory(ax.transData, ax.transAxes)
    up = lambda pt: S.matplotlib.transforms.offset_copy(blend, fig=fig, y=pt, units='points')  # noqa: E731
    for x, r in zip(xs, regs):                                         # column headers, drawn so that the region
        ax.text(x, 1.0, f'{r.rsid}\ntier {r.tier}\n{status(r)}', transform=up(5), ha='center', va='bottom',
                fontsize=9.5, color=S.INK_SOFT, linespacing=1.2)     # name can be bold italic
        ax.text(x, 1.0, r.gene, transform=up(44), ha='center', va='bottom', fontsize=11.5, fontweight='bold',
                style='italic', color=S.INK)
    ax.set_xticks([])
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows.label, fontsize=11)
    ax.set_ylim(len(rows) - 0.5, -0.6)
    ax.set_xlim(-0.6, xs[-1] + 0.5)
    ax.tick_params(length=0)
    for meta, text in ((True, 'G · CTEPH meta-analysis P < 5 × 10⁻⁸'), (False, 'G* · only after leveraging PE')):
        sel = [x for x, r in zip(xs, regs) if (r.evidence == 'meta') == meta]
        if not sel:
            continue
        ax.plot([min(sel) - 0.4, max(sel) + 0.4], [1.0, 1.0], transform=up(64), color=S.INK, lw=1.0, clip_on=False)
        ax.text((min(sel) + max(sel)) / 2, 1.0, text, transform=up(67), ha='center', va='bottom', fontsize=11,
                fontweight='bold', color=S.INK)
    prev = None                                                        # category bands
    for i, x in enumerate(rows.itertuples()):
        if x.category != prev:
            k = (rows.category == x.category).sum()
            ax.add_patch(S.matplotlib.patches.Rectangle((-0.6, i - 0.5), 0.12, k, color=P['cat'][x.category],
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
               label='Colocalizes with the\nCTEPH credible set'),
             L([0], [0], marker='o', ls='', markersize=10, markerfacecolor='white', markeredgecolor=S.DATA_DARK,
               markeredgewidth=1.4, label=f'Its credible set holds the\nCTEPH lead variant\n(PIP ≥ {clpp_min}; weaker)')]
    b0, t0 = 0.2 / h, 1 - 1.25 / h                                     # legends at figure level (a tight bounding
    fy = lambda ya: b0 + ya * (t0 - b0)  # noqa: E731                    # box drops axes-level extra legends)
    fig.legend(handles=kinds, loc='upper left', bbox_to_anchor=(0.815, fy(0.68)), frameon=False, fontsize=10.5,
               labelspacing=1.1, borderpad=0.2, title='Shared region', title_fontsize=11, alignment='left')
    sizes = [L([0], [0], marker='o', ls='', markersize=math.sqrt(size(k)), color=S.NEUTRAL_D, label=f'{k}')
             for k in (1, 5, 20)]
    fig.legend(handles=sizes, title='GWAS signals', loc='upper left', bbox_to_anchor=(0.815, fy(0.43)),
               frameon=False, fontsize=10.5, title_fontsize=11, labelspacing=1.3, borderpad=0.2, alignment='left')
    handles = [S.matplotlib.patches.Patch(color=P['cat'][c], label=c) for c in CATEGORIES if (rows.category == c).any()]
    fig.legend(handles=handles, title='Category (bar at left)', loc='upper left', bbox_to_anchor=(0.815, fy(0.17)),
               frameon=False, fontsize=10.5, title_fontsize=11, alignment='left', borderpad=0.2)
    save(fig, out, 'shared_regions', S)
    return {'conditions': len(rows), 'regions': [r.gene for r in regs],
            'colocalisation_cells': int(sum((n, r.gene) in ca.index for n in rows.node_id for r in regs)),
            'overlap_cells': int(sum((n, r.gene) not in ca.index and (n, r.gene) in oa.index
                                     for n in rows.node_id for r in regs))}


# ── gene_map ─────────────────────────────────────────────────────────────────
def gene_map(T, out, S, P, l2g_min, curated_min):
    """Every gene studied genetically in CTEPH (the evidence ladder), by level, x CTEPH and the related conditions
    with GWAS or curated evidence for any of them; each cell lists the kinds of evidence (G, G*, C, L), with how
    each kind is derived beside the map. Genes shared with no related condition stay in the map (an empty row)."""
    plt = S.plt
    gc, cond, L = T['gene_condition'], T['conditions'], T['ladder'].copy()
    L['level'] = L.level.astype(int)
    gc = gc.assign(code=gc.evidence.map(KIND))
    rel = gc[(gc.node_id != CORE) & gc.symbol.isin(set(L.symbol))]
    L['nrel'] = L.symbol.map(rel.groupby('symbol').node_id.nunique()).fillna(0).astype(int)
    L = L.sort_values(['level', 'nrel', 'symbol'], ascending=[True, False, True])
    genes, level, nrel = list(L.symbol), dict(zip(L.symbol, L.level)), dict(zip(L.symbol, L.nrel))
    cat = dict(zip(cond.node_id, cond.category))
    names = dict(zip(cond.node_id, cond.name))
    cols = (rel.groupby('node_id').symbol.nunique().rename('n').reset_index()
            .assign(_o=lambda d: d.node_id.map(lambda n: CATEGORIES.index(cat[n]))))
    cols = cols.sort_values(['_o', 'n'], ascending=[True, False])
    nodes = [CORE] + list(cols.node_id)
    cell = gc[gc.symbol.isin(genes) & gc.node_id.isin(nodes)].groupby(['symbol', 'node_id']).code.apply(
        lambda x: ''.join(sorted(set(x), key=CODES.index))).to_dict()
    with_genes = set(gc[gc.node_id != CORE].node_id)
    n_none = int((cond.genetics == 'none').sum())
    n_coloc = int((cond.genetics == 'coloc_only').sum())
    n_other = len(with_genes - set(cols.node_id))
    ng, nn = len(genes), len(nodes)
    H = 0.3 * ng + 3.6
    fig = plt.figure(figsize=(17.5, H))
    ax = fig.add_axes([0.11, 0.04, 0.52, 0.78])
    ax.set_xlim(-1.5, nn + 0.6)
    ax.set_ylim(ng - 0.5, -0.5)
    for i, g in enumerate(genes):
        lv = level[g]
        ax.add_patch(S.matplotlib.patches.Rectangle((-1.45, i - 0.42), 0.9, 0.84, facecolor=P['level'][lv],
                                                    edgecolor=S.NEUTRAL_D, lw=0.4))
        ax.text(-1.0, i, str(lv), ha='center', va='center', fontsize=8.5, color=S.INK)
        for j, n in enumerate(nodes):
            k = cell.get((g, n), '')
            if not k:
                if j == 0:                                             # levels 6-7: studied, no positive result
                    ax.text(j, i, '–', ha='center', va='center', fontsize=9, color=S.INK_SOFT)
                continue
            fc = P['kind'].get(k, P['multi'])
            ax.add_patch(S.matplotlib.patches.Rectangle((j - 0.45, i - 0.42), 0.9, 0.84, color=fc, lw=0))
            ax.text(j, i, k, ha='center', va='center', fontsize=8.5 if len(k) < 3 else 7.5,
                    color=S.INK if k in P['kind'] else 'white', fontweight='bold')
        ax.text(nn + 0.15, i, str(nrel[g]), ha='center', va='center', fontsize=9,
                color=S.INK if nrel[g] else S.INK_SOFT)
    ax.set_xticks([-1] + list(range(nn)) + [nn + 0.15])
    ax.set_xticklabels(['Evidence level'] + ['CTEPH'] + [label(names[n]) for n in cols.node_id] +
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
        c = cat.get(n) if n else None
        if c != prev and prev is not None:
            ax.add_patch(S.matplotlib.patches.Rectangle((j0 - 0.5, ng - 0.5 + 0.25), j - j0, 0.35,
                                                        color=P['cat'][prev], clip_on=False, lw=0))
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
    lines = [('How each cell is derived', 'h'),
             ('G   GWAS. CTEPH: the gene the publication names for a genome-wide significant region (meta-analysis '
              'P < 5 × 10⁻⁸; Liley 2024). Related condition: the top locus-to-gene gene (score ≥ {l2g}) of at least '
              'one of its GWAS credible sets in Open Targets.'.format(l2g=l2g_min), 'G'),
             ('G*  GWAS, CTEPH only: a region genome-wide significant only after leveraging pulmonary embolism '
              '(V < 5 × 10⁻⁸); marked apart because the signal is partly that of PE.', 'G*'),
             ('C   Curated: curated gene–disease evidence for the condition or its descendants, from a source whose '
              'Open Targets score for the gene is ≥ {c} (ClinGen, PanelApp, Orphanet, ClinVar, UniProt, '
              'Gene2Phenotype, gene burden).'.format(c=curated_min), 'C'),
             ('L   Literature, CTEPH only: a genetic study of CTEPH in PubMed reports an association, an enrichment '
              'or a variant in the gene; the sentence is checked against the abstract.', 'L'),
             ('Gray: two or more kinds (letters shown). –: studied in CTEPH without a positive result '
              '(levels 6–7).', ''),
             ('', ''),
             ('Rows and columns', 'h'),
             (f'Rows: all {ng} genes studied genetically in CTEPH, by evidence level (square at the left: 1 G, 2 G*, '
              '3 associated in two or more studies, 4 in one, 5 variant identified without a test, 6 tested, not '
              'associated, 7 other result), then by the number of related conditions in which the gene has G or C '
              f'evidence (right). {sum(1 for g in genes if not nrel[g])} genes have no such condition: they are '
              'specific to CTEPH in this map.', ''),
             (f'Columns: CTEPH, then the {nn - 1} related conditions with G or C evidence for at least one of these '
              f'genes, by category (bar below) and number of genes. Not shown: {n_coloc} traits related by a shared '
              f'GWAS signal only, whose genes are not fetched by design (their shared regions are in the regions '
              f'figure); {n_none} conditions without an Open Targets record; {n_other} conditions with gene evidence '
              'for none of these genes.', ''),
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
        if k in CODES:
            w = textwrap.fill(t[4:], 54)
            tx.add_patch(S.matplotlib.patches.Rectangle((0, y - 0.26 * inch), 0.045, 0.26 * inch, color=P['kind'][k],
                                                        transform=tx.transAxes, clip_on=False))
            tx.text(0.0225, y - 0.13 * inch, k, fontsize=10, color=S.INK, fontweight='bold', ha='center',
                    va='center')
            tx.text(0.065, y, w, fontsize=11, va='top', color=S.INK)
        else:
            tx.text(0, y, w, fontsize=11, va='top', color=S.INK_SOFT)
        y -= 0.2 * inch * (w.count('\n') + 1) + 0.18 * inch
    handles = [S.matplotlib.patches.Patch(color=P['cat'][c], label=c) for c in CATEGORIES
               if c in {cat[n] for n in cols.node_id}]
    tx.legend(handles=handles, loc='upper left', bbox_to_anchor=(0, y - 0.2 * inch), frameon=False, fontsize=10.5,
              title='Category of the related condition', title_fontsize=11, alignment='left')
    save(fig, out, 'gene_map', S)
    return {'genes': ng, 'genes_shared': sum(1 for g in genes if nrel[g]), 'conditions': nn - 1,
            'not_shown': {'shared_signal_only': n_coloc, 'no_record': n_none, 'no_cteph_gene': n_other},
            'cells': {k: int(v) for k, v in pd.Series(list(cell.values())).value_counts().items()}}


# ── literature_genes ─────────────────────────────────────────────────────────
def literature_genes(T, out, S, P, themes):
    plt = S.plt
    fx = assign_themes(T['findings'], themes)
    gt = gene_tallies(fx)
    gt = gt[(gt.n_studies >= 2) & (gt.theme != NO_GENE)].copy()
    order = {t['name']: i for i, t in enumerate(themes)}
    short = {t['name']: t['short'] for t in themes}
    gt['_t'] = gt.theme.map(lambda t: order.get(t, len(order)))
    gt = gt.sort_values(['_t', 'n_studies', 'label'], ascending=[True, False, True]).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(10.5, 0.42 * len(gt) + 1.8))
    S.despine(ax, grid_axis='x')
    cols = {'associated': S.ACCENT, 'within_cteph': S.ACCENT_LT, 'not_associated': S.NEUTRAL_D, 'identified': S.DATA,
            'other': S.NEUTRAL}
    left = np.zeros(len(gt))
    for key, lab in RESULTS:
        v = gt[key].to_numpy(dtype=float)
        ax.barh(range(len(gt)), v, left=left, color=cols[key], label=lab, height=0.68, edgecolor='white', lw=0.6)
        left += v
    ax.set_yticks(range(len(gt)))
    ax.set_yticklabels(gt.label, fontsize=11.5, fontstyle='italic')
    ax.invert_yaxis()
    ax.set_xlabel('Genetic studies of CTEPH (a study can report several results)', fontsize=11.5)
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


# design classes, top to bottom (config/curation.yaml, design_classes), and the place groups of a study
DESIGN_ROWS = ['Genome-wide studies (GWAS, genome sequencing, Mendelian randomization)',
               'Clonal hematopoiesis and somatic variants', 'Exome and gene-panel sequencing', 'HLA typing',
               'Candidate-gene association', 'Other designs', 'Case reports']
PLACE_GROUPS = [('Japan', ['Japan']), ('China', ['China']), ('Korea, Taiwan', ['Korea', 'Taiwan', 'Hong Kong']),
                ('Europe', ['United Kingdom', 'France', 'Germany', 'Austria', 'Switzerland', 'Czech Republic',
                            'Hungary', 'Italy', 'Spain', 'Netherlands', 'Belgium', 'Denmark', 'Sweden', 'Norway',
                            'Finland', 'Poland', 'Ireland', 'Portugal', 'Greece']),
                ('North America', ['United States', 'Canada']), ('Other', None), ('Not known', None)]


def genetic_studies(T, out, S, P):
    """Every genetic study of CTEPH: year (x) x design class (rows); color, the place; size, the cases."""
    plt = S.plt
    st = T['studies'].copy()
    st['year_n'] = pd.to_numeric(st.year)
    st['cases'] = pd.to_numeric(st.n_cases, errors='coerce')

    def group(c):
        if not c:
            return 'Not known'
        return next((g for g, cs in PLACE_GROUPS if cs and c in cs), 'Other')
    st['group'] = st.place.map(group)
    color = {'Japan': S.ACCENT, 'China': '#B08A3C', 'Korea, Taiwan': '#D9A03C', 'Europe': S.DATA,
              'North America': '#3F8A6E', 'Other': '#7A5C9E', 'Not known': S.NEUTRAL_D}
    rows = [r for r in DESIGN_ROWS if (st.design_class_label == r).any()]
    short = {DESIGN_ROWS[0]: 'Genome-wide (GWAS, genome\nsequencing, MR)', DESIGN_ROWS[1]: 'Clonal hematopoiesis',
             DESIGN_ROWS[2]: 'Exome or gene panel', DESIGN_ROWS[3]: 'HLA typing', DESIGN_ROWS[4]: 'Candidate gene',
             DESIGN_ROWS[5]: 'Other designs', DESIGN_ROWS[6]: 'Case reports'}

    def size(n):
        return 30 + 9 * math.sqrt(n) if n == n else 30
    fig, ax = plt.subplots(figsize=(12.5, 0.85 * len(rows) + 2.6))
    S.despine(ax, grid_axis='x')
    lo, hi = int(st.year_n.min()), int(st.year_n.max())
    pts = []
    for i, r in enumerate(rows):
        sub = st[st.design_class_label == r].sort_values(['year_n', 'pmid'])
        for _y, g in sub.groupby('year_n'):
            offs = (np.arange(len(g)) - (len(g) - 1) / 2) * 0.22
            for x, o in zip(g.itertuples(), offs):
                ax.scatter(x.year_n, i + o, s=size(x.cases), color=color[x.group], alpha=0.85, edgecolor='white',
                           lw=0.6, zorder=3)
                if r == DESIGN_ROWS[0]:
                    pts.append((x, i + o))
    # the genome-wide studies cluster in a few recent years: their labels stand in a column to the left of the
    # cluster, in the order of the points' heights, so that the leader lines do not cross
    pts.sort(key=lambda q: (q[1], q[0].year_n))
    x0 = min(x.year_n for x, _y in pts) - 2.5 if pts else hi
    for k, (x, y) in enumerate(pts):
        lab = f'{x.first_author.split()[0]} {x.year}' + (f', {int(x.cases):,} cases' if x.cases == x.cases else '')
        ty = -0.62 + k * (1.05 / max(len(pts) - 1, 1))
        ax.annotate(lab, (x.year_n, y), xytext=(x0, ty), textcoords='data', ha='right', va='center', fontsize=9.5,
                    color=S.INK_SOFT, arrowprops=dict(arrowstyle='-', color=S.NEUTRAL_D, lw=0.7, shrinkA=3,
                                                      shrinkB=4))
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f'{short[r]} ({int((st.design_class_label == r).sum())})' for r in rows], fontsize=11)
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
            'designs': {short[r].replace(chr(10), ' '): int((st.design_class_label == r).sum()) for r in rows}}


README = """# Figures of the snapshot ({date})

Generated by `scripts/draw_figures.py` from the tables of this snapshot (Open Targets {release}); do not edit.
Each figure is written as PNG (300 dpi) and PDF. Every threshold named here has its basis in sheet 11 of the
workbook (Choices and basis).

## evidence_grid — how strong is the genetic evidence for each CTEPH gene, and what does the gene do
Rows: the evidence ladder, strongest first. Levels 1–2 come from the published CTEPH GWAS (Liley 2024,
Table 1): G, meta-analysis P < 5 × 10⁻⁸; G*, genome-wide significant only after leveraging pulmonary embolism
(V < 5 × 10⁻⁸). Levels 3–7 come from the {studies} genetic studies of CTEPH in PubMed: associated in two or more
studies; in one; a variant identified without an association test; tested, not associated; other result only.
Levels 3–6 count only findings about CTEPH itself: a gene associated only with an outcome within CTEPH (e.g.
GNB3) stays at level 7. Columns: the function of the gene (gene theme), each theme with its source. A gene sits at its strongest
level, in one cell; East Asian evidence is in sheet 9 of the workbook. Shown: {grid_genes} genes in {grid_themes}
columns ({by_level}).

## condition_map — conditions related to CTEPH
CTEPH at the center; each related condition sits in the sector of its category. The line to CTEPH shows the
kind of evidence: solid, clinical evidence quoted from PubMed (line width grows with the number of items);
dashed, the disease ontology only; dotted, a shared GWAS signal only. Node size grows with the amount of
evidence (clinical items, CTEPH GWAS regions shared, ontology). A red outline marks a condition sharing two or
more CTEPH GWAS regions; a small diamond, clinical evidence from an East Asian population. Traits sharing a
single CTEPH region ({single} traits, {single_abo} of them the pleiotropic ABO region) are shown as one gray node:
they are weak evidence of a relation. Shown: {shown} conditions.

## shared_regions — which conditions share which CTEPH GWAS region
Rows: conditions sharing two or more of the nine published CTEPH regions, and the clinically related conditions
sharing one. Columns: the regions, named as published, with their lead variant, tier and status in Open Targets.
A filled dot: a GWAS signal of the condition colocalizes with the CTEPH credible set. The colocalization is
computed by Open Targets (coloc H4 and eCAVIAR CLPP for each pair of credible sets) and kept here at CLPP ≥ {clpp}
or H4 ≥ 0.8; color, the strongest CLPP. It is possible only where Open Targets holds a CTEPH credible set. An open
dot, checked by this pipeline: a GWAS credible set of the condition contains the CTEPH lead variant
(PIP ≥ {clpp}), for regions whose lead variant Open Targets holds without a credible set; weaker, since the CTEPH
signal is not fine-mapped. Dot size: number of GWAS signals, the GWAS studies of the condition with such a
credible set at the region (one per study). Gray columns: lead variant not in Open Targets, nothing to compare. Rows: {reg_rows};
{reg_coloc} filled and {reg_overlap} open cells.

## gene_map — genes shared by CTEPH and its related conditions
Rows: all {genes} genes studied genetically in CTEPH, by evidence level (the square at the left), then by the
number of related conditions in which the gene has GWAS or curated evidence (the column at the right; {shared} genes
have at least one, the others are specific to CTEPH in this map). Columns: CTEPH, then the {conds} related
conditions with such evidence for at least one of these genes, by category; the traits related by a shared GWAS
signal only are not columns, since their genes are not fetched (their shared regions are in shared_regions). A
cell lists the kinds of evidence: G, for CTEPH the gene the publication names for a genome-wide significant
region, for a related condition the top locus-to-gene gene (score ≥ {l2g}) of a GWAS credible set; G*, CTEPH only,
a region significant only after leveraging PE; C, curated gene–disease evidence from a source whose score is
≥ {cur}; L (CTEPH only), a genetic study of CTEPH reports an association, an enrichment or a variant in the gene;
– (CTEPH only), studied without a positive result (levels 6–7). Gray cells hold two or more kinds. The derivation
is printed beside the map; the evidence behind each cell is in sheet 6 of the workbook.

## literature_genes — what the genetic studies of CTEPH report, gene by gene
The {lit_genes} genes reported by two or more genetic studies of CTEPH, grouped by gene theme. Each bar counts
studies by result: associated (or enriched), associated with an outcome within CTEPH (not counted as evidence for
CTEPH itself), not associated, a variant identified without a test (case reports, sequencing), other. A study
that reports two results for one gene counts in both. Studies are not independent samples: some share cases.

## genetic_studies — every genetic study of CTEPH, by year, design and place
The {gs_n} genetic studies of CTEPH that the curation keeps (retracted articles and conference material are left
out), by year of publication (x) and design class (rows: the first class that applies, case reports last). Color:
the place of the study, the country of its first author's institution ({gs_places}). Size: the number of cases
genotyped, where the abstract states it. The genome-wide studies are labeled with their first author, year and
cases.
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tables', required=True)
    ap.add_argument('--themes', required=True, help='gene themes of the CTEPH literature (YAML)')
    ap.add_argument('--plot-style', required=True, help='scripts/plot_style.py')
    ap.add_argument('--l2g-min', type=float, default=0.5)
    ap.add_argument('--curated-min', type=float, default=0.5)
    ap.add_argument('--clpp-min', type=float, default=0.01)
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
    # the fills of the workbook, so that one color means one thing in every output: kinds of evidence G, G*, C, L
    # (gray: several kinds) and the levels of the evidence ladder (GWAS blue, literature association orange, the
    # weaker levels gray); categories of related conditions take hues of their own
    P = {'cat': {'Precursor / spectrum': '#C2627A', 'Risk factor': '#B08A3C', 'Comorbidity': '#7A5C9E',
                 'Intermediate trait': '#2E8B8B', 'PH group': '#4D5B6A', 'Shared genetic signal': S.REFERENCE,
                 'one region': S.NEUTRAL_D},
         'kind': {'G': '#C6DBEF', 'G*': '#E3EEF8', 'C': '#C7E9C0', 'L': '#FCBBA1'}, 'multi': '#8C8C8C',
         'level': {1: '#C6DBEF', 2: '#E3EEF8', 3: '#FCE4CC', 4: '#FDF2E6', 5: '#F2F2F2', 6: '#F2F2F2', 7: '#F2F2F2'},
         'dark': set()}
    T = Path(a.tables)
    D = {'conditions': read_tsv(T / 'conditions.tsv'), 'relations': read_tsv(T / 'relations.tsv'),
         'coloc': read_tsv(T / 'colocalisation.tsv'), 'overlap': read_tsv(T / 'overlap.tsv'),
         'regions': read_tsv(T / 'cteph_regions.tsv'), 'ladder': read_tsv(T / 'cteph_genes.tsv'),
         'gene_condition': read_tsv(T / 'gene_condition.tsv'), 'findings': read_tsv(T / 'cteph_findings.tsv'),
         'studies': read_tsv(T / 'cteph_studies.tsv')}
    themes = load_themes(a.themes)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    st = {'evidence_grid': evidence_grid(D, out, S, P, themes), 'condition_map': condition_map(D, out, S, P),
          'shared_regions': shared_regions(D, out, S, P, a.clpp_min),
          'gene_map': gene_map(D, out, S, P, a.l2g_min, a.curated_min),
          'literature_genes': literature_genes(D, out, S, P, themes),
          'genetic_studies': genetic_studies(D, out, S, P)}
    tidy = json.loads((T / 'tidy.json').read_text())
    eg, cm, sr = st['evidence_grid'], st['condition_map'], st['shared_regions']
    (out / 'README.md').write_text(README.format(
        date=tidy['fetch_date'], release=tidy['release'], studies=eg['studies'], grid_genes=eg['genes'],
        grid_themes=eg['themes'], by_level=', '.join(f'level {k}: {v}' for k, v in eg['by_level'].items()),
        single=cm['single_region_traits'], single_abo=cm['single_region_abo'], shown=cm['conditions_shown'],
        clpp=a.clpp_min, reg_rows=sr['conditions'], reg_coloc=sr['colocalisation_cells'],
        reg_overlap=sr['overlap_cells'], genes=st['gene_map']['genes'], conds=st['gene_map']['conditions'],
        shared=st['gene_map']['genes_shared'],
        l2g=a.l2g_min, cur=a.curated_min, lit_genes=st['literature_genes']['genes'],
        gs_n=st['genetic_studies']['studies'],
        gs_places=', '.join(f'{k} {v}' for k, v in st['genetic_studies']['places'].items())))
    (out / 'figures.json').write_text(json.dumps(st, indent=1))
    log(f'figures: {", ".join(st)}')


if __name__ == '__main__':
    main()
