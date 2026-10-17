#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The snapshot workbook of the CTEPH knowledge graph, organized by
#           question. Each question sheet opens with its answer in brief, then
#           a nested table: group > record > source, opened with the outline
#           buttons (+ / -, levels 1 2 3). The appendix holds the flat tables,
#           with filters. House style of the CTEPH project's integration tables
#           booktabs rules, italic gene symbols, P as
#           m x 10^-e, lettered footnotes, the reason for every "—", one
#           Abbreviations line. Presentation only: every value is taken as is
#           from the tables; generated summaries are labeled.
#             Guide                    the questions by part, how to read, the
#                                      snapshot and its sources
#             A  CTEPH genetics
#               1 CTEPH genes at a glance evidence level x function grid
#               2 CTEPH gene evidence  level > gene > finding (source sentence)
#               3 Genetic studies      design > study, its methods by step
#             B  Related conditions
#               4 Related conditions   category > condition > evidence
#               5 Shared GWAS regions  condition x the nine CTEPH regions
#               6 Shared genes         gene x condition map of genetic evidence
#               7 Condition genetics   category > condition > gene
#             C  Context
#               8 Mechanisms           Reactome theme > target > publication
#               9 East Asian           kind of evidence > item
#             D  Provenance
#               10 Curation            what was left out or corrected, and why
#               11 Choices and basis   every choice, its value and basis
#               12 Changes             since the previous snapshot
#             A1-A6                    GWAS studies, GWAS loci, shared signals,
#                                      curated genes, CTEPH studies, screening
#           The workbook names no file, script, parameter or pipeline.
# Inputs  : the tidy tables, the curation outputs, the graph schema, summaries,
#           changes, knowledge JSON, gene themes
# Outputs : the workbook, workbook.json (layout and rows per sheet)
# ---------------------------------------------------------------------------
import argparse
import json
import math
import re
import sys
import warnings
import zipfile
from pathlib import Path

import pandas as pd
import xlsxwriter

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import (CURATED, LADDER, MISSING, OTHER_GENES, ancestry_label, assign_themes,  # noqa: E402
                       gene_tallies, is_east_asian, is_japanese, load_choices, load_themes, log, num, pub_label,
                       pub_url,
                       read_tsv, ref_short, study_population, tally_text, variant_label)

FORBIDDEN_TOKENS = ['.tsv', '.py', '.nf', '.yaml', '/Users', '/home', 'nextflow', 'europepmc', 'clinical_precedence',
                    'gwas_credible_sets', '/snapshots', 'FetchDate', 'ensembl_id', 'l2g1', 'datasource', '_mainland',
                    'node_id', 'custom:']
FONT = 'Arial'
INK, INK2, GREY, RULE_GREY, NOTE = '#1A1A1A', '#595959', '#7F7F7F', '#BFBFBF', '#404040'
LINK = '#1F4E79'
GROUP0, GROUP1 = '#DCE6F1', '#F2F2F2'        # rows of outline level 0 / rows that open
KIND = {'GWAS': 'G', 'GWAS (PE-adjusted)': 'G*', 'Mendelian': 'C', 'Literature': 'L'}
CODE_ORDER = ['G', 'G*', 'C', 'L']
KIND_FILL = {'G': '#C6DBEF', 'G*': '#E3EEF8', 'C': '#C7E9C0', 'L': '#FCBBA1'}
MULTI_FILL = '#8C8C8C'
SOURCE_FILL = {'full text': '#C7E9C0', 'abstract': '#EDEDED'}
THIN, MEDIUM = 1, 2
CORE = 'EFO_0022954'

# The question sheets, in four parts: A, CTEPH's own genetics (its genes by strength of evidence and function,
# what each study reports, the studies themselves); B, the conditions related to it and what they share with it
# (GWAS regions, genes), and their own genetics; C, context (mechanisms, East Asian evidence); D, provenance (what
# the curation left out or corrected, the basis of every choice, the changes since the previous snapshot).
SHEETS = [('gene_map', '1 CTEPH genes at a glance',
           'Which genes have genetic evidence for CTEPH, how strong is it, and what do they do?'),
          ('genes', '2 CTEPH gene evidence', 'What does each study report for each CTEPH gene?'),
          ('methods', '3 Genetic studies', 'Which genetic studies of CTEPH are there, and how were they done?'),
          ('conditions', '4 Related conditions', 'Which conditions are related to CTEPH, and on what evidence?'),
          ('regions', '5 Shared GWAS regions', 'Which CTEPH GWAS regions do the related conditions share?'),
          ('shared', '6 Shared genes', 'Which genes do CTEPH and its related conditions share?'),
          ('cond_gen', '7 Condition genetics', 'What are the genetics of each related condition?'),
          ('mech', '8 Mechanisms', 'Which mechanisms does the CTEPH literature describe?'),
          ('eas', '9 East Asian', 'What evidence comes from East Asian and Japanese populations?'),
          ('curation', '10 Curation', 'What was left out or corrected, and why?'),
          ('choices', '11 Choices and basis', 'On what basis was each item and value chosen?'),
          ('changes', '12 Changes', 'What changed since the previous snapshot?')]
SN = {k: n for k, n, _q in SHEETS}
# the four parts of the line of argument, and the tab color of their sheets
PARTS = [('A  CTEPH genetics', ['gene_map', 'genes', 'methods'], '#1F4E79'),
         ('B  Related conditions', ['conditions', 'regions', 'shared', 'cond_gen'], '#2E7D7D'),
         ('C  Context', ['mech', 'eas'], '#6A4C93'),
         ('D  Provenance', ['curation', 'choices', 'changes'], '#595959')]
PART_OF = {k: (name, color) for name, keys, color in PARTS for k in keys}
QUESTIONS = [(n, q) for _k, n, q in SHEETS]


def QT(key):
    """A question sheet's title: its number and question."""
    n, q = next((n, q) for k, n, q in SHEETS if k == key)
    return f'{n.split()[0]}  {q}'


APPENDIX = [('A1 GWAS studies', 'Every GWAS study of CTEPH and the related conditions: design, size, ancestry'),
            ('A2 GWAS loci', 'Every credible set: lead variant, P, β, fine-mapping, locus-to-gene genes'),
            ('A3 Shared signals', 'CTEPH GWAS regions × GWAS signals of other traits: colocalization and lead-variant '
             'overlap'),
            ('A4 Curated genes', 'Curated gene evidence of the related conditions, by source'),
            ('A5 CTEPH studies', 'The genetic studies of CTEPH, one row each'),
            ('A6 Screening', 'Every PubMed record screened for the CTEPH literature, and the decision')]
SH = ['Guide'] + [s for s, _q in QUESTIONS] + [s for s, _d in APPENDIX]
CATEGORY_ORDER = ['Precursor / spectrum', 'Risk factor', 'Comorbidity', 'Intermediate trait', 'PH group']
CATEGORY_DESC = {
    'Precursor / spectrum': 'Venous thromboembolism and the chronic thromboembolic disease spectrum CTEPH develops '
                            'from.',
    'Risk factor': 'Conditions and procedures that raise the risk of CTEPH (mostly against non-thromboembolic PH '
                   'controls).',
    'Comorbidity': 'Conditions more frequent in CTEPH, or reported with it.',
    'Intermediate trait': 'Blood traits between genotype and thrombosis.',
    'PH group': 'Pulmonary hypertension and its other forms: the neighbors of CTEPH in the disease ontology.',
    'multi': 'Traits without a clinical or ontology relation that share two or more CTEPH GWAS regions '
             '(colocalization or lead-variant overlap).',
    'single': 'Traits that share one CTEPH GWAS region only, most often the pleiotropic ABO region: weak evidence '
              'of a relation.'}
SHORT = {'chronic thromboembolic pulmonary hypertension': 'CTEPH', 'Venous thromboembolism': 'VTE',
         'Deep vein thrombosis': 'DVT', 'Pulmonary embolism': 'PE', 'Pulmonary arterial hypertension': 'PAH',
         'Pulmonary arterial hypertension associated with another disease': 'PAH, associated',
         'heritable pulmonary arterial hypertension': 'Heritable PAH',
         'idiopathic pulmonary arterial hypertension': 'Idiopathic PAH',
         'pulmonary veno-occlusive disease and/or pulmonary capillary hemangiomatosis': 'PVOD / PCH',
         'hyperuricemia-pulmonary hypertension-renal failure-alkalosis syndrome': 'HUPRA syndrome',
         'Hypothyroidism / thyroid replacement therapy': 'Hypothyroidism',
         'Sickle cell disease and hemoglobinopathies': 'Sickle cell disease',
         'Fibrinogen (level and variants)': 'Fibrinogen', 'Factor XI level (F11 locus)': 'Factor XI level',
         'Antiphospholipid syndrome / antiphospholipid antibodies': 'Antiphospholipid syndrome',
         'von Willebrand factor level': 'VWF level', 'pulmonary hypertension': 'Pulmonary hypertension'}
METHOD_STEPS = [('Samples and diagnosis', ['cohort / recruitment', 'diagnosis of CTEPH', 'controls', 'DNA source']),
                ('Genotyping and sequencing', ['genotyping platform', 'sequencing', 'HLA typing', 'imputation',
                                               'variant calling / quality control']),
                ('Analysis', ['statistical test', 'covariates', 'multiple testing', 'power', 'software']),
                ('Replication', ['replication'])]
DESIGN_GROUPS = ['Genome-wide studies (GWAS, genome sequencing, Mendelian randomization)',
                 'Exome and gene-panel sequencing', 'Clonal hematopoiesis and somatic variants', 'HLA typing',
                 'Candidate-gene association', 'Other designs', 'Case reports']
ABBR = {'CTEPH': 'chronic thromboembolic pulmonary hypertension', 'GWAS': 'genome-wide association study',
        'L2G': 'locus-to-gene', 'PMID': 'PubMed identifier', 'CLPP': 'colocalization posterior probability',
        'PICS': 'probabilistic identification of causal SNPs', 'EAS': 'East Asian', 'PE': 'pulmonary embolism',
        'VTE': 'venous thromboembolism', 'DVT': 'deep vein thrombosis', 'PH': 'pulmonary hypertension',
        'PAH': 'pulmonary arterial hypertension', 'β': 'effect of the variant allele as reported by the study',
        'PMC': 'PubMed Central', 'PPR': 'Europe PMC preprint record', 'HLA': 'human leukocyte antigen',
        'ClinGen': 'Clinical Genome Resource', 'PanelApp': 'Genomics England PanelApp',
        'PVOD': 'pulmonary veno-occlusive disease', 'PCH': 'pulmonary capillary hemangiomatosis',
        'HUPRA': 'hyperuricemia, pulmonary hypertension, renal failure and alkalosis', 'VWF': 'von Willebrand factor',
        'CH': 'clonal hematopoiesis',
        'PIP': 'posterior inclusion probability', 'OR': 'odds ratio', 'MAF': 'minor allele frequency',
        'MPN': 'myeloproliferative neoplasm', 'OT': 'Open Targets Platform', 'rsID': 'dbSNP reference SNP identifier'}
# what a finding is about (curated), as shown in the narrow column
ABOUT_SHORT = {'CTEPH': 'CTEPH', 'a phenotype or outcome within CTEPH': 'outcome within CTEPH',
               'another trait, with CTEPH as the exposure': 'another trait',
               'a functional experiment': 'functional', 'a background statement, not a finding of the study':
               'background'}


# ── text helpers ─────────────────────────────────────────────────────────────
def p_parts(mant, exp):
    """P from mantissa / exponent: 7.00 × 10⁻²⁷ (no float underflow)."""
    m, e = num(pd.Series([mant]))[0], num(pd.Series([exp]))[0]
    if pd.isna(m) or pd.isna(e) or m <= 0:
        return None
    m, e = float(m), int(e)
    while m >= 9.995:
        m, e = m / 10, e + 1
    while m < 1:
        m, e = m * 10, e - 1
    if e >= -2:
        return [f'{m * 10 ** e:.2f}']
    return [f'{m:.2f} × 10', ('sup', f'−{-e}')]


def genes_parts(pairs):
    """[(symbol, score), ...] -> italic symbols with scores."""
    out = []
    for k, (g, s) in enumerate(pairs):
        if k:
            out.append(', ')
        out += [('i', g)] + ([f' {s:.2f}'] if s is not None and not pd.isna(s) else [])
    return out


def italic_list(genes, sep=', '):
    out = []
    for k, g in enumerate(genes):
        out += ([sep] if k else []) + [('i', g)]
    return out


def md_parts(text):
    """'**Lead.** text' -> rich parts with the bold runs."""
    out = []
    for k, seg in enumerate(re.split(r'\*\*(.+?)\*\*', text)):
        if seg:
            out.append(('b', seg) if k % 2 else seg)
    return out


def text_of(parts):
    return ''.join(p if isinstance(p, str) else p[1] for p in parts)


def text_len(parts):
    return len(text_of(parts))


def blank(x):
    return x is None or (isinstance(x, float) and math.isnan(x)) or (isinstance(x, str) and x == '')


def plain(v, sup=False):
    """Text of a cell read back with rich_text=True, without its superscript runs (sup=True: only those)."""
    from openpyxl.cell.rich_text import CellRichText
    if isinstance(v, CellRichText):
        return ''.join((b if isinstance(b, str) else b.text) for b in v
                       if (not isinstance(b, str) and b.font is not None and b.font.vertAlign == 'superscript') == sup)
    return '' if (v is None or sup) else str(v)


def wrap_parts(parts, width):
    """Split rich parts (str / (style, str)) into lines of <= width characters, at spaces and newlines."""
    tokens = []
    for p in parts:
        st, tx = ('', p) if isinstance(p, str) else p
        for li, line in enumerate(tx.split('\n')):
            if li:
                tokens.append(('', '\n'))
            words = line.split(' ')
            for k, wd in enumerate(words):
                tokens.append((st, wd + (' ' if k < len(words) - 1 else '')))
    lines, cur, n = [], [], 0
    for st, tx in tokens:
        if tx == '\n':
            lines.append(cur)
            cur, n = [], 0
            continue
        if tx == '':
            continue
        if n + len(tx.rstrip()) > width and cur and st not in ('sup', 'sub') and cur[-1][1].endswith(' '):
            lines.append(cur)
            cur, n = [], 0
        cur.append((st, tx))
        n += len(tx)
    if cur:
        lines.append(cur)
    out = []
    for ln in lines:
        if ln:
            ln[-1] = (ln[-1][0], ln[-1][1].rstrip())
        out.append([tx if st == '' else (st, tx) for st, tx in ln])
    return out


def year_span(a, b):
    a, b = num(pd.Series([a]))[0], num(pd.Series([b]))[0]
    if pd.isna(a):
        return None
    return f'{int(a)}' if a == b else f'{int(a)}–{int(b)}'


def fmt2(x):
    v = num(pd.Series([x]))[0]
    return None if pd.isna(v) else f'{v:.2f}'


def intx(x):
    v = num(pd.Series([x]))[0]
    return None if pd.isna(v) else f'{int(v):,}'


def count(x):
    """A count as a number cell (None when not a number)."""
    v = num(pd.Series([x]))[0]
    return None if pd.isna(v) else int(v)


def short(name):
    s = SHORT.get(name, name)
    return s if len(s) <= 30 else s[:29] + '…'


def b(x):
    return str(x) == 'True'


def cite(r):
    return ' '.join(y for y in (str(r.first_author).split()[0] if r.first_author else '', str(r.year)) if y)


def lines_needed(text, width):
    """Wrapped lines of a cell (Arial 9 pt: about 1.3 characters per unit of column width, less what wrapping at
    spaces loses)."""
    return sum(max(1, math.ceil(len(seg) / max(width * 1.3 - 4, 1))) for seg in text.split('\n'))


def plural(n, word, words=None):
    return f'{n:,} {word if n == 1 else (words or word + "s")}'


# ── the writer ───────────────────────────────────────────────────────────────
class Sheet:
    """An XlsxWriter worksheet with a format cache and a rich-text helper."""

    def __init__(self, wb, ws):
        self.wb, self.ws, self._fmt = wb, ws, {}

    def fmt(self, **kw):
        props = {'font_name': FONT, 'font_size': 9, 'font_color': INK, 'valign': 'top'}
        props.update({k: v for k, v in kw.items() if v is not None})
        key = tuple(sorted(props.items()))
        if key not in self._fmt:
            self._fmt[key] = self.wb.add_format(props)
        return self._fmt[key]

    def font(self, style, b, i, color, sz):
        props = {'font_name': FONT, 'font_size': sz, 'font_color': color, 'bold': bool(b), 'italic': bool(i)}
        if style == 'i':
            props['italic'] = True
        elif style == 'b':
            props['bold'] = True
        elif style == 'bi':
            props['bold'], props['italic'] = True, True
        elif style == 'sup':
            props['font_script'] = 1
        elif style == 'g':
            props['font_color'] = INK2
        elif style == 'gi':
            props['font_color'], props['italic'] = INK2, True
        key = ('font',) + tuple(sorted(props.items()))
        if key not in self._fmt:
            self._fmt[key] = self.wb.add_format(props)
        return self._fmt[key]

    def rich(self, r, c, parts, cell, b=False, i=False, color=INK, sz=9):
        frags = []
        for p in parts:
            st, tx = ('', p) if isinstance(p, str) else p
            if tx == '':
                continue
            if frags and frags[-1][0] == st:
                frags[-1][1] += tx
            else:
                frags.append([st, tx])
        if not frags:
            self.ws.write_blank(r, c, None, self.fmt(**cell))
        elif len(frags) == 1:
            st, tx = frags[0]
            f = dict(cell, font_size=sz, font_color=INK2 if st in ('g', 'gi') else color,
                     bold=True if (b or st in ('b', 'bi')) else None,
                     italic=True if (i or st in ('i', 'bi', 'gi')) else None, font_script=1 if st == 'sup' else None)
            self.ws.write_string(r, c, tx, self.fmt(**f))
        else:
            args = []
            for st, tx in frags:
                args += [self.font(st, b, i, color, sz), tx]
            self.ws.write_rich_string(r, c, *args, self.fmt(**cell))


def new_sheet(wb, name, tab=None):
    sh = Sheet(wb, wb.add_worksheet(name))
    sh.ws.hide_gridlines(2)
    sh.ws.set_zoom(90)
    if tab:
        sh.ws.set_tab_color(tab)
    return sh


def title(sh, text, sub):
    sh.rich(0, 0, [('b', text)], {'align': 'left', 'valign': 'vcenter'}, sz=11)
    sh.rich(1, 0, sub if isinstance(sub, list) else [sub], {'align': 'left'}, i=True, color=INK2)
    sh.ws.set_row(0, 18)


def header_row(sh, top, cols):
    rot = any(c.get('rot') for c in cols)
    for j, c in enumerate(cols):
        sh.ws.set_column(j, j, c['w'])
        h = c['h'] if isinstance(c['h'], list) else [('b', c['h'])]
        cell = {'align': c.get('a', 'left'), 'valign': 'bottom', 'text_wrap': not c.get('rot'), 'bold': True,
                'top': MEDIUM, 'top_color': INK, 'bottom': THIN, 'bottom_color': INK}
        if c.get('rot'):
            cell.update({'rotation': 90, 'align': 'center'})
        sh.rich(top, j, h + [('sup', c.get('n') or '')], cell, b=True)
        if c.get('tip'):
            sh.ws.write_comment(top, j, c['tip'], {'font_name': FONT, 'font_size': 9, 'x_scale': 2})
    if rot:
        sh.ws.set_row(top, min(150, max(5.3 * max(text_len(c['h'] if isinstance(c['h'], list) else [c['h']])
                                                  for c in cols if c.get('rot')) + 12, 40)))
    else:
        hl = max(math.ceil(text_len(c['h'] if isinstance(c['h'], list) else [c['h']]) / max(c['w'] * 1.05, 1))
                 for c in cols)
        sh.ws.set_row(top, max(26, 12.5 * hl + 2))


def table(sh, top, cols, rows, links=None, freeze=1, autofilter=False):
    """A flat booktabs table. cols: [{'h': header, 'w': width, 'n': footnote letter, 'a': align, 'wrap': bool}];
    rows: lists of cells (a str, rich parts, None for '—', '' for blank); links: {(i, j): url}.
    Returns the next free row."""
    links = links or {}
    header_row(sh, top, cols)
    for i, row in enumerate(rows):
        r = top + 1 + i
        last = i == len(rows) - 1
        lines = 1
        for j, (c, v) in enumerate(zip(cols, row)):
            cell = {'align': c.get('a', 'left'), 'text_wrap': bool(c.get('wrap'))}
            if last:
                cell.update({'bottom': MEDIUM, 'bottom_color': INK})
            parts = [] if isinstance(v, str) and v == '' else [MISSING] if blank(v) else \
                (v if isinstance(v, list) else [str(v)])
            if (i, j) in links and not blank(v):
                sh.ws.write_url(r, j, links[(i, j)], sh.fmt(**cell, font_color=LINK, underline=1),
                                string=text_of(parts))
            else:
                sh.rich(r, j, parts, cell)
            if c.get('wrap'):
                lines = max(lines, lines_needed(text_of(parts), c['w']))
        if lines > 1:
            sh.ws.set_row(r, 12.0 * lines + 1)
    if not rows:
        sh.rich(top + 1, 0, ['No records.'], {'align': 'left', 'bottom': MEDIUM, 'bottom_color': INK})
    if freeze is not None:
        sh.ws.freeze_panes(top + 1, freeze)
    if autofilter and rows:
        sh.ws.autofilter(top, 0, top + len(rows), len(cols) - 1)
    return top + 1 + max(len(rows), 1)


def notes(sh, row, items, width=150):
    """items: [(letter, parts)]; letter '' for unlettered lines. A heading repeated right after its own block (two
    "Missing values" headings, one from the sheet and one from missing_block) is written once."""
    seen = set()
    for letter, parts in items:
        if parts and parts[0] == ('b', 'Missing values') and len(parts) == 1:
            if 'missing' in seen:
                continue
            seen.add('missing')
        for k, ln in enumerate(wrap_parts(parts, width - 3)):
            lead = [('sup', letter), ' '] if (letter and k == 0) else []
            sh.rich(row, 0, lead + ln, {'align': 'left', 'indent': 1 if k else None}, color=NOTE)
            row += 1
    return row


def missing_block(items):
    out = [('', [('b', 'Missing values')])]
    for where, text in items:
        out.append(('', ['   ', ('b', MISSING), f', in {where}: ', text]))
    return out


def abbrev(keys):
    return ('', ['Abbreviations: ' + '; '.join(f'{k}, {ABBR[k]}' for k in sorted(keys, key=str.lower)) + '.'])


def print_setup(ws, last_row, ncol, footer, header=3):
    ws.set_landscape()
    ws.set_paper(8)
    ws.fit_to_pages(1, 0)
    ws.repeat_rows(header, header)
    ws.print_area(0, 0, last_row, ncol - 1)
    ws.set_margins(left=0.4, right=0.4)
    ws.set_footer(f'&L{footer}&RPage &P of &N')


# ── question sheets: a brief, then a nested table ────────────────────────────
def R(lvl, kind, cells, span=None, links=None, tips=None, fills=None, hide=False, key=None, bold=None):
    """One row of a nested table. kind: group0 | group1 (a row that opens) | item | quote.
    cells: one per column (str, rich parts, number, None for '—', '' for blank); span: (j0, j1) merges those
    columns into one cell holding cells[j0]; links: {j: url or ('@', sheet, key)}; tips: {j: comment};
    fills: {j: (background, font color)}; key: an anchor other sheets can link to; bold: a bold label (default:
    group0 rows)."""
    return {'lvl': lvl, 'kind': kind, 'cells': cells, 'span': span, 'links': links or {}, 'tips': tips or {},
            'fills': fills or {}, 'hide': hide, 'key': key, 'bold': kind == 'group0' if bold is None else bold}


def brief_lines(spec):
    """The lines above the table: [(parts, style, outline level, hidden)]. A summary marked closed sits under its
    heading as a closed outline group."""
    width = min(int(sum(c['w'] for c in spec['cols']) * 1.05), 175)
    out = []
    if spec.get('brief'):
        out.append(([('b', 'In brief')], 'head', 0, False))
        for item in spec['brief']:
            for k, ln in enumerate(wrap_parts(item, width - 3)):
                out.append(((['• '] if k == 0 else ['   ']) + ln, 'text', 0, False))
    closed = bool(spec.get('summary_closed'))
    for head, paras in spec.get('summaries') or []:
        out.append(([('b', head + ('; open it with + at the left' if closed and paras else ''))], 'shead', 0, False))
        for p in paras:
            for ln in wrap_parts(p, width):
                out.append((ln, 'text', 1, closed))
            out.append(([], 'text', 1, closed))
    return out + [([], 'text', 0, False)] if out else out


def layout(spec):
    """Row of the header and of every body row (0-based)."""
    top = 3 + len(brief_lines(spec)) + (1 if spec.get('spanners') else 0)
    return top, {r['key']: top + 1 + i for i, r in enumerate(spec['rows']) if r['key'] is not None}


STYLE = {'group0': {'bg_color': GROUP0}, 'group1': {'bg_color': GROUP1}, 'item': {}, 'quote': {}}


def write_spec(wb, spec, anchors, foot):
    sh = new_sheet(wb, spec['name'], spec.get('tab'))
    ws, cols = sh.ws, spec['cols']
    ws.outline_settings(True, False, True, False)
    title(sh, spec['title'], spec['sub'])
    ws.write_url(0, len(cols) - 1, "internal:'Guide'!A1", sh.fmt(font_color=LINK, underline=1, align='right'),
                 string='← Guide')
    r = 3
    lines = brief_lines(spec)
    for k, (parts, style, lvl, hide) in enumerate(lines):
        sh.rich(r, 0, parts, {'align': 'left'}, color=INK if style != 'shead' else LINK)
        nxt = lines[k + 1] if k + 1 < len(lines) else None
        if lvl or style == 'shead':
            ws.set_row(r, None, None, {'level': lvl, 'hidden': hide,
                                       'collapsed': bool(style == 'shead' and nxt and nxt[3])})
        r += 1
    top, _ = layout(spec)
    if spec.get('spanners'):
        for j0, j1, text in spec['spanners']:
            f = sh.fmt(bold=True, align='center', valign='bottom', text_wrap=True, bottom=THIN, bottom_color=INK)
            if j1 > j0:
                ws.merge_range(top - 1, j0, top - 1, j1, text, f)
            else:
                ws.write_string(top - 1, j0, text, f)
        ws.set_row(top - 1, 24)
    header_row(sh, top, cols)
    if spec.get('header_height'):
        ws.set_row(top, spec['header_height'])
    rows = spec['rows']
    for i, row in enumerate(rows):
        rr = top + 1 + i
        last = i == len(rows) - 1
        style = STYLE[row['kind']]
        height = 1
        j = 0
        while j < len(cols):
            j1 = row['span'][1] if row['span'] and row['span'][0] == j else j
            c, v = cols[j], row['cells'][j]
            cell = {'align': c.get('a', 'left') if j1 == j else 'left', 'text_wrap': bool(c.get('wrap') or j1 > j),
                    **style}
            if j == 0:
                cell['indent'] = row['lvl']
            if last:
                cell.update({'bottom': MEDIUM, 'bottom_color': INK})
            color = INK2 if row['kind'] == 'quote' else INK
            bold = j == 0 and row['bold']
            if j in row['fills']:
                cell['bg_color'], color = row['fills'][j]
                bold = bold or color != INK
            if isinstance(v, (int, float)) and not isinstance(v, bool):
                ws.write_number(rr, j, v, sh.fmt(**cell, num_format='#,##0', font_color=color))
                text = f'{v:,}'
            else:
                parts = [] if isinstance(v, str) and v == '' else [MISSING] if blank(v) else \
                    (v if isinstance(v, list) else [str(v)])
                text = text_of(parts)
                if j1 > j:
                    ws.merge_range(rr, j, rr, j1, '', sh.fmt(**cell))
                link = row['links'].get(j)
                if link and not blank(v):
                    url = link if isinstance(link, str) else \
                        f"internal:'{link[1]}'!A{anchors[(link[1], link[2])] + 1}"
                    ital = all(not isinstance(x, str) and x[0] in ('i', 'bi') for x in parts)   # a gene symbol
                    ws.write_url(rr, j, url, sh.fmt(**cell, font_color=LINK, underline=1, bold=bold or None,
                                                    italic=ital or None), string=text)
                else:
                    sh.rich(rr, j, parts, cell, b=bold, color=color)
            if j in row['tips']:
                tip = row['tips'][j]
                ws.write_comment(rr, j, tip, {'font_name': FONT, 'font_size': 9, 'x_scale': 3.2,
                                              'y_scale': max(1.0, lines_needed(tip, 60) / 4.5)})
            if cell['text_wrap']:
                height = max(height, lines_needed(text, sum(cols[k]['w'] for k in range(j, j1 + 1))))
            j = j1 + 1
        nxt = rows[i + 1] if i + 1 < len(rows) else None
        opts = {'level': row['lvl'], 'hidden': row['hide'],
                'collapsed': bool(nxt and nxt['lvl'] > row['lvl'] and nxt['hide'])}
        ws.set_row(rr, 12.0 * height + 1 if height > 1 else None, None, opts)
    end = top + 1 + max(len(rows), 1)
    if not rows:
        sh.rich(top + 1, 0, ['No records.'], {'align': 'left', 'bottom': MEDIUM, 'bottom_color': INK})
    if spec.get('freeze'):
        ws.freeze_panes(top + 1, spec['freeze'])
    end = notes(sh, end + 1, spec.get('notes') or [], spec.get('note_width', 150))
    print_setup(ws, end, len(cols), foot, header=top)
    return {'header': top, 'rows': len(rows), 'spanner': bool(spec.get('spanners'))}


def spec(name, title_, sub, cols, rows, **kw):
    key = next((k for k, n, _q in SHEETS if n == name), None)
    tab = PART_OF[key][1] if key in PART_OF else '#1F4E79'
    return {'name': name, 'title': title_, 'sub': sub, 'cols': cols, 'rows': rows, 'tab': tab, **kw}


def summary_block(D, kind, head):
    s = D['summ'][D['summ'].kind == kind]
    if not len(s):
        return []
    x = s.iloc[0]
    if not b(x.grounded) and x.problems.startswith('not generated'):
        return [(f'{head}: not generated in this build (its source records changed with the curation; no language '
                 'model was called)', [])]
    if not b(x.grounded):
        return [(f'{head}: withheld, it failed its check', [])]
    return [(f'{head} (generated by a language model from the records below; it cites only their PMIDs)',
             [md_parts(p) for p in x.summary.split('\n') if p.strip()])]


# ── data shared by several sheets ────────────────────────────────────────────
def prepare(D):
    cond, gc = D['conditions'], D['gene_condition'].copy()
    gc['code'] = gc.evidence.map(KIND)
    assert gc.code.notna().all(), sorted(set(gc[gc.code.isna()].evidence))
    D['gc'] = gc
    D['cteph_genes'] = set(gc[gc.node_id == CORE].symbol)
    D['cteph_kinds'] = gc[gc.node_id == CORE].groupby('symbol').code.apply(codes).to_dict()
    D['rel_n'] = gc[gc.node_id != CORE].groupby('symbol').node_id.nunique().to_dict()
    cst = D['cst'].copy()
    cst['population'] = [study_population(r.country, r.ancestry) for r in cst.itertuples()]
    cst['place'] = cst.place.replace('', MISSING)
    cst['jp'] = [is_japanese(r) for r in cst.itertuples()]
    cst['eas'] = [is_east_asian(r) for r in cst.itertuples()]
    cst['cite'] = [f'{r.first_author} {r.year}' for r in cst.itertuples()]
    D['cst'] = cst
    fx = assign_themes(D['cfi'], D['themes']).merge(
        cst[['pmid', 'year', 'cite', 'design', 'design_class_label', 'population', 'place', 'jp', 'eas', 'n_cases',
             'n_controls']], on='pmid', how='left')
    D['fx'] = fx.assign(_y=num(fx.year)).sort_values(['_y', 'pmid'])
    D['gt'] = gene_tallies(D['fx'])
    other = sorted(set(D['gt'][D['gt'].theme == OTHER_GENES].label))
    if other:
        log(f'genes in no theme (listed under "{OTHER_GENES}"): {", ".join(other)}')
    D['in_sheet5'] = set(cond[cond.genetics.eq('full') & cond.ot_ids.ne('')].node_id)
    D['short'] = {t['name']: t.get('short', t['name']) for t in D['themes']} | {OTHER_GENES: OTHER_GENES}
    D['choice'] = {it['key']: it for g in D['choices'] for it in g['items']}


def codes(x):
    """'GWAS', 'Literature', ... codes of one gene x condition, in a fixed order: 'GL', 'G*L', 'GC' ..."""
    return ''.join(sorted(set(x), key=CODE_ORDER.index))


def basis(D, key):
    """'Basis: Mountjoy 2021 (sheet 10).' — the short form of a choice's basis, for footnotes."""
    it = D['choice'][key]
    who = '; '.join(r[0] for r in it['refs']) or {'design': 'a design decision', 'data check': 'checked on the data',
                                                  'user': 'decided with the user'}[it['basis_type']]
    return f'Basis: {who} (sheet {SN["choices"].split()[0]}).'


def rel_category(D, node):
    c = D['conditions'].set_index('node_id')
    cat = c.category[node]
    if cat != 'Shared genetic signal':
        return cat
    return 'multi' if int(c.n_shared_regions[node]) >= 2 else 'single'


OT_TARGET = 'https://platform.opentargets.org/target/'
OT_DISEASE = 'https://platform.opentargets.org/disease/'


def condition_id(c):
    """The ontology term of a condition row ('' for a curated grouping without one)."""
    if not str(c.node_id).startswith('custom'):
        return c.node_id
    ids = c.ot_ids.split()
    return ids[0] if len(ids) == 1 else ''


def p_text(x):
    """6.8e-27 -> rich parts 6.8 × 10⁻²⁷."""
    v = float(x)
    e = int(math.floor(math.log10(v)))
    m = v / 10 ** e
    return [f'{m:.1f} × 10', ('sup', f'−{-e}')] if e < -2 else [f'{v:.2g}']


# ── Guide ────────────────────────────────────────────────────────────────────
def answers(D):
    st, cond = D['tidy'], D['conditions']
    L, rg = D['ladder'], D['regions']
    bs = st['conditions_by_source']
    lv = lambda k: italic_list(list(L[L.level == k].symbol))  # noqa: E731
    shared = sorted(((D['rel_n'].get(g, 0), g) for g in L.symbol if D['rel_n'].get(g, 0)), reverse=True)
    n_multi = int(((cond.n_shared_regions.astype(int) >= 2)).sum())
    cst = D['cst']
    des = cst.design_class_label.value_counts()
    mech = D['summ'][D['summ'].kind == 'mechanism'].assign(n=lambda d: num(d.n_targets)).sort_values('n',
                                                                                                   ascending=False)
    eas_rel = int(D['relations'].east_asian.map(b).sum())
    full = cond[cond.node_id.isin(D['in_sheet5'])]
    ch = D['changes_json']
    n_choices = sum(len(g['items']) for g in D['choices'])
    by_type = pd.Series([it['basis_type'] for g in D['choices'] for it in g['items']]).value_counts()
    return {
        SN['gene_map']: ['Genome-wide significant in the CTEPH GWAS (G): '] + lv(1) +
        ['; only after leveraging PE (G*): '] + lv(2) + ['; associated in two or more studies: '] + lv(3) + ['.'],
        SN['genes']: [f'{len(cst)} genetic studies of CTEPH, {st["cteph_findings"]} findings, each with its source '
                      'sentence; the 9 published GWAS regions with their statistics.'],
        SN['conditions']: [f'{st["conditions"]} conditions: {bs["clinical"]} with clinical evidence quoted from PubMed, '
                           f'{bs["ontology"]} by the disease ontology, {bs["colocalization"]} by colocalization and '
                           f'{bs["lead-variant overlap"]} by lead-variant overlap with a CTEPH GWAS region.'],
        SN['regions']: [f'{len(rg)} published CTEPH regions; Open Targets has a CTEPH credible set at '
                        f'{int(rg.cteph_study_locus_id.ne("").sum())}. {n_multi} conditions share two or more regions.'],
        SN['shared']: [f'{len(L)} genes studied genetically in CTEPH ({len(D["cteph_genes"])} with a positive '
                       f'result); {sum(1 for g in L.symbol if D["rel_n"].get(g))} of them also have evidence in a '
                       'related condition, most widely '] + italic_list([g for _n, g in shared[:6]]) +
                      [f'; {sum(1 for g in L.symbol if not D["rel_n"].get(g))} only in CTEPH.'],
        SN['cond_gen']: [f'{len(full)} conditions with genetic data: {st["gwas_studies"]} GWAS studies '
                         f'({st["eas_studies"]} East Asian), {st["credible_sets"]:,} credible sets, '
                         f'{st["mendelian_links"]} curated gene links.'],
        SN['mech']: [f'{len(mech)} Reactome themes of the text-mined CTEPH targets; the largest: ' +
                     ', '.join(f'{r.title} ({int(r.n)})' for r in mech.head(4).itertuples()) + '.'],
        SN['eas']: [f'{eas_rel} clinical evidence items, {int(cst.eas.sum())} genetic studies of CTEPH '
                    f'({int(cst.jp.sum())} Japanese) and {st["eas_studies"]} GWAS of related conditions.'],
        SN['methods']: ['Designs: ' + ', '.join(f'{k} {v}' for k, v in des.head(4).items()) + '; methods from the full '
                        f'text for {int(cst.methods_source.eq("full text").sum())} studies, from the abstract for '
                        f'{int(cst.methods_source.eq("abstract").sum())}.'],
        SN['curation']: [f'{len(D["cur"]["excluded"])} publications left out (retracted or conference material), '
                         f'{D["cur"]["studies_excluded"]} of them genetic studies of CTEPH; '
                         f'{int((D["cfi"].association_with != "CTEPH").sum())} findings not about CTEPH itself; '
                         f'{D["cur"]["textmining_dropped"]} text-mined links dropped.'],
        SN['choices']: [f'{n_choices} choices, each with its basis: ' + ', '.join(
            f'{by_type.get(k, 0)} {k}' for k in ('source', 'data check', 'design', 'user')) + '.'],
        SN['changes']: ['First snapshot: nothing to compare with.' if ch['previous'] is None else
                        f'Since {ch["previous"]}: ' + ('; '.join(f'{k.replace("|", " ")} {v}' for k, v in
                                                             ch['counts'].items()) or 'no change') + '.' +
                        (' Not compared (new tables): ' + ', '.join(ch['not_compared']) + '.'
                         if ch.get('not_compared') else '')],
    }


def s_guide(wb, D, foot):
    sh = new_sheet(wb, 'Guide', '#1A1A1A')
    st, ch = D['tidy'], D['changes_json']
    title(sh, 'CTEPH knowledge graph: CTEPH, the conditions related to it, and their genetics',
          f'Open Targets Platform release {st["release"]} and PubMed, fetched {st["fetch_date"]}. External knowledge '
          'only. The sheets follow one line in four parts: A, CTEPH\'s own genetics; B, the conditions related to it '
          'and what they share; C, context; D, provenance (what was left out or corrected, the basis of every '
          'choice, the changes).')
    ans = answers(D)
    rows, links = [], {}
    prev = None
    for i, (k, s, q) in enumerate(SHEETS):
        part = PART_OF[k][0]
        rows.append([[('b', part)] if part != prev else part, s, q, ans[s]])
        prev = part
        links[(i, 1)] = f"internal:'{s}'!A1"
    r = table(sh, 3, [{'h': 'Part', 'w': 22}, {'h': 'Sheet', 'w': 24}, {'h': 'Question', 'w': 44, 'wrap': True},
                      {'h': 'Answer in brief', 'w': 100, 'wrap': True}], rows, links, freeze=None)
    rows, links = [], {}
    for i, (s, d) in enumerate(APPENDIX):
        rows.append([s, d, f'{D["appendix_rows"][s]:,}'])
        links[(i, 0)] = f"internal:'{s}'!A1"
    r = table(sh, r + 1, [{'h': 'Appendix', 'w': 24}, {'h': 'Flat table, with filters', 'w': 44, 'wrap': True},
                          {'h': 'Records', 'w': 100}], rows, links, freeze=None)
    items = [('', [('b', 'How to read')]),
             ('', ['   Evidence ladder of a CTEPH gene: 1 G, genome-wide significant in the CTEPH GWAS; 2 G*, '
                   'significant only after leveraging pulmonary embolism; 3 associated in two or more studies; '
                   '4 in one study; 5 a variant identified without a test; 6 tested, not associated; 7 other. Levels '
                   '3–6 count findings about CTEPH itself; an association with an outcome within CTEPH is shown but '
                   'does not raise the level.']),
             ('', ['   Rows are nested: group > record > source. Open or close them with + / − at the left, or show a '
                   'whole level with the buttons 1 2 3 above the row numbers.']),
             ('', ['   Shaded bold rows are groups; shaded rows open to their records; gray rows hold source quotes '
                   'and sentences. Blue underlined text is a link: to PubMed, or to the same gene or condition on '
                   'another sheet. A cell with a red corner has a note: hover over it.']),
             ('', ['   Gene codes: G, GWAS (for CTEPH, a published genome-wide significant region; for a related '
                   'condition, the top locus-to-gene gene of a credible set); G*, a CTEPH region significant only '
                   'after leveraging PE; C, curated gene evidence; L, a genetic study of CTEPH in the literature.']),
             ('', ['   Shared GWAS regions: ●, a GWAS signal of the condition colocalizes with the CTEPH credible '
                   'set; ○, a GWAS credible set of the condition contains the CTEPH lead variant (regions without '
                   'a CTEPH credible set; weaker).']),
             ('', [f'   Every threshold and grouping has its basis in sheet {SN["choices"]}. Paragraphs headed '
                   '"Generated summary" were written by a language model from the records listed with them and '
                   'checked against them.']),
             ('', ['   ', ('b', MISSING), ' marks a value that is not available; the reason is given under each '
                   'table.'])]
    r = notes(sh, r + 1, items, 168)
    srows = [['Open Targets release', st['release'], ''], ['Fetched', st['fetch_date'], ''],
             ['Previous snapshot', ch['previous'] or 'none', ''] if ch['previous'] else
             ['Previous snapshot', 'none', 'this is the first'],
             ['Published CTEPH GWAS regions', f'{st["published_regions"]}',
              f'Open Targets credible sets: {st["published_with_credible_set"]}'],
             ['Related conditions', f'{st["conditions"]}', ''], ['GWAS studies / credible sets',
                                                                 f'{st["gwas_studies"]} / {st["credible_sets"]:,}', ''],
             [f'{CURATED} gene links', f'{st["mendelian_links"]}', ''],
             ['PubMed records / CTEPH studies', f'{st["cteph_genetic_records"]} / {st["cteph_genetic_studies"]}',
              'records of the genetic-study search / genetic studies of CTEPH kept'],
             ['Knowledge-graph nodes / edges', f'{D["graph"]["n_nodes"]:,} / {D["graph"]["n_edges"]:,}',
              f'{len(D["graph"]["counts"]["nodes"])} node and {len(D["graph"]["counts"]["edges"])} edge types'],
             ['Publications checked / left out', f'{D["cur"]["publications"]:,} / {len(D["cur"]["excluded"])}',
              f'PubMed, {D["cur"]["checked_at"][:10]}']]
    r = table(sh, r + 1, [{'h': 'Snapshot', 'w': 24}, {'h': 'Value', 'w': 44}, {'h': 'Note', 'w': 100}], srows,
              freeze=None)
    rel, cst = D['relations'], D['cst']
    cl = rel[rel.evidence_type == 'clinical']
    src = [['Open Targets Platform', f'release {st["release"]}, read {st["fetch_date"]}',
            'GWAS credible sets and locus-to-gene genes, colocalization, curated gene evidence, text-mined CTEPH '
            'literature, the disease ontology'],
           ['PubMed', f'searched up to {st["fetch_date"]}; records checked {D["cur"]["checked_at"][:10]}',
            f'{st["cteph_genetic_records"]} records of the genetic-study search; the publication type, issue and '
            f'first-author affiliation of all {D["cur"]["publications"]:,} publications the graph uses'],
           ['Europe PMC', f'read {st["fetch_date"]}',
            f'the Methods sections of {int(cst.methods_source.eq("full text").sum())} open-access articles; '
            'affiliations PubMed lacks'],
           ['Liley 2024, Table 1', 'PMID 38470220', 'the nine genome-wide significant CTEPH GWAS regions, as published'],
           ['Curated clinical evidence', f'{len(cl)} quotes from {cl.pmid.nunique()} publications',
            'every quote checked against its PubMed abstract at every build (sheet ' + SN['conditions'].split()[0] +
            ')']]
    links = {(0, 0): 'https://platform.opentargets.org/', (1, 0): 'https://pubmed.ncbi.nlm.nih.gov/',
             (2, 0): 'https://europepmc.org/', (3, 0): pub_url('38470220')}
    r = table(sh, r + 1, [{'h': 'Source', 'w': 24}, {'h': 'Version or date', 'w': 44},
                          {'h': 'Used for', 'w': 100, 'wrap': True}], src, links, freeze=None)
    items = [('', [('b', 'Caveats')])] + [('', ['   ' + c]) for c in D['caveats']]
    items.append(abbrev(['CTEPH', 'GWAS', 'PE', 'PMID']))
    r = notes(sh, r + 1, items, 168)
    print_setup(sh.ws, r, 3, foot)
    return {'header': 3, 'rows': len(QUESTIONS), 'spanner': False}


# ── 1 CTEPH genes at a glance: evidence level x function ─────────────────────
LEVEL_FILL = {1: '#C6DBEF', 2: '#E3EEF8', 3: '#FCE4CC', 4: '#FDF2E6'}


def gene_tip(x):
    """The evidence behind one gene of the ladder, for a cell note."""
    out = [f'{x.symbol}: level {x.level}, {x.level_label.lower()}']
    if x.gwas_region:
        out.append(f'GWAS region {x.gwas_region} ({x.gwas_rsid}), tier {x.gwas_tier}: P {float(x.gwas_p):.1e}, '
                   f'V {float(x.gwas_v):.1e}, OR {x.gwas_or}')
    if int(x.n_studies):
        t = tally_text(x)
        out.append(f'Literature: {t} ({x.n_studies} studies'
                   + (f', {x.eas_studies} East Asian' if int(x.eas_studies) else '') + ')')
    return '\n'.join(out)


def spec_gene_map(D):
    L = D['ladder'].copy()
    order = [t['name'] for t in D['themes']] + [OTHER_GENES]
    present = [t for t in order if (L.theme == t).any()]
    cols = [{'h': 'Evidence level', 'w': 34, 'wrap': True}] + \
        [{'h': D['short'][t], 'w': 13, 'wrap': True, 'a': 'center'} for t in present]
    rows = []
    for lv in sorted(L.level.unique()):
        sub = L[L.level == lv]
        code, lab = LADDER[lv]
        cells, tips, fills = [[('b', f'{lv}  ')] + ([('b', f'{code}  ')] if code else []) + [lab]], {}, {}
        for j, t in enumerate(present, start=1):
            g = sub[sub.theme == t]
            cells.append(italic_list(list(g.symbol)) if len(g) else '')
            if len(g):
                tips[j] = '\n\n'.join(gene_tip(x) for x in g.itertuples())
        if lv in LEVEL_FILL:
            fills = {j: (LEVEL_FILL[lv], INK) for j in range(len(cols))}
        rows.append(R(0, 'item', cells, tips=tips, fills=fills, bold=False))
    rg = D['regions']
    lit = L[L.level >= 3]
    brief = [['Each gene stands at its strongest evidence (rows) and under its function (columns). Hover over a cell '
              'for the evidence behind each gene; sheet ' + SN['genes'].split()[0] + ' lists every study.'],
             [f'GWAS: the {len(rg)} genome-wide significant regions of the only CTEPH GWAS (Liley 2024), named by the '
              'gene the publication reports; regions significant only after leveraging pulmonary embolism are G*.'],
             [f'Literature: {len(lit)} further genes reported by {D["tidy"]["cteph_genetic_studies"]} genetic studies '
              'of CTEPH, ranked by what the studies found.']]
    nts = [('', [('b', 'Evidence levels')])]
    nts += [('', [f'   {k}  ' + (f'{v[0]}  ' if v[0] else '') + v[1] + '.']) for k, v in LADDER.items()]
    nts.append(('', ['   ' + basis(D, 'ladder') + ' Studies are counted, not independent samples. Levels 3–6 count '
                     'only findings about CTEPH itself; an association with an outcome within CTEPH (e.g. ACE, GNB3) '
                     'is shown on sheet ' + SN['genes'].split()[0] + ' but does not raise the level. ' +
                     basis(D, 'about')]))
    nts.append(('', [('b', 'Functions (gene themes) and their basis')]))
    for t in D['themes']:
        refs = ', '.join(ref_short(D['refs'][k]) for k in t.get('refs', []))
        nts.append(('', ['   ', ('b', D['short'][t['name']]), f': {t["basis"]} ({refs}).']))
    nts.append(('', ['   ', ('b', OTHER_GENES), ': genes no source above groups.']))
    nts.append(abbrev(['CTEPH', 'GWAS', 'PE', 'VWF', 'MPN', 'PAH']))
    return spec(SN['gene_map'], QT('gene_map'), 'Evidence level × function. Levels 1–2 from the CTEPH GWAS, 3–7 from '
                'the genetic studies in PubMed.', cols, rows, brief=brief, notes=nts)


# ── 2 CTEPH gene evidence: level > gene > study ──────────────────────────────
def spec_cteph_genes(D):
    L, fx, rg = D['ladder'].copy(), D['fx'], D['regions']
    cols = [{'h': 'Level › gene › study', 'w': 28, 'wrap': True}, {'h': 'Function', 'w': 16, 'wrap': True},
            {'h': 'Design', 'w': 17, 'wrap': True}, {'h': 'Place', 'w': 11, 'wrap': True, 'n': 'c'},
            {'h': 'Cases', 'w': 7, 'a': 'right'}, {'h': 'Controls', 'w': 8, 'a': 'right'},
            {'h': 'Variant', 'w': 15, 'wrap': True, 'n': 'd'}, {'h': 'Result', 'w': 11, 'n': 'a'},
            {'h': 'About', 'w': 11, 'wrap': True, 'n': 'e'}, {'h': 'Effect', 'w': 17, 'wrap': True},
            {'h': 'P', 'w': 10, 'wrap': True}, {'h': 'Source sentence', 'w': 54, 'wrap': True, 'n': 'b'},
            {'h': 'Link', 'w': 13}]
    nc = len(cols)
    rows = []
    gwas_rows = {r.gene: r for r in rg.itertuples()}

    def finding_rows(ff):
        return [R(2, 'quote', [y.cite or f'PMID {y.pmid}', '', y.design_class_label or None, y.place or None,
                               count(y.n_cases) if y.n_cases else '', count(y.n_controls) if y.n_controls else '',
                               y.variant_name or '', y.result, ABOUT_SHORT.get(y.association_with,
                                                                              y.association_with),
                               y.effect or '', y.p_value or '', [('gi', y.quote)], f'PMID {y.pmid}'],
                  links={nc - 1: pub_url(y.pmid)}, hide=True,
                  tips={6: f'As written: {y.variant}' + (f'; rsID {y.rsid}' if y.rsid else '')}
                  if (y.variant and y.variant != y.variant_name) or y.rsid else None)
                for y in ff.itertuples()]
    for lv in sorted(L.level.unique()):
        sub = L[L.level == lv]
        code, lab = LADDER[lv]
        head = f'Level {lv}' + (f' · {code}' if code else '') + f' · {lab} ({len(sub)})'
        what = {1: 'CTEPH meta-analysis P < 5 × 10⁻⁸ (Liley 2024, Table 1).',
                2: 'Genome-wide significant only in the co-analysis with pulmonary embolism, V < 5 × 10⁻⁸ (Liley 2024).',
                3: 'An association with CTEPH in two or more genetic studies.',
                4: 'An association with CTEPH in one genetic study.',
                5: 'Variants found in people with CTEPH (case reports, sequencing) without an association test.',
                6: 'Tested and not associated with CTEPH in every study that tested the gene.',
                7: 'Only other results: an association with an outcome within CTEPH, or a functional or descriptive '
                   'result.'}[lv]
        rows.append(R(0, 'group0', [head, what] + [''] * (nc - 2), span=(1, nc - 1)))
        for x in sub.itertuples():
            ff = fx[fx.symbol == x.symbol]
            txt = []
            if x.gwas_region:
                rr = gwas_rows[x.gwas_region]
                txt += [f'GWAS region {x.gwas_region} ({x.gwas_rsid}, tier {x.gwas_tier}): P '] + p_text(x.gwas_p) + \
                    [', V '] + p_text(x.gwas_v) + [f', OR {x.gwas_or}; Open Targets credible set: '
                                                   f'{"yes" if rr.cteph_study_locus_id else "no"}. ']
            if x.n_studies:
                tt = tally_text(x) or 'reported without a result on CTEPH itself'
                txt += [tt[0].upper() + tt[1:] + f' ({x.n_studies} '
                        f'{"study" if x.n_studies == 1 else "studies"}' +
                        (f', {x.eas_studies} East Asian' if x.eas_studies else '') + '). ']
            n_rel = D['rel_n'].get(x.symbol, 0)
            if n_rel:
                txt += [f'Genetic evidence in {n_rel} related {"condition" if n_rel == 1 else "conditions"}.']
            link = ('@', SN['shared'], x.symbol) if x.symbol in D['gene_map_genes'] else None
            lk = {0: OT_TARGET + x.ensembl_id} if x.ensembl_id else {}
            if link:
                lk[nc - 1] = link
            rows.append(R(1, 'group1', [[('bi', x.symbol)], D['short'].get(x.theme, x.theme), txt] + [''] * (nc - 4) +
                          ['map →' if link else ''], span=(2, nc - 2), links=lk, key=(SN['genes'], x.symbol)))
            rows += finding_rows(ff)
    rest = fx[fx.symbol.eq('')]
    if len(rest):
        rows.append(R(0, 'group0', [f'Gene sets and findings without a single gene ({len(rest)})',
                                    'Findings that name several genes at once, or no gene.'] + [''] * (nc - 2),
                      span=(1, nc - 1)))
        for lab, ff in rest.groupby('label', sort=False):
            rows.append(R(1, 'group1', [[('bi', lab)], D['short'].get(ff.theme.iloc[0], ff.theme.iloc[0]),
                                        f'{ff.pmid.nunique()} {"study" if ff.pmid.nunique() == 1 else "studies"}.']
                          + [''] * (nc - 3), span=(2, nc - 1)))
            rows += finding_rows(ff)
    brief = [[f'{len(D["cst"])} genetic studies of CTEPH ({D["tidy"]["cteph_genetic_records"]} PubMed records '
              f'screened), {len(fx)} findings, and the {len(rg)} published GWAS regions, by evidence level. Open a gene '
              '(+) for each finding with its source sentence; a gene symbol opens its Open Targets page, a PMID its '
              'PubMed record.'],
             ['Level definitions and the function of each gene: sheet ' + SN['gene_map'].split()[0] + '; their basis: '
              'sheet ' + SN['choices'].split()[0] + '.']]
    nts = [('a', ['As the study reports it: associated, not associated, enriched (a group of variants), identified (a '
                  'variant found, without a test), or other.']),
           ('b', ['Verbatim from the abstract or the open-access full text; every gene, variant, number and quote was '
                  'checked against it.']),
           ('c', ['Country of the first author\'s institution (PubMed affiliation). ' + basis(D, 'place')]),
           ('d', ['One name per variant, with its rsID where established; the cell note holds the name as written. '
                  + basis(D, 'variants')]),
           ('e', ['What the finding is about: CTEPH itself, or an outcome or phenotype within CTEPH, another trait, '
                  'a functional experiment or a background statement. Only findings about CTEPH set the level. '
                  + basis(D, 'about')]),
           ('', ['   ' + basis(D, 'g_codes') + ' ' + basis(D, 'ladder')])]
    nts += [('', [('b', 'Missing values')]),
            ('', ['   Empty cells: not stated in the source sentence or the abstract.'])]
    nts += missing_block([('Design', 'not stated in the abstract'),
                          ('Place', 'no affiliation in PubMed or Europe PMC')])
    nts.append(abbrev(['CTEPH', 'GWAS', 'PMID', 'HLA', 'CH', 'PE', 'OR', 'rsID']))
    return spec(SN['genes'], QT('genes'), 'Evidence level › gene › finding, with the GWAS statistics of each region.',
                cols, rows, brief=brief, notes=nts, summary_closed=True,
                summaries=summary_block(D, 'findings', 'Generated summary of the findings'))


# ── 3 related conditions ─────────────────────────────────────────────────────
def regions_parts(D, node):
    """The CTEPH regions a condition shares: colocalizing ones plain, lead-variant overlaps marked °."""
    rel = D['relations']
    r = rel[rel.node_id == node]
    co = sorted({g for v in r[r.evidence_type == 'colocalization'].loci for g in v.split()})
    ov = sorted({g for v in r[r.evidence_type == 'lead-variant overlap'].loci for g in v.split()} - set(co))
    out = []
    for k, g in enumerate(co + ov):
        out += ([', '] if k else []) + [('i', g)] + (['°'] if g in ov else [])
    return out


def spec_conditions(D):
    cond, rel = D['conditions'].copy(), D['relations'].copy()
    gc = D['gc']
    gene_sets = gc.groupby('node_id').symbol.apply(set).to_dict()
    cols = [{'h': 'Category › condition › evidence', 'w': 44, 'wrap': True},
            {'h': 'Clinical items', 'w': 8, 'a': 'center', 'n': 'a'},
            {'h': 'Ontology', 'w': 8, 'a': 'center', 'n': 'b'},
            {'h': 'CTEPH regions shared', 'w': 18, 'wrap': True, 'n': 'c'},
            {'h': 'East Asian evidence', 'w': 9, 'a': 'center'},
            {'h': 'GWAS studies', 'w': 8, 'a': 'right'},
            {'h': 'EAS studies', 'w': 7, 'a': 'right', 'n': 'd'},
            {'h': f'{CURATED} genes', 'w': 8, 'a': 'center', 'n': 'e'},
            {'h': 'Genes shared with CTEPH', 'w': 26, 'wrap': True, 'n': 'f'},
            {'h': 'Identifier', 'w': 15, 'n': 'h'},
            {'h': 'Link', 'w': 14, 'n': 'g'}]
    nc = len(cols)
    rel['_t'] = rel.evidence_type.map({'clinical': 0, 'ontology': 1, 'colocalization': 2, 'lead-variant overlap': 3})
    cond['_g'] = [rel_category(D, n) for n in cond.node_id]
    cond['_c'] = cond.n_clinical.astype(int)
    cond['_r'] = cond.n_shared_regions.astype(int)
    rows = []
    for g in CATEGORY_ORDER + ['multi', 'single']:
        sub = cond[cond._g == g].assign(_n=lambda d: d.name.str.lower()).sort_values(
            ['_c', '_r', '_n'], ascending=[False, False, True])
        if not len(sub):
            continue
        label = g if g in CATEGORY_ORDER else ('Shared GWAS signal: two or more CTEPH regions' if g == 'multi' else
                                               'Shared GWAS signal: one CTEPH region')
        rows.append(R(0, 'group0', [f'{label} ({len(sub)})', CATEGORY_DESC[g]] + [''] * (nc - 2), span=(1, nc - 1)))
        for c in sub.itertuples():
            ev = rel[rel.node_id == c.node_id].sort_values(['_t', 'pmid'])
            ea = ev.east_asian.map(b).any()
            shared = sorted(gene_sets.get(c.node_id, set()) & D['cteph_genes'])
            full = c.genetics == 'full'
            link = ('@', SN['cond_gen'], c.node_id) if c.node_id in D['in_sheet5'] else None
            oid = condition_id(c)
            cells = [c.name, count(c.n_clinical) or '', 'yes' if 'ontology' in c.sources else '',
                     regions_parts(D, c.node_id) or '', 'yes' if ea else '',
                     count(c.n_gwas_studies) if full else '', count(c.n_eas_studies) if full else '',
                     count(c.n_mendelian_genes) if full else '', italic_list(shared) if shared else '',
                     oid or None, 'genetics →' if link else '']
            lk = {nc - 2: OT_DISEASE + oid} if oid else {}
            if link:
                lk[nc - 1] = link
            rows.append(R(1, 'group1', cells, links=lk, hide=g == 'single', key=(SN['conditions'], c.node_id)))
            for x in ev.itertuples():
                if x.evidence_type == 'clinical':
                    lab = f'Clinical: {x.relation}' if x.relation else 'Clinical'
                    txt = [x.detail] + (['\n', ('gi', f'“{x.quote}”')] if x.quote else []) + \
                        ([('g', '  (' + '; '.join(y for y in (x.design, x.population) if y) + ')')]
                         if (x.design or x.population) else [])
                elif x.evidence_type == 'ontology':
                    lab, txt = 'Ontology', [x.detail]
                else:
                    ovl = x.evidence_type == 'lead-variant overlap'
                    lab = 'Lead-variant overlap' if ovl else 'Colocalization'
                    n = int(x.n_studies or 0)
                    txt = ([f'{n} GWAS credible {"set contains" if n == 1 else "sets contain"} the CTEPH lead '
                            f'variant of '] if ovl else
                           [f'{n} GWAS {"signal colocalizes" if n == 1 else "signals colocalize"} with the CTEPH '
                            f'{"region" if len(x.loci.split()) == 1 else "regions"} '])
                    for k, v in enumerate(x.loci.split()):
                        txt += ([', '] if k else []) + [('i', v)]
                    txt += [(f'; highest PIP of the lead {float(x.strength):.2f}.' if ovl else
                             f'; strongest CLPP {float(x.strength):.2f}.') if x.strength else '.']
                ref = (f'{x.citation.split(",")[0]}' if x.pmid else x.citation) or None
                rows.append(R(2, 'quote', [lab] + [txt] + [''] * (nc - 3) + [ref], span=(1, nc - 2),
                              links={nc - 1: pub_url(x.pmid)} if x.pmid else None, hide=True))
    bs = D['tidy']['conditions_by_source']
    brief = [[f'{len(cond)} conditions: {bs["clinical"]} with clinical evidence (quotes from PubMed abstracts, '
              f'checked at every build), {bs["ontology"]} by the disease ontology, {bs["colocalization"]} by '
              f'colocalization and {bs["lead-variant overlap"]} by lead-variant overlap with a CTEPH GWAS region. A '
              'condition can have several kinds.'],
             ['Within each category, conditions are ordered by their evidence: clinical items, then CTEPH regions '
              'shared. Conditions sharing a single CTEPH region only are closed (level 1 button to open).'],
             ['Contrary findings are kept: against survivors of pulmonary embolism, most classic risk factors were '
              'not enriched (open Risk factor and Comorbidity).']]
    nts = [('a', ['Clinical evidence items: guideline, review or study sentences quoted from PubMed abstracts. ' +
                  basis(D, 'clinical_rule')]),
           ('b', ['A relation in the disease ontology: the parent of CTEPH, its siblings, or a child of pulmonary '
                  'arterial hypertension.']),
           ('c', ['CTEPH GWAS regions, named by the gene the publication reports. Plain: a GWAS signal of the condition '
                  f'colocalizes with the CTEPH credible set ({D["choice"]["clpp"]["value_text"]}). Marked °: a GWAS '
                  'credible set of the condition contains the CTEPH lead variant '
                  f'({D["choice"]["overlap"]["value_text"]}), for regions without a CTEPH credible set. ' +
                  basis(D, 'clpp') + ' ' + basis(D, 'overlap')]),
           ('d', [f'GWAS studies with a mostly East Asian sample ({D["choice"]["eas"]["value_text"]}).']),
           ('e', ['Genes with curated gene-disease evidence for the condition or its descendants '
                  f'({D["choice"]["curated"]["value_text"]}). ' + basis(D, 'curated')]),
           ('f', ['Genes with genetic evidence in both the condition and CTEPH (sheet ' + SN['shared'].split()[0] +
                  ').']),
           ('g', ['Condition rows: the condition on sheet ' + SN['cond_gen'].split()[0] + '; evidence rows: the PubMed '
                  'record.']),
           ('h', ['The ontology term of the condition, linked to its Open Targets page (EFO, MONDO, HP, Orphanet or '
                  'OBA).']),
           ('', ['Categories: ' + basis(D, 'categories')])]
    nts += [('', [('b', 'Missing values')]),
            ('', ['   Empty cells: no evidence of that kind, or no genetic data fetched (a procedure, a phenotype too '
                  'broad, or a trait taken only for its shared signal).'])]
    nts += missing_block([('Link', 'no PubMed record (ontology and shared-signal evidence)'),
                          ('Identifier', 'a grouping of several conditions or a procedure, without one ontology term')])
    nts.append(abbrev(['CTEPH', 'GWAS', 'EAS', 'CLPP', 'PIP', 'PH']))
    return spec(SN['conditions'], QT('conditions'),
                'Category › condition › evidence. Open a condition (+) for its quotes, ontology relation and shared '
                'GWAS signals.', cols, rows, brief=brief, notes=nts)


# ── 4 shared GWAS regions ────────────────────────────────────────────────────
def spec_regions(D):
    rg = D['regions'].rename(columns={'or': 'odds_ratio'})
    cond, co, ov = D['conditions'].copy(), D['coloc'], D['overlap']
    cond['nreg'] = cond.n_shared_regions.astype(int)
    cond['_g'] = [rel_category(D, n) for n in cond.node_id]
    names = list(rg.gene)

    def head(r):
        status = ('OT credible set' if r.cteph_study_locus_id else 'OT lead only' if r.variant_id else 'not in OT')
        code = 'G' if r.evidence == 'meta' else 'G*'
        return [('bi', r.gene), '\n', r.rsid, '\n', f'tier {r.tier} · {code}', '\n', ('g', status)]
    cols = [{'h': 'Category › condition', 'w': 42, 'wrap': True}] + \
        [{'h': head(r), 'w': 12, 'a': 'center'} for r in rg.itertuples()] + \
        [{'h': 'Regions shared', 'w': 8, 'a': 'right', 'n': 'a'}]
    nc = len(cols)
    cagg = co.assign(clpp=num(co.clpp)).groupby(['node_id', 'cteph_gene']).agg(n=('study_id', 'nunique'),
                                                                                  best=('clpp', 'max'))
    oagg = ov.assign(pip=num(ov.pip)).groupby(['node_id', 'cteph_gene']).agg(n=('study_id', 'nunique'),
                                                                               best=('pip', 'max'))
    alias = {i: c.node_id for c in cond.itertuples() for i in c.ot_ids.split()}
    cagg = cagg.reset_index().assign(node_id=lambda d: d.node_id.map(lambda i: alias.get(i, i))).groupby(
        ['node_id', 'cteph_gene']).agg(n=('n', 'sum'), best=('best', 'max'))
    oagg = oagg.reset_index().assign(node_id=lambda d: d.node_id.map(lambda i: alias.get(i, i))).groupby(
        ['node_id', 'cteph_gene']).agg(n=('n', 'sum'), best=('best', 'max'))
    rows = []
    for g in CATEGORY_ORDER + ['multi', 'single']:
        sub = cond[(cond._g == g) & (cond.nreg >= 1)].assign(_n=lambda d: d.name.str.lower()).sort_values(
            ['nreg', '_n'], ascending=[False, True])
        if not len(sub):
            continue
        label = g if g in CATEGORY_ORDER else ('Shared GWAS signal only: two or more CTEPH regions' if g == 'multi'
                                               else 'Shared GWAS signal only: one CTEPH region')
        rows.append(R(0, 'group0', [f'{label} ({len(sub)})'] + [''] * (nc - 1), span=(0, nc - 1)))
        for c in sub.itertuples():
            cells, tips, fills = [c.name], {}, {}
            for j, gname in enumerate(names, start=1):
                if (c.node_id, gname) in cagg.index:
                    x = cagg.loc[(c.node_id, gname)]
                    cells.append('●')
                    fills[j] = ('#C6DBEF', INK)
                    tips[j] = f'{gname}: {int(x.n)} GWAS signals of {c.name} colocalize with the CTEPH credible set; ' \
                              f'strongest CLPP {x.best:.2f}.'
                elif (c.node_id, gname) in oagg.index:
                    x = oagg.loc[(c.node_id, gname)]
                    cells.append('○')
                    fills[j] = ('#EEF4FA', INK)
                    tips[j] = f'{gname}: {int(x.n)} GWAS credible sets of {c.name} contain the CTEPH lead variant; ' \
                              f'highest PIP {x.best:.2f} (no CTEPH credible set in Open Targets).'
                else:
                    cells.append('')
            cells.append(int(c.nreg))
            rows.append(R(1, 'item', cells, tips=tips, fills=fills, hide=g == 'single',
                          links={0: ('@', SN['conditions'], c.node_id)}))
    n_ot = int(rg.variant_id.ne('').sum())
    n_cs = int(rg.cteph_study_locus_id.ne('').sum())
    brief = [[f'The {len(rg)} genome-wide significant regions of the CTEPH GWAS (Liley 2024), with their statistics in '
              f'the notes below. Open Targets has a CTEPH credible set at {n_cs} of them, so only these can be '
              f'colocalized (●); for the {n_ot - n_cs} other lead variants it holds, a related trait shares the region '
              f'when one of its GWAS credible sets contains the CTEPH lead variant (○, weaker); '
              f'{len(rg) - n_ot} lead variants are not in Open Targets.'],
             [f'{int((cond.nreg >= 2).sum())} conditions share two or more CTEPH regions, '
              f'{int((cond.nreg == 1).sum())} share one (closed). Hover over a mark for its GWAS signals.']]
    nts = [('a', ['CTEPH regions the condition shares, by either route.']),
           ('', ['● Colocalization: ' + D['choice']['clpp']['value_text'] + ' (or ' + D['choice']['h4']['value_text'] +
                 '). ' + basis(D, 'clpp')]),
           ('', ['○ Lead-variant overlap: ' + D['choice']['overlap']['value_text'] + '. ' + basis(D, 'overlap')]),
           ('', ['Regions and their names: ' + basis(D, 'published_regions') + ' ' + basis(D, 'g_codes')]),
           ('', [('b', 'The regions, as published')])]
    for r in rg.itertuples():
        nts.append(('', ['   ', ('i', r.gene), f' ({r.rsid}, chr{r.chrom}:{r.pos_grch37} GRCh37): P '] + p_text(r.p) +
                    ['; P in PE '] + p_text(r.p_pe) + ['; V '] + p_text(r.v) +
                    [f'; OR {r.odds_ratio}; MAF {r.maf}; tier {r.tier}.']))
    nts.append(('', ['Column headers: the region (gene the publication reports), its lead variant, its tier and code '
                     '(G, CTEPH meta-analysis P < 5 × 10⁻⁸; G*, significant only after leveraging PE), and its status '
                     'in Open Targets: a CTEPH credible set (colocalization possible), the lead variant only (overlap '
                     'possible), or not in Open Targets (neither).']))
    nts.append(abbrev(['CTEPH', 'GWAS', 'CLPP', 'PIP', 'PE', 'OR', 'MAF', 'OT']))
    return spec(SN['regions'], QT('regions'), 'Related condition × CTEPH GWAS region: ● colocalization, ○ lead-variant '
                'overlap.', cols, rows, brief=brief, notes=nts, header_height=64, freeze=1)


# ── 5 shared genes ───────────────────────────────────────────────────────────
def spec_shared_genes(D):
    gc, cond = D['gc'], D['conditions']
    nodes = [CORE] + [n for cat in CATEGORY_ORDER for n in
                      cond[(cond.category == cat) & cond.node_id.isin(set(gc.node_id))]
                      .sort_values('name', key=lambda s: s.str.lower()).node_id]
    names = dict(zip(cond.node_id, cond.name))
    names[CORE] = 'chronic thromboembolic pulmonary hypertension'
    label = {'G': 'GWAS', 'G*': 'GWAS, significant only after leveraging PE', 'C': CURATED, 'L': 'CTEPH literature'}
    cell = {}
    for (g, n), x in gc.groupby(['symbol', 'node_id']):
        cs = codes(x.code)
        tip = '\n'.join(f'{label[k]}: {d}' for k, d in zip(x.code, x.details))
        cell[(g, n)] = (cs, f'{x.symbol.iloc[0]} in {names.get(n, n)}\n{tip}')
    rel_n = D['rel_n']
    lvl = dict(zip(D['ladder'].symbol, D['ladder'].level))
    by_level = lambda gs: sorted(gs, key=lambda g: (lvl[g], -rel_n.get(g, 0), g))  # noqa: E731
    pos, neg = by_level(g for g in lvl if lvl[g] <= 5), by_level(g for g in lvl if lvl[g] >= 6)
    groups = [('CTEPH genes with a positive result (levels 1–5)', pos, False),
              ('CTEPH genes studied without a positive result (levels 6–7)', neg, False),
              ('Other genes in three or more related conditions', sorted([g for g, k in rel_n.items() if k >= 3 and
                                                                          g not in lvl], key=lambda g: (-rel_n[g], g)),
               False),
              ('Other genes in two related conditions', sorted(g for g, k in rel_n.items() if k == 2 and g not in lvl),
               True)]
    cols = [{'h': 'Gene', 'w': 12}, {'h': 'Level', 'w': 4.6, 'a': 'center', 'rot': True, 'n': 'a'},
            {'h': 'Related conditions', 'w': 6, 'a': 'right', 'rot': True, 'n': 'b'}] + \
        [{'h': short(names.get(n, n)), 'w': 4.6, 'a': 'center', 'rot': True,
          'tip': names.get(n, n)} for n in nodes]
    rows = []
    for lab, genes, closed in groups:
        rows.append(R(0, 'group0', [f'{lab} ({len(genes)})'] + [''] * (len(cols) - 1), span=(0, len(cols) - 1)))
        for g in genes:
            cells, fills, tips = [[('i', g)], int(lvl[g]) if g in lvl else '', rel_n.get(g, 0)], {}, {}
            if lvl.get(g, 9) in LEVEL_FILL:
                fills[1] = (LEVEL_FILL[lvl[g]], INK)
            for k, n in enumerate(nodes, start=3):
                cs, tip = cell.get((g, n), ('', ''))
                if n == CORE and not cs and g in lvl:
                    cs = '–'                                           # studied in CTEPH, no positive result
                cells.append(cs)
                if cs and cs != '–':
                    one = cs in KIND_FILL
                    fills[k] = (KIND_FILL[cs], INK) if one else (MULTI_FILL, '#FFFFFF')
                    tips[k] = tip
            links = {0: ('@', SN['genes'], g)} if g in lvl else None
            rows.append(R(1, 'item', cells, links=links, tips=tips, fills=fills, hide=closed,
                          key=(SN['shared'], g)))
    spans, j = [], 4
    for cat in CATEGORY_ORDER:
        k = sum(1 for n in nodes[1:] if cond.set_index('node_id').category[n] == cat)
        if k:
            spans.append((j, j + k - 1, cat))
            j += k
    shared = sorted(((rel_n.get(g, 0), g) for g in lvl if rel_n.get(g, 0)), reverse=True)
    only = by_level(g for g in lvl if not rel_n.get(g))
    ch = D['choice']
    brief = [[f'{len(lvl)} genes were studied genetically in CTEPH: {len(pos)} with a positive result (levels 1–5), '
              f'{len(neg)} without. {len(shared)} of them also have evidence in a related condition, most widely '] +
             italic_list([f'{g} ({n})' for n, g in shared[:8]]) +
             [f'. {len(only)} have none and are specific to CTEPH in this map, among them '] +
             italic_list([g for g in only if lvl[g] <= 3]) + [' (levels 1–3).'],
             ['How a cell is derived. ', ('b', 'G'), ', GWAS: for CTEPH, the gene of a published genome-wide significant '
              'region; for a related condition, the top locus-to-gene gene of one of its GWAS credible sets ('
              f'{ch["l2g"]["value_text"]}). ', ('b', 'G*'), ': a CTEPH region significant only after leveraging '
              'pulmonary embolism. ', ('b', 'C'), f', curated gene–disease evidence ({ch["curated"]["value_text"]}). ',
              ('b', 'L'), ': a genetic study of CTEPH reports an association, an enrichment or a variant in the gene.'],
             ['Rows: every gene studied in CTEPH, by evidence level (– in the CTEPH column: no positive result), '
              'then other genes in three or more related conditions, then two (closed). Columns: CTEPH, then the '
              'related conditions with gene data, by category; traits related by a shared GWAS signal only have no '
              f'gene data by design (sheet {SN["regions"].split()[0]}). Fill: blue G, pale blue G*, green C, red L, gray two kinds. Hover over '
              'a code for the records behind it.']]
    nts = [('a', [f'Evidence level of the gene in CTEPH (sheet {SN["gene_map"].split()[0]}): 1 G; 2 G*; 3 associated in two or more studies; '
                  '4 in one; 5 variant identified without a test; 6 tested, not associated; 7 other result only.']),
           ('b', ['Related conditions (not counting CTEPH) in which the gene has genetic evidence.']),
           ('', ['G for a related condition: ' + basis(D, 'l2g') + ' G and G* for CTEPH: ' + basis(D, 'g_codes')]),
           ('', ['C: ' + basis(D, 'curated')]),
           ('', ['Columns are the conditions with genetic data, by category; hover over a column header for its full '
                 'name.'])]
    nts.append(abbrev(['CTEPH', 'GWAS', 'L2G', 'PE', 'PAH', 'PVOD', 'PCH', 'HUPRA', 'VWF', 'VTE', 'DVT']))
    return spec(SN['shared'], QT('shared'),
                'Gene × condition map of genetic evidence: CTEPH first, then the related conditions by category.',
                cols, rows, brief=brief, notes=nts, spanners=spans, freeze=4)


# ── 10 curation ──────────────────────────────────────────────────────────────
def spec_curation(D):
    cols = [{'h': 'Kind › item', 'w': 30, 'wrap': True}, {'h': 'What', 'w': 70, 'wrap': True},
            {'h': 'Decision', 'w': 40, 'wrap': True, 'n': 'a'}, {'h': 'Detail', 'w': 30, 'wrap': True},
            {'h': 'Link', 'w': 20}]
    pubs, lg, cur, cfi = D['pubs'], D['clog'], D['cur'], D['cfi']
    cite_of = dict(zip(D['cst'].pmid, D['cst'].cite))
    rows = []

    def group(text, desc='', closed=False):
        rows.append(R(0, 'group0', [text, desc + (' Open with + at the left.' if closed else ''), '', '', ''],
                      span=(1, 4) if desc else (0, 4)))
    ex = pubs[pubs.status == 'excluded']        # counted, never listed: they are not part of the graph
    fl = pubs[pubs.status == 'flagged']
    group(f'Publications flagged and kept ({len(fl)})', 'Kept, with the flag shown wherever they are used.', True)
    for x in fl.itertuples():
        rows.append(R(1, 'item', [f'{x.first_author} {x.year}'.strip() or pub_label(x.pmid), x.title or None,
                                  f'kept: {x.reason}', f'used in: {x.roles}', pub_label(x.pmid)],
                      links={4: pub_url(x.pmid)} if pub_url(x.pmid) else None, hide=True))
    fa = cfi[cfi.association_with != 'CTEPH']
    group(f'Findings not about CTEPH itself ({len(fa)})', 'Kept and marked; they do not raise the evidence level of '
          'the gene.')
    for x in fa.itertuples():
        rows.append(R(1, 'item', [cite_of.get(x.pmid, f'PMID {x.pmid}'),
                                  ([('bi', x.symbol), ': '] if x.symbol else []) + [('gi', x.quote)],
                                  f'{x.result}, about {x.association_with}', x.association_note or None,
                                  f'PMID {x.pmid}'], links={4: pub_url(x.pmid)}))
    vr = lg[(lg.step == 'findings') & (lg.rule == 'variant name')]
    group(f'Variants named ({len(vr)})', 'One name per variant, with its rsID where established.', True)
    for x in vr.itertuples():
        pm = x.item.split()[1]
        how = re.match(r'Matched in the (.+?);', x.basis)
        rows.append(R(1, 'item', [cite_of.get(pm, f'PMID {pm}'), [('i', x.item.split(' ', 2)[2] if
                                                                  x.item.count(' ') >= 2 else ''),
                                                                 f'  written: {x.before}'], f'named: {x.after}',
                                  f'matched in the {how.group(1)}' if how else None, f'PMID {pm}'],
                      links={4: pub_url(pm)}, hide=True))
    tm = lg[(lg.step == 'text mining') & ~lg.after.str.startswith('publication left out')]
    tmr = D['tmr']
    n_wrong = int(((tmr.decision == 'dropped') & tmr.reason.str.startswith('mapped to the wrong gene')).sum())
    group(f'Text-mined links mapped to the wrong gene ({len(tm)} genes)', 'Each matched word read against its '
          'sentences.', True)
    for x in tm.itertuples():
        dec = 'dropped' if x.rule == 'link dropped' else 'kept'
        pm = [p_ for p_ in tmr[(tmr.symbol == x.item) & (tmr.decision == dec) & tmr.reason.ne('')].pmid if p_]
        more = f' and {len(pm) - 5} more' if len(pm) > 5 else ''
        rows.append(R(1, 'item', [[('bi', x.item)], x.after, 'link dropped' if x.rule == 'link dropped' else
                                  'words removed, link kept',
                                  f'{x.before}: ' + ', '.join(pub_label(p_) for p_ in pm[:5]) + more,
                                  pub_label(pm[0]) if pm else ''],
                      links={4: pub_url(pm[0])} if pm and pub_url(pm[0]) else None, hide=True))
    pl = lg[(lg.step == 'studies') & lg.rule.isin(['place', 'Japanese flag'])]
    group(f'Place of the genetic studies, from the institution ({len(pl)})', 'The country of the first author\'s '
          'institution; it sets the Japanese and East Asian flags.')
    for x in pl.itertuples():
        pm = x.item.split()[1]
        rows.append(R(1, 'item', [cite_of.get(pm, x.item), {'place': f'Population stated: '
                                                                     f'{x.before.split(": ", 1)[-1]}',
                                                            'Japanese flag': 'Done in Japan (first author\'s '
                                                                             'institution)'}[x.rule],
                                  x.after if x.rule == 'place' else 'Japanese: yes',
                                  'the abstract states another population' if x.rule == 'place' else
                                  'the abstract names no population', f'PMID {pm}'],
                      links={4: pub_url(pm)}))
    dc = D['cst'].design_class_label.value_counts()
    dc = dc.reindex([g for g in DESIGN_GROUPS if g in dc.index])
    group(f'Design classes ({len(D["cst"])} studies)', 'One class per study, the first that applies.')
    for k, v in dc.items():
        rows.append(R(1, 'item', [k, f'{v} {"study" if v == 1 else "studies"}', '', '', '']))
    brief = [[f'Every publication the graph uses was looked up in PubMed on {cur["checked_at"][:10]}. Retracted '
              f'articles and conference material are not part of the graph ({len(ex)}, none of them listed); '
              f'{len(fl)} publications are flagged and kept.'],
             [f'{len(fa)} findings are about something other than CTEPH itself; {len(vr)} variants carry one name; '
              f'{n_wrong} text-mined links mapped to the wrong gene are dropped. No rule calls a language model.']]
    nts = [('a', ['The rule applied. ' + basis(D, 'exclusion') + ' ' + basis(D, 'about') + ' ' +
                  basis(D, 'textmining_review') + ' ' + basis(D, 'place')])]
    nts += missing_block([('What', 'the record has no title in PubMed (Europe PMC preprint identifiers)'),
                          ('Detail', 'nothing further recorded')])
    nts.append(abbrev(['CTEPH', 'PMID', 'PE', 'rsID']))
    return spec(SN['curation'], QT('curation'), 'Kind › item: what was left out or corrected, and the rule.', cols,
                rows, brief=brief, notes=nts)


# ── 11 choices and their basis ───────────────────────────────────────────────
BASIS_FILL = {'source': '#DEEAF6', 'data check': '#E2F0D9', 'design': '#F2F2F2', 'user': '#FCEFD2'}


def spec_choices(D):
    cols = [{'h': 'Group › choice', 'w': 34, 'wrap': True}, {'h': 'Value in this snapshot', 'w': 28, 'wrap': True},
            {'h': 'Basis type', 'w': 10, 'a': 'center', 'n': 'a'}, {'h': 'Basis', 'w': 70, 'wrap': True},
            {'h': 'References', 'w': 24, 'wrap': True}, {'h': 'Used in', 'w': 20, 'wrap': True}]
    rows, refs_used = [], {}
    for g in D['choices']:
        rows.append(R(0, 'group0', [f'{g["name"]} ({len(g["items"])})'] + [''] * 5, span=(0, 5)))
        for it in g['items']:
            txt = [it['basis']]
            if it['check_text']:
                txt += ['\n', ('b', 'Check passed. ' if it['check_ok'] else 'Check FAILED. '), it['check_text']]
            refs = [r[0] for r in it['refs']]
            for r in it['refs']:
                refs_used[r[0]] = r[2]
            link = it['refs'][0][1] if len(it['refs']) == 1 and it['refs'][0][1] else None   # one link per cell
            rows.append(R(1, 'item', [it['item'], it['value_text'], it['basis_type'], txt, '; '.join(refs),
                                      it['used_in']], links={4: link} if link else None,
                          fills={2: (BASIS_FILL[it['basis_type']], INK)}, key=(SN['choices'], it['key'])))
    n = sum(len(g['items']) for g in D['choices'])
    by = pd.Series([it['basis_type'] for g in D['choices'] for it in g['items']]).value_counts()
    brief = [[f'{n} choices decide what enters the knowledge graph. Each states its value in this snapshot (read from '
              'the run itself) and its basis: ' + ', '.join(f'{by.get(k, 0)} {k}' for k in BASIS_FILL) + '.'],
             ['Statements attributed to a source were checked against its full text or documentation; data checks are '
              'recomputed at every snapshot.']]
    nts = [('a', ['Source: a publication or the documentation of a data source. Data check: a check on this '
                  'snapshot\'s data, recomputed at every build. Design: a decision of this project, with its reason. '
                  'User: decided with the user, on the date given.']),
           ('', [('b', 'References')])]
    for short_, full in sorted(refs_used.items()):
        nts.append(('', ['   ', ('b', short_), f'. {full}']))
    nts.append(abbrev(['CTEPH', 'GWAS', 'PE', 'CLPP', 'PIP', 'L2G']))
    return spec(SN['choices'], QT('choices'), 'Group › choice: value, basis and references.', cols, rows, brief=brief,
                notes=nts)


# ── 4 study methods ──────────────────────────────────────────────────────────
def spec_methods(D):
    cst, cmf = D['cst'].copy(), D['cmf']
    known = {f for _s, fs in METHOD_STEPS for f in fs}
    cols = [{'h': 'Design › study', 'w': 24, 'wrap': True, 'n': 'c'}, {'h': 'Place', 'w': 11, 'wrap': True},
            {'h': 'Cases', 'w': 7, 'a': 'right'}, {'h': 'Controls', 'w': 8, 'a': 'right'},
            {'h': 'Methods from', 'w': 9, 'a': 'center', 'n': 'a'}] + \
        [{'h': s, 'w': 20 if s == 'Replication' else 40, 'wrap': True, 'n': 'b' if k == 0 else None}
         for k, (s, _f) in enumerate(METHOD_STEPS)] + [{'h': 'Link', 'w': 13}]
    nc = len(cols)
    cst['group'] = cst.design_class_label
    assert set(cst.group) <= set(DESIGN_GROUPS), sorted(set(cst.group) - set(DESIGN_GROUPS))
    cst['_y'] = num(cst.year)
    rows = []
    for gname in DESIGN_GROUPS:
        sub = cst[cst.group == gname].sort_values(['_y', 'pmid'])
        if not len(sub):
            continue
        rows.append(R(0, 'group0', [f'{gname} ({len(sub)})'] + [''] * (nc - 1), span=(0, nc - 1)))
        for s in sub.itertuples():
            m = cmf[cmf.pmid == s.pmid]
            cells, tips = [], {}
            for k, (step, fields) in enumerate(METHOD_STEPS):
                fs = fields + ([f for f in m.field.unique() if f not in known] if step == 'Analysis' else [])
                mm = m[m.field.isin(fs)].assign(_o=lambda d: d.field.map(lambda f: fs.index(f)))
                mm = mm.sort_values('_o')
                parts, tip = [], []
                for y in mm.itertuples():
                    parts += (['\n'] if parts else []) + ([('b', f'{y.field[0].upper() + y.field[1:]}: ')]
                                                          if step != 'Replication' else []) + [y.value]
                    tip.append(f'{y.field}: “{y.quote}”')
                cells.append(parts or None)
                if tip:
                    tips[5 + k] = f'Source ({s.methods_source}):\n' + '\n\n'.join(tip)
            src = s.methods_source or None
            rows.append(R(1, 'item', [s.cite, s.place,
                                      count(s.n_cases) if s.n_cases else '', count(s.n_controls) if s.n_controls
                                      else '', src] + cells + [f'PMID {s.pmid}'],
                          links={nc - 1: pub_url(s.pmid)}, tips=tips,
                          fills={4: (SOURCE_FILL[src], INK)} if src in SOURCE_FILL else None))
    des = cst.design.str.split('; ').explode().value_counts()
    pops = cst.place.replace(MISSING, 'not known').value_counts()
    nca = num(cst.n_cases)
    gw = D['studies']
    brief = [['Designs (a study can have several): ' + ', '.join(f'{k} {v}' for k, v in des.items()) + '.'],
             ['Places (country of the first author\'s institution): ' + ', '.join(f'{k} {v}' for k, v in pops.items())
              + '. Cases per study: median '
              f'{nca.median():.0f}, range {nca.min():.0f}–{nca.max():,.0f} (where stated).'],
             [f'Methods were read from the open-access Methods section for {int(cst.methods_source.eq("full text").sum())}'
              f' studies and from the abstract for {int(cst.methods_source.eq("abstract").sum())} (no open-access full '
              'text). Hover over a methods cell for its source passages.'],
             [f'GWAS in Open Targets (CTEPH and the related conditions): {len(gw)} studies; fine-mapping of their '
              'credible sets: ' + ', '.join(f'{k} {v:,}' for k, v in D['loci'].finemapping.replace('', 'not given')
                                            .value_counts().items()) + '.']]
    nts = [('a', ['Full text: the Methods section of the open-access article (Europe PMC); abstract: the PubMed '
                  'abstract, used when no full text could be fetched.']),
           ('b', ['Methods as the source states them, by step; every value was checked against its source passage '
                  '(in the cell note).']),
           ('c', ['One design class per study, the first that applies in the order of the groups (case reports '
                  'first). ' + basis(D, 'design_class')])]
    nts += [('', [('b', 'Missing values')]),
            ('', ['   Empty cells: the number of cases or controls is not stated.'])]
    nts += missing_block([('the methods columns', 'the step is not described in the source text'),
                          ('Place', 'no affiliation in PubMed or Europe PMC')])
    nts.append(abbrev(['CTEPH', 'GWAS', 'PMID', 'HLA', 'PICS']))
    return spec(SN['methods'], QT('methods'),
                'Design › study, with its methods by step. Hover over a methods cell for the source passage.',
                cols, rows, brief=brief, summaries=summary_block(D, 'methods', 'Generated summary of the methods'),
                notes=nts)


# ── 5 condition genetics ─────────────────────────────────────────────────────
def spec_condition_genetics(D):
    cond, lo, men, gc = D['conditions'], D['loci'].copy(), D['mendelian'], D['gc']
    lo['_s'], lo['_p'] = num(lo.l2g1_score), num(lo.log10_p)
    lo = lo[lo._s >= D['l2g_min']].sort_values('_p')
    best = lo.groupby(['node_id', 'l2g1_symbol']).agg(n=('study_locus_id', 'size'), m=('p_mantissa', 'first'),
                                                       e=('p_exponent', 'first'), rs=('rsid', 'first'),
                                                       vid=('variant_id', 'first'), s=('_s', 'first'),
                                                       p=('_p', 'first'), eas=('eas', lambda x: x.map(b).any()))
    best = best.to_dict('index')
    cur = {}
    for x in men.itertuples() if len(men) else []:
        cls = x.classifications or x.clinical_significance
        cur.setdefault((x.node_id, x.symbol), []).append(
            f'{x.source}: {cls}' + (f' ({x.review_status})' if getattr(x, 'review_status', '') else '')
            if cls else f'{x.source}: score {float(x.source_score):.2f}')
    cols = [{'h': 'Category › condition › gene', 'w': 40, 'wrap': True}, {'h': 'GWAS loci', 'w': 7, 'a': 'right',
                                                                         'n': 'a'},
            {'h': 'Best P', 'w': 12, 'a': 'right'}, {'h': 'Lead variant', 'w': 13}, {'h': 'L2G score', 'w': 7,
                                                                                    'a': 'right', 'n': 'b'},
            {'h': 'EAS', 'w': 5, 'a': 'center', 'n': 'c'}, {'h': f'{CURATED} evidence', 'w': 44, 'wrap': True,
                                                             'n': 'd'},
            {'h': 'In CTEPH', 'w': 7, 'a': 'center', 'n': 'e'}, {'h': 'Other related conditions', 'w': 9,
                                                                 'a': 'right', 'n': 'f'}]
    nc = len(cols)
    rows = []
    for cat in CATEGORY_ORDER:
        sub = cond[(cond.category == cat) & cond.node_id.isin(D['in_sheet5'])].sort_values(
            'name', key=lambda s: s.str.lower())
        if not len(sub):
            continue
        rows.append(R(0, 'group0', [f'{cat} ({len(sub)})', CATEGORY_DESC[cat]] + [''] * (nc - 2), span=(1, nc - 1)))
        for c in sub.itertuples():
            g = gc[(gc.node_id == c.node_id)]
            genes = sorted(set(g.symbol))
            shared = [x for x in genes if x in D['cteph_genes']]
            n_g, n_c = g[g.code == 'G'].symbol.nunique(), g[g.code == 'C'].symbol.nunique()
            if genes:
                txt = [f'{plural(int(c.n_gwas_studies), "GWAS study", "GWAS studies")} ({c.n_eas_studies} East '
                       f'Asian), {plural(int(c.n_credible_sets), "credible set")}; {plural(len(genes), "gene")} '
                       f'({n_g} by GWAS, {n_c} curated)'] + \
                    ([', shared with CTEPH: '] + italic_list(shared) if shared else []) + ['.']
            else:
                txt = [f'{plural(int(c.n_gwas_studies), "GWAS study", "GWAS studies")}, '
                       f'{plural(int(c.n_credible_sets), "credible set")}; no gene with a locus-to-gene score ≥ '
                       f'{D["l2g_min"]} or curated evidence.']
            rows.append(R(1, 'group1', [c.name, txt] + [''] * (nc - 2), span=(1, nc - 1),
                          key=(SN['cond_gen'], c.node_id)))
            order = sorted(genes, key=lambda x: (x not in D['cteph_genes'], (c.node_id, x) not in cur,
                                                 best[(c.node_id, x)]['p'] if (c.node_id, x) in best else 0, x))
            for x in order:
                bb = best.get((c.node_id, x))
                rel_other = D['rel_n'].get(x, 0) - 1
                cells = [[('i', x)], bb['n'] if bb else '', p_parts(bb['m'], bb['e']) if bb else '',
                         (bb['rs'] or variant_label(bb['vid'])) if bb else '', f'{bb["s"]:.2f}' if bb else '',
                         'yes' if bb and bb['eas'] else '', '\n'.join(cur.get((c.node_id, x), [])) or '',
                         D['cteph_kinds'].get(x, ''), rel_other or '']
                link = ('@', SN['shared'], x) if x in D['gene_map_genes'] else None
                rows.append(R(2, 'item', cells, links={0: link} if link else None, hide=True))
    st = D['tidy']
    brief = [[f'{len(D["in_sheet5"])} conditions with genetic data: {st["gwas_studies"]} GWAS studies ({st["eas_studies"]} '
              f'East Asian), {st["credible_sets"]:,} credible sets, {st["mendelian_links"]} curated gene links.'],
             ['A condition (+) opens its genes: those that are the top locus-to-gene gene of one of its credible sets, '
              'and those with curated evidence. Genes also found in CTEPH come first.'],
             ['The full GWAS records are in the appendix (A1 studies, A2 loci, A4 curated genes).']]
    nts = [('a', ['Credible sets of the condition whose top locus-to-gene gene is this gene; P, lead variant and score '
                  'are those of the strongest.']),
           ('b', [f'Locus-to-gene score (0–1), a prediction of the causal gene; genes from {D["l2g_min"]} upward. ' +
                  basis(D, 'l2g')]),
           ('c', [f'Any of those credible sets from a study with a mostly East Asian sample (≥ {D["eas_min"]:.0%}).']),
           ('d', ['Source: ClinGen or Gene2Phenotype class, PanelApp rating, Orphanet assessment, ClinVar clinical '
                  'significance (review status), UniProt confidence; or the source score when it gives none. Kept '
                  f'from {D["choice"]["curated"]["value_text"]}. ' + basis(D, 'curated')]),
           ('e', ['Genetic evidence in CTEPH: G, a published genome-wide significant GWAS region; G*, significant only '
                  'after leveraging pulmonary embolism; L, a genetic study of CTEPH. ' + basis(D, 'g_codes')]),
           ('f', ['Related conditions other than this one with genetic evidence for the gene.'])]
    nts += [('', [('b', 'Missing values')]), ('', ['   Empty cells: no evidence of that kind.'])]
    nts.append(abbrev(['CTEPH', 'GWAS', 'L2G', 'EAS', 'ClinGen', 'PanelApp']))
    return spec(SN['cond_gen'], QT('cond_gen'),
                'Category › condition › gene: GWAS loci and curated gene evidence of each related condition.',
                cols, rows, brief=brief, notes=nts)


# ── 6 mechanisms ─────────────────────────────────────────────────────────────
def spec_mechanisms(D):
    m = D['summ'][D['summ'].kind == 'mechanism'].copy()
    m = m.assign(_n=-num(m.n_targets)).sort_values('_n')
    tg = D['targets'].copy()
    tg['lit'] = num(tg.ds_europepmc)
    tg = tg[tg.lit.notna()].sort_values('lit', ascending=False)
    L = D['literature'].assign(_s=lambda d: num(d.evidence_score)).sort_values('_s', ascending=False)
    cols = [{'h': 'Theme › target › publication', 'w': 30, 'wrap': True},
            {'h': 'Count', 'w': 8, 'a': 'right', 'n': 'a'},
            {'h': 'Summary › target › sentences', 'w': 110, 'wrap': True, 'n': 'b'}, {'h': 'Link', 'w': 15}]
    rows = []
    for x in m.itertuples():
        src = set(x.genes.split('; '))
        sub = tg[tg.themes.str.split('; ').map(lambda xs: x.key in xs)]
        txt = md_parts(x.summary) if b(x.grounded) else [f'Summary withheld: {x.problems}.']
        rows.append(R(0, 'group1', [x.title, count(x.n_targets), [('b', 'Generated summary. ')] + txt
                                    if b(x.grounded) else txt, ''], bold=True))
        for t in sub.itertuples():
            ll = L[L.ensembl_id == t.ensembl_id]
            detail = [t.name + f'. Literature score {float(t.lit):.2f}; {len(ll)} publications'
                      + (f' ({year_span(t.first_year, t.last_year)})' if year_span(t.first_year, t.last_year) else '')
                      + '.' + (' A source of the summary.' if t.symbol in src else '')]
            if len(ll) and ll.matched_as.iloc[0]:
                detail += [('g', ' Matched as: ' + ', '.join(sorted({w for v in ll.matched_as for w in v.split(' | ')
                                                                     if w})[:6]) + '.')]
            rows.append(R(1, 'group1', [[('i', t.symbol)], len(ll), detail, ''], hide=True))
            for p in ll.itertuples():
                sents = [s for s in p.sentences.split(' || ') if s]
                body = ([('b', p.title)] if p.title else []) + \
                    [('g', ('\n' if (p.title or k) else '') + '• ' + s) for k, s in enumerate(sents)]
                rows.append(R(2, 'quote', [p.year or None, '', body or [('g', 'only the title names both')],
                                           (f'PMID {p.pmid}' if str(p.pmid).isdigit() else f'Europe PMC {p.pmid}')
                                           if p.pmid else None],
                              links={3: pub_url(p.pmid)} if pub_url(p.pmid) else None, hide=True))
    brief = [[f'{len(tg)} CTEPH targets with text-mined literature, in {len(m)} Reactome themes (a target can be in '
              'several). Each theme opens to its targets and their publications with the sentences that name the '
              'target and CTEPH.'],
             ['The summary of each theme was generated from up to two sentences of the best-scored publications of '
              'its best-scored targets, and checked: every PMID it cites is among those, every gene it names occurs '
              'in them.'],
             ['Text mining can map an abbreviation to the wrong gene (for example TF, tissue factor, mapped to '
              'transferrin): check the words a target was matched as.']]
    nts = [('a', ['Theme rows: targets in the theme; target rows: their publications.']),
           ('b', ['Theme rows: the generated summary; target rows: the target, its literature score and the words it '
                  'was matched as; publication rows: the title and up to three sentences, abstract first.'])]
    nts += missing_block([('the first column', 'the publication has no date'),
                          ('Link', 'the record has no identifier')])
    nts.append(abbrev(['CTEPH', 'PMID']))
    return spec(SN['mech'], QT('mech'),
                'Reactome theme › text-mined CTEPH target › publication and its sentences.', cols, rows, brief=brief,
                notes=nts)


# ── 7 East Asian ─────────────────────────────────────────────────────────────
def spec_eas(D):
    cols = [{'h': 'Kind › item', 'w': 34, 'wrap': True}, {'h': 'Evidence', 'w': 80, 'wrap': True},
            {'h': 'Population', 'w': 24, 'wrap': True}, {'h': 'Reference', 'w': 16}]
    names = D['names']
    rows = []
    rel = D['relations']
    ea = rel[rel.east_asian.map(b)]
    rows.append(R(0, 'group0', [f'Clinical evidence ({len(ea)})', '', '', ''], span=(0, 3)))
    for x in ea.itertuples():
        rows.append(R(1, 'item', [names.get(x.node_id, x.node_id), [x.detail] + (['\n', ('gi', f'“{x.quote}”')]
                                                                                 if x.quote else []),
                                  x.population or None, x.citation.split(',')[0] if x.pmid else None],
                      links={3: pub_url(x.pmid)} if x.pmid else None))
    cst = D['cst']
    s = cst[cst.eas].sort_values(['year', 'pmid'])
    rows.append(R(0, 'group0', [f'Genetic studies of CTEPH done in East Asia ({len(s)}; in Japan {int(s.jp.sum())})',
                                '', '', ''], span=(0, 3)))
    fx = D['fx']
    for x in s.itertuples():
        f = fx[fx.pmid == x.pmid]
        res = '; '.join(dict.fromkeys(
            (y.variant_name if y.variant_name.startswith(y.label) else f'{y.label} {y.variant_name}').strip() +
            f' {y.result}' + ('' if y.association_with == 'CTEPH' else f' (about {y.association_with})')
            for y in f.itertuples()))
        rows.append(R(1, 'item', [x.cite, [f'{x.design_class_label}. ', ('g', res or 'no gene-level finding')],
                                  x.place if x.population in (x.place, 'not stated') else f'{x.place} ({x.population})',
                                  f'PMID {x.pmid}'], links={3: pub_url(x.pmid)}))
    st = D['studies']
    e = st[st.eas.map(b)].assign(_n=lambda d: -num(d.n_samples)).sort_values(['condition', '_n'])
    rows.append(R(0, 'group0', [f'GWAS of the related conditions with a mostly East Asian sample ({len(e)})', '',
                                '', ''], span=(0, 3)))
    for x in e.itertuples():
        pub = cite(x) if x.first_author else (x.project or None)
        rows.append(R(1, 'item', [x.condition, f'{x.trait}; {intx(x.n_samples) or "sample size not given"} samples; '
                                  f'{x.n_credible_sets} credible sets ({x.study_id}).', ancestry_label(x.ancestry)
                                  or None, pub], links={3: pub_url(x.pmid)} if pub_url(x.pmid) else None))
    co = D['coloc']
    c = co[co.eas.map(b)].assign(_c=lambda d: -num(d.clpp)).sort_values(['cteph_locus', '_c'])
    rows.append(R(0, 'group0', [f'GWAS signals from East Asian studies that colocalize with a CTEPH locus ({len(c)})',
                                '', '', ''], span=(0, 3)))
    for x in c.itertuples():
        rows.append(R(1, 'item', [x.condition or x.trait, [f'Colocalizes with the CTEPH locus {x.cteph_locus} (',
                                                         ('i', x.cteph_gene), f'), CLPP {float(x.clpp):.2f}; '
                                                         f'{x.trait} ({x.study_id}).'],
                                  ancestry_label(x.ancestry) or None,
                                  ' '.join(y for y in (x.first_author, x.pub_date[:4]) if y) or x.project or None],
                      links={3: pub_url(x.pmid)} if pub_url(x.pmid) else None))
    brief = [['Evidence from East Asian populations, gathered from the other sheets: clinical evidence quoted from '
              'PubMed, genetic studies of CTEPH done in Japan or another East Asian country (the country of the first '
              'author\'s institution, or an East Asian population the abstract states), and GWAS whose sample is at '
              f'least {D["eas_min"]:.0%} East Asian.']]
    nts = missing_block([('Population', 'not stated'), ('Reference', 'no publication listed (consortium release)')])
    nts.append(abbrev(['CTEPH', 'GWAS', 'PMID', 'CLPP']))
    return spec(SN['eas'], QT('eas'), 'Kind of evidence › item.', cols, rows, brief=brief,
                notes=nts)


# ── 8 changes ────────────────────────────────────────────────────────────────
def s_changes(wb, D, foot):
    sh = new_sheet(wb, SN['changes'], '#1F4E79')
    c, ch = D['changes'], D['changes_json']
    sub = (f'Compared with the snapshot of {ch["previous"]} (release {ch["previous_release"]})' if ch['previous']
           else 'First snapshot: nothing to compare with')
    title(sh, QT('changes'), sub)
    sh.ws.write_url(0, 3, "internal:'Guide'!A1", sh.fmt(font_color=LINK, underline=1, align='right'),
                    string='← Guide')
    cols = [{'h': 'Category', 'w': 20}, {'h': 'Change', 'w': 10}, {'h': 'Item', 'w': 44, 'wrap': True},
            {'h': 'Detail', 'w': 80, 'wrap': True}]
    rows = [[x.category, x.change, x.item, x.detail or None] for x in c.itertuples()]
    r = table(sh, 3, cols, rows, freeze=0, autofilter=True)
    items = missing_block([('Detail', 'no further detail recorded')]) if len(c) and c.detail.eq('').any() else []
    r = notes(sh, r + 1, items or [('', [' '])])
    print_setup(sh.ws, r, len(cols), foot)
    return {'header': 3, 'rows': len(c), 'spanner': False}


# ── appendix: flat tables with filters ───────────────────────────────────────
def flat(wb, name, sub, cols, rows, links, nts, foot, freeze=1):
    sh = new_sheet(wb, name, '#A6A6A6')
    title(sh, dict(APPENDIX)[name], sub)
    sh.ws.write_url(0, len(cols) - 1, "internal:'Guide'!A1", sh.fmt(font_color=LINK, underline=1, align='right'),
                    string='← Guide')
    r = table(sh, 3, cols, rows, links, freeze=freeze, autofilter=True)
    r = notes(sh, r + 1, nts)
    print_setup(sh.ws, r, len(cols), foot)
    return {'header': 3, 'rows': len(rows), 'spanner': False}


def a_studies(wb, D, foot):
    s = D['studies'].copy()
    s['_c'] = (s.node_id != CORE).astype(int)
    s['_n'] = num(s.n_samples)
    s = s.sort_values(['_c', 'condition', '_n'], ascending=[True, True, False])
    cols = [{'h': 'Condition', 'w': 28, 'wrap': True}, {'h': 'Study', 'w': 22}, {'h': 'Trait', 'w': 30, 'wrap': True},
            {'h': 'Publication', 'w': 18}, {'h': 'PMID', 'w': 9}, {'h': 'Samples', 'w': 10, 'a': 'right'},
            {'h': 'Cases', 'w': 9, 'a': 'right'}, {'h': 'Controls', 'w': 10, 'a': 'right'},
            {'h': 'Cohorts', 'w': 22, 'wrap': True}, {'h': 'Ancestry', 'w': 22, 'wrap': True, 'n': 'a'},
            {'h': 'EAS', 'w': 6, 'a': 'center', 'n': 'b'}, {'h': 'Summary statistics', 'w': 10, 'a': 'center'},
            {'h': 'Credible sets', 'w': 9, 'a': 'right'}]
    rows, links = [], {}
    for i, x in enumerate(s.itertuples()):
        pub = ' '.join(y for y in (x.first_author, x.year) if y) or x.project or None
        rows.append([x.condition, x.study_id, x.trait or None, pub, x.pmid or None, intx(x.n_samples),
                     intx(x.n_cases), intx(x.n_controls), x.cohorts or None, ancestry_label(x.ancestry) or None,
                     'yes' if b(x.eas) else 'no', 'yes' if b(x.summary_statistics) else 'no', x.n_credible_sets])
        if pub_url(x.pmid):
            links[(i, 4)] = pub_url(x.pmid)
    nts = [('a', ['Population of the study\'s linkage-disequilibrium reference, with its share of the sample.']),
           ('b', [f'East Asian share of the sample ≥ {D["eas_min"]:.0%}.'])]
    nts += missing_block([('Trait, Publication, PMID, Samples, Cases, Controls, Cohorts, Ancestry',
                           'not given by Open Targets for the study')])
    nts.append(abbrev(['CTEPH', 'GWAS', 'EAS', 'PMID']))
    return flat(wb, 'A1 GWAS studies', f'{len(s)} studies, CTEPH first; use the filter buttons of the header row',
                cols, rows, links, nts, foot, freeze=2)


def a_loci(wb, D, foot):
    g = D['loci'].copy()
    g['_c'] = (g.node_id != CORE).astype(int)
    g['_p'] = num(g.log10_p)
    g = g.sort_values(['_c', 'condition', '_p'])
    cols = [{'h': 'Condition', 'w': 26, 'wrap': True}, {'h': 'Study', 'w': 20}, {'h': 'Publication', 'w': 16},
            {'h': 'Lead variant', 'w': 13}, {'h': 'Position (GRCh38)', 'w': 22},
            {'h': 'P', 'w': 12, 'a': 'right', 'n': 'a'}, {'h': 'β', 'w': 8, 'a': 'right', 'n': 'a'},
            {'h': 'Fine-mapping', 'w': 10}, {'h': 'EAS', 'w': 6, 'a': 'center'},
            {'h': 'Locus-to-gene genes (score)', 'w': 30, 'wrap': True, 'n': 'b'}]
    rows = []
    for x in g.itertuples():
        l2g = [(getattr(x, f'l2g{k}_symbol'), num(pd.Series([getattr(x, f'l2g{k}_score')]))[0]) for k in (1, 2, 3)
               if getattr(x, f'l2g{k}_symbol')]
        be = num(pd.Series([x.beta]))[0]
        pub = ' '.join(y for y in (x.first_author, x.year) if y) or x.project or None
        rows.append([x.condition, x.study_id, pub, x.rsid or None, variant_label(x.variant_id) or None,
                     p_parts(x.p_mantissa, x.p_exponent), None if pd.isna(be) else f'{be:.3f}', x.finemapping or None,
                     'yes' if b(x.eas) else 'no', genes_parts(l2g) or None])
    nts = [('a', ['P and β of the credible set\'s lead variant, as reported by the study.']),
           ('b', ['Up to three genes by Open Targets locus-to-gene score (0–1), a prediction of the causal gene.'])]
    nts += missing_block([('Publication', 'none listed (consortium release)'), ('Lead variant', 'no rsID'),
                          ('β, Fine-mapping', 'not given'), ('Locus-to-gene genes', 'no prediction')])
    nts.append(abbrev(['CTEPH', 'GWAS', 'L2G', 'EAS', 'β']))
    return flat(wb, 'A2 GWAS loci', f'{len(g):,} credible sets, by condition and P', cols, rows, {}, nts, foot)


def a_shared(wb, D, foot):
    c = D['coloc'].assign(kind='colocalization', score=lambda d: num(d.clpp))
    o = D['overlap'].assign(kind='lead-variant overlap', score=lambda d: num(d.pip), method='lead-variant PIP',
                            l2g_gene='', eas='False')
    x = pd.concat([c, o], ignore_index=True)
    order = {g: i for i, g in enumerate(D['regions'].gene)}
    x['_r'] = x.cteph_gene.map(order)
    x = x.sort_values(['_r', 'kind', 'score'], ascending=[True, True, False])
    cols = [{'h': 'CTEPH region', 'w': 11}, {'h': 'Lead variant (CTEPH)', 'w': 13}, {'h': 'Evidence', 'w': 15, 'n': 'a'},
            {'h': 'Condition / trait', 'w': 30, 'wrap': True}, {'h': 'Study', 'w': 20},
            {'h': 'Trait as reported', 'w': 30, 'wrap': True}, {'h': 'Publication', 'w': 16},
            {'h': 'Ancestry', 'w': 18, 'wrap': True}, {'h': 'EAS', 'w': 6, 'a': 'center'},
            {'h': 'Lead variant (trait)', 'w': 13}, {'h': 'P', 'w': 12, 'a': 'right'},
            {'h': 'CLPP or PIP', 'w': 8, 'a': 'right', 'n': 'b'}, {'h': 'Top locus-to-gene gene', 'w': 10, 'n': 'c'}]
    rows = []
    for r in x.itertuples():
        pub = ' '.join(y for y in (r.first_author, str(r.pub_date)[:4]) if y) or r.project or None
        rows.append([[('i', r.cteph_gene)], r.cteph_locus, r.kind, r.condition or None, r.study_id, r.trait or None,
                     pub, ancestry_label(r.ancestry) or None, 'yes' if b(r.eas) else 'no', r.other_lead or None,
                     p_parts(r.p_mantissa, r.p_exponent), f'{r.score:.3f}' if pd.notna(r.score) else None,
                     [('i', r.l2g_gene)] if r.l2g_gene else ''])
    nts = [('a', ['Colocalization: a GWAS credible set of the trait colocalizes with the CTEPH credible set '
                  f'({D["choice"]["clpp"]["value_text"]}); lead-variant overlap: a GWAS credible set of the trait '
                  f'contains the CTEPH lead variant ({D["choice"]["overlap"]["value_text"]}), for the regions without a '
                  'CTEPH credible set. ' + basis(D, 'overlap')]),
           ('b', ['CLPP: eCAVIAR colocalization posterior probability of the two credible sets; PIP: the posterior '
                  'inclusion probability of the CTEPH lead variant in the trait\'s credible set.']),
           ('c', ['Of the CTEPH credible set, for colocalizations; the region itself is named by the gene the '
                  'publication reports. ' + basis(D, 'region_names')]),
           ('', [('b', 'Missing values')]),
           ('', ['   Empty cells: no CTEPH credible set (lead-variant overlap rows).'])]
    nts += missing_block([('Condition / trait', 'the study has no disease mapping'),
                          ('Publication, Ancestry, Lead variant (trait)', 'not given')])
    nts.append(abbrev(['CTEPH', 'GWAS', 'CLPP', 'PIP', 'EAS']))
    return flat(wb, 'A3 Shared signals', f'{len(x)} shared signals: {len(c)} colocalizations, {len(o)} lead-variant '
                'overlaps', cols, rows, {}, nts, foot)


def a_curated(wb, D, foot):
    m = D['mendelian'].sort_values(['condition', 'symbol', 'source']) if len(D['mendelian']) else D['mendelian']
    cols = [{'h': 'Condition', 'w': 28, 'wrap': True}, {'h': 'Gene', 'w': 10}, {'h': 'Source', 'w': 16},
            {'h': 'Score', 'w': 7, 'a': 'right', 'n': 'a'}, {'h': 'Classification', 'w': 16, 'wrap': True, 'n': 'b'},
            {'h': 'ClinVar clinical significance', 'w': 24, 'wrap': True},
            {'h': 'ClinVar review status', 'w': 26, 'wrap': True}, {'h': 'Specific diseases', 'w': 40, 'wrap': True},
            {'h': 'Records', 'w': 8, 'a': 'right'}, {'h': 'PMIDs', 'w': 24, 'wrap': True}]
    rows = []
    for x in m.itertuples():
        cv = x.source == 'ClinVar'
        rows.append([x.condition, [('i', x.symbol)], x.source, fmt2(x.source_score),
                     '' if cv else (x.classifications or None), (x.clinical_significance or None) if cv else '',
                     (getattr(x, 'review_status', '') or None) if cv else '', x.specific_diseases or None,
                     x.n_records, ' '.join(x.pmids.split()[:6]) + (' …' if len(x.pmids.split()) > 6 else '') or None])
    nts = [('a', ['Open Targets score of the source for the condition and its descendants (0–1); genes from 0.5 '
                  'upward are listed.']),
           ('b', ['ClinGen or Gene2Phenotype class, PanelApp rating (green, amber), Orphanet assessment, UniProt '
                  'confidence, as given.'])]
    nts += [('', [('b', 'Missing values')]),
            ('', ['   Empty cells: the field belongs to another source (ClinVar fields for ClinVar rows only).'])]
    nts += missing_block([('Classification', 'the source gives none (gene burden)'),
                          ('Specific diseases', 'none given'), ('PMIDs', 'none given')])
    nts.append(abbrev(['ClinGen', 'PanelApp', 'PMID']))
    return flat(wb, 'A4 Curated genes', f'{len(m)} condition × gene × source', cols, rows, {}, nts, foot, freeze=2)


def a_cteph_studies(wb, D, foot):
    s = D['cst'].assign(_y=lambda d: num(d.year)).sort_values(['_y', 'pmid'])
    cols = [{'h': 'PMID', 'w': 9}, {'h': 'Year', 'w': 6}, {'h': 'First author', 'w': 14},
            {'h': 'Title', 'w': 44, 'wrap': True}, {'h': 'Study type', 'w': 12},
            {'h': 'Design class', 'w': 20, 'wrap': True}, {'h': 'Design', 'w': 22, 'wrap': True},
            {'h': 'Place', 'w': 12, 'wrap': True, 'n': 'c'}, {'h': 'Population', 'w': 12, 'wrap': True},
            {'h': 'Japanese', 'w': 8, 'a': 'center', 'n': 'a'},
            {'h': 'East Asian', 'w': 8, 'a': 'center'}, {'h': 'Cases', 'w': 7, 'a': 'right'},
            {'h': 'Controls', 'w': 8, 'a': 'right'}, {'h': 'Control type', 'w': 20, 'wrap': True},
            {'h': 'Genes', 'w': 22, 'wrap': True, 'n': 'b'}, {'h': 'Variants', 'w': 20, 'wrap': True},
            {'h': 'Replication', 'w': 9}, {'h': 'Methods from', 'w': 9, 'a': 'center'},
            {'h': 'Publication flag', 'w': 14, 'wrap': True}]
    rows, links = [], {}
    for i, x in enumerate(s.itertuples()):
        rows.append([x.pmid, x.year, x.first_author, x.title, x.study_type or None, x.design_class_label, x.design,
                     x.place, x.population,
                     'yes' if x.jp else 'no', 'yes' if x.eas else 'no', x.n_cases or '', x.n_controls or '',
                     x.control_type or '', italic_list(x.gene_symbols.split('; ')) if x.gene_symbols else
                     (x.genes or ''), x.variants or '', x.replication, x.methods_source or None,
                     x.publication_flag or ''])
        links[(i, 0)] = pub_url(x.pmid)
    nts = [('a', ['Done in Japan (first author\'s institution), or a Japanese population stated in the abstract.']),
           ('c', ['Country of the first author\'s institution (PubMed affiliation; Europe PMC when PubMed has none).']),
           ('b', ['Official gene symbols checked against Open Targets; the record\'s own names where none was '
                  'found.'])]
    nts += [('', [('b', 'Missing values')]),
            ('', ['   Empty cells: not stated in the abstract, or no gene named.'])]
    nts += missing_block([('Study type', 'not classified'), ('Place', 'no affiliation in PubMed or Europe PMC')])
    nts.append(abbrev(['CTEPH', 'PMID', 'HLA']))
    return flat(wb, 'A5 CTEPH studies', f'{len(s)} studies, by year (publications left out by the curation are in '
                f'sheet {SN["curation"].split()[0]})', cols, rows, links, nts, foot)


def a_screening(wb, D, foot):
    s = D['csc'][D['csc'].stage != 'publication status']     # retracted / conference records are not listed
    s = s.assign(_y=lambda d: num(d.year), _i=lambda d: d.include != 'True').sort_values(['_i', '_y', 'pmid'])
    cols = [{'h': 'PMID', 'w': 9}, {'h': 'Year', 'w': 6}, {'h': 'Title', 'w': 70, 'wrap': True},
            {'h': 'Included', 'w': 8, 'a': 'center'}, {'h': 'Record type', 'w': 16},
            {'h': 'Decided by', 'w': 12, 'n': 'a'}, {'h': 'Reason', 'w': 60, 'wrap': True}]
    rows, links = [], {}
    for i, x in enumerate(s.itertuples()):
        rows.append([x.pmid, x.year or None, x.title or None, 'yes' if x.include == 'True' else 'no',
                     x.study_type or '', {'rule': 'pre-screen', 'model': 'language model'}.get(x.stage, x.stage),
                     x.reason.replace('pre-screen: ', '') or None])
        links[(i, 0)] = pub_url(x.pmid)
    nts = [('a', ['Pre-screen: a deterministic rule (publication type, no genetic-study term in the title or '
                  'abstract); language model: the remaining records, read in full abstract. Retracted articles and '
                  'conference material are not listed (the rule: sheet ' + SN['curation'].split()[0] + ').'])]
    nts += [('', [('b', 'Missing values')]), ('', ['   Empty cells: the pre-screen does not classify the record.'])]
    nts += missing_block([('Year, Title', 'not in the PubMed record'), ('Reason', 'none recorded')])
    nts.append(abbrev(['CTEPH', 'PMID']))
    return flat(wb, 'A6 Screening', f'{len(s)} PubMed records: included first', cols, rows, links, nts, foot)


# ── checks ───────────────────────────────────────────────────────────────────
def verify(path, stats):
    """Re-open the workbook: sheets; body rows; header letters against footnotes; a reason for every '—'; no
    missing-value artifacts or internal names, in cells or notes; outline levels; internal links; the XML Excel is
    strict about."""
    from openpyxl import load_workbook
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        wb = load_workbook(path, rich_text=True)
    assert wb.sheetnames == SH, wb.sheetnames
    for sheet in wb.worksheets:
        L = stats['layout'][sheet.title]
        h = L['header'] + 1
        head = {ch for c in range(1, sheet.max_column + 1) for ch in plain(sheet.cell(h, c).value, sup=True)
                if ch.isalpha() and ch.islower()}           # footnote letters, not exponents
        n = max(L['rows'], 1)
        body = range(h + 1, h + 1 + n)
        assert all(sheet.cell(r, 1).value is not None for r in body), (sheet.title, 'a body row without a label')
        assert sheet.cell(h + 1 + n, 1).value is None, (sheet.title, 'more body rows than recorded')
        letters, notes_text, has_dash = set(), '', False
        for r in range(h + 1 + n, sheet.max_row + 1):
            v = sheet.cell(r, 1).value
            if v is None:
                continue
            sup, txt = plain(v, sup=True), plain(v)
            if sup and txt.startswith(' ') and len(sup) == 1:
                letters.add(sup)
            notes_text += txt
        for row in sheet.iter_rows(min_row=h + 1, max_row=h + n):
            for c in row:
                if c.value is not None and plain(c.value) == MISSING:
                    has_dash = True
        assert head <= letters, (sheet.title, sorted(head), sorted(letters))
        if has_dash:
            assert '—,' in notes_text or '—:' in notes_text, (sheet.title, 'a — without its reason')
        prev = 0
        for r in body:
            d = sheet.row_dimensions[r]
            lvl = d.outlineLevel or 0
            assert lvl <= prev + 1, (sheet.title, r, 'outline level jumps')
            assert not d.hidden or lvl >= 1, (sheet.title, r, 'a hidden row outside a group')
            prev = lvl
        for row in sheet.iter_rows():
            for c in row:
                texts = [] if c.value is None else [plain(c.value) + plain(c.value, sup=True)]
                if c.comment is not None:
                    texts.append(c.comment.text)
                for s in texts:
                    assert s not in ('nan', 'None', 'NaN') and 'nan' not in s.split() and 'NaN' not in s.split(), \
                        (sheet.title, c.coordinate, s[:80])
                    bad = [t for t in FORBIDDEN_TOKENS if t in s]
                    assert not bad, (sheet.title, c.coordinate, bad, s[:80])
                if c.hyperlink is not None and c.hyperlink.location:
                    m = re.match(r"'?([^'!]+)'?!A(\d+)$", c.hyperlink.location)
                    assert m and m.group(1) in SH, (sheet.title, c.coordinate, c.hyperlink.location)
                    tgt = wb[m.group(1)]
                    assert tgt.cell(int(m.group(2)), 1).value is not None, (sheet.title, c.coordinate,
                                                                          c.hyperlink.location, 'empty target')
    z = zipfile.ZipFile(path)
    for name in z.namelist():
        if name.startswith('xl/worksheets/sheet'):
            x = z.read(name).decode()
            assert 'inlineStr' not in x, f'inline strings in {name}'
            for pane in re.findall(r'<pane [^>]*>', x):
                assert 'xSplit=' in pane or 'ySplit=' in pane, (name, pane)
    ss = z.read('xl/sharedStrings.xml').decode()
    bare = [t for p, t in re.findall(r'<t( xml:space="preserve")?>([^<]*)</t>', ss) if not p and t != t.strip()]
    assert not bare, bare[:3]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tables', required=True)
    ap.add_argument('--curated', required=True, help='curate.py outputs')
    ap.add_argument('--graph-schema', required=True, help='build_graph.py schema.json')
    ap.add_argument('--summaries', required=True)
    ap.add_argument('--changes', required=True)
    ap.add_argument('--changes-json', required=True)
    ap.add_argument('--knowledge-json', required=True, help='for the caveats')
    ap.add_argument('--themes', required=True, help='gene themes of the CTEPH literature (YAML)')
    ap.add_argument('--choices', required=True, help='every choice and its basis (YAML)')
    ap.add_argument('--run-params', required=True, help="the run's parameters (JSON)")
    ap.add_argument('--l2g-min', type=float, default=0.5)
    ap.add_argument('--eas-min', type=float, default=0.5)
    ap.add_argument('--out', required=True)
    ap.add_argument('--out-stats', default='workbook.json')
    a = ap.parse_args()
    T = Path(a.tables)
    cond = read_tsv(T / 'conditions.tsv')
    D = {'tidy': json.loads((T / 'tidy.json').read_text()), 'conditions': cond, 'relations': read_tsv(T / 'relations.tsv'),
         'coloc': read_tsv(T / 'colocalisation.tsv'), 'studies': read_tsv(T / 'gwas_studies.tsv'),
         'loci': read_tsv(T / 'gwas_loci.tsv'), 'mendelian': read_tsv(T / 'mendelian.tsv'),
         'gene_condition': read_tsv(T / 'gene_condition.tsv'), 'literature': read_tsv(T / 'literature.tsv'),
         'targets': read_tsv(T / 'cteph_targets.tsv'), 'cst': read_tsv(T / 'cteph_studies.tsv'),
         'cfi': read_tsv(T / 'cteph_findings.tsv'), 'cmf': read_tsv(T / 'cteph_methods.tsv'),
         'csc': read_tsv(T / 'cteph_screening.tsv'), 'summ': read_tsv(a.summaries),
         'changes': read_tsv(a.changes), 'changes_json': json.loads(Path(a.changes_json).read_text()),
         'caveats': json.loads(Path(a.knowledge_json).read_text())['caveats'], 'themes': load_themes(a.themes),
         'names': dict(zip(cond.node_id, cond.name)), 'l2g_min': a.l2g_min, 'eas_min': a.eas_min,
         'regions': read_tsv(T / 'cteph_regions.tsv'), 'ladder': read_tsv(T / 'cteph_genes.tsv'),
         'overlap': read_tsv(T / 'overlap.tsv'), 'pubs': read_tsv(Path(a.curated) / 'publications.tsv'),
         'clog': read_tsv(Path(a.curated) / 'curation_log.tsv'),
         'tmr': read_tsv(Path(a.curated) / 'textmining_review.tsv'),
         'cur': json.loads((Path(a.curated) / 'curation.json').read_text()),
         'graph': json.loads(Path(a.graph_schema).read_text())['graph']}
    for c in ('level', 'n_studies', 'associated', 'not_associated', 'identified', 'other', 'within_cteph',
              'eas_studies'):
        D['ladder'][c] = D['ladder'][c].astype(int)
    D['choices'], D['refs'] = load_choices(a.choices, json.loads(Path(a.run_params).read_text()), T)
    prepare(D)
    st = D['tidy']
    foot = f'CTEPH knowledge graph · Open Targets {st["release"]} · fetched {st["fetch_date"]}'
    # every question sheet is laid out before any is written, so that sheets can link to each other's rows
    D['gene_map_genes'] = set()
    gmap = spec_shared_genes(D)
    D['gene_map_genes'] = {r['key'][1] for r in gmap['rows'] if r['key']}
    specs = [spec_gene_map(D), spec_cteph_genes(D), spec_methods(D), spec_conditions(D), spec_regions(D), gmap,
             spec_condition_genetics(D), spec_mechanisms(D), spec_eas(D), spec_curation(D), spec_choices(D)]
    assert [s['name'] for s in specs] == [n for _k, n, _q in SHEETS[:-1]], 'sheets out of order'
    anchors = {}
    for s in specs:
        top, keys = layout(s)
        assert not set(keys) & set(anchors), set(keys) & set(anchors)
        anchors.update(keys)
    D['appendix_rows'] = {'A1 GWAS studies': len(D['studies']), 'A2 GWAS loci': len(D['loci']),
                          'A3 Shared signals': len(D['coloc']) + len(D['overlap']),
                          'A4 Curated genes': len(D['mendelian']),
                          'A5 CTEPH studies': len(D['cst']),
                          'A6 Screening': int((D['csc'].stage != 'publication status').sum())}
    wb = xlsxwriter.Workbook(a.out, {'nan_inf_to_errors': False})
    wb.set_properties({'title': f'CTEPH knowledge graph ({st["fetch_date"]})',
                       'subject': 'Conditions related to CTEPH and their genetics'})
    lay = {'Guide': s_guide(wb, D, foot)}
    for s in specs:
        lay[s['name']] = write_spec(wb, s, anchors, foot)
    lay[SN['changes']] = s_changes(wb, D, foot)
    for (name, _d), f in zip(APPENDIX, (a_studies, a_loci, a_shared, a_curated, a_cteph_studies, a_screening)):
        lay[name] = f(wb, D, foot)
    wb.close()
    stats = {'sheets': SH, 'rows': {k: v['rows'] for k, v in lay.items()}, 'layout': lay, 'release': st['release'],
             'fetch_date': st['fetch_date']}
    verify(a.out, stats)
    Path(a.out_stats).write_text(json.dumps(stats, indent=1))
    log(f'workbook: {", ".join(f"{k} {v}" for k, v in stats["rows"].items())}')


if __name__ == '__main__':
    main()
