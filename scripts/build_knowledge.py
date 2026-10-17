#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The knowledge base of one snapshot (Markdown + JSON), in the order
#           of the workbook's sheets: CTEPH's genes first, then its related
#           conditions, what they share, the context, and the basis of every
#           choice. Every section is built from the tables; the generated
#           summaries (summarise.py) are shown only when grounded, always under
#           a "Generated summary" heading and followed by the source records
#           they rest on. The figures of the snapshot (09.figures) are embedded.
#             1  Scope (and the questions, by part)
#             A  CTEPH genetics
#             2  CTEPH genes: evidence level x function; the published GWAS
#                regions; what the genetic studies report, by gene theme
#             3  Genetic studies of CTEPH: designs, places and methods
#             B  Related conditions
#             4  Conditions related to CTEPH, and why
#             5  The CTEPH GWAS regions the related conditions share
#             6  Genes shared by CTEPH and the related conditions
#             7  Genetics of the related conditions (and how their GWAS were done)
#             C  Context
#             8  Mechanisms (text-mined CTEPH targets by Reactome theme)
#             9  East Asian evidence
#             D  Provenance
#             10 Curation: what was left out or corrected, and why
#             11 Choices and their basis
#             12 Changes since the previous snapshot
#             13 Caveats
#             References
# Inputs  : the tidy tables, the curation outputs, the graph schema,
#           summaries, changes, the seed and its PubMed records, the gene
#           themes, the choices registry and the run's parameters
# Outputs : cteph_knowledge_graph.md, cteph_knowledge_graph.json
# ---------------------------------------------------------------------------
import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import (LADDER, NO_GENE, OTHER_GENES, assign_themes, gene_tallies,  # noqa: E402
                       is_east_asian, is_japanese, load_choices, load_json, load_themes, log, num, read_tsv,
                       study_population, tally_text)

CORE = 'EFO_0022954'
CATEGORY_ORDER = ['Precursor / spectrum', 'Risk factor', 'Comorbidity', 'Intermediate trait', 'PH group',
                  'Shared genetic signal']
QUESTIONS = ['Which genes have genetic evidence for CTEPH, how strong is it, and what do they do?',
             'Which genetic studies of CTEPH are there, and how were they done?',
             'Which conditions are related to CTEPH, and on what evidence?',
             'Which CTEPH GWAS regions do the related conditions share?',
             'Which genes do CTEPH and its related conditions share?',
             'What are the genetics of each related condition?',
             'Which mechanisms does the CTEPH literature describe?',
             'What evidence comes from East Asian and Japanese populations?',
             'What was left out or corrected, and why?',
             'On what basis was each item and value chosen?']
SHEETS = ['1 and 2', '3', '4', '5', '6', '7', '8', '9', '10', '11']     # workbook sheets of each question
# the four parts of the line of argument (sections 2-12; the workbook's sheets follow the same order)
PARTS = [('A. CTEPH genetics', 'sections 2-3'), ('B. Related conditions', 'sections 4-7'),
         ('C. Context', 'sections 8-9'), ('D. Provenance', 'sections 10-12')]
KIND = {'GWAS': 'G', 'GWAS (PE-adjusted)': 'G*', 'Mendelian': 'C', 'Literature': 'L'}
CODES = ['G', 'G*', 'C', 'L']
SUP = str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
CAVEATS = [
    'Relations to CTEPH are evidence of association (clinical studies, ontology, shared genetic signal), not of '
    'causation. Many classic risk factors were found against non-thromboembolic PH controls; against PE survivors '
    'most were not enriched (Dodson 2022, PMID 35000431).',
    'The CTEPH GWAS regions are those of the publication (Liley 2024, Table 1), named as there. Open Targets holds no '
    'summary statistics for this study and fine-maps only some regions, so colocalization is possible at those '
    'regions only.',
    'G* regions are genome-wide significant only when the genetic similarity with pulmonary embolism is leveraged '
    '(conditional false discovery rate V); part of their signal is that of PE, so they are marked apart from G.',
    'Colocalization is eCAVIAR on fine-mapped credible sets (CLPP) or coloc (H4). Lead-variant overlap is weaker: '
    'it takes the CTEPH lead variant as a stand-in for an unfine-mapped CTEPH signal, and its PIP in the other '
    'trait\'s credible set bounds the CLPP only if that variant is causal.',
    'A trait sharing only the ABO region mostly reflects the pleiotropy of ABO, so single-region relations are '
    'weaker than multi-region ones.',
    'Locus-to-gene (L2G) scores of the related conditions are predictions of the causal gene, not proof.',
    'Curated gene evidence comes from curated resources (ClinGen, PanelApp, Orphanet, ClinVar, UniProt, '
    'Gene2Phenotype, gene burden) for the condition or its descendants in the ontology; it is mostly rare-variant, '
    'Mendelian evidence.',
    'The literature studies were screened and extracted by a language model from PubMed abstracts and open-access '
    'Methods sections; every gene, variant, number and quote was checked against the source text, and items that '
    'failed were dropped. Studies behind paywalls contribute their abstract only. Retracted articles and '
    'conference material are left out (section 10).',
    'The place of a study is the country of its first author\'s institution (PubMed affiliation); a multinational '
    'study takes the place of its lead institution, and the population the abstract states is kept beside it.',
    'The evidence ladder counts studies, not independent samples: some studies share cases. A gene with several '
    'functions sits in the first gene theme that lists it.',
    'Generated summaries are syntheses of the listed source records; they cite only those records, but they are not '
    'a systematic review.',
    'Text-mined target links can map an abbreviation to the wrong gene. The links of this snapshot were reviewed and '
    'the wrong mappings found were removed (section 10); others may remain, so the matched words are listed with '
    'every link.',
]


def md_table(cols, rows):
    esc = lambda x: str(x).replace('|', '/').replace('\n', ' ')  # noqa: E731
    return ['| ' + ' | '.join(cols) + ' |', '|' + '|'.join('---' for _ in cols) + '|'] + \
        ['| ' + ' | '.join(esc(x) for x in r) + ' |' for r in rows]


def counts(series):
    return ', '.join(f'{k} {v}' for k, v in series.value_counts().items())


def p_text(p):
    m, e = f'{float(p):.1e}'.split('e')
    return f'{m} × 10{str(int(e)).translate(SUP)}'


def question(k):
    return [f'Question: {QUESTIONS[k]} (workbook sheet{"s" if " " in SHEETS[k] else ""} {SHEETS[k]})', '']


def plain(text):
    """A generated summary without its bold markup (the paragraph labels stay as plain text)."""
    return re.sub(r'\*\*(.+?)\*\*', r'\1', text)


def figure(name, alt):
    return [f'![{alt}](../09.figures/{name}.png)', '']


def by_condition(d, cond, col):
    alias = {i: c.node_id for c in cond.itertuples() for i in c.ot_ids.split()}
    d = d.assign(node_id=d.node_id.map(lambda i: alias.get(i, i)), v=num(d[col]))
    return d.groupby(['node_id', 'cteph_gene']).agg(n=('study_id', 'nunique'), best=('v', 'max'))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tables', required=True)
    ap.add_argument('--curated', required=True, help='curate.py outputs')
    ap.add_argument('--graph-schema', required=True, help='build_graph.py schema.json')
    ap.add_argument('--summaries', required=True)
    ap.add_argument('--changes', required=True)
    ap.add_argument('--changes-json', required=True)
    ap.add_argument('--seed', required=True)
    ap.add_argument('--seed-pubmed', required=True)
    ap.add_argument('--themes', required=True, help='gene themes of the CTEPH literature (YAML)')
    ap.add_argument('--choices', required=True, help='the choices registry (YAML)')
    ap.add_argument('--run-params', required=True, help='the parameters of this run (JSON)')
    ap.add_argument('--l2g-min', type=float, default=0.5)
    ap.add_argument('--top-genes', type=int, default=8, help='GWAS genes listed per related condition')
    ap.add_argument('--max-changes', type=int, default=40)
    ap.add_argument('--out-md', default='cteph_knowledge_graph.md')
    ap.add_argument('--out-json', default='cteph_knowledge_graph.json')
    a = ap.parse_args()
    T = Path(a.tables)
    st = json.loads((T / 'tidy.json').read_text())
    cond, rel, co = read_tsv(T / 'conditions.tsv'), read_tsv(T / 'relations.tsv'), read_tsv(T / 'colocalisation.tsv')
    ov, rg, lad = read_tsv(T / 'overlap.tsv'), read_tsv(T / 'cteph_regions.tsv'), read_tsv(T / 'cteph_genes.tsv')
    studies, loci, men = read_tsv(T / 'gwas_studies.tsv'), read_tsv(T / 'gwas_loci.tsv'), read_tsv(T / 'mendelian.tsv')
    gc = read_tsv(T / 'gene_condition.tsv')
    cst, cfi = read_tsv(T / 'cteph_studies.tsv'), read_tsv(T / 'cteph_findings.tsv')
    cmf = read_tsv(T / 'cteph_methods.tsv')
    Q = Path(a.curated)
    pubs, clog = read_tsv(Q / 'publications.tsv'), read_tsv(Q / 'curation_log.tsv')
    cur = json.loads((Q / 'curation.json').read_text())
    graph = json.loads(Path(a.graph_schema).read_text())['graph']
    summ = read_tsv(a.summaries)
    changes, chj = read_tsv(a.changes), json.loads(Path(a.changes_json).read_text())
    seed = yaml.safe_load(Path(a.seed).read_text())
    themes = load_themes(a.themes)
    choices, _refs = load_choices(a.choices, json.loads(Path(a.run_params).read_text()), T)
    choice = {it['key']: it for g in choices for it in g['items']}
    pub = {r['pmid']: r for r in load_json(a.seed_pubmed)['records']}
    cite = lambda p: (f'{pub[p]["first_author"].split()[0]} {pub[p]["year"]}, PMID {p}'  # noqa: E731
                      if p in pub else f'PMID {p}')
    names = dict(zip(cond.node_id, cond.name))
    lad['level'] = lad.level.astype(int)
    cond['nreg'] = cond.n_shared_regions.astype(int)
    cst['population'] = [study_population(r.country, r.ancestry) for r in cst.itertuples()]
    cst['place'] = cst.place.replace('', 'not known')
    cst['jp'] = [is_japanese(r) for r in cst.itertuples()]
    cst['eas'] = [is_east_asian(r) for r in cst.itertuples()]
    short = {t['name']: t['short'] for t in themes}
    short[OTHER_GENES] = OTHER_GENES
    theme_cols = [t['name'] for t in themes] + [OTHER_GENES]
    K = {'release': st['release'], 'fetch_date': st['fetch_date'], 'questions': QUESTIONS,
         'counts': {k: v for k, v in st.items() if isinstance(v, (int, dict))}}
    md = [f'# CTEPH knowledge graph (Open Targets {st["release"]}, PubMed; fetched {st["fetch_date"]})', '',
          'Generated by the `cteph_knowledge_graph` pipeline; do not edit by hand. Sections are built from the '
          'snapshot tables, in the order of the workbook: CTEPH\'s genes, its related conditions, what they share, '
          'the context, then the basis of every choice. Text under a "Generated summary" heading was written by a '
          'language model from the records listed with it and checked against those records. Figures: '
          '`../09.figures/` (catalog in its README); the graph for programs: `../05.graph/` (schema in its '
          'README); the viewer: `../cteph_knowledge_graph.html`.', '']

    # 1 scope ───────────────────────────────────────────────────────────────
    bs = st['conditions_by_source']
    bylv = lad.level.value_counts()
    md += ['## 1. Scope', '',
           f'- CTEPH genes: {len(lad)} genes studied genetically in CTEPH ({int((lad.level <= 5).sum())} with a '
           'positive result), placed on an evidence ladder (levels 1–7: ' +
           ', '.join(f'{int(bylv.get(k, 0))}' for k in sorted(LADDER)) + ') and in a function (gene theme).',
           f'- CTEPH GWAS regions: the {len(rg)} genome-wide significant regions of the publication (Liley 2024); '
           f'Open Targets holds {int(rg.variant_id.ne("").sum())} lead variants and a CTEPH credible set at '
           f'{int(rg.cteph_study_locus_id.ne("").sum())}.',
           f'- Related conditions: {st["conditions"]} ({bs["clinical"]} by clinical evidence, PubMed-verified; '
           f'{bs["ontology"]} by the disease ontology; {bs["colocalization"]} by colocalization and '
           f'{bs["lead-variant overlap"]} by lead-variant overlap with a CTEPH GWAS region; a condition can have '
           'several).',
           f'- Genetics of the related conditions: {st["gwas_studies"]} GWAS studies ({st["eas_studies"]} East Asian) '
           f'with {st["credible_sets"]:,} credible sets; {st["mendelian_links"]} curated condition-gene links.',
           f'- CTEPH literature: {st["cteph_genetic_records"]} PubMed records screened, {st["cteph_genetic_studies"]} '
           f'genetic studies of CTEPH ({st["cteph_findings"]} findings, {st["methods_fields"]} methods fields), '
           f'after the curation left out {cur["studies_excluded"]} (section 10).',
           f'- Curation: {cur["publications"]:,} publications checked in PubMed, {len(cur["excluded"])} left out '
           f'(retracted or conference material); {cur["textmining_dropped"]} text-mined links dropped; every '
           'finding classified by what it is about (section 10).',
           f'- Choices: {sum(len(g["items"]) for g in choices)}, each with its value in this run and its basis '
           '(section 11).',
           f'- Graph: {graph["n_nodes"]:,} nodes, {graph["n_edges"]:,} edges of '
           f'{len(graph["counts"]["nodes"])} node and {len(graph["counts"]["edges"])} edge types (`05.graph`).', '',
           'Sections 2 to 11 answer these questions, in four parts and in the order of the workbook:', '']
    for name, secs in PARTS:
        md += [f'- {name} ({secs})']
    md += [''] + [f'{k}. {q}' for k, q in enumerate(QUESTIONS, 1)] + ['']
    K['graph'] = graph

    # 2 CTEPH genes ─────────────────────────────────────────────────────────
    md += ['## 2. CTEPH genes: strength and function of the evidence', ''] + question(0)
    md += ['Each gene sits at its strongest level of evidence (rows) and in the column of its function (gene '
           'theme). ' + choice['ladder']['basis'] + ' ' + choice['about']['basis'] + ' ' + choice['themes']['basis'],
           '']
    md += figure('evidence_grid', 'Evidence level x function of the CTEPH genes')
    rows = []
    for lv in sorted(LADDER):
        sub = lad[lad.level == lv]
        parts = [f'{short[t]}: ' + ', '.join(sub[sub.theme == t].symbol) for t in theme_cols if (sub.theme == t).any()]
        code, text = LADDER[lv]
        rows.append([lv, code or '—', text, len(sub), '; '.join(parts) or '—'])
    md += md_table(['Level', 'Code', 'Evidence', 'Genes', 'Genes by function'], rows) + ['']
    K['ladder'] = lad.to_dict('records')

    md += ['### 2.1 The published GWAS regions', '',
           'Liley et al. 2024 (PMID 38470220), Table 1: three tiers of significance on the meta-analysis P value (G), '
           'and four further regions significant only in the co-analysis with pulmonary embolism (G*, V < 5 × 10⁻⁸). '
           + choice['region_names']['basis'], '']
    l2g = {}
    for r in loci[loci.node_id == CORE].itertuples():
        l2g[r.rsid] = ', '.join(f'{getattr(r, f"l2g{k}_symbol")} {float(getattr(r, f"l2g{k}_score")):.2f}'
                                for k in (1, 2, 3) if getattr(r, f'l2g{k}_symbol'))
    rows = [[r.gene, r.rsid, f'chr{r.chrom}:{r.pos_grch37}', r.tier, 'G' if r.evidence == 'meta' else 'G*',
             p_text(r.p), p_text(r.p_pe), p_text(r.v), r['or'], r.maf,
             'credible set' if r.cteph_study_locus_id else 'lead variant only' if r.variant_id else 'not held',
             l2g.get(r.rsid, '—')] for _i, r in rg.iterrows()]
    md += md_table(['Region', 'Lead variant', 'Position (GRCh37)', 'Tier', 'Code', 'P', 'P in PE', 'V', 'OR', 'MAF',
                    'Open Targets', 'Open Targets locus-to-gene'], rows) + ['']
    K['regions'] = rg.to_dict('records')

    fx = assign_themes(cfi, themes).merge(cst[['pmid', 'year', 'first_author', 'eas']], on='pmid', how='left')
    gt = gene_tallies(fx)
    rep = gt[(gt.associated >= 2) & (gt.theme != NO_GENE)].sort_values('associated', ascending=False)
    md += ['### 2.2 What the genetic studies of CTEPH report, by gene theme', '',
           f'{len(cst)} studies, {len(cfi)} findings on {gt.label.nunique()} genes or gene sets. Associated in two or '
           'more studies: ' + ', '.join(f'{r.label} ({r.associated})' for r in rep.itertuples()) + '.', '']
    sf = summ[(summ.kind == 'findings') & (summ.grounded == 'True')]
    if len(sf):
        md += ['#### Generated summary of the findings', '', plain(sf.iloc[0].summary), '']
    desc = {t['name']: t['description'] for t in themes}
    level = dict(zip(lad.symbol, lad.level))
    K['themes'] = []
    for th in theme_cols + [NO_GENE]:
        g = gt[gt.theme == th].sort_values(['associated', 'n_studies', 'label'], ascending=[False, False, True])
        if not len(g):
            continue
        f = fx[fx.theme == th]
        md += [f'#### {th}', '', (desc.get(th, '') + ' ').lstrip() + f'{len(g)} genes, {f.pmid.nunique()} studies, '
               f'{len(f)} findings.', '']
        rows = []
        for r in g.itertuples():
            ff = f[f.label == r.label].assign(_y=lambda d: num(d.year)).sort_values('_y')
            refs = '; '.join(dict.fromkeys(
                f'{x.first_author} {x.year} ({x.variant_name + ", " if x.variant_name else ""}{x.result}'
                f'{"" if x.association_with == "CTEPH" else " with " + x.association_with}; PMID {x.pmid})'
                for x in ff.itertuples()))
            rows.append([r.label, level.get(r.label, '—'), tally_text(r), refs])
        md += md_table(['Gene', 'Level', 'Studies', 'Findings (variant, result)'], rows) + ['']
        K['themes'].append({'theme': th, 'genes': g.to_dict('records')})
    md += ['### 2.3 Curated facts', ''] + [f'- {f["topic"]}: {f["statement"]} ({cite(str(f["pmid"]))})'
                                           for f in seed.get('cteph_facts') or []] + ['']
    K['cteph_genetics'] = {'regions': len(rg), 'literature_studies': len(cst), 'findings': len(cfi),
                           'facts': seed.get('cteph_facts') or []}

    # 3 genetic studies of CTEPH ──────────────────────────────────────────
    sm = summ[(summ.kind == 'methods') & (summ.grounded == 'True')]
    md += ['## 3. Genetic studies of CTEPH: designs, places and methods', ''] + question(1)
    md += [choice['design_class']['basis'] + ' ' + choice['place']['basis'], '']
    md += figure('genetic_studies', 'The genetic studies of CTEPH by year, design class and place')
    if len(cst):
        md += ['- Design classes (one per study): ' + counts(cst.design_class_label) + '.',
               '- Designs (a study can have several): ' + counts(cst.design.str.split('; ').explode()) + '.',
               '- Places (country of the first author\'s institution): ' + counts(cst.place) + '.',
               '- Populations as the abstracts state them: ' + counts(cst.population) + '.',
               f'- Cases per study: median {num(cst.n_cases).median():.0f}, range {num(cst.n_cases).min():.0f}–'
               f'{num(cst.n_cases).max():.0f} (where stated).',
               f'- Methods fields ({len(cmf)}: platform, sequencing, quality control, statistics, covariates, '
               f'multiple testing): from the open-access Methods section for {cst.methods_source.eq("full text").sum()} '
               f'studies, from the abstract for {cst.methods_source.eq("abstract").sum()} (no full text).', '']
        rows = [[s.year, s.first_author, s.design_class_label, s.place, s.population, s.n_cases or '—',
                 s.n_controls or '—', s.methods_source or '—', s.pmid]
                for s in cst.assign(_y=num(cst.year)).sort_values('_y').itertuples()]
        md += md_table(['Year', 'First author', 'Design class', 'Place', 'Population', 'Cases', 'Controls',
                        'Methods from', 'PMID'], rows) + ['']
    if len(sm):
        md += ['#### Generated summary of the methods', '', plain(sm.iloc[0].summary), '']

    # 4 related conditions ──────────────────────────────────────────────────
    md += ['## 4. Conditions related to CTEPH, and why', ''] + question(2)
    md += [choice['categories']['basis'] + ' ' + choice['clinical_rule']['basis'], '']
    md += figure('condition_map', 'CTEPH and its related conditions')
    K['conditions'] = []
    for cat in CATEGORY_ORDER:
        sub = cond[cond.category == cat]
        if not len(sub):
            continue
        md += [f'### {cat}', '']
        if cat == 'Shared genetic signal':
            multi, single = sub[sub.nreg >= 2], sub[sub.nreg < 2]
            md += [f'Traits sharing a GWAS signal with a CTEPH region and with no clinical or ontology relation: '
                   f'{len(sub)}. ' + choice['loci_grading']['basis'], '',
                   f'Sharing two or more CTEPH regions ({len(multi)}):', '']
            md += md_table(['Trait', 'CTEPH regions shared', 'Regions'],
                           [[r.name, ', '.join(r.shared_regions.split()), r.nreg]
                            for r in multi.sort_values(['nreg', 'name'], ascending=[False, True]).itertuples()])
            md += ['', f'Sharing one CTEPH region only ({len(single)}), by region:', '']
            for reg, s in single.groupby('shared_regions'):
                md += [f'- {reg} ({len(s)}): ' + ', '.join(sorted(s.name, key=str.lower)) + '.']
            md += ['']
            continue
        rows = []
        for c in sub.itertuples():
            ev = rel[rel.node_id == c.node_id]
            cl = ev[ev.evidence_type == 'clinical']
            other = ev[ev.evidence_type != 'clinical']
            items = [f'{r.detail} ({cite(r.pmid)}{"; East Asian" if r.east_asian == "True" else ""})'
                     for r in cl.itertuples()]
            items += [r.detail for r in other.itertuples()]
            K['conditions'].append({'condition': c.name, 'id': c.node_id, 'category': cat, 'sources': c.sources,
                                    'shared_regions': c.shared_regions, 'evidence': items})
            rows.append([c.name, c.sources.replace(' ', ', '), ', '.join(c.shared_regions.split()) or '—',
                         '<br>'.join(items)])
        md += md_table(['Condition', 'Evidence types', 'CTEPH regions shared', 'Evidence'], rows) + ['']

    # 5 shared regions ──────────────────────────────────────────────────────
    ca, oa = by_condition(co, cond, 'clpp'), by_condition(ov, cond, 'pip')
    md += ['## 5. The CTEPH GWAS regions the related conditions share', ''] + question(3)
    md += [f'Colocalization: {choice["clpp"]["value_text"]} or {choice["h4"]["value_text"]}, possible only at a region '
           'with a CTEPH credible set in Open Targets. ' + choice['clpp']['basis'],
           '',
           f'Lead-variant overlap (weaker): {choice["overlap"]["value_text"]}, used only where Open Targets holds the '
           'CTEPH lead variant without a credible set. ' + choice['overlap']['basis'], '']
    md += figure('shared_regions', 'Related conditions x CTEPH GWAS regions')
    rows = []
    for r in rg.itertuples():
        c_n = sorted({n for n, g in ca.index if g == r.gene})
        o_n = sorted({n for n, g in oa.index if g == r.gene} - set(c_n))
        how = ('colocalization' if r.cteph_study_locus_id else 'lead-variant overlap' if r.variant_id
               else 'none (lead variant not in Open Targets)')
        n_sig = ca[ca.index.get_level_values(1) == r.gene].n.sum() + oa[oa.index.get_level_values(1) == r.gene].n.sum()
        rows.append([r.gene, 'G' if r.evidence == 'meta' else 'G*', how, len(c_n) + len(o_n), int(n_sig)])
    md += md_table(['Region', 'Code', 'How a condition can share it', 'Conditions sharing it', 'GWAS signals'], rows)
    md += ['', 'Conditions sharing two or more regions, or clinically related and sharing one (a ° marks a region '
           'shared by lead-variant overlap):', '']
    sh = cond[(cond.nreg >= 2) | ((cond.category != 'Shared genetic signal') & (cond.nreg >= 1))].copy()
    sh['_o'] = sh.category.map(CATEGORY_ORDER.index)
    rows = []
    for c in sh.sort_values(['_o', 'nreg', 'name'], ascending=[True, False, True]).itertuples():
        regs = [g + ('' if (c.node_id, g) in ca.index else '°') for g in rg.gene
                if (c.node_id, g) in ca.index or (c.node_id, g) in oa.index]
        rows.append([c.name, c.category, ', '.join(regs), c.nreg])
    md += md_table(['Condition', 'Category', 'CTEPH regions', 'Regions'], rows) + ['']
    K['shared_regions'] = rows

    # 6 shared genes ───────────────────────────────────────────────────────
    gc['code'] = gc.evidence.map(KIND)
    cell = gc.groupby(['symbol', 'node_id']).code.apply(lambda x: ''.join(sorted(set(x), key=CODES.index)))
    cteph_code = cell.xs(CORE, level='node_id') if CORE in cell.index.get_level_values(1) else pd.Series(dtype=str)
    relc = cell[(cell.index.get_level_values(1) != CORE) & cell.index.get_level_values(0).isin(set(lad.symbol))]
    nrel = relc.groupby(level=0).size()
    md += ['## 6. Genes shared by CTEPH and the related conditions', ''] + question(4)
    md += ['Codes: G, GWAS (for CTEPH, the gene the publication names for a genome-wide significant region; for a '
           f'related condition, the top locus-to-gene gene, score ≥ {a.l2g_min}, of one of its GWAS credible sets); '
           'G*, CTEPH only, a region significant only after leveraging PE; C, curated gene–disease evidence; L, '
           'CTEPH only, a genetic study in the literature; –, CTEPH only, studied without a positive result (levels '
           '6–7). ' + choice['g_codes']['basis'], '']
    md += figure('gene_map', 'Genes studied in CTEPH x the related conditions')
    only = lad[~lad.symbol.isin(nrel.index)]
    md += [f'{len(lad)} genes were studied genetically in CTEPH ({int((lad.level <= 5).sum())} with a positive result); '
           f'{len(nrel)} of them also have GWAS or curated evidence in at least one related condition. Traits related '
           'by a shared GWAS signal only have no gene data by design (section 5 shows their shared regions).', '']
    rows = []
    for x in lad[lad.symbol.isin(nrel.index)].assign(n=lambda d: d.symbol.map(nrel)).sort_values(
            ['level', 'n', 'symbol'], ascending=[True, False, True]).itertuples():
        cc = relc.xs(x.symbol, level=0)
        rows.append([x.symbol, x.level, cteph_code.get(x.symbol, '–'), int(x.n),
                     '; '.join(f'{names.get(n, n)} {k}' for n, k in sorted(cc.items(), key=lambda y: names.get(
                         y[0], y[0]).lower()))])
    md += md_table(['Gene', 'Level', 'CTEPH', 'Related conditions', 'Which (code)'], rows) + ['']
    md += [f'Specific to CTEPH in this map ({len(only)}: no GWAS or curated evidence in any related condition), by '
           'level: ' + '; '.join(f'{lv} ' + ', '.join(only[only.level == lv].symbol)
                                 for lv in sorted(only.level.unique())) + '.', '']
    gcs = gc.groupby('symbol').node_id.apply(lambda x: len(set(x) - {CORE}))
    other = gcs[~gcs.index.isin(set(lad.symbol)) & (gcs >= 3)].sort_values(ascending=False)
    md += [f'Not studied in CTEPH, {len(other)} genes have genetic evidence in three or more related conditions; '
           'the most widely shared: ' + ', '.join(f'{g} ({n})' for g, n in other.head(15).items()) + '.', '']
    K['shared_genes'] = rows

    # 7 genetics of related conditions ──────────────────────────────────────
    md += ['## 7. Genetics of the related conditions', ''] + question(5)
    md += [choice['genetics_mode']['basis'] + ' ' + choice['curated']['basis'], '']
    rows = []
    loci['l2g1_score_n'] = num(loci.l2g1_score)
    for c in cond[cond.genetics == 'full'].sort_values(['category', 'name']).itertuples():
        L = loci[(loci.node_id == c.node_id) & (loci.l2g1_score_n >= a.l2g_min)]
        top = L.groupby('l2g1_symbol').study_id.nunique().sort_values(ascending=False).head(a.top_genes)
        M = men[men.node_id == c.node_id] if len(men) else men
        mg = M.groupby('symbol').source.apply(lambda x: '/'.join(sorted(set(x)))) if len(M) else pd.Series(dtype=str)
        rows.append([c.name, c.category, f'{c.n_gwas_studies} ({c.n_eas_studies} EAS)', c.n_credible_sets,
                     ', '.join(f'{g} ({n})' for g, n in top.items()) or '—',
                     ', '.join(f'{g} [{s}]' for g, s in list(mg.items())[:12]) or '—'])
    md += md_table(['Condition', 'Category', 'GWAS studies', 'Credible sets', 'Top GWAS genes (studies)',
                    'Curated genes [source]'], rows) + ['']

    gw = studies.copy()
    gw['n_samples'] = num(gw.n_samples)
    md += ['### How the GWAS were done (Open Targets: CTEPH and the related conditions)', '',
           f'- {len(gw)} studies; sample size median {gw.n_samples.median():,.0f} (range {gw.n_samples.min():,.0f}–'
           f'{gw.n_samples.max():,.0f}); with summary statistics {gw.summary_statistics.eq("True").sum()}.',
           '- Fine-mapping of the credible sets: ' + counts(loci.finemapping) + '.',
           '- Largest ancestry group of each study: ' + counts(gw.ancestry.str.split('; ').str[0].str.split().str[0])
           + '.', '']


    # 8 mechanisms ──────────────────────────────────────────────────────────
    md += ['## 8. Mechanisms (text-mined CTEPH literature, by Reactome theme)', ''] + question(6)
    md += [choice['reactome']['basis'], '']
    mech = summ[summ.kind == 'mechanism'].copy()
    mech['n'] = num(mech.n_targets)
    for r in mech.sort_values('n', ascending=False).itertuples():
        md += [f'### {r.title} ({int(r.n)} targets)', '']
        if r.grounded == 'True':
            md += ['Generated summary: ' + plain(r.summary), '']
        elif r.problems.startswith('not generated'):
            md += ['(No generated summary in this build: the source sentences changed with the curation, and this '
                   'build calls no language model.)', '']
        else:
            md += [f'(Summary withheld: {r.problems}.)', '']
        md += [f'Source targets: {r.genes}. Source PMIDs: {r.pmids}.', '']
    K['mechanisms'] = summ[summ.kind == 'mechanism'].to_dict('records')

    # 9 East Asian evidence ─────────────────────────────────────────────────
    eas = studies[studies.eas == 'True']
    ea_rel = rel[rel.east_asian == 'True']
    ea_lit = cst[cst.eas]
    md += ['## 9. East Asian evidence', ''] + question(7)
    md += [f'- Clinical evidence from East Asian populations ({len(ea_rel)}): ' +
           '; '.join(f'{names[r.node_id]}: {r.detail} ({cite(r.pmid)})' for r in ea_rel.itertuples()) + '.',
           f'- Genetic studies of CTEPH done in East Asia ({len(ea_lit)}; in Japan {int(ea_lit.jp.sum())}): ' +
           '; '.join(f'{r.first_author} {r.year} ({r.place}; {r.design_class_label}; '
                     f'{r.gene_symbols or r.genes or "no gene"}; PMID {r.pmid})'
                     for r in ea_lit.sort_values('year').itertuples()) + '.',
           f'- GWAS of the related conditions with a mostly East Asian sample ({len(eas)}): ' +
           '; '.join(f'{r.condition}: {r.study_id} ({r.first_author or r.project} {r.year}, n = {r.n_samples})'
                     for r in eas.itertuples()) + '.', '']

    # 10 curation ───────────────────────────────────────────────────────────
    md += ['## 10. Curation: what was left out or corrected, and why', ''] + question(8)
    md += [f'Every publication the graph uses ({cur["publications"]:,}) was looked up in PubMed on '
           f'{cur["checked_at"][:10]}; the rules are in `config/curation.yaml`, and every rule applied is in '
           '`03.curated/curation_log.tsv`. No rule calls a language model.', '', '### Publications left out', '',
           choice['exclusion']['basis'] + f' {len(cur["excluded"])} publications meet the rule; they are part of no '
           'table, figure or count of this snapshot and are not listed.', '']
    fl = pubs[pubs.status == 'flagged']
    md += [f'Flagged, kept: {len(fl)} ({counts(fl.reason)}).', '', '### Findings that are not about CTEPH itself',
           '', choice['about']['basis'], '']
    fa = cfi[cfi.association_with != 'CTEPH'].merge(cst[['pmid', 'first_author', 'year']], on='pmid', how='left')
    md += md_table(['Study', 'Gene', 'Result', 'About', 'Note', 'PMID'],
                   [[f'{x.first_author} {x.year}', x.symbol or x.gene or '—', x.result, x.association_with,
                     x.association_note or '—', x.pmid] for x in fa.itertuples()]) + ['']
    renamed = int(((clog.step == 'findings') & (clog.rule == 'variant name')).sum())
    md += ['### Variants named', '', choice['variants']['basis'] + f' Renamed in this snapshot: {renamed} findings '
           '(e.g. factor V Leiden for "FVL", "1691G>A" and "factor V leidan"; HLA-B*52:01 for "B(*)5201").', '',
           '### Text-mined links mapped to the wrong gene', '', choice['textmining_review']['basis'], '']
    tm = clog[(clog.step == 'text mining') & (clog.rule == 'link dropped')]
    wrong = tm[tm.after.str.startswith('mapped to the wrong gene')]
    md += md_table(['Gene', 'Publications', 'What the word means'],
                   [[x.item, x.before.split()[0], x.after.split(': ', 1)[-1]] for x in wrong.itertuples()]) + ['']
    md += ['### Place and design class of the genetic studies', '',
           f'{choice["place"]["basis"]} Places: ' + counts(cst.place) + '. ' + choice['design_class']['basis'], '']
    K['curation'] = cur

    # 11 choices ────────────────────────────────────────────────────────────
    by = pd.Series([it['basis_type'] for g in choices for it in g['items']]).value_counts()
    md += ['## 11. Choices and their basis', ''] + question(9)
    md += [f'{int(by.sum())} choices decide what enters the knowledge graph: ' +
           ', '.join(f'{int(v)} {k}' for k, v in by.items()) + '. Source: a publication or the documentation of a '
           'data source, checked against its text. Data check: recomputed on this snapshot. Design: a decision of '
           'this project, with its reason. User: decided with the user, on the date given.', '']
    used = {}
    for g in choices:
        md += [f'### {g["name"]}', '']
        rows = []
        for it in g['items']:
            chk = (' Check ' + ('passed' if it['check_ok'] else 'FAILED') + ': ' + it['check_text']) \
                if it['check_text'] else ''
            rows.append([it['item'], it['value_text'], it['basis_type'], it['basis'] + chk,
                         '; '.join(r[0] for r in it['refs']) or '—'])
            used.update({r[0]: r[2] for r in it['refs']})
        md += md_table(['Choice', 'Value in this run', 'Basis type', 'Basis', 'References'], rows) + ['']
    K['choices'] = [{'group': g['name'], 'items': [{k: it[k] for k in ('key', 'item', 'value_text', 'basis_type',
                                                                         'basis', 'check_text', 'check_ok')}
                                                    for it in g['items']]} for g in choices]

    # 12 changes ────────────────────────────────────────────────────────────
    md += ['## 12. Changes since the previous snapshot', '']
    if chj['previous'] is None:
        md += ['First snapshot: nothing to compare with.', '']
    else:
        md += [f'Compared with {chj["previous"]} (release {chj["previous_release"]}): ' +
               ('; '.join(f'{k.replace("|", " ")} {v}' for k, v in chj['counts'].items()) or 'no change') + '.', '']
        if chj.get('not_compared'):
            md += ['Not compared (tables the earlier snapshot does not have): ' + ', '.join(chj['not_compared']) + '.',
                   '']
        if len(changes):
            md += md_table(['Category', 'Change', 'Item', 'Detail'],
                           [[r.category, r.change, r.item, r.detail] for r in changes.head(a.max_changes).itertuples()])
            md += ['']
    K['changes'] = chj

    md += ['## 13. Caveats', ''] + [f'- {c}' for c in CAVEATS] + ['']
    K['caveats'] = CAVEATS
    md += ['## References', ''] + [f'- {full}' for _s, full in sorted(used.items())] + ['']
    K['references'] = used
    Path(a.out_md).write_text('\n'.join(md))
    Path(a.out_json).write_text(json.dumps(K, indent=1, ensure_ascii=False, default=str))
    log(f'knowledge: {len(md)} lines')


if __name__ == '__main__':
    main()
