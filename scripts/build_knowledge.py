#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The knowledge base of one snapshot of one disease (Markdown +
#           JSON), in the order of the workbook's sheets: the disease's genes
#           first, then its related conditions, what they share, the context,
#           and the basis of every choice. Every section is built from the
#           tables; the generated summaries (summarise.py) are shown only when
#           grounded, always under a "Generated summary" heading and followed
#           by the source records they rest on. The figures of the snapshot
#           (09.figures) are embedded. Names come from the disease profile;
#           text specific to a disease from its text/knowledge.yaml, and a
#           generic sentence stands in when it has none.
#             1  Scope (and the questions, by part)
#             A  Genetics of the disease
#             2  Its genes: evidence level x function; its GWAS regions; what
#                its genetic studies report, by gene theme
#             3  Its genetic studies: designs, places and methods
#             B  Related conditions
#             4  Conditions related to the disease, and why
#             5  The disease's GWAS regions the related conditions share
#             6  Genes shared by the disease and the related conditions
#             7  Genetics of the related conditions (and how their GWAS were done)
#             C  Context
#             8  Mechanisms (text-mined targets of the disease, by Reactome theme)
#             9  East Asian evidence
#             D  Provenance
#             10 Curation: what was left out or corrected, and why; the
#                language-model readers and their agreement
#             11 Choices and their basis
#             12 Changes since the previous snapshot
#             13 Caveats
#             References
# Inputs  : the disease profile, the tidy tables (with the ladder and the gene
#           themes of the snapshot), the curation outputs, the graph schema,
#           summaries, changes, the PubMed records of the seed and the run's
#           parameters
# Outputs : <key>_knowledge_graph.md, <key>_knowledge_graph.json
# ---------------------------------------------------------------------------
import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import (MISSING, a_, NO_ASSOCIATION, NO_GENE, OTHER_GENES, SHARED_SIGNAL, assign_themes,  # noqa: E402
                       gene_tallies,
                       is_east_asian, is_japanese, load_choices, load_json, load_ladder, load_profile, load_themes,
                       log, not_in_gene_map, not_regions_of, not_regions_sentence, num, p_unicode, read_tsv, study_population, tally_text)

QUESTIONS = ['Which genes have genetic evidence for {short}, how strong is it, and what do they do?',
             'Which genetic studies of {short} are there, and how were they done?',
             'Which conditions are related to {short}, and on what evidence?',
             'Which {short} GWAS regions do the related conditions share?',
             'Which genes do {short} and its related conditions share?',
             'What are the genetics of each related condition?',
             'Which mechanisms does the {short} literature describe?',
             'What evidence comes from East Asian and Japanese populations?',
             'What was left out or corrected, and why?',
             'On what basis was each item and value chosen?']
SHEETS = ['1 and 2', '3', '4', '5', '6', '7', '8', '9', '10', '11']     # workbook sheets of each question
# the four parts of the line of argument (sections 2-12; the workbook's sheets follow the same order)
PARTS = [('A. {short} genetics', 'sections 2-3'), ('B. Related conditions', 'sections 4-7'),
         ('C. Context', 'sections 8-9'), ('D. Provenance', 'sections 10-12')]
# statistics of a GWAS region (disease_regions.tsv), in the order shown, with the header used when the profile's
# published table (gwas_regions.yaml, columns) names none
REGION_STATS = [('p', 'P'), ('p_pe', 'P (secondary analysis)'), ('v', 'V'), ('or', 'OR'), ('maf', 'MAF')]


def md_table(cols, rows):
    esc = lambda x: str(x).replace('|', '/').replace('\n', ' ')  # noqa: E731
    return ['| ' + ' | '.join(cols) + ' |', '|' + '|'.join('---' for _ in cols) + '|'] + \
        ['| ' + ' | '.join(esc(x) for x in r) + ' |' for r in rows]


def counts(series):
    return ', '.join(f'{k} {v}' for k, v in series.value_counts().items())


def plural(n, word, words=None):
    return f'{n} {word if n == 1 else words or word + "s"}'


def reader_count(n):
    return {1: 'one language-model reader', 2: 'two language-model readers'}.get(n, f'{n} language-model readers')


def p_text(p):
    """A P value as 6.8 × 10⁻²⁷; the em dash when the table gives none."""
    if pd.isna(num(pd.Series([p])).iloc[0]):
        return MISSING
    m, e = f'{float(p):.1e}'.split('e')
    return p_unicode(float(m), int(e))


def question(qs, k):
    return [f'Question: {qs[k]} (workbook sheet{"s" if " " in SHEETS[k] else ""} {SHEETS[k]})', '']


def plain(text):
    """A generated summary without its bold markup (the paragraph labels stay as plain text)."""
    return re.sub(r'\*\*(.+?)\*\*', r'\1', text)


def figure(name, alt):
    return [f'![{alt}](../09.figures/{name}.png)', '']


def sentence(*parts):
    """Sentences joined by a space, the empty ones left out."""
    return ' '.join(p for p in parts if p)


def by_condition(d, cond, col):
    alias = {i: c.node_id for c in cond.itertuples() for i in c.ot_ids.split()}
    d = d.assign(node_id=d.node_id.map(lambda i: alias.get(i, i)), v=num(d[col]))
    return d.groupby(['node_id', 'disease_gene']).agg(n=('study_id', 'nunique'), best=('v', 'max'))


def category_order(P, cond):
    """The categories of the related conditions in the order shown: the profile's, then any other the table holds."""
    order = list(P.category_order)
    return order + sorted(set(cond.category) - set(order))


def caveats(P, rg, single, n_readers, wrong_links, automatic):
    """The caveats of the snapshot. The first block concerns the relations and the GWAS regions of the disease: each
    applies only when the snapshot has what it is about, and a disease's text/knowledge.yaml (caveats: {key: text})
    replaces it by its own version; the profile's caveats and the other keys of that file follow. The second block
    is generic."""
    s, own = P.short, P.text('knowledge', 'caveats', {}) or {}
    pub = rg[rg.source.isin(['', 'published'])]
    lead_only = rg[rg.variant_id.ne('') & rg.disease_study_locus_id.eq('')]
    first = {'relations': f'Relations to {s} are evidence of association (clinical studies, ontology, shared genetic '
                          'signal), not of causation.'}
    if len(pub):
        first['regions'] = (f'The {s} GWAS regions are those of the publication, named as there. Colocalization is '
                            f'possible only at the regions where Open Targets holds {a_(s)} credible set.')
    elif len(rg):
        first['regions'] = (f'The {s} GWAS regions are the genome-wide significant credible sets of the {s} GWAS in '
                            'Open Targets, each named by its top locus-to-gene gene, or by its position when that '
                            'gene\'s score is below the threshold.')
    if rg.code.eq('G*').any():
        first['secondary'] = ('G* regions are genome-wide significant only in a secondary analysis of the published '
                              'GWAS, so they are marked apart from G.')
    if len(rg):
        first['colocalization'] = sentence(
            'Colocalization is eCAVIAR on fine-mapped credible sets (CLPP) or coloc (H4).',
            f'Lead-variant overlap is weaker: it takes the {s} lead variant as a stand-in for an unfine-mapped {s} '
            'signal, and its PIP in the other trait\'s credible set bounds the CLPP only if that variant is causal.'
            if len(lead_only) else '')
    if len(single):
        top = single.shared_regions.value_counts()
        top = sorted(top.items(), key=lambda x: (-x[1], x[0]))[0]
        first['single_region'] = (f'A trait sharing a single {s} region may reflect the pleiotropy of that region '
                                  f'(most often {top[0]}, {plural(top[1], "trait")}), so single-region relations are '
                                  'weaker than multi-region ones.')
    out = [own.get(k) or t for k, t in first.items()]
    out += list(P.caveats) + [t for k, t in own.items() if k not in first and k not in
                              ('relations', 'regions', 'secondary', 'colocalization', 'single_region') and t]
    reader = ('two independent language-model readers' if n_readers >= 2 else 'a language model')
    out += [
        'Locus-to-gene (L2G) scores of the related conditions are predictions of the causal gene, not proof.',
        'Curated gene evidence comes from curated resources (ClinGen, PanelApp, Orphanet, ClinVar, UniProt, '
        'Gene2Phenotype, gene burden) for the condition or its descendants in the ontology; it is mostly rare-variant, '
        'Mendelian evidence.',
        f'The literature studies were screened and extracted by {reader} from PubMed abstracts and open-access Methods '
        'sections; ' + ('an item enters when both readers report it or an adjudicating call confirms it (section '
                        '10); ' if n_readers >= 2 else '') +
        'every gene, variant, number and quote was checked against the source text, and items that failed were '
        'dropped. Studies behind paywalls contribute their abstract only. Retracted articles and conference material '
        'are left out (section 10).',
        'The place of a study is the country of its first author\'s institution (PubMed affiliation); a multinational '
        'study takes the place of its lead institution, and the population the abstract states is kept beside it.',
        'The evidence ladder counts studies, not independent samples: some studies share cases. ' +
        ('A gene with several Reactome top-level pathways sits in the one that holds the most of the disease\'s genes'
         + ('; a gene the curated themes list sits in the first that lists it.' if P.themes else '.')
         if automatic else 'A gene with several functions sits in the first gene theme that lists it.'),
        'Generated summaries are syntheses of the listed source records; they cite only those records, but they are '
        'not a systematic review.',
        'Text-mined target links can map an abbreviation to the wrong gene. ' +
        ('The links of this snapshot were reviewed and the wrong mappings found were removed (section 10); others may '
         'remain, so the matched words are listed with every link.' if wrong_links else
         'The matched words are listed with every link.')]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--profile', required=True, help='the disease profile (profile.py)')
    ap.add_argument('--tables', required=True)
    ap.add_argument('--curated', required=True, help='curate.py outputs')
    ap.add_argument('--graph-schema', required=True, help='build_graph.py schema.json')
    ap.add_argument('--summaries', required=True)
    ap.add_argument('--changes', required=True)
    ap.add_argument('--changes-json', required=True)
    ap.add_argument('--seed-pubmed', required=True, help='the PubMed records of the seed (citations of its facts)')
    ap.add_argument('--run-params', required=True, help='the parameters of this run (JSON)')
    ap.add_argument('--l2g-min', type=float, default=0.5)
    ap.add_argument('--top-genes', type=int, default=8, help='GWAS genes listed per related condition')
    ap.add_argument('--max-changes', type=int, default=40)
    ap.add_argument('--out-md', required=True)
    ap.add_argument('--out-json', required=True)
    a = ap.parse_args()
    P = load_profile(a.profile)
    s, core = P.short, P.id
    T = Path(a.tables)
    st = json.loads((T / 'tidy.json').read_text())
    # the columns each table must have, so that a table without rows (a disease without them) reads as empty
    cond = read_tsv(T / 'conditions.tsv', ['node_id', 'name', 'category', 'ot_ids', 'genetics', 'sources',
                                           'shared_regions', 'n_shared_regions', 'n_gwas_studies', 'n_eas_studies', 'genetics_note',
                                           'n_credible_sets'])
    rel = read_tsv(T / 'relations.tsv', ['node_id', 'evidence_type', 'relation', 'detail', 'east_asian', 'pmid'])
    co = read_tsv(T / 'colocalisation.tsv', ['node_id', 'disease_gene', 'study_id', 'clpp'])
    ov = read_tsv(T / 'overlap.tsv', ['node_id', 'disease_gene', 'study_id', 'pip'])
    lad = read_tsv(T / 'disease_genes.tsv', ['symbol', 'level', 'code', 'theme', 'curated_sources'])
    rg = read_tsv(T / 'disease_regions.tsv', ['gene', 'rsid', 'chrom', 'pos_grch37', 'pos_grch38', 'tier', 'code',
                                              'source', 'variant_id', 'disease_study_locus_id'])
    studies = read_tsv(T / 'gwas_studies.tsv', ['node_id', 'condition', 'study_id', 'project', 'first_author', 'year',
                                                'ancestry', 'eas', 'n_samples', 'summary_statistics'])
    loci = read_tsv(T / 'gwas_loci.tsv', ['node_id', 'study_id', 'rsid', 'variant_id', 'finemapping'] +
                    [f'l2g{k}_{c}' for k in (1, 2, 3) for c in ('symbol', 'score')])
    men = read_tsv(T / 'mendelian.tsv', ['node_id', 'symbol', 'source'])
    gc = read_tsv(T / 'gene_condition.tsv', ['symbol', 'node_id', 'evidence'])
    cst = read_tsv(T / 'disease_studies.tsv', ['pmid', 'year', 'first_author', 'design', 'country', 'ancestry',
                                               'n_cases', 'n_controls', 'gene_symbols', 'genes', 'methods_source',
                                               'place', 'japanese', 'east_asian', 'design_class_label'])
    cfi = read_tsv(T / 'disease_findings.tsv', ['pmid', 'gene', 'symbol', 'result', 'association_with',
                                                'association_note', 'variant_name'])
    cmf = read_tsv(T / 'disease_methods.tsv', ['pmid'])
    LADDER, codes, evidence_codes = load_ladder(T)
    CODES = list(codes)
    Q = Path(a.curated)
    pubs, clog = read_tsv(Q / 'publications.tsv'), read_tsv(Q / 'curation_log.tsv')
    mr = read_tsv(Q / 'model_review.tsv') if (Q / 'model_review.tsv').exists() else pd.DataFrame()
    cur = json.loads((Q / 'curation.json').read_text())
    graph = json.loads(Path(a.graph_schema).read_text())['graph']
    summ = read_tsv(a.summaries, ['kind', 'title', 'genes', 'pmids', 'n_targets', 'summary', 'grounded', 'problems',
                                  'check'])
    changes, chj = read_tsv(a.changes), json.loads(Path(a.changes_json).read_text())
    themes = load_themes(T)
    automatic = any(t.get('automatic') for t in themes)     # Reactome themes, alone or for the genes not curated
    choices, _refs = load_choices(P, json.loads(Path(a.run_params).read_text()), T)
    choice = {it['key']: it for g in choices for it in g['items']}
    basis = lambda k: (choice.get(k) or {}).get('basis', '')  # noqa: E731
    value = lambda k: (choice.get(k) or {}).get('value_text', '')  # noqa: E731
    pub = {r['pmid']: r for r in load_json(a.seed_pubmed)['records']}
    who = lambda p: f'{pub[p]["first_author"].split()[0]} {pub[p]["year"]}' if p in pub else ''  # noqa: E731
    cite = lambda p: f'{who(p)}, PMID {p}' if p in pub else f'PMID {p}'  # noqa: E731
    names = dict(zip(cond.node_id, cond.name))
    qs = [P.fill(q) for q in QUESTIONS]
    cats = category_order(P, cond)
    lad['level'] = lad.level.astype(int)
    cond['nreg'] = num(cond.n_shared_regions).fillna(0).astype(int)
    cst['population'] = [study_population(r.country, r.ancestry) for r in cst.itertuples()]
    cst['place'] = cst.place.replace('', 'not known')
    cst['jp'] = [is_japanese(r) for r in cst.itertuples()]
    cst['eas'] = [is_east_asian(r) for r in cst.itertuples()]
    short = {t['name']: t['short'] for t in themes}
    short[OTHER_GENES] = OTHER_GENES
    theme_cols = [t['name'] for t in themes] + [OTHER_GENES]
    pubr, otr = rg[rg.source.isin(['', 'published'])], rg[rg.source.eq('Open Targets')]
    src = ((P.gwas or {}).get('source') or {})
    src_pmid = str(src.get('pmid', '') or '')
    readers = (st.get('extract') or {}).get('readers') or []
    K = {'release': st['release'], 'fetch_date': st['fetch_date'], 'questions': qs,
         'core': {'id': core, 'name': P.name, 'short': s, 'key': P.key},
         'counts': {k: v for k, v in st.items() if isinstance(v, (int, dict))}}
    # the full name beside the short name, which can name other things too (KD is also Kennedy disease)
    title = s if P.name.lower() == s.lower() else f'{P.Name} ({s})'
    md = [f'# {title} knowledge graph (Open Targets {st["release"]}, PubMed; fetched {st["fetch_date"]})', '',
          'Generated by DiseaseLens; do not edit by hand. Sections are built from the '
          f'snapshot tables, in the order of the workbook: {s}\'s genes, its related conditions, what they share, '
          'the context, then the basis of every choice. Text under a "Generated summary" heading was written by a '
          'language model from the records listed with it and checked against those records. Figures: '
          '`../09.figures/` (catalog in its README); the graph for programs: `../05.graph/` (schema in its '
          f'README); the viewer: `../{P.stem}.html`.', '']

    # 1 scope ───────────────────────────────────────────────────────────────
    bs = st.get('conditions_by_source') or {}
    bylv = lad.level.value_counts()
    if len(rg):
        parts = []
        if len(pubr):
            parts.append(f'the {len(pubr)} genome-wide significant regions of the publication '
                         f'({who(src_pmid) or "PMID " + src_pmid}); Open Targets holds '
                         f'{int(pubr.variant_id.ne("").sum())} lead variants and {a_(s)} credible set at '
                         f'{int(pubr.disease_study_locus_id.ne("").sum())}')
        if len(otr):
            parts.append(f'{len(otr)} from the genome-wide significant credible sets of the {s} GWAS in Open Targets')
        regions_line = f'- {s} GWAS regions: ' + '; '.join(parts) + '.'
    else:
        regions_line = f'- {s} GWAS regions: none in this snapshot (section 2.1).'
    n_nr = len(not_regions_of(T))
    if n_nr:
        regions_line = regions_line.replace(' (section 2.1).', '.') + (
            f' {n_nr} genome-wide significant credible set{"s" if n_nr != 1 else ""} of its GWAS in Open Targets '
            f'{"are" if n_nr != 1 else "is"} no region, the lead variant not credible (section 2.1).')
    n_read = len(readers)
    md += ['## 1. Scope', '']
    if st.get('evidence_also'):
        md += [f'- Open Targets entries: {P.id} ({P.name}), with the evidence of '
               + '; '.join(f'{x["id"]} ({x["name"]})' for x in st['evidence_also'])
               + ' (the profile\'s evidence_also: targets, text-mined literature, GWAS and curated genes merged).']
    md += [
           f'- {s} genes: {len(lad)} genes studied genetically in {s} ({int((lad.level <= 5).sum())} with a '
           'positive result), placed on an evidence ladder (levels 1–7: ' +
           ', '.join(f'{int(bylv.get(k, 0))}' for k in sorted(LADDER)) + ') and in a function (gene theme).',
           regions_line,
           f'- Related conditions: {st["conditions"]} ({bs.get("clinical", 0)} by clinical evidence, PubMed-verified; '
           f'{bs.get("ontology", 0)} by the disease ontology; {bs.get("colocalization", 0)} by colocalization and '
           f'{bs.get("lead-variant overlap", 0)} by lead-variant overlap with {a_(s)} GWAS region; a condition can have '
           'several).',
           f'- Genetics of the related conditions: {st["gwas_studies"]} GWAS studies ({st["eas_studies"]} East Asian) '
           f'with {st["credible_sets"]:,} credible sets; {st["mendelian_links"]} curated condition-gene links.',
           f'- {s} literature: {st["disease_genetic_records"]} PubMed records screened, '
           f'{st["disease_genetic_studies"]} genetic studies of {s} ({st["disease_findings"]} findings, '
           f'{st["methods_fields"]} methods fields), after the curation left out {cur["studies_excluded"]}'
           + (f'; read by {reader_count(n_read)}' if n_read else '')
           + (f'; {(st.get("extract") or {}).get("not_read")} records not read in this build (no language model was '
              'called, and no answer was cached)' if (st.get("extract") or {}).get("not_read") else '')
           + ' (section 10).',
           f'- Curation: {cur["publications"]:,} publications checked in PubMed, {cur["excluded"]} left out '
           f'(retracted or conference material); {cur["textmining_dropped"]} text-mined links dropped; every '
           'finding classified by what it is about (section 10).',
           f'- Choices: {sum(len(g["items"]) for g in choices)}, each with its value in this run and its basis '
           '(section 11).',
           f'- Graph: {graph["n_nodes"]:,} nodes, {graph["n_edges"]:,} edges of '
           f'{len(graph["counts"]["nodes"])} node and {len(graph["counts"]["edges"])} edge types (`05.graph`).', '',
           'Sections 2 to 11 answer these questions, in four parts and in the order of the workbook:', '']
    for name, secs in PARTS:
        md += [f'- {P.fill(name)} ({secs})']
    md += [''] + [f'{k}. {q}' for k, q in enumerate(qs, 1)] + ['']
    K['graph'] = graph

    # 2 the disease's genes ─────────────────────────────────────────────────
    md += [f'## 2. {s} genes: strength and function of the evidence', ''] + question(qs, 0)
    md += [sentence('Each gene sits at its strongest level of evidence (rows) and in the column of its function (gene '
                    'theme).', basis('ladder'), basis('about'), basis('themes')), '']
    if len(lad):
        md += figure('evidence_grid', f'Evidence level x function of the {s} genes')
    rows = []
    for lv in sorted(LADDER):
        sub = lad[lad.level == lv]
        parts = [f'{short.get(t, t)}: ' + ', '.join(sub[sub.theme == t].symbol)
                 for t in theme_cols + sorted(set(sub.theme) - set(theme_cols)) if (sub.theme == t).any()]
        code, text = LADDER[lv]
        rows.append([lv, code or MISSING, text, len(sub), '; '.join(parts) or MISSING])
    md += md_table(['Level', 'Code', 'Evidence', 'Genes', 'Genes by function'], rows) + ['']
    # the curated gene-disease evidence of the disease itself (code C; a gene with G as well keeps G)
    curg = lad[lad.curated_sources.ne('')]
    cv = value('curated_core')
    head = f'Curated gene–disease evidence for {s} or its descendants in the ontology (code C' + \
        (f'; {cv}' if cv else '') + ')'
    if len(curg):
        md += [f'{head}: {plural(len(curg), "gene")}, ' + '; '.join(
            f'{r.symbol} ({r.curated_sources}{", level 1 by G as well" if r.code == "G" else ""})'
            for r in curg.itertuples()) + '.']
    else:
        md += [f'{head}: no gene in Open Targets in this snapshot.']
    if cur.get('core_curated_dropped'):
        md[-1] += (' Left out: ' + ', '.join(cur['core_curated_dropped']) + ' (the gene symbol is the abbreviation '
                   'of the disease or of an ontology neighbor; `03.curated/curation_log.tsv`).')
    md += ['']
    K['ladder'] = lad.to_dict('records')

    # 2.1 the GWAS regions of the disease
    head = ('The published GWAS regions' if len(rg) and len(pubr) == len(rg) else
            'The GWAS regions in Open Targets' if len(rg) and len(otr) == len(rg) else 'The GWAS regions')
    md += [f'### 2.1 {head}', '']
    if len(rg):
        intro = []
        if len(pubr):
            table = str(src.get('table', '') or '').split('.')[0]
            intro.append(P.text('knowledge', 'regions', '') or
                         f'The genome-wide significant regions of the publication ({cite(src_pmid)}'
                         + (f', {table}' if table else '') + '), named as there'
                         + ('; the regions significant only in a secondary analysis are marked G*'
                            if pubr.code.eq('G*').any() else '') + '.')
        if len(otr):
            intro.append(f'{plural(len(otr), "region")} from the genome-wide significant credible sets '
                         f'({value("genome_wide") or "lead P < 5 × 10⁻⁸"}) of the {s} GWAS in Open Targets'
                         + (f', merged when their {value("region_window")}' if value('region_window') else '') + '.')
        md += [sentence(*intro, basis('region_names')), '']
        l2g = {}
        for r in loci[loci.node_id == core].itertuples():
            txt = ', '.join(f'{getattr(r, f"l2g{k}_symbol")} {float(getattr(r, f"l2g{k}_score")):.2f}'
                            for k in (1, 2, 3) if getattr(r, f'l2g{k}_symbol'))
            for key in (r.rsid, r.variant_id):
                if key:
                    l2g.setdefault(key, txt)
        labels = ((P.gwas or {}).get('columns') or {})
        stats = [(k, (labels.get(k) or d)[:1].upper() + (labels.get(k) or d)[1:]) for k, d in REGION_STATS
                 if k in rg and rg[k].ne('').any()]
        tier = 'tier' in rg and rg.tier.ne('').any()
        b37, b38 = rg.pos_grch37.ne('').all(), rg.pos_grch38.ne('').all()
        pos_head = 'Position (GRCh37)' if b37 else 'Position (GRCh38)' if b38 else 'Position'
        mixed = rg.source.nunique() > 1

        def pos(r):
            if b37 or (not b38 and r.pos_grch37):
                return f'chr{r.chrom}:{r.pos_grch37}' + ('' if b37 else ' (GRCh37)')
            return f'chr{r.chrom}:{r.pos_grch38}' + ('' if b38 else ' (GRCh38)')
        rows = []
        for _i, r in rg.iterrows():
            vals = [p_text(r[k]) if k in ('p', 'p_pe', 'v') else (r[k] or MISSING) for k, _h in stats]
            rows.append([r.gene, r.rsid or MISSING, pos(r)] + ([r.tier or MISSING] if tier else []) + [r.code] + vals
                        + ['credible set' if r.disease_study_locus_id else 'lead variant only' if r.variant_id
                           else 'not held', l2g.get(r.rsid) or l2g.get(r.variant_id) or MISSING]
                        + ([r.source or 'published'] if mixed else []))
        md += md_table(['Region', 'Lead variant', pos_head] + (['Tier'] if tier else []) + ['Code']
                       + [h for _k, h in stats] + ['Open Targets', 'Open Targets locus-to-gene']
                       + (['Source'] if mixed else []), rows) + ['']
    else:
        md += [f'No GWAS region of {s} in this snapshot' +
               ('' if P.gwas else ': the disease profile gives no published table of regions, and Open Targets holds '
                f'no genome-wide significant credible set of a GWAS of {s}'
                + (' with a credible lead variant' if len(not_regions_of(T)) else '')) + '. The evidence ladder '
               'therefore starts from the curated gene evidence and the genetic studies in the literature.', '']
    nr = not_regions_of(T)
    if len(nr):
        md += [sentence(not_regions_sentence(nr, s), basis('region_mac')), '']
        md += md_table(['Lead variant', 'Study', 'P', 'MAF', 'Expected copies', 'Why not a region'],
                       [[r.rsid or r.variant_id, r.study_id, p_unicode(r.p_mantissa, r.p_exponent), r.maf,
                         r.expected_mac, r.reason] for r in nr.itertuples()]) + ['']
    K['regions'] = rg.to_dict('records')
    K['signals_not_regions'] = nr.to_dict('records')

    fx = assign_themes(cfi, themes).merge(cst[['pmid', 'year', 'first_author', 'eas']], on='pmid', how='left')
    gt = gene_tallies(fx, P.about, P.within)
    rep = gt[(gt.associated >= 2) & (gt.theme != NO_GENE)].sort_values('associated', ascending=False)
    md += [f'### 2.2 What the genetic studies of {s} report, by gene theme', '',
           f'{len(cst)} studies, {len(cfi)} findings on {gt.label.nunique()} genes or gene sets. Associated in two or '
           'more studies: ' + (', '.join(f'{r.label} ({r.associated})' for r in rep.itertuples()) or 'none') + '.',
           '']
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
        md += [f'#### {th}', '', (desc.get(th, '') + ' ').lstrip() + f'{plural(len(g), "gene")}, '
               f'{plural(f.pmid.nunique(), "study", "studies")}, {plural(len(f), "finding")}.', '']
        rows = []
        for r in g.itertuples():
            ff = f[f.label == r.label].assign(_y=lambda d: num(d.year)).sort_values('_y')
            refs = '; '.join(dict.fromkeys(
                f'{x.first_author} {x.year} ({x.variant_name + ", " if x.variant_name else ""}{x.result}'
                f'{"" if x.association_with == P.about else " with " + x.association_with}; PMID {x.pmid})'
                for x in ff.itertuples()))
            rows.append([r.label, level.get(r.label, MISSING), tally_text(r, s), refs])
        md += md_table(['Gene', 'Level', 'Studies', 'Findings (variant, result)'], rows) + ['']
        K['themes'].append({'theme': th, 'genes': g.to_dict('records')})
    if P.facts:
        md += ['### 2.3 Curated facts', ''] + [f'- {f["topic"] + ": " if f.get("topic") else ""}{f["statement"]} '
                                               f'({cite(str(f["pmid"]))})' for f in P.facts] + ['']
    K['disease_genetics'] = {'regions': len(rg), 'literature_studies': len(cst), 'findings': len(cfi),
                             'facts': P.facts}

    # 3 genetic studies of the disease ────────────────────────────────────
    sm = summ[(summ.kind == 'methods') & (summ.grounded == 'True')]
    md += [f'## 3. Genetic studies of {s}: designs, places and methods', ''] + question(qs, 1)
    md += [sentence(basis('design_class'), basis('place')), '']
    if len(cst):
        if num(cst.year).notna().any():         # the figure places each study by its year
            md += figure('genetic_studies', f'The genetic studies of {s} by year, design class and place')
        nc = num(cst.n_cases)
        md += ['- Design classes (one per study): ' + counts(cst.design_class_label) + '.',
               '- Designs (a study can have several): ' + counts(cst.design.str.split('; ').explode()) + '.',
               '- Places (country of the first author\'s institution): ' + counts(cst.place) + '.',
               '- Populations as the abstracts state them: ' + counts(cst.population) + '.',
               (f'- Cases per study: median {nc.median():.0f}, range {nc.min():.0f}–{nc.max():.0f} (where stated).'
                if nc.notna().any() else '- Cases per study: stated by no study.'),
               f'- Methods fields ({len(cmf)}: platform, sequencing, quality control, statistics, covariates, '
               'multiple testing): from the open-access Methods section for '
               f'{cst.methods_source.eq("full text").sum()} studies, from the abstract for '
               f'{cst.methods_source.eq("abstract").sum()} (no full text).', '']
        rows = [[r.year, r.first_author, r.design_class_label, r.place, r.population, r.n_cases or MISSING,
                 r.n_controls or MISSING, r.methods_source or MISSING, r.pmid]
                for r in cst.assign(_y=num(cst.year)).sort_values('_y').itertuples()]
        md += md_table(['Year', 'First author', 'Design class', 'Place', 'Population', 'Cases', 'Controls',
                        'Methods from', 'PMID'], rows) + ['']
    else:
        md += [f'No genetic study of {s} was kept in this snapshot.', '']
    if len(sm):
        md += ['#### Generated summary of the methods', '', plain(sm.iloc[0].summary), '']

    # 4 related conditions ──────────────────────────────────────────────────
    md += [f'## 4. Conditions related to {s}, and why', ''] + question(qs, 2)
    md += [sentence(basis('categories'), basis('clinical_rule')), '']
    md += figure('condition_map', f'{s} and its related conditions')
    K['conditions'] = []
    for cat in cats:
        sub = cond[cond.category == cat]
        if not len(sub):
            continue
        md += [f'### {cat}', '']
        if cat == SHARED_SIGNAL:
            multi, single = sub[sub.nreg >= 2], sub[sub.nreg == 1]
            md += [sentence(f'Traits sharing a GWAS signal with {a_(s)} region and with no clinical or ontology '
                            f'relation: {len(sub)}.', basis('loci_grading')), '',
                   f'Sharing two or more {s} regions ({len(multi)}):', '']
            if len(multi):
                md += md_table(['Trait', f'{s} regions shared', 'Regions'],
                               [[r.name, ', '.join(r.shared_regions.split()), r.nreg]
                                for r in multi.sort_values(['nreg', 'name'], ascending=[False, True]).itertuples()])
            else:
                md += ['None.']
            md += ['', f'Sharing one {s} region only ({len(single)}), by region:', '']
            for reg, g in single.groupby('shared_regions'):
                md += [f'- {reg} ({len(g)}): ' + ', '.join(sorted(g.name, key=str.lower)) + '.']
            none = sub[sub.nreg < 1]
            if len(none):
                md += [f'- No shared region recorded ({len(none)}): ' + ', '.join(sorted(none.name, key=str.lower))
                       + '.']
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
            rows.append([c.name, c.sources.replace(' ', ', '), ', '.join(c.shared_regions.split()) or MISSING,
                         '<br>'.join(items)])
        md += md_table(['Condition', 'Evidence types', f'{s} regions shared', 'Evidence'], rows) + ['']

    # 5 shared regions ──────────────────────────────────────────────────────
    ca, oa = by_condition(co, cond, 'clpp'), by_condition(ov, cond, 'pip')
    md += [f'## 5. The {s} GWAS regions the related conditions share', ''] + question(qs, 3)
    if len(rg):
        lead_only = rg[rg.variant_id.ne('') & rg.disease_study_locus_id.eq('')]
        md += [sentence(f'Colocalization: {value("clpp")} or {value("h4")}, possible only at a region with {a_(s)} '
                        'credible set in Open Targets.', basis('clpp')), '']
        if len(lead_only):
            md += [sentence(f'Lead-variant overlap (weaker): {value("overlap")}, used only where Open Targets holds '
                            f'the {s} lead variant without a credible set.', basis('overlap')), '']
        md += figure('shared_regions', f'Related conditions x {s} GWAS regions')
        rows = []
        for r in rg.itertuples():
            c_n = sorted({n for n, g in ca.index if g == r.gene})
            o_n = sorted({n for n, g in oa.index if g == r.gene} - set(c_n))
            how = ('colocalization' if r.disease_study_locus_id else 'lead-variant overlap' if r.variant_id
                   else 'none (lead variant not in Open Targets)')
            n_sig = ca[ca.index.get_level_values(1) == r.gene].n.sum() + \
                oa[oa.index.get_level_values(1) == r.gene].n.sum()
            rows.append([r.gene, r.code, how, len(c_n) + len(o_n), int(n_sig)])
        md += md_table(['Region', 'Code', 'How a condition can share it', 'Conditions sharing it', 'GWAS signals'],
                       rows)
        md += ['', 'Conditions sharing two or more regions, or clinically related and sharing one (a ° marks a region '
               'shared by lead-variant overlap):', '']
        sh = cond[(cond.nreg >= 2) | ((cond.category != SHARED_SIGNAL) & (cond.nreg >= 1))].copy()
        sh['_o'] = sh.category.map(cats.index)
        rows = []
        for c in sh.sort_values(['_o', 'nreg', 'name'], ascending=[True, False, True]).itertuples():
            regs = [g + ('' if (c.node_id, g) in ca.index else '°') for g in rg.gene
                    if (c.node_id, g) in ca.index or (c.node_id, g) in oa.index]
            rows.append([c.name, c.category, ', '.join(regs), c.nreg])
        md += (md_table(['Condition', 'Category', f'{s} regions', 'Regions'], rows) if rows else ['None.']) + ['']
    else:
        rows = []
        md += [f'{s} has no GWAS region in this snapshot (section 2.1), so no related condition can share one.', '']
    K['shared_regions'] = rows

    # 6 shared genes ───────────────────────────────────────────────────────
    gc['code'] = gc.evidence.map(evidence_codes).fillna('')
    gc = gc[gc.code.ne('')]
    rank = lambda c: CODES.index(c) if c in CODES else len(CODES)  # noqa: E731
    cell = gc.groupby(['symbol', 'node_id']).code.apply(lambda x: ''.join(sorted(set(x), key=rank)))
    core_code = cell.xs(core, level='node_id') if core in cell.index.get_level_values(1) else pd.Series(dtype=str)
    relc = cell[(cell.index.get_level_values(1) != core) & cell.index.get_level_values(0).isin(set(lad.symbol))]
    nrel = relc.groupby(level=0).size()
    md += [f'## 6. Genes shared by {s} and the related conditions', ''] + question(qs, 4)
    if len(pubr):
        g_core = 'the gene the publication names for a genome-wide significant region'
    elif len(otr):
        g_core = (f'the top locus-to-gene gene, score ≥ {a.l2g_min}, of a genome-wide significant region of its GWAS '
                  'in Open Targets')
    else:
        g_core = 'none in this snapshot, which has no GWAS region'
    gstar = (f'; G*, {s} only, ' + (P.text('knowledge', 'code_g_star', '') or
                                    'a region significant only in a secondary analysis of the published GWAS')
             if rg.code.eq('G*').any() else '')
    md += [sentence(f'Codes: G, GWAS (for {s}, {g_core}; for a related condition, the top locus-to-gene gene, score ≥ '
                    f'{a.l2g_min}, of one of its GWAS credible sets){gstar}; C, curated gene–disease evidence; L, {s} '
                    f'only, a genetic study in the literature; –, {s} only, studied without a positive result (levels '
                    '6–7).', basis('g_codes')), '']
    if len(lad):
        md += figure('gene_map', f'Genes studied in {s} x the related conditions')
    only = lad[~lad.symbol.isin(nrel.index)]
    md += [f'{len(lad)} genes were studied genetically in {s} ({int((lad.level <= 5).sum())} with a positive result); '
           f'{len(nrel)} of them also have GWAS or curated evidence in at least one related condition. '
           + not_in_gene_map(cond, set(relc.index.get_level_values(1)) if len(relc) else set(),
                             'section 5' if len(rg) else '')[1], '']
    rows = []
    for x in lad[lad.symbol.isin(nrel.index)].assign(n=lambda d: d.symbol.map(nrel)).sort_values(
            ['level', 'n', 'symbol'], ascending=[True, False, True]).itertuples():
        cc = relc.xs(x.symbol, level=0)
        rows.append([x.symbol, x.level, core_code.get(x.symbol, '–'), int(x.n),
                     '; '.join(f'{names.get(n, n)} {k}' for n, k in sorted(cc.items(), key=lambda y: names.get(
                         y[0], y[0]).lower()))])
    if rows:
        md += md_table(['Gene', 'Level', s, 'Related conditions', 'Which (code)'], rows) + ['']
    if len(only):
        md += [f'Specific to {s} in this map ({len(only)}: no GWAS or curated evidence in any related condition), by '
               'level: ' + '; '.join(f'{lv} ' + ', '.join(only[only.level == lv].symbol)
                                     for lv in sorted(only.level.unique())) + '.', '']
    gcs = gc.groupby('symbol').node_id.apply(lambda x: len(set(x) - {core}))
    other = gcs[~gcs.index.isin(set(lad.symbol)) & (gcs >= 3)].sort_values(ascending=False)
    md += [f'Not studied in {s}, {len(other)} genes have genetic evidence in three or more related conditions' +
           ('; the most widely shared: ' + ', '.join(f'{g} ({n})' for g, n in other.head(15).items()) if len(other)
            else '') + '.', '']
    K['shared_genes'] = rows

    # 7 genetics of related conditions ──────────────────────────────────────
    md += ['## 7. Genetics of the related conditions', ''] + question(qs, 5)
    md += [sentence(basis('genetics_mode'), basis('curated')), '']
    rows = []
    loci['l2g1_score_n'] = num(loci.l2g1_score)
    for c in cond[cond.genetics == 'full'].sort_values(['category', 'name']).itertuples():
        L = loci[(loci.node_id == c.node_id) & (loci.l2g1_score_n >= a.l2g_min)]
        top = L.groupby('l2g1_symbol').study_id.nunique().sort_values(ascending=False).head(a.top_genes)
        M = men[men.node_id == c.node_id] if len(men) else men
        mg = M.groupby('symbol').source.apply(lambda x: '/'.join(sorted(set(x)))) if len(M) else pd.Series(dtype=str)
        rows.append([c.name, c.category, f'{c.n_gwas_studies} ({c.n_eas_studies} EAS)', c.n_credible_sets,
                     ', '.join(f'{g} ({n})' for g, n in top.items()) or MISSING,
                     ', '.join(f'{g} [{s_}]' for g, s_ in list(mg.items())[:12]) or MISSING])
    md += md_table(['Condition', 'Category', 'GWAS studies', 'Credible sets', 'Top GWAS genes (studies)',
                    'Curated genes [source]'], rows) + ['']

    gw = studies.copy()
    gw['n_samples'] = num(gw.n_samples)
    md += [f'### How the GWAS were done (Open Targets: {s} and the related conditions)', '']
    if len(gw):
        md += [f'- {len(gw)} studies; sample size median {gw.n_samples.median():,.0f} (range '
               f'{gw.n_samples.min():,.0f}–{gw.n_samples.max():,.0f}); with summary statistics '
               f'{gw.summary_statistics.eq("True").sum()}.',
               '- Fine-mapping of the credible sets: ' + (counts(loci.finemapping) or MISSING) + '.',
               '- Largest ancestry group of each study: ' +
               counts(gw.ancestry.str.split('; ').str[0].str.split().str[0]) + '.', '']
    else:
        md += ['No GWAS study in Open Targets for these conditions.', '']

    # 8 mechanisms ──────────────────────────────────────────────────────────
    md += [f'## 8. Mechanisms (text-mined {s} literature, by Reactome theme)', ''] + question(qs, 6)
    md += [basis('reactome'), '']
    mech = summ[summ.kind == 'mechanism'].copy()
    mech['n'] = num(mech.n_targets)
    for r in mech.sort_values('n', ascending=False).itertuples():
        md += [f'### {r.title}' + (f' ({int(r.n)} targets)' if r.n == r.n else ''), '']
        if r.grounded == 'True':
            md += ['Generated summary: ' + plain(r.summary), '']
        elif r.problems.startswith('not generated'):
            md += ['(No generated summary in this build: none of these source sentences is cached, and this '
                   'build calls no language model.)', '']
        else:
            md += [f'(Summary withheld: {r.problems}.)', '']
        md += [f'Source targets: {r.genes}. Source PMIDs: {r.pmids}.', '']
    if not len(mech):
        md += ['No Reactome theme holds enough text-mined targets for a summary in this snapshot.', '']
    K['mechanisms'] = mech.drop(columns='n').to_dict('records')

    # 9 East Asian evidence ─────────────────────────────────────────────────
    eas_all = studies[studies.eas == 'True']
    eas_own, eas = eas_all[eas_all.node_id == P.id], eas_all[eas_all.node_id != P.id]
    # the relations whose own evidence is East Asian: clinical evidence that relates the condition, or a shared signal
    # whose two GWAS have mostly East Asian samples
    ea_rel = rel[(rel.east_asian == 'True') & (rel.relation.str.lower() != NO_ASSOCIATION)]
    ea_lit = cst[cst.eas.astype(bool)]
    md += ['## 9. East Asian evidence', ''] + question(qs, 7)
    md += [f'- Relations whose own evidence is East Asian ({len(ea_rel)}: clinical evidence from East Asia, or a '
           f'shared signal whose two GWAS have mostly East Asian samples)' + (': ' + '; '.join(
               f'{names.get(r.node_id, r.node_id)}: {r.detail} '
               + (f'({cite(r.pmid)})' if r.pmid else '(shared GWAS signal)') for r in ea_rel.itertuples())
               if len(ea_rel) else '') + '.',
           f'- Genetic studies of {s} done in East Asia ({len(ea_lit)}; in Japan {int(ea_lit.jp.sum())})' + (': ' +
               '; '.join(f'{r.first_author} {r.year} ({r.place}; {r.design_class_label}; '
                         f'{r.gene_symbols or r.genes or "no gene"}; PMID {r.pmid})'
                         for r in ea_lit.sort_values('year').itertuples()) if len(ea_lit) else '') + '.',
           *([f'- GWAS of {s} with a mostly East Asian sample ({eas_own.study_id.nunique()}): ' + '; '.join(
               f'{r.study_id} ({r.first_author or r.project} {r.year}, n = {r.n_samples})'
               for r in eas_own.drop_duplicates('study_id').itertuples()) + '.'] if len(eas_own) else []),
           f'- GWAS of the related conditions with a mostly East Asian sample ({eas.study_id.nunique()} studies)'
           + (': ' + '; '.join(f'{r.condition}: {r.study_id} ({r.first_author or r.project} {r.year}, '
                               f'n = {r.n_samples})' for r in eas.itertuples()) if len(eas) else '') + '.', '']

    # 10 curation ───────────────────────────────────────────────────────────
    own_rules = f'diseases/{P.key}/curation.yaml'
    md += ['## 10. Curation: what was left out or corrected, and why', ''] + question(qs, 8)
    md += [f'Every publication the graph uses ({cur["publications"]:,}) was looked up in PubMed on '
           f'{cur["checked_at"][:10]}; the rules are in `config/curation.yaml`'
           + (f' and `{own_rules}`' if own_rules in (P.d.get('files') or {}) else '') + ', and every rule applied is '
           'in `03.curated/curation_log.tsv`. No rule calls a language model.', '', '### Publications left out', '',
           sentence(basis('exclusion'), f'{cur["excluded"]} publications meet the rule; they are part of no '
                    'table, figure or count of this snapshot and are not listed.'), '']
    fl = pubs[pubs.status == 'flagged']
    md += [f'Flagged, kept: {len(fl)}' + (f' ({counts(fl.reason)})' if len(fl) else '') + '.', '',
           f'### Findings that are not about {s} itself', '', basis('about') or
           f'Findings about another trait, or about a phenotype or outcome within {s}, are kept and marked; they do '
           'not raise a gene\'s level.', '']
    fa = cfi[cfi.association_with != P.about].merge(cst[['pmid', 'first_author', 'year']], on='pmid', how='left')
    md += (md_table(['Study', 'Gene', 'Result', 'About', 'Note', 'PMID'],
                    [[f'{x.first_author} {x.year}', x.symbol or x.gene or MISSING, x.result, x.association_with,
                      x.association_note or MISSING, x.pmid] for x in fa.itertuples()])
           if len(fa) else ['None.']) + ['']
    vn = clog[(clog.step == 'findings') & (clog.rule == 'variant name')]
    example = P.text('knowledge', 'variant_examples', '')
    if not example and len(vn):     # the name written in the most ways, with up to three of its spellings
        spell = vn.groupby('after').before.apply(lambda x: list(dict.fromkeys(x)))
        top = sorted(spell.items(), key=lambda x: (-len(x[1]), x[0]))[0]
        quoted = [f'"{b}"' for b in top[1][:3]]
        example = f'e.g. {top[0]} for ' + (', '.join(quoted[:-1]) + ' and ' + quoted[-1] if len(quoted) > 1
                                          else quoted[0])
    md += ['### Variants named', '', sentence(basis('variants'), f'Renamed in this snapshot: {len(vn)} findings'
                                              + (f' ({example})' if example else '') + '.'), '',
           '### Text-mined links mapped to the wrong gene', '', basis('textmining_review'), '']
    tm = clog[(clog.step == 'text mining') & (clog.rule == 'link dropped')]
    wrong = tm[tm.after.str.startswith('mapped to the wrong gene')]
    md += (md_table(['Gene', 'Publications', 'What the word means'],
                    [[x.item, x.before.split()[0], x.after.split(': ', 1)[-1]] for x in wrong.itertuples()])
           if len(wrong) else ['None in this snapshot.']) + ['']
    md += ['### Place and design class of the genetic studies', '',
           sentence(basis('place'), ('Places: ' + counts(cst.place) + '.') if len(cst) else '',
                    basis('design_class')), '']
    # the language-model readers of the literature, their agreement and the items they disputed
    ex = st.get('extract') or {}
    agree, adj = ex.get('agreement') or {}, ex.get('adjudicator')
    desc_reader = lambda x: f'{x.get("backend", "")}' + (  # noqa: E731
        f', whose answers came from {" and ".join(x["models"])}' if x.get('models') else '')
    tab = [(lab, d) for lab, d in (('studies', cst), ('findings', cfi), ('methods fields', cmf))
           if 'agreement' in d and len(d)]
    shown = tab if len(readers) >= 2 else []        # with one reader, every item kept is "one reader"
    sc = summ[(summ.grounded == 'True') & summ.check.ne('')].check.str.split(':').str[0]
    flags = int(cfi.result_check.ne('').sum()) if 'result_check' in cfi else 0
    n_disputed = int(cur.get('model_review', len(mr)))
    md += ['### Language-model readers and their agreement', '',
           sentence(('The genetic studies were read by ' + reader_count(len(readers)) + ' (' + '; '.join(
                        f'reader {x["reader"]}: {desc_reader(x)}' for x in readers) + ').') if readers else
                    'The language-model readers of this snapshot are not recorded.',
                    f'Adjudicator: {desc_reader(adj) if isinstance(adj, dict) else adj}.' if adj else '',
                    ('Items both readers reported or the adjudicator confirmed, by step: ' + '; '.join(
                        f'{k}: {v.get("both", 0)} agreed, {v.get("adjudicated", 0)} adjudicated'
                        for k, v in agree.items()) + '.') if agree else
                    'With one reader, every item rests on that reader alone and is marked "one reader".'
                    if len(readers) == 1 else '',
                    ('Agreement of the items kept: ' + '; '.join(
                        f'{lab}, {counts(d.agreement.replace("", "not recorded"))}' for lab, d in shown) + '.')
                    if shown else '',
                    f'Disputed items: {n_disputed}' + (', each with every reader\'s answer, the final value and the '
                                                       'decision in `03.curated/model_review.tsv`' if n_disputed
                                                       else '') + '.',
                    f'Findings whose result disagrees with the wording of their source sentence: {flags}, kept and '
                    'flagged.',
                    ('Generated summaries shown, by check: ' + counts(sc) + '.') if len(sc) else ''), '']
    K['readers'] = {'readers': readers, 'adjudicator': adj, 'agreement': agree,
                    'kept_by_agreement': {lab: d.agreement.value_counts().to_dict() for lab, d in tab},
                    'disputed_items': n_disputed, 'result_flags': flags,
                    'summary_checks': sc.value_counts().to_dict()}
    K['curation'] = cur

    # 11 choices ────────────────────────────────────────────────────────────
    by = pd.Series([it['basis_type'] for g in choices for it in g['items']]).value_counts()
    md += ['## 11. Choices and their basis', ''] + question(qs, 9)
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
                         '; '.join(r[0] for r in it['refs']) or MISSING])
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

    cav = caveats(P, rg, cond[(cond.category == SHARED_SIGNAL) & (cond.nreg == 1)], len(readers), len(wrong),
                  automatic)
    md += ['## 13. Caveats', ''] + [f'- {c}' for c in cav] + ['']
    K['caveats'] = cav
    md += ['## References', ''] + [f'- {full}' for _s, full in sorted(used.items())] + ['']
    K['references'] = used
    Path(a.out_md).write_text('\n'.join(md))
    Path(a.out_json).write_text(json.dumps(K, indent=1, ensure_ascii=False, default=str))
    log(f'knowledge: {len(md)} lines')


if __name__ == '__main__':
    main()
