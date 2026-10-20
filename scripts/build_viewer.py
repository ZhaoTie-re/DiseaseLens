#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The viewer of one snapshot of one disease: a single HTML file that
#           opens in any browser, offline, and shows the whole knowledge graph.
#           No library, no network; the data are embedded.
#             Overview     counts, the evidence ladder (level x function), the
#                          GWAS regions of the disease x the conditions sharing
#                          them
#             Genes        the genes of the disease by level (or every gene),
#                          filterable
#             Conditions   the related conditions: an evidence matrix (one row
#                          per condition, grouped by category; publications,
#                          best design, direction, East Asian evidence, ontology
#                          relation, each shared GWAS region, its own genetics),
#                          filtered by category or name
#             Regions      the genome-wide significant GWAS regions of the
#                          disease (published, or from Open Targets); a sentence
#                          when there is none
#             Studies      the genetic studies of the disease, by place and
#                          design, with the agreement of the model readers
#             Curation     what was left out or corrected, and why
#             How it is built  the one place that answers how the graph is
#                          built: the general rules, each followed by how it
#                          applies to this disease, at a glance, step by step,
#                          and in sections reachable as #/build/<section>: the
#                          sources, the GWAS regions, the genetic studies, the
#                          evidence ladder and gene functions, the related
#                          conditions (routes, steps and the models' part,
#                          categories, the evidence matrix, what is left out),
#                          the associations (edge types), the language models,
#                          snapshots and files, and every rule with its basis
#           Every node has its own page (#/n/<id>): its facts, a diagram of its
#           neighbors, and every edge with its evidence and source link
#           (PubMed, Open Targets, dbSNP), so that each statement can be traced
#           to where it comes from. Names come from the profile; the texts
#           specific to a disease from its text/viewer.yaml (generic otherwise).
#           Made to be shared as a file, with a print style, so that printing a
#           page (or saving it as PDF) gives a clean document; it works from
#           file://.
# Inputs  : the profile, 05.graph (nodes_*.tsv, edges_*.tsv, schema.json with
#           the ladder), 03.curated, 04.tables/literature.tsv (one sentence per
#           text-mined link), 06.summaries (grounded summaries only),
# Outputs : <key>_knowledge_graph.html
# ---------------------------------------------------------------------------
import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import (EAS_COUNTRIES, NO_LEVEL, OTHER_GENES, RELATION_DESIGNS, SHARED_SIGNAL,  # noqa: E402
                       category_colors, disease_badge, load_choices, load_glosses, load_profile, load_status, log, not_regions_of,
                       not_regions_sentence, read_tsv, relates)

WORDS = ['no', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine', 'ten', 'eleven', 'twelve']


def compact(row, skip=()):
    return {k: v for k, v in row.items() if k not in skip and v not in ('', None)}


def word(n):
    return WORDS[n] if 0 <= n < len(WORDS) else str(n)


def texts(P, regions, ladder, rp, not_regions):
    """The prose of the viewer: the disease's text/viewer.yaml, else a generic text from the profile and the data.
    regions: the region nodes (with `source`); ladder: {level: [code, label]}; rp: the run's parameters."""
    s = P.short
    n_pub = int(regions.source.isin(['', 'published']).sum())
    n_ot = int(regions.source.eq('Open Targets').sum())
    # overview
    name = P.name if P.name.lower() == s.lower() else f'{P.name} ({s})'
    overview = (f'The knowledge graph of {name}: the genes with genetic evidence for {s} and how strong that evidence '
                f'is, the genetic studies it comes from, the conditions related to {s} and the genetics they share '
                'with it. Every statement links to its source.')
    # the note under the ladder: the same seven levels for every disease, each defined in the table above it
    note = (f'Every disease shows the same seven levels; an empty level says why. Levels 1–2 come from GWAS and '
            f'curated gene–disease evidence, levels 3–7 from the genetic studies of {s}; levels 3–6 count only '
            f'findings about {s} itself.')
    # the regions: the published ones (the profile's table), then those from Open Targets
    src = (P.gwas or {}).get('source') or {}
    cite = ', '.join(x for x in (f'PMID {src["pmid"]}' if src.get('pmid') else '',
                                 str(src.get('table') or '').split('.')[0].strip()) if x)
    cite = f' ({cite})' if cite else ''
    lead = []
    if n_pub:
        default = (f'The genome-wide significant region of the published {s} GWAS{cite}, named by the gene the '
                   'publication reports.' if n_pub == 1 else
                   f'The {{n}} genome-wide significant regions of the published {s} GWAS{cite}, named by the gene '
                   'the publication reports.')
        lead.append(P.text('viewer', 'regions_lead', default).replace('{n}', word(n_pub)))
    if n_ot:
        lead.append(f'{word(n_ot).capitalize()} region{"s" if n_ot != 1 else ""} from the genome-wide significant '
                    f'credible sets of the {s} GWAS in Open Targets' + (' that no published region holds' if n_pub
                                                                         else '') +
                    ', with a credible lead variant, named '
                    'by the top locus-to-gene (L2G) gene of the lead credible set when its score reaches '
                    f'{rp.get("L2gMin", 0.5)}, else by the position of the lead variant; credible sets within '
                    f'{rp.get("RegionWindowKb", 500)} kb are one region.')
    if not lead:
        lead.append(f'This snapshot has no GWAS region of {s}: its profile names no published region, and Open Targets '
                    f'holds no genome-wide significant credible set from a GWAS of {s}'
                    + (' with a credible lead variant.' if len(not_regions) else '.'))
    if len(not_regions):
        lead.append(not_regions_sentence(not_regions, s) + ' They are listed below.')
    if n_ot or len(not_regions):
        lead.append(f'A lead is credible with at least {rp.get("RegionMinMac", 25)} expected copies of its minor allele '
                    f'in the smaller group of its study, and, below a frequency of {rp.get("RegionMinMaf", 0.01)}, '
                    f'P < {rp.get("RareVariantP", 5e-9)}; or when another GWAS of {s} replicates it.')
    curation = (f'Retracted articles and conference material are not part of the graph; findings that are not about '
                f'{s} itself are kept and marked; words that text mining mapped to the wrong gene are removed. No '
                'rule calls a language model.')
    return {'overview_lead': P.text('viewer', 'overview_lead', overview),
            'ladder_note': P.text('viewer', 'ladder_note', note),
            'regions_lead': ' '.join(lead),
            'curation_lead': P.text('viewer', 'curation_lead', curation)}


def condition_evidence(T, G):
    """What the evidence matrix of the conditions page shows of each related condition: the publications of its
    clinical evidence that relate it, the best design class among them (RELATION_DESIGNS, by its index), the
    directions, whether a study found no association, East Asian evidence of the relation, its place in the ontology,
    the regions of the disease it shares a signal with, by method ('c' colocalization, 'o' lead-variant overlap), and
    why its own genetics were not read ('nf': none, shared, broad, cooccurrence, gwas; '' when they were)."""
    rel = read_tsv(T / 'relations.tsv')
    if 'design_class' not in rel:
        rel['design_class'] = ''
    sh = read_tsv(G / 'edges_shares_region.tsv', ['source', 'target', 'method', 'n_signals', 'best', 'measure'])
    out = {}
    for n, g in rel.groupby('node_id'):
        ok = [relates(t, r) for t, r in zip(g.evidence_type, g.relation)]
        clin = [t == 'clinical' for t in g.evidence_type]
        cl = g[[c and o for c, o in zip(clin, ok)]]
        onto = g[g.evidence_type.eq('ontology')].relation.tolist()
        rank = [RELATION_DESIGNS.index(d) for d in cl.design_class if d in RELATION_DESIGNS]
        out[n] = {'pubs': cl.pmid.nunique(), 'design': min(rank) if rank else None,
                  'dirs': sorted(set(cl.relation)), 'contrary': any(c and not o for c, o in zip(clin, ok)),
                  'ea': any(o and str(x) == 'True' for o, x in zip(ok, g.east_asian)),
                  'onto': (onto[0].split(' ')[0] if onto else ''), 'onto_text': '; '.join(onto)}
    # why its own genetics were not read: no Open Targets entry, or (coloc_only) the reason resolve.py wrote
    cond = read_tsv(T / 'conditions.tsv', ['node_id', 'genetics', 'genetics_note'])
    for x in cond.itertuples():
        note = x.genetics_note
        code = ('none' if x.genetics == 'none' else '' if x.genetics != 'coloc_only' else
                'broad' if 'broad class' in note else 'gwas' if 'GWAS studies, more than' in note else
                'cooccurrence' if 'publication' in note else 'shared')
        out.setdefault(x.node_id, {'pubs': 0, 'design': None, 'dirs': [], 'contrary': False, 'ea': False,
                                   'onto': '', 'onto_text': ''})['nf'] = code
    for x in sh.itertuples():
        r = out.setdefault(x.source, {'pubs': 0, 'design': None, 'dirs': [], 'contrary': False, 'ea': False,
                                      'onto': '', 'onto_text': ''}).setdefault('regions', {})
        m = 'c' if x.method == 'colocalization' else 'o'
        r[x.target] = 'c' if r.get(x.target) == 'c' else m
    return out


def model_use(P, T, Q, summ):
    """How the language models were used in this snapshot: the models, the counts of each step that uses them, and
    the disputed items (03.curated/model_review.tsv)."""
    tidy = json.loads((T / 'tidy.json').read_text()) if (T / 'tidy.json').exists() else {}
    ex = tidy.get('extract') or {}
    cur = json.loads((Q / 'curation.json').read_text()) if (Q / 'curation.json').exists() else {}
    rr = cur.get('relations_reading')
    role = {'A': 'Reader A', 'B': 'Reader B', 'adjudicator': 'Adjudicator'}
    models = []
    for r in (ex.get('readers') or []) + ([ex['adjudicator']] if ex.get('adjudicator') else []):
        backend, _, model = str(r.get('backend', '')).partition(' ')
        models.append({'role': role.get(r.get('reader'), r.get('reader', '')), 'backend': backend,
                       'model': model or backend})
    made = summ[~summ.problems.str.startswith('not generated')]
    review = read_tsv(Q / 'model_review.tsv', ['pmid', 'pass', 'item', 'A', 'B', 'final', 'decision'])
    keep = ('records', 'prescreened_out', 'read_by_model', 'not_read', 'offline', 'included', 'findings',
            'methods_fields', 'dropped_items',
            'disputed_items', 'agreement', 'result_flags')
    return {'models': models, 'themes': tidy.get('themes', ''), 'fetch_date': tidy.get('fetch_date', ''),
            'evidence_also': tidy.get('evidence_also') or [], 'not_related': tidy.get('not_related') or [],
            # where the clinical relations come from: the profile says (a curated list, or read from PubMed)
            'curated_relations': (P.d.get('relations') or {}).get('mode') == 'curated',
            'extract': {k: ex[k] for k in keep if k in ex},
            'relations': ({k: rr[k] for k in ('records', 'read', 'not_read', 'relations', 'conditions_read', 'split',
                                              'consolidation', 'conditions_consolidated', 'verification',
                                              'conditions', 'mapped_to_open_targets', 'one_publication',
                                              'agreement', 'disputed_items', 'dropped_items') if k in rr}
                          if rr else None),
            # the traits a shared signal alone would relate, left out: too few regions shared, a broad class
            'shared_left_out': {'too_few': tidy.get('shared_signal_too_few_regions') or {},
                                'broad': tidy.get('shared_signal_broad_classes') or {}},
            'summaries': {'n': len(summ), 'shown': int(summ.grounded.eq('True').sum()),
                          'withheld': int((~made.grounded.eq('True')).sum()), 'not_generated': len(summ) - len(made),
                          'checked': int(summ.check.ne('').sum()) if 'check' in summ else 0},
            'review': review[['pmid', 'pass', 'item', 'A', 'B', 'final', 'decision']].values.tolist()}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--profile', required=True, help='profile.json of the disease')
    ap.add_argument('--graph', required=True, help='05.graph')
    ap.add_argument('--curated', required=True, help='03.curated')
    ap.add_argument('--tables', required=True, help='04.tables')
    ap.add_argument('--summaries', required=True, help='summaries.tsv')
    ap.add_argument('--run-params', required=True, help='the parameters of this run (JSON): the values of the rules')
    ap.add_argument('--out', required=True, help='<key>_knowledge_graph.html')
    a = ap.parse_args()
    P = load_profile(a.profile)
    G, Q, T = Path(a.graph), Path(a.curated), Path(a.tables)
    schema = json.loads((G / 'schema.json').read_text())
    lit = read_tsv(T / 'literature.tsv', ['pmid', 'ensembl_id', 'sentences'])
    first = {(r.pmid, r.ensembl_id): r.sentences.split(' || ')[0] for r in lit.itertuples() if r.sentences}
    nodes, edges = [], []
    for k in schema['nodes']:
        for r in read_tsv(G / f'nodes_{k}.tsv').to_dict('records'):
            nodes.append([r['id'], k, r['label'], r['url'], compact(r, ('id', 'type', 'label', 'url'))])
    for k in schema['edges']:
        for r in read_tsv(G / f'edges_{k}.tsv').to_dict('records'):
            at = compact(r, ('source', 'target', 'type'))
            if k == 'mentions':
                pm = r['source'].split(':', 1)[1]
                if (pm, r['target']) in first:
                    at['sentence'] = first[(pm, r['target'])]
            edges.append([r['source'], r['target'], k, at])
    ladder = {str(k): [v['code'], v['label']] for k, v in sorted(schema['ladder'].items(), key=lambda x: int(x[0]))}
    ladder_def = {str(k): v.get('definition') or v['label'] for k, v in schema['ladder'].items()}
    # every level of the ladder, for every disease; what an empty or partly empty level holds, said in its row
    level_status = {str(k): v for k, v in load_status(T).items() if v.get('brief')}
    ladder_gloss = {str(k): v for k, v in load_glosses(T).items()}
    regions = read_tsv(G / 'nodes_region.tsv', ['id', 'source'])
    not_regions = not_regions_of(T)
    core = schema['graph']['core']
    # the categories of related conditions in the profile's order; a category the profile does not know is appended
    known = P.category_order
    seen = {n[4].get('category', '') for n in nodes if n[1] == 'condition' and n[0] != core}
    cats = known + sorted(c for c in seen if c and c not in known)
    cols = (P.gwas or {}).get('columns') or {}
    pubs = read_tsv(Q / 'publications.tsv')
    cur = json.loads((Q / 'curation.json').read_text())
    lg = read_tsv(Q / 'curation_log.tsv', ['step', 'item', 'rule', 'before', 'after'])
    summ = read_tsv(a.summaries, ['kind', 'grounded', 'title', 'summary', 'problems', 'check'])
    run_params = json.loads(Path(a.run_params).read_text())
    choices, _refs = load_choices(P, run_params, T)
    data = {
        'meta': {**schema['graph'], 'checked_at': cur['checked_at'][:10]},
        'disease': {'id': P.id, 'key': P.key, 'name': P.name, 'Name': P.Name, 'short': P.short,
                    **dict(zip(('mark', 'color'), disease_badge(P.short, P.id)))},
        'about': schema.get('about') or {'default': P.about, 'within': P.within},
        'schema': {'nodes': {k: v['meaning'] for k, v in schema['nodes'].items()},
                   'edges': {k: v['meaning'] for k, v in schema['edges'].items()}},
        'ladder': ladder, 'ladder_def': ladder_def, 'ladder_gloss': ladder_gloss, 'level_status': level_status,
        'no_level': P.fill(NO_LEVEL), 'categories': cats, 'shared_signal': SHARED_SIGNAL, 'other_genes': OTHER_GENES,
        # what each category means (the readers see the clinical ones), and the rule that makes the ontology neighbors
        'category_desc': P.category_desc, 'ontology_category': P.ontology_category,
        'cooccurrence': [c['name'] for c in P.categories if c.get('cooccurrence')],
        'relation_designs': RELATION_DESIGNS, 'cond_evidence': condition_evidence(T, G),
        'prescreen': P.d.get('prescreen') or {},
        'design_classes': [c.get('label', c.get('key', '')) for c in (P.curation.get('design_classes') or [])],
        'ontology': {k: P.ontology.get(k, True) for k in ('parents', 'siblings', 'children')}
                    | {'children_of': len(P.ontology.get('children_of') or [])},
        'eas': list(EAS_COUNTRIES),
        # the labels of the region statistics of a published table that only some diseases have
        'region_labels': {'p_pe': cols.get('p_pe') or 'P in the secondary analysis', 'v': cols.get('v') or 'V'},
        'text': texts(P, regions, ladder, run_params, not_regions), 'category_colors': category_colors(P, cats),
        'not_regions': not_regions.to_dict('records'),
        'llm': model_use(P, T, Q, summ), 'params': run_params,
        'choices': [{'name': g['name'], 'items': [{'key': it['key'], 'item': it['item'], 'value': it['value_text'],
                                                    'basis': it.get('basis', ''), 'refs': [[r[0], r[1]] for r in it['refs']]}
                                                   for it in g['items']]} for g in choices],
        'nodes': nodes, 'edges': edges,
        'curation': {
            'counts': {k: cur[k] for k in ('publications', 'studies', 'studies_excluded', 'textmining_dropped',
                                           'findings_by_subject')},
            'n_excluded': int((pubs.status == 'excluded').sum()),
            'flagged': [compact(r) for r in pubs[pubs.status == 'flagged'][
                ['pmid', 'first_author', 'year', 'title', 'reason', 'roles']].to_dict('records')],
            # the log without the rows of an excluded publication (counted, never listed)
            'log': lg[~lg.after.str.startswith('publication left out') & (lg.step != 'publications')
                      & (lg.rule != 'excluded publication')
                      & ~lg.item.map(lambda x: bool(set(re.findall(r'PMID (\d+)', x))
                                                    & set(pubs[pubs.status == 'excluded'].pmid)))][
                ['step', 'item', 'rule', 'before', 'after']].values.tolist()},
        'summaries': {r.title: r.summary for r in summ[(summ.kind == 'mechanism') & (summ.grounded == 'True')]
                      .itertuples()},
    }
    blob = json.dumps(data, ensure_ascii=False, separators=(',', ':')).replace('<', '\\u003c')
    title = f'{P.short} knowledge graph'.replace('&', '&amp;').replace('<', '&lt;')
    Path(a.out).write_text(TEMPLATE.replace('__DATA__', blob).replace('__TITLE__', title))
    log(f'viewer: {len(nodes)} nodes, {len(edges)} edges, {len(regions)} regions, '
        f'{Path(a.out).stat().st_size / 1e6:.1f} MB')


TEMPLATE = r'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
:root {
  --bg: #ffffff; --panel: #f6f7f9; --line: #d9dde3; --ink: #1f2328; --ink2: #57606a; --link: #1f4e79;
  --accent: #b0413e; --head: #1f3b5f; --chip: #eef2f7; --hi: #fff6d6;
  --t-condition: #b0413e; --t-gene: #3c6e9f; --t-region: #b08a3c; --t-publication: #6a4c93;
  --t-gwas_study: #8a9199; --t-theme: #3f8a6e; --t-pathway: #4f8fbf;
  --l1: #c6dbef; --l2: #e3eef8; --l3: #fce4cc; --l4: #fdf2e6; --l5: #f2f2f2; --navh: 42px;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--ink); font: 14px/1.45 -apple-system, "Segoe UI", Roboto,
       "Helvetica Neue", Arial, "Noto Sans", sans-serif; }
header { background: var(--head); color: #fff; padding: 12px 24px; display: flex; flex-wrap: wrap; gap: 10px 28px;
         align-items: center; }
header .brand { flex: 1 1 520px; min-width: 0; }
header h1 { font-size: 19px; margin: 0; font-weight: 650; letter-spacing: .2px; display: flex; align-items: center; }
/* the disease's mark, the same as on its browser tab */
.mark { display: inline-flex; align-items: center; justify-content: center; flex: none; width: 28px; height: 28px;
        border-radius: 6px; color: #fff; font-size: 12.5px; font-weight: 700; letter-spacing: .3px;
        box-shadow: 0 0 0 1.5px rgba(255,255,255,.55); margin-right: 10px; }
header .meta { font-size: 12.5px; opacity: .9; margin-top: 2px; }
header .chips { display: flex; flex-wrap: wrap; gap: 6px; margin-top: 7px; }
header .chips span { font-size: 11.5px; background: rgba(255,255,255,.13); border: 1px solid rgba(255,255,255,.22);
                     border-radius: 10px; padding: 1px 9px; white-space: nowrap; }
#search-wrap { position: relative; margin-left: auto; min-width: 280px; flex: 0 1 380px; }
#search { width: 100%; padding: 7px 10px; border-radius: 4px; border: 0; font-size: 14px; }
#suggest { position: absolute; top: 36px; left: 0; right: 0; background: #fff; color: var(--ink); z-index: 10;
           border: 1px solid var(--line); box-shadow: 0 4px 14px rgba(0,0,0,.15); display: none; max-height: 420px;
           overflow: auto; }
#suggest div { padding: 6px 10px; cursor: pointer; display: flex; gap: 8px; align-items: center; }
#suggest div.on, #suggest div:hover { background: var(--chip); }
nav { display: flex; gap: 2px; padding: 0 24px; background: var(--panel); border-bottom: 1px solid var(--line);
      overflow-x: auto; position: sticky; top: 0; z-index: 30; height: var(--navh); align-items: stretch; }
nav a { padding: 0 14px; color: var(--ink2); text-decoration: none; border-bottom: 3px solid transparent;
        white-space: nowrap; display: flex; align-items: center; gap: 6px; }
nav a:hover { color: var(--head); background: #eef1f5; text-decoration: none; }
nav a.on { color: var(--head); border-bottom-color: var(--head); font-weight: 600; }
nav a .cnt { font-size: 11px; font-weight: 600; color: var(--ink2); background: #e4e8ee; border-radius: 9px;
             padding: 0 7px; }
nav a.on .cnt { background: var(--head); color: #fff; }
main { padding: 20px 24px 40px; max-width: 1400px; margin: 0 auto; }
footer { border-top: 1px solid var(--line); color: var(--ink2); font-size: 12px; padding: 12px 24px 20px;
         max-width: 1400px; margin: 0 auto; }
h2 { font-size: 20px; margin: 4px 0 6px; }
table.ladder-def { border-collapse: collapse; margin: 6px 0; font-size: 13px; }
table.ladder-def td { padding: 4px 8px; border-bottom: 1px solid var(--line); vertical-align: top; }
table.ladder-def td.lvk { white-space: nowrap; } table.ladder-def td.lvn { font-weight: 600; white-space: nowrap; } h3 { font-size: 15px; margin: 22px 0 8px; color: var(--head); }
p.lead { color: var(--ink2); margin: 0 0 14px; max-width: 980px; }
a { color: var(--link); text-decoration: none; } a:hover { text-decoration: underline; }
a.ext { white-space: nowrap; } a.ext::after { content: " \2197"; font-size: 11px; }
.cards { display: grid; grid-template-columns: repeat(auto-fill, minmax(190px, 1fr)); gap: 10px; margin: 10px 0 6px; }
.card { background: var(--panel); border: 1px solid var(--line); border-radius: 6px; padding: 10px 12px;
        transition: border-color .15s, box-shadow .15s; }
a.card:hover { border-color: var(--link); box-shadow: 0 2px 8px rgba(31,59,95,.12); text-decoration: none; }
.card .n { font-size: 22px; font-weight: 650; color: var(--head); } .card .t { font-size: 12px; color: var(--ink2); }
table { border-collapse: collapse; width: 100%; font-size: 13px; }
th, td { border-bottom: 1px solid var(--line); padding: 5px 7px; text-align: left; vertical-align: top; }
th { background: var(--panel); font-weight: 600; position: sticky; top: var(--navh); cursor: pointer; white-space: nowrap;
     z-index: 2; box-shadow: inset 0 -1px 0 var(--line); }
th.nosort { cursor: default; }
tbody tr:nth-child(even) td { background: #fbfcfd; }
tr:hover td { background: #f1f5fa; }
td.num { text-align: right; font-variant-numeric: tabular-nums; }
.quote { color: var(--ink2); font-style: italic; }
.badge { display: inline-block; font-size: 11px; padding: 1px 7px; border-radius: 10px; color: #fff;
         white-space: nowrap; }
.badge.adj { background: #a8792e; font-size: 10px; padding: 0 6px; }
.badge.both { background: #4f8a68; font-size: 10px; padding: 0 6px; }
.badge.one { background: #8a9199; font-size: 10px; padding: 0 6px; }
.dot { display: inline-block; width: 9px; height: 9px; border-radius: 50%; flex: none; }
.gene { font-style: italic; }
.chip { display: inline-block; background: #fff; border: 1px solid var(--line); border-radius: 3px; padding: 0 5px;
        margin: 2px 3px 2px 0; font-style: italic; text-decoration: none; color: var(--ink); }
.chip:hover { border-color: var(--link); color: var(--link); }
.grid td { min-width: 80px; } .grid th { cursor: default; }
/* a table in a horizontally scrolled box keeps its header in place (a sticky header would sit below the navigation
   inside the box); its first column stays visible while the box scrolls (top: auto, or its header cell would be
   pushed under the first row); the function names wrap, so that more columns fit the window */
.grid th, .matrix th { position: static; top: auto; box-shadow: none; }
.grid td:first-child, .grid th:first-child { position: sticky; left: 0; z-index: 1; }
.grid th:first-child { background: var(--panel); }
.grid thead th { white-space: normal; vertical-align: bottom; }
.hscroll-hint { margin: 0 0 4px; text-align: right; }
.grid td:first-child { background: inherit; box-shadow: inset -1px 0 0 var(--line); }
.grid tbody tr:nth-child(even) td { background: inherit; }
.lv1 { background: var(--l1); } .lv2 { background: var(--l2); } .lv3 { background: var(--l3); }
.lv4 { background: var(--l4); } .lv5, .lv6, .lv7 { background: var(--l5); }
.matrix td.c { text-align: center; font-size: 15px; } .matrix th.rot { writing-mode: vertical-rl;
  transform: rotate(180deg); height: 96px; font-weight: 600; font-style: italic; }
.two { display: grid; grid-template-columns: minmax(300px, 1fr) minmax(320px, 560px); gap: 22px; align-items: start; }
@media (max-width: 900px) { .two { grid-template-columns: 1fr; } }
.facts { display: grid; grid-template-columns: max-content 1fr; gap: 4px 14px; font-size: 13px; }
.facts .k { color: var(--ink2); }
.tools { display: flex; gap: 10px; align-items: center; flex-wrap: wrap; margin: 8px 0; }
.tools input[type=search] { padding: 5px 8px; border: 1px solid var(--line); border-radius: 4px; min-width: 240px; }
.tools select { padding: 4px 6px; max-width: 340px; }
.grid td:first-child { min-width: 150px; }
.btn { border: 1px solid var(--line); background: #fff; border-radius: 4px; padding: 3px 9px; cursor: pointer;
       font-size: 12px; }
.btn.on { background: var(--head); color: #fff; border-color: var(--head); }
.note { font-size: 12px; color: var(--ink2); margin-top: 6px; }
svg text { font: 11px sans-serif; fill: var(--ink); }
.legend { display: flex; flex-wrap: wrap; gap: 4px 14px; font-size: 12px; color: var(--ink2); margin: 6px 0; }
.legend span { display: inline-flex; gap: 5px; align-items: center; }
.summary { background: var(--panel); border-left: 3px solid var(--t-theme); padding: 8px 12px; margin: 8px 0;
           font-size: 13px; }
.cmap-tip { position: absolute; pointer-events: none; background: #fff; border: 1px solid var(--line);
            box-shadow: 0 4px 14px rgba(0,0,0,.15); border-radius: 5px; padding: 7px 10px; font-size: 12.5px;
            max-width: 320px; display: none; z-index: 5; }
.cmap-tip b { color: var(--head); }
.cmap-tip .k { color: var(--ink2); }
.cmap-cats { display: flex; flex-wrap: wrap; gap: 6px; justify-content: center; margin: 4px 0; }
.cmap-cats .btn { display: inline-flex; align-items: center; gap: 6px; }
.cmap-cats .btn.on .dot { box-shadow: 0 0 0 2px #fff; }
.cmap-key { display: flex; flex-wrap: wrap; gap: 4px 18px; justify-content: center; font-size: 12px;
            color: var(--ink2); margin: 2px 0 6px; }
.cmap-key span { display: inline-flex; gap: 6px; align-items: center; }
.cmap-all { margin: 10px 0 4px; }
.cmap-all summary { cursor: pointer; color: var(--head); font-weight: 600; }
/* the evidence matrix of the related conditions: counts shaded by bins, the number in the cell */
table.mx { width: auto; min-width: 100%; font-size: 12.5px; }
table.mx th, table.mx td { padding: 3px 6px; border-bottom: 1px solid #eef0f3; }
table.mx tr.mx-g th { text-align: center; font-size: 11.5px; color: var(--ink2); border-bottom: 1px solid var(--line); }
table.mx th.mx-s { cursor: pointer; } table.mx th.mx-s.on { color: var(--link); text-decoration: underline; }
table.mx th { cursor: default; vertical-align: bottom; }
table.mx th.mx-r div { writing-mode: vertical-rl; transform: rotate(180deg); font-style: italic; max-height: 96px;
  white-space: nowrap; margin: 0 auto; }
table.mx .mx-name { position: sticky; left: 0; z-index: 1; background: #fff; min-width: 220px; max-width: 320px;
  white-space: normal; box-shadow: inset -1px 0 0 var(--line); }
table.mx th.mx-name { background: var(--panel); }
table.mx tr.mx-cat td { background: var(--panel); padding-top: 9px; font-size: 13px; }
table.mx tr.mx-cat .dot { margin-right: 6px; }
.mx-desc { color: var(--ink2); font-size: 12px; }
table.mx td.mx-n { text-align: center; min-width: 40px; font-variant-numeric: tabular-nums; }
table.mx td.q1, .mx-sw.q1 { background: #e3edf7; } table.mx td.q2, .mx-sw.q2 { background: #bcd3ea; }
table.mx td.q3, .mx-sw.q3 { background: #6f9fcf; color: #fff; } table.mx td.q4, .mx-sw.q4 { background: #2f5f8f;
  color: #fff; }
table.mx td.g1, .mx-sw.g1 { background: #f1f3f5; } table.mx td.g2, .mx-sw.g2 { background: #dfe3e8; }
table.mx td.g3, .mx-sw.g3 { background: #c3cad3; } table.mx td.g4, .mx-sw.g4 { background: #9aa5b1; color: #fff; }
table.mx td.mx-c { text-align: center; } table.mx td.mx-d, table.mx td.mx-o { color: var(--ink2); white-space: nowrap; }
table.mx td.mx-dir { white-space: nowrap; } table.mx td.mx-dir span { margin-right: 4px; font-weight: 700; }
.mx-x { color: var(--accent); } .mx-ea { color: #b08a3c; }
.mx-dot { display: inline-block; width: 11px; height: 11px; border-radius: 50%; border: 2px solid #3c6e9f;
  vertical-align: -1px; } .mx-dot.c { background: #3c6e9f; } .mx-dot.o { background: #fff; }
table.mx td.mx-nf { color: var(--ink2); font-style: italic; text-align: center; font-size: 12px; }
table.mx tbody tr:hover td.mx-name { background: var(--hi); }
.mx-key { display: flex; flex-wrap: wrap; gap: 6px 18px; font-size: 12px; color: var(--ink2); margin: 6px 0 2px; }
.mx-sw { display: inline-block; padding: 0 6px; margin: 0 2px; border-radius: 3px; }
.guide { font-size: 13px; max-width: 1100px; }
/* the build page: a methods document with a contents column */
.doc { display: grid; grid-template-columns: 250px minmax(0, 1fr); gap: 30px; margin-top: 18px; align-items: start; }
.doc-nav { position: sticky; top: calc(var(--navh) + 14px); max-height: calc(100vh - var(--navh) - 28px);
  overflow: auto; font-size: 13px; border-right: 1px solid var(--line); padding-right: 12px; }
.doc-nav-h { font-weight: 700; color: var(--head); margin-bottom: 6px; text-transform: uppercase; font-size: 11.5px;
  letter-spacing: .04em; }
.doc-nav ol { list-style: none; margin: 0; padding: 0; }
.doc-nav a { display: flex; gap: 8px; align-items: baseline; text-decoration: none; border-radius: 4px;
  line-height: 1.35; }
.doc-nav a:hover { background: var(--panel); color: var(--link); }
/* a section: bold, dark, its number aligned; its subsections nested under it, smaller and lighter, joined by a line */
.doc-nav .toc-l1 > li { margin: 0 0 8px; }
.doc-nav .toc-l1 > li > a { padding: 4px 6px; font-weight: 700; font-size: 13.5px; color: var(--head); }
.doc-nav .toc-l1 > li > a .doc-no { min-width: 20px; }
.doc-nav .toc-l2 { margin: 2px 0 4px 14px; padding-left: 8px; border-left: 2px solid var(--line); }
.doc-nav .toc-l2 a { padding: 2px 6px; font-size: 12.5px; font-weight: 400; color: var(--ink2); }
.doc-nav .toc-l2 .doc-no { min-width: 24px; color: #8a9199; font-weight: 400; }
.doc-nav a.on { background: #e3ecf6; color: var(--head); }
.doc-nav .toc-l2 a.on { font-weight: 600; } .doc-nav .toc-l2 li:has(a.on) { margin-left: -10px;
  border-left: 2px solid var(--head); padding-left: 8px; }
.doc-nav a.in { color: var(--head); text-decoration: underline; text-decoration-color: var(--line);
  text-underline-offset: 3px; }
.doc-no { display: inline-block; min-width: 30px; color: var(--head); font-variant-numeric: tabular-nums; }
.doc-body { min-width: 0; max-width: 1160px; }
.doc-sec { padding-bottom: 18px; margin-bottom: 18px; border-bottom: 1px solid var(--line); }
.doc-sec > h3 { font-size: 19px; margin: 6px 0 10px; } .doc-sec > h3 .doc-no { min-width: 38px; }
.doc-sub { margin: 18px 0 8px; } .doc-sub > h4 { font-size: 15px; color: var(--head); margin: 0 0 6px; }
.doc-sec, .doc-sub { scroll-margin-top: calc(var(--navh) + 14px); }
.insh2 { font-size: 13.5px; color: var(--ink); margin: 0 0 10px; font-style: italic; }
.doc-k { font-size: 11.5px; font-weight: 700; color: var(--ink2); text-transform: uppercase; letter-spacing: .04em;
  margin: 14px 0 4px; }
.doc-list { padding-left: 20px; margin: 6px 0 4px; font-size: 13px; } .doc-list li { margin: 4px 0; }
table.rules td.rule-k { width: 24%; font-weight: 600; } table.rules td.rule-v { width: 22%; }
table.rules td { vertical-align: top; } .rule-refs { margin-top: 3px; font-size: 12px; color: var(--ink2); }
table.proc td { vertical-align: top; } table.proc td.proc-no { width: 22px; color: var(--ink2); }
table.proc td.proc-step { width: 13%; font-weight: 600; } table.proc td.proc-by { width: 12%; }
table.proc td.proc-here { width: 22%; color: var(--ink); font-size: 12.5px; }
@media (max-width: 1000px) { .doc { grid-template-columns: 1fr; } .doc-nav { position: static; max-height: none;
  border-right: 0; padding-right: 0; } }
.toc { display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 6px 22px; margin: 6px 0 16px;
  padding: 10px 14px; background: var(--panel); border: 1px solid var(--line); border-radius: 6px; font-size: 13px; }
.toc .toc-h { grid-column: 1 / -1; font-weight: 700; color: var(--head); }
.toc a { font-weight: 600; } .toc .toc-s { display: block; color: var(--ink2); font-size: 12px; }
.toc .toc-q a { font-weight: 500; } .toc .toc-h + .toc-q, .toc .toc-q + .toc-h { margin-top: 0; }
.toc .toc-h:not(:first-child) { margin-top: 8px; padding-top: 8px; border-top: 1px solid var(--line); }
.insh { font-size: 14px; line-height: 1.5; color: var(--ink); border-left: 3px solid var(--head); padding: 2px 0 2px 10px;
  margin: 2px 0 12px; max-width: 1080px; }
.fb { display: flex; flex-wrap: wrap; align-items: stretch; gap: 6px; margin: 6px 0 14px; }
.fb-box { border: 1px solid var(--line); border-radius: 6px; padding: 7px 10px; min-width: 120px; max-width: 190px;
  background: #fff; }
.fb-n { font-size: 20px; font-weight: 700; color: var(--head); font-variant-numeric: tabular-nums; }
.fb-u { font-size: 12px; font-weight: 400; color: var(--ink2); }
.fb-l { font-size: 12px; line-height: 1.35; margin: 2px 0 4px; } .fb-by .badge { font-size: 10px; margin-right: 2px; }
.fb-arrow { align-self: center; color: var(--ink2); font-size: 16px; }
/* the matrix scrolls in its own box: its header and its first column stay in place */
.mx-scroll { max-height: calc(100vh - 170px); overflow: auto; border: 1px solid var(--line); border-radius: 4px; }
table.mx thead { position: sticky; top: 0; z-index: 3; }
table.mx thead th { background: var(--panel); }
table.mx th.mx-name { z-index: 4; }
.mx-catl { position: sticky; left: 8px; }
.mx-notes { display: grid; gap: 3px; font-size: 12.5px; color: var(--ink2); margin: 6px 0 4px; max-width: 1180px; }
.mx-notes b { color: var(--ink); }
.mx-tools { display: flex; align-items: center; gap: 12px; margin: 4px 0 6px; } .guide h4 { color: var(--head); font-size: 14px; margin: 16px 0 6px; }
.guide ul { margin: 4px 0; padding-left: 20px; } .guide li { margin: 3px 0; } .guide .stps { margin-top: 6px; }
/* the page on how the graph is built: at a glance, step by step, in detail */
.glance { display: grid; grid-template-columns: 1fr 34px 1fr 34px 1fr; align-items: stretch; margin: 14px 0 8px; }
.gl-col { background: #fff; border: 1px solid var(--line); border-radius: 8px; padding: 10px 12px 10px;
          box-shadow: 0 1px 3px rgba(0,0,0,.05); }
.gl-h { font-size: 12px; font-weight: 700; letter-spacing: .6px; text-transform: uppercase; color: var(--head);
        border-bottom: 2px solid var(--head); padding-bottom: 5px; margin-bottom: 8px; }
.gl-item { margin: 0 0 9px; }
.gl-t { font-size: 12.5px; font-weight: 600; color: var(--ink); }
.gl-t .badge { font-size: 11px; }
.gl-d { font-size: 12.5px; color: var(--ink2); margin-top: 2px; line-height: 1.4; }
.gl-n { display: grid; grid-template-columns: 14px 1fr auto; align-items: center; gap: 6px; font-size: 13px;
        padding: 3px 0; border-bottom: 1px dotted var(--line); }
.gl-n b { font-variant-numeric: tabular-nums; color: var(--head); }
.gl-foot { font-size: 12px; color: var(--ink2); margin-top: 8px; padding-top: 6px; border-top: 1px solid var(--line); }
.gl-arrow { display: flex; align-items: center; justify-content: center; font-size: 26px; color: var(--ink2); }
@media (max-width: 1000px) { .glance { grid-template-columns: 1fr; } .gl-arrow { transform: rotate(90deg);
  height: 30px; } }
.stps { position: relative; max-width: 1100px; margin: 6px 0 4px; }
.stps::before { content: ""; position: absolute; left: 15px; top: 8px; bottom: 8px; width: 2px; background: var(--line); }
.stp { display: grid; grid-template-columns: 32px 1fr; gap: 12px; margin: 0 0 12px; position: relative; }
.stp-no { width: 32px; height: 32px; border-radius: 50%; background: var(--head); color: #fff; font-weight: 700;
          display: flex; align-items: center; justify-content: center; font-size: 14px; z-index: 1; }
.stp-b { font-size: 13px; padding-top: 4px; }
.stp-h { margin-bottom: 2px; } .stp-h b { font-size: 14px; color: var(--head); margin-right: 4px; }
.stp-h .badge { font-size: 10.5px; margin: 0 2px 2px 0; }
details.bd { border: 1px solid var(--line); border-radius: 6px; margin: 8px 0; background: #fff; }
details.bd > summary { cursor: pointer; padding: 9px 14px; background: var(--panel); border-radius: 6px;
                       list-style-position: inside; }
details.bd[open] > summary { border-bottom: 1px solid var(--line); border-radius: 6px 6px 0 0; }
details.bd > summary b { color: var(--head); font-size: 14px; }
.bd-sub { color: var(--ink2); font-size: 12.5px; margin-left: 10px; }
.bd-body { padding: 8px 14px 12px; }
table.plain th { cursor: default; position: static; }
table.plain td { font-size: 13px; }
.stages { padding-left: 22px; max-width: 1080px; font-size: 13px; }
.stages li { margin: 0 0 8px; }
.stage-h { margin-bottom: 1px; }
.stage-h .badge, table.plain .badge { margin: 0 3px 2px 0; font-size: 10.5px; }
.llm-step h4, #main > h4 { color: var(--head); }
#main > h4 { font-size: 14px; margin: 18px 0 6px; }
.flow-wrap { max-width: 1130px; margin: 6px 0; }
.flow { display: block; width: 100%; height: auto; }
.swatch { display: inline-block; width: 22px; height: 12px; border: 1.4px solid transparent; border-radius: 3px; }
.llm-rules { max-width: 980px; font-size: 13px; }
.llm-rules p { margin: 4px 0 6px; }
.llm-steps { display: grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap: 12px; }
.llm-step { background: var(--panel); border: 1px solid var(--line); border-radius: 6px; padding: 4px 14px 10px;
            font-size: 13px; }
.llm-step h4 { margin: 8px 0 4px; color: var(--head); font-size: 14px; }
.llm-step p { margin: 3px 0; }
.here { margin-top: 5px; padding-left: 8px; border-left: 3px solid var(--line); color: var(--ink); font-size: 13px; }
.agree { max-width: 980px; }
.agree-row { display: grid; grid-template-columns: 170px minmax(120px, 1fr) 290px; gap: 10px; align-items: center;
             margin: 5px 0; font-size: 13px; }
.agree-bar { display: flex; height: 16px; background: var(--panel); border-radius: 3px; overflow: hidden; }
.agree-bar .both { background: #5b8c6e; } .agree-bar .adj { background: #c9a25b; }
.agree-v { color: var(--ink2); }
@media (max-width: 700px) { .agree-row { grid-template-columns: 1fr; } .llm-steps { grid-template-columns: 1fr; } }
.pager { display: flex; flex-wrap: wrap; align-items: center; gap: 8px 16px; margin: 8px 0 4px; font-size: 12.5px;
         color: var(--ink2); }
.pg-links { display: flex; flex-wrap: wrap; align-items: center; gap: 4px; }
.pg { border: 1px solid var(--line); background: #fff; color: var(--link); border-radius: 4px; padding: 2px 9px;
      min-width: 30px; cursor: pointer; font-size: 12.5px; }
.pg:hover:not([disabled]) { border-color: var(--link); }
.pg.on { background: var(--head); border-color: var(--head); color: #fff; font-weight: 600; cursor: default; }
.pg[disabled]:not(.on) { color: #b0b6bd; cursor: default; }
.pager .gap { padding: 0 2px; }
.pg-size { margin-left: auto; white-space: nowrap; }
.pg-size select { margin-left: 4px; padding: 1px 4px; font-size: 12.5px; }
.tbl { scroll-margin-top: calc(var(--navh) + 8px); }
/* print, or save as PDF: the current page as a document */
@media print {
  @page { margin: 12mm; }
  body { font-size: 9.5pt; }
  nav, #search-wrap, #suggest, .tools, .btn, .cmap-tip, .hscroll-hint, .pager { display: none !important; }
  th { top: 0; }
  header { background: none; color: var(--ink); padding: 0 0 4px; border-bottom: 2px solid var(--head); }
  header h1 { color: var(--head); }
  header .meta { opacity: 1; color: var(--ink2); }
  .cmap-all { display: block; } .cmap-all > summary { list-style: none; }
  .mx-scroll { max-height: none; overflow: visible; border: 0; }
  .doc { display: block; } .doc-nav { display: none; }
  .cmap-cats .btn { display: inline-flex !important; border: 0; padding: 0 4px; }
  main { padding: 8px 0 0; max-width: none; }
  table { font-size: 8.5pt; }
  th { position: static; background: var(--panel); white-space: normal; }
  thead { display: table-header-group; break-after: avoid; page-break-after: avoid; }
  tbody tr:first-child { break-before: avoid; page-break-before: avoid; }
  section.short { break-inside: avoid; page-break-inside: avoid; }
  tr, .card, .facts, .summary, svg { break-inside: avoid; page-break-inside: avoid; }
  h2, h3 { break-after: avoid; page-break-after: avoid; }
  td, th { overflow-wrap: break-word; }
  td.num { white-space: nowrap; }
  .grid td { min-width: 0; }
  .grid td:first-child { min-width: 30mm; }
  .chip { margin: 1px 2px 1px 0; padding: 0 3px; }
  div[style*="overflow-x"] { overflow: visible !important; }
  .two { display: block; }
  .two svg { display: block; width: 15cm; max-width: 100%; margin: 6px auto; }
  .cards { grid-template-columns: repeat(3, 1fr); }
  a { color: inherit; text-decoration: none; }
  a.ext::after { content: ""; }
  .facts a.ext::after { content: " <" attr(href) ">"; font-size: 8pt; color: var(--ink2); }
  .llm-step, .agree-row { break-inside: avoid; page-break-inside: avoid; }
  .lv1, .lv2, .lv3, .lv4, .lv5, .lv6, .lv7, .badge, .dot, .card, th, svg, .agree-bar span, .swatch {
    -webkit-print-color-adjust: exact; print-color-adjust: exact; }
}
</style>
</head>
<body>
<header>
  <div class="brand"><h1><span class="mark" id="mark"></span>__TITLE__</h1><div class="meta" id="meta"></div><div class="chips" id="chips"></div></div>
  <div id="search-wrap"><input id="search" type="search" placeholder="Search a gene, condition, PMID, region…"
       autocomplete="off"><div id="suggest"></div></div>
</header>
<nav id="nav"></nav>
<main id="main"></main>
<footer id="foot"></footer>
<script type="application/json" id="kg-data">__DATA__</script>
<script>
'use strict';
const D = JSON.parse(document.getElementById('kg-data').textContent);
const N = new Map(), OUT = new Map(), IN = new Map();
for (const [id, type, label, url, a] of D.nodes) N.set(id, {id, type, label, url, a: a || {}});
for (const [s, t, type, a] of D.edges) {
  const e = {s, t, type, a: a || {}};
  if (!OUT.has(s)) OUT.set(s, []); OUT.get(s).push(e);
  if (!IN.has(t)) IN.set(t, []); IN.get(t).push(e);
}
const CORE = D.meta.core, S = D.disease.short, ABOUT = D.about.default;
const TYPE_NAME = {condition: 'Condition', gene: 'Gene', region: `${S} GWAS region`, publication: 'Publication',
                   gwas_study: 'GWAS study', theme: 'Function', pathway: 'Reactome pathway'};
const COLOR = t => getComputedStyle(document.documentElement).getPropertyValue('--t-' + t).trim() || '#888';
// the categories of related conditions in the profile's order (unknown ones appended by build_viewer.py)
const CATS = D.categories;
const catRank = c => { const i = CATS.indexOf(c); return i < 0 ? CATS.length : i; };
// the published regions come from the profile's table; the others from the disease's credible sets in Open Targets
const published = r => !r.a.source || r.a.source === 'published';

// ── small helpers ─────────────────────────────────────────────────────────
function h(tag, attrs, ...kids) {
  const el = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v === null || v === undefined || v === false) continue;
    if (k === 'class') el.className = v; else if (k === 'html') el.innerHTML = v;
    else if (k.startsWith('on')) el.addEventListener(k.slice(2), v); else el.setAttribute(k, v);
  }
  for (const k of kids.flat()) if (k !== null && k !== undefined && k !== '')
    el.append(k instanceof Node ? k : document.createTextNode(String(k)));
  return el;
}
const href = id => '#/n/' + encodeURIComponent(id);
// a table wider than the window scrolls sideways in its box; a hint above the box says so while it does
const hscroll = (tbl, what) => [h('div', {class: 'note hscroll-hint', style: 'display:none'},
  `Wider than the window: scroll sideways to see all ${what} →`), h('div', {class: 'hscroll', style: 'overflow-x:auto'}, tbl)];
const fitHints = () => document.querySelectorAll('.hscroll').forEach(w => {
  w.previousElementSibling.style.display = w.scrollWidth > w.clientWidth + 1 ? '' : 'none'; });
// a link from a list (a table, the ladder, the map) opens the node in a new tab, so that the list, its filter and
// its place are kept; a link on a node's own page follows the graph in the same tab (the browser's Back returns)
let LIST_PAGE = true;
const newTab = () => LIST_PAGE ? {target: '_blank', rel: 'noopener'} : {};
const nodeLink = (id, text) => {
  const n = N.get(id);
  if (!n) return text || id;
  return h('a', {href: href(id), class: n.type === 'gene' ? 'gene' : null, ...newTab(),
                 title: TYPE_NAME[n.type] + ': ' + n.id + (LIST_PAGE ? ' (opens in a new tab)' : '')},
           text || n.label);
};
const openNode = id => { if (LIST_PAGE) window.open(href(id), '_blank', 'noopener'); else location.hash = href(id); };
const ext = (url, text) => url ? h('a', {href: url, target: '_blank', rel: 'noopener', class: 'ext'}, text) : text;
const pmid = p => {
  if (!p) return '';
  const id = String(p);
  return /^\d+$/.test(id) ? ext('https://pubmed.ncbi.nlm.nih.gov/' + id + '/', 'PMID ' + id)
                          : ext('https://europepmc.org/search?query=' + id, id);
};
const rsid = r => r ? ext('https://www.ncbi.nlm.nih.gov/snp/' + r, r) : '';
const badge = t => h('span', {class: 'badge', style: 'background:' + COLOR(t)}, TYPE_NAME[t] || t);
const dot = t => h('span', {class: 'dot', style: 'background:' + COLOR(t)});
const fmt = v => {
  if (v === undefined || v === null || v === '') return '';
  const x = Number(v);
  if (typeof v === 'string' && /^-?\d+(\.\d+)?(e-?\d+)?$/i.test(v) && !Number.isInteger(x)) {
    return Math.abs(x) < 1e-3 ? x.toExponential(1) : String(+x.toPrecision(3));
  }
  return String(v);
};
const num = v => v === undefined || v === '' ? -Infinity : Number(v);
const lvl = n => n.a.disease_level !== undefined ? Number(n.a.disease_level) : null;
// the agreement of the two model readers on a record: 'adjudicated' as a badge, the other values as text
const agree = v => v === 'adjudicated' ? h('span', {class: 'badge adj', title: 'the readers disagreed; an ' +
  'adjudicator decided'}, 'adjudicated') : v === 'both' ? h('span', {class: 'badge both', title: 'both readers gave ' +
  'the same answer'}, 'both') : v === 'one reader' ? h('span', {class: 'badge one'}, 'one reader') : (v || '');

// the rows per page of every long table: one choice for the whole viewer, kept in this browser (links open new tabs,
// which keep it); 0 shows every row. Printing shows every row.
const PAGE_SIZES = [10, 25, 50, 100, 0];
let PAGE_SIZE = 25, PRINTING = false;
try { const v = localStorage.getItem('diseaselens.pageSize');
      if (v !== null && PAGE_SIZES.includes(Number(v))) PAGE_SIZE = Number(v); }
catch (e) { /* no storage: the default */ }
const PAGED = [];   // the redraw of each table on the page

// a sortable, filterable, paged table: cols [{label, get(row) -> value for sort/filter, cell(row) -> node}]
function table(cols, rows, opts = {}) {
  let sortBy = opts.sort ?? null, dir = opts.dir ?? 1, q = '', page = 0;
  const wrap = h('div', {class: 'tbl'});
  const input = opts.filter === false || rows.length <= 8 ? null :
    h('input', {type: 'search', placeholder: 'Filter these rows…', oninput: e => { q = e.target.value.toLowerCase();
      page = 0; draw(); }});
  const count = h('span', {class: 'note'});
  if (input || opts.tools) wrap.append(h('div', {class: 'tools'}, input, opts.tools || null, count));
  const tbl = h('table'), body = h('tbody'), pager = h('div', {class: 'pager'});
  const head = h('tr');
  cols.forEach((c, i) => head.append(h('th', {class: c.get ? null : 'nosort', onclick: () => {
    if (!c.get) return; dir = sortBy === i ? -dir : (c.num ? -1 : 1); sortBy = i; page = 0; draw(); }}, c.label)));
  tbl.append(h('thead', {}, head), body); wrap.append(tbl, pager);
  // a new page starts at the top of the table, below the navigation
  const go = p => { page = p; draw();
    const top = wrap.getBoundingClientRect().top, nav = document.getElementById('nav').offsetHeight;
    if (top < nav) window.scrollBy(0, top - nav - 8); };
  function draw() {
    let rs = rows;
    if (q) rs = rs.filter(r => cols.some(c => c.get && String(c.get(r) ?? '').toLowerCase().includes(q)));
    if (sortBy !== null) {
      const c = cols[sortBy];
      rs = [...rs].sort((a, b) => {
        const x = c.get(a), y = c.get(b);
        return (c.num ? (num(x) - num(y)) : String(x ?? '').localeCompare(String(y ?? ''))) * dir;
      });
    }
    const size = PRINTING || !PAGE_SIZE ? Math.max(rs.length, 1) : PAGE_SIZE;
    const pages = Math.max(1, Math.ceil(rs.length / size));
    page = Math.min(page, pages - 1);
    const from = page * size, shown = rs.slice(from, from + size);
    body.replaceChildren(...shown.map(r => h('tr', {}, cols.map(c =>
      h('td', {class: c.num ? 'num' : null}, c.cell ? c.cell(r) : c.num ? fmt(c.get(r)) : linkIds(fmt(c.get(r))))))));
    count.textContent = rs.length === rows.length ? `${rows.length} rows` : `${rs.length} of ${rows.length} rows`;
    // the pager: only for a table longer than a page (or than the default page, when every row is shown)
    if (PRINTING || rs.length <= (PAGE_SIZE || 25)) { pager.replaceChildren(); return; }
    const btn = (text, p, attrs = {}) => h('button', {class: 'pg' + (p === page ? ' on' : ''), disabled:
      p < 0 || p >= pages || p === page ? '' : null, onclick: () => go(p), ...attrs}, text);
    const nums = [...new Set([0, page - 1, page, page + 1, pages - 1])].filter(p => p >= 0 && p < pages)
      .sort((a, b) => a - b);
    const links = [];
    nums.forEach((p, i) => { if (i && p - nums[i - 1] > 1) links.push(h('span', {class: 'gap'}, '…'));
      links.push(btn(String(p + 1), p, {'aria-label': `Page ${p + 1}`})); });
    const sel = h('select', {'aria-label': 'Rows per page', onchange: e => {
      const first = page * (PAGE_SIZE || rs.length);
      PAGE_SIZE = Number(e.target.value);
      try { localStorage.setItem('diseaselens.pageSize', String(PAGE_SIZE)); } catch (err) { /* kept for this page */ }
      page = PAGE_SIZE ? Math.floor(first / PAGE_SIZE) : 0;
      PAGED.forEach(f => f()); }},
      PAGE_SIZES.map(n => h('option', {value: n, selected: n === PAGE_SIZE ? '' : null}, n ? String(n) : 'All')));
    pager.replaceChildren(
      h('span', {class: 'pg-range'}, `${(from + 1).toLocaleString()}–${(from + shown.length).toLocaleString()} of ` +
        `${rs.length.toLocaleString()}`),
      h('span', {class: 'pg-links'}, pages > 1 ? [btn('‹ Previous', page - 1), ...links, btn('Next ›', page + 1)] : []),
      h('label', {class: 'pg-size'}, 'Rows per page ', sel));
  }
  PAGED.push(draw);
  draw();
  return wrap;
}

// ── the neighborhood diagram ─────────────────────────────────────────────
function diagram(id) {
  const edges = [...(OUT.get(id) || []).map(e => [e, e.t]), ...(IN.get(id) || []).map(e => [e, e.s])];
  const seen = new Map();
  for (const [e, o] of edges) if (o !== id && !seen.has(o)) seen.set(o, e.type);
  let others = [...seen.entries()];
  const total = others.length, MAX = 48;
  const order = ['region', 'gene', 'condition', 'publication', 'theme', 'pathway', 'gwas_study'];
  others.sort((a, b) => order.indexOf(N.get(a[0])?.type) - order.indexOf(N.get(b[0])?.type) ||
                        String(N.get(a[0])?.label).localeCompare(String(N.get(b[0])?.label)));
  if (others.length > MAX) {
    const keep = [], perType = {};
    for (const o of others) { const t = N.get(o[0])?.type; perType[t] = (perType[t] || 0) + 1;
      if (perType[t] <= Math.max(6, Math.floor(MAX / 4))) keep.push(o); }
    others = keep.slice(0, MAX);
  }
  const W = 680, H = 460, cx = W / 2, cy = H / 2, R = 160;
  const ns = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(ns, 'svg');
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`); svg.setAttribute('width', '100%');
  const el = (tag, at) => { const x = document.createElementNS(ns, tag);
    for (const [k, v] of Object.entries(at)) x.setAttribute(k, v); return x; };
  const me = N.get(id);
  // labels are shortened to fit; hovering a node (its dot or its label) shows its full name, kind and identifier, a
  // publication's title, and the edge that links it to this node
  const tip = h('div', {class: 'cmap-tip', role: 'tooltip'}), wrap = h('div', {style: 'position:relative'}, svg, tip);
  const place = e => {
    const r = wrap.getBoundingClientRect();
    let x = e.clientX - r.left + 14;
    if (x + tip.offsetWidth > r.width) x = Math.max(0, e.clientX - r.left - tip.offsetWidth - 14);
    tip.style.left = x + 'px'; tip.style.top = (e.clientY - r.top + 14) + 'px';
  };
  const hover = (target, n, via) => {
    const kids = () => [h('b', {class: n.type === 'gene' ? 'gene' : null}, n.label),
      h('div', {class: 'k'}, `${TYPE_NAME[n.type] || n.type} · ${n.id}`),
      n.a.title && n.a.title !== n.label ? h('div', {}, n.a.title) : null,
      via ? h('div', {class: 'k', style: 'margin-top:3px'}, `${via.replace(/_/g, ' ')} — ${D.schema.edges[via] || ''}`)
        : null].filter(Boolean);
    target.addEventListener('mouseenter', e => { tip.replaceChildren(...kids()); tip.style.display = 'block'; place(e); });
    target.addEventListener('mousemove', place);
    target.addEventListener('mouseleave', () => { tip.style.display = 'none'; });
  };
  const labeled = others.length <= 36, nR = Math.ceil(others.length / 2), nL = others.length - nR;
  others.forEach(([o, via], i) => {
    // labeled: right half top to bottom, then left half bottom to top, evenly spaced in height, so labels never
    // collide; unlabeled: a full circle
    const right = i < nR, m = right ? nR : nL, k = right ? i : nL - 1 - (i - nR);
    const ang = -Math.PI / 2 + 2 * Math.PI * i / others.length;
    const dy = labeled ? R * 0.94 * (2 * (k + 0.5) / m - 1) : R * Math.sin(ang);
    const dx = labeled ? (right ? 1 : -1) * Math.sqrt(R * R - dy * dy) : R * Math.cos(ang);
    const x = cx + dx, y = cy + dy, n = N.get(o);
    if (!n) return;
    svg.append(el('line', {x1: cx, y1: cy, x2: x, y2: y, stroke: '#c9ced6', 'stroke-width': 1}));
    const g = el('a', {href: href(o), 'aria-label': `${TYPE_NAME[n.type]}: ${n.label}`});
    g.append(el('circle', {cx: x, cy: y, r: 7, fill: COLOR(n.type), stroke: '#fff', 'stroke-width': 1.5}));
    hover(g, n, via);
    if (labeled) {
      const lab = n.label.length > 26 ? n.label.slice(0, 25) + '…' : n.label;
      const t = el('text', {x: x + (right ? 10 : -10), y: y + 4, 'text-anchor': right ? 'start' : 'end'});
      t.textContent = lab; if (n.type === 'gene') t.setAttribute('font-style', 'italic'); g.append(t);
    }
    svg.append(g);
  });
  const center = el('g', {'aria-label': `${TYPE_NAME[me.type]}: ${me.label}`});
  center.append(el('circle', {cx, cy, r: 13, fill: COLOR(me.type), stroke: '#fff', 'stroke-width': 2}));
  const t = el('text', {x: cx, y: cy + 30, 'text-anchor': 'middle', 'font-weight': 700});
  t.textContent = me.label.length > 40 ? me.label.slice(0, 39) + '…' : me.label; center.append(t);
  hover(center, me, null); svg.append(center);
  const types = [...new Set(others.map(([o]) => N.get(o)?.type).filter(Boolean))];
  return h('div', {}, wrap, h('div', {class: 'legend'}, types.map(x => h('span', {}, dot(x), TYPE_NAME[x]))),
           h('div', {class: 'note'}, total > others.length ? `${others.length} of ${total} neighbors shown; every ` +
             'neighbor is listed in the tables.' : `${total} neighbors. Click one to open it.`));
}

// ── edge sections of a node page ─────────────────────────────────────────
const ALL_PUBLISHED = [...N.values()].filter(n => n.type === 'region').every(published);
const SECTION = {
  'related_to:out': `Why it relates to ${S}`, 'related_to:in': 'Related conditions, one row per piece of evidence',
  'supports:in': 'Publications that give its clinical evidence', 'supports:out': 'Clinical evidence it gives',
  'has_genetic_evidence:out': 'Genetic evidence in conditions',
  'has_genetic_evidence:in': 'Genes with genetic evidence',
  'reports:in': `Findings of the genetic studies of ${S}`, 'reports:out': 'Findings it reports',
  'studies:out': 'Genetic study of', 'studies:in': `Genetic studies of ${S}`,
  'region_of:out': 'Region of', 'region_of:in': `${ALL_PUBLISHED ? 'Published ' : ''}${S} GWAS regions`,
  'names_gene:out': 'Gene named', 'names_gene:in': `${S} GWAS region`,
  'shares_region:out': `${S} GWAS regions it shares`,
  'shares_region:in': 'Conditions sharing this region', 'in_theme:out': 'Function (gene theme)',
  'in_theme:in': 'Genes of this function', 'in_pathway:out': 'Reactome pathways',
  'in_pathway:in': `Text-mined ${S} targets in this pathway`, 'mentions:in': `Text-mined ${S} publications`,
  'mentions:out': 'Genes it mentions (text mining)', 'gwas_of:in': 'GWAS studies', 'gwas_of:out': 'GWAS of'};
const ORDER = ['names_gene:in', 'region_of:in', 'reports:in', 'reports:out', 'studies:out', 'related_to:out',
               'related_to:in', 'shares_region:out', 'shares_region:in', 'names_gene:out', 'has_genetic_evidence:out',
               'has_genetic_evidence:in', 'supports:out', 'supports:in', 'in_theme:out', 'in_theme:in',
               'in_pathway:out', 'in_pathway:in', 'studies:in', 'mentions:in', 'mentions:out', 'gwas_of:in',
               'gwas_of:out', 'region_of:out'];
const A = (key, label, kind) => ({key, label, kind});
const EDGE_COLS = {
  related_to: [A('evidence_type', 'Evidence'), A('relation', 'Relation'), A('detail', 'Statement'),
               A('quote', 'Quote', 'quote'), A('pmid', 'Source', 'pmid'), A('strength', 'Strength', 'num'),
               A('east_asian', 'East Asian', 'bool'), A('agreement', 'Readers', 'agree')],
  supports: [A('relation', 'Relation'), A('detail', 'Statement')],
  has_genetic_evidence: [A('code', 'Code'), A('kind', 'Kind'), A('n', 'Records', 'num'), A('details', 'Details')],
  reports: [A('result', 'Result'), A('about', 'About'), A('variant_name', 'Variant'), A('rsid', 'rsID', 'rsid'),
            A('effect', 'Effect'), A('p_value', 'P'), A('gene_set', 'Gene set'), A('quote', 'Source sentence', 'quote'),
            A('agreement', 'Readers', 'agree')],
  studies: [A('design_class', 'Design class'), A('n_cases', 'Cases', 'num'), A('n_controls', 'Controls', 'num'),
            A('place', 'Place'), A('agreement', 'Readers', 'agree')],
  region_of: [A('code', 'Code'), A('tier', 'Tier'), A('p', 'P', 'num')], names_gene: [],
  shares_region: [A('method', 'Method'), A('n_signals', 'GWAS signals', 'num'), A('best', 'Best', 'num'),
                  A('measure', 'Measure')],
  in_theme: [], in_pathway: [],
  mentions: [A('matched_as', 'Matched as'), A('sentence', 'Sentence', 'quote'), A('n_sentences', 'Sentences', 'num'),
             A('flag', 'Flag')],
  gwas_of: []};
// columns read from the node at the other end of an edge
const OTHER_COLS = {
  publication: [A('year', 'Year', 'num'), A('place', 'Place'), A('design_class', 'Design class')],
  gene: [A('disease_level', `${S} level`, 'num')], condition: [A('category', 'Category')],
  gwas_study: [A('n_samples', 'Samples', 'num'), A('ancestry', 'Ancestry')]};

// the identifiers in a text, as links (the text itself unchanged): a GWAS study (GCST..., FINNGEN_..., UKB_PPP_...)
// to its page in the graph, else in Open Targets; an rsID to dbSNP; a PMID to its page in the graph, else PubMed
const ID_RX = /\b(GCST\d+|FINNGEN_R\d+_[A-Za-z0-9_]+|UKB_PPP_[A-Za-z0-9_]+)\b|\b(rs\d{3,})\b|\bPMID:?\s?(\d{4,9})\b/g;
function linkIds(text) {
  const s = String(text ?? '');
  if (!/GCST\d|FINNGEN_R|UKB_PPP_|rs\d|PMID/.test(s)) return s;
  const out = [];
  let last = 0, m;
  ID_RX.lastIndex = 0;
  while ((m = ID_RX.exec(s))) {
    out.push(s.slice(last, m.index));
    const [url, inside] = m[1] ? (N.has(m[1]) ? [href(m[1]), true] : ['https://platform.opentargets.org/study/' + m[1]])
      : m[2] ? ['https://www.ncbi.nlm.nih.gov/snp/' + m[2]]
      : N.has('PMID:' + m[3]) ? [href('PMID:' + m[3]), true] : ['https://pubmed.ncbi.nlm.nih.gov/' + m[3] + '/'];
    out.push(inside ? h('a', {href: url}, m[0]) : h('a', {href: url, target: '_blank', rel: 'noopener'}, m[0]));
    last = m.index + m[0].length;
  }
  out.push(s.slice(last));
  return h('span', {}, out);
}

// a flag as words: True / False -> yes / no
const yesNo = v => v === 'True' || v === true ? 'yes' : v === 'False' || v === false ? 'no' : v;

function cellOf(kind, v) {
  if (v === undefined || v === '') return '';
  if (kind === 'pmid') return pmid(v);
  if (kind === 'rsid') return rsid(v);
  if (kind === 'quote') return h('span', {class: 'quote'}, linkIds(v));
  if (kind === 'agree') return agree(v);
  if (kind === 'bool') return yesNo(v);
  return linkIds(fmt(v));
}

function edgeSection(id, key, list) {
  const [type, dir] = key.split(':');
  const otherOf = e => dir === 'out' ? e.t : e.s;
  const sample = N.get(otherOf(list[0]));
  const ocols = (OTHER_COLS[sample?.type] || []).filter(c => list.some(e => N.get(otherOf(e))?.a[c.key] !== undefined));
  const ecols = (EDGE_COLS[type] || []).filter(c => list.some(e => e.a[c.key] !== undefined));
  const single = list.length > 1 && new Set(list.map(otherOf)).size === 1;
  const cols = (single ? [] : [{label: sample ? TYPE_NAME[sample.type] : 'Node', get: e => N.get(otherOf(e))?.label,
                 cell: e => nodeLink(otherOf(e))}]
    .concat(ocols.map(c => ({label: c.label, num: c.kind === 'num', get: e => N.get(otherOf(e))?.a[c.key],
                             cell: e => fmt(N.get(otherOf(e))?.a[c.key])}))))
    .concat(ecols.map(c => ({label: c.label, num: c.kind === 'num', get: e => e.a[c.key],
                             cell: e => cellOf(c.kind, e.a[c.key])})));
  // a short section is kept on one printed page (heading, note and rows together)
  return h('section', {class: list.length <= 10 ? 'short' : null},
           h('h3', {}, `${SECTION[key] || key} (${list.length})`),
           h('div', {class: 'note'}, D.schema.edges[type] || ''),
           single ? h('div', {class: 'note'}, 'Every row: ', nodeLink(otherOf(list[0]))) : null,
           table(cols, list));
}

// ── facts of a node ───────────────────────────────────────────────────────
const FACTS = {
  gene: [['disease_level', `${S} evidence level`, v => `${v} — ${(D.ladder[v] || [])[1] || ''}` +
          ((D.ladder[v] || [])[0] ? ` (${D.ladder[v][0]})` : '')], ['theme', 'Function'],
         ['disease_level_note', 'Why this level'], ['n_related_conditions', 'Related conditions with genetic evidence'],
         ['n_text_mined_publications', `Text-mined ${S} publications`],
         ['disease_ot_score', `Open Targets association score with ${S}`], ['pathways', 'Reactome pathways']],
  condition: [['category', 'Category'], ['sources', 'Evidence types'], ['n_clinical', 'Clinical evidence items'],
              ['shared_regions', `${S} regions shared`], ['n_gwas_studies', 'GWAS studies'],
              ['n_eas_studies', 'East Asian GWAS studies'], ['n_credible_sets', 'Credible sets'],
              ['n_mendelian_genes', 'Curated genes'], ['genetics', 'Genetic data fetched'],
              ['genetics_note', 'Why not in full'], ['ot_ids', 'Ontology IDs']],
  // a statistic of the published table is shown only where the region has it (Open Targets regions have none)
  region: [['rsid', 'Lead variant', rsid], ['chrom', 'Chromosome'], ['pos_grch37', 'Position (GRCh37)'],
           ['pos_grch38', 'Position (GRCh38)'], ['tier', 'Tier'], ['code', 'Code'], ['p', 'P'],
           ['p_pe', D.region_labels.p_pe], ['v', D.region_labels.v], ['or', 'OR'], ['maf', 'MAF'],
           ['source', 'Source'], ['open_targets', 'In Open Targets']],
  publication: [['title', 'Title'], ['first_author', 'First author'], ['year', 'Year'], ['journal', 'Journal'],
                ['roles', 'Role in the graph'], ['place', "Place (first author's institution)"],
                ['design_class', 'Design class'], ['n_cases', 'Cases'], ['n_controls', 'Controls'],
                ['agreement', 'Readers', agree], ['status', 'Flag']],
  gwas_study: [['trait', 'Trait'], ['first_author', 'First author'], ['year', 'Year'], ['n_samples', 'Samples'],
               ['n_credible_sets', 'Credible sets'], ['ancestry', 'Ancestry (LD reference)'], ['eas', 'East Asian', yesNo],
               ['pmid', 'Publication', pmid]],
  theme: [['description', 'Description']], pathway: []};

function nodePage(id) {
  const n = N.get(id);
  if (!n) return h('div', {}, h('h2', {}, 'Not found'), h('p', {}, id));
  const facts = h('div', {class: 'facts'});
  for (const [k, label, f] of FACTS[n.type] || []) {
    const v = n.a[k]; if (v === undefined || v === '') continue;
    facts.append(h('div', {class: 'k'}, label), h('div', {}, f ? f(v) : linkIds(fmt(v))));
  }
  facts.append(h('div', {class: 'k'}, 'Identifier'), h('div', {}, n.url ? ext(n.url, n.id) : n.id));
  const left = h('div', {}, h('div', {class: 'tools'}, badge(n.type)),
                 h('h2', {class: n.type === 'gene' ? 'gene' : null}, n.label), facts);
  if (n.type === 'pathway' && D.summaries[n.label])
    left.append(h('h3', {}, 'Generated summary'), h('div', {class: 'summary'}, D.summaries[n.label]),
                h('div', {class: 'note'}, 'Written by a language model from the sentences of the listed publications ' +
                  'and checked against them: every cited PMID is among its sources.'));
  const groups = new Map();
  for (const e of OUT.get(id) || []) { const k = e.type + ':out'; if (!groups.has(k)) groups.set(k, []);
    groups.get(k).push(e); }
  for (const e of IN.get(id) || []) { const k = e.type + ':in'; if (!groups.has(k)) groups.set(k, []);
    groups.get(k).push(e); }
  const keys = [...groups.keys()].sort((a, b) => ORDER.indexOf(a) - ORDER.indexOf(b));
  return h('div', {}, h('div', {class: 'two'}, left, diagram(id)),
           keys.map(k => edgeSection(id, k, groups.get(k))));
}

// ── list pages ───────────────────────────────────────────────────────────
const byType = t => [...N.values()].filter(n => n.type === t);
const themes = () => byType('theme').map(t => t.label).concat([D.other_genes]);
const short = label => (byType('theme').find(t => t.label === label) || {a: {}}).a.short || label;
// the related conditions: every condition node but the disease itself
const related = () => byType('condition').filter(c => c.id !== CORE);
// '{value} {count}, ...' of a list of values ('' counted as 'not recorded')
const tally = vs => Object.entries(vs.reduce((m, v) => { m[v || 'not recorded'] = (m[v || 'not recorded'] || 0) + 1;
  return m; }, {})).map(([k, n]) => `${k} ${n}`).join(', ');

function overview() {
  const m = D.meta, c = D.curation.counts, genes = byType('gene').filter(n => lvl(n) !== null);
  const pubs = byType('publication').filter(p => p.a.design_class), regions = byType('region');
  const card = (n, t, link) => h('a', {href: link || null, class: 'card', style: 'text-decoration:none;color:inherit'},
                                  h('div', {class: 'n'}, n), h('div', {class: 't'}, t));
  const out = [h('h2', {}, 'Overview'), h('p', {class: 'lead'}, D.text.overview_lead),
    h('div', {class: 'cards'},
      card(genes.length, `${S} genes on the evidence ladder`, '#/genes'),
      card(pubs.length, `genetic studies of ${S} (${pubs.filter(p => p.a.japanese === 'True').length} in Japan)`,
           '#/studies'),
      card(related().length, 'related conditions', '#/conditions'),
      card(regions.length, `${regions.length && ALL_PUBLISHED ? 'published ' : ''}${S} GWAS regions`, '#/regions'),
      card(`${c.publications.toLocaleString()}`, 'publications checked in PubMed', '#/curation'),
      card(`${m.n_nodes.toLocaleString()} / ${m.n_edges.toLocaleString()}`, 'nodes / edges')),
    h('h3', {}, `Evidence ladder: how strong is the genetic evidence for each ${S} gene, and what does it do?`)];
  const th = themes();
  const grid = h('table', {class: 'grid'});
  grid.append(h('thead', {}, h('tr', {}, h('th', {}, 'Level'), th.map(t => h('th', {title: t}, short(t))))));
  const tb = h('tbody');
  for (const [k, [code, label]] of Object.entries(D.ladder)) {
    const row = h('tr', {class: 'lv' + k}, h('td', {title: D.ladder_def[k]}, h('b', {}, `${k} ${code}`), h('br'), label,
      D.ladder_gloss[k] ? h('div', {class: 'note', style: 'margin:2px 0 0'}, D.ladder_gloss[k]) : null,
      D.level_status[k] ? h('div', {class: 'note', style: 'margin:2px 0 0', title: D.level_status[k].full},
        D.level_status[k].brief) : null));
    for (const t of th) {
      const gs = genes.filter(g => lvl(g) === Number(k) && (g.a.theme || D.other_genes) === t)
                      .sort((a, b) => a.label.localeCompare(b.label));
      row.append(h('td', {}, gs.map(g => h('a', {class: 'chip', href: href(g.id), ...newTab(),
        title: 'Gene ' + g.label + ' (opens in a new tab)'}, g.label))));
    }
    tb.append(row);
  }
  grid.append(tb);
  // what each level means: the same definitions for every disease
  const defs = h('table', {class: 'ladder-def'}, h('tbody', {}, Object.entries(D.ladder).map(([k, [code, label]]) =>
    h('tr', {class: 'lv' + k}, h('td', {class: 'lvk'}, h('b', {}, `${k}`), code ? ` ${code}` : ''),
      h('td', {class: 'lvn'}, label), h('td', {}, D.ladder_def[k],
        D.level_status[k] ? h('div', {class: 'note'}, `${S}: ${D.level_status[k].full}`) : null)))));
  out.push(...hscroll(grid, `${th.length} function columns`),
    h('details', {class: 'cmap-all'}, h('summary', {}, 'What each level means'), defs),
    h('div', {class: 'note'}, D.text.ladder_note));
  out.push(h('h3', {}, `Which conditions share the ${S} GWAS regions?`), regionMatrix());
  return out;
}

function regionMatrix() {
  const regions = byType('region');
  if (!regions.length) {
    return h('div', {class: 'note'}, `This snapshot has no GWAS region of ${S}, so no condition is compared with ` +
             'one (GWAS regions page).');
  }
  const nShared = c => new Set((OUT.get(c.id) || []).filter(e => e.type === 'shares_region').map(e => e.t)).size;
  const sharing = related().filter(c => nShared(c) >= 1);
  const rows = sharing.filter(c => nShared(c) >= 2 || c.a.category !== D.shared_signal)
    .map(c => ({c, s: new Map((OUT.get(c.id) || []).filter(e => e.type === 'shares_region').map(e => [e.t, e]))}))
    .sort((a, b) => b.s.size - a.s.size || catRank(a.c.a.category) - catRank(b.c.a.category) ||
                    a.c.label.localeCompare(b.c.label));
  // the traits related only by a shared signal that share one region: not rows, but counted
  const single = sharing.length - rows.length;
  const singleNote = single === 1 ? '1 trait related only by a shared GWAS signal shares one region and is not a ' +
    'row here.' : single ? `${single} traits related only by a shared GWAS signal share one region each and are not ` +
    'rows here.' : '';
  if (!rows.length) {
    return h('div', {class: 'note'}, sharing.length ? `No related condition shares two or more GWAS regions of ${S}, ` +
      `and no clinically related condition shares one. ${singleNote}` : `No related condition shares a GWAS region ` +
      `of ${S}.`);
  }
  const tbl = h('table', {class: 'matrix'});
  tbl.append(h('thead', {}, h('tr', {}, h('th', {}, 'Condition'), h('th', {}, 'Category'),
    regions.map(r => h('th', {class: 'rot', title: [r.a.code, r.a.tier ? `tier ${r.a.tier}` : '',
                                                    published(r) ? '' : r.a.source].filter(Boolean).join(' · ')},
                        r.label.replace(' region', ''))),
    h('th', {}, 'Regions'))));
  tbl.append(h('tbody', {}, rows.map(({c, s}) => h('tr', {}, h('td', {}, nodeLink(c.id)), h('td', {}, c.a.category),
    regions.map(r => { const e = s.get(r.id);
      return h('td', {class: 'c', title: e ? `${e.a.method}: ${e.a.n_signals} GWAS signals, best ${e.a.measure} ` +
        fmt(e.a.best) : ''}, e ? (e.a.method === 'colocalization' ? '●' : '○') : ''); }),
    h('td', {class: 'num'}, s.size)))));
  return h('div', {}, ...hscroll(tbl, `${regions.length} region columns`),
           h('div', {class: 'note'}, `● a GWAS signal of the condition colocalizes with the ${S} credible set; ○ a ` +
             `GWAS credible set of the condition contains the ${S} lead variant (regions without a ${S} credible ` +
             'set; weaker). Rows: conditions sharing two or more regions, and clinically related conditions sharing ' +
             'one.' + (singleNote ? ' ' + singleNote : '')));
}

function genesPage() {
  let all = false;
  const box = h('div');
  const draw = () => {
    const rows = byType('gene').filter(g => all || lvl(g) !== null);
    box.replaceChildren(table([
      {label: 'Gene', get: g => g.label, cell: g => nodeLink(g.id)},
      {label: `${S} level`, num: true, get: g => g.a.disease_level, cell: g => g.a.disease_level === undefined ? '' :
        h('span', {class: 'lv' + g.a.disease_level, style: 'padding:0 6px'},
          `${g.a.disease_level} ${g.a.disease_code || ''}`)},
      {label: 'Evidence', get: g => lvl(g) === null ? `No level: text-mined in ${g.a.n_text_mined_publications || 0} ` +
        `${S} publications; genetic evidence in ${g.a.n_related_conditions || 0} related conditions` :
        g.a.disease_level_label + (g.a.disease_level_note ? ` (${g.a.disease_level_note})` : '')},
      {label: 'Function', get: g => g.a.theme},
      {label: 'Related conditions', num: true, get: g => g.a.n_related_conditions},
      {label: 'Text-mined publications', num: true, get: g => g.a.n_text_mined_publications}],
      rows, {sort: all ? 4 : 1, dir: all ? -1 : 1, tools: h('span', {},
        h('button', {class: 'btn' + (all ? '' : ' on'), onclick: () => { all = false; draw(); }}, `${S} genes`), ' ',
        h('button', {class: 'btn' + (all ? ' on' : ''), onclick: () => { all = true; draw(); }}, 'Every gene'))}));
  };
  draw();
  return [h('h2', {}, 'Genes'), h('p', {class: 'lead'}, `The genes studied genetically in ${S}, by evidence level ` +
    `(levels 1–7). Level 7 holds the genes a genetic study of ${S} names without a result about ${S} itself. ` +
    `“Every gene” adds the genes without a level, which appear in no genetic study of ${S}: the text-mined ${S} ` +
    'targets (Europe PMC, as Open Targets gives them) and the genes with genetic evidence in a related condition.'), box];
}

// ── SVG helpers ──
const SVGNS = 'http://www.w3.org/2000/svg';
function svgEl(tag, at, ...kids) {
  const x = document.createElementNS(SVGNS, tag);
  for (const [k, v] of Object.entries(at || {})) if (v !== null && v !== undefined) x.setAttribute(k, v);
  for (const k of kids.flat()) if (k !== null && k !== undefined) x.append(k);
  return x;
}

// the evidence matrix: one row per related condition, grouped by category, one column per kind of evidence; counts
// are shaded by bins (1, 2–3, 4–9, 10 or more) so that a condition with a hundred publications does not hide those
// with one, and the number is in the cell
const BIN = (v, q = 'q') => !v ? '' : q + (v === 1 ? 1 : v <= 3 ? 2 : v <= 9 ? 3 : 4);
const DESIGN_SHORT = ['guideline', 'meta-analysis', 'cohort / case-control', 'cross-sectional', 'case report',
                      'other'];
const DIR_MARK = {'increased risk': ['↑', 'increased risk'], 'decreased risk': ['↓', 'decreased risk'],
                  'associated': ['~', 'associated']};
const catColor = c => (D.category_colors || {})[c] || '#8a9199';
function evidenceMatrix(focus, sortBy, setSort, query, own) {
  const E = D.cond_evidence || {};
  const conds = related().filter(c => !query || c.label.toLowerCase().includes(query) ||
    c.id.toLowerCase().includes(query));
  const regions = byType('region'), perRegion = regions.length && regions.length <= 20;
  const ev = c => E[c.id] || {pubs: 0, dirs: [], regions: {}};
  const nReg = c => Object.keys(ev(c).regions || {}).length;
  const num = v => v === '' || v === undefined || v === null ? null : Number(v);
  const KEYS = {
    pubs: c => -ev(c).pubs, design: c => ev(c).design ?? 9, regions: c => -nReg(c),
    gwas: c => -(num(c.a.n_gwas_studies) ?? -1), genes: c => -(num(c.a.n_mendelian_genes) ?? -1)};
  if (!own && (sortBy === 'gwas' || sortBy === 'genes')) sortBy = 'pubs';
  const order = (a, b) => (KEYS[sortBy](a) - KEYS[sortBy](b)) || (ev(b).pubs - ev(a).pubs) ||
    ((ev(a).design ?? 9) - (ev(b).design ?? 9)) || (nReg(b) - nReg(a)) || a.label.localeCompare(b.label);
  const shown = CATS.filter(k => !focus || focus === k);
  // the header: column groups, then columns (a sortable column says so)
  const sortTh = (key, label, title) => h('th', {class: 'mx-s' + (sortBy === key ? ' on' : ''), title:
    (title ? title + '. ' : '') + 'Click to sort each category by this column.', onclick: () => setSort(key)}, label);
  const regTh = perRegion ? regions.map(r => h('th', {class: 'mx-r', title: r.label}, h('div', {},
    r.label.replace(/ region$/, '')))) : [sortTh('regions', `${S} regions`, `How many ${S} GWAS regions it shares`)];
  const head = h('thead', {},
    h('tr', {class: 'mx-g'}, h('th', {}), h('th', {colspan: 4}, 'Clinical evidence'), h('th', {}, 'Ontology'),
      regions.length ? h('th', {colspan: perRegion ? regions.length : 1}, `Shared ${S} GWAS regions`) : null,
      own ? h('th', {colspan: 2}, 'Its own genetics') : null),
    h('tr', {}, h('th', {class: 'mx-name'}, 'Condition'),
      sortTh('pubs', 'Publications', 'Publications whose clinical evidence relates it'),
      sortTh('design', 'Best design', 'The best design among those publications'),
      h('th', {title: 'Direction of the relation; ✕: a study found no association'}, 'Direction'),
      h('th', {title: 'East Asian evidence of the relation'}, 'EA'),
      h('th', {}, 'Relation'), ...(regions.length ? regTh : []),
      own ? sortTh('gwas', 'GWAS', 'Its GWAS studies in Open Targets') : null,
      own ? sortTh('genes', 'Curated genes', 'Genes with curated gene–disease evidence for it') : null));
  const nCols = (own ? 8 : 6) + (regions.length ? (perRegion ? regions.length : 1) : 0);
  const body = h('tbody');
  for (const k of shown) {
    const rows = conds.filter(c => c.a.category === k).sort(order);
    if (!rows.length) continue;
    body.append(h('tr', {class: 'mx-cat'}, h('td', {colspan: nCols}, h('span', {class: 'mx-catl'},
      h('span', {class: 'dot', style: 'background:' + catColor(k)}), h('b', {}, `${SHOWN_CAT(k)} (${rows.length})`),
      D.category_desc[k] ? h('span', {class: 'mx-desc'}, ' · ' + D.category_desc[k]) : null))));
    for (const c of rows) {
      const e = ev(c), g = c.a.genetics || 'none';
      // clinical counts in blue, the condition's own genetics in gray: the matrix is about the relation
      const cnt = (v, what, q, more) => v ? h('td', {class: 'mx-n ' + BIN(v, q), title: `${v} ${what}` +
        (more || '')}, String(v)) : h('td', {class: 'mx-n'});
      const dirs = Object.keys(DIR_MARK).filter(d => (e.dirs || []).includes(d)).map(d => h('span',
        {title: DIR_MARK[d][1]}, DIR_MARK[d][0]));
      if (e.contrary) dirs.push(h('span', {class: 'mx-x', title: 'a study found no association'}, '✕'));
      const regCells = !regions.length ? [] : perRegion ? regions.map(r => {
        const m = (e.regions || {})[r.id];
        return h('td', {class: 'mx-c', title: m ? `${r.label}: ${m === 'c' ? 'colocalization' :
          'lead-variant overlap'}` : ''}, m ? h('span', {class: 'mx-dot ' + m}) : null);
      }) : [cnt(nReg(c), `${S} regions shared`)];
      const ea = num(c.a.n_eas_studies);
      const ownCells = g === 'full' ? [cnt(num(c.a.n_gwas_studies), 'GWAS studies', 'g', ea ? `, ${ea} East Asian` : ''),
                                  cnt(num(c.a.n_mendelian_genes), 'curated genes', 'g')] :
        [h('td', {colspan: 2, class: 'mx-nf', title: c.a.genetics_note || ''},
           g === 'none' ? 'no Open Targets entry' : 'not fetched')];
      body.append(h('tr', {}, h('td', {class: 'mx-name'}, nodeLink(c.id)),
        cnt(e.pubs, 'publications'),
        h('td', {class: 'mx-d', title: e.design !== null && e.design !== undefined ?
          D.relation_designs[e.design] : ''}, e.design !== null && e.design !== undefined ?
          DESIGN_SHORT[e.design] : ''),
        h('td', {class: 'mx-dir'}, dirs), h('td', {class: 'mx-c'}, e.ea ? h('span', {class: 'mx-ea',
          title: 'East Asian evidence of the relation'}, '◆') : null),
        h('td', {class: 'mx-o', title: e.onto_text || ''}, e.onto || ''), ...regCells, ...(own ? ownCells : [])));
    }
  }
  const key = h('div', {class: 'mx-key'},
    h('span', {}, 'Publications: ', ...['q1', 'q2', 'q3', 'q4'].map((q, i) => h('span', {class: 'mx-sw ' + q},
      ['1', '2–3', '4–9', '≥ 10'][i]))),
    own ? h('span', {}, 'Its own genetics: ', ...['g1', 'g2', 'g3', 'g4'].map((q, i) => h('span',
      {class: 'mx-sw ' + q}, ['1', '2–3', '4–9', '≥ 10'][i]))) : null,
    h('span', {}, '↑ increased risk · ↓ decreased risk · ~ associated · ', h('span', {class: 'mx-x'}, '✕'),
      ' a study found no association'),
    h('span', {}, h('span', {class: 'mx-ea'}, '◆'), ' East Asian evidence'),
    regions.length ? h('span', {}, h('span', {class: 'mx-dot c'}), ' colocalization ', h('span', {class: 'mx-dot o'}),
      ' lead-variant overlap') : null,
    h('span', {}, (own ? 'Hover a GWAS count for its East Asian studies, “not fetched” for the reason; ' : '') +
      'click a column title to sort each category by it.'));
  if (!body.childNodes.length) return h('div', {class: 'note'}, 'No related condition matches.');
  // under the matrix: what a blank count of its own genetics means, with this disease's counts
  const P = D.params, nf = count(related(), c => (E[c.id] || {}).nf || 'full');
  const why = {shared: 'related only by a shared GWAS signal',
               broad: `a broad class (more than ${P.MaxDescendants ?? 100} subtypes in the ontology)`,
               cooccurrence: `related only as ${D.cooccurrence.join(' or ').toLowerCase()} by fewer than ` +
                 `${P.ComorbidityMinPublications ?? 2} publications`,
               gwas: `more than ${P.MaxGwasStudies ?? 100} GWAS studies (a general trait)`};
  const reasons = Object.keys(why).filter(k => nf[k]).map(k => `${why[k]} (${nf[k]})`);
  const notes = h('div', {class: 'mx-notes'},
    h('div', {}, h('b', {}, 'No Open Targets entry'), ' — the condition has no entry in Open Targets, so it has no ' +
      `genetics of its own to read (${nf.none || 0} conditions).`),
    h('div', {}, h('b', {}, 'Not fetched'), ` — its relation and its signals shared with ${S} stand, but its own GWAS ` +
      `studies and curated genes are not read: ${reasons.join('; ') || 'no condition here'}.`),
    h('div', {}, own ? 'Blank counts are not zeros. ' : 'Its own genetics are hidden: “Show its own genetics” ' +
      'above the matrix. ', 'Why, and every rule of this page: ', h('a', {href: '#/build/relations'},
      'How it is built › Related conditions'), '.'));
  return h('div', {}, h('div', {class: 'mx-scroll'}, h('table', {class: 'matrix mx'}, head, body)), key, notes);
}

function conditionsPage(selected) {
  // the selection is part of the address (#/conditions/<category>), so that it can be shared
  const conds = related(), cats = CATS.filter(k => conds.some(c => c.a.category === k));
  let focus = selected && cats.includes(selected) ? selected : null, sortBy = 'pubs', query = '';
  // the condition's own genetics (GWAS studies, curated genes): hidden unless the viewer asks, remembered here
  let own = false;
  try { own = localStorage.getItem('diseaselens.ownGenetics') === '1'; } catch (e) { /* no storage: hidden */ }
  const ownBtn = h('button', {class: 'btn'});
  const btnBox = h('div', {class: 'cmap-cats'}), matBox = h('div'), shown = h('span', {class: 'note'});
  const pick = f => {
    focus = focus === f ? null : f;
    try { history.replaceState(null, '', '#/conditions' + (focus ? '/' + encodeURIComponent(focus) : '')); }
    catch (e) { /* an address that cannot change: the selection still applies */ }
    draw();
  };
  const setSort = k => { sortBy = k; draw(); };
  const draw = () => {
    btnBox.replaceChildren(h('button', {class: 'btn' + (!focus ? ' on' : ''), onclick: () => focus && pick(focus)},
      `All (${conds.length})`), ...cats.map(k => h('button', {class: 'btn' + (focus === k ? ' on' : ''),
      onclick: () => pick(k)}, h('span', {class: 'dot', style: 'background:' + catColor(k)}),
      `${SHOWN_CAT(k)} (${conds.filter(c => c.a.category === k).length})`)));
    const n = conds.filter(c => (!focus || c.a.category === focus) && (!query ||
      c.label.toLowerCase().includes(query) || c.id.toLowerCase().includes(query))).length;
    shown.textContent = n === conds.length ? '' : `${n} of ${conds.length} conditions shown`;
    ownBtn.className = 'btn' + (own ? ' on' : '');
    ownBtn.textContent = own ? 'Hide its own genetics' : 'Show its own genetics (GWAS studies, curated genes)';
    matBox.replaceChildren(evidenceMatrix(focus, sortBy, setSort, query, own));
  };
  ownBtn.addEventListener('click', () => {
    own = !own;
    try { localStorage.setItem('diseaselens.ownGenetics', own ? '1' : '0'); } catch (e) { /* not remembered */ }
    draw();
  });
  const filter = h('input', {type: 'search', placeholder: 'Filter conditions…', oninput: e => {
    query = e.target.value.trim().toLowerCase(); draw(); }});
  draw();
  return [h('h2', {}, 'Related conditions'), h('p', {class: 'lead'}, `Every condition related to ${S}, with the ` +
    'kinds of evidence for the relation: one row per condition, grouped by category, one column per kind of ' +
    'evidence. Select a category or filter by name; click a condition to open its page, with every piece of its ' +
    'evidence. How the conditions are found and classified, and what each column means: ',
    h('a', {href: '#/build/relations'}, 'How it is built › Related conditions'), '.'),
    btnBox, h('div', {class: 'tools mx-tools'}, filter, ownBtn, shown), matBox];
}

// tables of the build page: a plain table; a procedure (one row per step: what happens, who acts, the result here)
const plainTable = (head, rows) => h('table', {class: 'plain'},
  h('thead', {}, h('tr', {}, head.map(x => h('th', {class: 'nosort'}, x)))),
  h('tbody', {}, rows.map(r => h('tr', {}, r.map(x => h('td', {}, x))))));
const procTable = rows => h('table', {class: 'plain proc'},
  h('thead', {}, h('tr', {}, ['', 'Step', 'By', 'What happens', S].map(x => h('th', {class: 'nosort'}, x)))),
  h('tbody', {}, rows.filter(Boolean).map(([no, step, by, what, mine]) => h('tr', {}, h('td', {class: 'proc-no'}, no),
    h('td', {class: 'proc-step'}, step), h('td', {class: 'proc-by'}, by.map(tag)), h('td', {}, what),
    h('td', {class: 'proc-here'}, mine || '')))));

// the related conditions, part by part (sections 4.1–4.7 of the build page): the routes into the graph, the clinical
// relations step by step with the models' part, the categories, the evidence matrix, the genetics of the conditions
// and what is left out; the rules are those of every disease, the counts this snapshot's
function relationsParts() {
  const L = D.llm, R = L.relations, P = D.params, conds = related(), cur = L.curated_relations;
  const n = v => (v ?? 0).toLocaleString(), pv = (k, d) => P[k] ?? d;
  const model = r => (L.models.find(x => x.role === r) || {}).model || '';
  const src = c => String(c.a.sources || '').split(' ');
  const nBy = f => conds.filter(f).length;
  const ag = o => o && (o.both || o.adjudicated) ? `${n(o.both)} agreed, ${n(o.adjudicated)} adjudicated` : '';
  const tags = (...ks) => h('span', {}, ks.map(tag));
  const code = (...xs) => h('span', {}, xs.map((x, i) => [i ? ', ' : '', h('code', {}, x)]).flat());
  const E = D.cond_evidence || {};
  const O = D.ontology || {}, onto = [O.parents !== false ? `a parent of ${S} (one with more than ` +
    `${pv('MaxOntologyChildren', 30)} children is too broad a class and gives none)` : '',
    O.parents !== false && O.siblings !== false ? 'another child of such a parent' : '',
    O.children !== false ? `a child of ${S}` : '', O.children_of ? 'a child of an entry the profile names' : '']
    .filter(Boolean);
  const parts = {};

  // 4.1 the routes into the graph
  parts.routes = [plainTable(['Route', 'Decided by', 'A condition enters when', 'Its category', 'Evidence type'], [
    ['Clinical evidence', tags(cur ? 'profile' : 'model', 'pubmed', 'rule'), `A PubMed abstract states, for ${S} ` +
      'itself, an association with the condition, quoted verbatim. A study that found no association is kept beside a ' +
      'relation as contrary evidence, never as one.', 'The most frequent category of its evidence rows (4.4).',
      code('clinical')],
    ['Disease ontology', tags('ot', 'rule'), `It is a neighbor of ${S} in the disease ontology: ${onto.join(', or ')}.`,
      `${D.ontology_category}, when it has no clinical evidence.`, code('ontology')],
    ['Shared GWAS signal', tags('ot', 'rule'), `One of its GWAS shares a signal with a ${S} GWAS region (4.3). Related ` +
      `by this alone, it must share at least ${pv('MinSharedRegions', 2)} regions and not be a broad class.`,
      `${SHOWN_CAT(D.shared_signal)}, when neither of the above relates it.`,
      code('colocalization', 'lead-variant overlap')]]),
    here(`${conds.length} related conditions: ${nBy(c => src(c).includes('clinical'))} by clinical evidence, ` +
      `${nBy(c => src(c).includes('ontology'))} by the ontology, ${nBy(c => src(c).includes('colocalization') ||
      src(c).includes('lead-variant'))} by a shared GWAS signal (a condition can have several).`)];

  // 4.2 the clinical relations, step by step
  if (cur) {
    parts.clinical = [h('p', {}, 'The clinical relations come from the curated list of the disease profile, each with ' +
      'its category and a verbatim quote; no language model reads them. Every quote is checked against its PubMed ' +
      'abstract at each build, and the build stops if one is not found.')];
  } else {
    const C = (R && R.consolidation) || {}, V = (R && R.verification) || {};
    const WHY = {not_about_disease: `not about ${S}`, no_relation_stated: 'no relation stated',
                 other_disease_comparison: 'a comparison with another disease only',
                 consequence_or_finding: 'a consequence or finding', other_form: `another form of ${S}'s parent class`};
    const dropped = Object.entries(V.rows_dropped || {}).map(([k, v]) => `${WHY[k] || k} ${v}`).join(', ');
    parts.clinical = [
      R ? flowBoxes([
        {n: R.records, unit: 'records', label: 'found in PubMed', by: ['pubmed', 'rule']},
        {n: R.relations, unit: 'relations', label: 'read by two models and checked', by: ['model', 'rule']},
        {n: R.conditions_read, unit: 'conditions', label: 'mapped to Open Targets' + (R.split ? `; ${R.split} ` +
          `joined name${R.split > 1 ? 's' : ''} split` : ''), by: ['ot', 'rule']},
        {n: C.conditions, unit: 'conditions', label: 'after the consolidation', by: ['model']},
        {n: V.conditions, unit: 'conditions', label: 'after the verification', by: ['model']},
        {n: nBy(c => src(c).includes('clinical')), unit: 'conditions', label: 'related by clinical evidence',
          by: ['rule']}]) : null,
      procTable([
        ['1', 'Search', ['pubmed', 'rule'], `The terms of ${S}, and a guideline, consensus statement, meta-analysis, ` +
          'systematic review, cohort or case-control study, and a relation word (risk, associat*, comorbid*, ' +
          'complicat*, predispos*, preval*, inciden*), to the snapshot date.', R ? `${n(R.records)} records` : ''],
        ['2', 'Reading', ['model'], `Two models read every abstract apart (${model('Reader A')}, ` +
          `${model('Reader B')}), with only the text given. For each relation: the condition as written and its usual ` +
          'name, a category (by name only at this step), the direction, design, population, a one-line statement and ' +
          `a verbatim quote. Not a relation: a consequence, complication or outcome of ${S} or of its treatment; a ` +
          `sign, symptom or finding among its patients; a treatment; a gene; a demographic attribute; ${S} itself or ` +
          `a subtype; a comparison between groups of ${S} patients; evidence from animals or cells.`, ''],
        ['3', 'Source-text check', ['rule'], 'The quote must occur in the abstract word for word and name the ' +
          'condition: its full name, an alias in parentheses, every word of its name, or an abbreviation for it (one ' +
          'the record defines, or one whose letters the words of the name spell, written out in the same record). ' +
          'Category and direction must be among the allowed ones. Any other item is dropped and logged.', ''],
        ['4', 'Agreement', ['rule', 'model'], 'A rule compares the two readings of each record by condition, ' +
          `direction and category; where they differ, a third model (${model('Adjudicator')}) reads the abstract with ` +
          'both readings and decides, and its answer passes the same check.',
          R ? `${n(R.relations)} relations (${ag(R.agreement)})` + (R.not_read ? `; ${n(R.not_read)} records not read` :
            '') : ''],
        ['5', 'Mapping', ['ot', 'rule'], 'A name becomes an Open Targets entry when it equals the entry\'s name or an ' +
          'exact synonym (a level written with its direction, such as “elevated factor VIII”, also as the measured ' +
          `trait). A name that maps to an ancestor of ${S} keeps no entry; ${S} itself and its subtypes leave. Names ` +
          'without an entry that differ only in words of use or exposure are one condition; a name joining other ' +
          'listed conditions (“A / B”) is split into them.',
          R && R.conditions_read ? `${n(R.conditions_read)} conditions` + (R.split ? `; ${n(R.split)} split` : '') :
            ''],
        ['6', 'Consolidation', ['model'], 'Both readers see the whole list and decide for each condition: keep; ' +
          'merge (the same condition under another name, or the marker that defines it); feature (a characteristic of ' +
          'another listed condition: severity, type, provocation, a finding at its diagnosis, its treatment; its ' +
          `evidence moves there); or drop (not a condition, or a consequence of ${S}). The adjudicator decides where ` +
          'they differ. Categories do not change here.',
          C.items ? `${n(C.items)} → ${n(C.conditions)} (` + tallyText({merged: C.merge, attached: C.feature,
            dropped: C.drop}) + `); ${ag(C.agreement)}` : ''],
        ['7', 'Verification', ['model'], 'Both readers judge every evidence row against its quote, with the category ' +
          `definitions (4.4) and the ontology neighbors of ${S}; the adjudicator decides where they differ. A row is ` +
          `kept with its category, or dropped as: not about ${S} (another disease, or a wider group that also holds ` +
          'patients without it); no relation stated (background, an aim, a resemblance); a comparison with another ' +
          `disease only; a consequence or finding among patients; a condition defining another form of ${S}'s parent ` +
          'class. Each condition is then kept, dropped, merged or attached; one left without a row leaves.',
          V.conditions_in ? `${n(V.conditions_in)} → ${n(V.conditions)}; rows dropped: ${dropped || 'none'}; ` +
            `${n(V.rows_recategorized)} rows recategorized; ${ag(V.agreement)}` : ''],
        ['8', 'Category and build check', ['rule'], 'A condition takes the most frequent category of its own rows that ' +
          'relate it (rows that found no association, and rows of an attached characteristic, do not vote; of tied ' +
          'categories, the one met first). At each build every quote is checked again against its abstract, and the ' +
          'build stops if one is not found. A condition whose clinical rows all found no association is not related.',
          (R ? `${n(R.conditions)} conditions` : '') + (L.not_related.length ? `; not related: ` +
            L.not_related.join('; ') : '')]]),
      h('p', {class: 'note'}, `The models read, adjudicate, consolidate and verify (steps 2, 4, 6, 7): readers ` +
        `${model('Reader A')} and ${model('Reader B')}, adjudicator ${model('Adjudicator')}. They never decide the ` +
        'search, the Open Targets entries, the ontology neighbors, the shared signals, the category vote or which ' +
        'genetics are read. Every answer is held to fixed lists (categories, directions, actions, reasons) and to its ' +
        'source text; an answer outside them is not used, and is logged. Every decision, with both readings and the ' +
        'adjudicator\'s reason, is in the snapshot\'s 01.relations folder; the disputed ones are listed in section 6.')];
  }

  // 4.4 the categories
  const tally = count(conds, c => c.a.category);
  parts.categories = [plainTable(['Category', 'Definition', S], CATS.map(c => [h('b', {}, SHOWN_CAT(c)),
    D.category_desc[c] || '', h('span', {class: 'num'}, n(tally[c] || 0))])),
    h('p', {class: 'note'}, 'A clinical category comes first; the ontology category is given only to a neighbor ' +
      'without clinical evidence, and the shared-signal category only to a condition related by nothing else. The ' +
      'readers see the clinical definitions when they verify each row (4.2, step 7).')];

  // 4.6 the evidence matrix
  const dz = count(conds.filter(c => (E[c.id] || {}).design != null), c => D.relation_designs[E[c.id].design]);
  parts.matrix = [plainTable(['Column', 'Meaning'], [
    ['Publications', 'The publications whose clinical evidence relates the condition; several statements of one ' +
      'publication count once, and a study that found no association is not counted. Shaded by bins (1, 2–3, 4–9, ' +
      '10 or more) so that a condition with many publications does not hide those with few; the number is in the cell.'],
    ['Best design', 'The best design among those publications, in this order: ' + D.relation_designs.join('; ') +
      '. From the PubMed publication types and MeSH design headings as NLM indexes each record (a case report stays ' +
      'one); the readers\' design words only for a record not yet indexed. A guideline, meta-analysis or systematic ' +
      'review counts as one only when it is about the relation: its title names the condition, or the risk factors, ' +
      `causes, comorbidities or associated conditions of ${S}; otherwise the quoted sentence is background from a ` +
      'synthesis of another question, and counts as other. The search asks for these designs, so most relations rest ' +
      'on cohort or case-control studies; the column sets apart those that a guideline or a meta-analysis is about.'],
    ['Direction', '↑ increased risk, ↓ decreased risk, ~ associated, as the publications state it; ✕ when a study ' +
      'found no association (shown on its page as contrary evidence).'],
    ['EA', 'East Asian evidence of the relation (section 5).'],
    ['Relation', 'Its place in the disease ontology: parent, sibling or child.'],
    [`Shared ${S} GWAS regions`, 'A filled dot: a credible set of the condition colocalizes with the region; an open ' +
      'dot: it holds the region\'s lead variant (lead-variant overlap). With more than 20 regions, a count.'],
    ['Its own genetics', 'Hidden unless shown with its button: its own GWAS studies in Open Targets (hover: how many ' +
      'are East Asian) and its genes with curated gene–disease evidence; “not fetched” and “no Open Targets entry” ' +
      'as in 4.5.']]),
    here('best design: ' + (tallyText(dz, D.relation_designs) || 'none') + '.')];

  // 4.5 its own genetics: read in full, not read, or nothing to read
  const nfc = count(conds, c => (E[c.id] || {}).nf || 'full'), NF = () => h('b', {}, 'Not fetched');
  parts.own = [plainTable(['', 'When', 'Rule in this run', S], [
    [h('b', {}, 'Fetched in full'), 'Every other condition with an Open Targets entry: its GWAS studies and credible ' +
      'sets, the top L2G gene of each, its curated gene–disease evidence. These bring its genes into the graph.', '—',
      n(nfc.full || 0)],
    [NF(), `Related only by a shared GWAS signal: the relation is the shared signal itself, and its other genes say ` +
      `nothing about ${S}.`, `shares at least ${pv('MinSharedRegions', 2)} regions (MinSharedRegions)`,
      n(nfc.shared || 0)],
    [NF(), 'A broad class: its own GWAS and curated genes are those of its many subtypes, each listed on its own ' +
      'when related.', `more than ${pv('MaxDescendants', 100)} subtypes (MaxDescendants)`, n(nfc.broad || 0)],
    ...(D.cooccurrence.length ? [[NF(), `Related only as ${D.cooccurrence.join(' or ').toLowerCase()}, not an ` +
      'ontology neighbor, on too few publications: co-occurrence says nothing of which came first, and its genes say ' +
      `little about the causes of ${S} until the relation is established.`,
      `fewer than ${pv('ComorbidityMinPublications', 2)} publications (ComorbidityMinPublications)`,
      n(nfc.cooccurrence || 0)]] : []),
    [NF(), 'A general trait (a blood count, body size): its genes are so many that they would flood the gene map.',
      `more than ${pv('MaxGwasStudies', 100)} GWAS studies (MaxGwasStudies)`, n(nfc.gwas || 0)],
    [h('b', {}, 'No Open Targets entry'), 'No entry in Open Targets: nothing to read; the relation stands on its ' +
      'clinical evidence.', '—', n(nfc.none || 0)]]),
    h('p', {class: 'note'}, `A condition whose own genetics are not fetched keeps its relation, its clinical evidence ` +
      `and its signals shared with ${S}, which are looked up from the regions of ${S}. Its counts are left blank, not ` +
      'zero: a zero would read as “looked up, none found”. When several reasons apply, the first in the order of ' +
      'this table is given, on hover in the matrix and on its page.')];

  // 4.7 what is left out
  const lo = L.shared_left_out || {}, few = lo.too_few || {}, broad = lo.broad || {};
  const names = (title, xs) => h('details', {class: 'cmap-all'}, h('summary', {}, title),
    h('div', {class: 'note'}, xs.join('; ')));
  parts.leftout = [h('p', {}, 'A condition whose clinical evidence all found no association is not related (4.2, ' +
    'step 8). A trait related only by a shared signal is left out when it shares too few regions, or when it is a ' +
    'broad class (each of its subtypes that shares a signal is listed on its own). Each is counted and named here.'),
    here(L.not_related.length ? `not related: ${L.not_related.join('; ')}.` : 'no condition left as not related.'),
    few.n ? names(`${S}: ${n(few.n)} trait${few.n > 1 ? 's' : ''} sharing fewer than ${few.min_regions} regions`,
      few.names || []) : null,
    broad.n ? names(`${S}: ${n(broad.n)} broad class${broad.n > 1 ? 'es' : ''}`, broad.names || []) : null];
  return parts;
}

function regionsPage() {
  const regions = byType('region');
  const out = [h('h2', {}, `${S} GWAS regions`), h('p', {class: 'lead'}, D.text.regions_lead)];
  const notRegions = () => D.not_regions.length ? [h('h3', {}, 'Genome-wide significant signals that are not regions'),
    table([{label: 'Lead variant', get: r => r.rsid || r.variant_id, cell: r => r.rsid ? rsid(r.rsid) : r.variant_id},
           {label: 'Study', get: r => r.study_id, cell: r => ext(`https://platform.opentargets.org/study/${r.study_id}`,
             r.study_id)},
           {label: 'P', num: true, get: r => Number(r.p_mantissa) * Math.pow(10, Number(r.p_exponent)),
            cell: r => `${Number(r.p_mantissa).toFixed(1)}e${r.p_exponent}`},
           {label: 'MAF', num: true, get: r => r.maf}, {label: 'Expected copies', num: true, get: r => r.expected_mac},
           {label: 'Why not a region', get: r => r.reason}], D.not_regions, {filter: false})] : [];
  if (!regions.length) return out.concat(notRegions());
  // a statistic of the published table is a column only when a region has it
  const has = k => regions.some(r => r.a[k] !== undefined);
  const cols = [{label: 'Region', get: r => r.label, cell: r => nodeLink(r.id)}, {label: 'Code', get: r => r.a.code}];
  if (has('tier')) cols.push({label: 'Tier', num: true, get: r => r.a.tier});
  cols.push({label: 'Lead variant', get: r => r.a.rsid, cell: r => rsid(r.a.rsid)},
            {label: 'P', num: true, get: r => r.a.p});
  if (has('v')) cols.push({label: D.region_labels.v, num: true, get: r => r.a.v});
  if (has('or')) cols.push({label: 'OR', num: true, get: r => r.a.or});
  cols.push({label: 'Source', get: r => r.a.source || 'published'},
            {label: 'In Open Targets', get: r => r.a.open_targets},
            {label: 'Conditions sharing it', num: true, get: r => (IN.get(r.id) || []).filter(e =>
              e.type === 'shares_region').length});
  return out.concat([table(cols, regions, {filter: false}), h('h3', {}, 'Conditions sharing them'), regionMatrix()],
                    notRegions());
}

function studiesPage() {
  const all = byType('publication').filter(p => p.a.design_class);
  let place = 'all', design = 'all';
  const box = h('div');
  const designs = [...new Set(all.map(p => p.a.design_class))].sort();
  const draw = () => {
    // the study's own flags (done there, or a population it names), as in the workbook and the knowledge base
    const rows = all.filter(p => (place === 'all' || (place === 'Japanese' ? p.a.japanese === 'True' :
      place === 'East Asian' ? p.a.east_asian === 'True' : p.a.east_asian !== 'True')) &&
      (design === 'all' || p.a.design_class === design));
    const genesOf = p => [...new Set((OUT.get(p.id) || []).filter(e => e.type === 'reports' && N.get(e.t)?.type ===
      'gene').map(e => e.t))];
    box.replaceChildren(table([
      {label: 'Year', num: true, get: p => p.a.year}, {label: 'Study', get: p => p.label, cell: p => nodeLink(p.id)},
      {label: 'Place', get: p => p.a.place}, {label: 'Design class', get: p => p.a.design_class},
      {label: 'Cases', num: true, get: p => p.a.n_cases}, {label: 'Controls', num: true, get: p => p.a.n_controls},
      {label: 'Genes reported', get: p => genesOf(p).map(g => N.get(g).label).join(' '),
       cell: p => h('span', {}, genesOf(p).map(g => h('a', {class: 'chip', href: href(g)}, N.get(g).label)))},
      {label: 'Readers', get: p => p.a.agreement, cell: p => agree(p.a.agreement)},
      {label: 'Title', get: p => p.a.title}, {label: 'Source', get: p => p.id, cell: p => pmid(p.id.split(':')[1])}],
      rows, {sort: 0, tools: h('span', {}, ['all', 'Japanese', 'East Asian', 'Not East Asian'].map(x =>
        h('button', {class: 'btn' + (place === x ? ' on' : ''), onclick: () => { place = x; draw(); }},
          x === 'all' ? 'All places' : x)), ' ',
        h('select', {onchange: e => { design = e.target.value; draw(); }},
          h('option', {value: 'all'}, 'All designs'), designs.map(d => h('option', {value: d,
            selected: d === design ? 'selected' : null}, d))))}));
  };
  draw();
  return [h('h2', {}, `Genetic studies of ${S}`), h('p', {class: 'lead'}, `The genetic studies of ${S} the curation ` +
    "keeps, with their place (country of the first author's institution), design class and the agreement of the " +
    'model readers (adjudicated: the readers disagreed and an adjudicator decided).'), box];
}

function curationPage() {
  const C = D.curation;
  const pubCols = [{label: 'Publication', get: r => `${r.first_author || ''} ${r.year || ''}`},
                   {label: 'Title', get: r => r.title}, {label: 'Reason', get: r => r.reason},
                   {label: 'Used in', get: r => r.roles}, {label: 'Source', get: r => r.pmid, cell: r => pmid(r.pmid)}];
  const findings = D.edges.filter(([, , t, a]) => t === 'reports' && a.about && a.about !== ABOUT)
    .map(([s, t, , a]) => ({s, t, a}));
  const studies = byType('publication').filter(p => p.a.design_class);
  const reports = D.edges.filter(([, , t]) => t === 'reports');
  const log = C.log.map(([step, item, rule, before, after]) => ({step, item, rule, before, after}));
  const tmCols = [{label: 'Gene', get: r => r.item, cell: r => h('span', {class: 'gene'}, r.item)},
                  {label: 'Decision', get: r => r.rule}, {label: 'Links', get: r => r.before},
                  {label: 'What the word means',
                   get: r => r.after.replace(/^(mapped to the wrong gene|words removed): /, '')}];
  const vnCols = [{label: 'Finding', get: r => r.item}, {label: 'As written', get: r => r.before},
                  {label: 'Name', get: r => r.after.replace(/^named: /, '')}];
  return [h('h2', {}, 'Curation'), h('p', {class: 'lead'}, `Every publication the graph uses ` +
    `(${C.counts.publications.toLocaleString()}) was looked up in PubMed on ${D.meta.checked_at}. ` +
    D.text.curation_lead),
    h('div', {class: 'note'}, `Agreement of the model readers: genetic studies ${tally(studies.map(p =>
      p.a.agreement))}; findings ${tally(reports.map(([, , , a]) => a.agreement))}.`),
    h('h3', {}, `Findings not about ${S} itself (${findings.length})`),
    table([{label: 'Study', get: r => N.get(r.s)?.label, cell: r => nodeLink(r.s)},
           {label: 'Gene', get: r => N.get(r.t)?.label, cell: r => N.get(r.t)?.type === 'gene' ? nodeLink(r.t) : ''},
           {label: 'Result', get: r => r.a.result}, {label: 'About', get: r => r.a.about},
           {label: 'Source sentence', get: r => r.a.quote, cell: r => h('span', {class: 'quote'}, r.a.quote)},
           {label: 'Readers', get: r => r.a.agreement, cell: r => agree(r.a.agreement)}],
          findings, {filter: false}),
    h('h3', {}, 'Text-mined links mapped to the wrong gene'), table(tmCols, log.filter(r => r.step === 'text mining')),
    h('h3', {}, 'Variants given one name'), table(vnCols, log.filter(r => r.rule === 'variant name')),
    h('h3', {}, `Publications flagged and kept (${C.flagged.length})`), table(pubCols, C.flagged)];
}

// ── how language models are used: the general rules, then this disease in this snapshot ──────────────────────
// one reading, from record to graph: the models (purple) and the rules (gray) it passes
function flowDiagram(models) {
  const name = r => (models.find(m => m.role === r) || {}).model || '';
  const svg = svgEl('svg', {viewBox: '0 0 1170 286', class: 'flow', role: 'img',
                            'aria-label': 'One reading, from record to graph'});
  svg.append(svgEl('defs', {}, svgEl('marker', {id: 'arr', viewBox: '0 0 10 10', refX: 9, refY: 5, markerWidth: 7,
    markerHeight: 7, orient: 'auto-start-reverse'}, svgEl('path', {d: 'M0,0L10,5L0,10z', fill: '#57606a'}))));
  const box = (x, y, kind, title, sub, w = 150, hgt = 58) => {
    const fill = {model: '#efe8f7', rule: '#f6f7f9', data: '#ffffff', out: '#1f3b5f'}[kind];
    const stroke = {model: '#6a4c93', rule: '#8a9199', data: '#8a9199', out: '#1f3b5f'}[kind];
    const g = svgEl('g', {}, svgEl('rect', {x, y, width: w, height: hgt, rx: 7, fill, stroke, 'stroke-width': 1.4,
                                            'stroke-dasharray': kind === 'data' ? '5 3' : null}));
    g.append(svgEl('text', {x: x + w / 2, y: y + (sub ? 23 : hgt / 2 + 5), 'text-anchor': 'middle',
                            style: `font-size:14px;font-weight:700;fill:${kind === 'out' ? '#fff' : '#1f2328'}`}, title));
    if (sub) g.append(svgEl('text', {x: x + w / 2, y: y + 42, 'text-anchor': 'middle',
                                     style: `font-size:11.5px;fill:${kind === 'out' ? '#fff' : '#57606a'}`}, sub));
    svg.append(g);
  };
  const arrow = (d, label, lx, ly) => {
    svg.append(svgEl('path', {d, fill: 'none', stroke: '#57606a', 'stroke-width': 1.4, 'marker-end': 'url(#arr)'}));
    if (label) svg.append(svgEl('text', {x: lx, y: ly, 'text-anchor': 'middle',
                                         style: 'font-size:11.5px;font-style:italic;fill:#57606a'}, label));
  };
  // the order of the code: each reading is checked against the source text, then the readings are compared;
  // the adjudicator's reading is checked the same way
  box(8, 112, 'data', 'Source text', 'abstract or methods');
  box(175, 112, 'rule', 'Pre-screen', 'record types, terms');
  box(342, 30, 'model', 'Reader A', name('Reader A'));
  box(342, 194, 'model', 'Reader B', name('Reader B'));
  box(509, 30, 'rule', 'Source-text check', 'not found: dropped');
  box(509, 194, 'rule', 'Source-text check', 'not found: dropped');
  box(676, 112, 'rule', 'Compare', 'the decisive items');
  box(860, 200, 'model', 'Adjudicator', name('Adjudicator'));
  box(1046, 112, 'out', 'Graph', 'with its quote', 110);
  arrow('M158,141 L173,141');
  arrow('M325,130 L340,76'); arrow('M325,152 L340,208');
  arrow('M492,59 L507,59'); arrow('M492,223 L507,223');
  arrow('M659,74 L674,128'); arrow('M659,208 L674,156');
  arrow('M826,141 L1044,141', 'agree', 935, 133);
  arrow('M751,170 L751,229 L858,229', 'differ', 800, 221);
  arrow('M1010,229 L1101,229 L1101,172', '', 0, 0);
  svg.append(svgEl('text', {x: 935, y: 274, 'text-anchor': 'middle',
                            style: 'font-size:11.5px;font-style:italic;fill:#57606a'},
                   'its reading checked against the text in the same way'));
  const key = h('div', {class: 'cmap-key'},
    h('span', {}, h('span', {class: 'swatch', style: 'background:#efe8f7;border-color:#6a4c93'}), 'a language model'),
    h('span', {}, h('span', {class: 'swatch', style: 'background:#f6f7f9;border-color:#8a9199'}), 'a rule (no model)'));
  return h('div', {class: 'flow-wrap'}, svg, key);
}

// the kinds of source of a step or an association
const SRC = {ot: ['Open Targets', '#3c6e9f'], pubmed: ['PubMed', '#2e7d6b'], epmc: ['Europe PMC', '#4f8fbf'],
             profile: ['disease profile', '#b08a3c'], model: ['language models', '#6a4c93'], rule: ['rules', '#6e7781']};
const tag = k => h('span', {class: 'badge', style: 'background:' + SRC[k][1]}, SRC[k][0]);
// a flow of counts: each box a number, its unit, what it is and who acts at that step (a rule, a database, the models)
const flowBoxes = items => h('div', {class: 'fb'}, items.map((x, i) => [i ? h('div', {class: 'fb-arrow'}, '→') : null,
  h('div', {class: 'fb-box'}, h('div', {class: 'fb-n'}, (x.n ?? 0).toLocaleString(), h('span', {class: 'fb-u'},
    ' ' + x.unit)), h('div', {class: 'fb-l'}, x.label), h('div', {class: 'fb-by'}, (x.by || []).map(tag)))]).flat());
const here = (...kids) => h('div', {class: 'here'}, h('b', {}, S + ': '), ...kids);
const count = (xs, f) => xs.reduce((m, x) => { const k = f(x); m[k] = (m[k] || 0) + 1; return m; }, {});
const tallyText = (m, order) => (order || Object.keys(m)).filter(k => m[k]).map(k => `${k} ${m[k]}`).join('; ');

// the build page: a methods document in numbered sections that follow the build, each in one form (its answer in
// short; what happens, step by step, by whom; its rules with their values and basis, from the choices registry; this
// snapshot's counts), with a contents column; every rule of the registry is in one section, and any other in the last
function buildPage(sel) {
  const P = D.params, L = D.llm, m = D.meta, cc = D.curation, R = L.relations, X = L.extract || {};
  const num0 = v => (v ?? 0).toLocaleString(), pv = (k, d) => P[k] ?? d;
  const nEdges = t => D.edges.filter(e => e[2] === t).length;
  const regions = byType('region'), conds = related(), genes = byType('gene').filter(g => lvl(g) !== null);
  const studies = byType('publication').filter(p => p.a.design_class);
  const nPub = regions.filter(published).length;
  const model = r => (L.models.find(x => x.role === r) || {}).model || '';
  const parts = relationsParts();

  // the rules of a section, from the choices registry by key; the keys used are kept, so that the last section lists
  // every rule no section names
  const RULES = new Map(D.choices.flatMap(g => g.items.map(it => [it.key, it])));
  const used = new Set();
  const refs = r => (r.refs || []).length ? h('div', {class: 'rule-refs'}, r.refs.map(([t, u], k) =>
    h('span', {}, k ? '; ' : '', u ? ext(u, t) : t))) : null;
  const ruleRows = items => h('table', {class: 'plain rules'},
    h('thead', {}, h('tr', {}, ['Rule', 'In this run', 'Basis'].map(x => h('th', {class: 'nosort'}, x)))),
    h('tbody', {}, items.map(r => h('tr', {}, h('td', {class: 'rule-k'}, r.item), h('td', {class: 'rule-v'}, r.value),
      h('td', {}, r.basis, refs(r))))));
  const rules = (...keys) => {
    const items = keys.map(k => RULES.get(k)).filter(Boolean);
    items.forEach(r => used.add(r.key));
    return items.length ? [h('div', {class: 'doc-k'}, 'Rules'), ruleRows(items)] : [];
  };
  // a section and a subsection: numbered, with an anchor (#/build/<key>) and its answer in short
  const TOC = [];
  const sec = (no, key, title, short, ...kids) => {
    TOC.push([no, key, title, 1]);
    return h('section', {class: 'doc-sec', id: 'b-' + key}, h('h3', {}, h('span', {class: 'doc-no'}, no), title),
      short ? h('p', {class: 'insh'}, short) : null, ...kids.flat());
  };
  const sub = (no, key, title, short, ...kids) => {
    TOC.push([no, key, title, 2]);
    return h('div', {class: 'doc-sub', id: 'b-' + key}, h('h4', {}, h('span', {class: 'doc-no'}, no), title),
      short ? h('p', {class: 'insh2'}, short) : null, ...kids.flat());
  };
  const k = text => h('div', {class: 'doc-k'}, text);

  const body = [];
  // 1 design and scope
  body.push(sec('1', 'design', 'Design and scope', 'Databases give the facts, two language models read the ' +
    'literature, and rules decide what the graph counts. The rules are the same for every disease; what is specific ' +
    `to ${S} is its profile and the counts given in each section.`,
    h('ul', {class: 'doc-list'},
      h('li', {}, h('b', {}, 'Every statement has a source: '), 'a database record with its identifier, or a verbatim ' +
        'quote of a PubMed abstract or open article, linked from the graph.'),
      h('li', {}, h('b', {}, 'Models read, rules decide: '), 'a language model\'s answer enters the graph only through ' +
        'checks against its source text; regions, levels, categories, counts and filters are rules (section 6).'),
      h('li', {}, h('b', {}, 'One set of rules for every disease: '), 'a disease is its Open Targets ID; optional ' +
        'curated files of its profile replace a default and are named where they do.'),
      h('li', {}, h('b', {}, 'Complete, with stated exclusions: '), 'what a rule leaves out is counted and named ' +
        '(4.7, the Curation page).'),
      h('li', {}, h('b', {}, 'Dated and rebuildable: '), 'a snapshot keeps what it downloaded and its model answers ' +
        '(section 9).')),
    rules('fetch_date', 'disease', 'evidence_scope', 'sources', 'external_only', 'excluded')));

  // 2 data sources
  const srcRow = (kk, gives, how) => [tag(kk), gives, how];
  body.push(sec('2', 'sources', 'Data sources', 'Four kinds of source, each read its own way: Open Targets by rules; ' +
    'PubMed and Europe PMC abstracts and methods by the language models, their records by rules; the disease profile ' +
    'by rules.', plainTable(['Source', 'What it gives', 'How it is read'], [
      srcRow('ot', 'The disease and its ontology (names, synonyms, parents, children); GWAS studies and credible sets ' +
        'with locus-to-gene (L2G) scores and colocalizations; curated gene–disease evidence (ClinGen, Gene2Phenotype, ' +
        'PanelApp, Orphanet, ClinVar, UniProt, gene burden); Reactome pathways; text-mined targets.',
        'By rules only; values are kept as given, with their identifiers.'),
      srcRow('pubmed', 'The genetic studies of the disease; the clinical relations; the record of every publication ' +
        'the graph uses (type, retraction, affiliation, MeSH).', 'Abstracts by the language models; types, ' +
        'retractions, affiliations and MeSH by rules.'),
      srcRow('epmc', 'The methods text of open-access articles; the sentences behind text-mined gene mentions.',
        'Methods by the language models; mentions by rules.'),
      srcRow('profile', 'The Open Targets ID (required); optional curated files: clinical relations, published GWAS ' +
        'regions, gene functions, curation entries.', 'By rules.')]),
    here(`Open Targets ${m.release}; PubMed to ${L.fetch_date || m.snapshot}; every publication checked in PubMed ` +
      `on ${m.checked_at} (${num0(cc.counts.publications)} publications).`)));

  // 3 the genetic evidence of the disease
  const nr = D.not_regions.length;
  const ps = D.prescreen || {};
  const terms = (ps.genetic_terms || []).map(t => t.replace(/\\b|\\d\{\d*,?\d*\}|[?\\]/g, '')).filter(Boolean);
  const lvCount = count(genes, g => String(lvl(g)));
  body.push(sec('3', 'genetics', `The genetic evidence of ${S}`, `Two lines of evidence rank the genes of ${S}: ` +
    'its GWAS regions (3.1) and its genetic studies in the literature (3.2), curated by rules (3.3). Each gene takes ' +
    'the strongest level it reaches on one seven-level ladder (3.4).',
    sub('3.1', 'regions', `GWAS regions of ${S}`, `A region is a genome-wide significant locus of ${S} with a ` +
      'credible lead variant, published or from Open Targets. Its genes are the strongest genetic evidence of the ' +
      'graph, and each related condition is tested against it for a shared signal (4.3).',
      h('p', {}, D.text.regions_lead),
      rules('published_regions', 'genome_wide', 'region_mac', 'region_maf', 'region_rare_p', 'region_window',
        'region_names'),
      here(`${regions.length} region${regions.length === 1 ? '' : 's'}` + (regions.length ? ` (${nPub} published, ` +
        `${regions.length - nPub} from Open Targets)` : '') + (nr ? `; ${nr} genome-wide significant ` +
        `signal${nr === 1 ? '' : 's'} not a region, each with its reason` : '') + '. Every one is on the ',
        h('a', {href: '#/regions'}, 'GWAS regions'), ' page.')),
    sub('3.2', 'studies', `Genetic studies of ${S} in the literature`, `A publication is a genetic study of ${S} when ` +
      'two models, reading it apart, include it, and every gene, variant and number they take from it is found in ' +
      'its text.',
      flowBoxes([
        {n: X.records, unit: 'records', label: 'found in PubMed', by: ['pubmed', 'rule']},
        {n: X.read_by_model, unit: 'records', label: 'read by the models after the pre-screen', by: ['rule']},
        {n: X.included, unit: 'studies', label: 'included by the readers', by: ['model']},
        {n: X.findings, unit: 'findings', label: 'extracted and found in the text', by: ['model', 'rule']},
        {n: cc.counts.studies, unit: 'studies', label: 'kept after the curation', by: ['rule']}]),
      procTable([
        ['1', 'Search', ['pubmed', 'rule'], `The name, short name and exact synonyms of ${S} in title or abstract, ` +
          'together with a genetic-study term, to the snapshot date.', `${num0(X.records)} records`],
        ['2', 'Pre-screen', ['rule'], 'Before any model reads a record: a record of a type that reports no original ' +
          `data (${(ps.excluded_types || []).join(', ').toLowerCase()}) is left out, and so is one whose title and ` +
          `abstract name no genetic-study term (${terms.length} terms, such as ${terms.slice(0, 8).join(', ')}).`,
          `${num0(X.prescreened_out)} left out; ${num0(X.read_by_model)} read`],
        ['3', 'Reading', ['model'], `Two models read each abstract apart and decide whether it is a genetic study of ` +
          `${S}; they extract its findings (gene, variant, result, effect, P), its sample sizes and whether it was ` +
          'done in Japan or East Asia, and quote its methods from the open full text where there is one. The ' +
          'adjudicator decides where they differ (section 6).',
          `${num0(X.included)} included; ${num0(X.findings)} findings`],
        ['4', 'Source-text check', ['rule'], 'Every gene, variant, number and quote must occur in its text, or the ' +
          'item is dropped and logged; a result label that the wording of its own sentence contradicts is flagged.',
          `${num0(X.dropped_items)} items dropped; ${num0(X.result_flags)} flagged`]]),
      rules('query', 'prescreen_types', 'prescreen_terms', 'inclusion', 'comparison', 'grounding', 'result_check',
        'methods_fields', 'fallback'),
      h('p', {}, 'Every study with its findings: ', h('a', {href: '#/studies'}, 'Genetic studies'), '.')),
    sub('3.3', 'curation', 'Curation of the literature', 'Rules check every publication in PubMed, give each study ' +
      'its design and place, name each variant once, mark what each finding is about, and correct text mining.',
      procTable([
        ['1', 'Publications', ['pubmed', 'rule'], 'Retracted articles and conference material leave the graph; ' +
          'preprints and errata are flagged and kept.', `${cc.n_excluded} left out; ${cc.flagged.length} flagged`],
        ['2', 'Design and place', ['rule'], 'The design class is the first of these whose rule the study meets, from ' +
          `the study type and design words the readers quote: ${D.design_classes.join('; ')}. The place is the ` +
          'country of the first author\'s institution (PubMed, else Europe PMC).',
          tallyText(count(studies, p => p.a.design_class))],
        ['3', 'Variants and findings', ['rule'], 'A variant gets one name; a finding is marked by what it is about ' +
          `(${S} itself, another trait, a functional result, or background), and only findings about ${S} itself ` +
          'count for levels 3–6.', `${num0(cc.counts.studies)} studies kept`],
        ['4', 'Text mining', ['epmc', 'rule'], 'A text-mined link to a gene its sentence does not name, or a word ' +
          `naming ${S} or a neighbor, is removed.`, `${num0(cc.counts.textmining_dropped)} links removed`]]),
      rules('exclusion', 'flags', 'design_class', 'place', 'variants', 'about', 'textmining_review', 'curated_symbols'),
      h('p', {}, 'Every curation decision: ', h('a', {href: '#/curation'}, 'Curation'), '.')),
    sub('3.4', 'ladder', 'Evidence ladder and gene functions', 'Each gene takes the strongest of seven levels it ' +
      `reaches, from a GWAS region or curated evidence (1) to being named without a result about ${S} (7); every ` +
      'disease has the same seven levels.',
      h('table', {class: 'ladder-def'}, h('tbody', {}, Object.entries(D.ladder).map(([kk, [code, label]]) =>
        h('tr', {class: 'lv' + kk}, h('td', {class: 'lvk'}, h('b', {}, kk), code ? ` ${code}` : ''),
          h('td', {class: 'lvn'}, label), h('td', {}, D.ladder_def[kk],
            D.level_status[kk] ? h('div', {class: 'note'}, `${S}: ${D.level_status[kk].full}`) : null),
          h('td', {class: 'num', style: 'white-space:nowrap'}, `${lvCount[kk] || 0} gene${lvCount[kk] === 1 ? '' :
            's'}`))))),
      h('p', {class: 'note'}, `${D.text.ladder_note} ${D.no_level}.`),
      rules('g_codes', 'curated_core', 'ladder', 'themes'),
      here(`${genes.length} genes on the ladder; ${byType('theme').length} functions ` +
        `(${L.themes || 'curated'}).`))));

  // 4 related conditions
  body.push(sec('4', 'relations', 'Related conditions', `A condition is related to ${S} by a quoted clinical ` +
    'relation, by the disease ontology or by GWAS signals shared with it. The models read and check the clinical ' +
    'relations; rules and databases decide the rest, including the category and which genetics are read.',
    sub('4.1', 'routes', 'Routes into the graph', 'Three routes, each with its own evidence; a condition can come by ' +
      'several.', parts.routes, rules('relations_source', 'clinical_rule', 'ontology', 'ontology_children')),
    sub('4.2', 'clinical', 'Clinical relations from the literature', 'Two models read the abstracts apart; rules ' +
      'check every quote; the models then consolidate the list and verify every row, and the adjudicator decides ' +
      'where the readers differ.', parts.clinical,
      rules('condition_mapping', 'joined_names', 'consolidation', 'verification', 'contrary', 'relation_design')),
    sub('4.3', 'shared', 'Shared GWAS signals', `A related condition shares a region of ${S} when one of its GWAS ` +
      'credible sets colocalizes with the region\'s, or holds its lead variant where the region has no credible set ' +
      `of ${S}. These are computed by Open Targets and selected here by rules.`,
      rules('clpp', 'h4', 'overlap', 'loci_grading', 'min_shared_regions', 'joint_traits', 'study_once',
        'duplicate_traits'),
      here(`${conds.filter(c => (D.cond_evidence[c.id] || {}).regions && Object.keys(D.cond_evidence[c.id].regions)
        .length).length} related conditions share at least one region.`)),
    sub('4.4', 'categories', 'Categories', 'A condition takes the most frequent category of its clinical evidence; ' +
      'the ontology and the shared signal give their own category only when nothing clinical relates it.',
      parts.categories, rules('categories')),
    sub('4.5', 'own', 'Genetics of the related conditions', 'A related condition\'s own GWAS and curated genes are ' +
      'read in full, unless they would mostly add genes unrelated to ' + S + '; then only its shared signals are kept.',
      parts.own, rules('genetics_mode', 'broad_class', 'comorbidity_genetics', 'max_gwas', 'l2g', 'curated',
        'gene_link')),
    sub('4.6', 'matrix', 'The evidence matrix', 'The Related conditions page shows each condition in one row, its ' +
      'evidence in columns.', parts.matrix),
    sub('4.7', 'leftout', 'What is left out', 'Everything a rule leaves out is counted and named.', parts.leftout)));

  // 5 East Asian evidence
  const eaStudies = studies.filter(p => p.a.east_asian === 'True'), jpStudies = studies.filter(p => p.a.japanese ===
    'True');
  body.push(sec('5', 'eastasia', 'East Asian evidence', 'A study, a relation and a GWAS are flagged East Asian by ' +
    'rules, each on its own evidence.', rules('eas_countries', 'eas_relation', 'eas'),
    here(`${eaStudies.length} genetic studies of ${S} done in East Asia (${jpStudies.length} in Japan); ` +
      `${conds.filter(c => (D.cond_evidence[c.id] || {}).ea).length} related conditions with East Asian evidence of ` +
      'the relation.')));

  // 6 the language models
  body.push(sec('6', 'models', 'Language models', 'Two models read every record apart, a third ' +
    'decides where they differ, and rules check every item against its source text. No model decides a count, a ' +
    'level, a category vote or a rule.', modelSections(), rules('model', 'readers', 'batch', 'isolation')));

  // 7 summaries
  const sm = L.summaries;
  body.push(sec('7', 'summaries', 'Summaries', 'One model writes each summary from the material it is given, the ' +
    'other checks every claim against that material, and the adjudicator rules on each challenged claim; a summary ' +
    'with an unsupported claim is withheld, with its reason.',
    rules('reactome', 'theme_size', 'theme_min', 'sentences', 'summary_checks'),
    here(`${sm.n} summaries: ${sm.shown} shown, ${sm.withheld} withheld` + (sm.not_generated ? `, ` +
      `${sm.not_generated} not generated in this build` : '') + '.')));

  // 8 the graph
  const EDGE = {
    related_to: [['profile', 'model', 'ot'], 'Clinical (4.2), ontology or shared signal (4.3); one edge per piece of ' +
      'evidence.'],
    supports: [['pubmed'], 'The publication a clinical relation quotes.'],
    has_genetic_evidence: [['ot', 'model', 'rule'], `G: the top L2G gene (≥ ${pv('L2gMin', 0.5)}) of a GWAS credible ` +
      'set (for the disease: the genes of its GWAS regions; G*: significant only in a secondary analysis). C: curated ' +
      `evidence (source score ≥ ${pv('MendelianMin', 0.5)}; ClinVar alone, a pathogenic record of two stars or ` +
      'more). L: a positive finding of a genetic study of the disease.'],
    reports: [['model', 'rule'], 'A finding the readers took from an abstract, kept when its gene, variant and ' +
      'numbers are in the text (3.2).'],
    studies: [['model', 'pubmed', 'rule'], 'A record the readers include as a genetic study of the disease (3.2).'],
    region_of: [['profile', 'ot', 'rule'], 'A GWAS region of the disease (3.1).'],
    names_gene: [['profile', 'ot'], `The gene the publication names, else the top L2G gene (≥ ${pv('L2gMin', 0.5)}), ` +
      'else the position.'],
    shares_region: [['ot', 'rule'], 'A signal shared with a region of the disease (4.3).'],
    in_theme: [['profile', 'ot', 'rule'], 'The function of a gene (3.4).'],
    in_pathway: [['ot'], 'The Reactome top-level pathways of a text-mined target.'],
    mentions: [['epmc', 'ot', 'rule'], 'Europe PMC text mining as Open Targets gives it, corrected by the curation ' +
      '(3.3).'],
    gwas_of: [['ot'], 'The GWAS studies of the disease and of each related condition whose genetics are read (4.5).']};
  body.push(sec('8', 'graph', 'The graph', 'Typed nodes and edges; every edge keeps its evidence and its source.',
    k('Nodes'), plainTable(['Type', 'Meaning', S], Object.keys(D.schema.nodes).map(t => [h('span', {}, dot(t), ' ',
      h('code', {}, t)), D.schema.nodes[t], h('span', {class: 'num'}, num0(byType(t).length + (t === 'condition' ?
        1 : 0)))])),
    k('Edges'), plainTable(['Type', 'Meaning', 'From', 'How', S], Object.keys(D.schema.edges).map(t => [
      h('code', {}, t), D.schema.edges[t], h('span', {}, (EDGE[t] || [[]])[0].map(tag)), (EDGE[t] || [, ''])[1],
      h('span', {class: 'num'}, num0(nEdges(t)))])),
    here(`${num0(m.n_nodes)} nodes, ${num0(m.n_edges)} edges.`)));

  // 9 reproducibility and files
  const file = (part, where, what) => [part, h('code', {}, where), what];
  body.push(sec('9', 'snapshots', 'Reproducibility and files', 'A snapshot is the graph of one date: it keeps what it ' +
    'downloaded, its model answers are cached, and every file has one job.',
    h('ul', {class: 'doc-list'},
      h('li', {}, `This is the graph of ${S} as of ${m.snapshot}: Open Targets release ${m.release}, PubMed to ` +
        `${L.fetch_date || m.snapshot}, every publication checked in PubMed on ${m.checked_at}.`),
      h('li', {}, 'Every download is saved with the snapshot (the Open Targets answers and the PubMed records), so ' +
        'that it can be rebuilt without fetching again.'),
      h('li', {}, 'Every model answer is cached by model and prompt, so a rebuild with the same cache gives the same ' +
        'graph. When a prompt changes (a new rule, or a list of conditions that changed), its answers are asked ' +
        'again and a few decisions can change; each is listed with its reason.'),
      h('li', {}, 'A later snapshot of the disease is compared with the one before it: genes, conditions and studies ' +
        'added, removed or changed.')),
    k('Files'), plainTable(['Part', 'Where', 'What it holds'], [
      file('Viewer', `${D.disease.key}_knowledge_graph.html`, 'this page: one file that works offline'),
      file('Workbook', `${D.disease.key}_knowledge_graph_${m.snapshot}.xlsx`, 'the tables by question, sources in ' +
        'cell notes; a summary version without the flat appendix tables'),
      file('Knowledge base', '08.knowledge', 'the graph as Markdown text, one section per question'),
      file('Graph', '05.graph', 'nodes and edges as TSV, graph.json (NetworkX) and SQLite, with the schema'),
      file('Tables', '04.tables', 'the tables the graph is built from'),
      file('Figures', '09.figures', 'the figures, with a catalogue'),
      file('Readings and decisions', '01.relations, 02.literature, 03.curated', 'the models\' readings, the disputed ' +
        'items, the consolidation and verification decisions, the curation log'),
      file('Changes', '07.changes', 'the comparison with the previous snapshot, when there is one'),
      file('Run record', '_run_info', 'the parameters, the models and the trace of the run')]),
    rules('rebuild', 'offline', 'smoke')));

  // 10 the parameter index: every rule no section names, then the whole registry by group
  const rest = [...RULES.values()].filter(r => !used.has(r.key));
  body.push(sec('10', 'rules', 'Parameter index', `Every rule of this build with its value and basis: those named in ` +
    'the sections above, by group, and any other.',
    rest.length ? [k('Rules no section above names'), ruleRows(rest)] : [],
    D.choices.map(g => h('details', {class: 'cmap-all'}, h('summary', {}, `${g.name} (${g.items.length})`),
      ruleRows(g.items)))));

  // the contents column (a div: the page's nav and the ladder's level classes have styles of their own), in the
  // order of the section numbers (a subsection is made before the section that holds it)
  const rank = no => no.split('.').map(Number);
  TOC.sort((a, b) => rank(a[0])[0] - rank(b[0])[0] || (rank(a[0])[1] || 0) - (rank(b[0])[1] || 0));
  // sections, each with its subsections nested under it
  const link = (no, key, title) => h('a', {href: '#/build/' + key, 'data-key': key}, h('span', {class: 'doc-no'}, no),
    h('span', {}, title));
  const nav = h('div', {class: 'doc-nav'}, h('div', {class: 'doc-nav-h'}, 'Contents'), h('ol', {class: 'toc-l1'},
    TOC.filter(t => t[3] === 1).map(([no, key, title]) => {
      const kids = TOC.filter(t => t[3] === 2 && t[0].split('.')[0] === no);
      return h('li', {}, link(no, key, title), kids.length ? h('ol', {class: 'toc-l2'}, kids.map(([n2, k2, t2]) =>
        h('li', {}, link(n2, k2, t2)))) : null);
    })));
  // the contents mark the section in view, and the section that holds it
  const spy = () => {
    let cur = null;
    for (const [, key] of TOC) {
      const el = document.getElementById('b-' + key);
      if (el && el.getBoundingClientRect().top < 140) cur = key;
    }
    const no = (TOC.find(t => t[1] === cur) || [''])[0], top = (TOC.find(t => t[0] === no.split('.')[0]) || [])[1];
    nav.querySelectorAll('a').forEach(a => {
      a.classList.toggle('on', a.dataset.key === cur);
      a.classList.toggle('in', a.dataset.key === top && top !== cur);
    });
  };
  let ticking = false;
  const onScroll = () => {
    if (!document.body.contains(nav)) { window.removeEventListener('scroll', onScroll); return; }
    if (!ticking) { ticking = true; requestAnimationFrame(() => { ticking = false; spy(); }); }
  };
  window.addEventListener('scroll', onScroll, {passive: true});
  // an address #/build/<section> scrolls to that section
  requestAnimationFrame(() => {
    const d = sel && document.getElementById('b-' + sel);
    if (d) d.scrollIntoView();
    spy();
  });
  return [h('h2', {}, 'How this knowledge graph is built'), h('p', {class: 'lead'}, 'The methods of this graph, ' +
    `section by section in the order of the build, each with its rules and the counts of ${S}.`),
    glance(), h('div', {class: 'doc'}, nav, h('div', {class: 'doc-body'}, ...body))];
}

// the build at a glance: sources → methods → graph
function glance() {
  const L = D.llm, m = D.meta, cc = D.curation, R = L.relations, X = L.extract || {};
  const num0 = v => (v ?? 0).toLocaleString();
  const regions = byType('region'), nPub = regions.filter(published).length, conds = related();
  const model = r => (L.models.find(x => x.role === r) || {}).model || '';
  const item = (kk, title, text) => h('div', {class: 'gl-item'}, h('div', {class: 'gl-t'}, tag(kk), ' ', title),
                                      h('div', {class: 'gl-d'}, text));
  const nodeRow = (t, label, v) => h('div', {class: 'gl-n'}, dot(t), h('span', {}, label),
                                     h('b', {}, (v ?? byType(t).length).toLocaleString()));
  return h('div', {class: 'glance'},
    h('div', {class: 'gl-col'}, h('div', {class: 'gl-h'}, 'Sources'),
      item('ot', `release ${m.release}`, 'disease and ontology · GWAS credible sets, L2G, colocalization · curated ' +
        'gene evidence · text-mined targets'),
      item('pubmed', `to ${L.fetch_date || m.snapshot}`, `${num0(X.records)} genetic-study records` +
        (R ? ` · ${num0(R.records)} relation records` : '') + ` · ${num0(cc.counts.publications)} publications checked`),
      item('epmc', '', 'methods texts of open articles · sentences of text-mined mentions'),
      item('profile', D.disease.id + (L.evidence_also.length ? ' + ' + L.evidence_also.map(x => x.id).join(', ') : ''),
        [L.evidence_also.length ? 'evidence also from ' + L.evidence_also.map(x => x.name).join(', ') : '',
         L.curated_relations ? 'curated relations' : '', nPub ? 'published GWAS regions' : '',
         L.themes === 'curated' ? 'gene functions' : ''].filter(Boolean).join(' · ') || 'ID only: every rule by default')),
    h('div', {class: 'gl-arrow'}, '→'),
    h('div', {class: 'gl-col'}, h('div', {class: 'gl-h'}, 'Methods'),
      item('rule', 'Rules', 'search and pre-screen · curation · thresholds · evidence ladder · GWAS regions · ' +
        'shared signals · categories'),
      item('model', 'Two readers, one adjudicator', `${model('Reader A')} and ${model('Reader B')} read each record ` +
        `apart; ${model('Adjudicator')} decides where they differ`),
      item('rule', 'Source-text check', 'every gene, variant, number and quote must be found in its text'),
      h('div', {class: 'gl-foot'}, `${num0(X.read_by_model)} records read` + (R ? `, ${num0(R.read ?? R.records)} ` +
        'for relations' : '') + (agreeShare() ? ` · ${agreeShare()} of decisive items agreed by both readers` : ''))),
    h('div', {class: 'gl-arrow'}, '→'),
    h('div', {class: 'gl-col'}, h('div', {class: 'gl-h'}, 'Knowledge graph'),
      nodeRow('condition', 'Conditions', conds.length + 1), nodeRow('gene', 'Genes'),
      nodeRow('region', `${S} GWAS regions`), nodeRow('publication', 'Publications'),
      nodeRow('gwas_study', 'GWAS studies'), nodeRow('theme', 'Gene functions'), nodeRow('pathway', 'Reactome pathways'),
      h('div', {class: 'gl-foot'}, `${num0(m.n_edges)} edges of ${Object.keys(D.schema.edges).length} types, each ` +
        'with its source · viewer, workbook, knowledge base, SQLite')));
}
// the share of the decisive items both readers gave the same answer (studies, findings, methods, relations)
function agreeShare() {
  const X = D.llm.extract || {}, A = X.agreement || {}, R = D.llm.relations;
  const os = Object.values(A).concat(R && R.agreement ? [R.agreement] : []);
  const both = os.reduce((s, o) => s + (o.both || 0), 0), tot = os.reduce((s, o) => s + (o.both || 0) +
    (o.adjudicated || 0), 0);
  return tot ? `${Math.round(100 * both / tot)}%` : '';
}
const SHOWN_CAT = c => c === D.shared_signal ? 'Shared GWAS signal only' : c;

function modelSections() {
  const L = D.llm, X = L.extract || {}, A = X.agreement || {}, R = L.relations, Sm = L.summaries;
  const two = L.models.length >= 3;
  const n = v => (v ?? 0).toLocaleString();
  const ag = o => o ? `${n(o.both)} agreed, ${n(o.adjudicated)} adjudicated` : '';
  const out = [h('p', {class: 'llm-rules'}, 'Language models read the literature; every rule of the graph is ' +
    'deterministic. Each reading is done twice, by two models that never see each other\'s answer; a third model ' +
    'decides where they disagree; and every item a model returns must be found in its source text.')];
  // the models
  out.push(h('h4', {}, two ? 'The three models' : 'The model'),
    h('div', {class: 'cards'}, L.models.map(m => h('div', {class: 'card'}, h('div', {class: 't'}, m.role),
      h('div', {class: 'n', style: 'font-size:17px'}, m.model), h('div', {class: 't'}, m.backend)))),
    h('div', {class: 'llm-rules'},
      h('p', {}, 'Each model is named rather than left to a backend\'s default, and the adjudicator is a third model: ' +
        'the pipeline refuses one that is also a reader. The backend is the user\'s choice (a command-line tool or ' +
        'an API); keys are read from the environment only.'),
      h('p', {}, 'Every answer is cached by model and prompt, so a rebuild with the same cache gives the same graph ' +
        'and a later snapshot asks only about new or changed records. A model sees only the record and its ' +
        'instructions: no file, memory or instruction of the user\'s; an e-mail address that is not in the source ' +
        'text is removed from every answer before it is cached.' +
        (D.params.LlmBatch ? ` Abstracts are read ${D.params.LlmBatch} to a call.` : '')),
      two ? null : here('one reader only, so no agreement is measured and nothing is adjudicated.')));
  // one reading
  out.push(h('h4', {}, 'One reading, from record to graph'), flowDiagram(L.models),
    h('div', {class: 'note'}, 'The two readers read the same record apart. Every gene, variant, number and quote of ' +
      'each reading must occur in the source text, or the item is dropped and logged. A rule then compares the items ' +
      'that decide what the graph counts; where they agree, the reading goes on; where they differ, the adjudicator ' +
      'reads the source text with both readings (shown in an order that alternates by record, so that neither reader ' +
      'is always first) and decides, and its reading is checked against the text in the same way. Summaries take ' +
      'another path: reader A writes, reader B checks each claim, the adjudicator rules on each challenged claim.'));
  // where models are used
  const step = (title, does, checks, mine) => h('div', {class: 'llm-step'}, h('h4', {}, title),
    h('p', {}, h('b', {}, 'The models: '), does), h('p', {}, h('b', {}, 'The checks: '), checks), here(...[].concat(mine)));
  const rel = L.curated_relations ? ['not used: the clinical relations come from the curated list of its profile.'] :
    !R ? ['read from PubMed; the counts of the reading are not recorded in this snapshot.'] :
    [`${n(R.read ?? R.records)} records read` + (R.not_read ? ` (${n(R.not_read)} not read: no language model ` +
      'was called, and no answer was cached)' : '') + `; ${n(R.relations)} relations kept (${ag(R.agreement)}) for ` +
    `${n(R.conditions)} conditions, ${n(R.mapped_to_open_targets)} of them in Open Targets; ${n(R.dropped_items)} items ` +
    `dropped by the source-text check; ${n(R.disputed_items)} disputed items.`];
  out.push(h('h4', {}, 'Where models are used, and how each answer is checked'), h('div', {class: 'llm-steps'},
    step('Genetic studies: screening and extraction',
      `decide from the abstract whether a record is a genetic study of the disease, and extract its findings (gene, ` +
      'variant, result, effect, P), its sample sizes and whether it was done in Japan or East Asia.',
      'the readers are compared on inclusion, each finding (gene and result), the sample sizes and the population ' +
      'flags; the adjudicator decides the rest. Every gene, variant and number must be in the abstract; a result ' +
      'label that the wording of its quote contradicts is flagged. A pre-screen by rules (record types, genetic ' +
      'terms) comes first.',
      `${n(X.records)} records found; ${n(X.prescreened_out)} left out by the pre-screen; ${n(X.read_by_model)} read` +
      (X.not_read ? ` (${n(X.not_read)} not read: no language model was called, and no answer was cached)` : '') + ', ' +
      `inclusion ${ag(A.inclusion)}; ${n(X.included)} studies included; ${n(X.findings)} findings ` +
      `(${ag(A.findings)}); ${n(X.dropped_items)} items dropped by the source-text check; ${n(X.result_flags)} result ` +
      'labels flagged.'),
    step('Genetic studies: methods',
      'quote, from the methods text where it is open (Europe PMC), else the abstract, how the study was done: the ' +
      'diagnosis, cohort, controls, genotyping or sequencing, quality control, statistical tests and replication.',
      'two quotes of a field agree when they overlap; the adjudicator decides the rest; every quote must be in the ' +
      'text.',
      `${n(X.methods_fields)} methods fields (${ag(A['methods fields'])}).`),
    step('Clinical relations, for a disease without a curated list',
      'read guidelines, consensus statements, meta-analyses, systematic reviews, cohort and case-control studies, ' +
      'and give each condition the record relates to the disease, with its category, direction and a verbatim quote.',
      'the readers are compared on the condition, its direction and its category; the adjudicator decides the rest. ' +
      'The quote must be in the abstract and name the condition; only evidence from human patients or populations ' +
      'counts (not a sign or finding of the disease, a mechanism, an animal or cell experiment, or a comparison ' +
      'between groups of patients). Each condition is mapped to Open Targets by its exact name.', rel),
    step('Summaries',
      'write a short summary of each mechanism theme, of the findings and of the methods of the genetic studies, ' +
      'citing only the publications they were given.',
      'reader B checks every claim against the material; the adjudicator rules on each claim reader B challenged; ' +
      'a summary with a confirmed unsupported claim is withheld, with its reason, and every cited PMID must be among ' +
      'its sources.',
      `${n(Sm.n)} summaries: ${n(Sm.shown)} shown, ${n(Sm.withheld)} withheld after the checks` +
      (Sm.not_generated ? `, ${n(Sm.not_generated)} not generated in this build` : '') + '.')));
  // how often the readers agreed
  if (two) {
    const rows = [['Inclusion of a record', A.inclusion], ['Findings', A.findings], ['Methods fields', A['methods fields']]]
      .concat(R ? [['Clinical relations', R.agreement]] : []).filter(([, o]) => o && (o.both || o.adjudicated));
    const bar = (label, o) => {
      const tot = (o.both || 0) + (o.adjudicated || 0), w = x => (100 * x / Math.max(tot, 1)).toFixed(1) + '%';
      return h('div', {class: 'agree-row'}, h('div', {class: 'agree-k'}, label),
        h('div', {class: 'agree-bar'}, h('span', {class: 'both', style: 'width:' + w(o.both || 0)}),
          h('span', {class: 'adj', style: 'width:' + w(o.adjudicated || 0)})),
        h('div', {class: 'agree-v'}, `${n(o.both)} agreed · ${n(o.adjudicated)} adjudicated ` +
          `(${Math.round(100 * (o.both || 0) / Math.max(tot, 1))}% agreed)`));
    };
    out.push(h('h4', {}, 'How often the two readers agreed'), h('div', {class: 'agree'}, rows.map(([l, o]) => bar(l, o))),
      h('div', {class: 'cmap-key', style: 'justify-content:flex-start'},
        h('span', {}, h('span', {class: 'swatch', style: 'background:#5b8c6e'}), 'both readers gave the same answer'),
        h('span', {}, h('span', {class: 'swatch', style: 'background:#c9a25b'}), 'they differed; the adjudicator decided')),
      h('div', {class: 'note'}, 'Agreement is counted on the items that decide what the graph counts, as the ' +
        'readers returned them (before the curation leaves out retracted articles and conference material; the ' +
        'workbook counts the items kept after it). Inclusion: a record counts as adjudicated when the readers ' +
        'differed on whether to include it, not on another of its items. A relation that one reader names ' +
        'differently from the other counts as a difference.'));
  }
  // what no model decides
  out.push(h('h4', {}, 'What no model decides'), h('p', {class: 'llm-rules'}, 'The pre-screen of the records; every ' +
    'curation rule (publications left out, the place and design of a study, variant names, what a finding is about, ' +
    'text-mining corrections); the evidence ladder; the GWAS regions, the ontology neighbors and the shared signals; ' +
    'the related conditions of a curated list and their categories; the graph, its figures, its tables and this ' +
    'viewer. A model\'s answer enters the graph only through the checks above.'));
  // every disputed item
  const disputed = L.review.filter(r => r[2] !== 'result check'), flagged = L.review.filter(r => r[2] === 'result check');
  if (disputed.length) {
    let pass = 'all';
    const box2 = h('div');
    const passes = [...new Set(disputed.map(r => r[1]))];
    const draw = () => {
      const rows = disputed.filter(r => pass === 'all' || r[1] === pass)
        .map(([pm, ps, item, a, b, fin, dec]) => ({pm, ps, item, a, b, fin, dec}));
      box2.replaceChildren(table([
        {label: 'Source', get: r => r.pm || 'Condition list', cell: r => r.pm ? pmid(r.pm) : 'Condition list'},
        {label: 'Pass', get: r => r.ps},
        {label: 'Item', get: r => r.item}, {label: 'Reader A', get: r => r.a}, {label: 'Reader B', get: r => r.b},
        {label: 'Final', get: r => r.fin}, {label: 'Decision', get: r => r.dec, cell: r => h('span', {class: 'quote'}, r.dec)}],
        rows, {tools: h('span', {}, ['all'].concat(passes).map(x => h('button', {class: 'btn' + (pass === x ?
          ' on' : ''), onclick: () => { pass = x; draw(); }}, x === 'all' ? 'Every pass' : x)))}));
    };
    draw();
    out.push(h('h4', {}, `Every disputed item (${disputed.length.toLocaleString()})`),
      h('div', {class: 'note'}, 'Each item on which the two readers differed, with both readings, the final reading and ' +
        'the adjudicator\'s reason (the workbook lists them too).'), box2);
  }
  if (flagged.length) {
    out.push(h('h4', {}, `Results flagged against their source sentence (${flagged.length.toLocaleString()})`),
      h('div', {class: 'note'}, 'A rule, not a reader: the result label of a finding that the wording of its own ' +
        'sentence contradicts. The finding is kept and marked.'),
      table([{label: 'Source', get: r => r[0], cell: r => pmid(r[0])}, {label: 'Finding', get: r => r[5]},
             {label: 'Why', get: r => r[6]}], flagged));
  }
  return out;
}

// ── routing, search ──────────────────────────────────────────────────────
const TABS = [['overview', 'Overview'], ['genes', 'Genes'], ['conditions', 'Related conditions'],
              ['regions', 'GWAS regions'], ['studies', 'Genetic studies'], ['curation', 'Curation'],
              ['build', 'How it is built']];
const PAGES = {overview, genes: genesPage, conditions: conditionsPage, regions: regionsPage, studies: studiesPage,
               curation: curationPage, build: buildPage};
// the tab: the disease first, then what the page shows (a list, or one node by its kind), so that the open tabs of
// several diseases and pages tell apart even when the browser shortens their titles; its icon is the disease's mark
const TITLE_KIND = {...TYPE_NAME, region: 'GWAS region'};
function pageTitle(hash, tab) {
  if (LIST_PAGE) return (TABS.find(x => x[0] === (tab || 'overview')) || [0, 'Overview'])[1];
  const n = N.get(hash.slice(2));
  if (!n) return 'Not found';
  const label = n.type === 'region' ? n.label.replace(/ region$/, '') : n.type === 'publication' &&
    n.id.startsWith('PMID:') ? `${n.label} (PMID ${n.id.slice(5)})` : n.label;
  return `${TITLE_KIND[n.type] || n.type}: ${label}`;
}
const MARK = document.getElementById('mark');
MARK.textContent = D.disease.mark; MARK.style.background = D.disease.color; MARK.title = D.disease.Name;
document.head.append(h('link', {rel: 'icon', href: 'data:image/svg+xml,' + encodeURIComponent(
  '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 32 32"><rect width="32" height="32" rx="7" fill="' +
  D.disease.color + '"/><text x="16" y="21.5" font-family="Arial,Helvetica,sans-serif" font-size="14" ' +
  'font-weight="700" fill="#fff" text-anchor="middle">' + D.disease.mark + '</text></svg>')}));
function render() {
  const hash = decodeURIComponent(location.hash.slice(2) || 'overview');
  const main = document.getElementById('main');
  PAGED.length = 0;
  let tab = hash.split('/')[0];
  LIST_PAGE = !hash.startsWith('n/');
  if (!LIST_PAGE) { main.replaceChildren(nodePage(hash.slice(2))); tab = ''; }
  else main.replaceChildren(...[].concat((PAGES[tab] || overview)(hash.split('/').slice(1).join('/'))));
  document.title = `${S} · ${pageTitle(hash, tab)} — DiseaseLens`;
  document.querySelectorAll('nav a').forEach(a => a.classList.toggle('on', a.dataset.tab === tab));
  window.scrollTo(0, 0);
  requestAnimationFrame(fitHints);
}
const TAB_COUNT = {genes: byType('gene').filter(g => lvl(g) !== null).length, conditions: related().length,
                   regions: byType('region').length, studies: byType('publication').filter(p => p.a.design_class).length};
document.getElementById('nav').append(...TABS.map(([k, l]) => h('a', {href: '#/' + k, 'data-tab': k}, l,
  TAB_COUNT[k] !== undefined ? h('span', {class: 'cnt'}, TAB_COUNT[k].toLocaleString()) : null)));
// the full name and the ID beside the short name, which can name other things too; the snapshot as chips
document.getElementById('meta').textContent = (D.disease.name.toLowerCase() === S.toLowerCase() ? '' :
  `${D.disease.Name} · `) + D.disease.id;
document.getElementById('chips').append(...[`Snapshot ${D.meta.snapshot}`, `Open Targets ${D.meta.release}`,
  `PubMed checked ${D.meta.checked_at}`, `${D.meta.n_nodes.toLocaleString()} nodes · ` +
  `${D.meta.n_edges.toLocaleString()} edges`].map(x => h('span', {}, x)));
document.getElementById('foot').append(`DiseaseLens · ${S} knowledge graph · snapshot ${D.meta.snapshot} · Open ` +
  `Targets ${D.meta.release}. Links to sources (↗) and to a gene, condition or study from a list open in a new tab. `,
  h('a', {href: '#/build'}, 'How it is built'), '.');
const INDEX = [...N.values()].map(n => ({n, l: n.label.toLowerCase(), i: n.id.toLowerCase()}));
const RANK = {gene: 0, condition: 1, region: 2, publication: 3, theme: 4, pathway: 5, gwas_study: 6};
const search = document.getElementById('search'), sug = document.getElementById('suggest');
let hits = [], cur = 0;
function showSuggest() {
  const q = search.value.trim().toLowerCase();
  if (!q) { sug.style.display = 'none'; return; }
  hits = INDEX.filter(x => x.l.includes(q) || x.i.includes(q)).sort((a, b) =>
    (b.l === q) - (a.l === q) || (b.l.startsWith(q)) - (a.l.startsWith(q)) || RANK[a.n.type] - RANK[b.n.type] ||
    a.l.length - b.l.length).slice(0, 14);
  cur = 0;
  sug.replaceChildren(...hits.map((x, k) => h('div', {class: k === cur ? 'on' : null, onmousedown: () => go(x.n.id)},
    dot(x.n.type), h('span', {class: x.n.type === 'gene' ? 'gene' : null}, x.n.label),
    h('span', {class: 'note', style: 'margin-left:auto'}, TYPE_NAME[x.n.type]))));
  sug.style.display = hits.length ? 'block' : 'none';
}
function go(id) { location.hash = href(id); search.value = ''; sug.style.display = 'none'; }
search.addEventListener('input', showSuggest);
search.addEventListener('keydown', e => {
  if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
    cur = (cur + (e.key === 'ArrowDown' ? 1 : hits.length - 1)) % Math.max(hits.length, 1);
    [...sug.children].forEach((c, k) => c.classList.toggle('on', k === cur)); e.preventDefault();
  } else if (e.key === 'Enter' && hits[cur]) go(hits[cur].n.id);
  else if (e.key === 'Escape') sug.style.display = 'none';
});
search.addEventListener('blur', () => setTimeout(() => { sug.style.display = 'none'; }, 150));
window.addEventListener('hashchange', render);
window.addEventListener('resize', fitHints);
// printing shows every rule of the condition map and every row of every table
window.addEventListener('beforeprint', () => { document.querySelectorAll('details.cmap-all, details.bd')
  .forEach(d => { d.open = true; }); PRINTING = true; PAGED.forEach(f => f()); });
window.addEventListener('afterprint', () => { PRINTING = false; PAGED.forEach(f => f()); });
render();
</script>
</body>
</html>
'''

if __name__ == '__main__':
    main()
