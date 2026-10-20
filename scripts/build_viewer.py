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
#             Conditions   the related conditions: an interactive map (the
#                          figure condition_map: each condition in the sector
#                          of its category, the line showing the kind of
#                          evidence; hover for its evidence, click to open it,
#                          a category to select it), with what it shows and leaves
#                          out, above the table, which shows the selection
#             Regions      the genome-wide significant GWAS regions of the
#                          disease (published, or from Open Targets); a sentence
#                          when there is none
#             Studies      the genetic studies of the disease, by place and
#                          design, with the agreement of the model readers
#             Curation     what was left out or corrected, and why
#             How it is built  the general rules, each followed by how it
#                          applies to this disease: the sources (databases), how
#                          the databases, the language models and the rules are
#                          combined, the associations (edge types: where each comes
#                          from, by which rule), the categories and levels, the
#                          language models (the three models, one reading from
#                          record to graph, where they are used and how each answer
#                          is checked, how often the readers agreed, every disputed
#                          item), and every rule with its value and basis
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
#           09.figures/condition_map.json (the condition map as drawn)
# Outputs : <key>_knowledge_graph.html
# ---------------------------------------------------------------------------
import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import (EAS_COUNTRIES, OTHER_GENES, SHARED_SIGNAL, load_choices,  # noqa: E402
                       load_profile, log, not_regions_of, not_regions_sentence, read_tsv)

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
    # the note under the ladder: the coded levels (GWAS, curated evidence), then the literature levels
    coded = [(k, c, lab) for k, (c, lab) in ladder.items() if c]
    lit = [k for k, (c, _lab) in ladder.items() if not c]
    parts = [f'level {k} ({c}): {lab[:1].lower()}{lab[1:]}' for k, c, lab in coded]
    note = '; '.join(parts)
    note = (note[:1].upper() + note[1:] + '. ') if note else ''
    if len(lit) >= 2:
        note += (f'Levels {lit[0]}–{lit[-1]}: the genetic studies of {s}; levels {lit[0]}–{lit[-2]} count only '
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
            'relations': ({k: rr[k] for k in ('records', 'read', 'not_read', 'relations', 'conditions',
                                              'mapped_to_open_targets',
                                              'agreement', 'disputed_items', 'dropped_items') if k in rr}
                          if rr else None),
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
    ap.add_argument('--figures', help='09.figures (its condition_map.json; without it, no condition map)')
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
    # the level of a secondary analysis of a published GWAS (G*) only for a disease that has one
    rg_codes = set(read_tsv(G / 'nodes_region.tsv', ['code']).code)
    if 'G*' not in rg_codes and not any(n[1] == 'gene' and str(n[4].get('disease_level')) == '2' for n in nodes):
        ladder = {k: v for k, v in ladder.items() if v[0] != 'G*'}
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
    cmap = Path(a.figures or '.') / 'condition_map.json'
    cmap = json.loads(cmap.read_text()) if a.figures and cmap.exists() else None
    run_params = json.loads(Path(a.run_params).read_text())
    choices, _refs = load_choices(P, run_params, T)
    data = {
        'meta': {**schema['graph'], 'checked_at': cur['checked_at'][:10]},
        'disease': {'id': P.id, 'key': P.key, 'name': P.name, 'Name': P.Name, 'short': P.short},
        'about': schema.get('about') or {'default': P.about, 'within': P.within},
        'schema': {'nodes': {k: v['meaning'] for k, v in schema['nodes'].items()},
                   'edges': {k: v['meaning'] for k, v in schema['edges'].items()}},
        'ladder': ladder, 'categories': cats, 'shared_signal': SHARED_SIGNAL, 'other_genes': OTHER_GENES,
        'eas': list(EAS_COUNTRIES),
        # the labels of the region statistics of a published table that only some diseases have
        'region_labels': {'p_pe': cols.get('p_pe') or 'P in the secondary analysis', 'v': cols.get('v') or 'V'},
        'text': texts(P, regions, ladder, run_params, not_regions), 'condition_map': cmap,
        'not_regions': not_regions.to_dict('records'),
        'llm': model_use(P, T, Q, summ), 'params': run_params,
        'choices': [{'name': g['name'], 'items': [{'item': it['item'], 'value': it['value_text'],
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
        f'condition map: {"yes" if cmap else "no"}, {Path(a.out).stat().st_size / 1e6:.1f} MB')


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
  --l1: #c6dbef; --l2: #e3eef8; --l3: #fce4cc; --l4: #fdf2e6; --l5: #f2f2f2;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--ink); font: 14px/1.45 -apple-system, "Segoe UI", Roboto,
       "Helvetica Neue", Arial, "Noto Sans", sans-serif; }
header { background: var(--head); color: #fff; padding: 12px 20px; display: flex; flex-wrap: wrap; gap: 12px 24px;
         align-items: center; }
header h1 { font-size: 18px; margin: 0; font-weight: 600; letter-spacing: .2px; }
header .meta { font-size: 12px; opacity: .85; }
#search-wrap { position: relative; margin-left: auto; min-width: 280px; flex: 0 1 380px; }
#search { width: 100%; padding: 7px 10px; border-radius: 4px; border: 0; font-size: 14px; }
#suggest { position: absolute; top: 36px; left: 0; right: 0; background: #fff; color: var(--ink); z-index: 10;
           border: 1px solid var(--line); box-shadow: 0 4px 14px rgba(0,0,0,.15); display: none; max-height: 420px;
           overflow: auto; }
#suggest div { padding: 6px 10px; cursor: pointer; display: flex; gap: 8px; align-items: center; }
#suggest div.on, #suggest div:hover { background: var(--chip); }
nav { display: flex; gap: 2px; padding: 0 20px; background: var(--panel); border-bottom: 1px solid var(--line);
      overflow-x: auto; }
nav a { padding: 10px 14px; color: var(--ink2); text-decoration: none; border-bottom: 3px solid transparent;
        white-space: nowrap; }
nav a.on { color: var(--head); border-bottom-color: var(--head); font-weight: 600; }
main { padding: 18px 20px 60px; max-width: 1400px; margin: 0 auto; }
h2 { font-size: 20px; margin: 4px 0 6px; } h3 { font-size: 15px; margin: 22px 0 8px; color: var(--head); }
p.lead { color: var(--ink2); margin: 0 0 14px; max-width: 980px; }
a { color: var(--link); } a.ext { white-space: nowrap; } a.ext::after { content: " \2197"; font-size: 11px; }
.cards { display: grid; grid-template-columns: repeat(auto-fill, minmax(190px, 1fr)); gap: 10px; margin: 10px 0 6px; }
.card { background: var(--panel); border: 1px solid var(--line); border-radius: 6px; padding: 10px 12px; }
.card .n { font-size: 22px; font-weight: 650; color: var(--head); } .card .t { font-size: 12px; color: var(--ink2); }
table { border-collapse: collapse; width: 100%; font-size: 13px; }
th, td { border-bottom: 1px solid var(--line); padding: 5px 7px; text-align: left; vertical-align: top; }
th { background: var(--panel); font-weight: 600; position: sticky; top: 0; cursor: pointer; white-space: nowrap; }
th.nosort { cursor: default; }
tr:hover td { background: #fafbfc; }
td.num { text-align: right; font-variant-numeric: tabular-nums; }
.quote { color: var(--ink2); font-style: italic; }
.badge { display: inline-block; font-size: 11px; padding: 1px 7px; border-radius: 10px; color: #fff;
         white-space: nowrap; }
.badge.adj { background: #8a6d3b; font-size: 10px; padding: 0 6px; }
.dot { display: inline-block; width: 9px; height: 9px; border-radius: 50%; flex: none; }
.gene { font-style: italic; }
.chip { display: inline-block; background: #fff; border: 1px solid var(--line); border-radius: 3px; padding: 0 5px;
        margin: 2px 3px 2px 0; font-style: italic; text-decoration: none; color: var(--ink); }
.chip:hover { border-color: var(--link); color: var(--link); }
.grid td { min-width: 80px; } .grid th { cursor: default; }
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
.more { margin-top: 6px; }
.summary { background: var(--panel); border-left: 3px solid var(--t-theme); padding: 8px 12px; margin: 8px 0;
           font-size: 13px; }
/* the condition map: hover a condition for its evidence, click to open it; a category selects its conditions */
.cmap-wrap { position: relative; max-width: 880px; margin: 4px auto 0; }
.cmap { display: block; width: 100%; height: auto; overflow: visible; }
.cmap text { font-size: 13px; }
.cmap .it { cursor: pointer; transition: opacity .15s; }
.cmap .it.dim { opacity: .15; }
.cmap .it:hover text { font-weight: 700; }
.cmap .it:hover line, .cmap .it.hl line { stroke-opacity: 1; }
.cmap .core { cursor: pointer; }
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
.toc { display: flex; flex-wrap: wrap; gap: 6px 14px; margin: 0 0 6px; font-size: 13px; }
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
.print-only { display: none; }
/* print, or save as PDF: the current page as a document */
@media print {
  @page { margin: 12mm; }
  body { font-size: 9.5pt; }
  nav, #search-wrap, #suggest, .tools, .btn, .cmap-tip { display: none !important; }
  header { background: none; color: var(--ink); padding: 0 0 4px; border-bottom: 2px solid var(--head); }
  header h1 { color: var(--head); }
  header .meta { opacity: 1; color: var(--ink2); }
  .cmap-wrap { max-width: 16cm; }
  .cmap-all { display: block; } .cmap-all > summary { list-style: none; }
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
  .print-only { display: block; }
  .llm-step, .agree-row { break-inside: avoid; page-break-inside: avoid; }
  .lv1, .lv2, .lv3, .lv4, .lv5, .lv6, .lv7, .badge, .dot, .card, th, svg, .agree-bar span, .swatch {
    -webkit-print-color-adjust: exact; print-color-adjust: exact; }
}
</style>
</head>
<body>
<header>
  <h1>__TITLE__</h1>
  <div class="meta" id="meta"></div>
  <div id="search-wrap"><input id="search" type="search" placeholder="Search a gene, condition, PMID, region…"
       autocomplete="off"><div id="suggest"></div></div>
</header>
<nav id="nav"></nav>
<main id="main"></main>
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
const nodeLink = (id, text) => {
  const n = N.get(id);
  if (!n) return text || id;
  return h('a', {href: href(id), class: n.type === 'gene' ? 'gene' : null, title: TYPE_NAME[n.type] + ': ' + n.id},
           text || n.label);
};
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
  'adjudicator decided'}, 'adjudicated') : (v || '');

// a sortable, filterable table: cols [{label, get(row) -> value for sort/filter, cell(row) -> node}]
function table(cols, rows, opts = {}) {
  let sortBy = opts.sort ?? null, dir = opts.dir ?? 1, q = '', limit = opts.limit || 200;
  const wrap = h('div');
  const input = opts.filter === false || rows.length <= 8 ? null :
    h('input', {type: 'search', placeholder: 'Filter these rows…', oninput: e => { q = e.target.value.toLowerCase();
      draw(); }});
  const count = h('span', {class: 'note'});
  if (input || opts.tools) wrap.append(h('div', {class: 'tools'}, input, opts.tools || null, count));
  const tbl = h('table'), body = h('tbody'), more = h('div', {class: 'more'});
  const head = h('tr');
  cols.forEach((c, i) => head.append(h('th', {class: c.get ? null : 'nosort', onclick: () => {
    if (!c.get) return; dir = sortBy === i ? -dir : (c.num ? -1 : 1); sortBy = i; draw(); }}, c.label)));
  tbl.append(h('thead', {}, head), body); wrap.append(tbl, more);
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
    body.replaceChildren(...rs.slice(0, limit).map(r => h('tr', {}, cols.map(c =>
      h('td', {class: c.num ? 'num' : null}, c.cell ? c.cell(r) : c.num ? fmt(c.get(r)) : linkIds(fmt(c.get(r))))))));
    count.textContent = rs.length === rows.length ? `${rows.length} rows` : `${rs.length} of ${rows.length} rows`;
    more.replaceChildren(...(rs.length > limit ? [h('button', {class: 'btn', onclick: () => { limit += 500; draw(); }},
      `Show more (${rs.length - limit} hidden)`), h('div', {class: 'note print-only'}, `The first ${limit} of ` +
      `${rs.length} rows; the viewer lists every row (“Show more”).`)] : []));
  }
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
  const labeled = others.length <= 36, nR = Math.ceil(others.length / 2), nL = others.length - nR;
  others.forEach(([o], i) => {
    // labeled: right half top to bottom, then left half bottom to top, evenly spaced in height, so labels never
    // collide; unlabeled: a full circle
    const right = i < nR, m = right ? nR : nL, k = right ? i : nL - 1 - (i - nR);
    const ang = -Math.PI / 2 + 2 * Math.PI * i / others.length;
    const dy = labeled ? R * 0.94 * (2 * (k + 0.5) / m - 1) : R * Math.sin(ang);
    const dx = labeled ? (right ? 1 : -1) * Math.sqrt(R * R - dy * dy) : R * Math.cos(ang);
    const x = cx + dx, y = cy + dy, n = N.get(o);
    if (!n) return;
    svg.append(el('line', {x1: cx, y1: cy, x2: x, y2: y, stroke: '#c9ced6', 'stroke-width': 1}));
    const g = el('a', {href: href(o)});
    const c = el('circle', {cx: x, cy: y, r: 7, fill: COLOR(n.type), stroke: '#fff', 'stroke-width': 1.5});
    const tt = el('title', {}); tt.textContent = `${TYPE_NAME[n.type]}: ${n.label}`; c.append(tt); g.append(c);
    if (labeled) {
      const lab = n.label.length > 26 ? n.label.slice(0, 25) + '…' : n.label;
      const t = el('text', {x: x + (right ? 10 : -10), y: y + 4, 'text-anchor': right ? 'start' : 'end'});
      t.textContent = lab; if (n.type === 'gene') t.setAttribute('font-style', 'italic'); g.append(t);
    }
    svg.append(g);
  });
  svg.append(el('circle', {cx, cy, r: 13, fill: COLOR(me.type), stroke: '#fff', 'stroke-width': 2}));
  const t = el('text', {x: cx, y: cy + 30, 'text-anchor': 'middle', 'font-weight': 700});
  t.textContent = me.label.length > 40 ? me.label.slice(0, 39) + '…' : me.label; svg.append(t);
  const types = [...new Set(others.map(([o]) => N.get(o)?.type).filter(Boolean))];
  return h('div', {}, svg, h('div', {class: 'legend'}, types.map(x => h('span', {}, dot(x), TYPE_NAME[x]))),
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
  const cols = [{label: sample ? TYPE_NAME[sample.type] : 'Node', get: e => N.get(otherOf(e))?.label,
                 cell: e => nodeLink(otherOf(e))}]
    .concat(ocols.map(c => ({label: c.label, num: c.kind === 'num', get: e => N.get(otherOf(e))?.a[c.key],
                             cell: e => fmt(N.get(otherOf(e))?.a[c.key])})))
    .concat(ecols.map(c => ({label: c.label, num: c.kind === 'num', get: e => e.a[c.key],
                             cell: e => cellOf(c.kind, e.a[c.key])})));
  // a short section is kept on one printed page (heading, note and rows together)
  return h('section', {class: list.length <= 10 ? 'short' : null},
           h('h3', {}, `${SECTION[key] || key} (${list.length})`),
           h('div', {class: 'note'}, D.schema.edges[type] || ''), table(cols, list, {limit: 100}));
}

// ── facts of a node ───────────────────────────────────────────────────────
const FACTS = {
  gene: [['disease_level', `${S} evidence level`, v => `${v} — ${(D.ladder[v] || [])[1] || ''}` +
          ((D.ladder[v] || [])[0] ? ` (${D.ladder[v][0]})` : '')], ['theme', 'Function'],
         ['n_related_conditions', 'Related conditions with genetic evidence'],
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
    const row = h('tr', {class: 'lv' + k}, h('td', {}, h('b', {}, `${k} ${code}`), h('br'), label));
    for (const t of th) {
      const gs = genes.filter(g => lvl(g) === Number(k) && (g.a.theme || D.other_genes) === t)
                      .sort((a, b) => a.label.localeCompare(b.label));
      row.append(h('td', {}, gs.map(g => h('a', {class: 'chip', href: href(g.id)}, g.label))));
    }
    tb.append(row);
  }
  grid.append(tb);
  out.push(h('div', {style: 'overflow-x:auto'}, grid), h('div', {class: 'note'}, D.text.ladder_note));
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
    'row here (Related conditions page; in the condition map, the gray node).' : single ? `${single} traits related ` +
    'only by a shared GWAS signal share one region each and are not rows here (Related conditions page; the ' +
    'condition map shows them as one gray node).' : '';
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
  return h('div', {}, h('div', {style: 'overflow-x:auto'}, tbl),
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
      {label: 'Evidence', get: g => g.a.disease_level_label},
      {label: 'Function', get: g => g.a.theme},
      {label: 'Related conditions', num: true, get: g => g.a.n_related_conditions},
      {label: 'Text-mined publications', num: true, get: g => g.a.n_text_mined_publications}],
      rows, {sort: all ? 4 : 1, dir: all ? -1 : 1, limit: 300, tools: h('span', {},
        h('button', {class: 'btn' + (all ? '' : ' on'), onclick: () => { all = false; draw(); }}, `${S} genes`), ' ',
        h('button', {class: 'btn' + (all ? ' on' : ''), onclick: () => { all = true; draw(); }}, 'Every gene'))}));
  };
  draw();
  return [h('h2', {}, 'Genes'), h('p', {class: 'lead'}, `The genes studied genetically in ${S}, by evidence level; ` +
    `“Every gene” adds the genes with genetic evidence in a related condition and the text-mined ${S} targets.`), box];
}

// ── the condition map: the figure condition_map as drawn (09.figures/condition_map.json), interactive ──────
const SVGNS = 'http://www.w3.org/2000/svg';
function svgEl(tag, at, ...kids) {
  const x = document.createElementNS(SVGNS, tag);
  for (const [k, v] of Object.entries(at || {})) if (v !== null && v !== undefined) x.setAttribute(k, v);
  for (const k of kids.flat()) if (k !== null && k !== undefined) x.append(k);
  return x;
}
// a label in at most `max` lines of about w characters (the full name is in the tooltip)
function labelLines(s, w = 30, max = 2) {
  if (max === 1) return [String(s).length > w + 4 ? String(s).slice(0, w + 3).trimEnd() + '…' : String(s)];
  const lines = [];
  for (const word of String(s).split(' ')) {
    const k = lines.length - 1;
    if (k >= 0 && (lines[k] + ' ' + word).length <= w) lines[k] += ' ' + word; else lines.push(word);
  }
  return lines.length <= 2 ? lines : [lines[0], (lines.slice(1).join(' ')).slice(0, w - 1) + '…'];
}
const sectorLabel = (M, c) => (M.sectors.find(s => s.category === c) || {label: c}).label;

function conditionMap(M, focus, pick) {
  const R = 250, V = 520;                                    // the radius of the nodes; half the view box
  const xy = (a, r) => [+(r * Math.cos(a)).toFixed(1), +(-r * Math.sin(a)).toFixed(1)];
  const color = Object.fromEntries(M.sectors.map(s => [s.category, s.color]));
  // labels fit the room each condition has on the circle (any number of conditions): the angle of one slot (the
  // figure leaves two empty slots between sectors) at the label radius gives the room across a label; the font
  // shrinks to it (8 to 13), and a label keeps two lines only where two fit, else one line cut with an ellipsis
  const nShown = M.items.length + (M.single ? 1 : 0), nGroups = M.sectors.length + (M.single ? 1 : 0);
  const room = 2 * Math.PI / Math.max(nShown + 2 * nGroups, 1) * 1.17 * R;
  const fs = Math.max(8, Math.min(13, room / 1.25)), lineMax = room >= 2.3 * fs ? 2 : 1;
  const svg = svgEl('svg', {viewBox: `${-V} ${-V} ${2 * V} ${2 * V}`, class: 'cmap', role: 'img',
                            'aria-label': `${M.short} and its related conditions`});
  // a map of many conditions takes more of the page width (880 to 1300 px)
  const tip = h('div', {class: 'cmap-tip'}), wrap = h('div', {class: 'cmap-wrap',
    style: `max-width:${Math.round(Math.min(1300, Math.max(880, 880 + (nShown - 40) * 9)))}px`}, svg, tip);
  let touch = false, shown = null;
  svg.addEventListener('pointerdown', e => { touch = e.pointerType !== 'mouse'; });
  const place = e => {
    const r = wrap.getBoundingClientRect();
    let x = e.clientX - r.left + 14;
    if (x + tip.offsetWidth > r.width) x = Math.max(0, e.clientX - r.left - tip.offsetWidth - 14);
    tip.style.left = x + 'px'; tip.style.top = (e.clientY - r.top + 14) + 'px';
  };
  const show = (e, key, kids) => { tip.replaceChildren(...kids); tip.style.display = 'block'; place(e); shown = key; };
  const hide = () => { tip.style.display = 'none'; shown = null; };
  // hover (mouse) shows the evidence, click opens; on a touch screen the first tap shows, the second opens
  const hover = (g, key, kids, open) => {
    g.addEventListener('pointerenter', e => { if (e.pointerType === 'mouse') show(e, key, kids()); });
    g.addEventListener('pointermove', e => { if (e.pointerType === 'mouse') place(e); });
    g.addEventListener('pointerleave', e => { if (e.pointerType === 'mouse') hide(); });
    g.addEventListener('click', e => { if (touch && shown !== key) { show(e, key, kids()); return; } hide(); open(); });
  };
  const hint = () => h('div', {class: 'note'}, touch ? 'Tap again to open it.' : 'Click to open it.');
  // a faint wedge behind each sector; a click selects its category
  for (const s of M.sectors) {
    const r1 = 1.07 * R, r0 = 0.27 * R, big = s.a1 - s.a0 > Math.PI ? 1 : 0;
    const [x0, y0] = xy(s.a0, r1), [x1, y1] = xy(s.a1, r1), [x2, y2] = xy(s.a1, r0), [x3, y3] = xy(s.a0, r0);
    const w = svgEl('path', {d: `M${x0},${y0} A${r1},${r1} 0 ${big} 0 ${x1},${y1} L${x2},${y2} ` +
      `A${r0},${r0} 0 ${big} 1 ${x3},${y3}Z`, fill: s.color, 'fill-opacity': focus === s.category ? 0.22 : 0.09,
      style: 'cursor:pointer'}, svgEl('title', {}, `${s.label} (${s.n}): select`));
    w.addEventListener('click', () => pick(s.category));
    svg.append(w);
  }
  for (const it of M.items) {
    const [x, y] = xy(it.angle, R), [lx, ly] = xy(it.angle, 1.17 * R), [ex, ey] = xy(it.angle, 1.095 * R);
    const deg = it.angle * 180 / Math.PI, d = ((deg % 360) + 360) % 360, flip = d > 90 && d < 270;
    const on = !focus || focus === it.category;
    const g = svgEl('g', {class: 'it' + (on ? (focus ? ' hl' : '') : ' dim')});
    g.append(svgEl('line', {x1: 0, y1: 0, x2: x, y2: y, stroke: M.colors.line, 'stroke-opacity': 0.55,
      'stroke-width': it.line === 'clinical' ? 1.1 + 0.6 * it.n_clinical : 1.3,
      'stroke-dasharray': it.line === 'ontology' ? '7 4' : it.line === 'shared' ? '1.5 3.5' : null}));
    const multi = it.n_shared_regions >= 2;
    g.append(svgEl('circle', {cx: x, cy: y, r: (0.62 * Math.sqrt(45 + 32 * it.evidence)).toFixed(1),
      fill: color[it.category] || '#9AA1A8', stroke: multi ? M.colors.multi : '#fff', 'stroke-width': multi ? 2.6 : 1}));
    if (it.east_asian.length)
      g.append(svgEl('path', {d: `M${ex},${ey - 4.5}l4.5,4.5l-4.5,4.5l-4.5,-4.5z`, fill: M.colors.east_asian}));
    const lines = labelLines(it.label, 30, lineMax);
    g.append(svgEl('text', {transform: `translate(${lx},${ly}) rotate(${flip ? 180 - deg : -deg})`,
      'text-anchor': flip ? 'end' : 'start', 'dominant-baseline': 'central', style: `font-size:${fs.toFixed(1)}px`},
      lines.map((l, k) => svgEl('tspan', {x: 0, dy: k ? '1.1em' : `${-(lines.length - 1) * 0.55}em`}, l))));
    hover(g, it.id, () => {
      const a = (N.get(it.id) || {a: {}}).a, ev = [];
      if (it.n_clinical) ev.push(`clinical evidence, ${it.n_clinical} item${it.n_clinical > 1 ? 's' : ''}`);
      if (it.ontology) ev.push('the disease ontology');
      if (it.n_shared_regions) ev.push(`${it.n_shared_regions} ${S} GWAS region${it.n_shared_regions > 1 ? 's' : ''} ` +
        `shared${a.shared_regions ? ' (' + a.shared_regions + ')' : ''}`);
      return [h('b', {}, it.name), h('div', {class: 'k'}, sectorLabel(M, it.category)),
              h('div', {}, h('span', {class: 'k'}, 'Evidence: '), ev.join('; ')),
              it.east_asian.length ? h('div', {}, h('span', {class: 'k'}, 'East Asian: '),
                it.east_asian.map(k => M.east_asian_kinds[k] || k).join('; ')) : null, hint()];
    }, () => { location.hash = href(it.id); });
    svg.append(g);
  }
  if (M.single) {                                            // the traits sharing a single region: one node
    const sg = M.single, [x, y] = xy(sg.angle, R), [lx, ly] = xy(sg.angle, 1.17 * R);
    const deg = sg.angle * 180 / Math.PI, d = ((deg % 360) + 360) % 360, flip = d > 90 && d < 270;
    const g = svgEl('g', {class: 'it' + (!focus || focus === 'single' ? (focus ? ' hl' : '') : ' dim')});
    g.append(svgEl('line', {x1: 0, y1: 0, x2: x, y2: y, stroke: M.colors.line, 'stroke-opacity': 0.55,
                            'stroke-width': 1.3, 'stroke-dasharray': '1.5 3.5'}),
             svgEl('circle', {cx: x, cy: y, r: 17, fill: M.colors.single, stroke: '#fff'}),
             svgEl('text', {x, y, 'text-anchor': 'middle', 'dominant-baseline': 'central', style: 'font-size:12px'},
                   String(sg.n)),
             svgEl('text', {transform: `translate(${lx},${ly}) rotate(${flip ? 180 - deg : -deg})`,
               'text-anchor': flip ? 'end' : 'start', 'dominant-baseline': 'central',
               style: 'font-style:italic;fill:#57606a'}, `One ${S} region only`));
    hover(g, '__single', () => [h('b', {}, `${sg.n} trait${sg.n === 1 ? ' shares' : 's share'} a single ${S} GWAS ` +
      'region'),
      h('div', {}, `${sg.n === 1 ? 'The' : sg.n_top + ' of them the'} ${sg.top} region: weak evidence of a relation.`),
      h('div', {class: 'note'}, touch ? 'Tap again to list them.' : 'Click to list them below.')],
      () => pick('single'));
    svg.append(g);
  }
  const core = svgEl('g', {class: 'core'}, svgEl('circle', {r: 44, fill: M.colors.core}),
    svgEl('text', {'text-anchor': 'middle', 'dominant-baseline': 'central',
                   style: `fill:#fff;font-weight:700;font-size:${S.length <= 6 ? 22 : Math.max(11, 132 / S.length)}px`}, S));
  hover(core, '__core', () => [h('b', {}, D.disease.Name), hint()], () => { location.hash = href(CORE); });
  svg.append(core);
  // the categories (a button selects one), and the key
  const cats = h('div', {class: 'cmap-cats'}, M.sectors.map(s => h('button', {class: 'btn' + (focus === s.category ?
      ' on' : ''), onclick: () => pick(s.category)}, h('span', {class: 'dot', style: 'background:' + s.color}),
      `${s.label} (${s.n})`)),
    M.single ? h('button', {class: 'btn' + (focus === 'single' ? ' on' : ''), onclick: () => pick('single')},
      h('span', {class: 'dot', style: 'background:' + M.colors.single}), `One ${S} region only (${M.single.n})`) : null);
  const sample = (at, extra) => svgEl('svg', {width: 30, height: 12, viewBox: '0 0 30 12'},
    svgEl('line', {x1: 1, y1: 6, x2: 29, y2: 6, stroke: M.colors.line, ...at}), extra || null);
  // the key names only what the map draws
  const has = f => M.items.some(f);
  const key = h('div', {class: 'cmap-key'},
    has(i => i.line === 'clinical') ? h('span', {}, sample({'stroke-width': 3}), 'Clinical evidence (width: items)') : null,
    has(i => i.line === 'ontology') ? h('span', {}, sample({'stroke-width': 1.5, 'stroke-dasharray': '7 4'}),
      'Disease ontology only') : null,
    has(i => i.line === 'shared') || M.single ? h('span', {}, sample({'stroke-width': 1.5,
      'stroke-dasharray': '1.5 3.5'}), 'Shared GWAS signal only') : null,
    has(i => i.n_shared_regions >= 2) ? h('span', {}, svgEl('svg', {width: 14, height: 14}, svgEl('circle', {cx: 7,
      cy: 7, r: 5, fill: '#fff', stroke: M.colors.multi, 'stroke-width': 2.4})), `Shares ≥ 2 ${S} GWAS regions`) : null,
    has(i => i.east_asian.length) ? h('span', {}, svgEl('svg', {width: 12, height: 12}, svgEl('path',
      {d: 'M6,1l5,5l-5,5l-5,-5z', fill: M.colors.east_asian})), 'East Asian evidence of the relation') : null,
    h('span', {}, 'Node size: amount of evidence'));
  // what the map shows and leaves out, of every related condition; the rules are on the build page
  const sg = M.single, ns = (M.not_shown || {}).shared_signal_without_region || 0;
  const coverage = h('div', {class: 'note', style: 'text-align:center'}, `Shown: ${M.items.length} conditions` +
    (sg ? ` and one node for the ${sg.n} trait${sg.n === 1 ? ' that shares' : 's that share'} a single ${S} region ` +
          `only (${sg.n === 1 ? '' : sg.n_top + ' of them '}` +
          `${sg.top}; weak evidence)` : '') +
    (ns ? `; not shown: ${ns} traits of the shared-signal category with no shared region recorded` : '') +
    `, of ${related().length} related conditions. How a condition enters the graph and how it is classified: `,
    h('a', {href: '#/build'}, 'How it is built'), '.');
  return h('div', {}, cats, wrap, key, coverage);
}

function conditionsPage(selected) {
  const M = D.condition_map && (D.condition_map.items.length || D.condition_map.single) ? D.condition_map : null;
  // the selection is part of the address (#/conditions/<category>), so that it can be shared
  let focus = M && selected && (selected === 'single' ? M.single : M.sectors.some(s => s.category === selected)) ?
    selected : null;
  const mapBox = h('div'), listBox = h('div');
  const members = f => new Set(f === 'single' ? M.single.ids : M.items.filter(i => i.category === f).map(i => i.id));
  const pick = f => {
    focus = focus === f ? null : f;
    try { history.replaceState(null, '', '#/conditions' + (focus ? '/' + encodeURIComponent(focus) : '')); }
    catch (e) { /* an address that cannot change: the selection still applies */ }
    draw();
  };
  const draw = () => {
    if (M) mapBox.replaceChildren(conditionMap(M, focus, pick));
    const keep = focus ? members(focus) : null, rows = related().filter(c => !keep || keep.has(c.id));
    const name = focus === 'single' ? `Traits sharing a single ${S} GWAS region` : focus ? sectorLabel(M, focus) : '';
    listBox.replaceChildren(h('h3', {}, focus ? `${name} (${rows.length})` : `Every related condition (${rows.length})`,
      focus ? h('button', {class: 'btn', style: 'margin-left:10px', onclick: () => pick(focus)},
                'Show every condition') : null),
      table([{label: 'Condition', get: c => c.label, cell: c => nodeLink(c.id)},
             {label: 'Category', get: c => String(catRank(c.a.category)).padStart(3, '0') + ' ' + c.a.category,
              cell: c => c.a.category},
             {label: 'Evidence types', get: c => c.a.sources},
             {label: 'Clinical items', num: true, get: c => c.a.n_clinical},
             {label: `${S} regions shared`, num: true, get: c => c.a.n_shared_regions,
              cell: c => c.a.shared_regions || ''},
             {label: 'GWAS studies', num: true, get: c => c.a.n_gwas_studies},
             {label: 'East Asian GWAS', num: true, get: c => c.a.n_eas_studies},
             {label: 'Curated genes', num: true, get: c => c.a.n_mendelian_genes}], rows, {sort: 1, limit: 200}));
  };
  draw();
  return [h('h2', {}, 'Related conditions'), h('p', {class: 'lead'}, `Every condition related to ${S}, with the ` +
    'kinds of evidence for the relation. In the map, hover over a condition for its evidence and click it to open its ' +
    'page; select a category (a button or its sector) to list its conditions in the table below.'),
    mapBox, listBox];
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
      rows, {sort: 0, limit: 200, tools: h('span', {}, ['all', 'Japanese', 'East Asian', 'Not East Asian'].map(x =>
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
const here = (...kids) => h('div', {class: 'here'}, h('b', {}, S + ': '), ...kids);
const count = (xs, f) => xs.reduce((m, x) => { const k = f(x); m[k] = (m[k] || 0) + 1; return m; }, {});
const tallyText = (m, order) => (order || Object.keys(m)).filter(k => m[k]).map(k => `${k} ${m[k]}`).join('; ');

function buildPage() {
  const P = D.params, L = D.llm, m = D.meta, cc = D.curation, R = L.relations, X = L.extract || {};
  const num0 = v => (v ?? 0).toLocaleString(), pv = (k, d) => P[k] ?? d;
  const n = t => D.edges.filter(e => e[2] === t).length;
  const regions = byType('region'), conds = related(), genes = byType('gene').filter(g => lvl(g) !== null);
  const studies = byType('publication').filter(p => p.a.design_class);
  const nPub = regions.filter(published).length;
  const kinds = {clinical: 0, ontology: 0, shared: 0};
  for (const c of conds) { const s = String(c.a.sources || '').split(' ');
    if (s.includes('clinical')) kinds.clinical++; if (s.includes('ontology')) kinds.ontology++;
    if (s.includes('colocalization') || s.includes('lead-variant')) kinds.shared++; }
  const model = r => (L.models.find(x => x.role === r) || {}).model || '';
  const out = [h('h2', {}, 'How this knowledge graph is built'), h('p', {class: 'lead'}, 'Databases give the facts, ' +
    'two language models read the literature, and rules decide what the graph counts. The rules are the same for ' +
    `every disease; the lines marked ${S} give this snapshot.`)];

  // ── at a glance: sources → methods → graph ──
  const item = (k, title, text) => h('div', {class: 'gl-item'}, h('div', {class: 'gl-t'}, tag(k), ' ', title),
                                     h('div', {class: 'gl-d'}, text));
  const tcount = t => byType(t).length;
  const nodeRow = (t, label, v) => h('div', {class: 'gl-n'}, dot(t), h('span', {}, label),
                                      h('b', {}, (v ?? tcount(t)).toLocaleString()));
  out.push(h('div', {class: 'glance'},
    h('div', {class: 'gl-col'}, h('div', {class: 'gl-h'}, '1 · Sources'),
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
    h('div', {class: 'gl-col'}, h('div', {class: 'gl-h'}, '2 · Methods'),
      item('rule', 'Rules', 'search and pre-screen · curation · thresholds · evidence ladder · GWAS regions · ' +
        'shared signals · categories'),
      item('model', 'Two readers, one adjudicator', `${model('Reader A')} and ${model('Reader B')} read each record ` +
        `apart; ${model('Adjudicator')} decides where they differ`),
      item('rule', 'Source-text check', 'every gene, variant, number and quote must be found in its text'),
      h('div', {class: 'gl-foot'}, `${num0(X.read_by_model)} records read` + (R ? `, ${num0(R.read ?? R.records)} ` +
        'for relations' : '') + ((X.not_read || (R && R.not_read)) ? ` · ${num0((X.not_read || 0) + ((R && R.not_read) ||
        0))} not read (no model called, none cached)` : '') +
        (agreeShare() ? ` · ${agreeShare()} of decisive items agreed by both readers` : ''))),
    h('div', {class: 'gl-arrow'}, '→'),
    h('div', {class: 'gl-col'}, h('div', {class: 'gl-h'}, '3 · Knowledge graph'),
      nodeRow('condition', 'Conditions', conds.length + 1), nodeRow('gene', 'Genes'),
      nodeRow('region', `${S} GWAS regions`), nodeRow('publication', 'Publications'),
      nodeRow('gwas_study', 'GWAS studies'), nodeRow('theme', 'Gene functions'), nodeRow('pathway', 'Reactome pathways'),
      h('div', {class: 'gl-foot'}, `${num0(m.n_edges)} edges of ${Object.keys(D.schema.edges).length} types, each ` +
        'with its source · viewer, workbook, knowledge base, SQLite'))));

  // ── the build, step by step ──
  const step = (no, title, kinds, text, mine) => h('div', {class: 'stp'}, h('div', {class: 'stp-no'}, no),
    h('div', {class: 'stp-b'}, h('div', {class: 'stp-h'}, h('b', {}, title), ' ', kinds.map(tag)),
      h('div', {}, text), mine ? here(mine) : null));
  const sm = L.summaries, ag = X.agreement || {};
  out.push(h('h3', {}, 'The build, step by step'), h('div', {class: 'stps'},
    step(1, 'The disease', ['profile', 'ot'], 'Identified by its Open Targets ID. The profile written from the ID gives ' +
      'its name, short name and PubMed terms, each with the rule that chose it; optional curated files replace a ' +
      'default. Other Open Targets entries can lend their evidence (a grouping name joins diseases Open Targets ' +
      'holds as entries of their own).', `${D.disease.Name} (${D.disease.id})` + (L.evidence_also.length ?
      '; evidence also from ' + L.evidence_also.map(x => `${x.name} (${x.id})`).join(', ') + ', whose targets, ' +
      'text-mined literature, GWAS and curated genes are merged with its own.' : '.')),
    step(2, 'Sources', ['ot', 'pubmed', 'epmc'], 'Fetched for the snapshot date and saved with it, so that the ' +
      'snapshot can be rebuilt.', `Open Targets ${m.release}; PubMed to ${L.fetch_date || m.snapshot}.`),
    step(3, 'Related conditions', L.curated_relations ? ['profile', 'ot', 'rule'] : ['model', 'ot', 'rule'],
      'A condition enters by ' +
      'clinical evidence that found an association (a quoted abstract: from a curated list, else read by the models; ' +
      'a study that found none is kept as contrary evidence beside a relation, never as one), by the disease ' +
      'ontology, or by a GWAS signal shared with the disease; its category follows its strongest evidence.',
      `${conds.length} conditions: ${kinds.clinical} by clinical evidence, ${kinds.ontology} by the ontology, ` +
      `${kinds.shared} by a shared signal (a condition can have several); relations ` +
      (L.curated_relations ? 'from the curated list.' : 'read from PubMed' + (R && R.not_read ? ` (${num0(R.not_read)}` +
        ' of its records not read: no language model was called, and no answer was cached).' : '.')) +
      (L.not_related.length ? ` Not related: ${L.not_related.length} condition${L.not_related.length > 1 ? 's' : ''}` +
        ` whose clinical evidence found no association (${L.not_related.join('; ')}).` : '')),
    step(4, 'Genetic studies of the disease', ['pubmed', 'epmc', 'model', 'rule'], 'A rule pre-screen, then two ' +
      'readers and an adjudicator decide inclusion and extract the findings, sample sizes and methods; every item is ' +
      'checked against its text.', `${num0(X.records)} records; ${num0(X.prescreened_out)} left out by the ` +
      `pre-screen; ${num0(X.read_by_model)} read` + (X.not_read ? `, ${num0(X.not_read)} not read (no language ` +
      'model was called, and no answer was cached)' : '') + `; ${num0(X.included)} included by the readers, ` +
      `${num0(X.findings)} findings.`),
    step(5, 'Curation', ['pubmed', 'rule'], 'Retracted articles and conference material leave the graph; preprints ' +
      'and errata are flagged; each study gets its place and design; a variant gets one name; a finding is marked by ' +
      'what it is about; text-mined links to the wrong gene are removed.',
      `${num0(cc.counts.publications)} publications checked; ${cc.n_excluded} left out, ${cc.flagged.length} ` +
      `flagged; ${num0(cc.counts.studies)} genetic studies kept` + (cc.counts.studies_excluded ?
      ` (${cc.counts.studies_excluded} left out)` : '') + `; ${num0(cc.counts.textmining_dropped)} text-mined links ` +
      'removed.'),
    step(6, 'Evidence and tables', ['ot', 'rule'], 'The GWAS regions of the disease, the evidence level and function ' +
      'of each gene, and the genes it shares with the related conditions.', `${regions.length} regions` +
      (regions.length ? ` (${nPub === regions.length ? 'published' : !nPub ? 'from Open Targets' :
        `${nPub} published, ${regions.length - nPub} from Open Targets`})` : '') + `; ${genes.length} genes on the ` +
      `ladder; ${byType('theme').length} functions.`),
    step(7, 'Graph', ['rule'], 'Typed nodes and edges; every edge keeps its evidence and its source.',
      `${num0(m.n_nodes)} nodes, ${num0(m.n_edges)} edges.`),
    step(8, 'Summaries', ['model', 'rule'], 'One model writes, the other checks every claim, the adjudicator rules; ' +
      'a summary with an unsupported claim is withheld with its reason.',
      `${sm.shown} shown, ${sm.withheld} withheld` + (sm.not_generated ? `, ${sm.not_generated} not generated` : '') +
      '.')));

  // ── in detail ──
  const det = (title, sub, ...kids) => h('details', {class: 'bd'}, h('summary', {}, h('b', {}, title),
    h('span', {class: 'bd-sub'}, sub)), h('div', {class: 'bd-body'}, ...kids));
  out.push(h('h3', {}, 'In detail'));
  // the sources
  const srcRow = (k, gives, how) => h('tr', {}, h('td', {}, tag(k)), h('td', {}, gives), h('td', {}, how));
  out.push(det('Sources', 'what each database gives, and how it is read', h('table', {class: 'plain'},
    h('thead', {}, h('tr', {}, ['Source', 'What it gives', 'How it is read'].map(x => h('th', {class: 'nosort'}, x)))),
    h('tbody', {},
      srcRow('ot', 'The disease and its ontology (names, synonyms, parents, children); GWAS studies and credible sets ' +
        'with locus-to-gene (L2G) scores and colocalizations; curated gene–disease evidence (ClinGen, Gene2Phenotype, ' +
        'PanelApp, Orphanet, ClinVar, UniProt, gene burden); Reactome pathways; text-mined targets.',
        'By rules only; values are kept as given, with their identifiers.'),
      srcRow('pubmed', 'The genetic studies of the disease (its name and synonyms with genetic terms); the clinical ' +
        'relations (guidelines, reviews, cohort and case-control studies with relation words); the record of every ' +
        'publication the graph uses (type, retraction, affiliation).',
        'Abstracts by the language models; types, retractions and affiliations by rules.'),
      srcRow('epmc', 'The methods text of open-access articles; the sentences behind text-mined gene mentions.',
        'Methods by the language models; mentions by rules.'),
      srcRow('profile', 'The Open Targets ID (required); optional curated files: clinical relations, published GWAS ' +
        'regions, gene functions, curation entries.', 'By rules.')))));
  // the associations
  const EDGE = {
    related_to: [['profile', 'model', 'ot'], 'Clinical: a PubMed abstract quoted verbatim (curated list, or read by ' +
      `the models). Ontology: an Open Targets neighbor (a parent with more than ${pv('MaxOntologyChildren', 30)} ` +
      'children left out). Shared signal: see shares_region.'],
    supports: [['pubmed'], 'The publication a clinical relation quotes.'],
    has_genetic_evidence: [['ot', 'model', 'rule'], `G: the top L2G gene (≥ ${pv('L2gMin', 0.5)}) of a GWAS credible ` +
      'set (for the disease: the genes of its GWAS regions; G*: significant only in a secondary analysis). C: curated ' +
      `evidence (source score ≥ ${pv('MendelianMin', 0.5)}; ClinVar alone, a pathogenic record of two stars or ` +
      'more). L: a positive finding of a genetic study of the disease.'],
    reports: [['model', 'rule'], 'A finding the readers took from an abstract, kept when its gene, variant and ' +
      'numbers are in the text.'],
    studies: [['model', 'pubmed', 'rule'], 'A record the readers include as a genetic study of the disease; design ' +
      'class and place by rule.'],
    region_of: [['profile', 'ot', 'rule'], 'The published regions of the profile, and the credible sets of the ' +
      `disease's GWAS in Open Targets that no published region holds, with lead P < ${pv('GenomeWideP', 5e-8)} and a ` +
      `credible lead: at least ${pv('RegionMinMac', 25)} expected copies of its minor allele in the smaller group of ` +
      `its study and, below a frequency of ${pv('RegionMinMaf', 0.01)}, P < ${pv('RareVariantP', 5e-9)}; or replicated ` +
      `by another GWAS of the disease. Merged within ${pv('RegionWindowKb', 500)} kb.`],
    names_gene: [['profile', 'ot'], `The gene the publication names, else the top L2G gene (≥ ${pv('L2gMin', 0.5)}), ` +
      'else the position.'],
    shares_region: [['ot', 'rule'], `Colocalization with the region (CLPP ≥ ${pv('ClppMin', 0.01)} or H4 ≥ ` +
      `${pv('H4Min', 0.8)}), or, without a credible set of the disease, its lead variant in a credible set of the ` +
      `condition (PIP ≥ ${pv('ClppMin', 0.01)}).`],
    in_theme: [['profile', 'ot', 'rule'], 'The gene functions of the profile, else the Reactome top-level pathway ' +
      'that holds most genes of the ladder.'],
    in_pathway: [['ot'], 'The Reactome top-level pathways of a text-mined target.'],
    mentions: [['epmc', 'ot', 'rule'], 'Europe PMC text mining as Open Targets gives it; a word mapped to the wrong ' +
      'gene, or naming the disease or a neighbor, is removed.'],
    gwas_of: [['ot'], 'The GWAS studies of the disease and of each related condition; a condition with more than ' +
      `${pv('MaxGwasStudies', 100)} GWAS studies keeps its shared signals only, and no study of it is listed.`]};
  out.push(det('Associations', `the ${Object.keys(D.schema.edges).length} edge types: where each comes from, by which ` +
    'rule', h('table', {class: 'plain'},
    h('thead', {}, h('tr', {}, ['Edge', 'Meaning', 'Comes from', 'Rule (this run)', S].map(x => h('th', {class: 'nosort'},
      x)))),
    h('tbody', {}, Object.keys(D.schema.edges).map(k => h('tr', {}, h('td', {}, h('code', {}, k)),
      h('td', {}, D.schema.edges[k]), h('td', {}, (EDGE[k] || [[]])[0].map(tag)), h('td', {}, (EDGE[k] || [, ''])[1]),
      h('td', {class: 'num'}, n(k).toLocaleString())))))));
  // categories and levels
  const ladder = Object.entries(D.ladder).map(([k, [code, label]]) => `${k}${code ? ' ' + code : ''} ${label}`);
  const card = (title, rule, mine) => h('div', {class: 'llm-step'}, h('h4', {}, title), h('p', {}, rule), here(mine));
  out.push(det('Categories and levels', 'how a condition, a gene and a study are classified', h('div',
    {class: 'llm-steps'},
    card('Category of a related condition', 'From its strongest evidence: a clinical relation gives its category ' +
      '(curated list, or as the readers gave it and the adjudicator confirmed); ontology only joins the disease group; ' +
      'a shared GWAS signal only joins “Shared GWAS signal only”.',
      tallyText(count(conds, x => SHOWN_CAT(x.a.category)), CATS.map(SHOWN_CAT)) + '.'),
    card('Evidence level of a gene', 'The strongest level it reaches: ' + ladder.join('; ') + '. Levels 3–6 count ' +
      'only findings about the disease itself.',
      tallyText(count(genes, g => 'level ' + lvl(g)), Object.keys(D.ladder).map(k => 'level ' + k)) + '.'),
    card('Function of a gene', 'The functions of the profile; otherwise the Reactome top-level pathway holding most ' +
      `genes of the ladder (“Disease” excluded), the rest “${D.other_genes}”.`,
      `${byType('theme').length} functions (${L.themes || 'curated'}).`),
    card('Design and place of a study', 'The design class from the design the readers quote; the place from the ' +
      'first author\'s institution (PubMed, else Europe PMC); East Asian when done there or stated by the abstract.',
      tallyText(count(studies, p => p.a.design_class)) + '.'))));
  // the language models, every rule
  out.push(det('Language models', 'the three models, one reading from record to graph, the checks, agreement, ' +
    'every disputed item', ...modelSections()));
  out.push(det('Every rule', 'the choices registry: each rule with its value in this run and its basis',
    ...D.choices.map(g => h('details', {class: 'cmap-all'}, h('summary', {}, `${g.name} (${g.items.length})`),
      table([{label: 'Rule', get: r => r.item}, {label: 'Value', get: r => r.value}, {label: 'Basis', get: r => r.basis},
             {label: 'References', cell: r => h('span', {}, r.refs.map(([tt, u], k) =>
               h('span', {}, k ? '; ' : '', u ? ext(u, tt) : tt)))}], g.items, {filter: false, limit: 200})))));
  return out;
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
        {label: 'Source', get: r => r.pm, cell: r => pmid(r.pm)}, {label: 'Pass', get: r => r.ps},
        {label: 'Item', get: r => r.item}, {label: 'Reader A', get: r => r.a}, {label: 'Reader B', get: r => r.b},
        {label: 'Final', get: r => r.fin}, {label: 'Decision', get: r => r.dec, cell: r => h('span', {class: 'quote'}, r.dec)}],
        rows, {limit: 100, tools: h('span', {}, ['all'].concat(passes).map(x => h('button', {class: 'btn' + (pass === x ?
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
             {label: 'Why', get: r => r[6]}], flagged, {limit: 100}));
  }
  return out;
}

// ── routing, search ──────────────────────────────────────────────────────
const TABS = [['overview', 'Overview'], ['genes', 'Genes'], ['conditions', 'Related conditions'],
              ['regions', 'GWAS regions'], ['studies', 'Genetic studies'], ['curation', 'Curation'],
              ['build', 'How it is built']];
const PAGES = {overview, genes: genesPage, conditions: conditionsPage, regions: regionsPage, studies: studiesPage,
               curation: curationPage, build: buildPage};
function render() {
  const hash = decodeURIComponent(location.hash.slice(2) || 'overview');
  const main = document.getElementById('main');
  let tab = hash.split('/')[0];
  if (hash.startsWith('n/')) { main.replaceChildren(nodePage(hash.slice(2))); tab = ''; }
  else main.replaceChildren(...[].concat((PAGES[tab] || overview)(hash.split('/').slice(1).join('/'))));
  document.querySelectorAll('nav a').forEach(a => a.classList.toggle('on', a.dataset.tab === tab));
  window.scrollTo(0, 0);
}
document.getElementById('nav').append(...TABS.map(([k, l]) => h('a', {href: '#/' + k, 'data-tab': k}, l)));
// the full name and the ID beside the short name, which can name other things too
document.getElementById('meta').textContent = (D.disease.name.toLowerCase() === S.toLowerCase() ? '' :
  `${D.disease.Name} · `) + `${D.disease.id} · Snapshot ${D.meta.snapshot} · Open Targets ${D.meta.release} · ` +
  `PubMed checked ${D.meta.checked_at} · ${D.meta.n_nodes.toLocaleString()} nodes, ` +
  `${D.meta.n_edges.toLocaleString()} edges`;
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
// printing shows every rule of the condition map
window.addEventListener('beforeprint', () => document.querySelectorAll('details.cmap-all, details.bd')
  .forEach(d => { d.open = true; }));
render();
</script>
</body>
</html>
'''

if __name__ == '__main__':
    main()
