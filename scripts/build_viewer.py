#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The viewer of one snapshot: a single HTML file that opens in any
#           browser, offline, and shows the whole knowledge graph. No library,
#           no network; the data are embedded.
#             Overview     counts, the evidence ladder (level x function), the
#                          CTEPH GWAS regions x the conditions sharing them
#             Genes        the CTEPH genes by level (or every gene), filterable
#             Conditions   the related conditions by category
#             Regions      the nine published CTEPH GWAS regions
#             Studies      the genetic studies of CTEPH, by place and design
#             Curation     what was left out or corrected, and why
#           Every node has its own page (#/n/<id>): its facts, a diagram of its
#           neighbors, and every edge with its evidence and source link
#           (PubMed, Open Targets, dbSNP), so that each statement can be traced
#           to where it comes from.
# Inputs  : 05.graph (nodes_*.tsv, edges_*.tsv, schema.json), 03.curated,
#           04.tables/literature.tsv (one sentence per text-mined link),
#           06.summaries (grounded summaries only)
# Outputs : cteph_knowledge_graph.html
# ---------------------------------------------------------------------------
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import LADDER, log, read_tsv  # noqa: E402


def compact(row, skip=()):
    return {k: v for k, v in row.items() if k not in skip and v not in ('', None)}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--graph', required=True, help='05.graph')
    ap.add_argument('--curated', required=True, help='03.curated')
    ap.add_argument('--tables', required=True, help='04.tables')
    ap.add_argument('--summaries', required=True, help='summaries.tsv')
    ap.add_argument('--out', default='cteph_knowledge_graph.html')
    a = ap.parse_args()
    G, Q, T = Path(a.graph), Path(a.curated), Path(a.tables)
    schema = json.loads((G / 'schema.json').read_text())
    lit = read_tsv(T / 'literature.tsv')
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
    pubs = read_tsv(Q / 'publications.tsv')
    cur = json.loads((Q / 'curation.json').read_text())
    lg = read_tsv(Q / 'curation_log.tsv')
    summ = read_tsv(a.summaries)
    data = {
        'meta': {**schema['graph'], 'checked_at': cur['checked_at'][:10]},
        'schema': {'nodes': {k: v['meaning'] for k, v in schema['nodes'].items()},
                   'edges': {k: v['meaning'] for k, v in schema['edges'].items()}},
        'ladder': {str(k): [c, lab] for k, (c, lab) in LADDER.items()},
        'nodes': nodes, 'edges': edges,
        'curation': {
            'counts': {k: cur[k] for k in ('publications', 'studies', 'studies_excluded', 'textmining_dropped',
                                           'findings_by_subject')},
            'n_excluded': int((pubs.status == 'excluded').sum()),
            'flagged': [compact(r) for r in pubs[pubs.status == 'flagged'][
                ['pmid', 'first_author', 'year', 'title', 'reason', 'roles']].to_dict('records')],
            'log': lg[~lg.after.str.startswith('publication left out') & (lg.step != 'publications')][
                ['step', 'item', 'rule', 'before', 'after']].values.tolist()},
        'summaries': {r.title: r.summary for r in summ[(summ.kind == 'mechanism') & (summ.grounded == 'True')]
                      .itertuples()},
    }
    blob = json.dumps(data, ensure_ascii=False, separators=(',', ':')).replace('<', '\\u003c')
    Path(a.out).write_text(TEMPLATE.replace('__DATA__', blob).replace('__TITLE__', 'CTEPH knowledge graph'))
    log(f'viewer: {len(nodes)} nodes, {len(edges)} edges, {Path(a.out).stat().st_size / 1e6:.1f} MB')


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
const CORE = D.meta.core;
const TYPE_NAME = {condition: 'Condition', gene: 'Gene', region: 'CTEPH GWAS region', publication: 'Publication',
                   gwas_study: 'GWAS study', theme: 'Function', pathway: 'Reactome pathway'};
const COLOR = t => getComputedStyle(document.documentElement).getPropertyValue('--t-' + t).trim() || '#888';
const CATS = ['Precursor / spectrum', 'Risk factor', 'Comorbidity', 'Intermediate trait', 'PH group',
              'Shared genetic signal'];

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
const lvl = n => n.a.cteph_level !== undefined ? Number(n.a.cteph_level) : null;

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
      h('td', {class: c.num ? 'num' : null}, c.cell ? c.cell(r) : fmt(c.get(r)))))));
    count.textContent = rs.length === rows.length ? `${rows.length} rows` : `${rs.length} of ${rows.length} rows`;
    more.replaceChildren(rs.length > limit ? h('button', {class: 'btn', onclick: () => { limit += 500; draw(); }},
      `Show more (${rs.length - limit} hidden)`) : '');
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
const SECTION = {
  'related_to:out': 'Why it relates to CTEPH', 'related_to:in': 'Related conditions, one row per piece of evidence',
  'supports:in': 'Publications that give its clinical evidence', 'supports:out': 'Clinical evidence it gives',
  'has_genetic_evidence:out': 'Genetic evidence in conditions', 'has_genetic_evidence:in': 'Genes with genetic evidence',
  'reports:in': 'Findings of the genetic studies of CTEPH', 'reports:out': 'Findings it reports',
  'studies:out': 'Genetic study of', 'studies:in': 'Genetic studies of CTEPH',
  'region_of:out': 'Region of', 'region_of:in': 'Published CTEPH GWAS regions', 'names_gene:out': 'Gene named',
  'names_gene:in': 'CTEPH GWAS region', 'shares_region:out': 'CTEPH GWAS regions it shares',
  'shares_region:in': 'Conditions sharing this region', 'in_theme:out': 'Function (gene theme)',
  'in_theme:in': 'Genes of this function', 'in_pathway:out': 'Reactome pathways',
  'in_pathway:in': 'Text-mined CTEPH targets in this pathway', 'mentions:in': 'Text-mined CTEPH publications',
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
               A('east_asian', 'East Asian')],
  supports: [A('relation', 'Relation'), A('detail', 'Statement')],
  has_genetic_evidence: [A('code', 'Code'), A('kind', 'Kind'), A('n', 'Records', 'num'), A('details', 'Details')],
  reports: [A('result', 'Result'), A('about', 'About'), A('variant_name', 'Variant'), A('rsid', 'rsID', 'rsid'),
            A('effect', 'Effect'), A('p_value', 'P'), A('gene_set', 'Gene set'), A('quote', 'Source sentence', 'quote')],
  studies: [A('design_class', 'Design class'), A('n_cases', 'Cases', 'num'), A('n_controls', 'Controls', 'num'),
            A('place', 'Place')],
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
  gene: [A('cteph_level', 'CTEPH level', 'num')], condition: [A('category', 'Category')],
  gwas_study: [A('n_samples', 'Samples', 'num'), A('ancestry', 'Ancestry')]};

function cellOf(kind, v) {
  if (v === undefined || v === '') return '';
  if (kind === 'pmid') return pmid(v);
  if (kind === 'rsid') return rsid(v);
  if (kind === 'quote') return h('span', {class: 'quote'}, v);
  return fmt(v);
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
  return h('section', {}, h('h3', {}, `${SECTION[key] || key} (${list.length})`),
           h('div', {class: 'note'}, D.schema.edges[type] || ''), table(cols, list, {limit: 100}));
}

// ── facts of a node ───────────────────────────────────────────────────────
const FACTS = {
  gene: [['cteph_level', 'CTEPH evidence level', v => `${v} — ${(D.ladder[v] || [])[1] || ''}` +
          ((D.ladder[v] || [])[0] ? ` (${D.ladder[v][0]})` : '')], ['theme', 'Function'],
         ['n_related_conditions', 'Related conditions with genetic evidence'],
         ['n_text_mined_publications', 'Text-mined CTEPH publications'],
         ['cteph_ot_score', 'Open Targets association score with CTEPH'], ['pathways', 'Reactome pathways']],
  condition: [['category', 'Category'], ['sources', 'Evidence types'], ['n_clinical', 'Clinical evidence items'],
              ['shared_regions', 'CTEPH regions shared'], ['n_gwas_studies', 'GWAS studies'],
              ['n_eas_studies', 'East Asian GWAS studies'], ['n_credible_sets', 'Credible sets'],
              ['n_mendelian_genes', 'Curated genes'], ['genetics', 'Genetic data fetched'], ['ot_ids', 'Ontology IDs']],
  region: [['rsid', 'Lead variant', rsid], ['chrom', 'Chromosome'], ['pos_grch37', 'Position (GRCh37)'],
           ['tier', 'Tier'], ['code', 'Code'], ['p', 'P'], ['p_pe', 'P in PE'], ['v', 'V (with PE)'], ['or', 'OR'],
           ['maf', 'MAF'], ['open_targets', 'In Open Targets']],
  publication: [['title', 'Title'], ['first_author', 'First author'], ['year', 'Year'], ['journal', 'Journal'],
                ['roles', 'Role in the graph'], ['place', "Place (first author's institution)"],
                ['design_class', 'Design class'], ['n_cases', 'Cases'], ['n_controls', 'Controls'],
                ['status', 'Flag']],
  gwas_study: [['trait', 'Trait'], ['first_author', 'First author'], ['year', 'Year'], ['n_samples', 'Samples'],
               ['n_credible_sets', 'Credible sets'], ['ancestry', 'Ancestry (LD reference)'], ['eas', 'East Asian'],
               ['pmid', 'Publication', pmid]],
  theme: [['description', 'Description']], pathway: []};

function nodePage(id) {
  const n = N.get(id);
  if (!n) return h('div', {}, h('h2', {}, 'Not found'), h('p', {}, id));
  const facts = h('div', {class: 'facts'});
  for (const [k, label, f] of FACTS[n.type] || []) {
    const v = n.a[k]; if (v === undefined || v === '') continue;
    facts.append(h('div', {class: 'k'}, label), h('div', {}, f ? f(v) : fmt(v)));
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
const themes = () => byType('theme').map(t => t.label).concat(['Other genes']);
const short = label => (byType('theme').find(t => t.label === label) || {a: {}}).a.short || label;

function overview() {
  const m = D.meta, c = D.curation.counts, genes = byType('gene').filter(n => lvl(n) !== null);
  const pubs = byType('publication').filter(p => p.a.design_class);
  const card = (n, t, link) => h('a', {href: link || null, class: 'card', style: 'text-decoration:none;color:inherit'},
                                  h('div', {class: 'n'}, n), h('div', {class: 't'}, t));
  const out = [h('h2', {}, 'Overview'),
    h('p', {class: 'lead'}, 'The knowledge graph of chronic thromboembolic pulmonary hypertension (CTEPH): the genes ' +
      'with genetic evidence for CTEPH and how strong that evidence is, the genetic studies it comes from, the ' +
      'conditions related to CTEPH and the genetics they share with it. Every statement links to its source.'),
    h('div', {class: 'cards'},
      card(genes.length, 'CTEPH genes on the evidence ladder', '#/genes'),
      card(pubs.length, `genetic studies of CTEPH (${pubs.filter(p => p.a.japanese === 'True').length} in Japan)`,
           '#/studies'),
      card(byType('condition').length - 1, 'related conditions', '#/conditions'),
      card(byType('region').length, 'published CTEPH GWAS regions', '#/regions'),
      card(`${c.publications.toLocaleString()}`, 'publications checked in PubMed', '#/curation'),
      card(`${m.n_nodes.toLocaleString()} / ${m.n_edges.toLocaleString()}`, 'nodes / edges')),
    h('h3', {}, 'Evidence ladder: how strong is the genetic evidence for each CTEPH gene, and what does it do?')];
  const th = themes();
  const grid = h('table', {class: 'grid'});
  grid.append(h('thead', {}, h('tr', {}, h('th', {}, 'Level'), th.map(t => h('th', {title: t}, short(t))))));
  const tb = h('tbody');
  for (const [k, [code, label]] of Object.entries(D.ladder)) {
    const row = h('tr', {class: 'lv' + k}, h('td', {}, h('b', {}, `${k} ${code}`), h('br'), label));
    for (const t of th) {
      const gs = genes.filter(g => lvl(g) === Number(k) && (g.a.theme || 'Other genes') === t)
                      .sort((a, b) => a.label.localeCompare(b.label));
      row.append(h('td', {}, gs.map(g => h('a', {class: 'chip', href: href(g.id)}, g.label))));
    }
    tb.append(row);
  }
  grid.append(tb);
  out.push(h('div', {style: 'overflow-x:auto'}, grid),
           h('div', {class: 'note'}, 'Levels 1–2: the CTEPH GWAS (G; G* only after leveraging pulmonary embolism). ' +
             'Levels 3–7: the genetic studies of CTEPH; levels 3–6 count only findings about CTEPH itself.'));
  out.push(h('h3', {}, 'Which conditions share the CTEPH GWAS regions?'), regionMatrix());
  return out;
}

function regionMatrix() {
  const regions = byType('region');
  const rows = byType('condition').filter(c => {
    const s = (OUT.get(c.id) || []).filter(e => e.type === 'shares_region');
    return s.length >= 2 || (s.length >= 1 && c.a.category !== 'Shared genetic signal');
  }).map(c => ({c, s: new Map((OUT.get(c.id) || []).filter(e => e.type === 'shares_region').map(e => [e.t, e]))}))
    .sort((a, b) => b.s.size - a.s.size || CATS.indexOf(a.c.a.category) - CATS.indexOf(b.c.a.category) ||
                    a.c.label.localeCompare(b.c.label));
  const tbl = h('table', {class: 'matrix'});
  tbl.append(h('thead', {}, h('tr', {}, h('th', {}, 'Condition'), h('th', {}, 'Category'),
    regions.map(r => h('th', {class: 'rot', title: `${r.a.code} · tier ${r.a.tier}`}, r.label.replace(' region', ''))),
    h('th', {}, 'Regions'))));
  tbl.append(h('tbody', {}, rows.map(({c, s}) => h('tr', {}, h('td', {}, nodeLink(c.id)), h('td', {}, c.a.category),
    regions.map(r => { const e = s.get(r.id);
      return h('td', {class: 'c', title: e ? `${e.a.method}: ${e.a.n_signals} GWAS signals, best ${e.a.measure} ` +
        fmt(e.a.best) : ''}, e ? (e.a.method === 'colocalization' ? '●' : '○') : ''); }),
    h('td', {class: 'num'}, s.size)))));
  return h('div', {}, h('div', {style: 'overflow-x:auto'}, tbl),
           h('div', {class: 'note'}, '● a GWAS signal of the condition colocalizes with the CTEPH credible set; ○ a GWAS ' +
             'credible set of the condition contains the CTEPH lead variant (regions without a CTEPH credible set; ' +
             'weaker). Rows: conditions sharing two or more regions, and clinically related conditions sharing one.'));
}

function genesPage() {
  let all = false;
  const box = h('div');
  const draw = () => {
    const rows = byType('gene').filter(g => all || lvl(g) !== null);
    box.replaceChildren(table([
      {label: 'Gene', get: g => g.label, cell: g => nodeLink(g.id)},
      {label: 'CTEPH level', num: true, get: g => g.a.cteph_level, cell: g => g.a.cteph_level === undefined ? '' :
        h('span', {class: 'lv' + g.a.cteph_level, style: 'padding:0 6px'}, `${g.a.cteph_level} ${g.a.cteph_code || ''}`)},
      {label: 'Evidence', get: g => g.a.cteph_level_label},
      {label: 'Function', get: g => g.a.theme},
      {label: 'Related conditions', num: true, get: g => g.a.n_related_conditions},
      {label: 'Text-mined publications', num: true, get: g => g.a.n_text_mined_publications}],
      rows, {sort: all ? 4 : 1, dir: all ? -1 : 1, limit: 300, tools: h('span', {},
        h('button', {class: 'btn' + (all ? '' : ' on'), onclick: () => { all = false; draw(); }}, 'CTEPH genes'), ' ',
        h('button', {class: 'btn' + (all ? ' on' : ''), onclick: () => { all = true; draw(); }}, 'Every gene'))}));
  };
  draw();
  return [h('h2', {}, 'Genes'), h('p', {class: 'lead'}, 'The genes studied genetically in CTEPH, by evidence level; ' +
    '“Every gene” adds the genes with genetic evidence in a related condition and the text-mined CTEPH targets.'), box];
}

function conditionsPage() {
  const rows = byType('condition').filter(c => c.id !== CORE);
  return [h('h2', {}, 'Related conditions'), h('p', {class: 'lead'}, 'Every condition related to CTEPH, with the kinds ' +
    'of evidence for the relation. Open a condition for each piece of evidence and its source.'),
    table([{label: 'Condition', get: c => c.label, cell: c => nodeLink(c.id)},
           {label: 'Category', get: c => CATS.indexOf(c.a.category) + ' ' + c.a.category, cell: c => c.a.category},
           {label: 'Evidence types', get: c => c.a.sources},
           {label: 'Clinical items', num: true, get: c => c.a.n_clinical},
           {label: 'CTEPH regions shared', num: true, get: c => c.a.n_shared_regions,
            cell: c => c.a.shared_regions || ''},
           {label: 'GWAS studies', num: true, get: c => c.a.n_gwas_studies},
           {label: 'Curated genes', num: true, get: c => c.a.n_mendelian_genes}], rows, {sort: 1, limit: 200})];
}

function regionsPage() {
  return [h('h2', {}, 'CTEPH GWAS regions'), h('p', {class: 'lead'}, 'The nine genome-wide significant regions of the ' +
    'CTEPH GWAS (Liley 2024, Table 1), named by the gene the publication reports.'),
    table([{label: 'Region', get: r => r.label, cell: r => nodeLink(r.id)}, {label: 'Code', get: r => r.a.code},
           {label: 'Tier', num: true, get: r => r.a.tier}, {label: 'Lead variant', get: r => r.a.rsid,
            cell: r => rsid(r.a.rsid)}, {label: 'P', num: true, get: r => r.a.p},
           {label: 'V (with PE)', num: true, get: r => r.a.v}, {label: 'OR', num: true, get: r => r.a.or},
           {label: 'In Open Targets', get: r => r.a.open_targets},
           {label: 'Conditions sharing it', num: true, get: r => (IN.get(r.id) || []).filter(e =>
             e.type === 'shares_region').length}], byType('region'), {filter: false}),
    h('h3', {}, 'Conditions sharing them'), regionMatrix()];
}

function studiesPage() {
  const all = byType('publication').filter(p => p.a.design_class);
  const EAS = ['Japan', 'China', 'Korea', 'Taiwan', 'Hong Kong', 'Mongolia'];
  let place = 'all', design = 'all';
  const box = h('div');
  const designs = [...new Set(all.map(p => p.a.design_class))].sort();
  const draw = () => {
    const rows = all.filter(p => (place === 'all' || (place === 'Japan' ? p.a.place === 'Japan' : place === 'East Asia'
      ? EAS.includes(p.a.place) : !EAS.includes(p.a.place))) && (design === 'all' || p.a.design_class === design));
    const genesOf = p => [...new Set((OUT.get(p.id) || []).filter(e => e.type === 'reports' && N.get(e.t)?.type ===
      'gene').map(e => e.t))];
    box.replaceChildren(table([
      {label: 'Year', num: true, get: p => p.a.year}, {label: 'Study', get: p => p.label, cell: p => nodeLink(p.id)},
      {label: 'Place', get: p => p.a.place}, {label: 'Design class', get: p => p.a.design_class},
      {label: 'Cases', num: true, get: p => p.a.n_cases}, {label: 'Controls', num: true, get: p => p.a.n_controls},
      {label: 'Genes reported', get: p => genesOf(p).map(g => N.get(g).label).join(' '),
       cell: p => h('span', {}, genesOf(p).map(g => h('a', {class: 'chip', href: href(g)}, N.get(g).label)))},
      {label: 'Title', get: p => p.a.title}, {label: 'Source', get: p => p.id, cell: p => pmid(p.id.split(':')[1])}],
      rows, {sort: 0, limit: 200, tools: h('span', {}, ['all', 'Japan', 'East Asia', 'Elsewhere'].map(x =>
        h('button', {class: 'btn' + (place === x ? ' on' : ''), onclick: () => { place = x; draw(); }},
          x === 'all' ? 'All places' : x)), ' ',
        h('select', {onchange: e => { design = e.target.value; draw(); }},
          h('option', {value: 'all'}, 'All designs'), designs.map(d => h('option', {value: d,
            selected: d === design ? 'selected' : null}, d))))}));
  };
  draw();
  return [h('h2', {}, 'Genetic studies of CTEPH'), h('p', {class: 'lead'}, 'The genetic studies of CTEPH the curation ' +
    "keeps, with their place (country of the first author's institution) and design class."), box];
}

function curationPage() {
  const C = D.curation;
  const pubCols = [{label: 'Publication', get: r => `${r.first_author || ''} ${r.year || ''}`},
                   {label: 'Title', get: r => r.title}, {label: 'Reason', get: r => r.reason},
                   {label: 'Used in', get: r => r.roles}, {label: 'Source', get: r => r.pmid, cell: r => pmid(r.pmid)}];
  const findings = D.edges.filter(([, , t, a]) => t === 'reports' && a.about && a.about !== 'CTEPH')
    .map(([s, t, , a]) => ({s, t, a}));
  const log = C.log.map(([step, item, rule, before, after]) => ({step, item, rule, before, after}));
  const tmCols = [{label: 'Gene', get: r => r.item, cell: r => h('span', {class: 'gene'}, r.item)},
                  {label: 'Decision', get: r => r.rule}, {label: 'Links', get: r => r.before},
                  {label: 'What the word means', get: r => r.after.replace(/^(mapped to the wrong gene|words removed): /, '')}];
  const vnCols = [{label: 'Finding', get: r => r.item}, {label: 'As written', get: r => r.before},
                  {label: 'Name', get: r => r.after.replace(/^named: /, '')}];
  return [h('h2', {}, 'Curation'), h('p', {class: 'lead'}, `Every publication the graph uses ` +
    `(${C.counts.publications.toLocaleString()}) was looked up in PubMed on ${D.meta.checked_at}. Retracted articles ` +
    'and conference material are not part of the graph; findings that are not about CTEPH itself are kept and ' +
    'marked; words that text mining mapped to the wrong gene are removed. No rule calls a language model.'),
    h('h3', {}, `Findings not about CTEPH itself (${findings.length})`),
    table([{label: 'Study', get: r => N.get(r.s)?.label, cell: r => nodeLink(r.s)},
           {label: 'Gene', get: r => N.get(r.t)?.label, cell: r => N.get(r.t)?.type === 'gene' ? nodeLink(r.t) : ''},
           {label: 'Result', get: r => r.a.result}, {label: 'About', get: r => r.a.about},
           {label: 'Source sentence', get: r => r.a.quote, cell: r => h('span', {class: 'quote'}, r.a.quote)}],
          findings, {filter: false}),
    h('h3', {}, 'Text-mined links mapped to the wrong gene'), table(tmCols, log.filter(r => r.step === 'text mining')),
    h('h3', {}, 'Variants given one name'), table(vnCols, log.filter(r => r.rule === 'variant name')),
    h('h3', {}, `Publications flagged and kept (${C.flagged.length})`), table(pubCols, C.flagged)];
}

// ── routing, search ──────────────────────────────────────────────────────
const TABS = [['overview', 'Overview'], ['genes', 'Genes'], ['conditions', 'Related conditions'],
              ['regions', 'GWAS regions'], ['studies', 'Genetic studies'], ['curation', 'Curation']];
const PAGES = {overview, genes: genesPage, conditions: conditionsPage, regions: regionsPage, studies: studiesPage,
               curation: curationPage};
function render() {
  const hash = decodeURIComponent(location.hash.slice(2) || 'overview');
  const main = document.getElementById('main');
  let tab = hash.split('/')[0];
  if (hash.startsWith('n/')) { main.replaceChildren(nodePage(hash.slice(2))); tab = ''; }
  else main.replaceChildren(...[].concat((PAGES[tab] || overview)()));
  document.querySelectorAll('nav a').forEach(a => a.classList.toggle('on', a.dataset.tab === tab));
  window.scrollTo(0, 0);
}
document.getElementById('nav').append(...TABS.map(([k, l]) => h('a', {href: '#/' + k, 'data-tab': k}, l)));
document.getElementById('meta').textContent = `Snapshot ${D.meta.snapshot} · Open Targets ${D.meta.release} · ` +
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
render();
</script>
</body>
</html>
'''

if __name__ == '__main__':
    main()
