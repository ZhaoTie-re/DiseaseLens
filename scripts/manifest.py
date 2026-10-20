#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Provenance of one snapshot of one disease, written once its
#           workbook exists:
#             run_manifest.json  the disease, the run's parameters and script md5s
#                                (from the .nf), the Open Targets release, the record
#                                counts, the workbook's rows per sheet and the change
#                                summary
#             latest.json        the newest complete snapshot of the disease (by fetch
#                                date), with the path of each product, relative to
#                                results/<key>/ (where it is published)
#             README.md          the index of the snapshot: where to start, what each
#                                folder holds, the headline counts, how to query
#           Project mode (--project <results root>): the index of every disease
#           under the results root, from each results/<key>/latest.json and the
#           tidy.json of its newest snapshot:
#             index.json         each disease (ID, names), its newest snapshot, the
#                                release, the counts and the product paths, relative
#                                to the results root
#             README.md          the same as a short table
#             index.html         the same as one self-contained page (no library,
#                                no network), each disease linked by relative paths
#                                to its viewer, workbook and knowledge base, so that
#                                a copied or zipped results folder works from disk
# Inputs  : the profile, the run json from the .nf, tidy.json, workbook.json,
#           changes.json, curation.json, the graph schema; or the results root
# Outputs : run_manifest.json, latest.json, README.md; or index.json, README.md,
#           index.html
# ---------------------------------------------------------------------------
import argparse
import datetime as dt
import html
import json
import sys
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import REPOSITORY, MISSING, csv_list, load_profile, read_tsv  # noqa: E402

# the paths of latest.json, in the order written (relative to results/<key>/)
PRODUCTS = ['dir', 'index', 'viewer', 'knowledge', 'graph', 'database', 'workbook']
# tidy.json keys of the snapshots built before the disease profile
LEGACY_TIDY = {'disease_genes_by_level': 'cteph_genes_by_level', 'disease_genetic_studies': 'cteph_genetic_studies',
               'disease_findings': 'cteph_findings', 'disease_genetic_records': 'cteph_genetic_records',
               'disease_targets': 'cteph_targets'}


def tget(tidy, key, default=None):
    """A count of tidy.json under its name, or under its name before the disease profile."""
    return tidy.get(key, tidy.get(LEGACY_TIDY.get(key, key), default))


def identity(P):
    return {'id': P.id, 'name': P.name, 'short': P.short, 'key': P.key}


def levels_of(lv, ladder):
    """The levels of the ladder in order: those of the graph schema, else those holding a gene."""
    return sorted(int(k) for k in (ladder or lv))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--profile', help='profile.json of the disease')
    ap.add_argument('--run-json', help='parameters and script md5s, from the .nf')
    ap.add_argument('--tidy')
    ap.add_argument('--workbook-stats')
    ap.add_argument('--changes-json')
    ap.add_argument('--snapshot', help="this snapshot's name (its fetch date)")
    ap.add_argument('--earlier', default='', help='the other snapshots that exist (comma-separated)')
    ap.add_argument('--snapshots-root', help='results/<key>/snapshots')
    ap.add_argument('--workbook', help="the workbook's file name")
    ap.add_argument('--curation', help='03.curated/curation.json')
    ap.add_argument('--graph-schema', help='05.graph/schema.json')
    ap.add_argument('--out', default='run_manifest.json')
    ap.add_argument('--out-latest', default='latest.json')
    ap.add_argument('--out-readme', default='README.md')
    ap.add_argument('--project', help='project mode: the results root (results/<key>/latest.json of every disease)')
    ap.add_argument('--out-index', default='index.json', help='project mode: the index of the diseases')
    ap.add_argument('--out-html', default='index.html', help='project mode: the index as a page (relative links)')
    ap.add_argument('--latest', nargs='*', default=[], help='project mode: KEY=latest.json of the diseases of this '
                    'run (read instead of results/<key>/latest.json, which may not be published yet)')
    a = ap.parse_args()
    if a.project:
        return project(Path(a.project), a.out_index, a.out_readme,
                       dict(x.split('=', 1) for x in a.latest), a.out_html)
    need = ['profile', 'run_json', 'tidy', 'workbook_stats', 'changes_json', 'snapshot', 'snapshots_root', 'workbook',
            'curation', 'graph_schema']
    missing = [f'--{k.replace("_", "-")}' for k in need if getattr(a, k) is None]
    if missing:
        ap.error('the snapshot mode needs ' + ', '.join(missing) + ' (or --project for the project index)')
    P = load_profile(a.profile)
    run = json.loads(Path(a.run_json).read_text())
    tidy = json.loads(Path(a.tidy).read_text())
    wb = json.loads(Path(a.workbook_stats).read_text())
    ch = json.loads(Path(a.changes_json).read_text())
    cur = json.loads(Path(a.curation).read_text())
    schema = json.loads(Path(a.graph_schema).read_text())
    graph = schema['graph']
    man = {'component': 'DiseaseLens', 'disease': identity(P), 'snapshot': a.snapshot, 'written_at':
           dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'), 'open_targets': {
               'data_version': tidy['release'], 'api_version': tidy['api_version'], 'fetched_at': tidy['fetched_at']},
           'counts': {k: v for k, v in tidy.items() if k not in ('release', 'api_version', 'fetched_at')},
           'workbook': {'file': a.workbook, 'rows': wb['rows']}, 'changes': ch,
           'curation': {k: cur[k] for k in ('checked_at', 'by_status', 'studies', 'studies_excluded',
                                             'findings_by_subject', 'textmining_dropped') if k in cur},
           'graph': graph, **run}
    Path(a.out).write_text(json.dumps(man, indent=1))
    names = sorted(set(csv_list(a.earlier)) | {a.snapshot})
    newest = names[-1]
    root = Path(a.snapshots_root)
    d = Path(root.name) / newest      # relative to results/<key>/, where latest.json is published
    latest = {'snapshot': newest, 'dir': str(d), 'index': str(d / 'README.md'), 'viewer': str(d / f'{P.stem}.html'),
              'knowledge': str(d / '08.knowledge' / f'{P.stem}.md'), 'graph': str(d / '05.graph'),
              'database': str(d / '05.graph' / f'{P.stem}.sqlite'), 'snapshots': names}
    if newest == a.snapshot:
        lv = tget(tidy, 'disease_genes_by_level', {}) or {}
        latest.update({'release': tidy['release'], 'workbook': str(d / a.workbook),
                       'counts': {'genes': sum(lv.values()) if lv else 0, 'genes_by_level': lv,
                                  'genetic_studies': tget(tidy, 'disease_genetic_studies'),
                                  'conditions': tidy.get('conditions'), 'regions': tidy.get('regions')}})
    else:   # an older date was re-run: the pointer stays on the newer snapshot, as its own manifest records it
        other = json.loads((root / newest / '_run_info' / 'run_manifest.json').read_text())
        latest.update({'release': other['open_targets']['data_version'],
                       'workbook': str(d / other['workbook']['file'])})
    latest['disease'] = identity(P)
    Path(a.out_latest).write_text(json.dumps(latest, indent=1, ensure_ascii=False))
    genes = read_tsv(Path(a.tidy).parent / 'disease_genes.tsv', ['symbol'])
    example = P.checks.get('example_gene') or (genes.symbol.iloc[0] if len(genes) else '')
    Path(a.out_readme).write_text(index(P, a.snapshot, a.workbook, tidy, cur, schema, run, wb, example))


def folders(P):
    s = P.short
    published = ' and published regions' if P.gwas else ''
    return [('00.raw/', f'Open Targets responses as fetched ({s}, its GWAS{published}; one file per related '
             'condition)'),
            ('00.pubmed/', 'PubMed and Europe PMC records as fetched: the genetic-study search, the seed, the Methods '
             'texts, and the PubMed record of every publication the graph uses (publications.json.gz)'),
            ('01.relations/', 'when the profile has no curated list of related conditions: the clinical relations '
             'read from PubMed by the model readers (the seed, the disputed and the dropped items)'),
            ('01.conditions/', 'the related conditions and the evidence for each relation'),
            ('02.literature/', f'the genetic studies of {s} as extracted (screening, studies, findings, methods)'),
            ('03.curated/', 'the curation: publication status, studies, findings and text-mined links as the graph '
             'uses them, and curation_log.tsv (every rule applied)'),
            ('04.tables/', 'the long tables: conditions, regions, the evidence ladder, GWAS, curated genes, '
             'literature'),
            ('05.graph/', 'the knowledge graph: typed nodes and edges, graph.json, the SQLite database, the schema '
             '(README.md)'),
            ('06.summaries/', 'the generated summaries (language model; cached, checked against their sources)'),
            ('07.changes/', 'the changes since the previous snapshot'),
            ('08.knowledge/', 'the knowledge base (Markdown and JSON)'),
            ('09.figures/', 'the figures (PNG, PDF) and their catalog'),
            ('_run_info/', 'the resolved profile of the disease (profile.json) and the run manifest '
             '(run_manifest.json)')]


def index(P, snap, workbook, tidy, cur, schema, run, wb, example):
    """The README of a snapshot: where to start, what each folder holds, the headline counts."""
    s, graph = P.short, schema['graph']
    lv = tget(tidy, 'disease_genes_by_level', {})
    levels = levels_of(lv, schema.get('ladder'))
    span = f'level {levels[0]}' if len(levels) == 1 else f'levels {levels[0]}–{levels[-1]}' if levels else 'levels'
    query = f'gene {example}' if example else 'info'
    choices = next((x.split(' ')[0] for x in wb.get('sheets') or [] if 'Choices' in x), None)
    md = [f'# {s} knowledge graph — snapshot {snap}', '',
          f'Open Targets {tidy["release"]} (read {str(tidy.get("fetched_at") or tidy["fetch_date"])[:10]}); PubMed '
          f'searched up to {tidy["fetch_date"]} (entry date); '
          f'every publication checked in PubMed on {cur["checked_at"][:10]}. Generated by '
          'DiseaseLens' + (' (rebuilt from the saved sources)' if run['params'].get('Rebuild')
                           else '') + '; do not edit by hand.', '',
          '## Start here', '',
          '| To | Open |', '|---|---|',
          '| browse the graph (search, a page per gene, condition, region and study, every source linked) | '
          f'`{P.stem}.html` (any browser, offline) |',
          f'| read the evidence by question, in four parts (A {s} genetics, B related conditions, C context, '
          f'D provenance) | `{workbook}` |',
          f'| read the knowledge base | `08.knowledge/{P.stem}.md` |',
          '| see the figures | `09.figures/` (catalog in its README) |',
          f'| query from code | `05.graph/{P.stem}.sqlite` (schema: `05.graph/README.md`), or '
          f'`python3 scripts/kg.py --disease {P.key} {query}` |', '',
          '## In this snapshot', '',
          f'- {s} genes on the evidence ladder: {sum(lv.values())} ({span}: ' +
          ', '.join(str(lv.get(str(k), lv.get(k, 0))) for k in levels) + ').',
          f'- Genetic studies of {s}: {tget(tidy, "disease_genetic_studies")} ({tget(tidy, "disease_findings")} '
          f'findings), from {tget(tidy, "disease_genetic_records")} PubMed records; {cur["studies_excluded"]} left '
          'out by the curation.']
    if 'regions' in tidy:
        n_pub, n_ot = tidy.get('published_regions', 0), tidy.get('regions_from_open_targets', 0)
        md += [f'- GWAS regions of {s}: {tidy["regions"]} ({n_pub} published, {n_ot} from the genome-wide '
               f'significant credible sets of the {s} GWAS in Open Targets).' if tidy['regions'] else
               f'- GWAS regions of {s}: none (no published region in the profile, and no genome-wide significant '
               f'credible set from a GWAS of {s} in Open Targets).']
    md += [f'- Related conditions: {tidy["conditions"]}; GWAS studies {tidy["gwas_studies"]}, credible sets '
           f'{tidy["credible_sets"]:,}, curated gene links {tidy["mendelian_links"]}.',
           f'- Publications checked: {cur["publications"]:,} ('
           f'{", ".join(f"{k} {v}" for k, v in cur["by_status"].items())}); text-mined links dropped: '
           f'{cur["textmining_dropped"]}.',
           f'- Graph: {graph["n_nodes"]:,} nodes and {graph["n_edges"]:,} edges (' +
           ', '.join(f'{k} {v}' for k, v in graph['counts']['nodes'].items()) + ').', '',
           '## Folders, in the order the pipeline writes them', '', '| Folder | Content |', '|---|---|']
    md += [f'| `{f}` | {c} |' for f, c in folders(P)]
    md += ['', 'How the graph is built: `docs/METHODS.md` of the DiseaseLens repository; every choice, with its '
           'value and basis, is also in ' +
           (f'sheet {choices} of the workbook.' if choices else 'the workbook.'), '']
    return '\n'.join(md)


# ── project mode: every disease under the results root ──────────────────────
def load(path):
    path = Path(path)
    return json.loads(path.read_text()) if path.exists() else {}


def disease_of(key, latest, snap_dir, tidy):
    """The disease of a results folder: from latest.json, else the snapshot's tidy.json, else its graph (a snapshot
    built before the disease profile: the disease's own node, whose category is its short name)."""
    d = dict(latest.get('disease') or tidy.get('core') or {})
    g = load(snap_dir / '05.graph' / 'schema.json').get('graph') or {}
    d.setdefault('id', g.get('core', ''))
    d.setdefault('name', g.get('core_name', ''))
    d.setdefault('short', g.get('short', ''))
    if (not d['name'] or not d['short']) and (snap_dir / '05.graph' / 'nodes_condition.tsv').exists():
        cond = read_tsv(snap_dir / '05.graph' / 'nodes_condition.tsv', ['id', 'label', 'category'])
        own = cond[cond.id == d['id']]
        if len(own):
            d['name'] = d['name'] or own.label.iloc[0]
            d['short'] = d['short'] or own.category.iloc[0]
    d['key'] = key
    return {k: d.get(k) or '' for k in ('id', 'name', 'short', 'key')}


def project(root, out_index, out_readme, given=None, out_html=None):
    """index.json, README.md and index.html of the results root: one entry per results/<key>/latest.json (given: the
    latest.json of the diseases of this run, read in their place)."""
    rows = []
    files = {lp.parent.name: lp for lp in root.glob('*/latest.json')} | {k: Path(v) for k, v in (given or {}).items()}
    for key, lp in sorted(files.items()):
        latest = load(lp)
        snap_dir = root / key / latest.get('dir', f'snapshots/{latest.get("snapshot", "")}')
        tidy = {} if latest.get('counts') else load(snap_dir / '04.tables' / 'tidy.json')
        lv = tget(tidy, 'disease_genes_by_level', {}) or {}
        counts = latest.get('counts') or {
            'genes': sum(lv.values()) if lv else None, 'genes_by_level': lv,
            'genetic_studies': tget(tidy, 'disease_genetic_studies'), 'conditions': tidy.get('conditions'),
            'regions': tidy.get('regions', tidy.get('published_regions'))}
        paths = {'latest': str(Path(key) / 'latest.json')}
        paths.update({k: str(Path(key) / latest[k]) for k in PRODUCTS if latest.get(k)})
        rows.append({**disease_of(key, latest, snap_dir, tidy), 'snapshot': latest.get('snapshot'),
                     'release': latest.get('release'), 'snapshots': latest.get('snapshots') or [], 'counts': counts,
                     'paths': paths})
    Path(out_index).write_text(json.dumps({'generated_by': 'DiseaseLens', 'diseases': rows}, indent=1,
                                          ensure_ascii=False))
    v = lambda x: MISSING if x is None or x == '' else f'{x:,}' if isinstance(x, int) else str(x)  # noqa: E731
    md = ['# DiseaseLens results', '',
          'One folder per disease. `<key>/latest.json` points to the newest snapshot of the disease '
          '(`<key>/snapshots/<date>/`), whose README is the place to start; `index.json` lists the same, with the '
          'path of every product. Generated by DiseaseLens; do not edit by hand.', '']
    if rows:
        md += ['| Disease | Open Targets ID | Newest snapshot | Release | Genes on the ladder | Genetic studies | '
               'Related conditions | GWAS regions | Start here |', '|---|---|---|---|---|---|---|---|---|']
        for r in rows:
            c, n = r['counts'], len(r['snapshots'])
            name = r['short'] + (f' ({r["name"]})' if r['name'] and r['name'].lower() != r['short'].lower() else '')
            md.append(f'| {name or r["key"]} | {v(r["id"])} | {v(r["snapshot"])} ({n} snapshot{"s" if n != 1 else ""})'
                      f' | {v(r["release"])} | {v(c["genes"])} | {v(c["genetic_studies"])} | {v(c["conditions"])} | '
                      f'{v(c["regions"])} | `{r["paths"].get("index", r["paths"]["latest"])}` |')
        gap = any(f'| {MISSING} ' in x for x in md[-len(rows):])
        md += ['', (f'{MISSING}: not in the snapshot\'s tidy.json. ' if gap else '') +
               'Query a disease from the shell: `python3 scripts/kg.py --disease <key> gene <symbol>` (repository '
               'README).', '']
    else:
        md += ['No disease has a snapshot under this folder yet.', '']
    Path(out_readme).write_text('\n'.join(md))
    if out_html:
        Path(out_html).write_text(index_html(rows, v))


def index_html(rows, v):
    """The project index as one page in the viewer's style: a row per disease (its snapshot, its headline counts, its
    viewer and other products), every link relative to the results root, where the page is written, so that the
    folder works from disk (file://), copied or unzipped. Every link opens a new tab, so that the index stays open."""
    e = html.escape

    def link(r, k, text, cls=None):
        p = r['paths'].get(k)
        return (f'<a{f" class={chr(34)}{cls}{chr(34)}" if cls else ""} href="{e(urllib.parse.quote(p))}" '
                f'target="_blank" rel="noopener">{text}</a>') if p else ''
    body = []
    for r in rows:
        c, n = r['counts'], len(r['snapshots'])
        short = r['short'] or r['key']
        name = r['name'] if r['name'] and r['name'].lower() != short.lower() else ''
        name = name[:1].upper() + name[1:]
        ot = (f'<a class="ext" href="https://platform.opentargets.org/disease/{e(r["id"])}" target="_blank" '
              f'rel="noopener">{e(r["id"])}</a>' if r['id'] else MISSING)
        lv = c.get('genes_by_level') or {}
        title = ', '.join(f'level {k}: {lv[k]}' for k in sorted(lv, key=int))
        files = ' · '.join(x for x in (link(r, 'workbook', 'Workbook'), link(r, 'knowledge', 'Knowledge base'),
                                       link(r, 'index', 'README')) if x)
        facts = ' · '.join(f'<span>{x}</span>' for x in [ot, f'snapshot {e(v(r["snapshot"]))}',
                                                          f'Open Targets {e(v(r["release"]))}']
                           + ([f'{n} snapshots'] if n > 1 else []))
        body.append(f'<tr><td class="dis"><div class="dz">{e(short)}</div><div class="name">{e(name)}</div>'
                    f'<div class="sub">{facts}</div></td>'
                    f'<td class="num"{f" title={chr(34)}{e(title)}{chr(34)}" if title else ""}>'
                    f'{e(v(c.get("genes")))}</td><td class="num">{e(v(c.get("genetic_studies")))}</td>'
                    f'<td class="num">{e(v(c.get("conditions")))}</td><td class="num">{e(v(c.get("regions")))}</td>'
                    f'<td class="open">{link(r, "viewer", "Open viewer", "btn")}<div class="files">{files}</div>'
                    '</td></tr>')
    table = ('<div class="panel"><table><thead><tr><th>Disease</th><th class="num">Genes on the ladder</th><th class="num">Genetic studies</th>'
             '<th class="num">Related conditions</th><th class="num">GWAS regions</th><th></th></tr></thead><tbody>'
             + ''.join(body) + '</tbody></table></div>') if rows else \
        '<p class="empty">No disease has a snapshot under this folder yet.</p>'
    return INDEX_HTML.replace('__TABLE__', table).replace('__REPO__', e(REPOSITORY))


INDEX_HTML = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>DiseaseLens</title>
<style>
:root { --panel: #f6f7f9; --line: #d9dde3; --ink: #1f2328; --ink2: #57606a; --link: #1f4e79; --head: #1f3b5f; }
* { box-sizing: border-box; }
body { margin: 0; background: #fff; color: var(--ink); font: 14px/1.45 -apple-system, "Segoe UI", Roboto,
       "Helvetica Neue", Arial, "Noto Sans", sans-serif; }
header { background: var(--head); color: #fff; padding: 14px 24px; }
header h1 { font-size: 19px; margin: 0; font-weight: 650; letter-spacing: .2px; }
main { padding: 22px 24px 32px; max-width: 1280px; margin: 0 auto; }
a { color: var(--link); text-decoration: none; } a:hover { text-decoration: underline; }
a.ext { white-space: nowrap; } a.ext::after { content: " \\2197"; font-size: 11px; }
.panel { border: 1px solid var(--line); border-radius: 6px; overflow-x: auto; }
table { border-collapse: collapse; width: 100%; font-size: 13.5px; }
th, td { padding: 10px 12px; text-align: left; vertical-align: middle; border-bottom: 1px solid var(--line); }
th { background: var(--panel); font-weight: 600; font-size: 12.5px; color: var(--ink2); vertical-align: bottom; }
th.num { white-space: nowrap; }
td.dis { min-width: 240px; }
.name { color: var(--ink); font-size: 13px; margin-top: 1px; }
tbody tr:last-child td { border-bottom: 0; }
tbody tr:hover td { background: #f8fafc; }
.num { text-align: right; font-variant-numeric: tabular-nums; }
td.num { font-size: 15px; font-weight: 600; color: var(--head); }
.dz { font-weight: 650; font-size: 15px; color: var(--head); }
.sub { color: var(--ink2); font-size: 12px; margin-top: 3px; } .sub span { white-space: nowrap; }
td.open { text-align: right; }
a.btn { display: inline-block; border: 1px solid var(--head); color: var(--head); border-radius: 4px; padding: 4px 12px;
        font-weight: 600; font-size: 12.5px; }
a.btn:hover { background: var(--head); color: #fff; text-decoration: none; }
.files { font-size: 12px; margin-top: 5px; color: var(--ink2); white-space: nowrap; }
.empty { color: var(--ink2); }
footer { color: var(--ink2); font-size: 12px; margin-top: 18px; }
/* a narrow window: the count headers wrap, so that the disease keeps its width */
@media (max-width: 1150px) { th.num { white-space: normal; } main { padding: 18px 14px 28px; } }
@media print {
  header { background: none; color: var(--ink); padding: 0 0 4px; border-bottom: 2px solid var(--head); }
  main { padding: 8px 0 0; max-width: none; }
  .panel { overflow: visible; border: 0; }
  a { color: inherit; } a.ext::after { content: ""; } a.btn { border: 0; padding: 0; }
  tr { break-inside: avoid; page-break-inside: avoid; }
}
</style>
</head>
<body>
<header><h1>DiseaseLens</h1></header>
<main>
__TABLE__
<footer>Generated by <a class="ext" href="__REPO__" target="_blank" rel="noopener">DiseaseLens</a></footer>
</main>
</body>
</html>
"""


if __name__ == '__main__':
    main()
