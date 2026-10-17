#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Provenance of one snapshot, written once its workbook exists:
#             run_manifest.json  the run's parameters and script md5s (from the .nf),
#                                the Open Targets release, the record counts, the
#                                workbook's rows per sheet and the change summary
#             latest.json        the newest complete snapshot under the results root
#                                (by fetch date), with the path of each product,
#                                relative to the results root
#             README.md          the index of the snapshot: where to start, what each
#                                folder holds, the headline counts, how to query
# Inputs  : the run json from the .nf, tidy.json, workbook.json, changes.json,
#           curation.json, the graph schema
# Outputs : run_manifest.json, latest.json, README.md
# ---------------------------------------------------------------------------
import argparse
import datetime as dt
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import csv_list  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-json', required=True, help='parameters and script md5s, from the .nf')
    ap.add_argument('--tidy', required=True)
    ap.add_argument('--workbook-stats', required=True)
    ap.add_argument('--changes-json', required=True)
    ap.add_argument('--snapshot', required=True, help="this snapshot's name (its fetch date)")
    ap.add_argument('--earlier', default='', help='the other snapshots that exist (comma-separated)')
    ap.add_argument('--snapshots-root', required=True)
    ap.add_argument('--workbook', required=True, help="the workbook's file name")
    ap.add_argument('--curation', required=True, help='03.curated/curation.json')
    ap.add_argument('--graph-schema', required=True, help='05.graph/schema.json')
    ap.add_argument('--out', default='run_manifest.json')
    ap.add_argument('--out-latest', default='latest.json')
    ap.add_argument('--out-readme', default='README.md')
    a = ap.parse_args()
    run = json.loads(Path(a.run_json).read_text())
    tidy = json.loads(Path(a.tidy).read_text())
    wb = json.loads(Path(a.workbook_stats).read_text())
    ch = json.loads(Path(a.changes_json).read_text())
    cur = json.loads(Path(a.curation).read_text())
    graph = json.loads(Path(a.graph_schema).read_text())['graph']
    man = {'component': 'cteph_knowledge_graph', 'snapshot': a.snapshot, 'written_at':
           dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'), 'open_targets': {
               'data_version': tidy['release'], 'api_version': tidy['api_version'], 'fetched_at': tidy['fetched_at']},
           'counts': {k: v for k, v in tidy.items() if k not in ('release', 'api_version', 'fetched_at')},
           'workbook': {'file': a.workbook, 'rows': wb['rows']}, 'changes': ch,
           'curation': {k: cur[k] for k in ('checked_at', 'by_status', 'studies', 'studies_excluded',
                                             'findings_by_subject', 'textmining_dropped')},
           'graph': graph, **run}
    Path(a.out).write_text(json.dumps(man, indent=1))
    names = sorted(set(csv_list(a.earlier)) | {a.snapshot})
    newest = names[-1]
    root = Path(a.snapshots_root)
    d = Path(root.name) / newest      # relative to the results root, where latest.json lives
    latest = {'snapshot': newest, 'dir': str(d), 'index': str(d / 'README.md'),
              'viewer': str(d / 'cteph_knowledge_graph.html'),
              'knowledge': str(d / '08.knowledge' / 'cteph_knowledge_graph.md'),
              'graph': str(d / '05.graph'), 'database': str(d / '05.graph' / 'cteph_knowledge_graph.sqlite'),
              'snapshots': names}
    if newest == a.snapshot:
        latest.update({'release': tidy['release'], 'workbook': str(d / a.workbook)})
    else:   # an older date was re-run: the pointer stays on the newer snapshot, as its own manifest records it
        other = json.loads((root / newest / '_run_info' / 'run_manifest.json').read_text())
        latest.update({'release': other['open_targets']['data_version'],
                       'workbook': str(d / other['workbook']['file'])})
    Path(a.out_latest).write_text(json.dumps(latest, indent=1))
    Path(a.out_readme).write_text(index(a.snapshot, a.workbook, tidy, cur, graph, run))


FOLDERS = [('00.raw/', 'Open Targets responses as fetched (CTEPH, its GWAS and published regions; one file per related '
            'condition)'),
           ('00.pubmed/', 'PubMed and Europe PMC records as fetched: the genetic-study search, the seed, the Methods '
            'texts, and the PubMed record of every publication the graph uses (publications.json.gz)'),
           ('01.conditions/', 'the related conditions and the evidence for each relation'),
           ('02.literature/', 'the genetic studies of CTEPH as extracted (screening, studies, findings, methods)'),
           ('03.curated/', 'the curation: publication status, studies, findings and text-mined links as the graph '
            'uses them, and curation_log.tsv (every rule applied)'),
           ('04.tables/', 'the long tables: conditions, regions, the evidence ladder, GWAS, curated genes, literature'),
           ('05.graph/', 'the knowledge graph: typed nodes and edges, graph.json, the SQLite database, the schema '
            '(README.md)'),
           ('06.summaries/', 'the generated summaries (language model; cached, checked against their sources)'),
           ('07.changes/', 'the changes since the previous snapshot'),
           ('08.knowledge/', 'the knowledge base (Markdown and JSON)'),
           ('09.figures/', 'the figures (PNG, PDF) and their catalog'),
           ('_run_info/', 'trace, report, timeline, DAG and the run manifest')]


def index(snap, workbook, tidy, cur, graph, run):
    """The README of a snapshot: where to start, what each folder holds, the headline counts."""
    lv = tidy['cteph_genes_by_level']
    md = [f'# CTEPH knowledge graph — snapshot {snap}', '',
          f'Open Targets {tidy["release"]} (read {tidy["fetch_date"]}); PubMed searched up to {tidy["fetch_date"]}; '
          f'every publication checked in PubMed on {cur["checked_at"][:10]}. Generated by '
          'DiseaseLens' + (' (rebuilt from the saved sources)' if run['params'].get('Rebuild')
                                                   else '') + '; do not edit by hand.', '',
          '## Start here', '',
          '| To | Open |', '|---|---|',
          '| browse the graph (search, a page per gene, condition, region and study, every source linked) | '
          '`cteph_knowledge_graph.html` (any browser, offline) |',
          f'| read the evidence by question, in four parts (A CTEPH genetics, B related conditions, C context, '
          f'D provenance) | `{workbook}` |',
          '| read the knowledge base | `08.knowledge/cteph_knowledge_graph.md` |',
          '| see the figures | `09.figures/` (catalog in its README) |',
          '| query from code | `05.graph/cteph_knowledge_graph.sqlite` (schema: `05.graph/README.md`), or '
          '`python3 scripts/kg.py gene FGG` |', '',
          '## In this snapshot', '',
          f'- CTEPH genes on the evidence ladder: {sum(lv.values())} (levels 1–7: ' +
          ', '.join(str(lv.get(str(k), lv.get(k, 0))) for k in range(1, 8)) + ').',
          f'- Genetic studies of CTEPH: {tidy["cteph_genetic_studies"]} ({tidy["cteph_findings"]} findings), from '
          f'{tidy["cteph_genetic_records"]} PubMed records; {cur["studies_excluded"]} left out by the curation.',
          f'- Related conditions: {tidy["conditions"]}; GWAS studies {tidy["gwas_studies"]}, credible sets '
          f'{tidy["credible_sets"]:,}, curated gene links {tidy["mendelian_links"]}.',
          f'- Publications checked: {cur["publications"]:,} ({", ".join(f"{k} {v}" for k, v in cur["by_status"].items())});'
          f' text-mined links dropped: {cur["textmining_dropped"]}.',
          f'- Graph: {graph["n_nodes"]:,} nodes and {graph["n_edges"]:,} edges (' +
          ', '.join(f'{k} {v}' for k, v in graph['counts']['nodes'].items()) + ').', '',
          '## Folders, in the order the pipeline writes them', '', '| Folder | Content |', '|---|---|']
    md += [f'| `{f}` | {c} |' for f, c in FOLDERS]
    md += ['', 'Methods, outputs and the basis of every choice: the component README and `docs/` '
           '(METHODS.md, OUTPUTS.md); every choice is also in sheet 11 of the workbook.', '']
    return '\n'.join(md)


if __name__ == '__main__':
    main()
