#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Query the knowledge graph of a disease's snapshot (its SQLite file,
#           results/<key>/snapshots/<date>/05.graph/<key>_knowledge_graph.sqlite)
#           from the shell or from Python. Not a pipeline step: a reader of its
#           outputs. Python standard library only.
#             info                     the snapshot, its counts and files
#             gene <symbol | Ensembl>  a gene: evidence level, function, every
#                                      finding of the disease's genetic studies
#                                      with its study and sentence, the related
#                                      conditions with evidence, text-mined
#                                      publications
#             condition <name | ID>    a condition: why it relates to the disease
#                                      (each piece of evidence), the disease's
#                                      GWAS regions it shares, its genes
#             region <gene>            a GWAS region of the disease and the
#                                      conditions sharing it
#             publication <PMID>       a publication: status, study design and
#                                      place, findings, what it supports
#             studies [--place X] [--design X]
#                                      the genetic studies of the disease
#             search <text>            nodes whose label or ID contains the text
#             neighbors <node> [--edge TYPE]
#                                      every edge of a node
#             sql "<SELECT ...>"       any read-only query (tables: see
#                                      05.graph/README.md)
#           The disease: --disease <key> (default: the only disease under
#           results/, else the first in results/index.json); its newest snapshot
#           (results/<key>/latest.json) unless --snapshot names one; or any
#           database with --db. Snapshots built before the disease profile
#           (results/latest.json, columns cteph_level ...) are read as well.
#           Output: Markdown (default) or JSON (--json).
#           Python: from kg import KG; kg = KG(disease='cteph'); kg.gene('FGG')
#                   -> dict; KG(db='path/to/<key>_knowledge_graph.sqlite')
# Usage   : python3 scripts/kg.py [--disease KEY] [--snapshot YYYY-MM-DD |
#           --db PATH] [--json] <command> ...
# ---------------------------------------------------------------------------
import argparse
import json
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / 'results'
# gene columns of the databases built before the disease profile, read under their new names
LEGACY = {'cteph_level': 'disease_level', 'cteph_code': 'disease_code', 'cteph_level_label': 'disease_level_label',
          'cteph_ot_score': 'disease_ot_score'}


def diseases(root=ROOT):
    """The diseases with a snapshot under the results root (results/<key>/latest.json)."""
    return sorted(p.parent.name for p in Path(root).glob('*/latest.json'))


def locate(disease=None, snapshot=None, root=ROOT):
    """(snapshot directory, database, disease key) of a disease's snapshot: the disease named, else the only one
    under the results root, else the first of results/index.json; its newest snapshot (latest.json) unless one is
    named. The layout before several diseases (results/latest.json, results/snapshots/) is read when it is the only
    one there."""
    root = Path(root)
    if disease is None:
        keys = diseases(root)
        if len(keys) == 1:
            disease = keys[0]
        elif (root / 'index.json').exists() and json.loads((root / 'index.json').read_text()).get('diseases'):
            disease = json.loads((root / 'index.json').read_text())['diseases'][0]['key']
        elif keys:
            raise SystemExit(f'several diseases under {root}: {", ".join(keys)}; choose one with --disease')
        elif not (root / 'latest.json').exists():
            raise SystemExit(f'no snapshot under {root}: run the pipeline, or give --db')
    base = root / disease if disease else root
    latest = json.loads((base / 'latest.json').read_text()) if (base / 'latest.json').exists() else {}
    if disease and not latest and not snapshot:
        raise SystemExit(f'no snapshot of "{disease}": {base / "latest.json"} does not exist (diseases: '
                         f'{", ".join(diseases(root)) or "none"})')
    snap = snapshot or latest['snapshot']
    d = base / 'snapshots' / snap
    # the database: its name in latest.json (the same in every snapshot of a disease), <key>_knowledge_graph.sqlite,
    # else the only SQLite file of the snapshot's graph
    names = [Path(latest['database']).name] if latest.get('database') else []
    names += [f'{disease}_knowledge_graph.sqlite'] if disease else []
    found = [d / '05.graph' / n for n in names if (d / '05.graph' / n).exists()]
    if not found:
        found = sorted((d / '05.graph').glob('*.sqlite'))
        if len(found) != 1:
            raise SystemExit(f'no single graph database for snapshot {snap} in {d / "05.graph"}')
    key = disease or (latest.get('disease') or {}).get('key') or found[0].stem.replace('_knowledge_graph', '')
    return d, found[0], key


class KG:
    """The knowledge graph of one snapshot of one disease: KG(disease='cteph') (default: the only disease under
    results/, its newest snapshot), KG(disease='cteph', snapshot='2026-10-16'), or KG(db='path/to/<file>.sqlite')."""

    def __init__(self, disease=None, snapshot=None, db=None, root=ROOT):
        if db:
            db = Path(db)
            if not db.exists():
                raise SystemExit(f'no such database: {db}')
            self.dir = db.parent.parent if db.parent.name == '05.graph' else db.parent
            key = disease
        else:
            self.dir, db, key = locate(disease, snapshot, root)
        self.db = db
        self.con = sqlite3.connect(f'file:{db}?mode=ro', uri=True)
        self.con.row_factory = sqlite3.Row
        self.meta = {r['key']: json.loads(r['value']) for r in self.sql('SELECT key, value FROM meta')}
        g = self.meta.get('graph') or {}
        self.snapshot = g.get('snapshot') or snapshot or self.dir.name
        self.core = g.get('core', '')
        core = self.one('SELECT label, category FROM node_condition WHERE id = ?', (self.core,)) or {}
        # the short name: the graph's, else the category of the disease's own node (databases before the profile)
        self.short = g.get('short') or core.get('category') or self.core
        self.name = g.get('core_name') or core.get('label') or self.short
        self.key = g.get('key') or key or db.stem.replace('_knowledge_graph', '')
        self.about = ((self.meta.get('schema') or {}).get('about') or {}).get('default', self.short)
        cols = {r['name'] for r in self.sql("SELECT name FROM pragma_table_info('node_gene')")}
        self.legacy = 'disease_level' not in cols and 'cteph_level' in cols
        self.level = 'cteph_level' if self.legacy else 'disease_level'

    # ── plumbing ────────────────────────────────────────────────────────────
    def sql(self, query, params=()):
        return [dict(r) for r in self.con.execute(query, params).fetchall()]

    def one(self, query, params=()):
        r = self.con.execute(query, params).fetchone()
        return dict(r) if r else None

    def renamed(self, row):
        """A gene row of an old database under the new column names."""
        return {LEGACY.get(k, k): v for k, v in row.items()} if self.legacy and row else row

    def find(self, text, ntype=None):
        """A node by ID, by label (case-insensitive), else the nodes whose label or ID contains the text."""
        t = (' AND type = ?', (ntype,)) if ntype else ('', ())
        for q in ('id = ?', 'lower(label) = lower(?)'):
            rows = self.sql(f'SELECT * FROM nodes WHERE {q}{t[0]}', (text,) + t[1])
            if len(rows) == 1:
                return rows[0], []
            if rows:
                return None, rows
        rows = self.sql(f"SELECT * FROM nodes WHERE (lower(label) LIKE lower(?) OR lower(id) LIKE lower(?)){t[0]} "
                        'ORDER BY length(label) LIMIT 25', (f'%{text}%', f'%{text}%') + t[1])
        return (rows[0], []) if len(rows) == 1 else (None, rows)

    def node(self, ntype, text):
        n, cands = self.find(text, ntype)
        if n is None:
            raise LookupError(f'no single {ntype} for "{text}"' + (': ' + '; '.join(
                f'{c["label"]} ({c["id"]})' for c in cands[:12]) if cands else ''))
        return self.renamed(self.one(f'SELECT * FROM node_{ntype} WHERE id = ?', (n['id'],)))

    # ── cards ───────────────────────────────────────────────────────────────
    def info(self):
        cur = self.dir / '03.curated' / 'curation.json'
        cur = json.loads(cur.read_text()) if cur.exists() else {}
        stem = self.db.stem
        return {'disease': {'id': self.core, 'name': self.name, 'short': self.short, 'key': self.key},
                'snapshot': self.snapshot, 'dir': str(self.dir), 'graph': self.meta.get('graph'),
                'curation': {k: cur[k] for k in ('by_status', 'studies', 'studies_excluded', 'findings_by_subject',
                                                  'textmining_dropped') if k in cur},
                'files': {'viewer': str(self.dir / f'{stem}.html'),
                          'knowledge': str(self.dir / '08.knowledge' / f'{stem}.md'),
                          'workbook': str(next(self.dir.glob(f'{stem}_*.xlsx'), '')),
                          'graph': str(self.db.parent), 'database': str(self.db)}}

    def gene(self, text):
        g = self.node('gene', text)
        findings = self.sql(
            'SELECT p.label AS study, p.id AS publication, p.place, p.design_class, r.result, r.about, '
            'r.variant_name, r.rsid, r.effect, r.p_value, r.quote FROM edge_reports r '
            'JOIN node_publication p ON p.id = r.source WHERE r.target = ? ORDER BY p.year', (g['id'],))
        conds = self.sql(
            'SELECT c.label AS condition, c.id, c.category, GROUP_CONCAT(e.code, \'\') AS codes, '
            'GROUP_CONCAT(e.details, \' | \') AS details FROM edge_has_genetic_evidence e '
            'JOIN node_condition c ON c.id = e.target WHERE e.source = ? GROUP BY c.id ORDER BY c.id = ? DESC, '
            'c.category, c.label', (g['id'], self.core))
        regions = self.sql('SELECT r.* FROM edge_names_gene n JOIN node_region r ON r.id = n.source '
                           'WHERE n.target = ?', (g['id'],))
        lit = self.sql('SELECT p.label, p.id, p.title, m.matched_as FROM edge_mentions m JOIN node_publication p '
                       'ON p.id = m.source WHERE m.target = ? ORDER BY p.year DESC LIMIT 15', (g['id'],))
        return {'gene': g, 'regions': regions, 'findings': findings, 'conditions': conds,
                'text_mined': {'n': g['n_text_mined_publications'], 'latest': lit}}

    def condition(self, text):
        c = self.node('condition', text)
        ev = self.sql('SELECT evidence_type, relation, detail, quote, pmid, strength, east_asian FROM edge_related_to '
                      'WHERE source = ? ORDER BY evidence_type', (c['id'],))
        shared = self.sql('SELECT REPLACE(target, \'region:\', \'\') AS region, method, n_signals, best, measure '
                          'FROM edge_shares_region WHERE source = ? ORDER BY method, region', (c['id'],))
        lv = self.level
        genes = self.sql(
            f'SELECT g.label AS gene, GROUP_CONCAT(e.code, \'\') AS codes, MAX(e.n) AS n, g.{lv} AS disease_level '
            'FROM edge_has_genetic_evidence e JOIN node_gene g ON g.id = e.source WHERE e.target = ? '
            f'GROUP BY g.id ORDER BY g.{lv} IS NULL, g.{lv}, n DESC LIMIT 40', (c['id'],))
        gwas = self.one('SELECT COUNT(*) AS n FROM edge_gwas_of WHERE target = ?', (c['id'],))['n']
        return {'condition': c, 'evidence': ev, 'shared_regions': shared, 'genes': genes, 'gwas_studies': gwas}

    def region(self, text):
        r = self.node('region', text if text.startswith('region:') else f'region:{text}')
        conds = self.sql('SELECT c.label AS condition, c.category, s.method, s.n_signals, s.best, s.measure '
                         'FROM edge_shares_region s JOIN node_condition c ON c.id = s.source WHERE s.target = ? '
                         'ORDER BY s.method, s.best DESC', (r['id'],))
        genes = [x['label'] for x in self.sql('SELECT g.label FROM edge_names_gene n JOIN node_gene g '
                                              'ON g.id = n.target WHERE n.source = ?', (r['id'],))]
        labels = {'p': 'P', 'p_pe': 'P in the secondary analysis', 'v': 'V', 'or': 'OR', 'maf': 'MAF',
                  **((self.meta.get('schema') or {}).get('region_labels') or {})}
        return {'region': r, 'genes': genes, 'conditions': conds, 'labels': labels}

    def publication(self, text):
        p_ = text if ':' in text else f'PMID:{text}'
        pub = self.one('SELECT * FROM curated_publications WHERE pmid = ?', (p_.split(':', 1)[1],))
        node = self.one('SELECT * FROM node_publication WHERE id = ?', (p_,))
        if pub is None and node is None:
            raise LookupError(f'{text}: not a publication of this snapshot')
        study = self.one('SELECT * FROM curated_studies WHERE pmid = ?', (p_.split(':', 1)[1],))
        findings = self.sql('SELECT COALESCE(g.label, \'\') AS gene, r.result, r.about, r.variant_name, r.gene_set, '
                            'r.quote FROM edge_reports r LEFT JOIN node_gene g ON g.id = r.target WHERE r.source = ?',
                            (p_,))
        supports = self.sql('SELECT c.label AS condition, s.relation, s.detail FROM edge_supports s '
                            'JOIN node_condition c ON c.id = s.target WHERE s.source = ?', (p_,))
        mentions = [x['label'] for x in self.sql('SELECT g.label FROM edge_mentions m JOIN node_gene g '
                                                 'ON g.id = m.target WHERE m.source = ?', (p_,))]
        return {'publication': pub or node, 'node': node, 'study': study, 'findings': findings,
                'supports': supports, 'mentions': mentions}

    def studies(self, place=None, design=None):
        q, prm = 'SELECT pmid, year, first_author, place, design_class_label AS design, n_cases, n_controls, title ' \
                 'FROM curated_studies WHERE 1', []
        if place:
            q += ' AND lower(place) LIKE lower(?)'
            prm.append(f'%{place}%')
        if design:
            q += ' AND lower(design_class_label) LIKE lower(?)'
            prm.append(f'%{design}%')
        return self.sql(q + ' ORDER BY CAST(year AS INTEGER), pmid', prm)

    def search(self, text):
        return self.sql('SELECT id, type, label FROM nodes WHERE lower(label) LIKE lower(?) OR lower(id) LIKE lower(?) '
                        'ORDER BY type, length(label) LIMIT 60', (f'%{text}%', f'%{text}%'))

    def neighbors(self, text, edge=None):
        n, cands = self.find(text)
        if n is None:
            raise LookupError(f'no single node for "{text}"' + (': ' + '; '.join(
                f'{c["label"]} ({c["type"]} {c["id"]})' for c in cands[:12]) if cands else ''))
        e = (' AND e.type = ?', (edge,)) if edge else ('', ())
        out = self.sql('SELECT e.type, \'→\' AS dir, m.type AS node_type, m.label, m.id FROM edges e JOIN nodes m '
                       f'ON m.id = e.target WHERE e.source = ?{e[0]} ORDER BY e.type, m.label', (n['id'],) + e[1])
        inn = self.sql('SELECT e.type, \'←\' AS dir, m.type AS node_type, m.label, m.id FROM edges e JOIN nodes m '
                       f'ON m.id = e.source WHERE e.target = ?{e[0]} ORDER BY e.type, m.label', (n['id'],) + e[1])
        return {'node': n, 'edges': out + inn}


# ── Markdown ────────────────────────────────────────────────────────────────
def table(rows, cols=None, width=90):
    if not rows:
        return ['(none)']
    cols = cols or list(rows[0])
    def cut(v):
        v = '' if v is None else f'{v:.3g}' if isinstance(v, float) else str(v)
        return v.replace('|', '/').replace('\n', ' ')[:width]
    return ['| ' + ' | '.join(cols) + ' |', '|' + '---|' * len(cols)] + \
        ['| ' + ' | '.join(cut(r.get(c)) for c in cols) + ' |' for r in rows]


def given(rows, cols):
    """The columns that have a value in at least one row (a statistic of a published GWAS table that a region from
    Open Targets does not have is left out)."""
    return [c for c in cols if any(r.get(c) not in (None, '') for r in rows)]


def md(cmd, d, short='the disease'):
    s = short
    if cmd == 'gene':
        g = d['gene']
        lv = f'level {g["disease_level"]} ({g["disease_level_label"]})' if g['disease_level'] is not None else \
            f'not studied genetically in {s}'
        out = [f'# {g["label"]} ({g["id"]})', '', f'- {s} evidence: {lv}' +
               (f'; code {g["disease_code"]}' if g['disease_code'] else ''), f'- Function: {g["theme"] or "—"}',
               f'- Related conditions with genetic evidence: {g["n_related_conditions"]}',
               f'- Text-mined {s} publications: {g["n_text_mined_publications"]}' +
               (f'; Reactome: {g["pathways"]}' if g['pathways'] else ''), f'- {g["url"]}', '']
        if d['regions']:
            out += [f'## {s} GWAS region', ''] + table(d['regions'], given(
                d['regions'], ['label', 'rsid', 'tier', 'code', 'p', 'v', 'or', 'open_targets', 'source'])) + ['']
        out += [f'## Findings of the genetic studies of {s}', ''] + table(
            d['findings'], ['study', 'place', 'design_class', 'result', 'about', 'variant_name', 'effect', 'p_value',
                            'quote'], 70) + ['']
        out += [f'## Conditions with genetic evidence (G GWAS, G* {s} region significant only in a secondary '
                f'analysis, C curated, L {s} literature)', ''] + \
            table(d['conditions'], ['condition', 'category', 'codes', 'details'], 80) + ['']
        if d['text_mined']['latest']:
            out += ['## Latest text-mined publications', ''] + table(d['text_mined']['latest'],
                                                                     ['label', 'id', 'matched_as', 'title'], 70)
        return '\n'.join(out)
    if cmd == 'condition':
        c = d['condition']
        out = [f'# {c["label"]} ({c["id"]})', '', f'- Category: {c["category"]}', f'- Evidence types: '
               f'{c.get("sources") or "—"}', f'- GWAS studies: {d["gwas_studies"]}', f'- {c["url"]}', '',
               f'## Why it relates to {s}', ''] + table(d['evidence'], ['evidence_type', 'relation', 'detail',
                                                                       'pmid', 'strength'], 100)
        out += ['', f'## {s} GWAS regions shared', ''] + table(d['shared_regions'])
        out += ['', f'## Genes (codes as for genes; level = {s} evidence level)', ''] + table(d['genes'])
        return '\n'.join(out)
    if cmd == 'region':
        r = d['region']
        pos = (f'chr{r["chrom"]}:{r["pos_grch37"]}' if r.get('pos_grch37') not in (None, '') else
               f'chr{r["chrom"]}:{r["pos_grch38"]} (GRCh38)' if r.get('pos_grch38') not in (None, '') else
               f'chr{r["chrom"]}')
        # the statistics the region has: those of a published table (labels from the graph schema) where it gives them
        labels = d.get('labels') or {'p': 'P', 'p_pe': 'P in the secondary analysis', 'v': 'V', 'or': 'OR', 'maf': 'MAF'}
        num = lambda v: f'{v:.2g}' if isinstance(v, float) else v  # noqa: E731   6.8e-27, not 6.800000000000001e-27
        stats = ', '.join(f'{labels[k]} {num(r[k])}' for k in ('p', 'p_pe', 'v', 'or', 'maf')
                          if r.get(k) not in (None, ''))
        head = f'- Code {r["code"]}' + (f', tier {r["tier"]}' if r.get('tier') not in (None, '') else '') + \
            (f'; {stats}' if stats else '') + f'; Open Targets: {r["open_targets"]}' + \
            (f'; source: {r["source"]}' if r.get('source') else '')
        return '\n'.join([f'# {r["label"]} ({r["rsid"] or "no rsID"}, {pos})', '', head,
                          f'- Genes named: {", ".join(d["genes"]) or "—"}', '', '## Conditions sharing it', ''] +
                         table(d['conditions']))
    if cmd == 'publication':
        p, st = d['publication'], d['study']
        out = [f'# {p.get("first_author", "")} {p.get("year", "")}: {p.get("title", "")}', '',
               f'- Status: {p.get("status", "")}' + (f' ({p["reason"]})' if p.get('reason') else ''),
               f'- Journal: {p.get("journal", "")}; roles: {p.get("roles", "")}',
               f'- Affiliation: {p.get("affiliation", "") or "—"}']
        if st:
            out += [f'- Genetic study of {s}: {st["design_class_label"]} ({st["design"]}); place {st["place"] or "—"}; '
                    f'cases {st["n_cases"] or "—"}, controls {st["n_controls"] or "—"}' +
                    (f'; readers: {st["agreement"]}' if st.get('agreement') else '')]
        out += ['', '## Findings', ''] + table(d['findings'], None, 80)
        if d['supports']:
            out += ['', '## Clinical evidence it gives', ''] + table(d['supports'], None, 100)
        if d['mentions']:
            out += ['', f'Text-mined genes: {", ".join(d["mentions"])}']
        return '\n'.join(out)
    if cmd == 'neighbors':
        n = d['node']
        return '\n'.join([f'# {n["label"]} ({n["type"]} {n["id"]}): {len(d["edges"])} edges', ''] +
                         table(d['edges']))
    if cmd == 'info':
        return json.dumps(d, indent=1, ensure_ascii=False)
    return '\n'.join(table(d, None, 90))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--disease', help='the disease key, a folder of results/ (default: the only one, else the first '
                                      'in results/index.json)')
    ap.add_argument('--snapshot', help='YYYY-MM-DD (default: the newest)')
    ap.add_argument('--db', help='a graph database (05.graph/<key>_knowledge_graph.sqlite) instead of results/')
    ap.add_argument('--root', default=str(ROOT), help='the results root (default: results/ of the repository)')
    ap.add_argument('--json', action='store_true', help='JSON instead of Markdown')
    sub = ap.add_subparsers(dest='cmd', required=True)
    sub.add_parser('info')
    for c in ('gene', 'condition', 'region', 'publication', 'search'):
        sub.add_parser(c).add_argument('text')
    s = sub.add_parser('studies')
    s.add_argument('--place')
    s.add_argument('--design')
    s = sub.add_parser('neighbors')
    s.add_argument('text')
    s.add_argument('--edge', help='an edge type (05.graph/README.md)')
    sub.add_parser('sql').add_argument('query')
    a = ap.parse_args()
    kg = KG(disease=a.disease, snapshot=a.snapshot, db=a.db, root=a.root)
    try:
        if a.cmd == 'sql':
            if not a.query.lstrip().lower().startswith(('select', 'with', 'pragma')):
                raise SystemExit('read-only: SELECT / WITH / PRAGMA only')
            d = kg.sql(a.query)
        elif a.cmd == 'studies':
            d = kg.studies(a.place, a.design)
        elif a.cmd == 'neighbors':
            d = kg.neighbors(a.text, a.edge)
        elif a.cmd == 'info':
            d = kg.info()
        else:
            d = getattr(kg, a.cmd)(a.text)
    except LookupError as e:
        raise SystemExit(str(e))
    print(json.dumps(d, indent=1, ensure_ascii=False, default=str) if a.json else md(a.cmd, d, kg.short))


if __name__ == '__main__':
    main()
