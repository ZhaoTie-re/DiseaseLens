#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Query the CTEPH knowledge graph of a snapshot (its SQLite file,
#           05.graph/cteph_knowledge_graph.sqlite) from the shell or from Python.
#           Not a pipeline step: a reader of its outputs.
#             info                     the snapshot, its counts and files
#             gene <symbol | Ensembl>  a gene: evidence level, function, every
#                                      finding of the CTEPH literature with its
#                                      study and sentence, the related conditions
#                                      with evidence, text-mined publications
#             condition <name | ID>    a condition: why it relates to CTEPH (each
#                                      piece of evidence), the CTEPH regions it
#                                      shares, its genes
#             region <gene>            a published CTEPH GWAS region and the
#                                      conditions sharing it
#             publication <PMID>       a publication: status, study design and
#                                      place, findings, what it supports
#             studies [--place X] [--design X]
#                                      the genetic studies of CTEPH
#             search <text>            nodes whose label or ID contains the text
#             neighbors <node> [--edge TYPE]
#                                      every edge of a node
#             sql "<SELECT ...>"       any read-only query (tables: see
#                                      05.graph/README.md)
#           Output: Markdown (default) or JSON (--json).
#           Python: from kg import KG; kg = KG(); kg.gene('FGG') -> dict
# Usage   : python3 scripts/kg.py [--snapshot YYYY-MM-DD] [--json] <command> ...
# ---------------------------------------------------------------------------
import argparse
import json
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / 'results'
CORE = 'EFO_0022954'


class KG:
    """The knowledge graph of one snapshot (default: the newest, results/latest.json)."""

    def __init__(self, snapshot=None, root=ROOT):
        root = Path(root)
        snap = snapshot or json.loads((root / 'latest.json').read_text())['snapshot']
        self.dir = root / 'snapshots' / snap
        db = self.dir / '05.graph' / 'cteph_knowledge_graph.sqlite'
        if not db.exists():
            raise SystemExit(f'no graph database for snapshot {snap}: {db}')
        self.snapshot = snap
        self.con = sqlite3.connect(f'file:{db}?mode=ro', uri=True)
        self.con.row_factory = sqlite3.Row

    # ── plumbing ────────────────────────────────────────────────────────────
    def sql(self, query, params=()):
        return [dict(r) for r in self.con.execute(query, params).fetchall()]

    def one(self, query, params=()):
        r = self.con.execute(query, params).fetchone()
        return dict(r) if r else None

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
        return self.one(f'SELECT * FROM node_{ntype} WHERE id = ?', (n['id'],))

    # ── cards ───────────────────────────────────────────────────────────────
    def info(self):
        meta = {r['key']: json.loads(r['value']) for r in self.sql('SELECT * FROM meta')}
        cur = json.loads((self.dir / '03.curated' / 'curation.json').read_text())
        return {'snapshot': self.snapshot, 'dir': str(self.dir), 'graph': meta['graph'],
                'curation': {k: cur[k] for k in ('by_status', 'studies', 'studies_excluded', 'findings_by_subject',
                                                  'textmining_dropped')},
                'files': {'viewer': str(self.dir / 'cteph_knowledge_graph.html'),
                          'knowledge': str(self.dir / '08.knowledge' / 'cteph_knowledge_graph.md'),
                          'workbook': str(next(self.dir.glob('cteph_knowledge_graph_*.xlsx'), '')),
                          'graph': str(self.dir / '05.graph')}}

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
            'c.category, c.label', (g['id'], CORE))
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
        genes = self.sql(
            'SELECT g.label AS gene, GROUP_CONCAT(e.code, \'\') AS codes, MAX(e.n) AS n, g.cteph_level '
            'FROM edge_has_genetic_evidence e JOIN node_gene g ON g.id = e.source WHERE e.target = ? '
            'GROUP BY g.id ORDER BY g.cteph_level IS NULL, g.cteph_level, n DESC LIMIT 40', (c['id'],))
        gwas = self.one('SELECT COUNT(*) AS n FROM edge_gwas_of WHERE target = ?', (c['id'],))['n']
        return {'condition': c, 'evidence': ev, 'shared_regions': shared, 'genes': genes, 'gwas_studies': gwas}

    def region(self, text):
        r = self.node('region', text if text.startswith('region:') else f'region:{text}')
        conds = self.sql('SELECT c.label AS condition, c.category, s.method, s.n_signals, s.best, s.measure '
                         'FROM edge_shares_region s JOIN node_condition c ON c.id = s.source WHERE s.target = ? '
                         'ORDER BY s.method, s.best DESC', (r['id'],))
        genes = [x['label'] for x in self.sql('SELECT g.label FROM edge_names_gene n JOIN node_gene g '
                                              'ON g.id = n.target WHERE n.source = ?', (r['id'],))]
        return {'region': r, 'genes': genes, 'conditions': conds}

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


def md(cmd, d):
    if cmd == 'gene':
        g = d['gene']
        lv = f'level {g["cteph_level"]} ({g["cteph_level_label"]})' if g['cteph_level'] is not None else \
            'not studied genetically in CTEPH'
        out = [f'# {g["label"]} ({g["id"]})', '', f'- CTEPH evidence: {lv}' +
               (f'; code {g["cteph_code"]}' if g['cteph_code'] else ''), f'- Function: {g["theme"] or "—"}',
               f'- Related conditions with genetic evidence: {g["n_related_conditions"]}',
               f'- Text-mined CTEPH publications: {g["n_text_mined_publications"]}' +
               (f'; Reactome: {g["pathways"]}' if g['pathways'] else ''), f'- {g["url"]}', '']
        if d['regions']:
            out += ['## CTEPH GWAS region', ''] + table(d['regions'], ['label', 'rsid', 'tier', 'code', 'p', 'v', 'or',
                                                                       'open_targets']) + ['']
        out += ['## Findings of the genetic studies of CTEPH', ''] + table(
            d['findings'], ['study', 'place', 'design_class', 'result', 'about', 'variant_name', 'effect', 'p_value',
                            'quote'], 70) + ['']
        out += ['## Conditions with genetic evidence (G GWAS, G* PE-adjusted CTEPH region, C curated, L CTEPH '
                'literature)', ''] + table(d['conditions'], ['condition', 'category', 'codes', 'details'], 80) + ['']
        if d['text_mined']['latest']:
            out += ['## Latest text-mined publications', ''] + table(d['text_mined']['latest'],
                                                                     ['label', 'id', 'matched_as', 'title'], 70)
        return '\n'.join(out)
    if cmd == 'condition':
        c = d['condition']
        out = [f'# {c["label"]} ({c["id"]})', '', f'- Category: {c["category"]}', f'- Evidence types: '
               f'{c.get("sources") or "—"}', f'- GWAS studies: {d["gwas_studies"]}', f'- {c["url"]}', '',
               '## Why it relates to CTEPH', ''] + table(d['evidence'], ['evidence_type', 'relation', 'detail',
                                                                       'pmid', 'strength'], 100)
        out += ['', '## CTEPH GWAS regions shared', ''] + table(d['shared_regions'])
        out += ['', '## Genes (codes as for genes; level = CTEPH evidence level)', ''] + table(d['genes'])
        return '\n'.join(out)
    if cmd == 'region':
        r = d['region']
        return '\n'.join([f'# {r["label"]} ({r["rsid"]}, chr{r["chrom"]}:{r["pos_grch37"]})', '',
                          f'- Code {r["code"]}, tier {r["tier"]}; P {r["p"]}, P in PE {r["p_pe"]}, V {r["v"]}, '
                          f'OR {r["or"]}, MAF {r["maf"]}; Open Targets: {r["open_targets"]}',
                          f'- Genes named: {", ".join(d["genes"])}', '', '## Conditions sharing it', ''] +
                         table(d['conditions']))
    if cmd == 'publication':
        p, s = d['publication'], d['study']
        out = [f'# {p.get("first_author", "")} {p.get("year", "")}: {p.get("title", "")}', '',
               f'- Status: {p.get("status", "")}' + (f' ({p["reason"]})' if p.get('reason') else ''),
               f'- Journal: {p.get("journal", "")}; roles: {p.get("roles", "")}',
               f'- Affiliation: {p.get("affiliation", "") or "—"}']
        if s:
            out += [f'- Genetic study of CTEPH: {s["design_class_label"]} ({s["design"]}); place {s["place"] or "—"}; '
                    f'cases {s["n_cases"] or "—"}, controls {s["n_controls"] or "—"}']
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
    ap.add_argument('--snapshot', help='YYYY-MM-DD (default: the newest)')
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
    kg = KG(a.snapshot)
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
    print(json.dumps(d, indent=1, ensure_ascii=False, default=str) if a.json else md(a.cmd, d))


if __name__ == '__main__':
    main()
