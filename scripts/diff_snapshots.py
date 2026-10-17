#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : What changed since the previous snapshot: the newest earlier
#           snapshot (by fetch date) whose tables are complete is compared with
#           this one by stable keys (records added / removed):
#             Condition            node ID
#             Condition regions    node ID x the CTEPH GWAS regions it shares
#             Relation evidence    condition x evidence type x PMID / relation
#             CTEPH GWAS region    region x lead variant (with its Open Targets status)
#             CTEPH gene level     gene x evidence level
#             Colocalization       CTEPH region x condition x credible set
#             Lead-variant overlap CTEPH region x condition x credible set
#             GWAS study           condition x study ID
#             GWAS credible set    condition x study-locus ID
#             Mendelian gene       condition x source x gene
#             Gene-condition link  gene x condition x evidence kind
#             CTEPH genetic study  PMID (literature extraction)
#             CTEPH finding        PMID x gene x variant x result x what it is about
#             Text-mined literature target x publication
#           A first snapshot records that there is nothing to compare with; a
#           table missing from the earlier snapshot (added to the component
#           later) is listed as not compared.
# Inputs  : this snapshot's tables; the snapshot root and the earlier snapshot names
# Outputs : changes.tsv (category, change, key, item, detail), changes.json (summary)
# ---------------------------------------------------------------------------
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import csv_list, log, read_tsv, write_tsv  # noqa: E402

COLS = ['category', 'change', 'key', 'item', 'detail']
# table, category, key columns, item, detail
SPECS = [
    ('conditions', 'Condition', ['node_id'], lambda r: r['name'], lambda r: f'{r["category"]}; {r["sources"]}'),
    ('conditions', 'Condition regions', ['node_id', 'shared_regions'], lambda r: r['name'],
     lambda r: r['shared_regions'] or 'no CTEPH region'),
    ('relations', 'Relation evidence', ['node_id', 'evidence_type', 'pmid', 'relation'],
     lambda r: r['node_id'], lambda r: r['detail']),
    ('cteph_regions', 'CTEPH GWAS region', ['gene', 'rsid', 'cteph_study_locus_id'], lambda r: r['gene'],
     lambda r: f'{r["rsid"]}; tier {r["tier"]}; Open Targets: ' +
     ('credible set' if r['cteph_study_locus_id'] else 'lead variant only' if r['variant_id'] else 'not held')),
    ('cteph_genes', 'CTEPH gene level', ['symbol', 'level'], lambda r: r['symbol'],
     lambda r: f'level {r["level"]}: {r["level_label"]}; {r["theme"]}'),
    ('colocalisation', 'Colocalisation', ['cteph_gene', 'node_id', 'other_study_locus_id'],
     lambda r: f'{r["condition"]}: {r["cteph_gene"]}', lambda r: f'{r["study_id"]}; CLPP {r["clpp"]}'),
    ('overlap', 'Lead-variant overlap', ['cteph_gene', 'node_id', 'other_study_locus_id'],
     lambda r: f'{r["condition"]}: {r["cteph_gene"]}', lambda r: f'{r["study_id"]}; PIP {r["pip"]}'),
    ('gwas_studies', 'GWAS study', ['node_id', 'study_id'], lambda r: f'{r["condition"]}: {r["study_id"]}',
     lambda r: f'{r["trait"]}; {r["first_author"]} {r["year"]}'.strip('; ')),
    ('gwas_loci', 'GWAS credible set', ['node_id', 'study_locus_id'],
     lambda r: f'{r["condition"]}: {r["rsid"] or r["variant_id"]}', lambda r: f'{r["study_id"]}; L2G {r["l2g1_symbol"]}'),
    ('mendelian', 'Mendelian gene', ['node_id', 'source', 'ensembl_id'], lambda r: f'{r["condition"]}: {r["symbol"]}',
     lambda r: r['source']),
    ('gene_condition', 'Gene-condition link', ['ensembl_id', 'node_id', 'evidence'],
     lambda r: f'{r["symbol"]} – {r["condition"]}', lambda r: r['evidence']),
    ('cteph_studies', 'CTEPH genetic study', ['pmid'], lambda r: f'PMID {r["pmid"]}',
     lambda r: f'{r["first_author"]} {r["year"]}; {r["design"]}'),
    ('cteph_findings', 'CTEPH finding', ['pmid', 'symbol', 'variant', 'result', 'association_with'],
     lambda r: f'PMID {r["pmid"]}',
     lambda r: f'{r["symbol"] or r["gene"]} {r["variant_name"]} {r["result"]} (about {r["association_with"]})'.strip()),
    ('literature', 'Text-mined literature', ['ensembl_id', 'pmid'], lambda r: f'{r["symbol"]} · {r["pmid"]}',
     lambda r: str(r['year'])),
]


def previous(root, names, current):
    for n in sorted((x for x in names if x < current), reverse=True):
        d = Path(root) / n / '04.tables'
        if (d / 'conditions.tsv').exists() and (d / 'tidy.json').exists():
            return n, d
    return None, None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tables', required=True, help="this snapshot's tables")
    ap.add_argument('--current', required=True, help="this snapshot's name (its fetch date)")
    ap.add_argument('--snapshots-root', required=True)
    ap.add_argument('--earlier', default='', help='names of the snapshots that exist (comma-separated)')
    ap.add_argument('--out', default='changes.tsv')
    ap.add_argument('--out-json', default='changes.json')
    a = ap.parse_args()
    cur = Path(a.tables)
    prev_name, prev = previous(a.snapshots_root, csv_list(a.earlier), a.current)
    summary = {'current': a.current, 'current_release': json.loads((cur / 'tidy.json').read_text())['release'],
               'previous': prev_name, 'previous_release': None, 'counts': {}, 'not_compared': []}
    rows = []
    if prev_name is None:
        log('no earlier snapshot: nothing to compare with')
    else:
        summary['previous_release'] = json.loads((prev / 'tidy.json').read_text())['release']
        for table, cat, keys, item, detail in SPECS:
            if not (prev / f'{table}.tsv').exists() or not set(keys) <= set(read_tsv(prev / f'{table}.tsv')):
                summary['not_compared'].append(cat)
                continue
            C, P = read_tsv(cur / f'{table}.tsv'), read_tsv(prev / f'{table}.tsv')
            for df in (C, P):
                df['_k'] = df[keys].agg('|'.join, axis=1) if len(df) else []
            for change, df, other in (('added', C, P), ('removed', P, C)):
                for r in df[~df._k.isin(set(other._k))].to_dict('records'):
                    rows.append({'category': cat, 'change': change, 'key': r['_k'], 'item': item(r),
                                 'detail': detail(r)})
        log(f'compared with {prev_name}: {len(rows)} changes')
    df = pd.DataFrame(rows, columns=COLS)
    summary['counts'] = {f'{c}|{ch}': int(n) for (c, ch), n in df.groupby(['category', 'change']).size().items()}
    write_tsv(df, a.out)
    Path(a.out_json).write_text(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
