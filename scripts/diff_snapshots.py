#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : What changed since the previous snapshot of a disease: the newest
#           earlier snapshot (by fetch date) whose tables are complete is
#           compared with this one by stable keys (records added / removed):
#             Condition            node ID
#             Condition regions    node ID x the disease's GWAS regions it shares
#             Relation evidence    condition x evidence type x PMID / relation
#             <disease> GWAS region region x lead variant x the disease's credible
#                                  set in Open Targets
#             <disease> gene level gene x evidence level
#             Colocalization       disease region x condition x credible set
#             Lead-variant overlap disease region x condition x credible set
#             GWAS study           condition x study ID
#             GWAS credible set    condition x study-locus ID
#             Mendelian gene       condition x source x gene
#             Gene-condition link  gene x condition x evidence kind
#             <disease> genetic study  PMID (literature extraction)
#             <disease> finding    PMID x gene x variant x result x what it is about
#             Text-mined literature target x publication
#           <disease> is the short name of the disease (tidy.json of this
#           snapshot). An earlier snapshot built before the disease profile is
#           read under its old names (cteph_regions.tsv -> disease_regions.tsv,
#           cteph_gene -> disease_gene, ...). A first snapshot records that there
#           is nothing to compare with; a table, or a key column, missing on
#           either side (a table added to the pipeline later) is listed as not
#           compared, with the reason.
# Inputs  : this snapshot's tables; the snapshot root and the earlier snapshot names
# Outputs : changes.tsv (category, change, key, item, detail), changes.json (summary)
# ---------------------------------------------------------------------------
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import csv_list, legacy_path, log, read_tsv, write_tsv  # noqa: E402

COLS = ['category', 'change', 'key', 'item', 'detail']
# the tables of the snapshots built before the disease profile, by their new names (their columns are renamed by
# read_tsv: cteph_gene -> disease_gene, cteph_study_locus_id -> disease_study_locus_id, ...)
OLD_TABLES = {'disease_regions': 'cteph_regions', 'disease_genes': 'cteph_genes', 'disease_studies': 'cteph_studies',
              'disease_findings': 'cteph_findings', 'disease_targets': 'cteph_targets',
              'disease_methods': 'cteph_methods', 'disease_screening': 'cteph_screening'}


def region_detail(r):
    ot = 'credible set' if r['disease_study_locus_id'] else 'lead variant only' if r['variant_id'] else 'not held'
    return '; '.join(x for x in (r['rsid'], f'tier {r["tier"]}' if r['tier'] else '', f'Open Targets: {ot}') if x)


def specs(s):
    """table, category, key columns, item, detail; s: the short name of the disease."""
    return [
        ('conditions', 'Condition', ['node_id'], lambda r: r['name'], lambda r: f'{r["category"]}; {r["sources"]}'),
        ('conditions', 'Condition regions', ['node_id', 'shared_regions'], lambda r: r['name'],
         lambda r: r['shared_regions'] or f'no {s} region'),
        ('relations', 'Relation evidence', ['node_id', 'evidence_type', 'pmid', 'relation'],
         lambda r: r['node_id'], lambda r: r['detail']),
        ('disease_regions', f'{s} GWAS region', ['gene', 'rsid', 'disease_study_locus_id'], lambda r: r['gene'],
         region_detail),
        ('disease_genes', f'{s} gene level', ['symbol', 'level'], lambda r: r['symbol'],
         lambda r: f'level {r["level"]}: {r["level_label"]}; {r["theme"]}'),
        ('colocalisation', 'Colocalization', ['disease_gene', 'node_id', 'other_study_locus_id'],
         lambda r: f'{r["condition"]}: {r["disease_gene"]}', lambda r: f'{r["study_id"]}; CLPP {r["clpp"]}'),
        ('overlap', 'Lead-variant overlap', ['disease_gene', 'node_id', 'other_study_locus_id'],
         lambda r: f'{r["condition"]}: {r["disease_gene"]}', lambda r: f'{r["study_id"]}; PIP {r["pip"]}'),
        ('gwas_studies', 'GWAS study', ['node_id', 'study_id'], lambda r: f'{r["condition"]}: {r["study_id"]}',
         lambda r: f'{r["trait"]}; {r["first_author"]} {r["year"]}'.strip('; ')),
        ('gwas_loci', 'GWAS credible set', ['node_id', 'study_locus_id'],
         lambda r: f'{r["condition"]}: {r["rsid"] or r["variant_id"]}',
         lambda r: f'{r["study_id"]}; L2G {r["l2g1_symbol"]}'),
        ('mendelian', 'Mendelian gene', ['node_id', 'source', 'ensembl_id'],
         lambda r: f'{r["condition"]}: {r["symbol"]}', lambda r: r['source']),
        ('gene_condition', 'Gene-condition link', ['ensembl_id', 'node_id', 'evidence'],
         lambda r: f'{r["symbol"]} – {r["condition"]}', lambda r: r['evidence']),
        ('disease_studies', f'{s} genetic study', ['pmid'], lambda r: f'PMID {r["pmid"]}',
         lambda r: f'{r["first_author"]} {r["year"]}; {r["design"]}'),
        ('disease_findings', f'{s} finding', ['pmid', 'symbol', 'variant', 'result', 'association_with'],
         lambda r: f'PMID {r["pmid"]}',
         lambda r: f'{r["symbol"] or r["gene"]} {r["variant_name"]} {r["result"]} (about {r["association_with"]})'
         .strip()),
        ('literature', 'Text-mined literature', ['ensembl_id', 'pmid'], lambda r: f'{r["symbol"]} · {r["pmid"]}',
         lambda r: str(r['year'])),
    ]


# the columns the items and details read (filled with '' where a side does not have them; keys must exist on both)
SHOWN = {'conditions': ['name', 'category', 'sources', 'shared_regions'], 'relations': ['detail'],
         'disease_regions': ['tier', 'variant_id'], 'disease_genes': ['level_label', 'theme'],
         'colocalisation': ['condition', 'study_id', 'clpp'], 'overlap': ['condition', 'study_id', 'pip'],
         'gwas_studies': ['condition', 'trait', 'first_author', 'year'],
         'gwas_loci': ['condition', 'rsid', 'variant_id', 'study_id', 'l2g1_symbol'],
         'mendelian': ['condition', 'symbol'], 'gene_condition': ['symbol', 'condition'],
         'disease_studies': ['first_author', 'year', 'design'],
         'disease_findings': ['gene', 'variant_name'], 'literature': ['symbol', 'year']}


def table_path(d, table):
    """A table of a snapshot under its name, or under its name before the disease profile."""
    return legacy_path(d, f'{table}.tsv', f'{OLD_TABLES[table]}.tsv') if table in OLD_TABLES else d / f'{table}.tsv'


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
    tidy = json.loads((cur / 'tidy.json').read_text())
    short = (tidy.get('core') or {}).get('short') or 'Disease'
    prev_name, prev = previous(a.snapshots_root, csv_list(a.earlier), a.current)
    summary = {'current': a.current, 'current_release': tidy['release'],
               'previous': prev_name, 'previous_release': None, 'counts': {}, 'not_compared': []}
    rows, why, hidden = [], {}, {}
    excluded = set(tidy.get('excluded_pmids') or [])
    if prev_name is None:
        log('no earlier snapshot: nothing to compare with')
    else:
        summary['previous_release'] = json.loads((prev / 'tidy.json').read_text())['release']
        for table, cat, keys, item, detail in specs(short):
            paths = {'this': table_path(cur, table), 'earlier': table_path(prev, table)}
            gone = [side for side, p in paths.items() if not p.exists()]
            if gone:
                why[cat] = f'no {table} table in the {" and the ".join(gone)} snapshot'
                continue
            C, P = (read_tsv(paths[s]) for s in ('this', 'earlier'))
            lost = {side: [k for k in keys if k not in df] for side, df in (('this', C), ('earlier', P))}
            if any(lost.values()):
                why[cat] = '; '.join(f'no column {", ".join(v)} in the {side} snapshot'
                                     for side, v in lost.items() if v)
                continue
            for df in (C, P):
                for c in SHOWN.get(table, []):
                    if c not in df:
                        df[c] = ''
                df['_k'] = df[keys].agg('|'.join, axis=1) if len(df) else []
            for change, df, other in (('added', C, P), ('removed', P, C)):
                for r in df[~df._k.isin(set(other._k))].to_dict('records'):
                    if r.get('pmid') in excluded:       # left out by the curation since: counted, never listed
                        hidden[cat] = hidden.get(cat, 0) + 1
                        continue
                    rows.append({'category': cat, 'change': change, 'key': r['_k'], 'item': item(r),
                                 'detail': detail(r)})
        for cat, n in hidden.items():
            rows.append({'category': cat, 'change': 'removed', 'key': 'excluded publications',
                         'item': f'{n} rows of publications the curation now leaves out', 'detail': 'not listed'})
        summary['not_compared'] = list(why)
        summary['not_compared_reasons'] = why
        log(f'compared with {prev_name}: {len(rows)} changes' +
            (f'; not compared: {", ".join(why)}' if why else ''))
    df = pd.DataFrame(rows, columns=COLS)
    summary['counts'] = {f'{c}|{ch}': int(n) for (c, ch), n in df.groupby(['category', 'change']).size().items()}
    write_tsv(df, a.out)
    Path(a.out_json).write_text(json.dumps(summary, indent=1, ensure_ascii=False))


if __name__ == '__main__':
    main()
