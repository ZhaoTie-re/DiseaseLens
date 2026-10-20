#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The knowledge graph of one snapshot of one disease, typed and
#           documented, for programs and for the viewer: every node and edge comes
#           from the curated literature (03.curated) and the long tables
#           (04.tables).
#             nodes  condition · gene · region (a GWAS region of the disease) ·
#                    publication · gwas_study · theme (function of a gene of the
#                    disease) · pathway (Reactome top level of a text-mined target)
#             edges  related_to · supports · has_genetic_evidence · reports ·
#                    studies · region_of · names_gene · shares_region ·
#                    in_theme · in_pathway · mentions · gwas_of
#           (schema_of() below: what each type means, its identifier and columns).
# Inputs  : 03.curated, 04.tables (with ladder.json and themes.json), the profile
# Outputs : nodes_<type>.tsv, edges_<type>.tsv, nodes.tsv / edges.tsv (all, with
#           the common columns), graph.json (node-link, networkx-compatible),
#           schema.json, README.md (the schema and example queries), and
#           <disease>_knowledge_graph.sqlite (the node and edge tables, every
#           curated and long table, and the schema)
# ---------------------------------------------------------------------------
import argparse
import json
import re
import sqlite3
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import a_, load_ladder, load_profile, load_themes, log, num, read_tsv, write_tsv  # noqa: E402


def schema_of(P):
    """The node and edge types of the graph of a disease: {type: (identifier, meaning, columns)} and
    {type: (source type, target type, meaning, columns)}."""
    s = P.short
    nodes = {
        'condition': ('the condition\'s Open Targets ID (EFO_, MONDO_, HP_, Orphanet_); custom:<name> for a condition '
                      'without one (a curated grouping, or a condition the readers found that Open Targets does not '
                      'name)', f'{s} ({P.id}) and the conditions related to it',
                      ['category', 'sources', 'n_clinical', 'shared_regions', 'n_shared_regions', 'n_gwas_studies',
                       'n_eas_studies', 'n_credible_sets', 'n_mendelian_genes', 'genetics', 'genetics_note',
                       'ot_ids']),
        'gene': ('Ensembl gene ID', f'a gene with genetic evidence for {s} or a related condition, or a text-mined '
                 f'{s} target',
                 ['disease_level', 'disease_code', 'disease_level_label', 'disease_level_note', 'theme',
                  'n_related_conditions',
                  'disease_ot_score', 'n_text_mined_publications', 'pathways']),
        'region': ('region:<its name: the gene the publication names, else its top locus-to-gene gene, else the '
                   'position of its lead variant>', f'a genome-wide significant GWAS region of {s} (published, '
                   'or from Open Targets)', ['rsid', 'chrom', 'pos_grch37', 'pos_grch38', 'tier', 'code', 'p', 'p_pe',
                                             'v', 'or', 'maf', 'open_targets', 'source']),
        'publication': ('PMID:<n> (Europe PMC IDs as EPMC:<id>)', f'a genetic study of {s}, a publication the clinical '
                        f'evidence cites, or a text-mined {s} publication',
                        ['year', 'first_author', 'title', 'journal', 'roles', 'status', 'place', 'japanese',
                         'east_asian', 'design_class', 'n_cases', 'n_controls', 'agreement']),
        'gwas_study': ('the GWAS study ID (GCST..., FINNGEN_...)', f'a GWAS of {s} or of a related condition (a study '
                       'of two conditions has one node and two gwas_of edges)',
                       ['condition', 'trait', 'first_author', 'year', 'pmid', 'n_samples', 'ancestry', 'eas',
                        'n_credible_sets']),
        'theme': ('theme:<key>', f'the function of a gene of {s}', ['short', 'description']),
        'pathway': ('reactome:<top-level pathway>', f'a Reactome top-level pathway of text-mined {s} targets', []),
    }
    edges = {
        'related_to': ('condition', 'condition', f'a condition is related to {s}; one edge per piece of evidence',
                       ['evidence_type', 'relation', 'detail', 'quote', 'pmid', 'strength', 'east_asian', 'agreement']),
        'supports': ('publication', 'condition', 'a publication is the clinical evidence of a relation',
                     ['relation', 'detail']),
        'has_genetic_evidence': ('gene', 'condition', f'genetic evidence links a gene to a condition: G (GWAS; for {s} '
                                 'a gene of its GWAS regions), G* (a region significant only in a secondary analysis), '
                                 f'C (curated), L (a positive finding about {s} in a genetic study)',
                                 ['code', 'kind', 'n', 'details', 'eas']),
        'reports': ('publication', 'gene|condition', f'a genetic study of {s} reports a finding on a gene (findings '
                    f'without one gene point to {s})', ['result', 'about', 'variant', 'variant_name', 'rsid', 'effect',
                                                        'p_value', 'gene_set', 'quote', 'agreement']),
        'studies': ('publication', 'condition', f'a genetic study of {s}', ['design_class', 'design', 'n_cases',
                                                                           'n_controls', 'place', 'agreement']),
        'region_of': ('region', 'condition', f'a GWAS region of {s}', ['code', 'tier', 'p']),
        'names_gene': ('region', 'gene', 'the gene that names a region (published, or its top locus-to-gene gene)', []),
        'shares_region': ('condition', 'region', f'a GWAS signal of a condition is shared with a region of {s}',
                          ['method', 'n_signals', 'best', 'measure']),
        'in_theme': ('gene', 'theme', f'the function of a gene of {s}', []),
        'in_pathway': ('gene', 'pathway', f'a Reactome top-level pathway of a text-mined {s} target', []),
        'mentions': ('publication', 'gene', f'{a_(s)} publication mentions a gene (Europe PMC text mining, curated)',
                     ['matched_as', 'n_sentences', 'evidence_score', 'flag']),
        'gwas_of': ('gwas_study', 'condition', 'a GWAS of the condition', []),
    }
    return nodes, edges


def typed(df):
    """Columns whose every non-empty value is a number become numeric ('' becomes NULL), for SQL comparisons;
    identifiers that look like numbers (PMIDs, positions) stay text."""
    df = df.copy()
    for c in df.columns:
        if c in ('id', 'source', 'target', 'pmid', 'pmids', 'study_id', 'chrom', 'pos_grch37') or c.endswith('_id'):
            continue
        s = df[c].astype(str).str.strip()
        v = pd.to_numeric(s.where(s.ne(''), None), errors='coerce')
        if s.ne('').any() and v[s.ne('')].notna().all():
            df[c] = v.astype('Int64') if (v.dropna() % 1 == 0).all() and v.dropna().abs().max() < 2 ** 53 else v
    return df


def pid(p):
    p = str(p)
    return f'PMID:{p}' if p.isdigit() else f'EPMC:{p}'


def url(kind, i):
    if kind == 'publication':
        n = i.split(':', 1)[1]
        return f'https://pubmed.ncbi.nlm.nih.gov/{n}/' if i.startswith('PMID:') else f'https://europepmc.org/search?query={n}'
    if kind == 'gene':
        return f'https://platform.opentargets.org/target/{i}'
    if kind == 'condition' and not i.startswith('custom:'):
        return f'https://platform.opentargets.org/disease/{i}'
    if kind == 'gwas_study':
        return f'https://platform.opentargets.org/study/{i}'
    return ''


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--profile', required=True, help='profile.json of the disease')
    ap.add_argument('--curated', required=True)
    ap.add_argument('--tables', required=True)
    ap.add_argument('--out-dir', default='05.graph')
    a = ap.parse_args()
    P = load_profile(a.profile)
    NODES, EDGES = schema_of(P)
    Q, T, out = Path(a.curated), Path(a.tables), Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    r = lambda d, n, cols=None: read_tsv(d / f'{n}.tsv', cols)  # noqa: E731
    cond, rel, co, ov = r(T, 'conditions'), r(T, 'relations'), r(T, 'colocalisation'), r(T, 'overlap')
    rg, lad, gc = r(T, 'disease_regions'), r(T, 'disease_genes'), r(T, 'gene_condition')
    gw, tg, lit = r(T, 'gwas_studies'), r(T, 'disease_targets'), r(T, 'literature')
    st, fi, pubs = r(Q, 'studies', ['agreement']), r(Q, 'findings', ['agreement']), r(Q, 'publications')
    rel = rel if 'agreement' in rel else rel.assign(agreement='')
    tidy = json.loads((T / 'tidy.json').read_text())
    LADDER, CODES, KIND = load_ladder(T)
    themes = load_themes(T)
    CORE, short = P.id, P.short
    core_name = (tidy.get('core') or {}).get('name') or P.name
    N = {k: [] for k in NODES}
    E = {k: [] for k in EDGES}
    seen = set()

    def node(ntype, i, label, **kw):
        if i in seen:
            return
        seen.add(i)
        N[ntype].append({'id': i, 'type': ntype, 'label': label, 'url': url(ntype, i), **kw})

    def edge(etype, s, t, **kw):
        E[etype].append({'source': s, 'target': t, 'type': etype, **kw})

    # ── conditions ──────────────────────────────────────────────────────────
    node('condition', CORE, core_name, category=short)
    for c in cond.itertuples():
        node('condition', c.node_id, c.name, **{k: getattr(c, k) for k in NODES['condition'][2]})
    alias = {i: c.node_id for c in cond.itertuples() for i in c.ot_ids.split()} | {CORE: CORE}
    alias.update({c.node_id: c.node_id for c in cond.itertuples()})

    # ── publications: genetic studies, seed evidence, text-mined ────────────
    pi = pubs.set_index('pmid')
    roles = {}
    for p in st.pmid:
        roles.setdefault(p, set()).add(f'genetic study of {short}')
    for p in rel.pmid[rel.pmid.ne('')]:
        roles.setdefault(p, set()).add('clinical evidence')
    for p in lit.pmid[lit.pmid.ne('')]:
        roles.setdefault(p, set()).add('text-mined')
    sti = st.set_index('pmid')
    for p, ro in roles.items():
        x = pi.loc[p] if p in pi.index else None
        s = sti.loc[p] if p in sti.index else None
        assert x is None or x.status != 'excluded', f'an excluded publication in the graph: {p}'
        get = lambda k: '' if x is None else x[k]  # noqa: E731
        surname = re.sub(r'\s+[A-Z]{1,3}$', '', get('first_author'))     # 'Arthur Ataam J' -> 'Arthur Ataam'
        label = (f'{surname} {get("year")}'.strip() if x is not None else p)
        node('publication', pid(p), label or p, year=get('year'), first_author=get('first_author'),
             title=get('title'), journal=get('journal'), roles='; '.join(sorted(ro)),
             status=(get('reason') if x is not None and x.status == 'flagged' else ''),
             place=(s.place if s is not None else get('institution_country')),
             japanese=(s.japanese if s is not None else ''), east_asian=(s.east_asian if s is not None else ''),
             design_class=(s.design_class_label if s is not None else ''),
             n_cases=(s.n_cases if s is not None else ''), n_controls=(s.n_controls if s is not None else ''),
             agreement=(s.agreement if s is not None else ''))

    # ── genes ───────────────────────────────────────────────────────────────
    L = lad.set_index('symbol')
    nrel = gc[gc.node_id != CORE].groupby('ensembl_id').node_id.nunique()
    tgi = tg.set_index('ensembl_id')
    nlit = lit.groupby('ensembl_id').pmid.nunique()
    sym = dict(zip(gc.ensembl_id, gc.symbol)) | dict(zip(lit.ensembl_id, lit.symbol)) | \
        {e: s for e, s in zip(lad.ensembl_id, lad.symbol) if e}
    for e, s in sorted(sym.items(), key=lambda x: x[1]):
        lv = L.loc[s] if s in L.index else None
        t = tgi.loc[e] if e in tgi.index else None
        node('gene', e, s, disease_level=(int(lv.level) if lv is not None else ''),
             disease_code=(lv.code if lv is not None else ''),
             disease_level_label=(lv.level_label if lv is not None else ''),
             disease_level_note=(lv.get('level_note', '') if lv is not None else ''),
             theme=(lv.theme if lv is not None else ''), n_related_conditions=int(nrel.get(e, 0)),
             disease_ot_score=(t.score if t is not None else ''), n_text_mined_publications=int(nlit.get(e, 0)),
             pathways=(t.themes if t is not None and nlit.get(e, 0) else ''))
    gid = {s: e for e, s in sym.items()}

    # ── regions ─────────────────────────────────────────────────────────────
    for x in rg.to_dict('records'):
        i = f'region:{x["gene"]}'
        node('region', i, f'{x["gene"]} region', code=x['code'], open_targets=(
            'credible set' if x['disease_study_locus_id'] else 'lead variant only' if x['variant_id'] else 'not held'),
            **{k: x.get(k, '') for k in ('rsid', 'chrom', 'pos_grch37', 'pos_grch38', 'tier', 'p', 'p_pe', 'v', 'or',
                                         'maf', 'source')})
        edge('region_of', i, CORE, code=x['code'], tier=x.get('tier', ''), p=x['p'])
        for g in x['genes'].split():
            if g in gid:
                edge('names_gene', i, gid[g])
    shared = pd.concat([co.assign(method='colocalization', measure='CLPP', val=num(co.clpp)),
                        ov.assign(method='lead-variant overlap', measure='PIP', val=num(ov.pip))])
    if len(shared):
        for (cid, g), d in shared.assign(cid=lambda d: d.node_id.map(alias)).dropna(subset=['cid']).groupby(
                ['cid', 'disease_gene']):
            m = 'colocalization' if (d.method == 'colocalization').any() else 'lead-variant overlap'
            dd = d[d.method == m]
            edge('shares_region', cid, f'region:{g}', method=m, n_signals=int(dd.study_id.nunique()),
                 best=float(dd.val.max()), measure=dd.measure.iloc[0])

    # ── relations and their publications ────────────────────────────────────
    for x in rel.itertuples():
        edge('related_to', x.node_id, CORE, evidence_type=x.evidence_type, relation=x.relation, detail=x.detail,
             quote=x.quote, pmid=x.pmid, strength=x.strength, east_asian=x.east_asian, agreement=x.agreement)
        if x.pmid:
            edge('supports', pid(x.pmid), x.node_id, relation=x.relation, detail=x.detail)

    # ── genetic evidence, studies, findings ─────────────────────────────────
    for x in gc.itertuples():
        assert x.evidence in KIND, f'an evidence label without a code: {x.evidence}'
        edge('has_genetic_evidence', x.ensembl_id, x.node_id, code=KIND[x.evidence], kind=x.evidence, n=int(x.n),
             details=x.details, eas=x.eas)
    for s in st.itertuples():
        edge('studies', pid(s.pmid), CORE, design_class=s.design_class_label, design=s.design, n_cases=s.n_cases,
             n_controls=s.n_controls, place=s.place, agreement=s.agreement)
    for f in fi.itertuples():
        tgt = f.ensembl_id if f.symbol and f.ensembl_id else CORE
        edge('reports', pid(f.pmid), tgt, result=f.result, about=f.association_with, variant=f.variant,
             variant_name=f.variant_name, rsid=f.rsid, effect=f.effect, p_value=f.p_value,
             gene_set=f.gene_set or ('' if f.symbol else f.gene), quote=f.quote, agreement=f.agreement)

    # ── function, pathways, text mining, GWAS studies ───────────────────────
    for t in themes:
        node('theme', f'theme:{t["key"]}', t['name'], short=t.get('short', t['name']), description=t['description'])
    tkey = {t['name']: t['key'] for t in themes}
    for x in lad.itertuples():
        if x.theme in tkey and x.symbol in gid:
            edge('in_theme', gid[x.symbol], f'theme:{tkey[x.theme]}')
    for x in tg[tg.ensembl_id.isin(set(lit.ensembl_id))].itertuples():
        for th in [t for t in x.themes.split('; ') if t]:
            node('pathway', f'reactome:{th}', th)
            edge('in_pathway', x.ensembl_id, f'reactome:{th}')
    for x in lit.itertuples():
        edge('mentions', pid(x.pmid), x.ensembl_id, matched_as=x.matched_as, n_sentences=x.n_sentences,
             evidence_score=x.evidence_score, flag=x.publication_flag)
    for x in gw.itertuples():
        node('gwas_study', x.study_id, f'{x.trait} ({x.first_author or x.project} {x.year})'.replace(' )', ')'),
             condition=x.node_id, trait=x.trait, first_author=x.first_author, year=x.year, pmid=x.pmid,
             n_samples=x.n_samples, ancestry=x.ancestry, eas=x.eas, n_credible_sets=x.n_credible_sets)
        edge('gwas_of', x.study_id, alias.get(x.node_id, x.node_id))

    # ── checks, tables, JSON, SQLite ────────────────────────────────────────
    nodes = {k: pd.DataFrame(v, columns=['id', 'type', 'label', 'url'] + NODES[k][2]) for k, v in N.items()}
    edges = {k: pd.DataFrame(v, columns=['source', 'target', 'type'] + EDGES[k][3]) for k, v in E.items()}
    types = {i: k for k, df in nodes.items() for i in df.id}
    for k, df in edges.items():
        s_t, t_t = EDGES[k][0].split('|'), EDGES[k][1].split('|')
        bad = df[~df.source.map(types).isin(s_t) | ~df.target.map(types).isin(t_t)]
        assert not len(bad), (k, bad.head(3).to_dict('records'))
    for k, df in nodes.items():
        write_tsv(df, out / f'nodes_{k}.tsv')
    for k, df in edges.items():
        write_tsv(df, out / f'edges_{k}.tsv')
    allN = pd.concat([df[['id', 'type', 'label', 'url']] for df in nodes.values()], ignore_index=True)
    allE = pd.concat([df[['source', 'target', 'type']] for df in edges.values()], ignore_index=True)
    write_tsv(allN, out / 'nodes.tsv')
    write_tsv(allE, out / 'edges.tsv')
    counts = {'nodes': {k: len(v) for k, v in nodes.items()}, 'edges': {k: len(v) for k, v in edges.items()}}
    schema = {'nodes': {k: {'id': v[0], 'meaning': v[1], 'columns': ['id', 'type', 'label', 'url'] + v[2],
                            'count': counts['nodes'][k]} for k, v in NODES.items()},
              'edges': {k: {'source': v[0], 'target': v[1], 'meaning': v[2],
                            'columns': ['source', 'target', 'type'] + v[3], 'count': counts['edges'][k]}
                        for k, v in EDGES.items()},
              'ladder': {str(k): {'code': c, 'label': lab} for k, (c, lab) in LADDER.items()},
              'codes': CODES, 'about': {'default': P.about, 'within': P.within},
              'region_labels': {'p': 'P', 'p_pe': 'P in the secondary analysis', 'v': 'V', 'or': 'OR', 'maf': 'MAF',
                                'tier': 'tier', **((P.gwas or {}).get('columns') or {})}}
    meta = {'snapshot': tidy['fetch_date'], 'release': tidy['release'], 'core': CORE, 'core_name': core_name,
            'short': short, 'key': P.key, 'n_nodes': len(allN), 'n_edges': len(allE), 'counts': counts}

    def recs(df):
        return json.loads(df.to_json(orient='records', force_ascii=False))
    g = {'directed': True, 'multigraph': True, 'graph': meta,
         'nodes': [x for df in nodes.values() for x in recs(df)],
         'edges': [{k: v for k, v in x.items() if k != 'quote'} for df in edges.values() for x in recs(df)]}
    (out / 'graph.json').write_text(json.dumps(g, ensure_ascii=False))
    (out / 'schema.json').write_text(json.dumps({**schema, 'graph': meta}, indent=1, ensure_ascii=False))
    db = out / f'{P.stem}.sqlite'
    db.unlink(missing_ok=True)
    con = sqlite3.connect(db)
    for k, df in nodes.items():
        typed(df).to_sql(f'node_{k}', con, index=False)
    for k, df in edges.items():
        typed(df).to_sql(f'edge_{k}', con, index=False)
    allN.to_sql('nodes', con, index=False)
    allE.to_sql('edges', con, index=False)
    tables = {}
    # the publications the curation left out are counted, never listed: their rows leave every table here (they
    # stay in 03.curated, the record of the curation)
    excluded = set(pubs[pubs.status == 'excluded'].pmid)

    def listed(df):
        if 'pmid' in df:
            df = df[~df.pmid.isin(excluded)]
        if 'item' in df:
            df = df[[not (set(re.findall(r'PMID (\d+)', str(x))) & excluded) for x in df.item]]
        return df
    for d, prefix in ((Q, 'curated_'), (T, 'table_')):
        for f in sorted(d.glob('*.tsv')):
            typed(listed(read_tsv(f))).to_sql(prefix + f.stem, con, index=False)
            tables[prefix + f.stem] = f'{d.name}/{f.name}'
    pd.DataFrame([{'key': 'schema', 'value': json.dumps(schema)}, {'key': 'graph', 'value': json.dumps(meta)},
                  {'key': 'tables', 'value': json.dumps(tables)}]).to_sql('meta', con, index=False)
    for sql in ('CREATE INDEX i_nodes ON nodes(id)', 'CREATE INDEX i_label ON nodes(label)',
                'CREATE INDEX i_esrc ON edges(source)', 'CREATE INDEX i_etgt ON edges(target)'):
        con.execute(sql)
    con.commit()
    con.close()
    example = P.checks.get('example_gene') or (lad.symbol.iloc[0] if len(lad) else 'BRCA1')
    (out / 'README.md').write_text(readme(schema, meta, tables, P, example))
    log(f'graph: {len(allN)} nodes ({counts["nodes"]}), {len(allE)} edges ({counts["edges"]}); SQLite {len(tables)} '
        'tables besides the graph')


def readme(schema, meta, tables, P, example):
    s = P.short
    md = [f'# {s} knowledge graph — graph files (snapshot {meta["snapshot"]}, Open Targets {meta["release"]})', '',
          'Generated by DiseaseLens (build_graph.py); do not edit by hand. '
          f'{meta["n_nodes"]:,} nodes and {meta["n_edges"]:,} edges, typed as below. Every node and edge comes from '
          'the curated literature (`../03.curated`) and the long tables (`../04.tables`). In the SQLite file a '
          'column whose values are all numbers is numeric (an empty value is NULL); true / false are the text '
          '"True" / "False". The publications the curation left out (retracted, conference material) are counted '
          'and listed nowhere here; `../03.curated` keeps the record of the curation.', '',
          '| File | Content |', '|---|---|',
          '| `nodes.tsv`, `edges.tsv` | every node (id, type, label, url) and edge (source, target, type) |',
          '| `nodes_<type>.tsv`, `edges_<type>.tsv` | one table per type, with the columns below |',
          '| `graph.json` | node-link JSON (`networkx.node_link_graph(json.load(f), edges="edges")`); edges without '
          'their source quotes |',
          '| `schema.json` | the types below, their counts and the evidence ladder |',
          f'| `{P.stem}.sqlite` | the node and edge tables (`node_<type>`, `edge_<type>`, `nodes`, '
          '`edges`), every curated table (`curated_<name>`) and long table (`table_<name>`), and `meta` |', '',
          f'Query from the shell: `python3 scripts/kg.py --disease {P.key}` (gene, condition, region, publication, '
          'studies, search, neighbors, sql); see the repository README.', '', '## Node types', '',
          '| Type | Identifier | Meaning | Count | Columns |', '|---|---|---|---|---|']
    for k, v in schema['nodes'].items():
        md.append(f'| `{k}` | {v["id"]} | {v["meaning"]} | {v["count"]:,} | {", ".join(v["columns"][4:]) or "—"} |')
    md += ['', '## Edge types', '', '| Type | From → to | Meaning | Count | Columns |', '|---|---|---|---|---|']
    for k, v in schema['edges'].items():
        md.append(f'| `{k}` | {v["source"]} → {v["target"].replace("|", " or ")} | {v["meaning"]} | '
                  f'{v["count"]:,} | '
                  f'{", ".join(v["columns"][3:]) or "—"} |')
    md += ['', f'## Evidence ladder of a gene of {s} (`node_gene.disease_level`)', '',
           '| Level | Code | Evidence |', '|---|---|---|']
    md += [f'| {k} | {v["code"] or "—"} | {v["label"]} |' for k, v in schema['ladder'].items()]
    md += ['', f'Levels 3–6 count only findings about {s} itself (`edge_reports.about` = "{P.about}"); an association '
           f'with a phenotype or outcome within {s} is kept on its edge but does not raise the level.', '',
           '## Example queries (SQLite)', '', '```sql',
           f'-- the genes of {s} by level, with their function and how many related conditions share them',
           "SELECT label, disease_level, disease_code, theme, n_related_conditions FROM node_gene",
           "  WHERE disease_level IS NOT NULL ORDER BY disease_level, n_related_conditions DESC;", '',
           f'-- every finding on {example}, with its study, place and source sentence',
           'SELECT p.label, p.place, r.result, r.about, r.variant_name, r.quote FROM edge_reports r',
           f"  JOIN node_publication p ON p.id = r.source JOIN node_gene g ON g.id = r.target WHERE g.label = '{example}';",
           '', f'-- conditions sharing two or more GWAS regions of {s}',
           'SELECT c.label, COUNT(*) AS regions, GROUP_CONCAT(REPLACE(s.target, \'region:\', \'\')) FROM edge_shares_region s',
           '  JOIN node_condition c ON c.id = s.source GROUP BY c.id HAVING regions >= 2 ORDER BY regions DESC;', '',
           f'-- genetic studies of {s} done in Japan',
           "SELECT label, design_class, n_cases, title FROM node_publication WHERE japanese = 'True' ORDER BY year;",
           '```', '', '## Tables in the SQLite file', '', '| Table | Source |', '|---|---|']
    md += [f'| `{k}` | `../{v}` |' for k, v in tables.items()]
    return '\n'.join(md) + '\n'


if __name__ == '__main__':
    main()
