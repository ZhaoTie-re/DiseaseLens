#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The long tables of one snapshot of one disease, every count asserted
#           against its source. The literature tables come curated (curate.py:
#           excluded publications left out, variants named, findings classified
#           by what they are about, text-mining mismatches dropped); the
#           publications the curation leaves out are also removed from the GWAS
#           studies and from the references of the curated gene evidence.
#             conditions.tsv         the conditions, why they relate to the disease, their
#                                    genetic data (GWAS studies, credible sets, East
#                                    Asian studies, Mendelian genes)
#             relations.tsv          one row per piece of relation evidence
#             colocalisation.tsv     disease region x colocalizing GWAS credible set
#             gwas_studies.tsv       every GWAS study of the disease and the conditions:
#                                    design (n, cases, controls, cohorts), ancestry,
#                                    East Asian flag, summary statistics
#             gwas_loci.tsv          every credible set: lead variant, its own P / beta,
#                                    fine-mapping, top three locus-to-gene genes
#             mendelian.tsv          condition x gene x source (the disease itself
#                                    included): score, classifications (ClinGen /
#                                    Gene2Phenotype class, PanelApp rating, ...), ClinVar
#                                    review status and clinical significance, specific
#                                    diseases, PMIDs
#             mendelian_records.tsv  every curated / Mendelian / burden evidence row
#             disease_targets.tsv    the disease's associated targets, their literature
#                                    score, curated text-mined publications, Reactome themes
#             literature.tsv         text-mined target x publication with sentences (curated)
#             disease_studies / findings / methods / screening .tsv
#                                    the genetic studies of the disease (curated)
#             overlap.tsv            published region x GWAS credible set containing its
#                                    lead variant (regions without a credible set of the
#                                    disease)
#             disease_regions.tsv    the GWAS regions of the disease (published, or from
#                                    Open Targets) and their status in Open Targets
#             disease_genes.tsv      the evidence ladder of the disease's genes (level,
#                                    code, function, GWAS region, curated sources,
#                                    literature tallies)
#             gene_condition.tsv     gene x condition: the kinds of genetic evidence; for
#                                    the disease, GWAS = a gene of a genome-wide significant
#                                    region (the secondary label when significant only in a
#                                    secondary analysis), Mendelian = curated evidence,
#                                    Literature = a positive finding about the disease in a
#                                    genetic study
#             ladder.json            the evidence ladder of this snapshot: levels, codes
#             themes.json            the gene themes used (curated, or Reactome; the
#                                    automatic ones also as themes_auto.yaml, a draft)
#             tidy.json              counts and release
#           The graph is built from these tables by build_graph.py.
# Inputs  : the fetch, resolve, extraction and curation outputs, the profile
# Outputs : as above, in --out-dir
# ---------------------------------------------------------------------------
import argparse
import json
import math
import shutil
import sys
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import (EAS_COUNTRIES, NOT_REGION_COLS, SHARED_SIGNAL, SHARED_SIGNAL_OLD, auto_themes, check_themes,  # noqa: E402
                       curated_links, disease_regions, gene_ladder, ladder_spec, mining_links, p_unicode, relates,
                       level_status, themes_with_fallback,
                       legacy_path, load_json, load_profile, log, num, read_tsv, shared_relations, us_english,
                       write_json, write_tsv)

SOURCE_LABEL = {'clingen': 'ClinGen', 'gene2phenotype': 'Gene2Phenotype', 'genomics_england': 'PanelApp',
                'orphanet': 'Orphanet', 'eva': 'ClinVar', 'uniprot_literature': 'UniProt literature',
                'uniprot_variants': 'UniProt variants', 'gene_burden': 'Gene burden'}
STUDY_COLS = ['node_id', 'condition', 'study_id', 'project', 'trait', 'first_author', 'year', 'journal', 'pmid',
              'ancestry', 'eas', 'n_samples', 'n_cases', 'n_controls', 'cohorts', 'summary_statistics',
              'n_credible_sets']
LOCUS_COLS = ['node_id', 'condition', 'study_id', 'project', 'trait', 'first_author', 'year', 'journal', 'pmid',
              'ancestry', 'eas', 'study_locus_id', 'variant_id', 'rsid', 'p_mantissa', 'p_exponent', 'log10_p', 'beta',
              'finemapping', 'confidence'] + [f'l2g{k}_{x}' for k in (1, 2, 3) for x in ('symbol', 'ensembl', 'score')]
MREC_COLS = ['node_id', 'condition', 'source', 'ensembl_id', 'symbol', 'specific_disease', 'disease_from_source',
             'score', 'confidence', 'review_status', 'clinical_significance', 'allele_origin', 'variant', 'pmids',
             'cohort', 'statistical_method', 'p_mantissa', 'p_exponent', 'odds_ratio', 'ancestry', 'sample_size']
M_COLS = ['node_id', 'condition', 'source', 'ensembl_id', 'symbol', 'source_score', 'n_records', 'classifications',
          'review_status', 'clinical_significance', 'specific_diseases', 'pmids']
GC_COLS = ['ensembl_id', 'symbol', 'node_id', 'condition', 'evidence', 'detail', 'eas']
REGION_COLS = ['gene', 'genes', 'rsid', 'chrom', 'pos_grch37', 'pos_grch38', 'maf', 'or', 'p', 'p_mantissa',
               'p_exponent', 'p_pe', 'v', 'tier',
               'evidence', 'code', 'source', 'variant_id', 'lookup', 'disease_study_locus_id', 'l2g_gene', 'l2g_score',
               'n_overlap_sets']


def ancestry(pops):
    return '; '.join(f'{p["ldPopulation"]} {100 * p["relativeSampleSize"]:.0f}%'
                     for p in sorted(pops or [], key=lambda p: -p['relativeSampleSize']))


def eas_share(pops):
    return sum(p['relativeSampleSize'] for p in pops or [] if p['ldPopulation'] == 'eas')


def log10p(m, e):
    return math.log10(m) + e if m and m > 0 and e is not None else None


def study_rows(node, name, studies, eas_min):
    srows, lrows = [], []
    for s in studies:
        pops = s.get('ldPopulationStructure') or []
        base = {'node_id': node, 'condition': name, 'study_id': s['id'], 'project': s.get('projectId') or '',
                'trait': s.get('traitFromSource') or '', 'first_author': s.get('publicationFirstAuthor') or '',
                'year': (s.get('publicationDate') or '')[:4], 'journal': s.get('publicationJournal') or '',
                'pmid': s.get('pubmedId') or '', 'ancestry': ancestry(pops),
                'eas': eas_share(pops) >= eas_min}
        srows.append({**base, 'n_samples': s.get('nSamples'), 'n_cases': s.get('nCases'),
                      'n_controls': s.get('nControls'), 'cohorts': '; '.join(s.get('cohorts') or []),
                      'summary_statistics': bool(s.get('hasSumstats')),
                      'n_credible_sets': len(s['credible_sets'])})
        for c in s['credible_sets']:
            v = c.get('variant') or {}
            l2g = (c.get('l2GPredictions') or {}).get('rows') or []
            row = {**base, 'study_locus_id': c['studyLocusId'], 'variant_id': v.get('id', ''),
                   'rsid': ' '.join(v.get('rsIds') or []), 'p_mantissa': c.get('pValueMantissa'),
                   'p_exponent': c.get('pValueExponent'),
                   'log10_p': log10p(c.get('pValueMantissa'), c.get('pValueExponent')), 'beta': c.get('beta'),
                   'finemapping': c.get('finemappingMethod') or '', 'confidence': c.get('confidence') or ''}
            for k in range(3):
                row[f'l2g{k + 1}_symbol'] = l2g[k]['target']['approvedSymbol'] if k < len(l2g) else ''
                row[f'l2g{k + 1}_ensembl'] = l2g[k]['target']['id'] if k < len(l2g) else ''
                row[f'l2g{k + 1}_score'] = l2g[k]['score'] if k < len(l2g) else None
            lrows.append(row)
    return srows, lrows


def curated_rows(node, name, d):
    """The curated / Mendelian evidence of one disease file: (records, links)."""
    rec, links = [], []
    for e in d['mendelian_evidence']:
        src = SOURCE_LABEL.get(e['datasourceId'], e['datasourceId'])
        conf = e.get('confidence') or ''   # ClinVar: the review status, not a classification
        rec.append({'node_id': node, 'condition': name, 'source': src,
                    'ensembl_id': e['target']['id'], 'symbol': e['target']['approvedSymbol'],
                    'specific_disease': (e.get('disease') or {}).get('name', ''),
                    'disease_from_source': e.get('diseaseFromSource') or '', 'score': e['score'],
                    'confidence': '' if src == 'ClinVar' else conf,
                    'review_status': conf if src == 'ClinVar' else '',
                    'clinical_significance': '; '.join(e.get('clinicalSignificances') or []),
                    'allele_origin': '; '.join(e.get('alleleOrigins') or []),
                    'variant': ((e.get('variant') or {}).get('rsIds') or [''])[0] or e.get('variantRsId') or
                    (e.get('variant') or {}).get('id', ''),
                    'pmids': ' '.join(e.get('literature') or []), 'cohort': e.get('cohortShortName') or '',
                    'statistical_method': e.get('statisticalMethod') or '',
                    'p_mantissa': e.get('pValueMantissa'), 'p_exponent': e.get('pValueExponent'),
                    'odds_ratio': e.get('oddsRatio'), 'ancestry': e.get('ancestry') or '',
                    'sample_size': e.get('studySampleSize')})
    for ds, ts in d['mendelian_targets'].items():
        for t in ts:
            links.append({'node_id': node, 'condition': name, 'source': SOURCE_LABEL.get(ds, ds),
                          'ensembl_id': t['id'], 'symbol': t['symbol'], 'source_score': t['score']})
    return rec, links


def region_detail(r):
    """The evidence text of a region gene: lead variant, region, P (and V and tier where published). P from Open
    Targets is written from its mantissa and exponent (a float cannot hold P below 1e-308)."""
    if r['source'] == 'published':
        parts = [f'P {r["p"]:.1e}'] + ([f'V {r["v"]:.1e}'] if r.get('v') is not None else []) + \
            ([f'tier {r["tier"]}'] if r.get('tier') is not None else [])
        return f'{r["rsid"]} ({r["gene"]} region; {", ".join(parts)})'
    return (f'{r["rsid"] or r["variant_id"]} ({r["study_id"]}; P {p_unicode(r["p_mantissa"], r["p_exponent"])}, '
            f'L2G {float(r["l2g_score"]):.2f})')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--profile', required=True, help='profile.json of the disease')
    ap.add_argument('--raw', required=True, help='00.raw of the disease')
    ap.add_argument('--conditions', required=True, help='resolve.py outputs')
    ap.add_argument('--diseases', nargs='*', default=[], help='disease.<id>.json.gz')
    ap.add_argument('--extract', required=True, help='extract_genetics.py outputs (for their counts)')
    ap.add_argument('--curated', required=True, help='curate.py outputs')
    ap.add_argument('--eas-min', type=float, default=0.5, help='East Asian share of a study to flag it')
    ap.add_argument('--region-window-kb', type=int, default=500, help='lead variants this close are one region')
    ap.add_argument('--genome-wide-p', type=float, default=5e-8, help='P of a GWAS region of the disease')
    ap.add_argument('--region-min-maf', type=float, default=0.01,
                    help='a lead of a region from Open Targets rarer than this needs P below --rare-variant-p')
    ap.add_argument('--region-min-mac', type=float, default=25,
                    help='expected copies of the minor allele of such a lead in the smaller group of its study')
    ap.add_argument('--rare-variant-p', type=float, default=5e-9, help='P of a rare lead of such a region')
    ap.add_argument('--l2g-min', type=float, default=0.5, help='top locus-to-gene score for a gene-condition link')
    ap.add_argument('--min-shared-regions', type=int, default=2,
                    help='regions a condition related only by shared GWAS signals shares with the disease, at least')
    ap.add_argument('--out-dir', default='.')
    a = ap.parse_args()
    P = load_profile(a.profile)
    CORE = P.id
    raw, R, X, out = Path(a.raw), Path(a.conditions), Path(a.extract), Path(a.out_dir)
    Q = Path(a.curated)
    pubs = read_tsv(Q / 'publications.tsv')
    excluded = set(pubs[pubs.status == 'excluded'].pmid)
    curation = json.loads((Q / 'curation.json').read_text())
    out.mkdir(parents=True, exist_ok=True)
    meta, fc = load_json(raw / 'meta.json'), load_json(legacy_path(raw, 'fetch_core.json', 'fetch_cteph.json'))
    assoc, ev = load_json(raw / 'associations.json.gz'), load_json(raw / 'evidence.json.gz')
    tg, cg = load_json(raw / 'targets.json.gz'), load_json(raw / 'gwas.json.gz')
    cond = read_tsv(R / 'conditions.tsv')
    rel = read_tsv(Q / 'relations.tsv') if (Q / 'relations.tsv').exists() else read_tsv(R / 'relations.tsv')
    co, ov = read_tsv(R / 'colocalisation.tsv'), read_tsv(R / 'overlap.tsv')
    # the saved sources in American English (names, evidence types, statements; quotes stay verbatim)
    for c in ('name', 'sources', 'category'):
        cond[c] = cond[c].map(us_english)
    cond['category'] = cond.category.replace({SHARED_SIGNAL_OLD: SHARED_SIGNAL})   # a saved snapshot's old name
    for c in ('evidence_type', 'relation', 'detail', 'design', 'population'):
        rel[c] = rel[c].map(us_english)
    rel['relation'] = rel.relation.str.replace('shared genetic signal (', 'shared GWAS signal (', regex=False)  # old
    co['condition'], ov['condition'] = co.condition.map(us_english), ov.condition.map(us_english)
    alias = {i: n for n, ids in zip(cond.node_id, cond.ot_ids) for i in str(ids).split()}

    # publications left out by the curation: a GWAS study they report defines no region of the disease and no shared
    # signal (the clinical relations they support were removed by the curation); the shared-signal relations are
    # then rebuilt from what remains
    gone_core = [s for s in cg['studies'] if str(s.get('pubmedId') or '') in excluded]
    if gone_core:
        cg = {**cg, 'studies': [s for s in cg['studies'] if s not in gone_core]}
    gone_sets = {c['studyLocusId'] for s in gone_core for c in s['credible_sets']}
    n_co, n_ov = len(co), len(ov)
    co = co[~co.pmid.isin(excluded) & ~co.disease_study_locus_id.isin(gone_sets)].copy()
    ov = ov[~ov.pmid.isin(excluded)].copy()
    if len(co) < n_co or len(ov) < n_ov or gone_core:
        log(f'excluded publications: {len(gone_core)} GWAS studies of {P.short}, {n_co - len(co)} colocalizations, '
            f'{n_ov - len(ov)} lead-variant overlaps left out; shared-signal relations rebuilt')
        rows, _shared, _names = shared_relations(co, ov, alias, P.short)
        new = pd.DataFrame(rows, columns=rel.columns).fillna('').astype(str)
        for c in ('evidence_type', 'relation', 'detail'):
            new[c] = new[c].map(us_english)
        rel = pd.concat([rel[~rel.evidence_type.isin(['colocalization', 'lead-variant overlap'])], new],
                        ignore_index=True)

    # a condition is related only by evidence that relates it (not clinical evidence that found no association, and
    # not evidence the curation removed); a condition left without such evidence leaves every table here
    relating = pd.Series([relates(t_, r_) for t_, r_ in zip(rel.evidence_type, rel.relation)], index=rel.index,
                         dtype=bool)
    kept = set(rel[relating].node_id)
    # a condition related only by shared GWAS signals shares at least --min-shared-regions of the disease's regions:
    # one region alone, often a pleiotropic one, does not relate a trait by itself (its signal still counts as the
    # evidence of a condition related otherwise)
    shared_types = {'colocalization', 'colocalisation', 'lead-variant overlap'}
    by_node = rel[relating].groupby('node_id')
    types_of = by_node.evidence_type.agg(set)
    regions_of = by_node.loci.agg(lambda x: {g for v in x for g in str(v).split() if g})
    too_few = sorted(n for n in kept if types_of[n] <= shared_types and len(regions_of[n]) < a.min_shared_regions)
    # a broad class (resolve.py: more descendants in the disease ontology than --max-descendants) related only by
    # shared signals leaves too: its signals are those of its subtypes, each listed on its own when it shares them
    broad = set(cond[cond.get('broad_class', pd.Series('', index=cond.index)).astype(str).eq('True')].node_id)
    too_broad = sorted(n for n in kept if n in broad and types_of[n] <= shared_types and n not in too_few)
    names_tf = sorted(cond.set_index('node_id').name.reindex(too_few).fillna('')) if too_few else []
    if too_few:
        log(f'{len(too_few)} conditions related only by shared GWAS signals in fewer than {a.min_shared_regions} '
            f'{P.short} regions left out')
        kept -= set(too_few)
    names_tb = sorted(cond.set_index('node_id').name.reindex(too_broad).fillna('')) if too_broad else []
    if too_broad:
        log(f'{len(too_broad)} broad classes related only by shared GWAS signals left out: {", ".join(names_tb)}')
        kept -= set(too_broad)
    unrelated = cond[~cond.node_id.isin(kept) & ~cond.node_id.isin(too_few) & ~cond.node_id.isin(too_broad)]
    if len(unrelated):
        log(f'not related (no evidence that relates them): {", ".join(unrelated.name)}')
    cond, rel = cond[cond.node_id.isin(kept)].copy(), rel[rel.node_id.isin(kept)].copy()
    relating = relating[rel.index]
    co = co[co.node_id.map(lambda n: alias.get(n, n)).isin(kept)].copy()
    ov = ov[ov.node_id.map(lambda n: alias.get(n, n)).isin(kept)].copy()
    # the sources, category and counts of each condition, from the relations that relate it (after the curation)
    types = rel[relating].groupby('node_id').evidence_type.agg(lambda x: ' '.join(sorted(set(x))))
    cond['sources'] = cond.node_id.map(types).fillna('')
    clinical_cats = {c['name'] for c in P.categories}
    cond['category'] = [c if 'clinical' in s.split() or c not in clinical_cats else
                        P.ontology_category if 'ontology' in s.split() else SHARED_SIGNAL
                        for c, s in zip(cond.category, cond.sources)]
    cond['n_clinical'] = cond.node_id.map(rel[(rel.evidence_type == 'clinical') & relating]
                                          .groupby('node_id').size()).fillna(0).astype(int).astype(str)
    _rows, shared, _names = shared_relations(co, ov, alias, P.short)
    for c, kind in (('coloc_loci', 'colocalization'), ('overlap_loci', 'lead-variant overlap')):
        cond[c] = cond.node_id.map(rel[rel.evidence_type == kind].set_index('node_id').n_loci).fillna('0')
    cond['shared_regions'] = cond.node_id.map(lambda n: ' '.join(sorted(shared.get(n, ()))))
    cond['n_shared_regions'] = cond.node_id.map(lambda n: str(len(shared.get(n, ()))))
    published = load_json(raw / 'published_loci.json.gz') if (raw / 'published_loci.json.gz').exists() else None
    assert len(assoc) == fc['associations'] and len(ev) == fc['evidence']
    core_name = load_json(raw / 'disease.json')['name']

    # ── the disease's targets, themes, text-mined literature ─────────────────
    pathways = {t['approvedSymbol']: sorted({p['topLevelTerm'] for p in t.get('pathways') or [] if p.get('topLevelTerm')})
                for t in tg}
    themes_of = {t['id']: pathways[t['approvedSymbol']] for t in tg}
    gp = read_tsv(X / 'gene_pathways.tsv', ['symbol', 'top_level_pathways'])
    for r in gp.itertuples():
        pathways.setdefault(r.symbol, [x for x in r.top_level_pathways.split('; ') if x])
    lit, review = read_tsv(Q / 'literature.tsv'), read_tsv(Q / 'textmining_review.tsv')
    assert len(review) == sum(len(e['literature'] or ['']) for e in mining_links(ev))
    lit['year'] = num(lit.year)
    dropped = review[review.decision == 'dropped'].groupby('ensembl_id').size()
    ln = lit.groupby('ensembl_id').agg(n_pmids=('pmid', 'nunique'), first_year=('year', 'min'),
                                       last_year=('year', 'max'))
    trows = []
    for k, r in enumerate(sorted(assoc, key=lambda r: (-r['score'], r['target']['approvedSymbol'])), 1):
        t, ds = r['target'], {x['id']: x['score'] for x in r['datasourceScores']}
        trows.append({'rank': k, 'ensembl_id': t['id'], 'symbol': t['approvedSymbol'], 'name': t.get('approvedName')
                      or '', 'score': r['score'], 'ds_europepmc': ds.get('europepmc'),
                      'ds_gwas_credible_sets': ds.get('gwas_credible_sets'),
                      'n_pmids': int(ln.n_pmids.get(t['id'], 0)),
                      'n_pmids_dropped': int(dropped.get(t['id'], 0)), 'first_year': ln.first_year.get(t['id']),
                      'last_year': ln.last_year.get(t['id']), 'themes': '; '.join(themes_of.get(t['id'], []))})
    targets = pd.DataFrame(trows, columns=['rank', 'ensembl_id', 'symbol', 'name', 'score', 'ds_europepmc',
                                           'ds_gwas_credible_sets', 'n_pmids', 'n_pmids_dropped', 'first_year',
                                           'last_year', 'themes'])

    # ── GWAS and curated evidence of the disease and of the conditions ───────
    sr, lr = study_rows(CORE, core_name, cg['studies'], a.eas_min)
    assert len(lr) + len(gone_sets) == fc['credible_sets']
    mrows, mrec = [], []
    core_men = load_json(raw / 'mendelian.json.gz') if (raw / 'mendelian.json.gz').exists() else None
    if core_men:
        drop = set(curation.get('core_curated_dropped') or [])
        core_men = {**core_men, 'mendelian_targets': {k: [t for t in v if t['symbol'] not in drop]
                                                      for k, v in core_men['mendelian_targets'].items()},
                    'mendelian_evidence': [e for e in core_men['mendelian_evidence']
                                           if e['target']['approvedSymbol'] not in drop]}
        r_, l_ = curated_rows(CORE, core_name, core_men)
        mrec, mrows = mrec + r_, mrows + l_
    names = dict(zip(cond.node_id, cond.name))
    fetched = set()
    for f in sorted(a.diseases):
        d = load_json(f)
        node = d['disease']['id']
        assert d['release']['data_version'] == meta['data_version'], f'{f}: another release'
        if node not in kept:                                  # fetched for a condition that is not related
            continue
        fetched.add(node)
        s1, l1 = study_rows(node, names.get(node, d['disease']['name']), d['studies'], a.eas_min)
        assert len(l1) == d['n_credible_sets']
        sr += s1
        lr += l1
        r_, l_ = curated_rows(node, names.get(node, d['disease']['name']), d)
        mrec, mrows = mrec + r_, mrows + l_
    want = set(cond[cond.genetics.eq('full') & cond.ot_ids.ne('')].node_id)
    assert fetched == want, (sorted(want - fetched), sorted(fetched - want))
    studies, loci = pd.DataFrame(sr, columns=STUDY_COLS), pd.DataFrame(lr, columns=LOCUS_COLS)
    # publications left out by the curation (retracted, conference material): their GWAS studies go, and they leave
    # the reference lists of the curated gene evidence
    gone = studies[studies.pmid.isin(excluded)]
    studies, loci = studies[~studies.pmid.isin(excluded)], loci[~loci.study_id.isin(set(gone.study_id))]
    mrec = pd.DataFrame(mrec, columns=MREC_COLS)
    refs_before = mrec.pmids.str.split().str.len().sum()
    mrec['pmids'] = mrec.pmids.map(lambda v: ' '.join(x for x in v.split() if x not in excluded))
    excluded_use = {'gwas_studies': len(gone) + len(gone_core), 'curated_evidence_references':
                    int(refs_before - mrec.pmids.str.split().str.len().sum())}
    M = pd.DataFrame(mrows)
    # ClinVar alone counts only with a reviewed pathogenic record; a gene needs an approved symbol (the same rule for
    # the disease and its related conditions)
    if len(M):
        M, mrec, curated_dropped = curated_links(M, mrec)
    else:
        curated_dropped = {'clinvar_review': 0, 'no_symbol': 0}
    if len(M):
        agg = mrec.groupby(['node_id', 'source', 'ensembl_id']).agg(
            n_records=('score', 'size'),
            classifications=('confidence', lambda x: '; '.join(sorted({v for v in x if v}))),
            review_status=('review_status', lambda x: '; '.join(sorted({v for v in x if v}))),
            clinical_significance=('clinical_significance',
                                   lambda x: '; '.join(sorted({s for v in x for s in v.split('; ') if s}))),
            specific_diseases=('specific_disease', lambda x: '; '.join(sorted({v for v in x if v}))),
            pmids=('pmids', lambda x: ' '.join(sorted({p for v in x for p in v.split() if p})))).reset_index()
        M = M.merge(agg, on=['node_id', 'source', 'ensembl_id'], how='left')
        M['n_records'] = M.n_records.fillna(0).astype(int)
        M = M.sort_values(['condition', 'symbol', 'source'])
    M = M.reindex(columns=M_COLS) if len(M) else pd.DataFrame(columns=M_COLS)

    # ── colocalization: East Asian flag of the colocalizing study ────────────
    co['eas'] = co.ancestry.map(lambda s: sum(float(x.split()[1].rstrip('%')) / 100 for x in s.split('; ')
                                              if x.startswith('eas ')) >= a.eas_min if s else False)

    # ── the East Asian evidence of a relation: the population of the evidence that relates the condition to the
    # disease, not of the condition's other studies. Clinical evidence: a population its source states, or its study
    # done in East Asia (below); a shared signal counts when both GWAS, the disease's and the condition's, have mostly
    # East Asian samples;
    # a lead-variant overlap has no credible set of the disease (its sample is not recorded): never; ontology: never.
    core_loci = loci[loci.node_id == CORE]
    core_eas = dict(zip(core_loci.study_locus_id, core_loci.eas.astype(bool)))
    co['eas_both'] = [bool(e) and core_eas.get(s, False) for e, s in zip(co.eas, co.disease_study_locus_id)]
    shared_eas = {alias.get(n, n) for n, both in zip(co.node_id, co.eas_both) if both}   # right for an empty table
    # clinical evidence is East Asian when its source states an East Asian population, or when its study was done in
    # East Asia (the country of the first author's institution): the rule of the genetic studies of the disease
    eas_place = {p for p, c in zip(pubs.pmid, pubs.get('institution_country', pd.Series(dtype=str)))
                 if c in EAS_COUNTRIES}
    rel['east_asian'] = [str(e == 'True' or p in eas_place) if t == 'clinical' else
                         str(t == 'colocalization' and n in shared_eas)
                         for t, e, p, n in zip(rel.evidence_type, rel.east_asian.astype(str), rel.pmid, rel.node_id)]

    # ── genetic studies of the disease (extracted, then curated) ─────────────
    cst, cfi = read_tsv(Q / 'studies.tsv', ['pmid', 'japanese']), read_tsv(Q / 'findings.tsv', [
        'pmid', 'gene', 'symbol', 'variant', 'result', 'ensembl_id', 'association_with', 'variant_name'])
    cmf, csc = read_tsv(Q / 'methods.tsv'), read_tsv(Q / 'screening.tsv')

    # ── the GWAS regions of the disease ──────────────────────────────────────
    not_credible = []
    regions = disease_regions(cg, published, P, a.l2g_min, a.region_window_kb, a.genome_wide_p, a.region_min_maf,
                              not_credible, a.region_min_mac, a.rare_variant_p)
    if not_credible:
        log(f'not regions (lead not credible, not replicated): {len(not_credible)} credible sets of '
            + ', '.join(sorted({x["study_id"] for x in not_credible})))
    write_tsv(pd.DataFrame(not_credible, columns=NOT_REGION_COLS), out / 'disease_signals_not_regions.tsv')

    # ── gene x condition evidence ────────────────────────────────────────────
    gc = []
    for r in loci.itertuples():
        if r.node_id == CORE:          # the disease's GWAS genes are the genes of its regions (below)
            continue
        if r.l2g1_symbol and num(pd.Series([r.l2g1_score]))[0] >= a.l2g_min:
            gc.append({'ensembl_id': r.l2g1_ensembl, 'symbol': r.l2g1_symbol, 'node_id': r.node_id,
                       'condition': r.condition, 'evidence': 'GWAS', 'detail': f'{r.rsid or r.variant_id} '
                       f'({r.study_id}, L2G {float(r.l2g1_score):.2f})', 'eas': r.eas})
    for r in M.itertuples():
        gc.append({'ensembl_id': r.ensembl_id, 'symbol': r.symbol, 'node_id': r.node_id, 'condition': r.condition,
                   'evidence': 'Mendelian', 'detail': f'{r.source} ({r.classifications or r.clinical_significance or "score " + format(float(r.source_score), ".2f")})',  # noqa: E501
                   'eas': False})
    for r in regions:
        for g in r['genes']:
            gc.append({'ensembl_id': r['ensembl'][g], 'symbol': g, 'node_id': CORE, 'condition': core_name,
                       'evidence': 'GWAS' if r['code'] == 'G' else P.secondary_label,
                       'detail': region_detail(r), 'eas': core_eas.get(r.get('disease_study_locus_id') or '', False)})
    eas_study = dict(zip(cst.pmid, cst.east_asian)) if len(cst) else {}
    for r in cfi.itertuples():
        if r.ensembl_id and r.association_with == P.about and r.result in ('associated', 'identified', 'enriched'):
            gc.append({'ensembl_id': r.ensembl_id, 'symbol': r.symbol, 'node_id': CORE, 'condition': core_name,
                       'evidence': 'Literature', 'detail': f'PMID {r.pmid}: {r.result} {r.variant_name}'.strip(),
                       'eas': eas_study.get(r.pmid) == 'True'})
    GC = pd.DataFrame(gc, columns=GC_COLS)
    GCs = (GC.groupby(['ensembl_id', 'symbol', 'node_id', 'condition', 'evidence'])
           .agg(n=('detail', 'size'), details=('detail', lambda x: '; '.join(list(dict.fromkeys(x))[:5])),
                eas=('eas', 'any')).reset_index())

    # ── the ladder, its themes and its definition ────────────────────────────
    cm = M[M.node_id == CORE]
    curated = {s: {'ensembl_id': g.ensembl_id.iloc[0], 'sources': sorted(set(g.source))}
               for s, g in cm.groupby('symbol')} if len(cm) else {}
    # the curated themes of the profile, with each gene they do not list in its Reactome top-level pathway; without
    # curated themes, Reactome for every gene (written as a draft to curate)
    genes = sorted({s for s in cfi.symbol if s} | {g for r in regions for g in r['genes']} | set(curated))
    if P.themes:
        themes = themes_with_fallback(check_themes(P.themes), genes, pathways)
    else:
        themes = auto_themes(genes, pathways)
        (out / 'themes_auto.yaml').write_text(yaml.safe_dump({'themes': themes}, sort_keys=False, allow_unicode=True,
                                                             width=118))
    spec = ladder_spec(P, with_curated=bool(curated))
    ladder = gene_ladder(cfi, cst, themes, regions, P.about, P.within, curated, spec)

    # ── write ────────────────────────────────────────────────────────────────
    cond['n_gwas_studies'] = cond.node_id.map(studies.groupby('node_id').size()).fillna(0).astype(int)
    cond['n_credible_sets'] = cond.node_id.map(loci.groupby('node_id').size()).fillna(0).astype(int)
    cond['n_eas_studies'] = cond.node_id.map(studies[studies.eas.astype(bool)].groupby('node_id').size()
                                             ).fillna(0).astype(int)
    cond['n_mendelian_genes'] = cond.node_id.map(M.groupby('node_id').ensembl_id.nunique() if len(M)
                                                 else pd.Series(dtype=int)).fillna(0).astype(int)
    # a condition whose genetics were not fetched (none: no Open Targets entry; coloc_only: its shared signals only)
    # has no count of its own: blank, not 0 (genetics_note says why)
    for c in ('n_gwas_studies', 'n_credible_sets', 'n_eas_studies', 'n_mendelian_genes'):
        cond[c] = [str(v) if g == 'full' else '' for v, g in zip(cond[c], cond.genetics)]
    cond['genetics_note'] = [n or ('no Open Targets entry' if g == 'none' else 'related only by a shared GWAS signal'
                                   if g == 'coloc_only' else '') for n, g in zip(cond.genetics_note, cond.genetics)]
    # the credible sets of other GWAS that hold the lead variant of a region without a credible set of the disease
    # (PIP >= ClppMin; overlap.tsv)
    n_ov_sets = ov.groupby('disease_gene').other_study_locus_id.nunique() if len(ov) else pd.Series(dtype=int)
    reg = pd.DataFrame([{**{k: r.get(k) for k in REGION_COLS if k not in ('genes', 'n_overlap_sets')},
                         'genes': ' '.join(r['genes']), 'n_overlap_sets': int(n_ov_sets.get(r['gene'], 0))}
                        for r in regions], columns=REGION_COLS)
    for name, df in [('disease_regions', reg), ('disease_genes', ladder), ('overlap', ov),
                     ('conditions', cond), ('relations', rel), ('colocalisation', co), ('gwas_studies', studies),
                     ('gwas_loci', loci), ('mendelian', M), ('mendelian_records', mrec), ('disease_targets', targets),
                     ('literature', lit), ('gene_condition', GCs)]:
        write_tsv(df, out / f'{name}.tsv')
    for f in ('studies', 'findings', 'methods', 'screening'):
        shutil.copy(Q / f'{f}.tsv', out / f'disease_{f}.tsv')
    spec['status'] = {str(k): dict(zip(('brief', 'full'), level_status(k, P, ladder))) for k in range(1, 8)}
    write_json(spec, out / 'ladder.json')
    write_json({'source': 'curated, and Reactome for the genes they do not list' if P.themes else 'Reactome (automatic)',
                'themes': themes}, out / 'themes.json')
    stats = {'release': meta['data_version'], 'api_version': meta['api_version'], 'fetch_date': meta['fetch_date'],
             # other Open Targets entries whose evidence the disease takes (its profile's evidence_also)
             'evidence_also': fc.get('evidence_also') or [],
             # conditions whose clinical evidence all found no association: not related (seed, or here)
             # genome-wide significant credible sets of the disease whose lead is not credible and not replicated
             # (disease_signals_not_regions.tsv)
             'signals_not_regions': len(not_credible),
             'not_related': sorted(set(read_tsv(R / 'not_related.tsv', ['name']).name if (R / 'not_related.tsv').exists()
                                       else []) | set(unrelated.name)),
             # traits related only by shared GWAS signals in fewer regions than the minimum: left out, named
             'shared_signal_too_few_regions': {'min_regions': a.min_shared_regions, 'n': len(too_few),
                                               'names': names_tf},
             # broad classes (more ontology descendants than the limit) related only by shared signals: left out
             'shared_signal_broad_classes': {'n': len(too_broad), 'names': names_tb},
             'fetched_at': meta['fetched_at'],
             'core': {'id': CORE, 'name': core_name, 'short': P.short, 'key': P.key},
             'conditions': len(cond),
             'conditions_by_source': {s: int(cond.sources.str.contains(s).sum())
                                      for s in ('clinical', 'ontology', 'colocalization', 'lead-variant overlap')},
             'relations': len(rel), 'colocalizations': len(co), 'overlaps': len(ov),
             # the GWAS of the related conditions (a study of two conditions counted once), and the disease's own
             'gwas_studies': int(studies[studies.node_id != CORE].study_id.nunique()),
             'disease_gwas_studies': int(studies[studies.node_id == CORE].study_id.nunique()),
             'regions': len(reg), 'published_regions': int(reg.source.eq('published').sum()),
             'published_in_ot': int((reg.source.eq('published') & reg.variant_id.fillna('').ne('')).sum()),
             'published_with_credible_set': int((reg.source.eq('published')
                                                 & reg.disease_study_locus_id.fillna('').ne('')).sum()),
             'regions_from_open_targets': int(reg.source.eq('Open Targets').sum()),
             'disease_genes_by_level': {int(k): int(v) for k, v in ladder.level.value_counts().sort_index().items()},
             'disease_genes_by_code': {k: int(v) for k, v in ladder.code.value_counts().items() if k},
             'curated_genes': len(curated), 'themes': 'curated, and Reactome for the genes they do not list' if P.themes else 'Reactome (automatic)',
             'credible_sets': int(loci[loci.node_id != CORE].study_locus_id.nunique()),
             'disease_credible_sets': int(loci[loci.node_id == CORE].study_locus_id.nunique()),
             'eas_studies': int(studies[(studies.node_id != CORE) & studies.eas.astype(bool)].study_id.nunique()),
             'disease_eas_studies': int(studies[(studies.node_id == CORE) & studies.eas.astype(bool)].study_id.nunique()),
             'mendelian_links': int((M.node_id != CORE).sum()) if len(M) else 0,
             'mendelian_records': len(mrec), 'disease_targets': len(targets), 'literature_rows': len(lit),
             'literature_pmids': int(lit.pmid.nunique()), 'disease_genetic_records': len(csc),
             'disease_genetic_studies': len(cst), 'disease_findings': len(cfi), 'methods_fields': len(cmf),
             'gene_condition': len(GCs), 'excluded_publications_removed': excluded_use,
             'curated_links_dropped': curated_dropped,
             'extract': json.loads((X / 'extract.json').read_text()), 'curation': curation,
             # for the comparison with an earlier snapshot only (diff_snapshots.py: counted, never listed)
             'excluded_pmids': sorted(excluded)}
    (out / 'tidy.json').write_text(json.dumps(stats, indent=1))
    log(f'{len(cond)} conditions, {len(studies)} GWAS studies, {len(loci)} credible sets, {len(M)} curated links, '
        f'{len(reg)} regions, {len(curated)} curated genes of {P.short}, {len(cst)} genetic studies (curated); '
        f'excluded publications removed: {excluded_use}; curated links left out: {curated_dropped}')


if __name__ == '__main__':
    main()
