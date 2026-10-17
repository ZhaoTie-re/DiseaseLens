#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The long tables of one snapshot, every count asserted against its
#           source. The literature tables come curated (curate.py: excluded
#           publications left out, variants named, findings classified by what
#           they are about, text-mining mismatches dropped); the publications
#           the curation leaves out are also removed from the GWAS studies and
#           from the references of the curated gene evidence.
#             conditions.tsv         the conditions, why they relate to CTEPH, their
#                                    genetic data (GWAS studies, credible sets, East
#                                    Asian studies, Mendelian genes)
#             relations.tsv          one row per piece of relation evidence
#             colocalisation.tsv     CTEPH locus x colocalizing GWAS credible set
#             gwas_studies.tsv       every GWAS study of CTEPH and the conditions:
#                                    design (n, cases, controls, cohorts), ancestry,
#                                    East Asian flag, summary statistics
#             gwas_loci.tsv          every credible set: lead variant, its own P / beta,
#                                    fine-mapping, top three locus-to-gene genes
#             mendelian.tsv          condition x gene x source: score, classifications
#                                    (ClinGen / Gene2Phenotype class, PanelApp rating, ...),
#                                    ClinVar review status and clinical significance,
#                                    specific diseases, PMIDs
#             mendelian_records.tsv  every curated / Mendelian / burden evidence row
#             cteph_targets.tsv      CTEPH's associated targets, their literature score,
#                                    curated text-mined publications and Reactome themes
#             literature.tsv         text-mined target x publication with sentences (curated)
#             cteph_studies / findings / methods / screening .tsv
#                                    the genetic studies of CTEPH (curated)
#             overlap.tsv            published CTEPH region x GWAS credible set containing
#                                    its lead variant (regions without a CTEPH credible set)
#             cteph_regions.tsv      the published CTEPH GWAS regions and their status in
#                                    Open Targets
#             cteph_genes.tsv        the evidence ladder of the CTEPH genes (level, function,
#                                    GWAS region, literature tallies)
#             gene_condition.tsv     gene x condition: the kinds of genetic evidence; for
#                                    CTEPH, GWAS = a published genome-wide significant region
#                                    (GWAS (PE-adjusted) when significant only after
#                                    leveraging PE), Literature = a positive finding about
#                                    CTEPH in a genetic study
#             tidy.json              counts and release
#           The graph is built from these tables by build_graph.py.
# Inputs  : the fetch, resolve, extraction and curation outputs
# Outputs : as above, in --out-dir
# ---------------------------------------------------------------------------
import argparse
import json
import math
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import (ABOUT_CTEPH, gene_ladder, load_json, load_themes, log, num, read_tsv, us_english,  # noqa: E402
                       write_tsv)

CORE = 'EFO_0022954'
SOURCE_LABEL = {'clingen': 'ClinGen', 'gene2phenotype': 'Gene2Phenotype', 'genomics_england': 'PanelApp',
                'orphanet': 'Orphanet', 'eva': 'ClinVar', 'uniprot_literature': 'UniProt literature',
                'uniprot_variants': 'UniProt variants', 'gene_burden': 'Gene burden'}


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


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cteph-raw', required=True)
    ap.add_argument('--conditions', required=True, help='resolve.py outputs')
    ap.add_argument('--diseases', nargs='*', default=[], help='disease.<id>.json.gz')
    ap.add_argument('--extract', required=True, help='extract_genetics.py outputs (for their counts)')
    ap.add_argument('--curated', required=True, help='curate.py outputs')
    ap.add_argument('--themes', required=True, help='gene themes of the CTEPH literature (YAML)')
    ap.add_argument('--eas-min', type=float, default=0.5, help='East Asian share of a study to flag it')
    ap.add_argument('--l2g-min', type=float, default=0.5, help='top locus-to-gene score for a gene-condition link')
    ap.add_argument('--out-dir', default='.')
    a = ap.parse_args()
    raw, R, X, out = Path(a.cteph_raw), Path(a.conditions), Path(a.extract), Path(a.out_dir)
    Q = Path(a.curated)
    pubs = read_tsv(Q / 'publications.tsv')
    excluded = set(pubs[pubs.status == 'excluded'].pmid)
    out.mkdir(parents=True, exist_ok=True)
    meta, fc = load_json(raw / 'meta.json'), load_json(raw / 'fetch_cteph.json')
    assoc, ev = load_json(raw / 'associations.json.gz'), load_json(raw / 'evidence.json.gz')
    tg, cg = load_json(raw / 'targets.json.gz'), load_json(raw / 'gwas.json.gz')
    cond, rel, co = read_tsv(R / 'conditions.tsv'), read_tsv(R / 'relations.tsv'), read_tsv(R / 'colocalisation.tsv')
    ov = read_tsv(R / 'overlap.tsv')
    # the saved sources in American English (names, evidence types, statements; quotes stay verbatim)
    for c in ('name', 'sources', 'category'):
        cond[c] = cond[c].map(us_english)
    for c in ('evidence_type', 'relation', 'detail', 'design', 'population'):
        rel[c] = rel[c].map(us_english)
    co['condition'], ov['condition'] = co.condition.map(us_english), ov.condition.map(us_english)
    pubr = load_json(raw / 'published_loci.json.gz')['regions']
    assert len(assoc) == fc['associations'] and len(ev) == fc['evidence']
    cteph_name = load_json(raw / 'disease.json')['name']

    # ── CTEPH targets, themes, text-mined literature ─────────────────────────
    themes = {t['id']: sorted({p['topLevelTerm'] for p in t.get('pathways') or [] if p.get('topLevelTerm')})
              for t in tg}
    lit, review = read_tsv(Q / 'literature.tsv'), read_tsv(Q / 'textmining_review.tsv')
    assert len(review) == sum(len(e['literature'] or ['']) for e in ev if e['datasourceId'] == 'europepmc')
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
                      'last_year': ln.last_year.get(t['id']), 'themes': '; '.join(themes.get(t['id'], []))})
    targets = pd.DataFrame(trows)

    # ── GWAS of CTEPH and of the conditions ─────────────────────────────────
    sr, lr = study_rows(CORE, cteph_name, cg['studies'], a.eas_min)
    assert len(lr) == fc['credible_sets']
    mrows, mrec = [], []
    names = dict(zip(cond.node_id, cond.name))
    fetched = set()
    for f in sorted(a.diseases):
        d = load_json(f)
        node = d['disease']['id']
        assert d['release']['data_version'] == meta['data_version'], f'{f}: another release'
        fetched.add(node)
        s1, l1 = study_rows(node, names.get(node, d['disease']['name']), d['studies'], a.eas_min)
        assert len(l1) == d['n_credible_sets']
        sr += s1
        lr += l1
        for e in d['mendelian_evidence']:
            src = SOURCE_LABEL.get(e['datasourceId'], e['datasourceId'])
            conf = e.get('confidence') or ''   # ClinVar: the review status, not a classification
            mrec.append({'node_id': node, 'condition': names.get(node, d['disease']['name']), 'source': src,
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
                mrows.append({'node_id': node, 'condition': names.get(node, d['disease']['name']),
                              'source': SOURCE_LABEL.get(ds, ds), 'ensembl_id': t['id'], 'symbol': t['symbol'],
                              'source_score': t['score']})
    want = set(cond[cond.genetics.eq('full') & cond.ot_ids.ne('')].node_id)
    assert fetched == want, (sorted(want - fetched), sorted(fetched - want))
    studies, loci = pd.DataFrame(sr), pd.DataFrame(lr)
    # publications left out by the curation (retracted, conference material): their GWAS studies go, and they leave
    # the reference lists of the curated gene evidence
    gone = studies[studies.pmid.isin(excluded)]
    studies, loci = studies[~studies.pmid.isin(excluded)], loci[~loci.study_id.isin(set(gone.study_id))]
    mrec = pd.DataFrame(mrec)
    refs_before = mrec.pmids.str.split().str.len().sum()
    mrec['pmids'] = mrec.pmids.map(lambda v: ' '.join(x for x in v.split() if x not in excluded))
    excluded_use = {'gwas_studies': len(gone), 'curated_evidence_references':
                    int(refs_before - mrec.pmids.str.split().str.len().sum())}
    M = pd.DataFrame(mrows)
    if len(M):
        agg = mrec.groupby(['node_id', 'source', 'ensembl_id']).agg(
            n_records=('score', 'size'),
            classifications=('confidence', lambda x: '; '.join(sorted({v for v in x if v}))),
            review_status=('review_status', lambda x: '; '.join(sorted({v for v in x if v}))),
            clinical_significance=('clinical_significance',
                                   lambda x: '; '.join(sorted({s for v in x for s in v.split('; ') if s}))),
            specific_diseases=('specific_disease', lambda x: '; '.join(sorted({v for v in x if v}))),
            pmids=('pmids', lambda x: ' '.join(sorted({p for v in x for p in v.split() if p}))),
            best_log10_p=('p_mantissa', lambda x: None)).reset_index()
        M = M.merge(agg, on=['node_id', 'source', 'ensembl_id'], how='left')
        M['n_records'] = M.n_records.fillna(0).astype(int)
        M = M.sort_values(['condition', 'symbol', 'source']).drop(columns='best_log10_p')

    # ── colocalization: East Asian flag of the colocalizing study ────────────
    co['eas'] = co.ancestry.map(lambda s: sum(float(x.split()[1].rstrip('%')) / 100 for x in s.split('; ')
                                              if x.startswith('eas ')) >= a.eas_min if s else False)

    # ── genetic studies of CTEPH (extracted, then curated) ──────────────────
    cst, cfi = read_tsv(Q / 'studies.tsv'), read_tsv(Q / 'findings.tsv')
    cmf, csc = read_tsv(Q / 'methods.tsv'), read_tsv(Q / 'screening.tsv')

    # ── gene x condition evidence ────────────────────────────────────────────
    gc = []
    for r in loci.itertuples():
        if r.node_id == CORE:          # CTEPH's GWAS genes are the published regions' genes (below)
            continue
        if r.l2g1_symbol and num(pd.Series([r.l2g1_score]))[0] >= a.l2g_min:
            gc.append({'ensembl_id': r.l2g1_ensembl, 'symbol': r.l2g1_symbol, 'node_id': r.node_id,
                       'condition': r.condition, 'evidence': 'GWAS', 'detail': f'{r.rsid or r.variant_id} '
                       f'({r.study_id}, L2G {float(r.l2g1_score):.2f})', 'eas': r.eas})
    for r in M.itertuples() if len(M) else []:
        gc.append({'ensembl_id': r.ensembl_id, 'symbol': r.symbol, 'node_id': r.node_id, 'condition': r.condition,
                   'evidence': 'Mendelian', 'detail': f'{r.source} ({r.classifications or r.clinical_significance or "score " + format(float(r.source_score), ".2f")})',
                   'eas': False})
    for r in pubr:
        for g in r['genes']:
            gc.append({'ensembl_id': r['ensembl'][g], 'symbol': g, 'node_id': CORE, 'condition': cteph_name,
                       'evidence': 'GWAS' if r['evidence'] == 'meta' else 'GWAS (PE-adjusted)',
                       'detail': f'{r["rsid"]} ({r["gene"]} region; P {r["p"]:.1e}, V {r["v"]:.1e}, tier '
                                 f'{r["tier"]})', 'eas': False})
    jp = dict(zip(cst.pmid, cst.japanese)) if len(cst) else {}
    for r in cfi.itertuples():
        if r.ensembl_id and r.association_with == ABOUT_CTEPH and r.result in ('associated', 'identified',
                                                                                'enriched'):
            gc.append({'ensembl_id': r.ensembl_id, 'symbol': r.symbol, 'node_id': CORE, 'condition': cteph_name,
                       'evidence': 'Literature', 'detail': f'PMID {r.pmid}: {r.result} {r.variant_name}'.strip(),
                       'eas': jp.get(r.pmid) == 'True'})
    GC = pd.DataFrame(gc)
    GCs = (GC.groupby(['ensembl_id', 'symbol', 'node_id', 'condition', 'evidence'])
           .agg(n=('detail', 'size'), details=('detail', lambda x: '; '.join(list(dict.fromkeys(x))[:5])),
                eas=('eas', 'any')).reset_index())

    # ── write ────────────────────────────────────────────────────────────────
    cond['n_gwas_studies'] = cond.node_id.map(studies.groupby('node_id').size()).fillna(0).astype(int)
    cond['n_credible_sets'] = cond.node_id.map(loci.groupby('node_id').size()).fillna(0).astype(int)
    cond['n_eas_studies'] = cond.node_id.map(studies[studies.eas].groupby('node_id').size()).fillna(0).astype(int)
    cond['n_mendelian_genes'] = cond.node_id.map(M.groupby('node_id').ensembl_id.nunique() if len(M)
                                                 else pd.Series(dtype=int)).fillna(0).astype(int)
    regions = pd.DataFrame([{k: r[k] for k in ('gene', 'rsid', 'chrom', 'pos_grch37', 'maf', 'or', 'p', 'p_pe', 'v',
                                              'tier', 'evidence', 'variant_id', 'lookup', 'cteph_study_locus_id')}
                            | {'genes': ' '.join(r['genes']),
                               'n_overlap_sets': int(sum(1 for x in r['overlaps'] if (x.get('pip') or 0) > 0))}
                            for r in pubr])
    ladder = gene_ladder(cfi, cst, load_themes(a.themes), pubr)
    for name, df in [('cteph_regions', regions), ('cteph_genes', ladder), ('overlap', ov),
                     ('conditions', cond), ('relations', rel), ('colocalisation', co), ('gwas_studies', studies),
                     ('gwas_loci', loci), ('mendelian', M), ('mendelian_records', mrec), ('cteph_targets', targets),
                     ('literature', lit), ('gene_condition', GCs)]:
        write_tsv(df, out / f'{name}.tsv')
    for f in ('studies', 'findings', 'methods', 'screening'):
        shutil.copy(Q / f'{f}.tsv', out / f'cteph_{f}.tsv')
    stats = {'release': meta['data_version'], 'api_version': meta['api_version'], 'fetch_date': meta['fetch_date'],
             'fetched_at': meta['fetched_at'], 'conditions': len(cond),
             'conditions_by_source': {s: int(cond.sources.str.contains(s).sum())
                                      for s in ('clinical', 'ontology', 'colocalization', 'lead-variant overlap')},
             'relations': len(rel), 'colocalizations': len(co), 'overlaps': len(ov), 'gwas_studies': len(studies),
             'published_regions': len(regions), 'published_in_ot': int(regions.variant_id.ne('').sum()),
             'published_with_credible_set': int(regions.cteph_study_locus_id.ne('').sum()),
             'cteph_genes_by_level': {int(k): int(v) for k, v in ladder.level.value_counts().sort_index().items()},
             'credible_sets': len(loci), 'eas_studies': int(studies.eas.sum()), 'mendelian_links': len(M),
             'mendelian_records': len(mrec), 'cteph_targets': len(targets), 'literature_rows': len(lit),
             'literature_pmids': int(lit.pmid.nunique()), 'cteph_genetic_records': len(csc),
             'cteph_genetic_studies': len(cst), 'cteph_findings': len(cfi), 'methods_fields': len(cmf),
             'gene_condition': len(GCs), 'excluded_publications_removed': excluded_use,
             'extract': json.loads((X / 'extract.json').read_text()),
             'curation': json.loads((Q / 'curation.json').read_text())}
    (out / 'tidy.json').write_text(json.dumps(stats, indent=1))
    log(f'{len(cond)} conditions, {len(studies)} GWAS studies, {len(loci)} credible sets, {len(M)} Mendelian links, '
        f'{len(cst)} CTEPH genetic studies (curated); excluded publications removed: {excluded_use}')


if __name__ == '__main__':
    main()
