#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Query the Open Targets Platform GraphQL API and keep the responses
#           as they came (gzipped JSON), every list's length and unique IDs
#           asserted against the count the API reports.
#             --mode cteph     the core disease:
#                                meta.json            release, fetch time
#                                disease.json         description, synonyms, parents
#                                associations.json.gz every associated target, datasource scores
#                                evidence.json.gz     text-mined literature and GWAS evidence
#                                targets.json.gz      the evidence targets' Reactome
#                                                     top-level pathways (mechanism themes)
#                                gwas.json.gz         the core disease's GWAS studies and
#                                                     credible sets, each with its
#                                                     colocalisations (GWAS and molecular
#                                                     QTL) and top locus-to-gene genes
#                                published_loci.json.gz  the published CTEPH GWAS regions
#                                                     (config/cteph_gwas.yaml): each lead
#                                                     variant's Open Targets ID, whether a
#                                                     CTEPH credible set exists there, and
#                                                     every GWAS credible set of another
#                                                     study that contains it, with the
#                                                     variant's PIP in that set
#             --mode disease   one related disease or trait:
#                                disease.<id>.json.gz record; GWAS studies and credible sets
#                                                     (--gwas full); genetic evidence of the
#                                                     curated / Mendelian / burden sources
#                                                     (--mendelian true)
# Inputs  : the API URL and disease IDs (command line)
# Outputs : as above, in --out-dir
# ---------------------------------------------------------------------------
import argparse
import datetime as dt
import sys
import time
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import get_json, gql, log, release, write_json  # noqa: E402

# genetic evidence that is not GWAS: curated gene-disease validity (ClinGen, Gene2Phenotype), rare-disease
# panels (Genomics England PanelApp, Orphanet), germline ClinVar (eva), UniProt, and gene burden tests
MENDELIAN = ['clingen', 'gene2phenotype', 'genomics_england', 'orphanet', 'eva', 'uniprot_literature',
             'uniprot_variants', 'gene_burden']
DISEASE = '''query($id: String!) { disease(efoId: $id) { id name description dbXRefs
  synonyms { relation terms } parents { id name } children { id name } therapeuticAreas { id name } } }'''
ASSOC = '''query($id: String!, $index: Int!, $size: Int!) { disease(efoId: $id) {
  associatedTargets(page: {index: $index, size: $size}) { count rows { score
    target { id approvedSymbol approvedName biotype }
    datasourceScores { id score } } } } }'''
EVIDENCE = '''query($id: String!, $ids: [String!]!, $ds: [String!], $size: Int!, $cursor: String) {
  disease(efoId: $id) { evidences(ensemblIds: $ids, datasourceIds: $ds, size: $size, cursor: $cursor) {
    count cursor rows { id datasourceId score target { id approvedSymbol }
      literature publicationYear textMiningSentences { section text tStart tEnd dStart dEnd }
      credibleSet { studyLocusId } studyId } } } }'''
TARGETS = '''query($ids: [String!]!) { targets(ensemblIds: $ids) { id approvedSymbol approvedName
  pathways { pathwayId pathway topLevelTerm } } }'''
STUDY_FIELDS = '''id projectId studyType traitFromSource publicationFirstAuthor publicationDate publicationJournal
    pubmedId nSamples nCases nControls cohorts hasSumstats analysisFlags
    ldPopulationStructure { ldPopulation relativeSampleSize } diseases { id name }'''
STUDIES = '''query($ids: [String!], $index: Int!, $size: Int!) {
  studies(diseaseIds: $ids, page: {index: $index, size: $size}) { count rows { %s
    credibleSets(page: {index: 0, size: 1}) { count } } } }''' % STUDY_FIELDS
# pValueMantissa / pValueExponent / beta of a credible set are its LEAD variant's own statistics
CS_FIELDS = '''studyLocusId variant { id rsIds } pValueMantissa pValueExponent beta finemappingMethod
    credibleSetIndex confidence l2GPredictions(page: {index: 0, size: 3}) { rows { score target { id approvedSymbol } } }'''
CREDIBLE = '''query($id: String!, $index: Int!, $size: Int!) { study(studyId: $id) {
  credibleSets(page: {index: $index, size: $size}) { count rows { %s } } } }''' % CS_FIELDS
COLOC = '''query($id: String!, $index: Int!, $size: Int!) { credibleSet(studyLocusId: $id) {
  colocalisation(page: {index: $index, size: $size}) { count rows {
    rightStudyType colocalisationMethod h3 h4 clpp numberColocalisingVariants betaRatioSignAverage
    otherStudyLocus { studyLocusId variant { id rsIds } pValueMantissa pValueExponent beta
      study { id studyType traitFromSource projectId publicationFirstAuthor publicationDate pubmedId nSamples
              ldPopulationStructure { ldPopulation relativeSampleSize } diseases { id name }
              target { id approvedSymbol } biosample { biosampleName } } } } } } }'''
SEARCH_VARIANT = '''query($q: String!) { search(queryString: $q, entityNames: ["variant"], page: {index: 0, size: 5}) {
  hits { id } } }'''
VARIANT_CS = '''query($id: String!, $index: Int!, $size: Int!) { variant(variantId: $id) { id rsIds
  credibleSets(page: {index: $index, size: $size}) { count rows { studyLocusId studyType pValueMantissa
    pValueExponent finemappingMethod variant { id rsIds }
    study { id studyType traitFromSource projectId publicationFirstAuthor publicationDate pubmedId nSamples
            ldPopulationStructure { ldPopulation relativeSampleSize } diseases { id name } }
    locus(variantIds: [$id]) { rows { posteriorProbability } } } } } }'''
VARIANT = '''query($id: String!) { variant(variantId: $id) { id rsIds } }'''
SEARCH_TARGET = '''query($q: String!) { search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 5}) {
  hits { id name } } }'''
BY_SOURCE = '''query($id: String!, $order: String!, $size: Int!) { disease(efoId: $id) {
  associatedTargets(enableIndirect: true, orderByScore: $order, page: {index: 0, size: $size}) { count
    rows { target { id approvedSymbol } datasourceScores { id score } } } } }'''
GEN_EVIDENCE = '''query($id: String!, $ids: [String!]!, $ds: [String!], $cursor: String) { disease(efoId: $id) {
  evidences(ensemblIds: $ids, datasourceIds: $ds, enableIndirect: true, size: 500, cursor: $cursor) {
    count cursor rows { id datasourceId score target { id approvedSymbol } disease { id name }
      confidence clinicalSignificances alleleOrigins variant { id rsIds } variantRsId studyId literature
      diseaseFromSource cohortShortName statisticalMethod pValueMantissa pValueExponent oddsRatio beta
      ancestry studySampleSize } } } }'''


def paged(url, query, variables, path, size, key):
    """All rows of an index / size paged list; path: keys down to the paged object, key: a row's unique ID.
    Successive pages are NOT stable where the sort key ties (the 615 CTEPH targets paged by 100 differed from one
    full page in 10 targets, 2026-10-16), so a list longer than one page is re-read as a single page of its
    full length. The rows must be unique and as many as the API's count."""
    def page(n):
        d = gql(url, query, {**variables, 'index': 0, 'size': n})
        for k in path:
            d = d[k]
        return d
    d = page(size)
    count, rows = d['count'], d['rows']
    if count > len(rows):
        rows = page(count)['rows']
    keys = [key(r) for r in rows]
    assert len(rows) == count and len(set(keys)) == count, (path, len(rows), len(set(keys)), count)
    return rows, count


def cursored(url, query, variables, path, key):
    rows, count, cursor = [], None, None
    while True:
        d = gql(url, query, {**variables, 'cursor': cursor})
        for k in path:
            d = d[k]
        if count is None:
            count = d['count']
        rows += d['rows']
        cursor = d.get('cursor')
        if not cursor or not d['rows'] or len(rows) >= count:
            break
    keys = [key(r) for r in rows]
    assert len(rows) == count and len(set(keys)) == count, (path, len(rows), len(set(keys)), count)
    return rows, count


def gwas(url, disease, size, coloc=False):
    """The disease's GWAS studies, each with its credible sets (and, for the core disease, their colocalisations)."""
    studies, n_st = paged(url, STUDIES, {'ids': [disease]}, ['studies'], size, lambda r: r['id'])
    n_cs = n_co = 0
    for s in studies:
        k = s.pop('credibleSets')['count']
        s['credible_sets'] = []
        if k:
            s['credible_sets'], _ = paged(url, CREDIBLE, {'id': s['id']}, ['study', 'credibleSets'], size,
                                          lambda r: r['studyLocusId'])
        n_cs += len(s['credible_sets'])
        if coloc:
            for c in s['credible_sets']:
                c['colocalisation'], _ = paged(url, COLOC, {'id': c['studyLocusId']},
                                               ['credibleSet', 'colocalisation'], size,
                                               lambda r: (r['otherStudyLocus']['studyLocusId'],
                                                          r['colocalisationMethod']))
                n_co += len(c['colocalisation'])
    return studies, n_st, n_cs, n_co


def variant_ids(url, rsid):
    """Open Targets variant IDs of an rsID: its search first; else the GRCh38 positions and alleles from Ensembl
    (the search does not index every rsID), each checked against Open Targets."""
    hits = [h['id'] for h in gql(url, SEARCH_VARIANT, {'q': rsid})['search']['hits']]
    if hits:
        return hits, 'Open Targets search'
    time.sleep(0.2)
    j = get_json(f'https://rest.ensembl.org/variation/human/{rsid}?content-type=application/json')
    cand = []
    for m in j.get('mappings') or []:
        if m.get('assembly_name') == 'GRCh38' and '_' not in m['seq_region_name']:
            al = m['allele_string'].split('/')
            cand += [f'{m["seq_region_name"]}_{m["start"]}_{al[0]}_{x}' for x in al[1:]]
    return [c for c in cand if gql(url, VARIANT, {'id': c})['variant']], 'Ensembl GRCh38 position'


def published_loci(url, cfg, study, page_size):
    """Each published CTEPH GWAS region: its lead variant in Open Targets, the CTEPH credible set there (if any),
    and every GWAS credible set of another study that contains the lead, with the lead's PIP in that set."""
    out = []
    for r in cfg['regions']:
        ids, how = variant_ids(url, r['rsid'])
        ens = {}
        for gname in r['genes']:
            hit = [h['id'] for h in gql(url, SEARCH_TARGET, {'q': gname})['search']['hits'] if h['name'] == gname]
            assert hit, f'{gname}: no Open Targets target of that symbol'
            ens[gname] = hit[0]
        rec = {**r, 'ensembl': ens, 'variant_id': '', 'lookup': how, 'cteph_study_locus_id': '', 'overlaps': []}
        if not ids:
            rec['lookup'] = 'not in Open Targets'
            out.append(rec)
            continue
        rec['variant_id'] = ids[0]
        rows, _ = paged(url, VARIANT_CS, {'id': ids[0]}, ['variant', 'credibleSets'], page_size,
                        lambda x: x['studyLocusId'])
        for x in rows:
            if x['studyType'] != 'gwas':
                continue
            if x['study']['id'] == study:
                rec['cteph_study_locus_id'] = x['studyLocusId']
                continue
            x['pip'] = ((x.pop('locus') or {}).get('rows') or [{}])[0].get('posteriorProbability')
            rec['overlaps'].append(x)
        out.append(rec)
    return out


def mendelian(url, disease, threshold, page_size):
    """Per curated / Mendelian / burden source: the targets it scores >= threshold for the disease (and its
    descendants), then every evidence row of that source for them."""
    out, targets = [], {}
    for ds in MENDELIAN:
        d = gql(url, BY_SOURCE, {'id': disease, 'order': f'{ds} desc', 'size': page_size})['disease']['associatedTargets']
        sc = [(r['target'], next((x['score'] for x in r['datasourceScores'] if x['id'] == ds), 0)) for r in d['rows']]
        keep = [t for t, s in sc if s >= threshold]
        assert len(keep) < len(sc) or len(sc) == d['count'], f'{disease} {ds}: more targets above the threshold'
        if not keep:
            continue
        targets[ds] = [{'id': t['id'], 'symbol': t['approvedSymbol'], 'score': s} for t, s in sc if s >= threshold]
        rows, _ = cursored(url, GEN_EVIDENCE, {'id': disease, 'ids': [t['id'] for t in keep], 'ds': [ds]},
                           ['disease', 'evidences'], lambda r: r['id'])
        out += rows
    return out, targets


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mode', choices=['cteph', 'disease'], required=True)
    ap.add_argument('--api', required=True)
    ap.add_argument('--disease', required=True, help='EFO / MONDO / HP / OBA ID')
    ap.add_argument('--fetch-date', default='', help='cteph: the snapshot date, recorded in meta.json')
    ap.add_argument('--published', default='', help='cteph: the published CTEPH GWAS regions (YAML)')
    ap.add_argument('--gwas', choices=['full', 'none'], default='full', help='disease: fetch its GWAS')
    ap.add_argument('--mendelian', default='true', help='disease: fetch its curated / Mendelian evidence')
    ap.add_argument('--threshold', type=float, default=0.5, help='disease: smallest datasource score kept')
    ap.add_argument('--page-size', type=int, default=100, help='first page of an index / size list')
    ap.add_argument('--evidence-page', type=int, default=500)
    ap.add_argument('--out-dir', default='.')
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rel = release(a.api)
    dis = gql(a.api, DISEASE, {'id': a.disease})['disease']
    assert dis, f'{a.disease}: not in Open Targets'
    if a.mode == 'cteph':
        log(f'Open Targets {rel["data_version"]} (API {rel["api_version"]})')
        assoc, n_assoc = paged(a.api, ASSOC, {'id': a.disease}, ['disease', 'associatedTargets'], a.page_size,
                               lambda r: r['target']['id'])
        ids = [r['target']['id'] for r in assoc]
        evid, n_ev = cursored(a.api, EVIDENCE, {'id': a.disease, 'ids': ids, 'size': a.evidence_page,
                                                'ds': ['europepmc', 'gwas_credible_sets']},
                              ['disease', 'evidences'], lambda r: r['id'])
        tids = sorted({e['target']['id'] for e in evid})
        targets = []
        for i in range(0, len(tids), 200):
            targets += gql(a.api, TARGETS, {'ids': tids[i:i + 200]})['targets']
        assert len(targets) == len(tids), 'a target record is missing'
        studies, n_st, n_cs, n_co = gwas(a.api, a.disease, a.page_size, coloc=True)
        log(f'{dis["name"]}: {n_assoc} targets, {n_ev} literature / GWAS evidence rows, {n_st} GWAS studies, '
            f'{n_cs} credible sets, {n_co} colocalisations')
        meta = {**rel, 'api': a.api, 'fetch_date': a.fetch_date, 'disease': a.disease,
                'fetched_at': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')}
        write_json(meta, out / 'meta.json')
        write_json(dis, out / 'disease.json')
        write_json(assoc, out / 'associations.json.gz')
        write_json(evid, out / 'evidence.json.gz')
        write_json(targets, out / 'targets.json.gz')
        write_json({'studies': studies, 'n_studies': n_st, 'n_credible_sets': n_cs, 'n_colocalisations': n_co},
                   out / 'gwas.json.gz')
        pub = []
        if a.published:
            cfg = yaml.safe_load(Path(a.published).read_text())
            gw = {s['id'] for s in studies if s.get('pubmedId') == cfg['source']['pmid']}
            assert len(gw) == 1, f'the published CTEPH GWAS (PMID {cfg["source"]["pmid"]}) is not one study: {gw}'
            pub = published_loci(a.api, cfg, gw.pop(), a.page_size)
            write_json({'source': cfg['source'], 'regions': pub}, out / 'published_loci.json.gz')
            log('published regions: ' + '; '.join(f'{r["gene"]} {r["lookup"]}, CTEPH credible set '
                                                   f'{"yes" if r["cteph_study_locus_id"] else "no"}, '
                                                   f'{len(r["overlaps"])} GWAS credible sets' for r in pub))
        write_json({'release': rel, 'associations': n_assoc, 'evidence': n_ev, 'targets': len(targets),
                    'studies': n_st, 'credible_sets': n_cs, 'colocalisations': n_co,
                    'published_regions': len(pub),
                    'published_in_ot': sum(1 for r in pub if r['variant_id']),
                    'published_with_cteph_credible_set': sum(1 for r in pub if r['cteph_study_locus_id'])},
                   out / 'fetch_cteph.json')
    else:
        studies, n_st, n_cs, _ = gwas(a.api, a.disease, a.page_size) if a.gwas == 'full' else ([], 0, 0, 0)
        ev, tg = (mendelian(a.api, a.disease, a.threshold, 500) if a.mendelian.lower() in ('true', '1', 'yes')
                  else ([], {}))
        log(f'{dis["name"]}: {n_st} GWAS studies, {n_cs} credible sets; genetic evidence rows {len(ev)} '
            f'({", ".join(f"{k} {len(v)}" for k, v in tg.items()) or "none"})')
        write_json({'release': rel, 'disease': dis, 'gwas_mode': a.gwas, 'studies': studies, 'n_studies': n_st,
                    'n_credible_sets': n_cs, 'mendelian_threshold': a.threshold, 'mendelian_targets': tg,
                    'mendelian_evidence': ev}, out / f'disease.{a.disease}.json.gz')


if __name__ == '__main__':
    main()
