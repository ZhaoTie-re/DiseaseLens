#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Query the Open Targets Platform GraphQL API and keep the responses
#           as they came (gzipped JSON), every list's length and unique IDs
#           asserted against the count the API reports.
#             --mode core      the disease itself:
#                                meta.json            release, fetch time
#                                disease.json         description, synonyms, parents
#                                associations.json.gz every associated target, datasource scores
#                                evidence.json.gz     text-mined literature and GWAS evidence
#                                targets.json.gz      the evidence targets' Reactome
#                                                     top-level pathways (mechanism themes)
#                                gwas.json.gz         the disease's GWAS studies and
#                                                     credible sets, each with its
#                                                     colocalisations (GWAS and molecular
#                                                     QTL) and top locus-to-gene genes
#                                mendelian.json.gz    its curated / Mendelian / burden
#                                                     evidence (as for a related disease)
#                                published_loci.json.gz  the published GWAS regions of the
#                                                     profile (gwas_regions.yaml), if any:
#                                                     each lead variant's Open Targets ID,
#                                                     whether a credible set of the disease
#                                                     exists there, and every GWAS credible
#                                                     set of another study that contains
#                                                     it, with the variant's PIP in that set
#                                fetch_core.json      the counts
#             --mode related   every condition of a fetch list (resolve.py), three at a
#                              time: for each, disease.<id>.json.gz: its record, GWAS
#                              studies and credible sets, and the genetic evidence of
#                              the curated / Mendelian / burden sources (when the list
#                              says so)
#             --mode disease   one related disease or trait, as above
# Inputs  : the API URL, the profile (core) or disease IDs
# Outputs : as above, in --out-dir
# ---------------------------------------------------------------------------
import argparse
import datetime as dt
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import get_json, gql, load_profile, log, release, write_json  # noqa: E402

# genetic evidence that is not GWAS: curated gene-disease validity (ClinGen, Gene2Phenotype), rare-disease
# panels (Genomics England PanelApp, Orphanet), germline ClinVar (eva), UniProt, and gene burden tests
MENDELIAN = ['clingen', 'gene2phenotype', 'genomics_england', 'orphanet', 'eva', 'uniprot_literature',
             'uniprot_variants', 'gene_burden']
DISEASE = '''query($id: String!) { disease(efoId: $id) { id name description dbXRefs
  synonyms { relation terms } parents { id name synonyms { relation terms } }
  children { id name synonyms { relation terms } } therapeuticAreas { id name } } }'''
ASSOC = '''query($id: String!, $index: Int!, $size: Int!) { disease(efoId: $id) {
  associatedTargets(page: {index: $index, size: $size}) { count rows { score
    target { id approvedSymbol approvedName biotype }
    datasourceScores { id score } } } } }'''
EVIDENCE = '''query($id: String!, $ids: [String!]!, $ds: [String!], $size: Int!, $cursor: String) {
  disease(efoId: $id) { evidences(ensemblIds: $ids, datasourceIds: $ds, enableIndirect: true, size: $size,
                                  cursor: $cursor) {
    count cursor rows { id datasourceId score target { id approvedSymbol }
      literature publicationYear textMiningSentences { section text tStart tEnd dStart dEnd }
      credibleSet { studyLocusId } studyId } } } }'''
TARGETS = '''query($ids: [String!]!) { targets(ensemblIds: $ids) { id approvedSymbol approvedName
  pathways { pathwayId pathway topLevelTerm } } }'''
STUDY_FIELDS = '''id projectId studyType traitFromSource publicationFirstAuthor publicationDate publicationJournal
    pubmedId nSamples nCases nControls cohorts hasSumstats analysisFlags
    ldPopulationStructure { ldPopulation relativeSampleSize } diseases { id name }'''
STUDIES = '''query($ids: [String!], $indirect: Boolean, $index: Int!, $size: Int!) {
  studies(diseaseIds: $ids, enableIndirect: $indirect, page: {index: $index, size: $size}) { count rows { %s
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
VARIANT_AF = '''query($id: String!) { variant(variantId: $id) { id alleleFrequencies { populationName alleleFrequency } } }'''
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


MAX_PAGE = 3000    # the largest page the Open Targets API serves


def paged(url, query, variables, path, size, key):
    """All rows of an index / size paged list; path: keys down to the paged object, key: a row's unique ID.
    Successive pages are NOT stable where the sort key ties (the 615 CTEPH targets paged by 100 differed from one
    full page in 10 targets, 2026-10-16), so a list longer than one page is re-read as a single page of its
    full length; a list longer than the API's largest page (MAX_PAGE) is read in pages of that size, again and
    again, until the union of the rows (by key) holds them all. The rows must be unique and as many as the API's
    count."""
    def page(i, n):
        d = gql(url, query, {**variables, 'index': i, 'size': n})
        for k in path:
            d = d[k]
        return d
    d = page(0, min(size, MAX_PAGE))
    count, rows = d['count'], d['rows']
    if count > len(rows) and count <= MAX_PAGE:
        rows = page(0, count)['rows']
    elif count > len(rows):
        got = {}
        for _ in range(10):
            for i in range(-(-count // MAX_PAGE)):
                for r in page(i, MAX_PAGE)['rows']:
                    got.setdefault(key(r), r)
            if len(got) >= count:
                break
        rows = list(got.values())
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


def gwas(url, disease, size, coloc=False, indirect=False):
    """The GWAS studies of a disease (an ID, or several), each with its credible sets (and, for the core disease,
    their colocalisations); indirect: with the studies of its descendants (the core disease)."""
    ids = [disease] if isinstance(disease, str) else list(disease)
    studies, n_st = paged(url, STUDIES, {'ids': ids, 'indirect': indirect}, ['studies'], size, lambda r: r['id'])
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
    """Open Targets variant IDs of an rsID: the hits of its search whose rsIDs include it (the search also matches
    longer rsIDs that begin with it: rs6025 finds rs60254775); else the GRCh38 positions and alleles from Ensembl
    (the search does not index every rsID), each checked against Open Targets."""
    hits = [h['id'] for h in gql(url, SEARCH_VARIANT, {'q': rsid})['search']['hits']]
    hits = [i for i in hits if rsid in ((gql(url, VARIANT, {'id': i})['variant'] or {}).get('rsIds') or [])]
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
    """Each published GWAS region of the disease: its lead variant in Open Targets, the disease's credible set
    there (if any), and every GWAS credible set of another study that contains the lead, with the lead's PIP in
    that set."""
    out = []
    for r in cfg['regions']:
        ids, how = variant_ids(url, r['rsid'])
        ens = {}
        for gname in r['genes']:
            hit = [h['id'] for h in gql(url, SEARCH_TARGET, {'q': gname})['search']['hits'] if h['name'] == gname]
            assert hit, f'{gname}: no Open Targets target of that symbol'
            ens[gname] = hit[0]
        rec = {**r, 'ensembl': ens, 'variant_id': '', 'lookup': how, 'disease_study_locus_id': '', 'overlaps': []}
        if not ids:
            rec['lookup'] = 'not in Open Targets'
            out.append(rec)
            continue
        # a site with several alleles: the first allele that holds credible sets (else the first allele)
        rows, rec['variant_id'] = [], ids[0]
        for vid in ids:
            got, _ = paged(url, VARIANT_CS, {'id': vid}, ['variant', 'credibleSets'], page_size,
                           lambda x: x['studyLocusId'])
            if got:
                rows, rec['variant_id'] = got, vid
                break
        for x in rows:
            if x['studyType'] != 'gwas':
                continue
            if x['study']['id'] == study:
                rec['disease_study_locus_id'] = x['studyLocusId']
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


def fetch_disease(url, disease, mendelian_on, threshold, page_size, gwas_mode='full'):
    """One related disease: its record, GWAS studies and credible sets, and curated / Mendelian evidence."""
    dis = gql(url, DISEASE, {'id': disease})['disease']
    assert dis, f'{disease}: not in Open Targets'
    studies, n_st, n_cs, _ = gwas(url, disease, page_size) if gwas_mode == 'full' else ([], 0, 0, 0)
    ev, tg = mendelian(url, disease, threshold, 500) if mendelian_on else ([], {})
    log(f'{dis["name"]}: {n_st} GWAS studies, {n_cs} credible sets; genetic evidence rows {len(ev)} '
        f'({", ".join(f"{k} {len(v)}" for k, v in tg.items()) or "none"})')
    return {'release': release(url), 'disease': dis, 'gwas_mode': gwas_mode, 'studies': studies, 'n_studies': n_st,
            'n_credible_sets': n_cs, 'mendelian_threshold': threshold, 'mendelian_targets': tg,
            'mendelian_evidence': ev}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mode', choices=['core', 'related', 'disease'], required=True)
    ap.add_argument('--api', required=True)
    ap.add_argument('--profile', help='core: profile.json of the disease')
    ap.add_argument('--disease', help='disease: an EFO / MONDO / HP / OBA ID')
    ap.add_argument('--list', help='related: fetch_list.tsv (ID, curated evidence true/false)')
    ap.add_argument('--fetch-date', default='', help='core: the snapshot date, recorded in meta.json')
    ap.add_argument('--gwas', choices=['full', 'none'], default='full', help='disease: fetch its GWAS')
    ap.add_argument('--mendelian', default='true', help='disease: fetch its curated / Mendelian evidence')
    ap.add_argument('--threshold', type=float, default=0.5, help='smallest datasource score kept (curated)')
    ap.add_argument('--parallel', type=int, default=3, help='related: conditions fetched at a time')
    ap.add_argument('--page-size', type=int, default=100, help='first page of an index / size list')
    ap.add_argument('--evidence-page', type=int, default=500)
    ap.add_argument('--out-dir', default='.')
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rel = release(a.api)
    on = lambda v: str(v).lower() in ('true', '1', 'yes')  # noqa: E731
    if a.mode == 'related':
        rows = [ln.split('\t') for ln in Path(a.list).read_text().splitlines() if ln.strip()]

        def one(r):
            d = fetch_disease(a.api, r[0], on(r[1]) if len(r) > 1 else True, a.threshold, a.page_size)
            write_json(d, out / f'disease.{r[0]}.json.gz')
        with ThreadPoolExecutor(a.parallel) as ex:
            list(ex.map(one, rows))
        log(f'{len(rows)} related conditions fetched')
        return
    if a.mode == 'disease':
        write_json(fetch_disease(a.api, a.disease, on(a.mendelian), a.threshold, a.page_size, a.gwas),
                   out / f'disease.{a.disease}.json.gz')
        return
    P = load_profile(a.profile)
    dis = gql(a.api, DISEASE, {'id': P.id})['disease']
    assert dis, f'{P.id}: not in Open Targets'
    log(f'Open Targets {rel["data_version"]} (API {rel["api_version"]}); {P.id} {dis["name"]}')
    # the disease and the other entries its profile names (evidence_also), each with its descendants; merged by
    # target (best score per source) and by evidence record
    also = [x for x in P.d.get('evidence_also') or [] if x != P.id]
    also_recs = []
    for x in also:
        r = gql(a.api, DISEASE, {'id': x})['disease']
        assert r, f'{x} (evidence_also): not in Open Targets'
        also_recs.append(r)
    by_target, evid, men, men_targets, per_entry = {}, {}, {}, {}, []
    for did in [P.id] + also:
        rows, n = paged(a.api, ASSOC, {'id': did}, ['disease', 'associatedTargets'], a.page_size,
                        lambda r: r['target']['id'])
        for r in rows:
            m = by_target.setdefault(r['target']['id'], {**r, 'datasourceScores': []})
            m['score'] = max(m['score'], r['score'])
            ds = {s['id']: s['score'] for s in m['datasourceScores']}
            for s in r['datasourceScores']:
                ds[s['id']] = max(ds.get(s['id'], 0), s['score'])
            m['datasourceScores'] = [{'id': k, 'score': v} for k, v in ds.items()]
        ids = [r['target']['id'] for r in rows]
        ev, n_e = cursored(a.api, EVIDENCE, {'id': did, 'ids': ids, 'size': a.evidence_page,
                                             'ds': ['europepmc', 'gwas_credible_sets']},
                           ['disease', 'evidences'], lambda r: r['id']) if ids else ([], 0)
        evid.update({e['id']: e for e in ev})
        mv, mt = mendelian(a.api, did, a.threshold, 500)
        men.update({e['id']: e for e in mv})
        for ds, ts in mt.items():
            cur = {t_['id']: t_ for t_ in men_targets.get(ds, [])}
            for t_ in ts:
                if t_['id'] not in cur or t_['score'] > cur[t_['id']]['score']:
                    cur[t_['id']] = t_
            men_targets[ds] = list(cur.values())
        per_entry.append({'id': did, 'targets': n, 'evidence': n_e, 'curated_evidence': len(mv)})
    assoc = sorted(by_target.values(), key=lambda r: (-r['score'], r['target']['id']))
    n_assoc, evid, n_ev, men = len(assoc), list(evid.values()), len(evid), list(men.values())
    if also:
        log('evidence also from ' + '; '.join(f'{x["id"]} {x["name"]}' for x in also_recs) + ': '
            + '; '.join(f'{e["id"]} {e["targets"]} targets, {e["evidence"]} literature / GWAS rows, '
                        f'{e["curated_evidence"]} curated rows' for e in per_entry))
    tids = sorted({e['target']['id'] for e in evid} | {t['id'] for ts in men_targets.values() for t in ts})
    targets = []
    for i in range(0, len(tids), 200):
        targets += gql(a.api, TARGETS, {'ids': tids[i:i + 200]})['targets']
    assert len(targets) == len(tids), 'a target record is missing'
    studies, n_st, n_cs, n_co = gwas(a.api, [P.id] + also, a.page_size, coloc=True, indirect=True)
    # the allele frequencies (gnomAD, by population) of each lead variant: a region needs a common lead (RegionMinMaf)
    af = {}
    for s in studies:
        for c in s['credible_sets']:
            vid = (c.get('variant') or {}).get('id')
            if vid and vid not in af:
                v = gql(a.api, VARIANT_AF, {'id': vid})['variant'] or {}
                af[vid] = {x['populationName']: x['alleleFrequency'] for x in v.get('alleleFrequencies') or []}
            c['lead_af'] = af.get(vid, {})
    log(f'{dis["name"]}: {n_assoc} targets, {n_ev} literature / GWAS evidence rows, {n_st} GWAS studies, '
        f'{n_cs} credible sets, {n_co} colocalisations; curated evidence rows {len(men)} '
        f'({", ".join(f"{k} {len(v)}" for k, v in men_targets.items()) or "none"})')
    meta = {**rel, 'api': a.api, 'fetch_date': a.fetch_date, 'disease': P.id,
            'fetched_at': dt.datetime.now().astimezone().isoformat(timespec='seconds')}
    write_json(meta, out / 'meta.json')
    # the other entries (evidence_also): their names and synonyms are words of the disease for the curation
    write_json({**dis, 'evidence_also': also_recs}, out / 'disease.json')
    write_json(assoc, out / 'associations.json.gz')
    write_json(evid, out / 'evidence.json.gz')
    write_json(targets, out / 'targets.json.gz')
    write_json({'studies': studies, 'n_studies': n_st, 'n_credible_sets': n_cs, 'n_colocalisations': n_co},
               out / 'gwas.json.gz')
    write_json({'release': rel, 'mendelian_threshold': a.threshold, 'mendelian_targets': men_targets,
                'mendelian_evidence': men}, out / 'mendelian.json.gz')
    pub = []
    if P.gwas:
        cfg = P.gwas
        gw = {s['id'] for s in studies if s.get('pubmedId') == str(cfg['source'].get('pmid', ''))}
        assert len(gw) <= 1, f'the published GWAS (PMID {cfg["source"].get("pmid")}) is several studies: {gw}'
        pub = published_loci(a.api, cfg, gw.pop() if gw else '', a.page_size)
        write_json({'source': cfg['source'], 'regions': pub}, out / 'published_loci.json.gz')
        log('published regions: ' + '; '.join(f'{r["gene"]} {r["lookup"]}, credible set of the disease '
                                               f'{"yes" if r["disease_study_locus_id"] else "no"}, '
                                               f'{len(r["overlaps"])} GWAS credible sets' for r in pub))
    write_json({'release': rel, 'evidence_also': [{'id': x['id'], 'name': x['name']} for x in also_recs],
                'per_entry': per_entry,
                'associations': n_assoc, 'evidence': n_ev, 'targets': len(targets),
                'studies': n_st, 'credible_sets': n_cs, 'colocalisations': n_co,
                'curated_evidence': len(men), 'curated_targets': len({t['id'] for ts in men_targets.values()
                                                                      for t in ts}),
                'published_regions': len(pub),
                'published_in_ot': sum(1 for r in pub if r['variant_id']),
                'published_with_credible_set': sum(1 for r in pub if r['disease_study_locus_id'])},
               out / 'fetch_core.json')


if __name__ == '__main__':
    main()
