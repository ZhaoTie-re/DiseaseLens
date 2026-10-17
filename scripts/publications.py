#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The PubMed record of every publication the knowledge graph uses, for
#           its status and its place: the genetic-study search and the curated
#           seed (00.pubmed), the text-mined CTEPH literature, the GWAS
#           publications and the references of the curated gene evidence
#           (00.raw), and the relation evidence (01.conditions). A record:
#           publication types, journal issue and pagination, whether PubMed holds
#           an abstract, title, journal, year, DOI, PMCID, the first author and
#           the first author's affiliation (from Europe PMC when PubMed has
#           none), and the comments and corrections
#           PubMed links to it (RetractionIn, ErratumIn, ...). curate.py decides
#           from these which publications are left out (retracted, conference
#           material) and the place of each study.
#           Europe PMC identifiers that are not PMIDs (PMC..., PPR...) cannot be
#           looked up in PubMed and are listed apart.
# Inputs  : the snapshot's PubMed records, 00.raw, its disease files, relations.tsv
# Outputs : publications.json.gz {checked_at, n, roles {pmid: [where used]},
#           records [...], not_found [...], other_ids {id: [where used]}}
# Usage   : publications.py --pubmed <search> <seed> --raw 00.raw --diseases <disease.*.json.gz>
#                           --relations relations.tsv --out publications.json.gz
# ---------------------------------------------------------------------------
import argparse
import datetime as dt
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import load_json, log, read_tsv, write_json  # noqa: E402

EFETCH = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
EPMC = 'https://www.ebi.ac.uk/europepmc/webservices/rest/searchPOST'
_last = [0.0]


def post(url, data, tries=5):
    """POST with a 0.4 s spacing (NCBI's limit without a key) and backoff on server errors / 429."""
    body = urllib.parse.urlencode(data).encode()
    err = None
    for k in range(tries):
        wait = 0.4 - (time.time() - _last[0])
        if wait > 0:
            time.sleep(wait)
        _last[0] = time.time()
        try:
            with urllib.request.urlopen(url, body, timeout=180) as r:
                return r.read().decode('utf-8')
        except urllib.error.HTTPError as e:
            if e.code < 500 and e.code != 429:
                raise RuntimeError(f'HTTP {e.code}: {url}')
            err = e
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            err = e
        time.sleep(2 ** k * 2)
    raise RuntimeError(f'unreachable after {tries} tries: {url}: {err}')


def text_of(el):
    return ' '.join(''.join(el.itertext()).split()) if el is not None else ''


def parse(xml):
    out = []
    root = ET.fromstring(xml)
    for art in root.findall('.//PubmedArticle'):
        mc = art.find('MedlineCitation')
        a = mc.find('Article')
        ji = a.find('Journal/JournalIssue')
        auth = a.findall('AuthorList/Author')
        first, aff = '', ''
        if auth:
            f = auth[0]
            first = (f.findtext('LastName') or f.findtext('CollectiveName') or '') + \
                (' ' + f.findtext('Initials') if f.findtext('Initials') else '')
            aff = text_of(f.find('AffiliationInfo/Affiliation'))
        aff_source = 'first author' if aff else ''
        if not aff:   # a collective first author: the first affiliation of the list
            aff = next((text_of(x) for x in a.findall('AuthorList/Author/AffiliationInfo/Affiliation')), '')
            aff_source = 'first listed affiliation' if aff else ''
        ids = {i.get('IdType'): i.text for i in art.findall('PubmedData/ArticleIdList/ArticleId')}
        year = (ji.findtext('PubDate/Year') or (ji.findtext('PubDate/MedlineDate') or '')[:4]
                or (mc.findtext('DateCompleted/Year') or ''))
        out.append({'pmid': mc.findtext('PMID'), 'title': text_of(a.find('ArticleTitle')),
                    'journal': a.findtext('Journal/ISOAbbreviation') or a.findtext('Journal/Title') or '',
                    'year': year, 'volume': ji.findtext('Volume') or '', 'issue': ji.findtext('Issue') or '',
                    'pages': a.findtext('Pagination/MedlinePgn') or '',
                    'pub_types': [text_of(p) for p in a.findall('PublicationTypeList/PublicationType')],
                    'has_abstract': a.find('Abstract/AbstractText') is not None,
                    'first_author': first, 'n_authors': len(auth), 'affiliation': aff, 'affiliation_source': aff_source,
                    'comments': sorted({c.get('RefType') for c in mc.findall('CommentsCorrectionsList/'
                                                                             'CommentsCorrections')}),
                    'doi': ids.get('doi', ''), 'pmcid': ids.get('pmc', '')})
    for art in root.findall('.//PubmedBookArticle'):
        bd = art.find('BookDocument')
        auth = bd.findall('AuthorList/Author')
        out.append({'pmid': bd.findtext('PMID'),
                    'title': text_of(bd.find('ArticleTitle')) or text_of(bd.find('Book/BookTitle')),
                    'journal': text_of(bd.find('Book/BookTitle')), 'year': bd.findtext('Book/PubDate/Year') or '',
                    'volume': '', 'issue': '', 'pages': '',
                    'pub_types': [text_of(p) for p in bd.findall('PublicationType')] or ['Book'],
                    'has_abstract': bd.find('Abstract/AbstractText') is not None,
                    'first_author': (auth[0].findtext('LastName') or '') if auth else '', 'n_authors': len(auth),
                    'affiliation': text_of(bd.find('AuthorList/Author/AffiliationInfo/Affiliation')),
                    'affiliation_source': 'first author', 'comments': [], 'doi': '', 'pmcid': ''})
    return out


def epmc_affiliations(pmids, batch=50):
    """{pmid: affiliation} from Europe PMC, for records whose PubMed entry has none (older records)."""
    import json
    out = {}
    for i in range(0, len(pmids), batch):
        q = ' OR '.join(f'EXT_ID:{p}' for p in pmids[i:i + batch])
        d = json.loads(post(EPMC, {'query': f'({q}) AND SRC:MED', 'resultType': 'core', 'format': 'json',
                                   'pageSize': 1000}))
        for r in d['resultList']['result']:
            au = ((r.get('authorList') or {}).get('author') or [{}])[0]
            aff = (r.get('affiliation') or (au.get('authorAffiliationDetailsList') or {}).get(
                'authorAffiliation', [{}])[0].get('affiliation') or '').strip()
            if aff and r.get('pmid'):
                out[r['pmid']] = aff
    return out


def collect(a):
    """{identifier: set of roles}: every publication the snapshot's sources name."""
    roles = {}

    def add(pid, role):
        pid = str(pid or '').strip()
        if pid:
            roles.setdefault(pid, set()).add(role)
    for f in a.pubmed:
        d = load_json(f)
        role = 'genetic-study search' if d.get('mode') == 'search' else 'curated seed'
        for r in d['records']:
            add(r['pmid'], role)
    raw = Path(a.raw)
    for e in load_json(raw / 'evidence.json.gz'):
        if e['datasourceId'] == 'europepmc':
            for p in e.get('literature') or []:
                add(p, 'text-mined literature')
    for s in load_json(raw / 'gwas.json.gz')['studies']:
        add(s.get('pubmedId'), 'GWAS publication')
    for f in a.diseases:
        d = load_json(f)
        for s in d['studies']:
            add(s.get('pubmedId'), 'GWAS publication')
        for e in d['mendelian_evidence']:
            for p in e.get('literature') or []:
                add(p, 'curated gene evidence')
    if a.relations:
        for p in read_tsv(a.relations).pmid:
            add(p, 'relation evidence')
    return roles


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pubmed', nargs='+', required=True, help='pubmed.py outputs (search, seed)')
    ap.add_argument('--raw', required=True, help='00.raw of the snapshot')
    ap.add_argument('--diseases', nargs='*', default=[], help='disease.<id>.json.gz')
    ap.add_argument('--relations', help='01.conditions/relations.tsv')
    ap.add_argument('--batch', type=int, default=200)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    roles = collect(a)
    pmids = sorted((p for p in roles if p.isdigit()), key=int)
    other = {p: sorted(r) for p, r in roles.items() if not p.isdigit()}
    recs = []
    for i in range(0, len(pmids), a.batch):
        recs += parse(post(EFETCH, {'db': 'pubmed', 'id': ','.join(pmids[i:i + a.batch]), 'retmode': 'xml'}))
    got = {r['pmid'] for r in recs}
    missing = [p for p in pmids if p not in got]
    need = [r['pmid'] for r in recs if not r['affiliation']]
    epmc = epmc_affiliations(need) if need else {}
    for r in recs:
        if not r['affiliation'] and r['pmid'] in epmc:
            r['affiliation'], r['affiliation_source'] = epmc[r['pmid']], 'Europe PMC'
    log(f'affiliations: {len(need)} records without one in PubMed, {len(epmc)} found in Europe PMC')
    write_json({'checked_at': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'), 'n': len(recs),
                'roles': {p: sorted(roles[p]) for p in pmids}, 'records': sorted(recs, key=lambda r: int(r['pmid'])),
                'not_found': missing, 'other_ids': other}, a.out)
    log(f'{len(pmids)} PMIDs: {len(recs)} PubMed records, {len(missing)} not found; {len(other)} other identifiers')


if __name__ == '__main__':
    main()
