#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : PubMed and Europe PMC records, fetched from the public E-utilities /
#           Europe PMC REST APIs (no key: at most 3 requests per second).
#             --mode ids       the records of given PMIDs
#             --mode seed      the records of every PMID the curated seed cites
#             --mode search    the records a PubMed query finds, published up to
#                              --max-date (the snapshot's fetch date)
#             --mode fulltext  the Methods text of open-access full texts
#                              (Europe PMC JATS XML), for given PMID=PMCID pairs
#           A record: PMID, title, abstract (labelled sections joined), journal,
#           year, first author, author count, publication types, MeSH headings,
#           PMCID, DOI.
# Inputs  : PMIDs (file, one per line) or a query; PMID=PMCID pairs
# Outputs : records.json.gz (ids / search) with the query and its count;
#           fulltext.json.gz {pmid: {pmcid, methods, sections}}
# ---------------------------------------------------------------------------
import argparse
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import log, write_json  # noqa: E402

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
EPMC = 'https://www.ebi.ac.uk/europepmc/webservices/rest/'
METHOD_HEAD = re.compile(r'method|material|patient|subject|participant|population|genotyp|sequenc|statistic|'
                         r'study design|cohort|sample', re.I)
_last = [0.0]


def get(url, tries=5, raw=False):
    """GET with a 0.4 s spacing (NCBI's limit without a key) and backoff on server errors / 429."""
    err = None
    for k in range(tries):
        wait = 0.4 - (time.time() - _last[0])
        if wait > 0:
            time.sleep(wait)
        _last[0] = time.time()
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                b = r.read()
                return b if raw else b.decode('utf-8')
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
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
    for art in ET.fromstring(xml).findall('.//PubmedArticle'):
        mc = art.find('MedlineCitation')
        a = mc.find('Article')
        pmid = mc.findtext('PMID')
        parts = []
        for ab in a.findall('Abstract/AbstractText'):
            lab = ab.get('Label')
            t = text_of(ab)
            parts.append(f'{lab}: {t}' if lab and t else t)
        year = (a.findtext('Journal/JournalIssue/PubDate/Year')
                or (a.findtext('Journal/JournalIssue/PubDate/MedlineDate') or '')[:4]
                or (mc.findtext('DateCompleted/Year') or ''))
        auth = a.findall('AuthorList/Author')
        first = ''
        if auth:
            f = auth[0]
            first = (f.findtext('LastName') or f.findtext('CollectiveName') or '') + \
                (' ' + f.findtext('Initials') if f.findtext('Initials') else '')
        ids = {i.get('IdType'): i.text for i in art.findall('PubmedData/ArticleIdList/ArticleId')}
        out.append({'pmid': pmid, 'title': text_of(a.find('ArticleTitle')), 'abstract': ' '.join(p for p in parts if p),
                    'journal': a.findtext('Journal/ISOAbbreviation') or a.findtext('Journal/Title') or '',
                    'year': year, 'first_author': first, 'n_authors': len(auth),
                    'pub_types': [text_of(p) for p in a.findall('PublicationTypeList/PublicationType')],
                    'mesh': [text_of(m.find('DescriptorName')) for m in mc.findall('MeshHeadingList/MeshHeading')],
                    'pmcid': ids.get('pmc', ''), 'doi': ids.get('doi', '')})
    # books and book chapters (e.g. GeneReviews) come as PubmedBookArticle
    for art in ET.fromstring(xml).findall('.//PubmedBookArticle'):
        bd = art.find('BookDocument')
        parts = []
        for ab in bd.findall('Abstract/AbstractText'):
            lab, t = ab.get('Label'), text_of(ab)
            parts.append(f'{lab}: {t}' if lab and t else t)
        title = text_of(bd.find('ArticleTitle')) or text_of(bd.find('Book/BookTitle'))
        auth = bd.findall('AuthorList/Author')
        first = ''
        if auth:
            first = (auth[0].findtext('LastName') or auth[0].findtext('CollectiveName') or '') + \
                (' ' + auth[0].findtext('Initials') if auth[0].findtext('Initials') else '')
        out.append({'pmid': bd.findtext('PMID'), 'title': title, 'abstract': ' '.join(p for p in parts if p),
                    'journal': text_of(bd.find('Book/BookTitle')), 'year': bd.findtext('Book/PubDate/Year') or '',
                    'first_author': first, 'n_authors': len(auth),
                    'pub_types': [text_of(p) for p in bd.findall('PublicationType')] or ['Book'],
                    'mesh': [], 'pmcid': '', 'doi': ''})
    return out


def fetch_ids(ids, batch=200):
    recs = []
    for i in range(0, len(ids), batch):
        q = urllib.parse.urlencode({'db': 'pubmed', 'id': ','.join(ids[i:i + batch]), 'retmode': 'xml'})
        recs += parse(get(EUTILS + 'efetch.fcgi?' + q))
    got = {r['pmid'] for r in recs}
    missing = [x for x in ids if x not in got]
    assert not missing, f'PMIDs PubMed does not return: {missing}'
    return sorted(recs, key=lambda r: int(r['pmid']))


def search(term, max_date):
    q = urllib.parse.urlencode({'db': 'pubmed', 'term': term, 'retmode': 'json', 'retmax': 10000,
                                'datetype': 'pdat', 'mindate': '1800/01/01', 'maxdate': max_date.replace('-', '/')})
    import json
    d = json.loads(get(EUTILS + 'esearch.fcgi?' + q))['esearchresult']
    ids = d['idlist']
    assert int(d['count']) == len(ids), f'query returns {d["count"]} records, more than one page'
    return ids


def methods_text(pmcid):
    """The Methods text of an open-access full text: sections whose title names methods, patients, genotyping,
    sequencing or statistics (with their subsections); '' when Europe PMC has no open full text. Europe PMC answers
    some non-open-access PMCIDs with a persistent HTTP 500 (PMC1729855, 2026-10-16): that is 'no full text' too."""
    try:
        xml = get(f'{EPMC}{pmcid}/fullTextXML', tries=3)
    except RuntimeError as e:
        log(f'{pmcid}: no full text ({str(e)[-60:]})')
        return '', ['(not available: ' + str(e)[-40:] + ')']
    if not xml:
        return '', []
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        return '', []
    body = root.find('body')
    if body is None:
        return '', []
    keep, titles = [], []
    for sec in body.findall('sec'):
        t = text_of(sec.find('title'))
        if METHOD_HEAD.search(t) or sec.get('sec-type', '').startswith(('methods', 'materials', 'subjects')):
            titles.append(t)
            keep.append(f'{t}. ' + ' '.join(text_of(p) for p in sec.iter('p')))
    return '\n'.join(keep), titles


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mode', choices=['ids', 'seed', 'search', 'fulltext'], required=True)
    ap.add_argument('--seed', help='seed: config/related_conditions.yaml (every PMID it cites)')
    ap.add_argument('--ids', help='ids: a file of PMIDs; fulltext: a file of PMID=PMCID lines')
    ap.add_argument('--query', help='search: the PubMed query')
    ap.add_argument('--max-date', help='search: YYYY-MM-DD, the newest publication date')
    ap.add_argument('--max-chars', type=int, default=20000, help='fulltext: Methods text kept per article')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    if a.mode in ('ids', 'seed'):
        if a.mode == 'seed':
            import yaml
            s = yaml.safe_load(Path(a.seed).read_text())
            ids = sorted({str(e['pmid']) for c in s['conditions'] for e in (c.get('evidence') or []) +
                          (c.get('east_asian') or [])} | {str(e['pmid']) for e in s.get('cteph_facts') or []}, key=int)
        else:
            ids = sorted({x.strip() for x in Path(a.ids).read_text().split() if x.strip()}, key=int)
        recs = fetch_ids(ids)
        write_json({'mode': a.mode, 'n': len(recs), 'records': recs}, a.out)
        log(f'{len(recs)} PubMed records')
    elif a.mode == 'search':
        ids = search(a.query, a.max_date)
        recs = fetch_ids(ids)
        write_json({'mode': 'search', 'query': a.query, 'max_date': a.max_date, 'count': len(ids), 'n': len(recs),
                    'records': recs}, a.out)
        log(f'query: {len(ids)} records up to {a.max_date}')
    else:
        pairs = [x.split('=') for x in Path(a.ids).read_text().split() if '=' in x]
        out = {}
        for pmid, pmcid in pairs:
            txt, titles = methods_text(pmcid)
            out[pmid] = {'pmcid': pmcid, 'methods': txt[:a.max_chars], 'sections': titles,
                         'truncated': len(txt) > a.max_chars}
        write_json(out, a.out)
        log(f'full texts: {sum(1 for v in out.values() if v["methods"])} of {len(out)} with a Methods text')


if __name__ == '__main__':
    main()
