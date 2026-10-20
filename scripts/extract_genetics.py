#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The genetic studies of a disease in PubMed, screened and extracted
#           by isolated language models (llm.py), then checked against their
#           source text. Two passes:
#             A  abstracts, --batch per call: is the record an original genetic
#                study (or case report) of the disease? If so: design, population,
#                sample sizes, genes, variants, findings, replication and the
#                methods the abstract states
#             B  the methods of every included record (platform, genotyping /
#                sequencing, quality control, statistics, covariates, multiple
#                testing, replication): from the Methods section of its open-
#                access full text, or else from its abstract
#           Two readers: every pass is read by two models that never see each
#           other's answer. Per record, the decisive items are compared:
#           inclusion, sample sizes, the Japanese / East Asian flags, the findings
#           (gene and result) and the methods fields. Where the readers disagree,
#           an adjudicating call reads the record with both readings and gives
#           the final reading. Only agreed or adjudicated items are kept, each
#           with `agreement` (both / adjudicated); every disputed item is listed
#           in review.tsv with both readings and the decision. With one reader,
#           its reading is kept (agreement: one reader).
#           Grounding (an item that fails is dropped and logged, never kept):
#             every gene, variant and finding quote must occur in the abstract
#             (a name of at most 4 characters as a whole word, in its case);
#             sample sizes must occur as numbers; a Japanese / East Asian
#             population or place must be named (whole words; 'Asian' alone
#             is not East Asian); every methods value must come with a quote that occurs in
#             the text it was read from (a field stated twice keeps both).
#           Pre-screen (before any model reads a record): the publication types
#           and genetic-study terms of the shared search rules (config/search.yaml).
#           Result check: a finding whose result disagrees with the wording of
#           its own sentence is flagged (review.tsv; patterns: curation rules).
# Inputs  : the PubMed search records, the Methods texts (pubmed.py --mode fulltext),
#           the disease profile, the model configuration
# Outputs : screening.tsv, studies.tsv, findings.tsv, methods.tsv (field, value, quote, source),
#           review.tsv, gene_pathways.tsv, extract_log.tsv, extract.json (counts); with
#           --fulltext-list only: fulltext_ids.txt (PMID=PMCID of the included records with a PMCID)
# ---------------------------------------------------------------------------
import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from llm import (EAST_ASIA_PLACES, in_symbol, in_text, in_words, load_llm, norm, numbers_in, overlap,  # noqa: E402
                 per_record, read_all, reconcile, record_text)
from ot_common import EAS_WORDS, gql, load_json, load_profile, log, release, write_tsv  # noqa: E402

MAPIDS = '''query($q: [String!]!) { mapIds(queryTerms: $q, entityNames: ["target"]) {
  mappings { term hits { id name } } } }'''
PATHWAYS = '''query($ids: [String!]!) { targets(ensemblIds: $ids) { id approvedSymbol
  pathways { topLevelTerm } } }'''


def map_symbols(api, symbols):
    """{symbol: Ensembl ID} for the symbols Open Targets knows by that exact name."""
    out = {}
    for i in range(0, len(symbols), 100):
        for m in gql(api, MAPIDS, {'q': symbols[i:i + 100]})['mapIds']['mappings']:
            hit = next((h for h in m['hits'] if h['name'] == m['term']), None)
            if hit:
                out[m['term']] = hit['id']
    return out


def pathways_of(api, ens):
    """{symbol: [Reactome top-level pathways]} of the given Ensembl IDs (the automatic gene themes)."""
    out = {}
    ids = sorted(set(ens.values()))
    for i in range(0, len(ids), 200):
        for t in gql(api, PATHWAYS, {'ids': ids[i:i + 200]})['targets']:
            out[t['approvedSymbol']] = sorted({p['topLevelTerm'] for p in t.get('pathways') or []
                                               if p.get('topLevelTerm')})
    return out


STUDY_COLUMNS = ['pmid', 'year', 'first_author', 'journal', 'title', 'pmcid', 'study_type', 'design', 'country',
                 'ancestry', 'japanese', 'east_asian', 'n_cases', 'n_controls', 'control_type', 'genes', 'gene_symbols',
                 'variants', 'replication', 'm_genotyping', 'm_statistics', 'm_covariates', 'm_multiple_testing',
                 'methods_quote', 'agreement']
DESIGNS = ['GWAS', 'candidate gene association', 'whole-exome sequencing', 'whole-genome sequencing',
           'targeted gene panel', 'HLA typing', 'copy-number analysis', 'clonal haematopoiesis / somatic',
           'mitochondrial DNA', 'Mendelian randomization', 'family / linkage', 'other']
def prompts(P):
    """The system prompts of the two passes for the disease. For CTEPH they are, character for character, the
    prompts of the first build (2026-10-16), so that its cached answers are reused."""
    name, short = P.name, P.short
    a = f'''You read PubMed records and extract facts about the GENETICS of {P.d.get('prompt_name') or name} ({short}). Use ONLY the text given; never add outside knowledge. Copy gene symbols, variant names, numbers
and quotes exactly as they are written in the record. Answer with JSON only (no prose).

For every record return one object:
{{"pmid": "<PMID>",
  "include": true or false,   // true only for an ORIGINAL study or case report that reports genetic data (germline
                              // or somatic variants, genotypes, HLA alleles, sequencing, copy number, Mendelian
                              // randomization) in people with {short}. Reviews, guidelines, gene-EXPRESSION-only
                              // studies, animal-only studies and studies without {short} patients: false.
  "reason": "<at most 15 words>",
  "study_type": "original research" | "case report" | "review" | "guideline" | "other",
  "design": [one or more of {json.dumps(DESIGNS)}],
  "population": {{"country": "", "ancestry": "", "japanese": true/false, "east_asian": true/false}},
  "n_cases": <{short} patients genotyped, integer or null>, "n_controls": <integer or null>, "control_type": "",
  "genes": [{{"as_written": "<the gene exactly as the record writes it>", "symbol": "<its official HGNC symbol,
              e.g. F5 for factor V; empty if unsure>"}}],
  "variants": ["<variant / allele as written>"],
  "findings": [{{"gene": "<as written>", "symbol": "<HGNC symbol or empty>", "variant": "", "result": "associated" | "not associated" | "identified" | "enriched" |
                "other", "effect": "<OR / frequency as written, or empty>", "p_value": "<as written, or empty>",
                "quote": "<one exact sentence or clause of the record that states this finding>"}}],
  "replication": "yes" | "no" | "not stated",
  "methods": {{"genotyping": "", "statistics": "", "covariates": "", "multiple_testing": "",
              "quote": "<exact text of the record that states these methods, or empty>"}}}}
For excluded records, give pmid, include, reason and study_type only. Return a JSON array, one object per record,
in the order given. One finding per gene (or per variant); quotes are copied whole, never shortened with '...'.'''
    if P.d.get('scope'):
        a += f'\nWhich studies count as studies of {short}: {P.d["scope"]}'
    b = f'''You read the Methods section of a genetic study of {name} and
extract how the study was done. Use ONLY the text given; copy values exactly. Answer with one JSON object (no prose):
{{"fields": [{{"field": one of "cohort / recruitment", "diagnosis of {short}", "controls", "DNA source",
             "genotyping platform", "sequencing", "HLA typing", "variant calling / quality control",
             "imputation", "statistical test", "covariates", "multiple testing", "replication", "power",
             "software",
             "value": "<short summary in the text's own words>",
             "quote": "<an exact passage of the text that states it>"}}]}}
List only fields the text states; omit the others. Quotes are copied whole, never shortened with '...'.'''
    # the same fields from the abstract, for the included studies without an open-access Methods section (the user,
    # 2026-10-16: "for those whose full text cannot be fetched, fall back to the abstract")
    b_abs = b.replace('the Methods section of a genetic study', 'the title and abstract of a genetic study').replace(
        'extract how the study was done.', 'extract how the study was done (no full text is available).')
    return a, b, b_abs


def ground_a(ans, r, log_rows, reader):
    """The answer for one record, with every ungrounded item removed (and logged)."""
    text = r['title'] + ' ' + r['abstract']
    pm = r['pmid']

    def drop(what, value):
        log_rows.append({'pmid': pm, 'pass': 'abstract', 'item': what, 'value': str(value)[:200], 'reader': reader})

    out = {k: ans.get(k) for k in ('include', 'reason', 'study_type')}
    if out['include'] is not True:          # only a JSON true includes a record (not "false" given as text)
        out['include'] = False
        return out
    design = ans.get('design')
    design = [design] if isinstance(design, str) else design if isinstance(design, list) else []
    out['design'] = [d for d in design if d in DESIGNS] or ['other']
    pop = ans.get('population') if isinstance(ans.get('population'), dict) else {}
    jp = pop.get('japanese') is True and (in_words('Japan', text) or in_words('Japanese', text))
    ea = pop.get('east_asian') is True and any(in_words(w, text) for w in EAS_WORDS + EAST_ASIA_PLACES)
    if pop.get('japanese') is True and not jp:
        drop('japanese', 'not stated in the abstract')
    if pop.get('east_asian') is True and not (ea or jp):
        drop('east_asian', 'not stated in the abstract')
    out.update({'country': pop.get('country') or '', 'ancestry': pop.get('ancestry') or '', 'japanese': jp,
                'east_asian': ea or jp, 'control_type': ans.get('control_type') or '',
                'replication': ans.get('replication') or 'not stated'})
    for k in ('n_cases', 'n_controls'):
        v = ans.get(k)
        if v not in (None, '') and not numbers_in(v, text):
            drop(k, v)
            v = None
        out[k] = v
    def items(k):
        v = ans.get(k)
        return v if isinstance(v, list) else []

    def strs(x, keys):
        return {k: '' if x.get(k) is None else str(x.get(k)) for k in keys}
    genes = []
    for g in items('genes'):
        g = strs(g, ('as_written', 'symbol')) if isinstance(g, dict) else {'as_written': str(g), 'symbol': ''}
        if in_symbol(g['as_written'], text) or (g['symbol'] and in_symbol(g['symbol'], text)):
            genes.append(g)
        else:
            drop('gene', g['as_written'])
    out['genes'] = genes
    out['variants'] = [str(v) for v in items('variants') if in_symbol(str(v), text) or drop('variant', v)]
    fs = []
    for f in items('findings'):
        if not isinstance(f, dict):
            drop('finding', f)
            continue
        f = strs(f, ('gene', 'symbol', 'variant', 'result', 'effect', 'p_value', 'quote'))
        if not in_text(f['quote'], text):
            drop('finding', f['quote'])
            continue
        if f['gene'] and not (in_symbol(f['gene'], text) or (f['symbol'] and in_symbol(f['symbol'], text))):
            drop('finding gene', f['gene'])
            continue
        if f['variant'] and not in_symbol(f['variant'], text):
            drop('finding variant', f['variant'])
            f['variant'] = ''
        fs.append(f)
    out['findings'] = fs
    m = ans.get('methods') if isinstance(ans.get('methods'), dict) else {}
    if isinstance(m.get('quote'), str) and m['quote'] and in_text(m['quote'], text):
        out['methods'] = {k: m.get(k) or '' for k in ('genotyping', 'statistics', 'covariates', 'multiple_testing')}
        out['methods_quote'] = m['quote']
    else:
        if any(m.get(k) for k in ('genotyping', 'statistics', 'covariates', 'multiple_testing')):
            drop('methods', m.get('quote', ''))
        out['methods'], out['methods_quote'] = {}, ''
    if ans.get('decision'):
        out['decision'] = str(ans['decision'])[:300]
    return out


# ── comparing two readings ───────────────────────────────────────────────────
def finding_key(f):
    """A finding as the readers are compared on it: the gene (symbol, else the gene as written; a finding without a
    gene: its variant) and the result."""
    g = (f.get('symbol') or '').strip().upper() or re.sub(r'[\s\-*:()]', '', norm(f.get('gene') or '')) or \
        ('variant ' + re.sub(r'[\s\-*:()]', '', norm(f.get('variant') or ''))).strip()
    res = f.get('result') or ''
    return g, 'associated' if res == 'enriched' else res


def as_int(v):
    try:
        return int(float(str(v).replace(',', '')))
    except (TypeError, ValueError):
        return None


def differences(ga, gb):
    """The decisive items on which two grounded readings of one record differ: [(item, value A, value B)]."""
    if bool(ga.get('include')) != bool(gb.get('include')):
        return [('inclusion', str(bool(ga.get('include'))), str(bool(gb.get('include'))))]
    if not ga.get('include'):
        return []
    out = []
    for k in ('n_cases', 'n_controls'):
        if as_int(ga.get(k)) != as_int(gb.get(k)):
            out.append((k, str(ga.get(k)), str(gb.get(k))))
    for k in ('japanese', 'east_asian'):
        if bool(ga.get(k)) != bool(gb.get(k)):
            out.append((k, str(ga.get(k)), str(gb.get(k))))
    ka, kb = {finding_key(f) for f in ga['findings']}, {finding_key(f) for f in gb['findings']}
    for g, res in sorted(ka - kb):
        out.append(('finding', f'{g}: {res}', 'not reported'))
    for g, res in sorted(kb - ka):
        out.append(('finding', 'not reported', f'{g}: {res}'))
    return out


def result_check(f, rules):
    """'' or why the result of a finding disagrees with the wording of its source sentence (rules: the negative and
    positive patterns of the curation rules, result_check)."""
    q, res = f.get('quote') or '', f.get('result') or ''
    NEGATIVE = re.compile('|'.join(rules['negative']), re.I)
    POSITIVE_WORDS = re.compile('|'.join(rules['positive']), re.I)
    if res == 'not associated' and not NEGATIVE.search(q):
        return 'result "not associated", but the sentence states no negative or non-significant result'
    if res in ('associated', 'enriched') and NEGATIVE.search(q) and not POSITIVE_WORDS.search(q):
        return f'result "{res}", but the sentence reads as negative'
    return ''


def method_fields(ans, text, src, pmid, logs, reader):
    """The grounded methods fields of one answer, in the schema of pass B: {'fields': [{field, value, quote}]}. A
    field stated twice (two cohorts, two platforms) keeps both entries; an entry repeated word for word is kept once."""
    out, seen = [], set()
    fields = ans.get('fields') if isinstance(ans, dict) else None
    for f in fields if isinstance(fields, list) else []:
        if not isinstance(f, dict):
            continue
        f = {k: '' if f.get(k) is None else str(f.get(k)) for k in ('field', 'value', 'quote')}
        if not in_text(f['quote'], text):
            logs.append({'pmid': pmid, 'pass': src, 'item': f['field'], 'value': f['quote'][:200], 'reader': reader})
            continue
        if (f['field'], norm(f['quote'])) not in seen:
            seen.add((f['field'], norm(f['quote'])))
            out.append(f)
    return {'fields': out}


def method_differences(fa, fb):
    """The methods fields on which two readings differ: a field one reading states and the other does not, or whose
    quotes do not overlap (a field stated twice agrees when any quote of one reading overlaps one of the other)."""
    qa, qb = {}, {}
    for q, fs in ((qa, fa), (qb, fb)):
        for f in fs['fields']:
            q.setdefault(f['field'], []).append(f['quote'])
    return [(f'methods: {f}', ' | '.join(qa.get(f, ['not stated']))[:200], ' | '.join(qb.get(f, ['not stated']))[:200])
            for f in sorted(set(qa) | set(qb))
            if not (f in qa and f in qb and any(overlap(x, y) for x in qa[f] for y in qb[f]))]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--records', required=True, help='pubmed.py --mode search output')
    ap.add_argument('--fulltext', help='pubmed.py --mode fulltext output (pass B); omit for pass A only')
    ap.add_argument('--fulltext-list', help='pass A only: write PMID=PMCID of the included records here')
    ap.add_argument('--profile', required=True, help='profile.json of the disease')
    ap.add_argument('--llm', required=True, help='llm.json: the readers and the adjudicator')
    ap.add_argument('--offline', default='false', help='true: cached answers only, never call a model')
    ap.add_argument('--batch', type=int, default=5)
    ap.add_argument('--limit', type=int, default=0, help='smoke run: the first N records only')
    ap.add_argument('--prescreen', default='true', help='drop reviews and records without a genetic-study term '
                    'before the model reads them')
    ap.add_argument('--api', required=True, help='Open Targets GraphQL (symbol check)')
    ap.add_argument('--expect-release', default='', help='rebuild: the Open Targets release of the snapshot; the '
                    'symbol lookups stop if the API serves another')
    ap.add_argument('--out-dir', default='.')
    a = ap.parse_args()
    P = load_profile(a.profile)
    SYS_A, SYS_B, SYS_B_ABS = prompts(P)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    recs = load_json(a.records)['records']
    recs = [r for r in recs if r['abstract'] or r['title']]
    pre = []
    EXCLUDED_TYPES = set(P.d['prescreen']['excluded_types'])
    GENETIC_TERMS = re.compile('|'.join(P.d['prescreen']['genetic_terms']), re.I)
    if a.prescreen.lower() in ('true', '1', 'yes'):
        keep = []
        for r in recs:
            types = set(r['pub_types']) & EXCLUDED_TYPES
            if types or not GENETIC_TERMS.search(r['title'] + ' ' + r['abstract']):
                pre.append({'pmid': r['pmid'], 'year': r['year'], 'title': r['title'], 'include': False,
                            'study_type': 'review' if types else '', 'stage': 'rule',
                            'reason': f'pre-screen: {", ".join(sorted(types))}' if types
                            else 'pre-screen: no genetic-study term in title or abstract', 'agreement': ''})
            else:
                keep.append(r)
        recs = keep
    if a.limit:
        recs = recs[:a.limit]
    readers, adj = load_llm(a.llm, offline=a.offline.lower() in ('true', '1', 'yes'))
    two = len(readers) == 2
    names = [r.name for r in readers]

    # ── pass A: abstracts, read by each reader, reconciled ───────────────────
    logs = []
    by_pmid = {r['pmid']: r for r in recs}
    readings = per_record(readers, recs, SYS_A, a.batch, lambda x, r, name: ground_a(x, r, logs, name), logs)
    finals, review = reconcile(readers, adj, [(r['pmid'], 'Record:\n\n' + record_text(r), readings[r['pmid']])
                                              for r in recs], SYS_A, differences,
                               lambda ans, key, _t: ground_a(ans, by_pmid[key], logs, 'adjudicator'))
    for row in review:
        row.update({'pmid': row.pop('key'), 'pass': 'abstract'})
        if row['pmid'] in finals:
            row['final'] = _final_value(row['item'], finals[row['pmid']][0])
    screen, studies, findings = [], [], []
    for r in recs:
        if r['pmid'] not in finals:
            continue
        g, why = finals[r['pmid']]
        rs = readings[r['pmid']]
        screen.append({'pmid': r['pmid'], 'year': r['year'], 'title': r['title'], 'include': bool(g['include']),
                       'study_type': g.get('study_type') or '', 'stage': 'model', 'reason': g.get('reason') or '',
                       'agreement': why})
        if not g['include']:
            continue
        sym = [x.get('symbol') or '' for x in g['genes']]
        studies.append({'pmid': r['pmid'], 'year': r['year'], 'first_author': r['first_author'],
                        'journal': r['journal'], 'title': r['title'], 'pmcid': r['pmcid'],
                        'study_type': g['study_type'], 'design': '; '.join(g['design']),
                        'country': g['country'], 'ancestry': g['ancestry'], 'japanese': g['japanese'],
                        'east_asian': g['east_asian'], 'n_cases': g['n_cases'], 'n_controls': g['n_controls'],
                        'control_type': g['control_type'], 'genes': '; '.join(x['as_written'] for x in g['genes']),
                        'gene_symbols': '; '.join(dict.fromkeys(s for s in sym if s)),
                        'variants': '; '.join(g['variants']), 'replication': g['replication'],
                        **{f'm_{k}': v for k, v in g['methods'].items()}, 'methods_quote': g['methods_quote'],
                        'agreement': why})
        # an adjudicated record: a finding both readers reported is 'both', the others 'adjudicated'
        both = set.intersection(*[{finding_key(f) for f in x.get('findings') or []} for x in rs]) \
            if why == 'adjudicated' and None not in rs else set()
        for f in g['findings']:
            ag = why if why != 'adjudicated' else ('both' if finding_key(f) in both else 'adjudicated')
            chk = result_check(f, P.curation['result_check'])
            if chk:
                review.append({'pmid': r['pmid'], 'pass': 'abstract', 'item': 'result check',
                               'final': f'{f.get("symbol") or f.get("gene")}: {f.get("result")}', 'decision': chk})
            findings.append({'pmid': r['pmid'], 'gene': f.get('gene') or '', 'symbol': f.get('symbol') or '',
                             'variant': f.get('variant') or '',
                             'result': f.get('result') or '', 'effect': f.get('effect') or '',
                             'p_value': f.get('p_value') or '', 'quote': f['quote'], 'agreement': ag,
                             'result_check': chk})
    # a record without a final reading is listed too: not read (no reader answered: an offline run without a cached
    # answer, or no answer), or left out (the readers differed and the adjudicator gave no answer)
    offline = a.offline.lower() in ('true', '1', 'yes')
    not_read = 0
    for r in recs:
        if r['pmid'] in finals:
            continue
        none = all(x is None for x in readings[r['pmid']])
        some = any(x is None for x in readings[r['pmid']])
        not_read += none
        screen.append({'pmid': r['pmid'], 'year': r['year'], 'title': r['title'], 'include': False, 'study_type': '',
                       'stage': 'not read' if none else 'model',
                       'reason': ('not read: no reader answered' + (' (this run calls no language model, and no '
                                                                     'answer is cached)' if offline else ''))
                       if none else 'not read by both readers (this run calls no language model, and one reader\'s '
                       'answer is not cached)' if offline and some else
                       'left out: a reader gave no answer, and the adjudicator gave none' if some else
                       'left out: the readers differed and the adjudicator gave no answer',
                       'agreement': ''})
    for x in logs:
        x.setdefault('pass', 'abstract')
    S = pd.DataFrame(studies, columns=STUDY_COLUMNS) if not studies else pd.DataFrame(studies)
    for k in STUDY_COLUMNS:
        if k not in S:
            S[k] = ''
    if a.fulltext_list:
        Path(a.fulltext_list).write_text('\n'.join(f'{r.pmid}={r.pmcid}' for r in S.itertuples() if r.pmcid) + '\n')

    # ── pass B: methods from the open-access Methods section, else from the abstract ──
    mrows = []
    if a.fulltext:
        ft = load_json(a.fulltext)
        todo = [(p, v['methods'], 'full text') for p, v in ft.items() if v['methods'] and p in set(S.pmid)]
        done = {p for p, _, _ in todo}
        todo += [(p, f'{by_pmid[p]["title"]}\n{by_pmid[p]["abstract"]}', 'abstract') for p in S.pmid if p not in done]
        head = {'full text': 'Methods', 'abstract': 'Title and abstract'}
        system = {'full text': SYS_B, 'abstract': SYS_B_ABS}
        mans = read_all(readers, [(system[src], f'PMID: {p}\n{head[src]}:\n{t}') for p, t, src in todo])
        for src in ('full text', 'abstract'):
            part = [(k, p, t) for k, (p, t, s) in enumerate(todo) if s == src]
            items = [(p, f'PMID: {p}\n{head[src]}:\n{t}',
                      [None if mans[ri][k] is None else method_fields(mans[ri][k], t, src, p, logs, rd.name)
                       for ri, rd in enumerate(readers)]) for k, p, t in part]
            text_of = {p: t for _k, p, t in part}
            fin, rv = reconcile(readers, adj, items, system[src], method_differences,
                                lambda ans, key, _t: method_fields(ans, text_of[key], src, key, logs, 'adjudicator'))
            agreed = {p: {f['field'] for f in rs[0]['fields']} - {d[0][9:] for d in method_differences(*rs)}
                      for p, _t, rs in items if None not in rs and len(rs) == 2}
            for row in rv:
                row.update({'pmid': row.pop('key'), 'pass': src})
                f = row['item'][9:]
                q = [x['quote'] for x in fin[row['pmid']][0]['fields'] if x['field'] == f] if row['pmid'] in fin else []
                row['final'] = (' | '.join(q) or 'not stated')[:200]
            review += rv
            for p, (reading, why) in fin.items():
                mrows += [{'pmid': p, 'field': f['field'], 'value': f['value'], 'quote': f['quote'], 'source': src,
                           'agreement': why if why != 'adjudicated' else
                           ('both' if f['field'] in agreed.get(p, set()) else 'adjudicated')}
                          for f in reading['fields']]
        S['methods_source'] = S.pmid.map(lambda p: 'full text' if p in done else 'abstract')
    else:
        S['methods_source'] = ''
    write_tsv(pd.DataFrame(pre + screen, columns=['pmid', 'year', 'title', 'include', 'study_type', 'stage', 'reason',
                                                  'agreement']), out / 'screening.tsv')
    write_tsv(S, out / 'studies.tsv')
    F = pd.DataFrame(findings, columns=['pmid', 'gene', 'symbol', 'variant', 'result', 'effect', 'p_value', 'quote',
                                        'agreement', 'result_check'])
    # official symbols: kept only when Open Targets knows them (exact name match), with the Ensembl ID
    syms = sorted({s for s in list(F.symbol) + [x for v in S.get('gene_symbols', []) for x in str(v).split('; ')] if s})
    if syms and a.expect_release:
        rel = release(a.api)['data_version']
        assert rel == a.expect_release, f'Open Targets serves {rel}, the snapshot is {a.expect_release}'
    ens = map_symbols(a.api, syms) if syms else {}
    for s in syms:
        if s not in ens:
            logs.append({'pmid': '', 'pass': 'symbols', 'item': 'symbol not in Open Targets', 'value': s,
                         'reader': ''})
    F['symbol'] = F.symbol.where(F.symbol.isin(ens.keys()), '')
    F['ensembl_id'] = F.symbol.map(ens).fillna('')
    if len(S):
        S['gene_symbols'] = S.gene_symbols.map(lambda v: '; '.join(x for x in str(v).split('; ') if x in ens))
    pw = pathways_of(a.api, {s: ens[s] for s in set(F.symbol) if s}) if ens and len(F) else {}
    write_tsv(F, out / 'findings.tsv')
    write_tsv(pd.DataFrame([{'symbol': s, 'ensembl_id': ens.get(s, ''), 'top_level_pathways': '; '.join(v)}
                            for s, v in sorted(pw.items())], columns=['symbol', 'ensembl_id', 'top_level_pathways']),
              out / 'gene_pathways.tsv')
    write_tsv(pd.DataFrame(mrows, columns=['pmid', 'field', 'value', 'quote', 'source', 'agreement']),
              out / 'methods.tsv')
    write_tsv(pd.DataFrame(logs, columns=['pmid', 'pass', 'item', 'value', 'reader']), out / 'extract_log.tsv')
    write_tsv(pd.DataFrame(review, columns=['pmid', 'pass', 'item'] + names + ['final', 'decision']),
              out / 'review.tsv')

    def tally(df, col='agreement'):
        return {k: int(v) for k, v in df[col].value_counts().items()} if len(df) else {}
    # inclusion: adjudicated only when the readers differed on it (or one gave no answer); a record adjudicated for
    # another item was included by both
    on_inclusion = {x['pmid'] for x in review if x['item'] in ('inclusion', 'reading')}
    Sc = pd.DataFrame([{'agreement': 'adjudicated' if x['pmid'] in on_inclusion else 'both'}
                       for x in screen if x['agreement'] in ('both', 'adjudicated')], columns=['agreement'])
    agreement = {'inclusion': tally(Sc), 'findings': tally(F),
                 'methods fields': tally(pd.DataFrame(mrows, columns=['agreement']))} if two else {}
    agreement = {k: {'both': v.get('both', 0), 'adjudicated': v.get('adjudicated', 0)} for k, v in agreement.items()}
    unanswered = [x for x in logs if x.get('item') == 'record' and x.get('value') == 'no answer']
    stats = {'records': len(recs) + len(pre), 'prescreened_out': len(pre), 'read_by_model': len(recs) - not_read,
             'not_read': not_read, 'offline': offline,
             'included': len(S), 'findings': len(findings), 'methods_fields': len(mrows),
             'dropped_items': len(logs) - len(unanswered),
             'disputed_items': sum(x['item'] != 'result check' for x in review), 'agreement': agreement,
             'result_flags': int(F.result_check.ne('').sum()) if len(F) else 0,
             'readers': [r.stats() for r in readers], 'adjudicator': adj.stats() if adj else None,
             'llm_calls': sum(r.calls for r in readers) + (adj.calls if adj else 0),
             'llm_cached': sum(r.cached for r in readers) + (adj.cached if adj else 0),
             'model': '; '.join(r.describe() for r in readers)}
    (out / 'extract.json').write_text(json.dumps(stats, indent=1))
    log(f'{len(recs) + len(pre)} records, {len(pre)} pre-screened out, {len(recs)} read: {len(S)} genetic studies, '
        f'{len(findings)} findings, {len(mrows)} methods fields; {len(logs)} items dropped by grounding; '
        f'{stats["disputed_items"]} disputed, {stats["result_flags"]} result flags; agreement {agreement}; '
        + '; '.join(f'{x["reader"]} {x["calls"]} calls, {x["cached"]} cached' for x in stats['readers'])
        + (f'; adjudicator {adj.calls} calls, {adj.cached} cached' if adj else ''))


def _final_value(item, g):
    """The adjudicated value of a disputed item of pass A."""
    if item in ('n_cases', 'n_controls', 'japanese', 'east_asian'):
        return str(g.get(item))
    if item == 'finding':
        return '; '.join(f'{k[0]}: {k[1]}' for k in sorted({finding_key(f) for f in g.get('findings') or []}))
    return 'included' if g.get('include') else 'excluded'


if __name__ == '__main__':
    main()
