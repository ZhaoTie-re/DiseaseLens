#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The genetic studies of CTEPH in PubMed, screened and extracted by an
#           isolated language model (llm.py), then checked against their source
#           text. Two passes:
#             A  abstracts, --batch per call: is the record an original genetic
#                study (or case report) of CTEPH? If so: design, population,
#                sample sizes, genes, variants, findings, replication and the
#                methods the abstract states
#             B  the methods of every included record (platform, genotyping /
#                sequencing, quality control, statistics, covariates, multiple
#                testing, replication): from the Methods section of its open-
#                access full text, or else from its abstract
#           Grounding (an item that fails is dropped and logged, never kept):
#             every gene, variant and finding quote must occur in the abstract;
#             sample sizes must occur as numbers; 'Japanese' / 'East Asian' must
#             be stated; every methods value must come with a quote that occurs
#             in the text it was read from.
# Inputs  : the PubMed search records, the Methods texts (pubmed.py --mode fulltext)
# Outputs : screening.tsv, studies.tsv, findings.tsv, methods.tsv (field, value, quote, source),
#           extract_log.tsv, extract.json (counts); with --fulltext-list only:
#           fulltext_ids.txt (PMID=PMCID of the included records with a PMCID)
# ---------------------------------------------------------------------------
import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from llm import LLM, in_text, numbers_in  # noqa: E402
from ot_common import gql, load_json, log, write_tsv  # noqa: E402

MAPIDS = '''query($q: [String!]!) { mapIds(queryTerms: $q, entityNames: ["target"]) {
  mappings { term hits { id name } } } }'''


def map_symbols(api, symbols):
    """{symbol: Ensembl ID} for the symbols Open Targets knows by that exact name."""
    out = {}
    for i in range(0, len(symbols), 100):
        for m in gql(api, MAPIDS, {'q': symbols[i:i + 100]})['mapIds']['mappings']:
            hit = next((h for h in m['hits'] if h['name'] == m['term']), None)
            if hit:
                out[m['term']] = hit['id']
    return out


# deterministic pre-screen (before any model call): publication types that are never original data, and the
# terms an original genetic study states in its title or abstract (kept 16 of 16 known CTEPH genetic studies,
# 98 of 471 records, 2026-10-16)
EXCLUDED_TYPES = {'Review', 'Systematic Review', 'Practice Guideline', 'Guideline', 'Editorial', 'Comment', 'News',
                  'Consensus Development Conference'}
GENETIC_TERMS = re.compile(r'polymorphism|genotyp|allele|mutation|sequencing|exome|genome-wide|GWAS|\bHLA\b|'
                           r'\bSNPs?\b|\brs\d{3,}|germline|somatic|clonal hematopoiesis|clonal haematopoiesis|CHIP|'
                           r'heritab|Mendelian randomi|copy number|pathogenic variant|rare variant|missense|'
                           r'haplotype|linkage', re.I)
DESIGNS = ['GWAS', 'candidate gene association', 'whole-exome sequencing', 'whole-genome sequencing',
           'targeted gene panel', 'HLA typing', 'copy-number analysis', 'clonal haematopoiesis / somatic',
           'mitochondrial DNA', 'Mendelian randomization', 'family / linkage', 'other']
SYSTEM_A = f'''You read PubMed records and extract facts about the GENETICS of chronic thromboembolic pulmonary
hypertension (CTEPH). Use ONLY the text given; never add outside knowledge. Copy gene symbols, variant names, numbers
and quotes exactly as they are written in the record. Answer with JSON only (no prose).

For every record return one object:
{{"pmid": "<PMID>",
  "include": true or false,   // true only for an ORIGINAL study or case report that reports genetic data (germline
                              // or somatic variants, genotypes, HLA alleles, sequencing, copy number, Mendelian
                              // randomization) in people with CTEPH. Reviews, guidelines, gene-EXPRESSION-only
                              // studies, animal-only studies and studies without CTEPH patients: false.
  "reason": "<at most 15 words>",
  "study_type": "original research" | "case report" | "review" | "guideline" | "other",
  "design": [one or more of {json.dumps(DESIGNS)}],
  "population": {{"country": "", "ancestry": "", "japanese": true/false, "east_asian": true/false}},
  "n_cases": <CTEPH patients genotyped, integer or null>, "n_controls": <integer or null>, "control_type": "",
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
SYSTEM_B = '''You read the Methods section of a genetic study of chronic thromboembolic pulmonary hypertension and
extract how the study was done. Use ONLY the text given; copy values exactly. Answer with one JSON object (no prose):
{"fields": [{"field": one of "cohort / recruitment", "diagnosis of CTEPH", "controls", "DNA source",
             "genotyping platform", "sequencing", "HLA typing", "variant calling / quality control",
             "imputation", "statistical test", "covariates", "multiple testing", "replication", "power",
             "software",
             "value": "<short summary in the text's own words>",
             "quote": "<an exact passage of the text that states it>"}]}
List only fields the text states; omit the others. Quotes are copied whole, never shortened with '...'.'''
# the same fields from the abstract, for the included studies without an open-access Methods section (the user,
# 2026-10-16: "for those whose full text cannot be fetched, fall back to the abstract")
SYSTEM_B_ABSTRACT = SYSTEM_B.replace('the Methods section of a genetic study',
                                     'the title and abstract of a genetic study').replace(
    'extract how the study was done.', 'extract how the study was done (no full text is available).')


def batches(xs, n):
    return [xs[i:i + n] for i in range(0, len(xs), n)]


def prompt_a(recs):
    parts = [f'PMID: {r["pmid"]}\nTitle: {r["title"]}\nPublication types: {"; ".join(r["pub_types"])}\n'
             f'Abstract: {r["abstract"] or "(no abstract)"}' for r in recs]
    return 'Records:\n\n' + '\n\n---\n\n'.join(parts)


def ground_a(ans, r, log_rows):
    """The answer for one record, with every ungrounded item removed (and logged)."""
    text = r['title'] + ' ' + r['abstract']
    pm = r['pmid']

    def drop(what, value):
        log_rows.append({'pmid': pm, 'pass': 'abstract', 'item': what, 'value': str(value)[:200]})

    out = {k: ans.get(k) for k in ('include', 'reason', 'study_type')}
    if not ans.get('include'):
        return out
    out['design'] = [d for d in ans.get('design') or [] if d in DESIGNS] or ['other']
    pop = ans.get('population') or {}
    jp = bool(pop.get('japanese')) and (in_text('Japan', text) or in_text('Japanese', text))
    ea = bool(pop.get('east_asian')) and any(in_text(w, text) for w in
                                             ('Japan', 'Chin', 'Korea', 'Taiwan', 'East Asian', 'Asian'))
    if pop.get('japanese') and not jp:
        drop('japanese', 'not stated in the abstract')
    out.update({'country': pop.get('country') or '', 'ancestry': pop.get('ancestry') or '', 'japanese': jp,
                'east_asian': ea or jp, 'control_type': ans.get('control_type') or '',
                'replication': ans.get('replication') or 'not stated'})
    for k in ('n_cases', 'n_controls'):
        v = ans.get(k)
        if v not in (None, '') and not numbers_in(v, text):
            drop(k, v)
            v = None
        out[k] = v
    genes = []
    for g in ans.get('genes') or []:
        g = g if isinstance(g, dict) else {'as_written': g, 'symbol': ''}
        if in_text(g.get('as_written', ''), text) or (g.get('symbol') and in_text(g['symbol'], text)):
            genes.append(g)
        else:
            drop('gene', g.get('as_written', ''))
    out['genes'] = genes
    out['variants'] = [v for v in ans.get('variants') or [] if in_text(v, text) or drop('variant', v)]
    fs = []
    for f in ans.get('findings') or []:
        if not in_text(f.get('quote', ''), text):
            drop('finding', f.get('quote', ''))
            continue
        if f.get('gene') and not (in_text(f['gene'], text) or (f.get('symbol') and in_text(f['symbol'], text))):
            drop('finding gene', f['gene'])
            continue
        if f.get('variant') and not in_text(f['variant'], text):
            drop('finding variant', f['variant'])
            f['variant'] = ''
        fs.append(f)
    out['findings'] = fs
    m = ans.get('methods') or {}
    if m.get('quote') and in_text(m['quote'], text):
        out['methods'] = {k: m.get(k) or '' for k in ('genotyping', 'statistics', 'covariates', 'multiple_testing')}
        out['methods_quote'] = m['quote']
    else:
        if any(m.get(k) for k in ('genotyping', 'statistics', 'covariates', 'multiple_testing')):
            drop('methods', m.get('quote', ''))
        out['methods'], out['methods_quote'] = {}, ''
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--records', required=True, help='pubmed.py --mode search output')
    ap.add_argument('--fulltext', help='pubmed.py --mode fulltext output (pass B); omit for pass A only')
    ap.add_argument('--fulltext-list', help='pass A only: write PMID=PMCID of the included records here')
    ap.add_argument('--cache-dir', required=True)
    ap.add_argument('--workdir', required=True, help='scratch directory the model runs in')
    ap.add_argument('--model', default='')
    ap.add_argument('--effort', default='')
    ap.add_argument('--parallel', type=int, default=3)
    ap.add_argument('--batch', type=int, default=5)
    ap.add_argument('--limit', type=int, default=0, help='smoke run: the first N records only')
    ap.add_argument('--prescreen', default='true', help='drop reviews and records without a genetic-study term '
                    'before the model reads them')
    ap.add_argument('--claude', default='claude')
    ap.add_argument('--api', required=True, help='Open Targets GraphQL (symbol check)')
    ap.add_argument('--out-dir', default='.')
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    recs = load_json(a.records)['records']
    recs = [r for r in recs if r['abstract'] or r['title']]
    pre = []
    if a.prescreen.lower() in ('true', '1', 'yes'):
        keep = []
        for r in recs:
            types = set(r['pub_types']) & EXCLUDED_TYPES
            if types or not GENETIC_TERMS.search(r['title'] + ' ' + r['abstract']):
                pre.append({'pmid': r['pmid'], 'year': r['year'], 'title': r['title'], 'include': False,
                            'study_type': 'review' if types else '', 'stage': 'rule',
                            'reason': f'pre-screen: {", ".join(sorted(types))}' if types
                            else 'pre-screen: no genetic-study term in title or abstract'})
            else:
                keep.append(r)
        recs = keep
    if a.limit:
        recs = recs[:a.limit]
    llm = LLM(a.cache_dir, a.workdir, a.model, a.effort, parallel=a.parallel, claude=a.claude)

    # ── pass A: abstracts ────────────────────────────────────────────────────
    jobs = [(SYSTEM_A, prompt_a(b)) for b in batches(recs, a.batch)]
    answers = llm.map(jobs)
    logs, screen, studies, findings = [], [], [], []
    for b, ans in zip(batches(recs, a.batch), answers):
        ans = ans if isinstance(ans, list) else [ans] if isinstance(ans, dict) else []
        by = {str(x.get('pmid')): x for x in ans}
        for r in b:
            x = by.get(r['pmid'])
            if x is None:
                logs.append({'pmid': r['pmid'], 'pass': 'abstract', 'item': 'record', 'value': 'no answer'})
                continue
            g = ground_a(x, r, logs)
            screen.append({'pmid': r['pmid'], 'year': r['year'], 'title': r['title'], 'include': bool(g['include']),
                           'study_type': g.get('study_type') or '', 'stage': 'model', 'reason': g.get('reason') or ''})
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
                            **{f'm_{k}': v for k, v in g['methods'].items()}, 'methods_quote': g['methods_quote']})
            for f in g['findings']:
                findings.append({'pmid': r['pmid'], 'gene': f.get('gene') or '', 'symbol': f.get('symbol') or '',
                                 'variant': f.get('variant') or '',
                                 'result': f.get('result') or '', 'effect': f.get('effect') or '',
                                 'p_value': f.get('p_value') or '', 'quote': f['quote']})
    S = pd.DataFrame(studies)
    for k in ('m_genotyping', 'm_statistics', 'm_covariates', 'm_multiple_testing'):
        if k not in S:
            S[k] = ''
    if a.fulltext_list:
        Path(a.fulltext_list).write_text('\n'.join(f'{r.pmid}={r.pmcid}' for r in S.itertuples() if r.pmcid) + '\n')

    # ── pass B: methods from the open-access Methods section, else from the abstract ──
    mrows = []
    if a.fulltext:
        ft = load_json(a.fulltext)
        by = {r['pmid']: r for r in recs}
        todo = [(p, v['methods'], 'full text') for p, v in ft.items() if v['methods'] and p in set(S.pmid)]
        done = {p for p, _, _ in todo}
        todo += [(p, f'{by[p]["title"]}\n{by[p]["abstract"]}', 'abstract') for p in S.pmid if p not in done]
        jobs = [(SYSTEM_B, f'PMID: {p}\nMethods:\n{t}') if src == 'full text' else
                (SYSTEM_B_ABSTRACT, f'PMID: {p}\nTitle and abstract:\n{t}') for p, t, src in todo]
        for (p, t, src), ans in zip(todo, llm.map(jobs)):
            for f in (ans or {}).get('fields') or []:
                if not in_text(f.get('quote', ''), t):
                    logs.append({'pmid': p, 'pass': src, 'item': f.get('field', ''),
                                 'value': f.get('quote', '')[:200]})
                    continue
                mrows.append({'pmid': p, 'field': f.get('field', ''), 'value': f.get('value', ''), 'quote': f['quote'],
                              'source': src})
        S['methods_source'] = S.pmid.map(lambda p: 'full text' if p in done else 'abstract')
    else:
        S['methods_source'] = ''
    write_tsv(pd.DataFrame(pre + screen, columns=['pmid', 'year', 'title', 'include', 'study_type', 'stage', 'reason']),
              out / 'screening.tsv')
    write_tsv(S, out / 'studies.tsv')
    F = pd.DataFrame(findings, columns=['pmid', 'gene', 'symbol', 'variant', 'result', 'effect', 'p_value', 'quote'])
    # official symbols: kept only when Open Targets knows them (exact name match), with the Ensembl ID
    syms = sorted({s for s in list(F.symbol) + [x for v in S.get('gene_symbols', []) for x in str(v).split('; ')] if s})
    ens = map_symbols(a.api, syms) if syms else {}
    for s in syms:
        if s not in ens:
            logs.append({'pmid': '', 'pass': 'symbols', 'item': 'symbol not in Open Targets', 'value': s})
    F['symbol'] = F.symbol.where(F.symbol.isin(ens.keys()), '')
    F['ensembl_id'] = F.symbol.map(ens).fillna('')
    if len(S):
        S['gene_symbols'] = S.gene_symbols.map(lambda v: '; '.join(x for x in str(v).split('; ') if x in ens))
    write_tsv(F,
              out / 'findings.tsv')
    write_tsv(pd.DataFrame(mrows, columns=['pmid', 'field', 'value', 'quote', 'source']), out / 'methods.tsv')
    write_tsv(pd.DataFrame(logs, columns=['pmid', 'pass', 'item', 'value']), out / 'extract_log.tsv')
    stats = {'records': len(recs) + len(pre), 'prescreened_out': len(pre), 'read_by_model': len(recs),
             'included': len(S), 'findings': len(findings), 'methods_fields': len(mrows),
             'dropped_items': len(logs), 'llm_calls': llm.calls, 'llm_cached': llm.cached,
             'model': a.model or 'CLI default', 'effort': a.effort or 'CLI default'}
    (out / 'extract.json').write_text(json.dumps(stats, indent=1))
    log(f'{len(recs) + len(pre)} records, {len(pre)} pre-screened out, {len(recs)} read: {len(S)} genetic studies, {len(findings)} findings, {len(mrows)} methods fields; '
        f'{len(logs)} items dropped by grounding; {llm.calls} model calls, {llm.cached} cached')


if __name__ == '__main__':
    main()
