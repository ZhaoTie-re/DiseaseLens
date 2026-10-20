#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The generated summaries of the knowledge base, each written by an
#           isolated language model (llm.py, the first reader) from the source
#           text it is given and checked against it:
#             mechanism.<theme>   per Reactome top-level pathway with >= --min-targets
#                                 text-mined targets of the disease: what the
#                                 sentences of the --top-targets best-scored
#                                 targets say
#             methods             how the genetic studies of the disease were done
#                                 (designs, populations, sample sizes, platforms,
#                                 analyses, limits), from the extracted study and
#                                 methods tables and their tallies
#             findings            what the genetic studies found (replicated or
#                                 single-study, population-specific, clonal)
#           Grounding: every [PMID] cited must be among the summary's inputs, and
#           every gene symbol of the input universe named in the text must be one
#           of its inputs. Second reader: the other model reads each summary with
#           its source material and lists the claims the material does not
#           support; the adjudicator rules on each, and a summary with a confirmed
#           unsupported claim is withheld. A summary that fails is kept with
#           grounded = false and left out of the knowledge base. The inputs are
#           the curated tables: only targets with a curated text-mined
#           publication enter a theme.
#           --offline: cached answers only; a summary whose inputs changed since it
#           was cached is not generated (problems: "not generated ..."). Summaries
#           are shown in American spelling (citations and gene symbols unchanged).
# Inputs  : the tidy tables, the profile, the model configuration
# Outputs : summaries.tsv (kind, key, title, summary, grounded, problems, check, inputs),
#           summaries.json
# ---------------------------------------------------------------------------
import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from llm import load_llm  # noqa: E402
from ot_common import legacy_path, load_profile, log, num, read_tsv, us_english, write_tsv  # noqa: E402

RULES = '''Use ONLY the material given; never add outside knowledge, numbers or references. After every claim, cite
the PubMed IDs it rests on in square brackets, e.g. [12345678] or [12345678, 23456789]. Name only genes that occur in
the material. Write plain, precise scientific English. Answer with one JSON object: {"summary": "<text>"}.'''
CHECK = '''You check a summary against the material it was written from. Read every claim of the summary (a sentence or
clause that states a fact) and list the claims the material does not support: a fact, number, gene or citation that is
not in the material, or that the material contradicts. Do not list a claim only because it is worded differently.
Answer with one JSON object: {"unsupported": [{"claim": "<the claim as written>", "reason": "<at most 20 words>"}]};
an empty list when every claim is supported.'''
RULE = '''You rule on claims of a summary that a reader found unsupported by the material the summary was written from.
For each claim, decide from the material alone whether it is supported. Answer with one JSON object:
{"rulings": [{"claim": "<the claim as given>", "supported": true or false, "reason": "<at most 20 words>"}]}.'''


def systems(P):
    """The system prompts of the three kinds of summary (for CTEPH, those of its first build, character for
    character)."""
    d = f'{P.name} ({P.short}). '
    return ('You summarise published findings on the mechanisms of ' + d + RULES,
            'You summarise the METHODS of the genetic studies of ' + d + RULES,
            'You summarise the FINDINGS of the genetic studies of ' + d + RULES)


def dehyphen(s):
    """VCAM-1 -> VCAM1, IL-1β -> IL1β: the official symbol and the written name compare equal."""
    return re.sub(r'(?<=[A-Za-z])-(?=\d)', '', s)


def check(summary, pmids, source, universe, synonyms):
    problems = []
    cited = set(re.findall(r'\b(\d{6,9})\b', ' '.join(re.findall(r'\[([^\]]+)\]', summary))))
    bad = cited - set(pmids)
    if bad:
        problems.append(f'cites PMIDs outside its inputs: {sorted(bad)}')
    if not cited:
        problems.append('cites nothing')
    word = lambda g, s: re.search(rf'(?<![A-Za-z0-9-]){re.escape(dehyphen(g))}(?![A-Za-z0-9-])', dehyphen(s))  # noqa: E731
    # a gene is grounded when its symbol, or a word Open Targets text mining matched to it, is in the source
    extra = {g for g in universe if word(g, summary) and not word(g, source)
             and not any(w and dehyphen(w).lower() in dehyphen(source).lower() for w in synonyms.get(g, ()))}
    if extra:
        problems.append(f'names genes absent from its source text: {sorted(extra)}')
    return not problems, '; '.join(problems)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tables', required=True)
    ap.add_argument('--profile', required=True, help='profile.json of the disease')
    ap.add_argument('--llm', required=True, help='llm.json: the readers and the adjudicator')
    ap.add_argument('--min-targets', type=int, default=3)
    ap.add_argument('--top-targets', type=int, default=12)
    ap.add_argument('--max-themes', type=int, default=0, help='smoke run: the N largest themes only')
    ap.add_argument('--offline', default='false', help='true: cached answers only, never call the model')
    ap.add_argument('--out', default='summaries.tsv')
    ap.add_argument('--out-json', default='summaries.json')
    a = ap.parse_args()
    P = load_profile(a.profile)
    short = P.short
    SYS_MECH, SYS_METH, SYS_FIND = systems(P)
    T = Path(a.tables)
    tab = lambda name, cols: read_tsv(legacy_path(T, f'disease_{name}.tsv', f'cteph_{name}.tsv'), cols)  # noqa: E731
    targets, lit = tab('targets', ['symbol', 'ensembl_id', 'ds_europepmc', 'themes']), read_tsv(T / 'literature.tsv')
    studies = tab('studies', ['pmid', 'year', 'first_author', 'design', 'country', 'ancestry', 'japanese', 'n_cases',
                              'n_controls', 'control_type', 'gene_symbols', 'genes', 'replication', 'methods_source'])
    findings = tab('findings', ['pmid', 'gene', 'symbol', 'variant', 'result', 'effect', 'p_value', 'quote'])
    mft = tab('methods', ['pmid', 'field', 'value'])
    universe = set(targets.symbol) | set(findings.symbol) | {s for v in studies.gene_symbols for s in v.split('; ') if s}
    # the gene symbols checked: of three characters or more (shorter ones are words too), and not a name of the
    # disease (PAH is a gene symbol and the abbreviation of pulmonary arterial hypertension)
    own = {w.upper() for w in [P.short, *P.short.split('/'), *re.findall(r'\b[A-Z][A-Z0-9]{1,10}\b', P.d.get(
        'pubmed_terms') or '')]}
    universe = {g for g in universe if g and len(g) >= 3 and g.upper() not in own}
    synonyms = {}
    for r in lit.itertuples():
        synonyms.setdefault(r.symbol, set()).update(w for w in r.matched_as.split(' | ') if len(w) >= 3)
    offline = a.offline.lower() in ('true', '1', 'yes')
    readers, adj = load_llm(a.llm, offline=offline)
    jobs, meta = [], []

    # ── mechanism themes ─────────────────────────────────────────────────────
    targets['lit'] = num(targets.ds_europepmc)
    lit['evidence_score'] = num(lit.evidence_score)
    tt = targets[targets.lit.notna() & targets.ensembl_id.isin(set(lit.ensembl_id))].copy()
    themes = sorted({t for v in tt.themes for t in v.split('; ') if t})
    if a.max_themes:
        size = {th: tt.themes.str.split('; ').map(lambda xs: th in xs).sum() for th in themes}
        themes = sorted(themes, key=lambda th: -size[th])[:a.max_themes]
    for th in themes:
        sub = tt[tt.themes.str.split('; ').map(lambda xs: th in xs)].sort_values('lit', ascending=False)
        if len(sub) < a.min_targets:
            continue
        lines, pm, gs = [], set(), []
        for t in sub.head(a.top_targets).itertuples():
            L = lit[lit.ensembl_id == t.ensembl_id].sort_values('evidence_score', ascending=False).head(2)
            for r in L.itertuples():
                s = r.sentences.split(' || ')[:2] if r.sentences else [r.title]
                lines += [f'{t.symbol} [{r.pmid}]: {x}' for x in s if x]
                pm.add(r.pmid)
            gs.append(t.symbol)
        prompt = (f'Theme (Reactome top-level pathway): {th}\nTargets: {", ".join(gs)}\n\nSentences from the '
                  'literature (gene [PMID]: sentence):\n' + '\n'.join(lines) +
                  '\n\nIn at most 150 words (3 to 6 sentences), state what this material says about the role of '
                  f'this theme in {short}.')
        jobs.append((SYS_MECH, prompt))
        meta.append({'kind': 'mechanism', 'key': th, 'title': th, 'pmids': sorted(pm), 'genes': gs,
                     'n_targets': len(sub), 'source': prompt})

    # ── methods and findings of the genetic studies ──────────────────────────
    if len(studies) and not a.max_themes:
        st = studies.copy()
        tally = {'studies': len(st),
                 'designs': st.design.str.split('; ').explode().value_counts().to_dict(),
                 'populations': st.country.replace('', 'not stated').value_counts().to_dict(),
                 'japanese': int((st.japanese == 'True').sum()),
                 'cases_median': float(num(st.n_cases).median()) if num(st.n_cases).notna().any() else None,
                 'cases_range': [float(num(st.n_cases).min()), float(num(st.n_cases).max())]
                 if num(st.n_cases).notna().any() else None,
                 'methods_from_full_text': int((st.methods_source == 'full text').sum()),
                 'methods_from_abstract_only': int((st.methods_source == 'abstract').sum())}
        rows = []
        for r in st.itertuples():
            m = '; '.join(f'{k}: {getattr(r, "m_" + k)}' for k in ('genotyping', 'statistics', 'covariates',
                                                                   'multiple_testing') if getattr(r, 'm_' + k, ''))
            ft = mft[mft.pmid == r.pmid]
            ftm = '; '.join(f'{x.field}: {x.value}' for x in ft.itertuples())
            rows.append(f'[{r.pmid}] {r.year} {r.first_author}: design {r.design}; population {r.country or "?"} '
                        f'{r.ancestry}; cases {r.n_cases or "?"}, controls {r.n_controls or "?"} ({r.control_type}); '
                        f'genes {r.gene_symbols or r.genes}; replication {r.replication}. Methods (abstract): '
                        f'{m or "not stated"}. Methods ({r.methods_source or "none"}): {ftm or "none stated"}')
        pm = list(st.pmid)
        prompt = ('Tallies: ' + json.dumps(tally) + '\n\nStudies:\n' + '\n'.join(rows) +
                  f'\n\nWrite 3 to 5 short paragraphs on how these genetic studies of {short} were done: study designs '
                  'and how they changed over time; populations (Japanese and other East Asian studies first); '
                  'sample sizes and controls; genotyping / sequencing platforms and variant quality control; '
                  'statistical analyses, covariates, multiple-testing correction and replication; and the common '
                  'methodological limits the material shows. At most 350 words.')
        jobs.append((SYS_METH, prompt))
        meta.append({'kind': 'methods', 'key': 'methods', 'title': f'Methods of the genetic studies of {short}',
                     'source': prompt,
                     'pmids': pm, 'genes': sorted(({s for v in st.gene_symbols for s in v.split('; ') if s}
                                                   | set(findings.symbol)) - {''}),
                     'tally': tally})
        frows = []
        for f in findings.itertuples():
            s = st[st.pmid == f.pmid]
            pop = (s.country.iloc[0] or 'population not stated') if len(s) else ''
            frows.append(f'[{f.pmid}] ({pop}) {f.symbol or f.gene} {f.variant}: {f.result}; {f.effect} '
                         f'{f.p_value}. "{f.quote}"')
        prompt = ('Findings:\n' + '\n'.join(frows) + '\n\nWrite 3 to 5 short paragraphs: which genes or variants '
                  f'were reported as associated with {short} and by how many independent studies; which findings are '
                  'specific to a population (Japanese first); which concern acquired (somatic / clonal) variants; '
                  'which candidate genes showed no association; and how strong the evidence is overall. At most '
                  '350 words.')
        fg = sorted({f.symbol for f in findings.itertuples() if f.symbol})
        jobs.append((SYS_FIND, prompt))
        meta.append({'kind': 'findings', 'key': 'findings', 'title': f'Findings of the genetic studies of {short}',
                     'source': prompt,
                     'pmids': sorted(set(findings.pmid)), 'genes': fg})

    # ── written by the first reader, checked against the material ────────────
    answers = readers[0].map(jobs)
    out = []
    for m, ans in zip(meta, answers):
        src = m['source']
        if ans is None and offline:
            out.append({**m, 'summary': '', 'grounded': False, 'check': '',
                        'problems': 'not generated: no summary of these inputs is cached, and this run calls no '
                                    'language model'})
            continue
        text = (ans or {}).get('summary', '') if isinstance(ans, dict) else ''
        text = us_english(text)       # American spelling; citations and gene symbols are unchanged
        ok, problems = check(text, m['pmids'], src, universe, synonyms)
        out.append({**m, 'summary': text, 'grounded': ok, 'problems': problems,
                    'check': 'one reader' if len(readers) == 1 else ''})

    # ── the second reader lists unsupported claims; the adjudicator rules on them ──
    if len(readers) == 2:
        todo = [r for r in out if r['summary']]
        audits = readers[1].map([(CHECK, f'Material:\n{r["source"]}\n\n---\n\nSummary:\n{r["summary"]}')
                                 for r in todo])
        disputed = []
        for r, au in zip(todo, audits):
            # a summary is shown only when it was checked: an answer without a list of claims is no check
            claims = au.get('unsupported') if isinstance(au, dict) else None
            if not isinstance(claims, list):
                why = 'no answer from the second reader' if au is None else \
                    'the second reader\'s answer gives no list of unsupported claims'
                r['check'] = f'not checked ({why})'
                r['grounded'] = False
                r['problems'] = '; '.join(x for x in (r['problems'], f'not checked: {why}') if x)
            elif not claims:
                r['check'] = 'both'
            else:
                disputed.append((r, claims))
        rulings = adj.map([(RULE, f'Material:\n{r["source"]}\n\n---\n\nSummary:\n{r["summary"]}\n\n---\n\nClaims:\n'
                            + json.dumps(claims, ensure_ascii=False)) for r, claims in disputed]) if adj and \
            disputed else [None] * len(disputed)
        for (r, claims), ru in zip(disputed, rulings):
            # one ruling per disputed claim, each true or false; anything else is no ruling
            rs = ru.get('rulings') if isinstance(ru, dict) else None
            if not (isinstance(rs, list) and len(rs) == len(claims) and
                    all(isinstance(x, dict) and isinstance(x.get('supported'), bool) for x in rs)):
                r['check'] = f'{len(claims)} claims disputed, no ruling on each'
                r['grounded'] = False
                r['problems'] = '; '.join(x for x in (r['problems'], 'claims disputed by the second reader, not '
                                                      'ruled on') if x)
                continue
            upheld = [x for x in rs if x['supported'] is False]
            r['check'] = f'adjudicated: {len(upheld)} of {len(claims)} disputed claims unsupported'
            if upheld:
                r['grounded'] = False
                r['problems'] = '; '.join(x for x in (r['problems'], 'unsupported claims: ' + ' | '.join(
                    f'{str(x.get("claim") or "")[:120]} ({x.get("reason") or ""})' for x in upheld)) if x)
    for r in out:
        r.pop('source', None)
    df = pd.DataFrame([{k: (('; '.join(v) if isinstance(v, list) else json.dumps(v) if isinstance(v, dict) else v))
                        for k, v in r.items()} for r in out],
                      columns=['kind', 'key', 'title', 'pmids', 'genes', 'n_targets', 'tally', 'summary', 'grounded',
                               'problems', 'check'])
    write_tsv(df, a.out)
    calls = sum(r.calls for r in readers) + (adj.calls if adj else 0)
    cached = sum(r.cached for r in readers) + (adj.cached if adj else 0)
    missed = sum(r.missed for r in readers) + (adj.missed if adj else 0)
    Path(a.out_json).write_text(json.dumps({'summaries': out, 'llm_calls': calls, 'llm_cached': cached,
                                            'not_generated': missed, 'offline': offline,
                                            'readers': [r.stats() for r in readers],
                                            'adjudicator': adj.stats() if adj else None},
                                           indent=1, ensure_ascii=False))
    log(f'{len(out)} summaries ({sum(r["grounded"] for r in out)} grounded); {calls} model calls, '
        f'{cached} cached, {missed} not generated (offline)')


if __name__ == '__main__':
    main()
