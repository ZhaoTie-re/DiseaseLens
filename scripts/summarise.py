#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The generated summaries of the knowledge base, each written by an
#           isolated language model (llm.py) from the source text it is given
#           and checked against it:
#             mechanism.<theme>   per Reactome top-level pathway with >= --min-targets
#                                 text-mined CTEPH targets: what the sentences of
#                                 the --top-targets best-scored targets say
#             methods             how the genetic studies of CTEPH were done
#                                 (designs, populations, sample sizes, platforms,
#                                 analyses, limits), from the extracted study and
#                                 methods tables and their tallies
#             findings            what the genetic studies of CTEPH found (replicated
#                                 or single-study, population-specific, clonal)
#           Grounding: every [PMID] cited must be among the summary's inputs, and
#           every gene symbol of the input universe named in the text must be one
#           of its inputs; a summary that fails is kept with grounded = false and
#           left out of the knowledge base. The inputs are the curated tables:
#           only targets with a curated text-mined publication enter a theme.
#           --offline: cached answers only; a summary whose inputs changed since it
#           was cached is not generated (problems: "not generated ..."). Summaries
#           are shown in American spelling (citations and gene symbols unchanged).
# Inputs  : the tidy tables
# Outputs : summaries.tsv (kind, key, title, summary, grounded, problems, inputs),
#           summaries.json
# ---------------------------------------------------------------------------
import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from llm import LLM  # noqa: E402
from ot_common import log, num, read_tsv, us_english, write_tsv  # noqa: E402

RULES = '''Use ONLY the material given; never add outside knowledge, numbers or references. After every claim, cite
the PubMed IDs it rests on in square brackets, e.g. [12345678] or [12345678, 23456789]. Name only genes that occur in
the material. Write plain, precise scientific English. Answer with one JSON object: {"summary": "<text>"}.'''
SYS_MECH = 'You summarise published findings on the mechanisms of chronic thromboembolic pulmonary hypertension ' \
           '(CTEPH). ' + RULES
SYS_METH = 'You summarise the METHODS of the genetic studies of chronic thromboembolic pulmonary hypertension ' \
           '(CTEPH). ' + RULES
SYS_FIND = 'You summarise the FINDINGS of the genetic studies of chronic thromboembolic pulmonary hypertension ' \
           '(CTEPH). ' + RULES


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
    ap.add_argument('--cache-dir', required=True)
    ap.add_argument('--workdir', required=True)
    ap.add_argument('--model', default='')
    ap.add_argument('--effort', default='')
    ap.add_argument('--parallel', type=int, default=3)
    ap.add_argument('--min-targets', type=int, default=3)
    ap.add_argument('--top-targets', type=int, default=12)
    ap.add_argument('--claude', default='claude')
    ap.add_argument('--max-themes', type=int, default=0, help='smoke run: the N largest themes only')
    ap.add_argument('--offline', default='false', help='true: cached answers only, never call the model')
    ap.add_argument('--out', default='summaries.tsv')
    ap.add_argument('--out-json', default='summaries.json')
    a = ap.parse_args()
    T = Path(a.tables)
    targets, lit = read_tsv(T / 'cteph_targets.tsv'), read_tsv(T / 'literature.tsv')
    studies, findings = read_tsv(T / 'cteph_studies.tsv'), read_tsv(T / 'cteph_findings.tsv')
    mft = read_tsv(T / 'cteph_methods.tsv')
    universe = set(targets.symbol) | set(findings.symbol) | {s for v in studies.gene_symbols for s in v.split('; ') if s}
    universe = {g for g in universe if g and len(g) >= 3}
    synonyms = {}
    for r in lit.itertuples():
        synonyms.setdefault(r.symbol, set()).update(w for w in r.matched_as.split(' | ') if len(w) >= 3)
    offline = a.offline.lower() in ('true', '1', 'yes')
    llm = LLM(a.cache_dir, a.workdir, a.model, a.effort, parallel=a.parallel, claude=a.claude, offline=offline)
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
                  'this theme in CTEPH.')
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
                  '\n\nWrite 3 to 5 short paragraphs on how these genetic studies of CTEPH were done: study designs '
                  'and how they changed over time; populations (Japanese and other East Asian studies first); '
                  'sample sizes and controls; genotyping / sequencing platforms and variant quality control; '
                  'statistical analyses, covariates, multiple-testing correction and replication; and the common '
                  'methodological limits the material shows. At most 350 words.')
        jobs.append((SYS_METH, prompt))
        meta.append({'kind': 'methods', 'key': 'methods', 'title': 'Methods of the genetic studies of CTEPH',
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
                  'were reported as associated with CTEPH and by how many independent studies; which findings are '
                  'specific to a population (Japanese first); which concern acquired (somatic / clonal) variants; '
                  'which candidate genes showed no association; and how strong the evidence is overall. At most '
                  '350 words.')
        fg = sorted({f.symbol for f in findings.itertuples() if f.symbol})
        jobs.append((SYS_FIND, prompt))
        meta.append({'kind': 'findings', 'key': 'findings', 'title': 'Findings of the genetic studies of CTEPH',
                     'source': prompt,
                     'pmids': sorted(set(findings.pmid)), 'genes': fg})

    answers = llm.map(jobs)
    out = []
    for m, ans in zip(meta, answers):
        src = m.pop('source')
        if ans is None and offline:
            out.append({**m, 'summary': '', 'grounded': False,
                        'problems': 'not generated: its inputs changed since the cached summary, and this run calls '
                                    'no model'})
            continue
        text = (ans or {}).get('summary', '') if isinstance(ans, dict) else ''
        text = us_english(text)       # American spelling; citations and gene symbols are unchanged
        ok, problems = check(text, m['pmids'], src, universe, synonyms)
        out.append({**m, 'summary': text, 'grounded': ok, 'problems': problems})
    df = pd.DataFrame([{k: (('; '.join(v) if isinstance(v, list) else json.dumps(v) if isinstance(v, dict) else v))
                        for k, v in r.items()} for r in out])
    write_tsv(df, a.out)
    Path(a.out_json).write_text(json.dumps({'summaries': out, 'llm_calls': llm.calls, 'llm_cached': llm.cached,
                                            'not_generated': llm.missed, 'offline': offline},
                                           indent=1, ensure_ascii=False))
    log(f'{len(out)} summaries ({sum(r["grounded"] for r in out)} grounded); {llm.calls} model calls, '
        f'{llm.cached} cached, {llm.missed} not generated (offline)')


if __name__ == '__main__':
    main()
