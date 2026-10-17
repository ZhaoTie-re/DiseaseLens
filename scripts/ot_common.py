#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Helpers shared by the cteph_knowledge_graph step scripts: the Open
#           Targets GraphQL client (retries with backoff, as in the CTEPH
#           project's Open Targets client), gzipped JSON
#           I/O, the TSV conventions of the project's integration tables, and
#           the gene themes and population labels of the CTEPH literature
#           (shared by the knowledge base and the workbook).
#           No paths, no project constants: every script receives its inputs
#           on the command line from main.nf.
# ---------------------------------------------------------------------------
import datetime as dt
import gzip
import json
import re
import time
import urllib.error
import urllib.request
from collections import OrderedDict
from pathlib import Path

import pandas as pd
import yaml

MISSING = '—'   # em dash: "not available / not applicable"
META = '{ meta { apiVersion { x y z } dataVersion { year month iteration } } }'
# Open Targets datasource ids -> the labels printed in the workbook and the knowledge base
DATASOURCES = {'europepmc': 'Literature', 'impc': 'Animal model', 'clinical_precedence': 'Clinical',
               'gwas_credible_sets': 'GWAS'}


# Open Targets / gnomAD population codes of a study's LD reference
ANCESTRY = {'nfe': 'European', 'fin': 'Finnish', 'eas': 'East Asian', 'afr': 'African', 'amr': 'admixed American',
            'sas': 'South Asian', 'csa': 'Central / South Asian', 'mid': 'Middle Eastern', 'asj': 'Ashkenazi Jewish'}
STAGE = {'APPROVAL': 'Approval', 'PHASE_4': 'Phase 4', 'PHASE_3': 'Phase 3', 'PHASE_2_3': 'Phase 2/3',
         'PHASE_2': 'Phase 2', 'PHASE_1_2': 'Phase 1/2', 'PHASE_1': 'Phase 1', 'EARLY_PHASE_1': 'Early phase 1',
         'PRECLINICAL': 'Preclinical', 'UNKNOWN': 'Unknown'}
CAVEATS = [
    'Open Targets associations are scores that aggregate evidence; they are not causal claims.',
    'Literature evidence is text mining (Europe PMC): a gene symbol can be mapped to the wrong gene (for example, '
    '"TF" for tissue factor is mapped to transferrin, TF). The words a target was matched as are listed beside it.',
    'Locus-to-gene (L2G) scores are predictions of the causal gene at a GWAS locus, not proof.',
    'GWAS P values and effects are those of each credible set\'s lead variant.',
    'Mouse-model evidence links a knockout\'s phenotypes to the disease\'s clinical features; it is not a model of '
    'the disease itself.',
    'Related-disease studies overlap in their samples (for example, FinnGen and UK Biobank endpoints); genes are '
    'therefore counted by independent publications, not by studies.',
]
SUP = str.maketrans('0123456789-', '⁰¹²³⁴⁵⁶⁷⁸⁹⁻')


def stage_label(s):
    s = str(s or '')
    return STAGE.get(s, s.replace('_', ' ').capitalize())


def phase_label(phase, stage=''):
    """A trial's phase as registered ('PHASE2/PHASE3', 'EARLY_PHASE1', 'NA') -> 'Phase 2/3' ...; the harmonized
    stage when the registry gives none."""
    s = str(phase or '').upper()
    digits = re.findall(r'PHASE_?(\d)', s)
    if digits:
        return ('Early phase ' if s.startswith('EARLY') else 'Phase ') + '/'.join(digits)
    if s in ('NA', 'N/A', 'NOT_APPLICABLE'):
        return 'phase not applicable'
    return stage_label(stage) if stage and str(stage).upper() != 'UNKNOWN' else 'phase not given'


def ancestry_label(s):
    """'nfe 70%; eas 30%' -> 'European 70%, East Asian 30%'."""
    out = []
    for part in [x.strip() for x in str(s or '').split(';') if x.strip()]:
        code, _, pct = part.partition(' ')
        out.append(f'{ANCESTRY.get(code, code)} {pct}'.strip())
    return ', '.join(out)


def p_unicode(mant, exp):
    """P from Open Targets' mantissa / exponent as 7.0 × 10⁻²⁷ (no float underflow)."""
    if mant is None or exp is None or pd.isna(mant) or pd.isna(exp):
        return ''
    m, e = float(mant), int(exp)
    if m >= 9.95:
        m, e = m / 10, e + 1
    return f'{m:.1f} × 10{str(e).translate(SUP)}'


def pub_label(pid):
    """A literature ID as printed: PubMed IDs as 'PMID n'; Europe PMC preprints (PPR...) and full texts (PMC...)
    by their Europe PMC ID."""
    pid = str(pid or '')
    return f'PMID {pid}' if pid.isdigit() else f'Europe PMC {pid}' if pid else ''


def pub_url(pid):
    pid = str(pid or '')
    if pid.isdigit():
        return f'https://pubmed.ncbi.nlm.nih.gov/{pid}/'
    m = re.match(r'(PPR|PMC)\d+$', pid)
    return f'https://europepmc.org/article/{m.group(1)}/{pid}' if m else ''


def variant_label(vid):
    """Open Targets 9_133261703_A_G -> chr9:133261703:A:G."""
    return 'chr' + str(vid).replace('_', ':') if vid else ''


# sections of a text-mined article, best first, for the sentences shown with a target x publication link
SECTION_RANK = {'abstract': 0, 'results': 1, 'concl': 2, 'discuss': 3, 'intro': 4, 'other': 5, 'table': 6}


def sentences(ev, k, drop=()):
    """(title, words the target was matched as, n distinct sentences, the k best sentences) of one text-mined
    evidence row (Europe PMC). A sentence matched through a word in `drop` (a word mapped to the wrong gene,
    config/curation.yaml) is left out, and so is the word."""
    title, seen, matched, rows = '', set(), OrderedDict(), []
    for i, s in enumerate(ev.get('textMiningSentences') or []):
        text = ' '.join(re.sub(r'<[^>]+>', ' ', s['text']).split())
        term = ''
        if s.get('tStart') is not None and s.get('tEnd') is not None:
            term = s['text'][s['tStart']:s['tEnd']].strip()
        if s['section'] == 'title':
            title = title or text
        if term in drop:
            continue
        if term:
            matched[term] = None
        if s['section'] == 'title' or text in seen:
            continue
        seen.add(text)
        rows.append((SECTION_RANK.get(s['section'], 7), i, text))
    return title, list(matched), len(rows), [t for _r, _i, t in sorted(rows)[:k]]


# ── American spelling of the text the graph shows ───────────────────────────
# Text carried from the saved sources (condition names, statements, evidence types, design labels) is shown in
# American English; quotes from publications and reference titles are never changed.
IZE = ['summar', 'organ', 'recogn', 'priorit', 'minim', 'standard', 'character', 'optim', 'visual', 'util', 'categor',
       'harmon', 'normal', 'random', 'general', 'hospital', 'immun', 'special', 'final', 'initial', 'local', 'colocal',
       'emphas', 'real', 'stabil', 'maxim', 'critic', 'apolog', 'author', 'central', 'mobil', 'neutral', 'polar',
       'sensit', 'steril', 'synthes', 'hypothes', 'memor', 'rational', 'capital', 'fertil', 'familiar', 'digit',
       'custom', 'vocal', 'industrial', 'commercial', 'symbol', 'stigmat', 'sanit', 'tantal', 'jeopard', 'scrutin']
US_RULES = []
for s in IZE:
    US_RULES.append((re.compile(rf'\b([{s[0].upper()}{s[0]}]{s[1:]})is(e|ed|es|ing|ation|ations|er|ers)\b'), r'\1iz\2'))
US_RULES += [
    (re.compile(r'\b([Aa]nal|[Cc]atal|[Pp]aral|[Dd]ial|[Hh]ydrol)ys(e|ed|ing)\b'), r'\1yz\2'),
    (re.compile(r'\b([Cc]ol|[Ff]av|[Bb]ehavi|[Nn]eighb|[Ll]ab|[Hh]on|[Hh]um|[Ff]lav|[Tt]um|[Rr]um|[Vv]ap|[Hh]arb|'
                r'[Ee]ndeav)our(s|ed|ing|ful|ite|ites|hood|hoods|able)?\b'), r'\1or\2'),
    (re.compile(r'\b([Cc])entre(s)?\b'), r'\1enter\2'), (re.compile(r'\b([Cc])entred\b'), r'\1entered'),
    (re.compile(r'\b([Ff])ibre(s)?\b'), r'\1iber\2'), (re.compile(r'\b([Tt])heatre(s)?\b'), r'\1heater\2'),
    (re.compile(r'\b([Cc])atalogue(s|d)?\b'), r'\1atalog\2'), (re.compile(r'\b([Aa])nalogue(s)?\b'), r'\1nalog\2'),
    (re.compile(r'\b([Pp])rogramme(s)?\b'), r'\1rogram\2'),
    (re.compile(r'\b([Ll]abel|[Mm]odel|[Tt]ravel|[Cc]ancel|[Ll]evel|[Ss]ignal|[Cc]hannel|[Ff]uel|[Tt]otal|[Cc]ounsel)l'
                r'(ed|ing|er|ers)\b'), r'\1\2'),
    (re.compile(r'\b([Gg])rey(s|ish|ed)?\b'), r'\1ray\2'),
    (re.compile(r'([a-z])aemi(a|as|c)\b'), r'\1emi\2'),
    (re.compile(r'\b([Hh])aem'), r'\1em'),
    (re.compile(r'\b([Oo])edema'), lambda m: ('E' if m.group(1) == 'O' else 'e') + 'dema'),
    (re.compile(r'\b([Oo])esophag'), lambda m: ('E' if m.group(1) == 'O' else 'e') + 'sophag'),
    (re.compile(r'\b([Oo])estr'), lambda m: ('E' if m.group(1) == 'O' else 'e') + 'str'),
    (re.compile(r'\b([Dd])iarrhoea'), r'\1iarrhea'), (re.compile(r'\b([Ff])oet(al|us)\b'), r'\1et\2'),
    (re.compile(r'\b([Cc])oeliac'), r'\1eliac'), (re.compile(r'\b([Aa])etiolog'), lambda m: ('E' if m.group(1) == 'A'
                                                                                              else 'e') + 'tiolog'),
    (re.compile(r'\b([Aa])naesthe'), r'\1nesthe'), (re.compile(r'\b([Pp])aediatr'), r'\1ediatr'),
    (re.compile(r'\b([Oo])rthopaed'), r'\1rthoped'), (re.compile(r'\b([Gg])ynaecolog'), r'\1ynecolog'),
    (re.compile(r'\b([Ll])eucocyte'), r'\1eukocyte'), (re.compile(r'\b([Ss])ulph'), r'\1ulf'),
    (re.compile(r'\b([Aa])rtefact'), r'\1rtifact'), (re.compile(r'\b([Jj])udgement'), r'\1udgment'),
    (re.compile(r'\b([Aa])cknowledgement'), r'\1cknowledgment'), (re.compile(r'\b([Aa])geing\b'), r'\1ging'),
    (re.compile(r'\b([Ww])hilst\b'), r'\1hile'), (re.compile(r'\b([Aa])mongst\b'), r'\1mong'),
    (re.compile(r'\b([Tt])owards\b'), r'\1oward'), (re.compile(r'\b([Aa])fterwards\b'), r'\1fterward'),
    (re.compile(r'\b([Ff])ocuss(ed|ing)\b'), r'\1ocus\2'), (re.compile(r'\b([Ll])earnt\b'), r'\1earned'),
    (re.compile(r'\bper cent\b'), 'percent'), (re.compile(r'\b([Ll])icence\b'), r'\1icense'),
    (re.compile(r'\b([Dd])efence\b'), r'\1efense'), (re.compile(r'\b([Pp])ractis(e|ed|ing)\b'), r'\1ractic\2'),
    (re.compile(r'\b([Ee])nrol(ment|ments)?\b'), r'\1nroll\2'), (re.compile(r'\b([Ff])ulfil(ment)?\b'), r'\1ulfill\2'),
    (re.compile(r'\b([Ss])ceptic'), r'\1keptic'), (re.compile(r'\b([Ss])peciality\b'), r'\1pecialty'),
    (re.compile(r'\b([Mm])anoeuvre'), r'\1aneuver'), (re.compile(r'\b([Tt])umour'), r'\1umor'),
]


def us_english(text):
    """British spellings of a text in American English (colocalisation -> colocalization, haem- -> hem-, ...)."""
    if not isinstance(text, str):
        return text
    for rx, rep in US_RULES:
        text = rx.sub(rep, text)
    return text


def log(msg):
    print(f'[{dt.datetime.now():%H:%M:%S}] {msg}', flush=True)


def gql(url, query, variables=None, tries=5):
    """One GraphQL request. Retries server errors, 429 and network failures with exponential backoff;
    any other HTTP error or a GraphQL error raises at once."""
    body = json.dumps({'query': query, 'variables': variables or {}}).encode()
    err = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, data=body, headers={'Content-Type': 'application/json'})
            d = json.load(urllib.request.urlopen(req, timeout=180))
            if d.get('errors'):
                raise RuntimeError(f'GraphQL errors: {d["errors"]}')
            return d['data']
        except urllib.error.HTTPError as e:
            if e.code < 500 and e.code != 429:
                raise RuntimeError(f'HTTP {e.code}: {e.read().decode()[:800]}')
            err = e
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            err = e
        time.sleep(2 ** k * 3)
    raise RuntimeError(f'Open Targets API unreachable after {tries} tries: {err}')


def get_json(url, tries=6):
    """One REST GET (Ensembl). Retries server errors, 429 and network failures with exponential backoff."""
    err = None
    for k in range(tries):
        try:
            return json.load(urllib.request.urlopen(url, timeout=60))
        except urllib.error.HTTPError as e:
            if e.code < 500 and e.code != 429:
                raise RuntimeError(f'HTTP {e.code} for {url}: {e.read().decode()[:400]}')
            err = e
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            err = e
        time.sleep(2 ** k * 3)
    raise RuntimeError(f'{url} unreachable after {tries} tries: {err}')


def release(url):
    m = gql(url, META)['meta']
    dv, av = m['dataVersion'], m['apiVersion']
    data = f'{dv["year"]}.{dv["month"]}' + (f'.{dv["iteration"]}' if dv.get('iteration') else '')
    return {'data_version': data, 'api_version': f'{av["x"]}.{av["y"]}.{av["z"]}'}


def write_json(obj, path):
    path = Path(path)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'wt') as fh:
        json.dump(obj, fh, indent=None if path.suffix == '.gz' else 1, ensure_ascii=False)


def load_json(path):
    path = Path(path)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as fh:
        return json.load(fh)


def read_tsv(path):
    """Every column as text; '' stays '' (never NaN). Callers convert with num()."""
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, na_values=[])


def write_tsv(df, path):
    """P columns in scientific notation; floats with 10 significant digits."""
    out = df.copy()
    for c in out.columns:
        if (c == 'p' or c.endswith('_p')) and pd.api.types.is_numeric_dtype(out[c]):
            out[c] = out[c].map(lambda x: '' if pd.isna(x) else f'{float(x):.6e}')
    out.to_csv(path, sep='\t', index=False, float_format='%.10g')


def num(s):
    return pd.to_numeric(s, errors='coerce')


def csv_list(s):
    """'a,b,c' -> ['a', 'b', 'c'] (the form in which the .nf passes lists)."""
    return [x.strip() for x in str(s).split(',') if x.strip()]


# ── the CTEPH literature: populations and gene themes ───────────────────────
EAS_COUNTRIES = ('Japan', 'China', 'Korea', 'Taiwan', 'Hong Kong', 'Mongolia')
OTHER_GENES, NO_GENE = 'Other genes', 'No gene named'
POSITIVE = ('associated', 'enriched')
CURATED = 'Curated'   # display name of the curated / Mendelian gene evidence ("Mendelian" in the tables)


def study_population(country, ancestry=''):
    """The population a study names: its country, else its ancestry, else 'not stated'."""
    for v in (country, ancestry):
        v = str(v or '').strip()
        if v and v != 'not stated':
            return v
    return 'not stated'


def is_japanese(s):
    """A literature study of a Japanese population: flagged by the extraction, or done in Japan."""
    return s.japanese == 'True' or 'Japan' in str(s.country)


def is_east_asian(s):
    return is_japanese(s) or s.east_asian == 'True' or any(c in str(s.country) for c in EAS_COUNTRIES)


def load_themes(path):
    """The gene themes (config/gene_themes.yaml); a gene may be listed in one theme only."""
    themes = yaml.safe_load(Path(path).read_text())['themes']
    seen = {}
    for t in themes:
        assert {'key', 'name', 'description', 'genes'} <= set(t), t
        for g in t['genes']:
            assert g not in seen, f'{g} is in two themes: {seen[g]}, {t["key"]}'
            seen[g] = t['key']
    return themes


def theme_order(themes):
    return [t['name'] for t in themes] + [OTHER_GENES, NO_GENE]


def assign_themes(findings, themes):
    """The findings with `theme` and `label` (the gene symbol; for a finding without one, its gene text).
    A symbol joins the theme that lists it; a finding without a symbol joins the first theme whose genes or
    words occur in its gene text."""
    by_gene = {g: t['name'] for t in themes for g in t['genes']}

    def hit(word, text, case):
        return re.search(rf'(?<![A-Za-z0-9]){re.escape(word)}(?![A-Za-z0-9])', text, 0 if case else re.I)

    def one(r):
        if r.symbol:
            return by_gene.get(r.symbol, OTHER_GENES)
        for t in themes:
            if any(hit(g, r.gene, True) for g in t['genes']) or any(hit(w, r.gene, False)
                                                                    for w in t.get('words', [])):
                return t['name']
        return NO_GENE
    f = findings.copy()
    f['theme'] = [one(r) for r in f.itertuples()]
    f['label'] = [r.symbol or (r.gene if t != NO_GENE else 'no gene named') for r, t in zip(f.itertuples(), f.theme)]
    return f


def gene_tallies(f):
    """One row per theme x gene label of assign_themes() output: studies by result (findings about CTEPH itself),
    and studies with a finding about a phenotype or outcome within CTEPH."""
    rows = []
    for (th, lab), g in f.groupby(['theme', 'label'], sort=False):
        about = g.association_with if 'association_with' in g else pd.Series(ABOUT_CTEPH, index=g.index)
        gc = g[about == ABOUT_CTEPH]
        n = lambda res: int(gc[gc.result.isin(res)].pmid.nunique())  # noqa: E731
        rows.append({'theme': th, 'label': lab, 'symbol': g.symbol.iloc[0], 'n_studies': int(g.pmid.nunique()),
                     'n_findings': len(g), 'associated': n(POSITIVE), 'not_associated': n(('not associated',)),
                     'identified': n(('identified',)), 'other': n(('other',)),
                     'within_cteph': int(g[about == WITHIN_CTEPH].pmid.nunique())})
    return pd.DataFrame(rows, columns=['theme', 'label', 'symbol', 'n_studies', 'n_findings', 'associated',
                                       'not_associated', 'identified', 'other', 'within_cteph'])


def tally_text(r):
    """'associated in 3 studies, not associated in 1, variants identified in 2, associated with an outcome within
    CTEPH in 1'."""
    within = int(getattr(r, 'within_cteph', 0) or 0)
    ns = (r.associated, r.not_associated, r.identified, r.other, within)
    parts = [(r.associated, 'associated in'), (r.not_associated, 'not associated in'),
             (r.identified, 'variants identified in'), (r.other, 'other result in'),
             (within, 'associated with a phenotype or outcome within CTEPH in')]
    parts = [f'{t} {n}' for n, t in parts if n]
    if not parts:
        return ''
    word = 'study' if [n for n in ns if n][0] == 1 else 'studies'
    return f'{parts[0]} {word}' + ''.join(f', {p}' for p in parts[1:])


# ── choices and their basis (config/choices.yaml) ───────────────────────────
BASIS_TYPES = ('source', 'data check', 'design', 'user')


def ref_short(text):
    """'Hormozdiari, F., ... (2016). ...' -> 'Hormozdiari 2016'; a documentation entry keeps its label."""
    m = re.match(r'([^,(]+),.*?\((\d{4})\)', text)
    return f'{m.group(1).strip()} {m.group(2)}' if m else text.split('.')[0].strip()


def ref_url(text):
    m = re.search(r'PMID (\d+)', text)
    if m:
        return f'https://pubmed.ncbi.nlm.nih.gov/{m.group(1)}/'
    m = re.search(r'(https?://\S+?)[).,]?(\s|$)', text)
    return m.group(1) if m else ''


def choice_checks(tables):
    """The data checks that choices.yaml cites, computed on this snapshot's tables (a directory)."""
    T = Path(tables)
    cond = read_tsv(T / 'conditions.tsv')
    single = cond[cond.n_shared_regions == '1']
    return {'single_locus_abo': single.shared_regions.value_counts(), 'screening': read_tsv(T / 'cteph_screening.tsv'),
            'studies': read_tsv(T / 'cteph_studies.tsv'), 'regions': read_tsv(T / 'cteph_regions.tsv')}


def load_choices(path, run_params, tables=None):
    """The groups of choices, each item with `value_text` (from the run's parameters when the item names one),
    `refs` [(short, url, full)], and the result of its data check (`check_text`, `check_ok`)."""
    cfg = yaml.safe_load(Path(path).read_text())
    refs = cfg.get('references') or {}
    checks = choice_checks(tables) if tables else {}
    out = []
    for g in cfg['groups']:
        items = []
        for it in g['items']:
            it = dict(it)
            assert it.get('basis_type') in BASIS_TYPES, (it['key'], it.get('basis_type'))
            if 'param' in it:
                assert it['param'] in run_params, f'{it["key"]}: no parameter {it["param"]} in the run'
                v = run_params[it['param']]
                v = ('on' if v else 'off') if isinstance(v, bool) else ('the tool default' if v == '' else v)
                it['value_text'] = it.get('format', '{v}').format(v=v)
            else:
                it['value_text'] = str(it.get('value', ''))
            for k in it.get('refs') or []:
                assert k in refs, f'{it["key"]}: reference {k} is not in the references block'
            it['refs'] = [(ref_short(refs[k]), ref_url(refs[k]), refs[k]) for k in it.get('refs') or []]
            it['check_text'], it['check_ok'] = '', True
            if it.get('check') == 'single_locus_abo' and checks:
                vc = checks['single_locus_abo']
                it['check_text'] = (f'This snapshot: {int(vc.sum())} conditions share a single CTEPH region ('
                                    + ', '.join(f'{k} {v}' for k, v in vc.items()) + ').')
                it['check_ok'] = bool(len(vc)) and vc.index[0] == 'ABO'
            if it.get('check') == 'published_regions' and checks:
                rg = checks['regions']
                it['check_text'] = (f'This snapshot: Open Targets holds {int(rg.variant_id.ne("").sum())} of the '
                                    f'{len(rg)} lead variants and has a CTEPH credible set at '
                                    f'{int(rg.cteph_study_locus_id.ne("").sum())} ('
                                    + ', '.join(rg[rg.cteph_study_locus_id.ne('')].gene) + ').')
            if it.get('benchmark') and checks:
                b = [str(p) for p in it['benchmark']]
                sc, st = checks['screening'], checks['studies']
                found = sum(p in set(sc.pmid) for p in b)
                passed = sum(p in set(sc[sc.stage == 'model'].pmid) for p in b)
                incl = sum(p in set(st.pmid) for p in b)
                it['check_text'] = (f'This snapshot: of the {len(b)} original genetic studies cited by the curated '
                                    f'seed, the search finds {found}, the pre-screen keeps {passed}, and {incl} are '
                                    'included.')
                it['check_ok'] = found == passed == incl == len(b)
            items.append(it)
        out.append({'name': g['name'], 'items': items})
    return out, refs


# ── the evidence ladder of the CTEPH genes ───────────────────────────────────
# A gene's level is the strongest evidence it has; the definitions and their bases are in config/choices.yaml
# (key: ladder). Levels 1-2 come from the published GWAS regions (config/cteph_gwas.yaml), 3-7 from the genetic
# studies in the literature.
LADDER = {1: ('G', 'Genome-wide significant in the CTEPH GWAS'),
          2: ('G*', 'Genome-wide significant only after leveraging PE'),
          3: ('', 'Associated in two or more studies'),
          4: ('', 'Associated in one study'),
          5: ('', 'Variant identified, no association test'),
          6: ('', 'Tested, not associated'),
          7: ('', 'Other result only')}
# what a finding of the CTEPH literature is about (curate.py; config/curation.yaml, findings): only findings about
# CTEPH itself count toward levels 3-6
ABOUT_CTEPH = 'CTEPH'
WITHIN_CTEPH = 'a phenotype or outcome within CTEPH'


def gene_ladder(findings, studies, themes, regions):
    """One row per gene with genetic evidence for CTEPH: its ladder level, function (gene theme), GWAS region and
    literature tallies. findings / studies: the curated literature tables (curate.py); regions: the published
    regions (with their `ensembl` map), as fetched. Levels 3-6 count the findings about CTEPH itself; a finding
    about a phenotype or outcome within CTEPH (column within_cteph), a functional experiment or a background
    statement leaves a gene at level 7 when it has nothing else."""
    fx = assign_themes(findings, themes).merge(studies[['pmid']].assign(eas=[is_east_asian(r) for r in
                                                                                studies.itertuples()]),
                                               on='pmid', how='left')
    fx = fx[fx.symbol.ne('')]
    if 'association_with' not in fx:
        fx['association_with'] = ABOUT_CTEPH
    by_gene = {g: t['name'] for t in themes for g in t['genes']}
    reg = {g: r for r in regions for g in r['genes']}
    rows = []
    for g in sorted(set(fx.symbol) | set(reg)):
        f = fx[fx.symbol == g]
        fc = f[f.association_with == ABOUT_CTEPH]
        n = lambda res: int(fc[fc.result.isin(res)].pmid.nunique())  # noqa: E731
        assoc, nota, ident, oth = n(POSITIVE), n(('not associated',)), n(('identified',)), n(('other',))
        within = int(f[f.association_with == WITHIN_CTEPH].pmid.nunique())
        r = reg.get(g)
        if r is not None:
            level = 1 if r['evidence'] == 'meta' else 2
        else:
            level = 3 if assoc >= 2 else 4 if assoc == 1 else 5 if ident else 6 if nota else 7
        ens = (r or {}).get('ensembl', {}).get(g) or (f.ensembl_id.iloc[0] if len(f) and 'ensembl_id' in f else '')
        rows.append({'symbol': g, 'ensembl_id': ens, 'level': level, 'code': LADDER[level][0],
                     'level_label': LADDER[level][1], 'theme': by_gene.get(g, OTHER_GENES),
                     'n_studies': int(f.pmid.nunique()), 'associated': assoc, 'not_associated': nota,
                     'identified': ident, 'other': oth, 'within_cteph': within,
                     'eas_studies': int(f[f.eas.eq(True)].pmid.nunique()),
                     'gwas_region': (r or {}).get('gene', ''), 'gwas_rsid': (r or {}).get('rsid', ''),
                     'gwas_p': (r or {}).get('p'), 'gwas_v': (r or {}).get('v'), 'gwas_tier': (r or {}).get('tier'),
                     'gwas_or': (r or {}).get('or'),
                     'pmids': ' '.join(sorted(set(f.pmid), key=lambda x: int(x) if str(x).isdigit() else 0))})
    order = {t['name']: i for i, t in enumerate(themes)}
    out = pd.DataFrame(rows)
    out['_t'] = out.theme.map(lambda t: order.get(t, len(order)))
    out = out.sort_values(['level', '_t', 'associated', 'n_studies', 'symbol'],
                          ascending=[True, True, False, False, True]).drop(columns='_t')
    return out.reset_index(drop=True)
