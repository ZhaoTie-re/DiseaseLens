#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Helpers shared by the DiseaseLens step scripts: the Open Targets
#           GraphQL client (retries with backoff), gzipped JSON I/O, the TSV
#           conventions of the tables, the disease profile (profile.json, written
#           by profile.py), the gene themes, the evidence ladder, the choices
#           registry and the population labels of the literature.
#           No paths and no disease constants: every script receives its inputs,
#           and the profile of its disease, on the command line from main.nf.
# ---------------------------------------------------------------------------
import datetime as dt
import http.client
import math
import gzip
import hashlib
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
# where the products send a reader, and the note on the reuse of the data they show
REPOSITORY = 'https://github.com/ZhaoTie-re/DiseaseLens'
DATA_NOTE = ("Open Targets data are CC0; quotes and abstracts remain under their publishers' copyright, each linked to "
             'its PubMed record.')


def disease_badge(short, did):
    """The mark and the color of a disease, on its viewer's browser-tab icon and in its header, so that the tabs of
    several diseases tell apart: the first two letters or digits of the short name, on a dark hue derived from the Open Targets ID (stable across
    runs, white text readable on it)."""
    mark = ''.join(ch for ch in short if ch.isalnum())[:2].upper() or '?'
    hue = int(hashlib.sha1(did.encode()).hexdigest(), 16) % 360
    return mark, f'hsl({hue}, 55%, 32%)'


def citation(short, snapshot, release):
    """How to cite a snapshot of a disease's graph."""
    return f'DiseaseLens, {short} knowledge graph, snapshot {snapshot}, Open Targets {release}, {REPOSITORY}'
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


def mining_links(evidence):
    """The text-mined links of the disease (Europe PMC evidence), one per gene and publication: Open Targets can give
    a pair more than once (for the disease and for a descendant, or for two matches); the link keeps its best score
    and every sentence. -> [evidence record]"""
    links = {}
    for x in (x for x in evidence if x['datasourceId'] == 'europepmc'):
        k = (x['target']['id'], tuple(x.get('literature') or []))
        if k not in links:
            links[k] = dict(x, textMiningSentences=list(x.get('textMiningSentences') or []))
            continue
        m = links[k]
        m['score'] = max(float(m['score']), float(x['score']))
        seen = {json.dumps(s, sort_keys=True) for s in m['textMiningSentences']}
        m['textMiningSentences'] += [s for s in x.get('textMiningSentences') or []
                                     if json.dumps(s, sort_keys=True) not in seen]
    return list(links.values())


# a clinical relation that found no association: kept as contrary evidence beside a relation, never a relation itself
NO_ASSOCIATION = 'no association found'


def relates(evidence_type, relation):
    """True when a piece of relation evidence relates the condition to the disease: any evidence but clinical
    evidence that found no association (the direction the readers, or the curated list, give)."""
    return not (evidence_type == 'clinical' and str(relation).strip().lower() == NO_ASSOCIATION)


# the design of the publication behind a clinical relation, best first: from its PubMed publication types and MeSH
# design headings as NLM indexes them (a case report stays one whatever else it is indexed with); the readers' design
# words only for a record NLM has not indexed yet
RELATION_DESIGNS = ['Guideline or consensus', 'Meta-analysis or systematic review', 'Cohort or case-control study',
                    'Cross-sectional or other observational study', 'Case report', 'Other or not stated']
PT_GUIDELINE = {'Practice Guideline', 'Guideline', 'Consensus Development Conference',
                'Consensus Development Conference, NIH'}
PT_META = {'Meta-Analysis', 'Systematic Review'}
MESH_COHORT = {'Cohort Studies', 'Longitudinal Studies', 'Follow-Up Studies', 'Prospective Studies',
               'Retrospective Studies', 'Case-Control Studies', 'Registries'}


def relation_design(record, words=''):
    """The design class (RELATION_DESIGNS) of the record behind a clinical relation."""
    pt, mesh, w = set(record.get('pub_types') or []), set(record.get('mesh') or []), str(words).lower()
    if pt & PT_GUIDELINE:
        return RELATION_DESIGNS[0]
    if pt & PT_META:
        return RELATION_DESIGNS[1]
    if 'Case Reports' in pt:
        return RELATION_DESIGNS[4]
    if mesh & MESH_COHORT:
        return RELATION_DESIGNS[2]
    if 'Cross-Sectional Studies' in mesh or 'Observational Study' in pt:
        return RELATION_DESIGNS[3]
    if not mesh:
        for i, pat in enumerate((r'guideline|consensus', r'meta-analy|systematic review',
                                 r'cohort|case-control|case control|follow-up|prospective|retrospective|registry',
                                 r'cross-sectional|case series|observational', r'case report')):
            if re.search(pat, w):
                return RELATION_DESIGNS[i]
    return RELATION_DESIGNS[5]


# the failures a request retries: the network, a timeout, a dropped connection, a response cut short
# (http.client.IncompleteRead) or a JSON body cut short
NETWORK_ERRORS = (urllib.error.URLError, TimeoutError, ConnectionError, http.client.HTTPException)


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
        except NETWORK_ERRORS + (json.JSONDecodeError,) as e:
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
        except NETWORK_ERRORS + (json.JSONDecodeError,) as e:
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


# column names of the snapshots built before the disease profile (CTEPH, 2026-10-16): read under their new names
LEGACY_COLUMNS = {'cteph_gene': 'disease_gene', 'cteph_locus': 'disease_locus',
                  'cteph_study_locus_id': 'disease_study_locus_id', 'within_cteph': 'within_disease'}


def read_tsv(path, columns=None):
    """Every column as text; '' stays '' (never NaN). Callers convert with num(). An empty file gives an empty
    table (with `columns`, when given); legacy column names are read under their new names."""
    try:
        df = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, na_values=[])
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=columns or [])
    df = df.rename(columns={k: v for k, v in LEGACY_COLUMNS.items() if k in df.columns})
    for c in columns or []:
        if c not in df.columns:
            df[c] = ''
    return df


def legacy_path(directory, name, old):
    """A file of a snapshot under its name, or under the name it had before the disease profile."""
    d = Path(directory)
    return d / name if (d / name).exists() or not (d / old).exists() else d / old


def write_tsv(df, path, columns=None):
    """P columns in scientific notation; floats with 10 significant digits. A table without rows keeps its header
    (`columns` gives it one when it was built from an empty list)."""
    if columns is not None and not len(df.columns):
        df = pd.DataFrame(columns=columns)
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




def a_(word):
    """The word with its indefinite article: 'a GWAS', 'an IPAH', 'an HPAH' (an acronym is read letter by letter),
    'an asthma'."""
    if len(word) > 1 and word.replace('/', '').isupper():
        return ('an ' if word[0] in 'AEFHILMNORSX' else 'a ') + word
    return ('an ' if word[:1].lower() in 'aeiou' else 'a ') + word


# ── the disease profile (profile.json, written by profile.py) ───────────────
class Profile:
    """The resolved profile of one disease: its identity (Open Targets ID, name, short name), the PubMed queries,
    the ontology rule, the categories of related conditions, and the merged configuration (curation rules, choices,
    curated seed, published GWAS regions, gene themes, texts). Disease names enter text as {short}, {name} and
    {Name} (the name with a capital)."""

    def __init__(self, d):
        self.d = d
        self.key, self.id, self.name, self.short = d['key'], d['id'], d['name'], d['short']
        self.Name = d.get('Name') or self.name[:1].upper() + self.name[1:]
        self.stem = d.get('stem') or f'{self.key}_knowledge_graph'
        self.about = self.short                          # findings about the disease itself (association_with)
        self.curation = d.get('curation') or {}
        kinds = (self.curation.get('findings') or {})
        self.within = (kinds.get('within_disease') or {}).get('label', f'a phenotype or outcome within {self.short}')
        self.seed = d.get('seed')
        self.gwas = d.get('gwas')
        self.themes = d.get('themes')
        self.choices = d.get('choices') or {'groups': [], 'references': {}}
        self.ontology = d.get('ontology') or {}
        self.relations = d.get('relations') or {}
        self.checks = d.get('checks') or {}
        self.facts = d.get('facts') or []
        self.caveats = d.get('caveats') or []
        self.categories = d.get('categories') or []

    def fill(self, s):
        """{short}, {name} and {Name} in a text replaced by the disease's names (nothing else is formatted); 'a {short}'
        takes the article of the short name (an IPAH)."""
        if not isinstance(s, str):
            return s
        art = a_(self.short)
        s = s.replace('a {short}', art).replace('A {short}', art[:1].upper() + art[1:])
        return s.replace('{short}', self.short).replace('{name}', self.name).replace('{Name}', self.Name)

    def text(self, section, key=None, default=''):
        """A text of the disease's text/<section>.yaml (or the whole section), with its names filled in."""
        sec = (self.d.get('text') or {}).get(section) or {}
        if key is None:
            return sec
        v = sec.get(key, default)
        return [self.fill(x) for x in v] if isinstance(v, list) else self.fill(v)

    @property
    def ontology_category(self):
        return self.ontology.get('category') or 'Disease group'

    @property
    def category_order(self):
        """The categories of related conditions, in the order they are shown."""
        out = [c['name'] for c in self.categories]
        for c in (self.ontology_category, SHARED_SIGNAL):
            if c not in out:
                out.append(c)
        return out

    @property
    def category_desc(self):
        out = {c['name']: self.fill(c.get('description', '')) for c in self.categories}
        out.setdefault(self.ontology_category, self.fill(self.ontology.get('description', '')))
        out.setdefault(SHARED_SIGNAL, f'Related only by a genetic signal shared with {self.short}.')
        return out

    @property
    def secondary_label(self):
        """The evidence label of a gene genome-wide significant only in a secondary analysis (code G*)."""
        return ((self.gwas or {}).get('evidence_labels') or {}).get('secondary', 'GWAS (secondary analysis)')


SHARED_SIGNAL = 'Shared GWAS signal only'
SHARED_SIGNAL_OLD = 'Shared genetic signal'   # its name in snapshots built before 2026-10-20


def load_profile(path):
    return Profile(load_json(path))


# ── populations and gene themes of the literature ───────────────────────────
EAS_COUNTRIES = ('Japan', 'China', 'Korea', 'Taiwan', 'Hong Kong', 'Mongolia')
# the words that state an East Asian population in a text (whole words, in their own case: not 'Asian' alone, which
# also names South Asia)
EAS_WORDS = EAS_COUNTRIES + ('Japanese', 'Chinese', 'Korean', 'Taiwanese', 'Mongolian', 'East Asian', 'East Asia')


def countries_named(text, C):
    """The countries a text names as whole words (C: the curation rules' countries: names, aliases mapped, phrases
    that name two countries ignored)."""
    s = str(text or '')
    for w in C.get('ignore') or []:
        s = s.replace(w, ' ')
    names = {n: n for n in C.get('names') or []} | (C.get('aliases') or {})
    return {v for k, v in names.items() if re.search(rf'(?<![A-Za-z]){re.escape(k)}(?![A-Za-z])', s)}
# the group of a gene that neither the curated themes nor Reactome place; findings that name no gene
OTHER_GENES, NO_GENE = 'No Reactome pathway', 'No gene named'
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


def check_themes(themes):
    """The gene themes, checked: a gene may be listed in one theme only."""
    seen = {}
    for t in themes:
        assert {'key', 'name', 'description', 'genes'} <= set(t), t
        for g in t['genes']:
            assert g not in seen, f'{g} is in two themes: {seen[g]}, {t["key"]}'
            seen[g] = t['key']
    return themes


def load_themes(path):
    """The gene themes the snapshot used (04.tables/themes.json), or a gene-themes YAML file."""
    p = Path(path)
    if p.is_dir():
        p = p / 'themes.json'
    data = load_json(p) if p.suffix == '.json' else yaml.safe_load(p.read_text())
    return check_themes(data['themes'] if isinstance(data, dict) else data)


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


def auto_themes(genes, pathways):
    """Gene themes for a disease without curated ones: each gene joins its Reactome top-level pathway that holds the
    most of these genes (ties: alphabetical); the top level 'Disease' (disease variants of normal pathways) is not a
    function. Genes without a Reactome pathway are left to OTHER_GENES. pathways: {symbol: [top-level terms]}."""
    from collections import Counter
    tops = {g: sorted({t for t in pathways.get(g) or [] if t and t != 'Disease'}) for g in genes}
    count = Counter(t for ts in tops.values() for t in ts)
    members = {}
    for g in sorted(genes):
        if tops[g]:
            members.setdefault(sorted(tops[g], key=lambda t: (-count[t], t))[0], []).append(g)
    out = []
    for t, gs in sorted(members.items(), key=lambda x: (-len(x[1]), x[0])):
        out.append({'key': re.sub(r'[^a-z0-9]+', '_', t.lower()).strip('_'), 'name': t, 'short': t,
                    'description': f'Genes whose Reactome top-level pathway, among their own, is "{t}".',
                    'basis': 'Reactome top-level pathway of the gene (Open Targets target annotation), assigned '
                             'automatically: the top-level pathway that holds the most of the disease\'s genes.',
                    'refs': ['reactome'], 'genes': gs, 'automatic': True})
    return out


def themes_with_fallback(curated, genes, pathways):
    """The curated gene themes of a disease, then, for each gene they do not list, its Reactome top-level pathway as
    auto_themes() assigns it among those genes (a theme name a curated theme has is followed by "(Reactome)"); a
    gene with neither is left to OTHER_GENES."""
    listed = {g for th in curated for g in th['genes']}
    names = {th['name'] for th in curated}
    auto = auto_themes(sorted(set(genes) - listed), pathways)
    for th in auto:
        if th['name'] in names:
            th['name'] = th['short'] = f'{th["name"]} (Reactome)'
            th['key'] = 'reactome_' + th['key']
        th['description'] = th['description'] + ' (a gene the curated themes do not list)'
    return [dict(th, automatic=False) for th in curated] + auto


def gene_tallies(f, about, within):
    """One row per theme x gene label of assign_themes() output: studies by result (findings about the disease
    itself), and studies with a finding about a phenotype or outcome within the disease."""
    rows = []
    for (th, lab), g in f.groupby(['theme', 'label'], sort=False):
        ab = g.association_with if 'association_with' in g else pd.Series(about, index=g.index)
        gc = g[ab == about]
        n = lambda res: int(gc[gc.result.isin(res)].pmid.nunique())  # noqa: E731
        rows.append({'theme': th, 'label': lab, 'symbol': g.symbol.iloc[0], 'n_studies': int(g.pmid.nunique()),
                     'n_findings': len(g), 'associated': n(POSITIVE), 'not_associated': n(('not associated',)),
                     'identified': n(('identified',)), 'other': n(('other',)),
                     'within_disease': int(g[ab == within].pmid.nunique())})
    return pd.DataFrame(rows, columns=['theme', 'label', 'symbol', 'n_studies', 'n_findings', 'associated',
                                       'not_associated', 'identified', 'other', 'within_disease'])


def tally_text(r, short):
    """'associated in 3 studies, not associated in 1, variants identified in 2, associated with an outcome within
    <disease> in 1'."""
    within = int(getattr(r, 'within_disease', 0) or 0)
    ns = (r.associated, r.not_associated, r.identified, r.other, within)
    parts = [(r.associated, 'associated in'), (r.not_associated, 'not associated in'),
             (r.identified, 'variants identified in'), (r.other, 'other result in'),
             (within, f'associated with a phenotype or outcome within {short} in')]
    parts = [f'{t} {n}' for n, t in parts if n]
    if not parts:
        return ''
    word = 'study' if [n for n in ns if n][0] == 1 else 'studies'
    return f'{parts[0]} {word}' + ''.join(f', {p}' for p in parts[1:])


# ── choices and their basis (config/choices.yaml merged with the disease's) ──
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
    """The tables the data checks of the choices read, from this snapshot's 04.tables."""
    T = Path(tables)
    cond = read_tsv(T / 'conditions.tsv', ['n_shared_regions', 'shared_regions'])
    single = cond[cond.n_shared_regions == '1']
    tidy = json.loads((T / 'tidy.json').read_text()) if (T / 'tidy.json').exists() else {}
    return {'single_region': single.shared_regions.value_counts(),
            'screening': read_tsv(T / 'disease_screening.tsv', ['pmid', 'stage']),
            'studies': read_tsv(T / 'disease_studies.tsv', ['pmid']),
            'regions': read_tsv(T / 'disease_regions.tsv', ['gene', 'variant_id', 'disease_study_locus_id', 'source']),
            'not_regions': not_regions_of(T),
            'genes': read_tsv(T / 'disease_genes.tsv', ['symbol', 'code', 'n_studies']),
            'tidy': tidy}


def load_choices(P, run_params, tables=None):
    """The groups of choices of a disease (the profile's merged registry), each item with `value_text` (from the
    run's parameters when the item names one), `refs` [(short, url, full)], and the result of its data check
    (`check_text`, `check_ok`)."""
    cfg = P.choices
    refs = cfg.get('references') or {}
    checks = choice_checks(tables) if tables else {}
    # the disease's own values: its identifier and its query (from the profile), the snapshot date (from the tables)
    run_params = {'DiseaseCore': P.id, 'PubmedQuery': P.d.get('pubmed_query', ''), 'Rebuild': False,
                  **({'FetchDate': checks['tidy'].get('fetch_date')} if checks.get('tidy') else {}), **run_params}
    short = P.short
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
                it['value_text'] = it.get('format', '{v}').replace('{v}', str(v))
            else:
                it['value_text'] = str(it.get('value', ''))
            for k in it.get('refs') or []:
                assert k in refs, f'{it["key"]}: reference {k} is not in the references block'
            it['refs'] = [(ref_short(refs[k]), ref_url(refs[k]), refs[k]) for k in it.get('refs') or []]
            it['check_text'], it['check_ok'] = '', True
            if it.get('check') == 'single_region' and checks:
                vc = checks['single_region']
                it['check_text'] = (f'This snapshot: {int(vc.sum())} conditions share a single {short} region'
                                    + (' (' + ', '.join(f'{k} {v}' for k, v in vc.items()) + ').' if len(vc) else '.'))
                want = P.checks.get('single_region_top')
                it['check_ok'] = (bool(len(vc)) and vc.index[0] == want) if want else True
            if it.get('check') == 'regions' and checks:
                rg = checks['regions']
                pub = rg[rg.source.isin(['', 'published'])] if 'source' in rg else rg
                auto = rg[rg.source.eq('Open Targets')] if 'source' in rg else rg.iloc[0:0]
                parts = []
                if len(pub):
                    cs = pub[pub.disease_study_locus_id.ne('')]
                    parts.append(f'Open Targets holds {int(pub.variant_id.ne("").sum())} of the {len(pub)} published '
                                 f'lead variants and has {a_(short)} credible set at {len(cs)}'
                                 + (f' ({", ".join(cs.gene)})' if len(cs) else '') + '.')
                if len(auto):
                    parts.append(f'{len(auto)} region{"s" if len(auto) != 1 else ""} from the genome-wide '
                                 f'significant credible sets of the {short} GWAS in Open Targets.')
                nr = checks.get('not_regions')
                if nr is not None and len(nr):
                    parts.append(not_regions_sentence(nr, short))
                it['check_text'] = 'This snapshot: ' + (' '.join(parts) if parts else
                                                        f'no genome-wide significant GWAS region of {short}.')
            if it.get('check') == 'curated_genes' and checks:
                gn = checks['genes']
                # every gene with curated evidence, whatever its code (a G gene can have curated evidence too)
                cur = gn[gn.code.str.contains('C', regex=False) | gn.curated_sources.ne('')] if len(gn) else gn
                found = cur[pd.to_numeric(cur.n_studies, errors='coerce').fillna(0) > 0]
                it['check_text'] = (f'This snapshot: {len(found)} of the {len(cur)} genes with curated evidence for '
                                    f'{short} are reported by a genetic study kept by the extraction'
                                    + (f' (not found: {", ".join(sorted(set(cur.symbol) - set(found.symbol)))})'
                                       if len(found) < len(cur) else '') + '.' if len(cur) else
                                    f'This snapshot: no gene has curated evidence for {short} in Open Targets.')
                # judged only on a literature read in full: not in a run that reads the first records only
                # (LlmLimit), or that left records unread (offline, nothing cached)
                ex = checks['tidy'].get('extract') or {}
                partial = int(ex.get('not_read') or 0) or int(run_params.get('LlmLimit') or 0)
                if partial and len(found) < len(cur):
                    it['check_text'] += (' Not judged: the literature of this run was not read in full ('
                                         + (f'{ex["not_read"]} records not read' if ex.get('not_read') else
                                            f'the first {run_params["LlmLimit"]} records only') + ').')
                it['check_ok'] = len(found) == len(cur) or bool(partial)
            if it.get('check') == 'readers' and checks:
                ag = ((checks['tidy'].get('extract') or {}).get('agreement') or {})
                if ag:
                    it['check_text'] = 'This snapshot: ' + '; '.join(
                        f'{k}: {v["both"]} agreed, {v["adjudicated"]} adjudicated' for k, v in ag.items()) + '.'
            if it.get('benchmark') and checks:
                b = [str(p) for p in it['benchmark']]
                sc, st = checks['screening'], checks['studies']
                found = sum(p in set(sc.pmid) for p in b)
                passed = sum(p in set(sc[sc.stage == 'model'].pmid) for p in b)
                incl = sum(p in set(st.pmid) for p in b)
                it['check_text'] = (f'This snapshot: of the {len(b)} original genetic studies of the benchmark, the '
                                    f'search finds {found}, the pre-screen keeps {passed}, and {incl} are included.')
                it['check_ok'] = found == passed == incl == len(b)
            items.append(it)
        out.append({'name': g['name'], 'items': items})
    return out, refs


# ── the evidence ladder of the disease's genes ──────────────────────────────
# A gene's level is the strongest evidence it has; the definitions and their bases are in the choices registry
# (key: ladder). Level 1: a genome-wide significant GWAS region of the disease (G) or curated gene-disease evidence
# for it (C); level 2: a region significant only in a secondary analysis of a published GWAS (G*); levels 3-7: the
# genetic studies in the literature. tidy_graph.py writes the ladder a snapshot uses to 04.tables/ladder.json.
# the name of each level, the same for every disease (rows of the grids and figures); its definition below
LADDER_NAMES = {1: ('G / C', 'Genome-wide or curated evidence'),
                2: ('G*', 'Secondary GWAS analysis'),
                3: ('', 'Replicated association'),
                4: ('', 'Single-study association'),
                5: ('', 'Variant reported, not tested'),
                6: ('', 'Tested, not associated'),
                7: ('', 'No result about {short} itself')}      # each gene with its reason (level_note)
LITERATURE_LEVELS = {k: v for k, v in LADDER_NAMES.items() if k >= 3}
# the meaning shown in the row of a level whose name alone could be misread
LADDER_GLOSSES = {7: 'named in a genetic study of {short}, without a result about {short} itself: not evidence for '
                     'it'}
# a gene of the graph without a level: in no genetic study of the disease
NO_LEVEL = ('No level: in no genetic study of {short} (a text-mined {short} target, or a gene with genetic evidence in '
            'a related condition)')


def ladder_definitions(P):
    """{level: definition}: what each level of the ladder means, the same for every disease, followed by the
    disease's own wording where its profile gives one ('CTEPH: ...')."""
    s = P.short
    lab = (P.gwas or {}).get('labels') or {}
    own = lambda code: f' {s}: {P.fill(lab[code])}.' if lab.get(code) else ''  # noqa: E731
    return {
        1: f'A gene of a genome-wide significant GWAS region of {s} (G: a region of a published GWAS, or a credible '
           f'genome-wide significant signal of its GWAS in Open Targets), or with curated gene–disease evidence for '
           f'{s} (C: ClinGen, Gene2Phenotype, PanelApp, Orphanet, ClinVar with a reviewed pathogenic record, UniProt, '
           f'gene burden).' + own('G'),
        2: 'A gene of a region genome-wide significant only in a secondary analysis of a published GWAS (G*), such as '
           'an analysis that leverages a correlated trait.' + own('G*'),
        3: f'Associated with {s} in two or more genetic studies: replication makes an association credible.',
        4: f'Associated with {s} in one genetic study.',
        5: f'A variant reported in people with {s} (case reports, sequencing series), without an association test.',
        6: f'Tested in genetic studies of {s} and not associated in any.',
        7: f'Named in a genetic study of {s} without a result about {s} itself: a result about an outcome within it or '
           'about another trait, or no association, variant or negative test; the reason is given for each gene.'}


def level_status(lv, P, L):
    """(brief, full): what a level of the ladder holds in this snapshot, when that needs saying, so that an empty or
    partly empty level is never read as missing data. Every disease shows the same seven levels. L: the ladder table
    of the snapshot (level, code, curated_sources). ('', '') when the level speaks for itself."""
    s = P.short
    sub = L[L.level.astype(int) == lv] if len(L) else L
    n = len(sub)
    if lv == 1:
        g = int(sub.code.eq('G').sum()) if n else 0
        c = int((sub.code.eq('C') | sub.curated_sources.ne('')).sum()) if n else 0
        if not n:
            return ('Empty: no GWAS region, no curated evidence',
                    f'Empty: {s} has no credible genome-wide significant GWAS region and no curated gene–disease '
                    'evidence in this snapshot.')
        if not g:
            return ('C only: no GWAS region',
                    f'C only: {s} has no credible genome-wide significant GWAS region in this snapshot (see the GWAS '
                    'regions); the genes at this level have curated gene–disease evidence.')
        if not c:
            return ('G only: no curated evidence',
                    f'G only: Open Targets holds no curated gene–disease evidence for {s} that the rules keep; the '
                    'genes at this level come from its GWAS regions.')
        return '', ''
    if n:
        return '', ''
    if lv == 2:
        if 'G*' in set(((P.gwas or {}).get('codes') or {}).values()):
            return ('Empty in this snapshot',
                    f'Empty: no region of the published GWAS of {s} is significant only in its secondary analysis.')
        return ('Empty: no published secondary GWAS analysis',
                f'Empty: this level holds the genes of a region significant only in a secondary analysis of a '
                f'published GWAS, and none is known for {s} (its profile names none). The level is kept so that the '
                'ladder of every disease has the same seven levels.')
    return 'Empty in this snapshot', 'Empty: no gene of this snapshot meets this level.'


def ladder_spec(P, with_curated=True):
    """The ladder of a disease: {level: (code, name)}, the same seven levels for every disease (an empty level is
    said in its row, level_status); {level: definition}; {code: what the code of a gene's evidence means}."""
    lab = (P.gwas or {}).get('labels') or {}
    g = P.fill(lab.get('G', 'Genome-wide significant in a GWAS of {short}'))
    gs = P.fill(lab.get('G*', 'Genome-wide significant only in a secondary analysis of a published GWAS'))
    c = P.fill('Curated gene–disease evidence for {short}')
    levels = {k: (code, P.fill(name)) for k, (code, name) in LADDER_NAMES.items()}
    codes = {'G': g, 'G*': gs, 'C': c, 'L': P.fill('A positive finding in a genetic study of {short}')}
    evidence = {'GWAS': 'G', P.secondary_label: 'G*', 'Mendelian': 'C', 'Literature': 'L'}
    return {'levels': {str(k): list(v) for k, v in levels.items()},
            'definitions': {str(k): v for k, v in ladder_definitions(P).items()},
            'glosses': {str(k): P.fill(v) for k, v in LADDER_GLOSSES.items()},
            'codes': codes, 'evidence_codes': evidence}


def load_ladder(tables):
    """(LADDER {level: (code, name)}, codes {code: label}, evidence_codes {evidence label: code}) of a snapshot."""
    d = load_json(Path(tables) / 'ladder.json')
    return ({int(k): tuple(v) for k, v in d['levels'].items()}, d['codes'], d['evidence_codes'])


def load_status(tables):
    """{level: {'brief', 'full'}}: what each level holds in a snapshot when that needs saying (level_status)."""
    d = load_json(Path(tables) / 'ladder.json')
    return {int(k): v for k, v in (d.get('status') or {}).items()}


def load_glosses(tables):
    """{level: the meaning shown in its row} of a snapshot's ladder (LADDER_GLOSSES)."""
    d = load_json(Path(tables) / 'ladder.json')
    return {int(k): v for k, v in (d.get('glosses') or {}).items()}


def load_definitions(tables):
    """{level: definition} of a snapshot's ladder (its names, for a ladder written before the definitions)."""
    d = load_json(Path(tables) / 'ladder.json')
    return {int(k): v for k, v in (d.get('definitions') or {k: n for k, (_c, n) in d['levels'].items()}).items()}


LADDER_COLUMNS = ['symbol', 'ensembl_id', 'level', 'code', 'level_label', 'theme', 'n_studies', 'associated',
                  'not_associated', 'identified', 'other', 'within_disease', 'eas_studies', 'gwas_region',
                  'gwas_rsid', 'gwas_p', 'gwas_v', 'gwas_tier', 'gwas_or', 'curated_sources', 'pmids', 'level_note']


def gene_ladder(findings, studies, themes, regions, about, within, curated=None, spec=None):
    """One row per gene with genetic evidence for the disease: its ladder level, function (gene theme), GWAS region,
    curated sources and literature tallies. findings / studies: the curated literature tables (curate.py); regions:
    the GWAS regions of the disease, each with its `code` (G or G*), `genes` and `ensembl` map; curated: {symbol:
    {'ensembl_id', 'sources'}} of the curated gene-disease evidence. Levels 3-6 count the findings about the disease
    itself (association_with == about); a finding about a phenotype or outcome within it, a functional experiment
    or a background statement leaves a gene at level 7 when it has nothing else."""
    curated = curated or {}
    codes = (spec or {}).get('codes') or {}
    fx = assign_themes(findings, themes).merge(
        studies[['pmid']].assign(eas=[is_east_asian(r) for r in studies.itertuples()]), on='pmid', how='left')
    fx = fx[fx.symbol.ne('')]
    if 'association_with' not in fx:
        fx['association_with'] = about
    by_gene = {g: t['name'] for t in themes for g in t['genes']}
    reg = {}
    for r in regions:
        for g in r['genes']:
            if g not in reg or (reg[g]['code'] != 'G' and r['code'] == 'G'):
                reg[g] = r
    rows = []
    for g in sorted(set(fx.symbol) | set(reg) | set(curated)):
        f = fx[fx.symbol == g]
        fc = f[f.association_with == about]
        n = lambda res: int(fc[fc.result.isin(res)].pmid.nunique())  # noqa: E731
        assoc, nota, ident, oth = n(POSITIVE), n(('not associated',)), n(('identified',)), n(('other',))
        within_n = int(f[f.association_with == within].pmid.nunique())
        r, c = reg.get(g), curated.get(g)
        if r is not None and r['code'] == 'G':
            level, code = 1, 'G'
        elif c:
            level, code = 1, 'C'
        elif r is not None:
            level, code = 2, 'G*'
        else:
            level, code = (3 if assoc >= 2 else 4 if assoc == 1 else 5 if ident else 6 if nota else 7), ''
        label = codes.get(code) or ((spec or {}).get('levels') or {}).get(str(level), ('', ''))[1] or \
            LITERATURE_LEVELS.get(level, ('', ''))[1].replace('{short}', about)
        note = ''
        if level == 7:                     # why: the gene's findings are about something else, or name it only
            other_kinds = f[f.association_with != about].groupby('association_with').pmid.nunique()
            parts = [f'{n} {"study" if n == 1 else "studies"} with a result about {k}' for k, n in other_kinds.items()]
            if oth:
                parts.append(f'{oth} {"study" if oth == 1 else "studies"} naming it without an association, a variant '
                             f'or a negative test of {about}')
            note = '; '.join(parts)
        ens = ((r or {}).get('ensembl', {}).get(g) or (c or {}).get('ensembl_id')
               or (f.ensembl_id.iloc[0] if len(f) and 'ensembl_id' in f else ''))
        rows.append({'symbol': g, 'ensembl_id': ens, 'level': level, 'code': code, 'level_label': label,
                     'theme': by_gene.get(g, OTHER_GENES),
                     'n_studies': int(f.pmid.nunique()), 'associated': assoc, 'not_associated': nota,
                     'identified': ident, 'other': oth, 'within_disease': within_n,
                     'eas_studies': int(f[f.eas.eq(True)].pmid.nunique()),
                     'gwas_region': (r or {}).get('gene', ''), 'gwas_rsid': (r or {}).get('rsid', ''),
                     'gwas_p': (r or {}).get('p'), 'gwas_v': (r or {}).get('v'), 'gwas_tier': (r or {}).get('tier'),
                     'gwas_or': (r or {}).get('or'), 'curated_sources': '; '.join((c or {}).get('sources') or []),
                     'pmids': ' '.join(sorted(set(f.pmid), key=lambda x: int(x) if str(x).isdigit() else 0)),
                     'level_note': note})
    out = pd.DataFrame(rows, columns=LADDER_COLUMNS)
    if not len(out):
        return out
    order = {t['name']: i for i, t in enumerate(themes)}
    out['_t'] = out.theme.map(lambda t: order.get(t, len(order)))
    out = out.sort_values(['level', '_t', 'associated', 'n_studies', 'symbol'],
                          ascending=[True, True, False, False, True]).drop(columns='_t')
    return out.reset_index(drop=True)


# ── the GWAS regions of the disease ──────────────────────────────────────────
# curated gene evidence: a ClinVar link counts only with a pathogenic or likely pathogenic record whose review status
# has two stars or more (ClinVar's review status); a gene needs an approved symbol
CLINVAR_REVIEWED = ('criteria provided, multiple submitters, no conflicts', 'reviewed by expert panel',
                    'practice guideline')
CLINVAR_PATHOGENIC = ('pathogenic', 'likely pathogenic', 'pathogenic/likely pathogenic')


def curated_links(M, mrec):
    """The curated gene-condition links that count, and their records: a ClinVar link needs a pathogenic or likely
    pathogenic record with a review status of two stars or more, since a single submitter or a record without
    assertion criteria is not enough to make a gene a curated gene of a disease; a gene without an approved symbol
    (its Ensembl ID in place of a symbol) is left out. -> (links, records, {reason: links dropped})"""
    if not len(M):
        return M, mrec, {'clinvar_review': 0, 'no_symbol': 0}
    ok = mrec[mrec.source.eq('ClinVar') & mrec.review_status.isin(CLINVAR_REVIEWED) &
              mrec.clinical_significance.map(lambda v: bool(set(str(v).split('; ')) & set(CLINVAR_PATHOGENIC)))]
    good = set(zip(ok.node_id, ok.ensembl_id))
    weak = M.source.eq('ClinVar') & pd.Series([k not in good for k in zip(M.node_id, M.ensembl_id)], index=M.index)
    nosym = M.symbol.fillna('').eq('') | M.symbol.eq(M.ensembl_id)
    drop = weak | nosym
    gone = set(zip(M[drop].node_id, M[drop].source, M[drop].ensembl_id))
    keep = [k not in gone for k in zip(mrec.node_id, mrec.source, mrec.ensembl_id)]
    return (M[~drop].copy(), mrec[keep].copy(),
            {'clinvar_review': int((weak & ~nosym).sum()), 'no_symbol': int(nosym.sum())})


def genome_wide(mant, exp, p_max=5e-8):
    """P < p_max, from Open Targets' mantissa and exponent (no float underflow)."""
    if mant is None or exp is None or float(mant) <= 0:
        return False
    return math.log10(float(mant)) + int(exp) < math.log10(p_max)


def shared_relations(co, ov, alias, short):
    """The shared-signal relations of the conditions, from the colocalizations (co) and lead-variant overlaps (ov) of
    the disease's regions: per condition (its ID through alias), one row per kind with the regions shared, the number
    of its studies and the best value (CLPP or H4; PIP), graded by the number of the disease's regions it shares
    across both kinds. -> (rows, {node: regions shared}, {node: name})"""
    def val(x):
        try:
            v = float(x)
        except (TypeError, ValueError):
            return 0.0
        return 0.0 if math.isnan(v) else v

    def collect(df, score):
        agg = {}
        for r in df[df.node_id.ne('')].itertuples():
            x = agg.setdefault(alias.get(r.node_id, r.node_id), {'regions': set(), 'studies': set(), 'best': 0.0,
                                                                  'name': ''})
            x['regions'].add(r.disease_gene)
            x['studies'].add(r.study_id)
            x['best'] = max(x['best'], score(r))
            if not x['name'] or r.node_id == alias.get(r.node_id, r.node_id):     # a merged trait: its own name
                x['name'] = r.condition
        return agg
    cagg = collect(co, lambda r: max(val(r.clpp), val(r.h4)))
    oagg = collect(ov, lambda r: val(r.pip))
    shared = {}
    for agg in (cagg, oagg):
        for i, v in agg.items():
            shared.setdefault(i, set()).update(v['regions'])
    rows, names = [], {}
    for kind, agg, what in (('colocalisation', cagg, f'GWAS signals colocalise with the {short} region'),
                            ('lead-variant overlap', oagg, f'GWAS credible sets contain the {short} lead variant of')):
        for i, v in agg.items():
            names.setdefault(i, v['name'])
            grade = f'two or more {short} regions' if len(shared[i]) >= 2 else f'one {short} region'
            reg = sorted(v['regions'])
            rows.append({'node_id': i, 'evidence_type': kind, 'relation': f'shared GWAS signal ({grade})',
                         'detail': f'{len(v["studies"])} {what}{"s" if kind == "colocalisation" and len(reg) > 1 else ""} '
                                   f'{", ".join(reg)}',
                         'quote': '', 'design': '', 'population': '', 'east_asian': False, 'pmid': '',
                         'citation': 'Open Targets colocalisation' if kind == 'colocalisation' else
                                     'Open Targets credible sets', 'strength': f'{v["best"]:.3f}',
                         'loci': ' '.join(reg), 'n_loci': len(reg), 'n_studies': len(v['studies']), 'agreement': ''})
    return rows, shared, names


def not_in_gene_map(cond, shown, regions_at=''):
    """The related conditions that are no column of the gene map, by reason, each counted:
    ([(key, n, text)], 'Not shown: ...'). cond: conditions.tsv; shown: the node IDs of the columns; regions_at: where
    the shared regions are shown ('sheet 5'), '' when the disease has no GWAS region."""
    g, note = cond.genetics, cond.genetics_note.astype(str)
    many = g.eq('coloc_only') & note.str.contains('more than')
    caps = sorted({m.group(1) for m in (re.search(r'more than (\d+)', x) for x in note[many]) if m})
    where = f' (their shared regions: {regions_at})' if regions_at else ''

    def noun(n, one, more):
        return one if n == 1 else more
    counts = [('shared_signal_only', int((g.eq('coloc_only') & ~many).sum())), ('too_many_gwas', int(many.sum())),
              ('no_entry', int(g.eq('none').sum())),
              ('no_gene', int((g.eq('full') & ~cond.node_id.isin(set(shown))).sum()))]
    texts = {'shared_signal_only': lambda n: noun(n, 'trait', 'traits') + ' related only by a shared GWAS signal, '
                                             'whose genes are not fetched by design' + where,
             'too_many_gwas': lambda n: noun(n, 'condition', 'conditions') + (
                 f' with more than {caps[0]} GWAS studies' if caps else ' with too many GWAS studies') +
                 ', whose genes are not fetched (the other evidence is kept)',
             'no_entry': lambda n: noun(n, 'condition', 'conditions') + ' without an Open Targets entry (no genetic '
                                   'data)',
             'no_gene': lambda n: noun(n, 'condition whose genetic data give', 'conditions whose genetic data give')
                                  + ' GWAS or curated evidence for no gene of this map'}
    parts = [(k, n, texts[k](n)) for k, n in counts if n]
    return parts, ('Not shown: ' + '; '.join(f'{n} {x}' for _k, n, x in parts) + '.') if parts else ''


def lead_maf(cs, study):
    """The minor allele frequency of a credible set's lead variant (gnomAD, as Open Targets gives it) in the largest
    ancestry group of its study (the LD reference sample), else the highest across groups; 0 for a lead absent from
    gnomAD; None for a snapshot fetched before the frequencies were (not judged)."""
    if 'lead_af' not in cs:
        return None
    af = {k.replace('_adj', ''): min(v, 1 - v) for k, v in (cs.get('lead_af') or {}).items() if v is not None}
    if not af:
        return 0.0
    pops = sorted(study.get('ldPopulationStructure') or [], key=lambda p: -p['relativeSampleSize'])
    main = pops[0]['ldPopulation'] if pops else ''
    return af[main] if main in af else max(af.values())


NOT_REGION_COLS = ['variant_id', 'rsid', 'study_id', 'p_mantissa', 'p_exponent', 'maf', 'expected_mac',
                   'n_smaller_group', 'reason']


def not_regions_of(tables):
    """The genome-wide significant credible sets of the disease that are no region (tidy_graph.py); an empty table for
    tables written before the rule."""
    f = Path(tables) / 'disease_signals_not_regions.tsv'
    return read_tsv(f, NOT_REGION_COLS) if f.exists() else pd.DataFrame(columns=NOT_REGION_COLS)


def not_regions_sentence(nr, short):
    """The sentence on the genome-wide significant credible sets of the disease's GWAS that are no region (tidy_graph:
    disease_signals_not_regions.tsv); '' when there are none."""
    if not len(nr):
        return ''
    n, studies = len(nr), sorted(set(nr.study_id))
    return (f'{n} genome-wide significant credible set{"s" if n != 1 else ""} of the {short} GWAS in Open Targets '
            f'({", ".join(studies)}) {"are no regions" if n != 1 else "is no region"}: the lead variant is not '
            'credible (too few expected copies of its minor allele in the smaller group of its study, or a rare lead '
            f'without the P of a rare variant), and no other GWAS of {short} replicates it.')


def lead_credibility(c, s, min_maf, min_mac, rare_p):
    """('' or why not, details) for the lead variant of a genome-wide significant credible set c of study s: credible
    when its study carries enough copies of the minor allele, an expected count of at least min_mac in the smaller
    group (cases, or all samples without cases and controls; 2 x n x MAF, the MAF of the study's largest ancestry
    group; Neale lab 2017: about 25 copies are needed to reach genome-wide significance, fewer give miscalibrated
    tests), and, when the lead is rarer than min_maf, when its P is below rare_p (Pulit 2017: 5e-9 once rare variants
    are tested). A snapshot fetched before the frequencies were is not judged."""
    maf = lead_maf(c, s)
    if maf is None:
        return '', {}
    nc, nn = s.get('nCases') or 0, s.get('nControls') or 0
    n = min(nc, nn) if nc and nn else (s.get('nSamples') or 0)
    mac = 2 * n * maf
    logp = math.log10(float(c['pValueMantissa'])) + int(c['pValueExponent'])
    info = {'maf': round(maf, 5), 'expected_mac': round(mac, 1), 'n_smaller_group': n}
    if mac < min_mac:
        return (f'{mac:.0f} expected {"copy" if round(mac) == 1 else "copies"} of the minor allele (MAF {maf:.2g}) in '
                f'the {n:,} {"cases" if nc and nn and nc <= nn else "samples"}, fewer than {min_mac:g}'), info
    if maf < min_maf and logp >= math.log10(rare_p):
        return f'a rare lead (MAF {maf:.2g}) with P not below {rare_p:g}', info
    return '', info


def disease_regions(gwas, published, P, l2g_min, window_kb=500, p_max=5e-8, min_maf=0.01, dropped=None, min_mac=25,
                    rare_p=5e-9):
    """The GWAS regions of the disease: the published regions of the profile (gwas_regions.yaml, as looked up by
    fetch.py), each with the disease's credible sets at its lead variant; and the genome-wide significant credible
    sets (lead P < p_max) of the disease's own GWAS in Open Targets that no published region holds, whose lead variant
    is credible (lead_credibility: enough expected copies of the minor allele in the study, and P < rare_p for a lead
    rarer than min_maf), or, if not, is replicated: another GWAS study of the disease has a genome-wide significant
    credible set within window_kb. A lead neither credible nor replicated is listed in `dropped`, with the reason.
    Credible sets whose leads lie within window_kb are merged. A credible set within
    window_kb of a published lead belongs to that published region. A region from Open Targets is named by its top
    locus-to-gene gene when the score reaches l2g_min (that gene gets the code G), else by the position of its lead
    variant; a name that an earlier region has is followed by the position. Every region: gene, genes, ensembl
    {gene: id}, rsid, variant_id,
    chrom, code (G or G*), source, study_locus_ids, disease_study_locus_id, l2g_gene, l2g_score, and the published
    statistics where there are some."""
    codes = (P.gwas or {}).get('codes') or {}
    regions, taken = [], set()
    sets = [(s, c) for s in gwas.get('studies') or [] for c in s.get('credible_sets') or []]

    def position(variant_id):                       # '9_133257521_T_C' -> ('9', 133257521)
        parts = str(variant_id or '').split('_')
        return (parts[0], int(parts[1])) if len(parts) >= 2 and parts[1].isdigit() else None
    for r in (published or {}).get('regions') or []:
        r = dict(r)
        own = r.pop('cteph_study_locus_id', None) or r.get('disease_study_locus_id', '')
        ids = [c['studyLocusId'] for _s, c in sets if r['rsid'] in ((c.get('variant') or {}).get('rsIds') or [])
               or (r.get('variant_id') and (c.get('variant') or {}).get('id') == r['variant_id'])]
        # the genome-wide significant credible sets near the published lead are the same region
        at = position(r.get('variant_id'))
        near = [c['studyLocusId'] for _s, c in sets if at and c['studyLocusId'] not in ids
                and genome_wide(c.get('pValueMantissa'), c.get('pValueExponent'), p_max)
                and (q := position((c.get('variant') or {}).get('id'))) and q[0] == at[0]
                and abs(q[1] - at[1]) <= window_kb * 1000]
        ids += near
        taken |= set(ids)
        # the region's credible set of the disease: the published study's, else the first at the lead or near it
        own = own or (ids[0] if ids else '')
        regions.append({**r, 'code': codes.get(r.get('evidence'), 'G'), 'source': 'published',
                        'disease_study_locus_id': own, 'study_locus_ids': ids or ([own] if own else []),
                        'l2g_gene': '', 'l2g_score': None})
    hits = []
    for s, c in sets:
        if c['studyLocusId'] in taken or not genome_wide(c.get('pValueMantissa'), c.get('pValueExponent'), p_max):
            continue
        v = c.get('variant') or {}
        parts = (v.get('id') or '').split('_')
        if len(parts) < 2 or not parts[1].isdigit():
            continue
        why, info = lead_credibility(c, s, min_maf, min_mac, rare_p)
        l2g = ((c.get('l2GPredictions') or {}).get('rows') or [{}])[0]
        hits.append({'chrom': parts[0], 'pos': int(parts[1]), 'cs': c, 'study': s, 'l2g': l2g, 'why': why,
                     'info': info, 'logp': math.log10(float(c['pValueMantissa'])) + int(c['pValueExponent'])})
    # a lead that is not credible on its own counts when another GWAS study of the disease replicates it: a genome-wide
    # significant credible set within window_kb
    kept = []
    for h in hits:
        rep_ = [o['study']['id'] for o in hits if o['study']['id'] != h['study']['id'] and o['chrom'] == h['chrom']
                and abs(o['pos'] - h['pos']) <= window_kb * 1000]
        if not h['why'] or rep_:
            kept.append(h)
        elif dropped is not None:
            dropped.append({'variant_id': (h['cs'].get('variant') or {}).get('id', ''),
                            'rsid': ((h['cs'].get('variant') or {}).get('rsIds') or [''])[0],
                            'study_id': h['study']['id'], 'p_mantissa': h['cs']['pValueMantissa'],
                            'p_exponent': h['cs']['pValueExponent'], **h['info'],
                            'reason': h['why'] + '; not replicated by another GWAS of the disease'})
    hits = kept
    hits.sort(key=lambda h: (h['chrom'], h['pos']))
    groups = []
    for h in hits:
        if groups and groups[-1][-1]['chrom'] == h['chrom'] and h['pos'] - groups[-1][-1]['pos'] <= window_kb * 1000:
            groups[-1].append(h)
        else:
            groups.append([h])
    names = {r['gene'] for r in regions}
    for g in groups:
        best = min(g, key=lambda h: h['logp'])
        c, l2g = best['cs'], best['l2g']
        sym = (l2g.get('target') or {}).get('approvedSymbol', '')
        score = l2g.get('score')
        named = bool(sym) and score is not None and score >= l2g_min
        rsid = ((c.get('variant') or {}).get('rsIds') or [''])[0]
        name = sym if named else f'chr{best["chrom"]}:{best["pos"]}'
        if name in names:                           # two regions with one top gene: the later one by its position
            name = f'{sym} (chr{best["chrom"]}:{best["pos"]})'
        names.add(name)
        regions.append({'gene': name, 'genes': [sym] if named else [],
                        'ensembl': {sym: l2g['target']['id']} if named else {}, 'rsid': rsid,
                        'chrom': best['chrom'], 'pos_grch38': best['pos'], 'variant_id': c['variant']['id'],
                        # P as a number only where a float holds it; the text from mantissa and exponent
                        'p': float(c['pValueMantissa']) * 10 ** int(c['pValueExponent'])
                        if int(c['pValueExponent']) > -300 else None,
                        'p_mantissa': c['pValueMantissa'], 'p_exponent': c['pValueExponent'],
                        'evidence': 'gwas', 'code': 'G', 'source': 'Open Targets',
                        'lookup': 'Open Targets credible set',
                        'disease_study_locus_id': c['studyLocusId'],
                        'study_locus_ids': [h['cs']['studyLocusId'] for h in g], 'study_id': best['study']['id'],
                        'l2g_gene': sym, 'l2g_score': score, 'overlaps': []})
    return regions
