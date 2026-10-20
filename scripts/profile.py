#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The profile of one disease, resolved: diseases/<key>/profile.yaml
#           and the disease's optional files merged with the shared
#           configuration (config/), every default filled in, and checked.
#           Every later step reads this one file (profile.json).
#             required   id            the Open Targets disease ID
#             defaults   name, short   Open Targets name; the shortest exact
#                                      synonym in capitals (2-11 characters)
#                        pubmed        terms: the name, the short name and the
#                                      exact synonyms in title or abstract;
#                                      query: terms AND the genetic-study terms
#                        relations     curated (related_conditions.yaml) when
#                                      given, else extracted by the language-
#                                      model readers from PubMed guidelines,
#                                      meta-analyses, cohort and case-control
#                                      studies
#                        ontology      parents, their other children, children
#                        categories    four clinical categories
#                        gwas          published regions (gwas_regions.yaml),
#                                      else the disease's Open Targets credible sets
#                        themes        gene_themes.yaml, else Reactome (tidy_graph)
#                        curation      config/curation.yaml + the disease's entries
#                        choices       config/choices.yaml + the disease's items
#           The search blocks and the pre-screen come from config/search.yaml.
#           Text may name the disease as {short}, {name} or {Name}.
# Inputs  : the disease folder, the shared config folder; Open Targets only
#           when the name or short name must be read from it
# Outputs : profile.json
# ---------------------------------------------------------------------------
import argparse
import copy
import hashlib
import json
import re
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import a_, gql, load_json, log, us_english  # noqa: E402

CATEGORIES = [
    {'name': 'Precursor / spectrum', 'description': 'Conditions that precede {short} or form one disease spectrum with '
                                                    'it.'},
    {'name': 'Risk factor', 'description': 'Conditions and exposures that raise the risk of {short}.'},
    {'name': 'Comorbidity', 'description': 'Conditions more frequent in {short}, or reported with it.'},
    {'name': 'Intermediate trait', 'description': 'Measurable traits between genotype and {short}.'}]
DISEASE = '''query($id: String!) { disease(efoId: $id) { id name synonyms { relation terms } } }'''
KINDS = ('within_disease', 'other_trait', 'functional', 'background')


def read_yaml(p):
    return yaml.safe_load(Path(p).read_text()) if Path(p).exists() else None


def filler(name, short):
    Name = name[:1].upper() + name[1:]

    def fill(x):
        if isinstance(x, str):
            x = x.replace('a {short}', a_(short)).replace('A {short}', a_(short)[:1].upper() + a_(short)[1:])
            return x.replace('{short}', short).replace('{name}', name).replace('{Name}', Name)
        if isinstance(x, list):
            return [fill(v) for v in x]
        if isinstance(x, dict):
            return {k: fill(v) for k, v in x.items()}
        return x
    return fill


def merge_choices(gen, dis):
    """The shared registry with the disease's items: an item with a known key replaces that item's fields (null
    removes a field); a new item joins the group of its name (a new group goes last); references are merged."""
    out = copy.deepcopy(gen)
    where = {it['key']: (gi, ii) for gi, g in enumerate(out['groups']) for ii, it in enumerate(g['items'])}
    for g in (dis or {}).get('groups') or []:
        for it in g.get('items') or []:
            if it['key'] in where:
                gi, ii = where[it['key']]
                merged = {**out['groups'][gi]['items'][ii], **it}
                out['groups'][gi]['items'][ii] = {k: v for k, v in merged.items() if v is not None}
                continue
            grp = next((x for x in out['groups'] if x['name'] == g['name']), None)
            if grp is None:
                grp = {'name': g['name'], 'items': []}
                out['groups'].append(grp)
            grp['items'].append({k: v for k, v in it.items() if v is not None})
            where[it['key']] = (out['groups'].index(grp), len(grp['items']) - 1)
    out['references'] = {**(gen.get('references') or {}), **((dis or {}).get('references') or {})}
    return out


def merge_curation(gen, dis):
    """The shared rules with the disease's entries. The findings of each kind, the gene sets and the text-mining
    corrections come from the disease file (a stale entry there stops the build); variants and review words add up."""
    dis = dis or {}
    out = copy.deepcopy(gen)
    for k in ('publications', 'design_classes', 'countries'):
        if k in dis:
            out[k] = dis[k]
    out['variants'] = (gen.get('variants') or []) + (dis.get('variants') or [])
    F, DF = gen.get('findings') or {}, dis.get('findings') or {}
    if 'within_cteph' in DF and 'within_disease' not in DF:      # the key before the disease profile
        DF = {**DF, 'within_disease': DF['within_cteph']}
    fo = {}
    for k in KINDS:
        spec = {**(F.get(k) or {}), **{x: v for x, v in (DF.get(k) or {}).items() if x != 'items'}}
        spec['items'] = (DF.get(k) or {}).get('items') or []
        fo[k] = spec
    fo['gene_sets'] = DF.get('gene_sets') or []
    fo['review_words'] = (F.get('review_words') or []) + (DF.get('review_words') or [])
    out['findings'] = fo
    out['text_mining'] = dis.get('text_mining') or []
    out['shared_text_mining'] = gen.get('text_mining') or []
    return out


MAPIDS = '''query($q: [String!]!, $e: [String!]!) { mapIds(queryTerms: $q, entityNames: $e) {
  mappings { term hits { id name } } } }'''
STOPWORDS = {'of', 'and', 'or', 'the', 'with', 'to', 'in', 'for', 'a', 'an', 'by', 'due', 'on', 'at'}


def name_words(name):
    return [w for w in re.split(r'[\s,/-]+', re.sub(r'\(.*?\)', ' ', name)) if w and w.lower() not in STOPWORDS]


def fits(acr, words):
    """The acronym spells the name: each word gives its first letter, and may give further letters of itself, in
    order (CTEPH: Chronic ThromboEmbolic Pulmonary Hypertension)."""
    acr, words = acr.lower().replace('-', ''), [w.lower() for w in words]

    def go(i, w, j):                       # acronym position i; word w, letters from j on already free to use
        if i == len(acr):
            return w == len(words) - 1 or not words
        if w + 1 < len(words) and acr[i] == words[w + 1][0] and go(i + 1, w + 1, 1):
            return True
        k = words[w].find(acr[i], j) if w >= 0 else -1
        return k >= 0 and go(i + 1, w, k + 1)
    return bool(words) and acr[:1] == words[0][:1] and go(1, 0, 1)


# words of a disease name that name no disease in particular: not enough to tie an acronym to the disease
GENERIC_WORDS = {'disease', 'diseases', 'syndrome', 'disorder', 'type', 'familial', 'primary', 'secondary'}


def parts_of(name):
    """The diseases a grouping name joins ('A and/or B'), else the name itself."""
    return [x.strip() for x in re.split(r'\s+and/or\s+', name) if x.strip()]


MADE_UP = 'made up'      # a short name that is no name of the disease: shown in the products, never searched


def short_from(name, synonyms):
    """(short name, rule): the shortest acronym among the exact synonyms that spells the name (CTEPH, HPAH); else the
    name itself when it has at most 11 characters (Asthma); else the initials of the name's words, stop words left
    out, a short name made up here (MADE_UP in its rule): a label of the products, never searched, to be checked; ''
    when not even that works (the profile must give `short`). A grouping name ('A and/or B') takes the short names
    of its parts joined by '/' (PVOD/PCH)."""
    parts = parts_of(name)
    if len(parts) > 1:
        syn = [x for s in synonyms for x in parts_of(s)]
        shorts = [short_from(x, syn) for x in parts]
        joined = '/'.join(x for x, _r in shorts)
        if len(joined) <= 11 and all(x for x, _r in shorts):
            made = any(MADE_UP in r for _x, r in shorts)
            return joined, 'the short names of the diseases the name joins ("A and/or B"), joined by "/"' + (
                f' ({MADE_UP} from initials for a part: check it)' if made else '')
    words = name_words(name)
    acr = sorted({s for s in synonyms if re.fullmatch(r'[A-Z][A-Z0-9-]{1,10}', s) and fits(s, words)},
                 key=lambda s: (len(s), s))
    if acr:
        return acr[0], 'the shortest acronym among the exact synonyms in Open Targets that spells the name'
    if len(name) <= 11 and "'" not in name:
        return name[:1].upper() + name[1:], 'the name itself (at most 11 characters)'
    ini = ''.join(w[0] for w in words).upper()[:11]
    return (ini, f'{MADE_UP}: the initials of the name (no acronym among the exact synonyms in Open Targets spells '
                 'it); a label only, never searched: replace it with the abbreviation the literature uses') \
        if len(ini) >= 2 else ('', 'none: give `short`, the abbreviation the literature uses (at most 11 characters)')


def mapped(api, terms, entity):
    """{term: [(id, name)]}: what Open Targets maps each term to (names and synonyms, for genes also aliases)."""
    out = {}
    terms = sorted({t for t in terms if t})
    for i in range(0, len(terms), 100):
        for m in gql(api, MAPIDS, {'q': terms[i:i + 100], 'e': [entity]})['mapIds']['mappings']:
            out[m['term']] = [(h['id'], h['name']) for h in m['hits']]
    return out


def inverted(name):
    """The natural form of a name in database form 'A, B' (one comma, no parentheses): 'B A' ('hypertension,
    pulmonary' -> 'pulmonary hypertension'); '' for any other name."""
    parts = [x.strip() for x in name.split(',')]
    return f'{parts[1]} {parts[0]}' if len(parts) == 2 and all(parts) and not re.search(r'[()\[\]]', name) else ''


def search_terms(api, name, synonyms, short, own_id='', short_rule='', also=()):
    """(PubMed block, notes): the name, the exact synonyms and the short name in title or abstract, without the terms
    that would find other things: a one-word term Open Targets maps to a gene (a symbol or alias, e.g. PAH, KD), a
    term that names another disease (not the disease, one of its descendants, or an entry whose evidence it takes),
    an acronym of two characters or fewer, and database forms of a name (with parentheses, brackets, ' - ', NOS, a
    comma that does not list ('A, B and C' is a list), or in capitals with spaces). The disease's own name is always kept, whatever its form: a gene may list it
    among its names (Marfan syndrome, FBN1) without the name finding the gene; a name in database form 'A, B' is
    also searched as 'B A'. A short name made up from initials (short_rule) is never searched. A grouping name ('A
    and/or B') is searched by its parts, which name the diseases it joins; the names and exact synonyms of the
    entries whose evidence the disease takes (also: [{id, name, synonyms}]) are searched too; every term also in
    American spelling. An acronym kept is searched only beside a word of the name."""
    # a grouping name ('A and/or B') is searched by its parts, with their own short names; every term also in
    # American spelling (haemangiomatosis -> hemangiomatosis), which a phrase search does not map
    groups = [s for s in [name] + synonyms if len(parts_of(s)) > 1]
    parts = [x for s in [name] + synonyms for x in parts_of(s)]
    made_up = MADE_UP in short_rule
    shorts = [] if made_up else short.split('/') if '/' in short else [short]
    also_names = [x for r in also for x in [r['name']] + [t for s in r.get('synonyms') or []
                                                          if s.get('relation') == 'hasExactSynonym' for t in s['terms']]]
    cand = [name] + ([inverted(name)] if inverted(name) else []) + synonyms + parts + shorts + also_names
    cand = list(dict.fromkeys(cand + [us_english(x) for x in cand]))
    # the diseases a grouping name joins: the parts of the terms with 'and/or' and their short names, not every term
    gparts = [x for s in groups for x in parts_of(s)]
    joined = {x.lower() for x in gparts + [us_english(x) for x in gparts]} | (
        {x.lower() for x in short.split('/')} if groups and not made_up else set())
    own = {own_id, *(r['id'] for r in also)}
    if api and own_id:                          # the disease's descendants are the disease (its evidence includes them)
        own |= set((gql(api, RECORD, {'id': own_id})['disease'] or {}).get('descendants') or [])
    genes = mapped(api, cand, 'target') if api else {}
    dis = mapped(api, cand, 'disease') if api else {}
    seen, keep, notes = set(), [], []
    if made_up:
        notes.append(f'{short}: not searched, a short name made up from the initials of the name')
    for t in cand:
        k = t.lower()
        if k in seen:
            continue
        seen.add(k)
        other = [n for i, n in dis.get(t, []) if i not in own]
        named = t in (name, inverted(name))
        why = ('' if named else
               'also a gene in Open Targets (' + ', '.join(n for _i, n in genes[t]) + ')'
               if genes.get(t) and ' ' not in t.strip() else
               'also names another disease (' + ', '.join(other) + ')' if other and k not in joined else
               'too short to be specific' if re.fullmatch(r'[A-Za-z0-9]{1,2}', t) else
               'a database form of a name' if re.search(r'\(|\)|\[|\]| - |\bNOS\b', t) or
               (',' in t and not re.search(r',\s*(and|or)\s', t)) or (' ' in t and t == t.upper()) else '')
        if why:
            notes.append(f'{t}: left out, {why}')
        else:
            keep.append(t)
    # an acronym is searched only beside a word of the disease's name: alone it names other things too (HPAH is
    # also human phenylalanine hydroxylase, and halogenated polycyclic aromatic hydrocarbons)
    # (a generic word such as 'disease' does not tie an acronym to this disease: PCH is also pontocerebellar
    # hypoplasia, a disease); American spellings too
    nw = [w.lower() for w in name_words(name) if len(w) >= 3 and w.lower() not in GENERIC_WORDS]
    words = ' OR '.join(f'{w}[tiab]' for w in dict.fromkeys(nw + [us_english(w) for w in nw]))

    def term(t):
        if re.search(r'[\s\-/,\']', t):
            return f'"{t}"[tiab]'
        return f'({t}[tiab] AND ({words}))' if words and re.fullmatch(r'[A-Za-z0-9]+', t) and t.upper() == t \
            else f'{t}[tiab]'
    assert keep, f'{name}: no search term is left'
    block = '(' + ' OR '.join(term(t) for t in keep) + ')'
    if words and any(re.fullmatch(r'[A-Z0-9]+', t) for t in keep):
        notes.append('acronyms: searched only beside a word of the name (' + words.replace('[tiab]', '') + ')')
    return block, notes


RECORD = '''query($id: String!) { disease(efoId: $id) { id name synonyms { relation terms } descendants
  associatedTargets(page: {index: 0, size: 500}) { rows { target { id } } } } }'''
EVIDENCE_COUNT = '''query($id: String!, $ids: [String!]!) { disease(efoId: $id) {
  evidences(ensemblIds: $ids, datasourceIds: ["europepmc", "gwas_credible_sets"], enableIndirect: true, size: 1) {
    count } } }'''


def squash(text):
    """A name for comparison: American spelling, lower case, letters and digits only (veno-occlusive = venoocclusive)."""
    return re.sub(r'[^a-z0-9]', '', us_english(text).lower())


def evidence_also(api, rec):
    """(entries, notes): the Open Targets entries whose evidence is fetched with the disease's own. A grouping name
    ('A and/or B') joins diseases that Open Targets can hold as entries of their own, with their literature (PVOD/PCH:
    the text mining of pulmonary veno-occlusive disease is mapped to MONDO_0009937, not to the grouping): each part
    whose name or exact synonym is that of another entry, neither the disease nor one of its descendants (whose
    evidence is fetched anyway), adds that entry. A note says when Open Targets holds no text-mined literature or GWAS
    evidence for the disease itself (with its descendants). -> ([{id, name, part}], [note])"""
    own = gql(api, RECORD, {'id': rec['id']})['disease']
    ids = [r['target']['id'] for r in own['associatedTargets']['rows']]
    n_ev = gql(api, EVIDENCE_COUNT, {'id': rec['id'], 'ids': ids})['disease']['evidences']['count'] if ids else 0
    notes = [] if n_ev else [f'Open Targets holds no text-mined literature or GWAS evidence for {rec["id"]} (with its '
                             'descendants)']
    parts = parts_of(rec['name'])
    if len(parts) < 2:
        return [], notes
    skip = {rec['id'], *(own.get('descendants') or [])}
    found, recs = {}, {}
    for part, hits in mapped(api, parts + [us_english(x) for x in parts], 'disease').items():
        for i, _n in hits:
            if i in skip:
                continue
            r = recs.get(i) or gql(api, RECORD, {'id': i})['disease']
            recs[i] = r
            names = [r['name']] + [x for s in r.get('synonyms') or [] if s['relation'] == 'hasExactSynonym'
                                   for x in s['terms']]
            if squash(part) in {squash(x) for x in names}:
                found.setdefault(i, {})[squash(part)] = next(x for x in parts if squash(x) == squash(part))
    return [{'id': i, 'name': recs[i]['name'], 'part': '" and "'.join(sorted(v.values()))}
            for i, v in found.items()], notes


def md5(p):
    return hashlib.md5(Path(p).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--disease-dir', required=True, help='diseases/<key>')
    ap.add_argument('--config-dir', required=True, help='config/ (shared curation rules and choices)')
    ap.add_argument('--api', default='', help='Open Targets GraphQL, to read a name or short name not given')
    ap.add_argument('--saved-disease', default='', help='a saved 00.raw/disease.json (rebuild: no query)')
    ap.add_argument('--saved-search', default='', help='a saved 00.pubmed (rebuild: its queries, as they were run)')
    ap.add_argument('--out', default='profile.json')
    a = ap.parse_args()
    D, C = Path(a.disease_dir), Path(a.config_dir)
    key = D.name
    assert re.fullmatch(r'[a-z0-9_]+', key), f'{key}: a disease folder name is lower case letters, digits and _'
    prof = read_yaml(D / 'profile.yaml')
    assert prof and prof.get('id'), f'{D}/profile.yaml must give the Open Targets disease id'
    name, short = prof.get('name'), prof.get('short')
    synonyms = []
    if not (name and short) or not (prof.get('pubmed') or {}).get('terms') and not (prof.get('pubmed') or {}).get(
            'query'):
        rec = load_json(a.saved_disease) if a.saved_disease else gql(a.api, DISEASE, {'id': prof['id']})['disease'] \
            if a.api else None
        assert rec, f'{key}: the name, short name or query is not given and Open Targets cannot be read (--api)'
        synonyms = [t for s in rec.get('synonyms') or [] if s.get('relation') == 'hasExactSynonym'
                    for t in s['terms']]
        name = name or rec['name']
        short = short or short_from(rec['name'], synonyms)[0]
    assert short, f'{key}: no short name: give `short` in profile.yaml (the abbreviation the literature uses)'
    assert len(short) <= 11 and "'" not in short, f'{key}: short name {short!r}: at most 11 characters, no apostrophe'
    fill = filler(name, short)
    search = read_yaml(C / 'search.yaml')
    pm = prof.get('pubmed') or {}
    # a rebuild records the queries that were run (saved with the snapshot's searches), never new ones
    saved = {}
    for what, f in (('genetics', 'pubmed_search.json.gz'), ('relations', 'relations_search.json.gz')):
        if a.saved_search and (Path(a.saved_search) / f).exists():
            saved[what] = load_json(Path(a.saved_search) / f).get('query') or ''
    tail = f' AND {search["genetic_terms"]}'
    if saved.get('genetics') and not pm.get('terms'):
        terms = saved['genetics'][:-len(tail)] if saved['genetics'].endswith(tail) else saved['genetics']
    else:
        also = [gql(a.api, DISEASE, {'id': i})['disease'] for i in prof.get('evidence_also') or []] if a.api else []
        terms = pm.get('terms') or search_terms(a.api, name, synonyms, short, prof['id'],
                                                '' if prof.get('short') else short_from(name, synonyms)[1],
                                                [x for x in also if x])[0]
    query = pm.get('query') or saved.get('genetics') or f'{terms} AND {search["genetic_terms"]}'
    seed = read_yaml(D / 'related_conditions.yaml')
    if seed is not None:
        seed['facts'] = seed.pop('facts', None) or seed.pop('cteph_facts', None) or []
    rel = prof.get('relations') or {}
    relations = {'mode': 'curated' if seed else 'extracted',
                 'query': rel.get('query') or saved.get('relations') or
                 f'{terms} AND {search["relation_designs"]} AND {search["relation_words"]}',
                 'query_source': 'profile' if rel.get('query') else 'default'}
    onto = {'parents': True, 'siblings': True, 'children': True, 'children_of': [], 'category': 'Disease group',
            'description': 'The neighbors of {short} in the disease ontology.', **(prof.get('ontology') or {})}
    gwas = read_yaml(D / 'gwas_regions.yaml')
    th = read_yaml(D / 'gene_themes.yaml')
    texts = {p.stem: read_yaml(p) for p in sorted((D / 'text').glob('*.yaml'))} if (D / 'text').is_dir() else {}
    files = sorted([p for p in D.rglob('*.yaml')] + [C / 'curation.yaml', C / 'choices.yaml', C / 'search.yaml'])
    out = {'key': key, 'id': prof['id'], 'name': name, 'short': short, 'Name': name[:1].upper() + name[1:],
           'stem': f'{key}_knowledge_graph',
           'pubmed_terms': terms, 'pubmed_query': query,
           'query_source': 'profile' if pm.get('query') else 'default',
           'scope': fill(pm.get('scope') or ''), 'prompt_name': prof.get('prompt_name') or '',
           'relations': relations, 'ontology': fill(onto),
           # other Open Targets entries whose evidence is fetched with the disease's own (new_disease.py proposes them)
           'evidence_also': [str(x) for x in prof.get('evidence_also') or []],
           'categories': fill(prof.get('categories') or CATEGORIES),
           'labels': prof.get('labels') or {}, 'abbreviations': prof.get('abbreviations') or {},
           'seed': seed, 'facts': (seed or {}).get('facts') or [], 'prescreen': search['prescreen'],
           'gwas': gwas, 'themes': th['themes'] if th else None,
           'curation': fill(merge_curation(read_yaml(C / 'curation.yaml'), read_yaml(D / 'curation.yaml'))),
           'choices': fill(merge_choices(read_yaml(C / 'choices.yaml'), read_yaml(D / 'choices.yaml'))),
           'checks': prof.get('checks') or {}, 'caveats': fill(prof.get('caveats') or []), 'text': fill(texts),
           'files': {str(p.relative_to(D.parent.parent) if D.parent.parent in p.parents else p.name): md5(p)
                     for p in files}}
    Path(a.out).write_text(json.dumps(out, indent=1, ensure_ascii=False))
    log(f'{key}: {prof["id"]} {name} ({short}); relations {relations["mode"]}; GWAS regions '
        f'{"published" if gwas else "from Open Targets"}; themes {"curated" if th else "Reactome"}; '
        f'{sum(len(g["items"]) for g in out["choices"]["groups"])} choices')


if __name__ == '__main__':
    main()
