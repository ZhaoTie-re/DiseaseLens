#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The clinical relations of a disease without a curated list, read
#           from PubMed: the records of the profile's relation search
#           (guidelines, consensus statements, meta-analyses, systematic
#           reviews, cohort and case-control studies whose abstract names the
#           disease and a relation word) are read by the two language-model
#           readers, independently; per record, the relations (condition,
#           direction and category) are compared, and an adjudicating call
#           decides where they disagree. Only relations between a condition and
#           the disease in human patients or populations count (not a sign or
#           finding of the disease, a mechanism, an animal or cell experiment,
#           or a comparison between groups of patients). Every relation keeps a verbatim quote of the abstract that
#           names the condition; a relation whose quote or condition is not in
#           the abstract is dropped and logged. Each condition is mapped to an
#           Open Targets disease whose name matches exactly (else it has no
#           genetics); a condition with one supporting publication is marked.
#           The result is a seed in the schema of a curated list
#           (related_conditions.yaml), which resolve.py checks quote by quote.
# Inputs  : the relation-search records (pubmed.py --mode search), the profile,
#           the model configuration
# Outputs : extracted_relations.yaml (the seed), relations_pubmed.json.gz (the
#           records it cites), relations_review.tsv (disputed items),
#           relations_log.tsv (dropped items), relations.json (counts)
# ---------------------------------------------------------------------------
import argparse
import json
import re
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from llm import EAST_ASIA_PLACES, in_text, in_words, load_llm, norm, per_record, reconcile, record_text  # noqa: E402
from ot_common import (EAS_WORDS, countries_named, gql, load_json, load_profile, log, release, us_english,  # noqa: E402
                       write_json, write_tsv)

DIRECTIONS = ['increased risk', 'decreased risk', 'associated', 'no association found']
MAPIDS = '''query($q: [String!]!) { mapIds(queryTerms: $q, entityNames: ["disease"]) {
  mappings { term hits { id name } } } }'''


def system_prompt(P):
    cats = [c['name'] for c in P.categories]
    return f'''You read PubMed records about {P.name} ({P.short}) and extract the CONDITIONS each record relates to
{P.short}: diseases, traits, exposures or procedures that precede it or form one spectrum with it, raise its risk,
occur with it more often than expected, or lie between genotype and disease. Use ONLY the text given; never add outside
knowledge. Copy names, numbers and quotes exactly as written in the record. Answer with JSON only (no prose).

For every record return one object:
{{"pmid": "<PMID>",
  "relations": [{{"condition": "<the condition exactly as the record writes it>",
                 "standard_name": "<its usual disease or trait name, e.g. systemic sclerosis; empty if unsure>",
                 "category": one of {json.dumps(cats)},
                 "direction": one of {json.dumps(DIRECTIONS)},
                 "design": "<the design the record states: guideline, meta-analysis, cohort, case-control, ...>",
                 "population": "<country or setting and size, as stated; empty if not stated>",
                 "statement": "<one sentence: how the condition relates to {P.short}, with the numbers the record gives>",
                 "quote": "<one exact sentence or clause of the record that states the relation>"}}]}}
Only relations the record states for {P.short} itself, between the condition and having {P.short}: for a group of
patients only when {P.short} patients are most of it, or when the record gives the result for {P.short} patients
apart (not for a broader group in which {P.short} is a minority, such as all patients with a related disease). Not a
relation: a consequence or complication of {P.short} or of its treatment, an outcome, a prognostic marker or a measurement among patients, a sign, symptom, imaging, histological or laboratory
finding of {P.short} itself, a molecular or cellular process (such as inflammation or a pathway), a treatment or a
drug given as therapy, a gene, a demographic attribute (sex, age, race or ethnicity), {P.short} itself or one of its
own subtypes, and a comparison between groups of {P.short} patients (such as carriers and non-carriers of a variant).
Only evidence from human patients or populations: not from an animal model or a cell or tissue experiment. Return a JSON
array, one object per record, in the order given; an empty relations list when a record states none. Quotes are copied
whole, never shortened with '...'.''' + (
        f'\nWhich patients count as {P.short}: {P.d["scope"]}' if P.d.get('scope') else '')


STOP = {'of', 'the', 'a', 'an', 'and', 'or', 'in', 'with', 'to', 'for', 'by', 'on', 'at', 'as'}


def named_in(cond, quote):
    """The quote names the condition: its full name; its name without parentheses; an abbreviation it gives in
    parentheses (a whole word: 'mitomycin-C (MMC) therapy' is named by 'MMC therapy'); or every word of its name
    but those that say how it was met (EXPOSURE_WORDS) and short connecting words, American spelling, a plural
    allowed ('fenfluramine exposure' is named by 'exposure to fenfluramine'). A different term is not the condition
    ('thyroid disease' is not named by 'thyrotoxicosis')."""
    if not cond or not quote:
        return False
    if in_text(cond, quote):
        return True
    outer = ' '.join(re.sub(r'\([^()]*\)', ' ', cond).split())
    if outer and in_text(outer, quote):
        return True
    if any(in_words(x.strip(), quote, case=False) for x in re.findall(r'\(([^()]*)\)', cond) if x.strip()):
        return True
    words = [w for w in re.split(r'[^a-z0-9]+', us_english(norm(outer))) if w and w not in EXPOSURE_WORDS | STOP]
    q = us_english(norm(quote))
    return bool(words) and all(re.search(rf'(?<![a-z0-9]){re.escape(w)}', q) for w in words)


def ground(ans, r, cats, logs, reader):
    """The grounded reading of one answer, {'relations': [...]}: quote and condition in the abstract."""
    text = r['title'] + ' ' + r['abstract']
    out = []
    rels = (ans.get('relations') if isinstance(ans, dict) else None) or []
    for x in rels if isinstance(rels, list) else []:
        if not isinstance(x, dict):
            logs.append({'pmid': r['pmid'], 'item': str(x)[:80], 'why': 'not a relation object', 'value': '',
                         'reader': reader})
            continue
        x = {k: v if isinstance(v, str) else '' if v is None else str(v) for k, v in x.items()}
        cond = x.get('condition', '')
        why = ('quote not in the abstract' if not in_text(x.get('quote', ''), text) else
               'condition not in its quote' if not (named_in(cond, x.get('quote', '')) or
                                                    named_in(x.get('standard_name', ''), x.get('quote', ''))) else
               'category not allowed' if x.get('category') not in cats else
               'direction not allowed' if x.get('direction') not in DIRECTIONS else '')
        if why:
            logs.append({'pmid': r['pmid'], 'item': x.get('condition', ''), 'why': why,
                         'value': str(x.get('quote', ''))[:200], 'reader': reader})
            continue
        out.append({k: str(x.get(k) or '') for k in ('condition', 'standard_name', 'category', 'direction', 'design',
                                                     'population', 'statement', 'quote')})
    return {'relations': out}


def relation_differences(ra, rb):
    """The relations on which two readings of a record differ (a different category is a difference)."""
    ka, kb = {rkey(x) for x in ra['relations']}, {rkey(x) for x in rb['relations']}
    return [('relation', rtext(k), 'not read') for k in sorted(ka - kb)] + \
        [('relation', 'not read', rtext(k)) for k in sorted(kb - ka)]


def rkey(x):
    """A relation as the readers are compared on it: the condition (its usual name, else as written), the
    direction (any positive direction counts as one) and the category."""
    name = re.sub(r'[^a-z0-9]+', ' ', norm(x.get('standard_name') or x.get('condition'))).strip()
    d = x.get('direction')
    return name, 'positive' if d in ('increased risk', 'associated') else d, x.get('category') or ''


def rtext(k):
    """A compared relation as text: 'condition (direction; category)'."""
    return f'{k[0]} ({k[1]}; {k[2]})'


RECORDS = '''query($ids: [String!]!) { diseases(efoIds: $ids) { id name synonyms { relation terms } } }'''


def squash(text):
    """A name for comparison: American spelling, lower case, letters and digits only (Tobacco-smoking = tobacco
    smoking)."""
    return re.sub(r'[^a-z0-9]', '', us_english(norm(text)))


# words that say how a condition was met, not what it is: 'fenfluramine use' and 'fenfluramine exposure' are one
EXPOSURE_WORDS = {'use', 'uses', 'user', 'users', 'exposure', 'exposures', 'exposed', 'intake', 'consumption',
                  'therapy', 'treatment', 'administration', 'of'}


def name_forms(name):
    """The forms a condition name can take for merging: the name and each parenthetical alias, in American spelling,
    lower case, without the words of EXPOSURE_WORDS, letters and digits only. An alias that is an abbreviation (one
    word of at most 5 letters or digits, or in capitals) is not a form: one abbreviation names several conditions
    (PE: pulmonary embolism, pleural effusion)."""
    def form(s):
        words = [w for w in re.split(r'[^a-z0-9]+', us_english(norm(s))) if w and w not in EXPOSURE_WORDS]
        return ''.join(words)

    def abbreviation(s):
        s = s.strip()
        return (' ' not in s and len(re.sub(r'[^A-Za-z0-9]', '', s)) <= 5) or s == s.upper()
    outer = re.sub(r'\([^()]*\)', ' ', name)
    aliases = [x for x in re.findall(r'\(([^()]*)\)', name) if not abbreviation(x)]
    return {f for f in [form(outer)] + [form(x) for x in aliases] if f}


def merge_names(names):
    """{name: the name it is merged into}: names that share a form (name_forms) are one condition, named by its
    shortest name (ties: alphabetical)."""
    parent = {n: n for n in names}

    def find(n):
        while parent[n] != n:
            parent[n] = parent[parent[n]]
            n = parent[n]
        return n
    by_form = {}
    for n in sorted(names):
        for f in name_forms(n):
            if f in by_form:
                parent[find(n)] = find(by_form[f])
            else:
                by_form[f] = n
    groups = {}
    for n in names:
        groups.setdefault(find(n), []).append(n)
    return {n: min(g, key=lambda x: (len(x), x)) for g in groups.values() for n in g}


def map_conditions(api, terms):
    """{term: (id, name)} for terms that Open Targets maps to a disease whose name or exact synonym is that term
    (case, spaces, hyphens and spelling aside): connective tissue disease is an exact synonym of an entry named
    otherwise. Names mapped to one entry become one condition."""
    out = {}
    terms = sorted({t for t in terms if t})
    maps = []
    for i in range(0, len(terms), 100):
        maps += gql(api, MAPIDS, {'q': terms[i:i + 100]})['mapIds']['mappings']
    ids = sorted({h['id'] for m in maps for h in m['hits']})
    names = {}
    for i in range(0, len(ids), 100):
        for r in gql(api, RECORDS, {'ids': ids[i:i + 100]})['diseases']:
            names[r['id']] = {squash(r['name'])} | {squash(x) for s in r.get('synonyms') or []
                                                     if s['relation'] == 'hasExactSynonym' for x in s['terms']}
    for m in maps:
        hit = next((h for h in m['hits'] if squash(m['term']) in names.get(h['id'], {squash(h['name'])})), None)
        if hit:
            out[m['term']] = (hit['id'], hit['name'])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--records', required=True, help='pubmed.py --mode search output of the relation search')
    ap.add_argument('--profile', required=True)
    ap.add_argument('--llm', required=True)
    ap.add_argument('--offline', default='false')
    ap.add_argument('--batch', type=int, default=5)
    ap.add_argument('--limit', type=int, default=0, help='smoke run: the first N records only')
    ap.add_argument('--api', required=True, help='Open Targets GraphQL (condition names)')
    ap.add_argument('--out-dir', default='.')
    a = ap.parse_args()
    P = load_profile(a.profile)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    recs = [r for r in load_json(a.records)['records'] if r['abstract']]
    if a.limit:
        recs = recs[:a.limit]
    cats = [c['name'] for c in P.categories]
    C = P.curation.get('countries') or {}
    EAS_SET = {n for n in C.get('names') or [] if n in EAS_WORDS}
    SYS = system_prompt(P)
    readers, adj = load_llm(a.llm, offline=a.offline.lower() in ('true', '1', 'yes'))
    names = [r.name for r in readers]
    logs = []
    by_rec = {r['pmid']: r for r in recs}
    readings = per_record(readers, recs, SYS, a.batch, lambda x, r, name: ground(x, r, cats, logs, name), logs)
    fin, review = reconcile(readers, adj, [(r['pmid'], 'Record:\n\n' + record_text(r), readings[r['pmid']])
                                           for r in recs], SYS, relation_differences,
                            lambda ans, key, _t: ground(ans, by_rec[key], cats, logs, 'adjudicator'))
    final = {}
    for pm, (reading, why) in fin.items():
        rs = readings[pm]
        both = ({rkey(x) for x in rs[0]['relations']} & {rkey(x) for x in rs[1]['relations']}) \
            if why == 'adjudicated' and len(rs) == 2 and None not in rs else set()
        final[pm] = [dict(x, agreement=why if why != 'adjudicated' else ('both' if rkey(x) in both else 'adjudicated'))
                     for x in reading['relations']]
    for row in review:
        row['pmid'] = row.pop('key')
        if row['pmid'] in fin:
            row['final'] = '; '.join(rtext(k) for k in sorted({rkey(x) for x in fin[row['pmid']][0][
                'relations']})) or 'none'

    # ── conditions: mapped to Open Targets, grouped, the disease itself left out ──
    rows = [{'pmid': p, **x} for p, xs in final.items() for x in xs]
    terms = {t for x in rows for t in (x['standard_name'], x['condition'])}
    rel = release(a.api)
    mapped = map_conditions(a.api, terms) if terms else {}
    own = {norm(P.name), norm(P.short)}
    own_ids = {P.id, *(P.d.get('evidence_also') or [])}     # the disease, and the entries whose evidence it takes
    # conditions without an Open Targets entry: names of one condition (use / exposure, an alias in parentheses) are
    # merged, and the merge is logged
    custom = {x['standard_name'] or x['condition'] for x in rows
              if not (mapped.get(x['standard_name']) or mapped.get(x['condition']))}
    merged = merge_names(custom)
    for n, m in sorted(merged.items()):
        if n != m:
            logs.append({'pmid': '', 'item': n, 'why': f'merged with "{m}" (one condition)', 'value': '', 'reader': ''})
    conds = {}
    for x in rows:
        hit = mapped.get(x['standard_name']) or mapped.get(x['condition'])
        name = hit[1] if hit else merged[x['standard_name'] or x['condition']]
        if (hit and hit[0] in own_ids) or norm(name) in own or norm(x['condition']) in own:
            logs.append({'pmid': x['pmid'], 'item': x['condition'], 'why': f'{P.short} itself', 'value': '',
                         'reader': ''})
            continue
        key = hit[0] if hit else re.sub(r'[^a-z0-9]+', '_', norm(name)).strip('_')     # resolve: custom:<key>
        c = conds.setdefault(key, {'name': name[:1].upper() + name[1:], 'ot_ids': [hit[0]] if hit else [],
                                   'categories': Counter(), 'evidence': [], 'east_asian': []})
        c['categories'][x['category']] += 1
        r = by_rec[x['pmid']]
        e = {'pmid': x['pmid'], 'direction': x['direction'], 'design': x['design'], 'population': x['population'],
             'statement': x['statement'], 'quote': x['quote'], 'agreement': x['agreement']}
        # an East Asian population: the final reading's population names an East Asian country, people or place,
        # and so does the abstract (whole words; not necessarily the same: 'China' and 'Peking'), and neither names
        # a country outside East Asia (a multinational study is not East Asian)
        text = r['title'] + ' ' + r['abstract']
        words = EAS_WORDS + EAST_ASIA_PLACES
        ea = (any(in_words(w, x['population']) for w in words) and any(in_words(w, text) for w in words)
              and not (countries_named(x['population'], C) | countries_named(text, C)) - EAS_SET)
        (c['east_asian'] if ea else c['evidence']).append(e)
        c.setdefault('years', set()).add(r['year'])
    seed = {'version': load_json(a.records).get('max_date', ''),
            'note': 'Extracted: read from PubMed by two independent language-model readers, adjudicated where they '
                    'disagreed; every quote is checked against its abstract at each build.',
            'source': 'extracted', 'conditions': [], 'facts': []}
    for key, c in sorted(conds.items(), key=lambda kv: kv[1]['name'].lower()):
        n_pub = len({e['pmid'] for e in c['evidence'] + c['east_asian']})
        seed['conditions'].append({
            'key': key.lower(), 'name': c['name'], 'category': c['categories'].most_common(1)[0][0],
            'ot_ids': c['ot_ids'], 'genetics': 'full' if c['ot_ids'] else 'none', 'mendelian': bool(c['ot_ids']),
            'genes': [], 'support': 'one publication' if n_pub == 1 else f'{n_pub} publications',
            'evidence': c['evidence'], 'east_asian': c['east_asian']})
    (out / 'extracted_relations.yaml').write_text(yaml.safe_dump(seed, sort_keys=False, allow_unicode=True, width=118))
    used = sorted({e['pmid'] for c in seed['conditions'] for e in c['evidence'] + c['east_asian']}, key=int)
    write_json({'mode': 'seed', 'n': len(used), 'records': [by_rec[p] for p in used]}, out / 'relations_pubmed.json.gz')
    write_tsv(pd.DataFrame(review, columns=['pmid', 'item'] + names + ['final', 'decision']),
              out / 'relations_review.tsv')
    write_tsv(pd.DataFrame(logs, columns=['pmid', 'item', 'why', 'value', 'reader']).fillna(''),
              out / 'relations_log.tsv')
    ag = Counter(x['agreement'] for x in rows)
    not_read = sum(all(x is None for x in readings[r['pmid']]) for r in recs)
    unanswered = sum(1 for x in logs if x.get('item') == 'record' and x.get('value') == 'no answer')
    stats = {'records': len(recs), 'read': len(recs) - not_read, 'not_read': not_read,
             'relations': len(rows), 'conditions': len(seed['conditions']),
             'mapped_to_open_targets': sum(1 for c in seed['conditions'] if c['ot_ids']),
             'one_publication': sum(1 for c in seed['conditions'] if c['support'] == 'one publication'),
             'agreement': dict(ag), 'disputed_items': len(review), 'dropped_items': len(logs) - unanswered,
             'open_targets_release': rel['data_version'],
             'readers': [r.stats() for r in readers], 'adjudicator': adj.stats() if adj else None}
    (out / 'relations.json').write_text(json.dumps(stats, indent=1))
    log(f'{len(recs)} records: {len(rows)} relations ({dict(ag)}), {len(seed["conditions"])} conditions '
        f'({stats["mapped_to_open_targets"]} in Open Targets, {stats["one_publication"]} on one publication); '
        f'{len(review)} disputed, {len(logs) - unanswered} dropped' + (f', {not_read} not read' if not_read else ''))


if __name__ == '__main__':
    main()
