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
#           The conditions are then consolidated: the two readers decide, for
#           each condition of the list, to keep it, merge it with another of the
#           list (the same condition under another name), attach it to another
#           as a feature of that condition's patients (unprovoked or massive
#           pulmonary embolism), or drop it (not a condition, or a consequence of
#           the disease); the adjudicator decides where they differ.
#           The result is a seed in the schema of a curated list
#           (related_conditions.yaml), which resolve.py checks quote by quote.
# Inputs  : the relation-search records (pubmed.py --mode search), the profile,
#           the model configuration
# Outputs : extracted_relations.yaml (the seed), relations_pubmed.json.gz (the
#           records it cites), relations_review.tsv (disputed items),
#           relations_log.tsv (dropped items), relations_consolidation.tsv (the
#           decision on each condition), relations.json (counts)
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
from llm import (EAST_ASIA_PLACES, decide, in_text, in_words, load_llm, norm, per_record, reconcile,  # noqa: E402
                 record_text)
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


def words_of(text):
    """The words of a text for name matching: lower case, American spelling, split at anything not a letter or digit."""
    return [w for w in re.split(r'[^a-z0-9]+', us_english(norm(text))) if w]


def spells(abbr, words):
    """The words spell the abbreviation in order: each gives its first letter and, after it, none or more of its later
    letters in order; a short connecting word may give none ('FVIII' is spelled by factor VIII, 'VTE' by venous
    thromboembolism, 'IL6' by interleukin 6; 'CT' is not spelled by clotting factor)."""
    a = [c for c in abbr.lower() if c.isalnum()]

    def go(i, j):
        if j == len(words):
            return i == len(a)
        w = words[j]
        if w in STOP and go(i, j + 1):
            return True
        if i >= len(a) or w[0] != a[i]:
            return False
        ends, k, p = [i + 1], i + 1, 1
        while k < len(a) and (p := w.find(a[k], p) + 1):
            k += 1
            ends.append(k)
        return any(go(k, j + 1) for k in reversed(ends))
    return bool(a) and bool(words) and go(0, 0)


def run_at(small, big):
    """Where the words of small occur as a run in big; -1 if they do not."""
    n = len(small)
    return next((i for i in range(len(big) - n + 1) if n and big[i:i + n] == small), -1)


def defined_abbreviations(text):
    """{abbreviation: long form words} that a text defines as 'long form (ABBR)': an abbreviation of 2 to 10 letters
    or digits with two capitals or more, its long form the fewest words before it that spell it."""
    out = {}
    for m in re.finditer(r'\(\s*([A-Za-z][A-Za-z0-9-]{1,9})\s*[,;)]', text):
        abbr = m.group(1)
        if sum(c.isupper() for c in abbr) < 2 or abbr in out:
            continue
        before = words_of(text[max(0, m.start() - 300):m.start()])[-(len(abbr) + 5):]
        k = next((k for k in range(1, len(before) + 1) if spells(abbr, before[-k:])), 0)
        if k:
            out[abbr] = before[-k:]
    return out


def named_by_abbreviation(cond, quote, text, defined):
    """The quote names the condition by an abbreviation (a word with two capitals or more) for part of its name: the
    long form the record defines for it ('factor VIII (FVIII)'), or, when the record defines none, two or more words
    of the name that spell it and are written out elsewhere in the record ('FVIII' and a title on clotting factor
    VIII). The rest of the name must be in the record too: 'factor VIII deficiency' is not named by 'FVIII' in a
    record on high factor VIII."""
    c = words_of(cond)
    if not c or not quote:
        return False
    t = ' '.join(words_of(text))
    for a in {x for x in re.findall(r'(?<![A-Za-z0-9])[A-Za-z][A-Za-z0-9]{1,9}(?![A-Za-z0-9])', quote)
              if sum(ch.isupper() for ch in x) >= 2}:
        own = next((c[i:j] for n in range(2, len(c) + 1) for i in range(len(c) - n + 1) for j in [i + n]
                    if c[i] not in STOP and c[j - 1] not in STOP and spells(a, c[i:j])), None)
        if a in defined:
            long = defined[a] if run_at(defined[a], c) >= 0 else own if own and run_at(own, defined[a]) >= 0 else None
        else:
            phrase = re.escape(' '.join(own or []))
            long = own if own and re.search(rf'(?<![a-z0-9]){phrase}(?![a-z0-9])', t) else None
        if not long:
            continue
        i = run_at(long, c)
        rest = [w for w in c[:i] + c[i + len(long):] if w not in EXPOSURE_WORDS | STOP]
        if all(re.search(rf'(?<![a-z0-9]){re.escape(w)}', t) for w in rest):
            return True
    return False


def ground(ans, r, cats, logs, reader):
    """The grounded reading of one answer, {'relations': [...]}: quote and condition in the abstract."""
    text = r['title'] + ' ' + r['abstract']
    defined = defined_abbreviations(text)
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
               'condition not in its quote' if not any(
                   named_in(n, x.get('quote', '')) or named_by_abbreviation(n, x.get('quote', ''), text, defined)
                   for n in (cond, x.get('standard_name', '')) if n) else
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


def custom_key(name):
    """The key of a condition without an Open Targets entry: the words of its name in American spelling, lower case,
    without the words that say how it was met (EXPOSURE_WORDS), so that one condition has one key in every disease
    ('Fenfluramine use' and 'Fenfluramine exposure': fenfluramine)."""
    words = [w for w in re.split(r'[^a-z0-9]+', us_english(norm(name))) if w and w not in EXPOSURE_WORDS]
    return '_'.join(words) or re.sub(r'[^a-z0-9]+', '_', norm(name)).strip('_')


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


# a level as the readers write it ('Elevated factor VIII', 'Low plasma ADAMTS13 level'): the words of its direction
# and of the measure, around the trait Open Targets names as '<trait> measurement', '<trait> level' or '<trait> amount'
DIRECTION_WORDS = r'^(?:elevated|high|higher|increased|raised|low|lower|decreased|reduced|altered|abnormal)\s+'
MEASURE_WORDS = r'(?:\s+(?:plasma|serum|blood|circulating))?'
TRAIT_TAILS = r'\s+(?:levels?|concentrations?|positivity|activity|amount|measurements?)$'


def term_forms(term):
    """The forms a condition name is looked up under, in order: the name; and, for a level written with its direction
    or its measure, the trait with each of the endings Open Targets gives measured traits."""
    t = term.strip()
    core = re.sub(TRAIT_TAILS, '', re.sub(DIRECTION_WORDS + MEASURE_WORDS, '', t, flags=re.I), flags=re.I).strip()
    core = re.sub(r'^(?:plasma|serum|blood|circulating)\s+', '', core, flags=re.I)
    if not core or core.lower() == t.lower():
        return [t]
    return [t, core] + [f'{core} {x}' for x in ('measurement', 'level', 'amount')]


def map_conditions(api, terms):
    """{term: (id, name)} for terms that Open Targets maps to a disease whose name or exact synonym is that term, or
    one of its forms (term_forms; case, spaces, hyphens and spelling aside): connective tissue disease is an exact
    synonym of an entry named otherwise. Names mapped to one entry become one condition."""
    out = {}
    terms = sorted({t for t in terms if t})
    forms = {t: term_forms(t) for t in terms}
    queries = sorted({f for fs in forms.values() for f in fs})
    maps = []
    for i in range(0, len(queries), 100):
        maps += gql(api, MAPIDS, {'q': queries[i:i + 100]})['mapIds']['mappings']
    ids = sorted({h['id'] for m in maps for h in m['hits']})
    names = {}
    for i in range(0, len(ids), 100):
        for r in gql(api, RECORDS, {'ids': ids[i:i + 100]})['diseases']:
            names[r['id']] = {squash(r['name'])} | {squash(x) for s in r.get('synonyms') or []
                                                     if s['relation'] == 'hasExactSynonym' for x in s['terms']}
    hit_of = {}
    for m in maps:
        hit = next((h for h in m['hits'] if squash(m['term']) in names.get(h['id'], {squash(h['name'])})), None)
        if hit:
            hit_of[m['term']] = (hit['id'], hit['name'])
    for t, fs in forms.items():
        h = next((hit_of[f] for f in fs if f in hit_of), None)
        if h:
            out[t] = h
    return out


# ── consolidation: one decision per condition of the list ───────────────────
ACTIONS = ('keep', 'merge', 'feature', 'drop')
CONSOLIDATION_COLUMNS = ['condition', 'open_targets', 'publications', 'action', 'into', 'agreement', 'reason',
                         'reader_a', 'reader_b']


def consolidation_prompt(P):
    s = P.short
    return f"""You review the conditions that two readers found related to {P.name} ({s}) in PubMed records: one item per
condition, each with its id, its name, whether Open Targets holds an entry for it, its category, the number of
publications that relate it to {s}, and up to three of their statements. Decide one action for every item:
- "keep": a condition of its own (a disease, trait, exposure or procedure) related to {s}.
- "merge", with "into": the id of the item that is the same condition under another name: a synonym or other wording,
  or the antibody, laboratory or genetic marker that defines it (antiphospholipid antibody positivity is
  antiphospholipid syndrome). A narrower disease of its own is not merged into a broader one (breast cancer stays
  apart from cancer).
- "feature", with "into": the id of the item whose patients it describes, when it is not a condition of its own but a
  characteristic of that condition which the record relates to {s}: its severity, type, location, extent,
  provocation, recurrence or persistence, a finding at its diagnosis (imaging, laboratory, hemodynamic), or its
  treatment (unprovoked or massive pulmonary embolism, and right ventricular dysfunction at the acute embolism, are
  features of pulmonary embolism).
- "drop": not a condition (a group or placeholder, such as "comorbid diseases" or "risk factors"), or a consequence of
  {s} itself (a complication, symptom or effect among {s} patients, such as their anxiety).
Use only the names and statements given. When unsure, keep. Answer with JSON only:
{{"decisions": [{{"id": <id>, "action": "keep" | "merge" | "feature" | "drop", "into": <id, or null>,
"reason": "<at most 15 words>"}}]}}, one decision per item."""


def votes(c):
    """The categories of a condition's evidence rows, counted: its own rows that relate it (not a finding of no
    association, and not a row attached as a feature of its patients), else its own rows, else every row. The
    condition takes the most frequent (its category)."""
    rows = c['evidence'] + c['east_asian']
    own = [e for e in rows if not e.get('feature')]
    rel = [e for e in own if e.get('direction') != 'no association found'] or own or rows
    return Counter(e.get('category') or '' for e in rel)


def consolidation_list(conds):
    """The conditions in a fixed order, and the list the readers are given (id, name, Open Targets entry, category,
    publications, the first three statements)."""
    keys = sorted(conds, key=lambda k: (conds[k]['name'].lower(), k))
    lines = []
    for i, k in enumerate(keys, 1):
        c = conds[k]
        ev = c['evidence'] + c['east_asian']
        lines.append(f'[{i}] {c["name"]} | Open Targets: {"yes" if c["ot_ids"] else "no"} | category: '
                     f'{votes(c).most_common(1)[0][0]} | publications: {len({e["pmid"] for e in ev})}\n'
                     + '\n'.join(f'    - {e["statement"][:240]}' for e in ev[:3]))
    return keys, '\n'.join(lines)


def ground_decisions(ans, ids, who, n, logs, key='decisions'):
    """{id: (action, into, reason)} of an answer: one of ACTIONS; a merge or a feature names another item of the
    list (1..n). A decision that is not one of these is logged, and its item is left as it was."""
    out = {}
    ds = (ans.get(key) if isinstance(ans, dict) else None) or []
    for d in ds if isinstance(ds, list) else []:
        if not isinstance(d, dict):
            continue
        try:
            i = int(d.get('id'))
        except (TypeError, ValueError):
            continue
        if i not in ids:
            continue
        act = str(d.get('action') or '').strip().lower()
        try:
            into = int(d['into']) if d.get('into') not in (None, '', 'null') else None
        except (TypeError, ValueError):
            into = None
        if act not in ACTIONS or (act in ('merge', 'feature') and (into is None or not 1 <= into <= n or into == i)):
            logs.append({'pmid': '', 'item': f'list item {i}', 'why': 'decision not allowed: left as it was',
                         'value': json.dumps(d)[:200], 'reader': who})
            continue
        out[i] = (act, into if act in ('merge', 'feature') else None, str(d.get('reason') or '')[:300])
    return out


def apply_decisions(conds, keys, final, agree, logs, step):
    """The conditions after one decision per condition (keys: the conditions by id 1..n; final: {id: (action, into,
    reason)}, an id without a decision is kept). Each merge or feature follows its "into" to an item that is kept;
    one that ends at a dropped item or loops is not applied (the item is kept, and logged). Merged items are one
    condition, named by its member with an Open Targets entry and the most publications (else the kept item), with
    every member's Open Targets IDs; a feature's evidence joins the condition it describes. Every moved statement
    starts with the name it was read under. -> (conditions, final decisions as applied)"""
    n = len(keys)
    fin = {i: final.get(i, ('keep', None, '')) for i in range(1, n + 1)}
    name = lambda i: conds[keys[i - 1]]['name']  # noqa: E731

    def dest(i):
        seen, feature = {i}, False
        while fin[i][0] in ('merge', 'feature'):
            feature |= fin[i][0] == 'feature'
            i = fin[i][1]
            if i in seen or fin[i][0] == 'drop':
                return None, False
            seen.add(i)
        return i, feature
    groups, features = {}, {}
    for i in range(1, n + 1):
        act = fin[i][0]
        if act == 'drop':
            logs.append({'pmid': '', 'item': name(i), 'why': f'dropped by the {step}: {fin[i][2]}', 'value': '',
                         'reader': agree.get(i, '')})
            continue
        d, feat = dest(i) if act != 'keep' else (i, False)
        if d is None:
            logs.append({'pmid': '', 'item': name(i), 'why': f'{step} not applied ({act} into "{name(fin[i][1])}"): '
                         'its target is dropped or the decisions loop; kept', 'value': '', 'reader': agree.get(i, '')})
            fin[i], d, feat = ('keep', None, fin[i][2]), i, False
        (features if feat else groups).setdefault(d, []).append(i)
    pubs = lambda c: len({e['pmid'] for e in c['evidence'] + c['east_asian']})  # noqa: E731
    out = {}
    for d, members in groups.items():
        with_ot = [i for i in members if conds[keys[i - 1]]['ot_ids']]
        rep_i = max(with_ot, key=lambda i: (pubs(conds[keys[i - 1]]), i == d)) if with_ot else d
        r = conds[keys[rep_i - 1]]
        c = {'name': r['name'], 'ot_ids': list(dict.fromkeys(r['ot_ids'] + [x for i in members
                                                                             for x in conds[keys[i - 1]]['ot_ids']])),
             'evidence': [], 'east_asian': [], 'years': set()}
        for i in members + features.get(d, []):
            m = conds[keys[i - 1]]
            pre = '' if m['name'] == r['name'] else f'{m["name"]}: '
            feat = {'feature': m['name']} if i in features.get(d, []) else {}
            for part in ('evidence', 'east_asian'):
                c[part] += [dict(e, statement=pre + e['statement'], **feat) if pre and not e['statement'].startswith(pre)
                            else dict(e, **feat) for e in m[part]]
            c['years'] |= m.get('years') or set()
            if i != rep_i:
                how = 'attached to' if i in features.get(d, []) else 'merged into'
                logs.append({'pmid': '', 'item': m['name'], 'why': f'{how} "{r["name"]}" by the {step}: {fin[i][2]}',
                             'value': '', 'reader': agree.get(i, '')})
        out[keys[rep_i - 1]] = c
    return out, fin


def consolidate(conds, readers, adj, P, logs):
    """The conditions after the readers' decisions on the list (see consolidation_prompt, decide in llm.py, and
    apply_decisions): (conditions, review rows, table rows, counts)."""
    keys, text = consolidation_list(conds)
    n = len(keys)
    if n < 2:
        return conds, [], [], {'items': n}
    name = lambda i: conds[keys[i - 1]]['name']  # noqa: E731
    show = lambda d: d[0] + (f' into "{name(d[1])}"' if d[1] else '')  # noqa: E731
    final, agree, review, answered = decide(
        readers, adj, consolidation_prompt(P), [(f'The conditions related to {P.short}:\n\n{text}',
                                                 list(range(1, n + 1)))],
        lambda ans, ids, who: ground_decisions(ans, ids, who, n, logs), show)
    if not answered:
        log('consolidation: no reader answered (offline without a cached answer): the conditions are kept as read')
        return conds, [], [], {'items': n, 'not_read': True}
    out, fin = apply_decisions(conds, keys, final, agree, logs, 'consolidation')
    reads = {r['item']: r for r in review}
    names = [r.name for r in readers]
    table = [{'condition': name(i), 'open_targets': 'yes' if conds[keys[i - 1]]['ot_ids'] else 'no',
              'publications': len({e['pmid'] for e in conds[keys[i - 1]]['evidence']
                                   + conds[keys[i - 1]]['east_asian']}),
              'action': fin[i][0], 'into': name(fin[i][1]) if fin[i][1] else '', 'agreement': agree.get(i, 'both'),
              'reason': fin[i][2], 'reader_a': reads.get(i, {}).get(names[0], show(fin[i])),
              'reader_b': reads.get(i, {}).get(names[-1], show(fin[i])) if len(names) > 1 else ''}
             for i in range(1, n + 1)]
    review = [{'pmid': '', 'item': f'consolidation: {name(r["item"])}', **{k: v for k, v in r.items() if k != 'item'}}
              for r in review]
    counts = {'items': n, 'conditions': len(out), **Counter(fin[i][0] for i in range(1, n + 1)),
              'agreement': dict(Counter(agree.get(i, 'both') for i in range(1, n + 1)))}
    return out, review, table, counts


# ── verification: every evidence row and every condition, against the rules ───
DROP_WHY = ('not_about_disease', 'no_relation_stated', 'other_disease_comparison', 'consequence_or_finding',
            'other_form')
VERIFICATION_COLUMNS = ['level', 'condition', 'pmid', 'decision', 'detail', 'agreement', 'reason', 'reader_a',
                        'reader_b']
ONTOLOGY = 'query($id: String!) { disease(efoId: $id) { parents { id name } children { id name } ancestors descendants } }'


def ontology_forms(api, P, max_children):
    """The disease's place in the ontology, for the verification and the mapping: (parent names, sibling names, child
    names, ancestor IDs, descendant IDs). A parent with more than max_children children gives no siblings."""
    d = gql(api, ONTOLOGY, {'id': P.id})['disease'] or {}
    sibs = []
    for par in d.get('parents') or []:
        ch = (gql(api, ONTOLOGY, {'id': par['id']})['disease'] or {}).get('children') or []
        if len(ch) <= max_children:
            sibs += [c['name'] for c in ch if c['id'] != P.id and c['id'] not in (P.d.get('evidence_also') or [])]
    return ([x['name'] for x in d.get('parents') or []], sorted(set(sibs)), [x['name'] for x in d.get('children') or []],
            set(d.get('ancestors') or []), set(d.get('descendants') or []))


def verification_prompt(P, siblings):
    s = P.short
    cats = '\n'.join(f'- "{c["name"]}": {c["description"]}' for c in P.categories)
    forms = f' ({", ".join(siblings)})' if siblings else ''
    return f"""You verify the evidence that relates conditions to {P.name} ({s}). Each condition comes with its id, its
name, whether Open Targets holds an entry for it, its current category, and its evidence rows, each with an id, a PMID,
a direction, the design and population, a statement and the verbatim quote of the PubMed abstract. Judge every row and
every condition by these rules, and by the text given only.

A row counts when its quote, read in its record, states a relation between the condition and {s} itself. Drop a row,
with the code of the reason, when:
- "not_about_disease": its result is about another disease or outcome (a related but different disease, a suspected
  or composite outcome), or about a broader group in which {s} patients are not most;
- "no_relation_stated": the quote is background, a study aim or rationale, a design statement, or says only that two
  conditions share features, resemble each other or must be told apart;
- "other_disease_comparison": the quote compares {s} patients with patients of another disease only, not with people
  without {s};
- "consequence_or_finding": the condition is a consequence or complication of {s} or of its treatment, or a finding,
  sign or measurement among {s} patients (heart failure, hypoxemia, anxiety or depression reported among {s}
  patients are consequences, not related conditions);
- "other_form": the condition defines, by definition, another form of {s}'s parent class{forms}: an exposure that
  defines a drug- or toxin-induced form, a disease that defines a form associated with another disease, a family
  history that defines a heritable form. Such a row relates that form, not {s}, even when an older record calls its
  patients by a former, broader name of {s} (a name that once covered all these forms).
Give every row that counts its category:
{cats}
Then decide each condition: "keep"; "drop" (not a condition, or none of its rows counts); "merge" with "into" the id
of another listed condition that is the same condition under another name; or "feature" with "into" the id of another
listed condition whose patients it describes (a characteristic of that condition: its severity, type, location,
provocation, recurrence, a finding at its diagnosis, its treatment, or a comparison between groups of its patients,
such as {s} after one kind of that condition against {s} after another). A narrower disease of its own is never merged
into a broader one. When unsure, keep. Answer with JSON only:
{{"conditions": [{{"id": <id>, "action": "keep" | "drop" | "merge" | "feature", "into": <id, or null>,
"reason": "<at most 20 words>", "rows": [{{"row": "<row id>", "keep": true | false, "why": "<a code, for a dropped
row>", "category": "<a category, for a kept row>"}}]}}]}}, one entry per condition of this part and one row entry
per row."""


def verification_batches(conds, keys, max_rows=40):
    """The parts the readers verify: every condition with all its rows, about max_rows rows a part (a condition is
    never split), each part with the numbered list of every condition (for "into")."""
    names = '\n'.join(f'[{i}] {conds[k]["name"]}' for i, k in enumerate(keys, 1))
    parts, cur, rows = [], [], 0
    for i, k in enumerate(keys, 1):
        n = len(conds[k]['evidence']) + len(conds[k]['east_asian'])
        if cur and rows + n > max_rows:
            parts.append(cur)
            cur, rows = [], 0
        cur.append(i)
        rows += n
    if cur:
        parts.append(cur)
    out = []
    for part in parts:
        blocks, ids = [], []
        for i in part:
            c = conds[keys[i - 1]]
            blocks.append(f'[{i}] {c["name"]} | Open Targets: {"yes" if c["ot_ids"] else "no"} | category: '
                          f'{votes(c).most_common(1)[0][0]}')
            ids.append(i)
            for k, e in enumerate(c['evidence'] + c['east_asian'], 1):
                blocks.append(f'  ({i}.{k}) PMID {e["pmid"]} | {e["direction"]} | {e["design"]}; {e["population"]}\n'
                              f'      statement: {e["statement"]}\n      quote: "{e["quote"]}"')
                ids.append(f'{i}.{k}')
        out.append((f'Every condition of the list (ids for "into"):\n{names}\n\n---\n\nThe conditions to verify in '
                    f'this part, with their evidence:\n\n' + '\n'.join(blocks), ids))
    return out


def ground_verification(ans, ids, who, n, cats, logs):
    """{id: decision}: a condition (int id) -> (action, into, reason) as in ground_decisions; a row ('i.k') -> (keep,
    category or '', reason code). An invalid decision is logged and its item left as it was."""
    out = ground_decisions(ans, {i for i in ids if isinstance(i, int)}, who, n, logs, key='conditions')
    want = {i for i in ids if isinstance(i, str)}
    for c in (ans.get('conditions') if isinstance(ans, dict) else None) or []:
        for r in (c.get('rows') if isinstance(c, dict) else None) or []:
            if not isinstance(r, dict) or str(r.get('row')) not in want:
                continue
            rid, keep = str(r['row']), r.get('keep')
            if keep is True and r.get('category') in cats:
                out[rid] = (True, r['category'], '')
            elif keep is False and r.get('why') in DROP_WHY:
                out[rid] = (False, '', r['why'])
            else:
                logs.append({'pmid': '', 'item': f'verification row {rid}', 'why': 'decision not allowed: left as it '
                             'was', 'value': json.dumps(r)[:200], 'reader': who})
    return out


def verify(conds, readers, adj, P, logs, siblings):
    """The conditions after the verification of every row and every condition (verification_prompt): dropped rows
    leave, kept rows take the category given, then each condition's decision applies (apply_decisions) and a
    condition left without a row leaves. -> (conditions, review rows, table rows, counts)"""
    keys = sorted(conds, key=lambda k: (conds[k]['name'].lower(), k))
    n = len(keys)
    if not n:
        return conds, [], [], {'conditions': 0}
    cats = [c['name'] for c in P.categories]
    name = lambda i: conds[keys[i - 1]]['name']  # noqa: E731

    def show(d):
        if isinstance(d[0], bool):
            return f'keep as {d[1]}' if d[0] else f'drop ({d[2]})'
        return d[0] + (f' into "{name(d[1])}"' if d[1] else '')
    final, agree, review, answered = decide(
        readers, adj, verification_prompt(P, siblings), verification_batches(conds, keys),
        lambda ans, ids, who: ground_verification(ans, ids, who, n, cats, logs), show)
    if not answered:
        log('verification: no reader answered (offline without a cached answer): the conditions are kept as read')
        return conds, [], [], {'conditions': n, 'not_read': True}
    names = [r.name for r in readers]
    reads = {r['item']: r for r in review}
    table, dropped_rows, recat, work = [], Counter(), 0, {}
    for i, k in enumerate(keys, 1):
        c = conds[k]
        keep, k_row = {'evidence': [], 'east_asian': []}, 0
        for part in ('evidence', 'east_asian'):
            for e in c[part]:
                k_row += 1
                rid = f'{i}.{k_row}'
                d = final.get(rid)
                if d and not d[0]:
                    dropped_rows[d[2]] += 1
                    logs.append({'pmid': e['pmid'], 'item': c['name'], 'why': f'evidence row dropped by the '
                                 f'verification: {d[2]}', 'value': e['quote'][:200], 'reader': agree.get(rid, '')})
                else:
                    if d and d[1] != e.get('category'):
                        recat += 1
                    keep[part].append(dict(e, category=d[1]) if d else e)
                table.append({'level': 'row', 'condition': c['name'], 'pmid': e['pmid'],
                              'decision': 'keep' if not d or d[0] else 'drop',
                              'detail': d[1] if d and d[0] else d[2] if d else e.get('category', ''),
                              'agreement': agree.get(rid, 'not decided'), 'reason': '',
                              'reader_a': reads.get(rid, {}).get(names[0], show(d) if d else ''),
                              'reader_b': (reads.get(rid, {}).get(names[-1], show(d) if d else '')
                                           if len(names) > 1 else '')})
        work[k] = dict(c, **keep)
    # a condition left without rows leaves; the others take their decision
    for i, k in enumerate(keys, 1):
        if not (work[k]['evidence'] + work[k]['east_asian']) and final.get(i, ('keep',))[0] != 'drop':
            final[i] = ('drop', None, 'none of its evidence rows counts')
    out, fin = apply_decisions(work, keys, {i: final[i] for i in range(1, n + 1) if i in final}, agree, logs,
                               'verification')
    for i in range(1, n + 1):
        table.append({'level': 'condition', 'condition': name(i), 'pmid': '', 'decision': fin[i][0],
                      'detail': name(fin[i][1]) if fin[i][1] else '', 'agreement': agree.get(i, 'not decided'),
                      'reason': fin[i][2], 'reader_a': reads.get(i, {}).get(names[0], show(fin[i])),
                      'reader_b': reads.get(i, {}).get(names[-1], show(fin[i])) if len(names) > 1 else ''})
    review_rows = []
    for r in review:
        i = r['item']
        rest = {k2: v for k2, v in r.items() if k2 != 'item'}
        if isinstance(i, str):
            ci, kk = (int(x) for x in i.split('.'))
            e = (conds[keys[ci - 1]]['evidence'] + conds[keys[ci - 1]]['east_asian'])[kk - 1]
            review_rows.append({'pmid': e['pmid'], 'item': f'verification: {name(ci)} (evidence row)', **rest})
        else:
            review_rows.append({'pmid': '', 'item': f'verification: {name(i)}', **rest})
    counts = {'conditions_in': n, 'conditions': len(out), 'rows_dropped': dict(dropped_rows),
              'rows_recategorized': recat, **Counter(fin[i][0] for i in range(1, n + 1)),
              'agreement': dict(Counter(agree.values()))}
    return out, review_rows, table, counts


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--records', required=True, help='pubmed.py --mode search output of the relation search')
    ap.add_argument('--profile', required=True)
    ap.add_argument('--llm', required=True)
    ap.add_argument('--offline', default='false')
    ap.add_argument('--batch', type=int, default=5)
    ap.add_argument('--limit', type=int, default=0, help='smoke run: the first N records only')
    ap.add_argument('--api', required=True, help='Open Targets GraphQL (condition names)')
    ap.add_argument('--max-children', type=int, default=30,
                    help='a parent with more children gives no siblings (the ontology neighbors)')
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
    # the disease's place in the ontology: its siblings (for the verification), its ancestors (never a related
    # condition: a name mapped to one keeps no entry) and its descendants (its own subtypes, left out)
    _parents, siblings, _children, ancestors, descendants = ontology_forms(a.api, P, a.max_children)
    for t, (i, nm) in sorted(mapped.items()):
        if i in ancestors:
            logs.append({'pmid': '', 'item': t, 'why': f'maps to {nm} ({i}), an ancestor of {P.short} in the '
                         'ontology: not mapped', 'value': '', 'reader': ''})
    mapped = {t: h for t, h in mapped.items() if h[0] not in ancestors}
    own_ids |= descendants
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
            logs.append({'pmid': x['pmid'], 'item': x['condition'], 'why': f'{P.short} itself or one of its subtypes',
                         'value': '', 'reader': ''})
            continue
        key = hit[0] if hit else custom_key(name)     # resolve: custom:<key>
        c = conds.setdefault(key, {'name': name[:1].upper() + name[1:], 'ot_ids': [hit[0]] if hit else [],
                                   'evidence': [], 'east_asian': []})
        r = by_rec[x['pmid']]
        e = {'pmid': x['pmid'], 'direction': x['direction'], 'design': x['design'], 'population': x['population'],
             'statement': x['statement'], 'quote': x['quote'], 'agreement': x['agreement'], 'category': x['category']}
        # an East Asian population: the final reading's population names an East Asian country, people or place,
        # and so does the abstract (whole words; not necessarily the same: 'China' and 'Peking'), and neither names
        # a country outside East Asia (a multinational study is not East Asian)
        text = r['title'] + ' ' + r['abstract']
        words = EAS_WORDS + EAST_ASIA_PLACES
        ea = (any(in_words(w, x['population']) for w in words) and any(in_words(w, text) for w in words)
              and not (countries_named(x['population'], C) | countries_named(text, C)) - EAS_SET)
        (c['east_asian'] if ea else c['evidence']).append(e)
        c.setdefault('years', set()).add(r['year'])
    n_read = len(conds)
    conds, review_c, table_c, counts_c = consolidate(conds, readers, adj, P, logs)
    review += review_c
    n_consolidated = len(conds)
    conds, review_v, table_v, counts_v = verify(conds, readers, adj, P, logs, siblings)
    review += review_v
    seed = {'version': load_json(a.records).get('max_date', ''),
            'note': 'Extracted: read from PubMed by two independent language-model readers, adjudicated where they '
                    'disagreed; every quote is checked against its abstract at each build.',
            'source': 'extracted', 'conditions': [], 'facts': []}
    for key, c in sorted(conds.items(), key=lambda kv: kv[1]['name'].lower()):
        n_pub = len({e['pmid'] for e in c['evidence'] + c['east_asian']})
        seed['conditions'].append({
            'key': key.lower(), 'name': c['name'], 'category': votes(c).most_common(1)[0][0],
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
    write_tsv(pd.DataFrame(table_c, columns=CONSOLIDATION_COLUMNS), out / 'relations_consolidation.tsv')
    write_tsv(pd.DataFrame(table_v, columns=VERIFICATION_COLUMNS), out / 'relations_verification.tsv')
    ag = Counter(x['agreement'] for x in rows)
    not_read = sum(all(x is None for x in readings[r['pmid']]) for r in recs)
    unanswered = sum(1 for x in logs if x.get('item') == 'record' and x.get('value') == 'no answer')
    stats = {'records': len(recs), 'read': len(recs) - not_read, 'not_read': not_read,
             'relations': len(rows), 'conditions_read': n_read, 'consolidation': counts_c,
             'conditions_consolidated': n_consolidated, 'verification': counts_v, 'siblings': siblings,
             'conditions': len(seed['conditions']),
             'mapped_to_open_targets': sum(1 for c in seed['conditions'] if c['ot_ids']),
             'one_publication': sum(1 for c in seed['conditions'] if c['support'] == 'one publication'),
             'agreement': dict(ag), 'disputed_items': len(review), 'dropped_items': len(logs) - unanswered,
             'open_targets_release': rel['data_version'],
             'readers': [r.stats() for r in readers], 'adjudicator': adj.stats() if adj else None}
    (out / 'relations.json').write_text(json.dumps(stats, indent=1))
    log(f'{len(recs)} records: {len(rows)} relations ({dict(ag)}), {n_read} conditions read, consolidated to '
        f'{n_consolidated} ({ {k: v for k, v in counts_c.items() if k in ACTIONS} }), verified to '
        f'{len(seed["conditions"])} ({ {k: v for k, v in counts_v.items() if k in ACTIONS} }; rows dropped '
        f'{counts_v.get("rows_dropped", {})}, recategorized {counts_v.get("rows_recategorized", 0)}); '
        f'{stats["mapped_to_open_targets"]} in Open Targets, {stats["one_publication"]} on one publication; '
        f'{len(review)} disputed, {len(logs) - unanswered} dropped' + (f', {not_read} not read' if not_read else ''))


if __name__ == '__main__':
    main()
