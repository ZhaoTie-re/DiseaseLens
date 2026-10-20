#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The curation step between what was fetched and extracted and what
#           the knowledge graph holds. Deterministic: every rule is in the
#           curation rules of the profile (config/curation.yaml merged with the
#           disease's own entries), and every rule applied is logged row by row.
#             publications   the status of every publication the graph uses,
#                            from its PubMed record: excluded (retracted, or
#                            conference material), flagged (preprint, erratum,
#                            expression of concern) or kept; the country of
#                            the first author's institution
#             studies        the genetic studies of the disease without the excluded
#                            ones; place (country of the first author's
#                            institution, else the country the abstract states),
#                            Japanese / East Asian flags from the place, one
#                            design class each
#             findings       without the excluded studies; one name and an rsID
#                            per variant; what each finding is about (the
#                            disease itself, a phenotype or outcome within it,
#                            another trait, a functional experiment, a
#                            background statement); gene sets marked
#             screening      every PubMed record; an excluded publication is
#                            marked as left out by its publication status
#             literature     the text-mined literature of the disease (Europe PMC
#                            via Open Targets) without excluded publications and
#                            without words mapped to the wrong gene (their
#                            sentences go with them): the words the disease lists,
#                            and any word that names the disease or an ontology
#                            neighbor (its name, an exact synonym)
#             relations      the clinical relations without excluded publications:
#                            a curated list may not cite one (the build stops); an
#                            extracted relation that does is left out
#             curated genes  the curated gene evidence of the disease itself
#                            without genes whose symbol is a word naming the disease
#                            or an ontology neighbor (e.g. PAH, phenylalanine
#                            hydroxylase, for pulmonary arterial hypertension)
# Inputs  : publications.json.gz (publications.py), 02.literature, 00.raw,
#           relations.tsv, the profile
# Outputs : publications.tsv, studies.tsv, findings.tsv, methods.tsv,
#           screening.tsv, literature.tsv, textmining_review.tsv,
#           relations.tsv, model_review.tsv (the disputed items of the readers),
#           curation_log.tsv (step, item, rule, before, after, basis),
#           curation.json (counts)
# ---------------------------------------------------------------------------
import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import (load_json, load_profile, log, mining_links, read_tsv, sentences, us_english,  # noqa: E402
                       write_tsv)

LOG_COLS = ['step', 'item', 'rule', 'before', 'after', 'basis']


# ── publications ─────────────────────────────────────────────────────────────
def supplement(rec):
    """A journal supplement: 'Suppl' in the volume or issue, or in the DOI."""
    return bool(re.search(r'suppl', f'{rec["volume"]} {rec["issue"]} {rec["doi"]}', re.I))


def pub_status(rec, rules):
    """(status, rule key, label): excluded, flagged or kept, by the first rule that applies."""
    pts, com = set(rec['pub_types']), set(rec.get('comments') or [])
    for r in rules['exclude']:
        if pts & set(r.get('pub_types') or []) or com & set(r.get('comments') or []):
            return 'excluded', r['key'], r['label']
        if r.get('supplement_without_abstract') and not rec['has_abstract'] and supplement(rec):
            return 'excluded', r['key'], f'{r["label"]} (journal supplement, no abstract)'
    for r in rules['flag']:
        if pts & set(r.get('pub_types') or []):
            return 'flagged', r['key'], r['label']
    return 'kept', '', ''


def country_of(aff, C):
    """The country of an affiliation: the last country it names (aliases mapped), else a US state, else a city."""
    s = aff or ''
    for w in C['ignore']:
        s = s.replace(w, ' ')
    names = {n: n for n in C['names']} | C['aliases']
    hits = [(m.start(), v) for k, v in names.items()
            for m in re.finditer(rf'(?<![A-Za-z]){re.escape(k)}(?![A-Za-z])', s)]
    if hits:
        return max(hits)[1]
    abbrev = [x for x in C['us_states'] if len(x) == 2]
    full = [x for x in C['us_states'] if len(x) > 2]
    if re.search(r',\s*(' + '|'.join(abbrev) + r')\b', s) or any(re.search(rf'\b{x}\b', s) for x in full):
        return 'United States'
    for city, country in C['cities'].items():
        if re.search(rf'\b{re.escape(city)}\b', s):
            return country
    return ''


def disease_words(rec):
    """{word: what it names}: the name and exact synonyms of the disease, its parents and its children, and of
    the other entries whose evidence it takes (evidence_also) (Open Targets), as words a text-mining dictionary can
    mistake for a gene symbol."""
    out = {}
    for d in [rec] + (rec.get('parents') or []) + (rec.get('children') or []) + (rec.get('evidence_also') or []):
        terms = [d['name']] + [t for s in d.get('synonyms') or [] if s.get('relation') == 'hasExactSynonym'
                               for t in s['terms']]
        for t in terms:
            out.setdefault(t, f'the disease "{d["name"]}"')
    return out


# ── studies and findings ─────────────────────────────────────────────────────
def design_class(s, classes):
    designs = set(str(s.design).split('; '))
    for c in classes:
        if c.get('study_type') and s.study_type == c['study_type']:
            return c
        if set(c.get('designs') or []) & designs:
            return c
    return classes[-1]


def variant_name(f, defs):
    """(name, rsID, how): the curated name of the finding's variant, matched in its variant text, then its gene
    text, then its source sentence; else the variant as written, with an rsID it states."""
    own = [v for v in defs if v['gene'] == f.symbol]
    for field, how in ((f.variant, 'variant'), (f.gene, 'gene text'), (f.quote, 'source sentence')):
        for v in own:
            if field and any(re.search(p, field, re.I) for p in v['patterns']):
                return v['name'], v['rsid'] or '', how
    rs = re.findall(r'(?<![A-Za-z])rs\d{3,}\b', f'{f.variant} {f.gene}')
    rs = rs[0] if rs else ''
    bare = re.sub(rf'^{re.escape(f.symbol)}\s*', '', f.variant).strip() if f.symbol else f.variant
    return (rs if rs and bare in ('', rs) else f.variant), rs, ''


def matches(f, it):
    gene_ok = not it.get('gene') or f.symbol == it['gene'] or it['gene'] in f.gene
    return str(f.pmid) == str(it['pmid']) and gene_ok and it['match'].lower() in f.quote.lower()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--publications', required=True, help='publications.py output')
    ap.add_argument('--literature', required=True, help='02.literature (extract_genetics.py outputs)')
    ap.add_argument('--raw', required=True, help='00.raw (text-mined evidence)')
    ap.add_argument('--relations', required=True, help='01.conditions/relations.tsv')
    ap.add_argument('--profile', required=True, help='profile.json of the disease (its curation rules)')
    ap.add_argument('--sentences', type=int, default=3, help='sentences kept per text-mined publication')
    ap.add_argument('--out-dir', default='03.curated')
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    PR = load_profile(a.profile)
    cur, ABOUT, short = PR.curation, PR.about, PR.short
    P, X = load_json(a.publications), Path(a.literature)
    C = cur['countries']
    logs = []

    def note(step, item, rule, before, after, basis=''):
        logs.append({'step': step, 'item': item, 'rule': rule, 'before': before, 'after': after, 'basis': basis})

    # ── publications ─────────────────────────────────────────────────────────
    rules = cur['publications']
    basis = {r['key']: r.get('basis', '') for r in rules['exclude'] + rules['flag']}
    prow = []
    for rec in P['records']:
        status, key, label = pub_status(rec, rules)
        prow.append({'pmid': rec['pmid'], 'status': status, 'rule': key, 'reason': label,
                     'roles': '; '.join(P['roles'].get(rec['pmid'], [])), 'first_author': rec['first_author'],
                     'year': rec['year'], 'title': rec['title'], 'journal': rec['journal'], 'volume': rec['volume'],
                     'issue': rec['issue'], 'pages': rec['pages'], 'pub_types': '; '.join(rec['pub_types']),
                     'has_abstract': rec['has_abstract'], 'comments': '; '.join(rec['comments']),
                     'doi': rec['doi'], 'pmcid': rec['pmcid'], 'affiliation': rec['affiliation'],
                     'affiliation_source': rec['affiliation_source'],
                     'institution_country': country_of(rec['affiliation'], C)})
    prefix = {p: r for r in rules['flag'] for p in r.get('id_prefixes') or []}
    for pid, roles in P['other_ids'].items():
        r = next((prefix[p] for p in prefix if pid.startswith(p)), None)
        prow.append({'pmid': pid, 'status': 'flagged' if r else 'not checked', 'rule': r['key'] if r else '',
                     'reason': r['label'] if r else 'not a PubMed identifier', 'roles': '; '.join(roles)})
    pubs = pd.DataFrame(prow).fillna('')
    for c in ('has_abstract',):
        pubs[c] = pubs[c].astype(str)
    excluded = dict(zip(pubs[pubs.status == 'excluded'].pmid, pubs[pubs.status == 'excluded'].reason))
    flagged = dict(zip(pubs[pubs.status == 'flagged'].pmid, pubs[pubs.status == 'flagged'].reason))
    for r in pubs[pubs.status == 'excluded'].itertuples():
        note('publications', f'PMID {r.pmid} ({r.first_author} {r.year})', r.rule, f'used in: {r.roles}',
             f'left out: {r.reason}', basis[r.rule])
    rel = read_tsv(a.relations)
    bad = sorted(set(rel.pmid) & set(excluded))
    extracted = PR.relations.get('mode') == 'extracted'
    assert extracted or not bad, f'the curated seed cites excluded publications: {bad}; replace them in ' \
                                 'related_conditions.yaml'
    for r in rel[rel.pmid.isin(excluded)].itertuples():
        note('relations', f'PMID {r.pmid}', 'excluded publication', f'{r.evidence_type}: {r.detail}',
             f'left out: {excluded[r.pmid]}', basis[pubs.set_index('pmid').rule[r.pmid]])
    rel = rel[~rel.pmid.isin(excluded)]

    # ── genetic studies of the disease ──────────────────────────────────────
    st = read_tsv(X / 'studies.tsv', ['pmid', 'year', 'first_author', 'title', 'study_type', 'design', 'country',
                                      'japanese', 'east_asian'])
    fi = read_tsv(X / 'findings.tsv', ['pmid', 'gene', 'symbol', 'variant', 'result', 'quote', 'ensembl_id'])
    me = read_tsv(X / 'methods.tsv', ['pmid', 'field', 'value', 'quote', 'source'])
    sc = read_tsv(X / 'screening.tsv', ['pmid', 'year', 'title', 'include', 'study_type', 'stage', 'reason'])
    drop = st.pmid.isin(excluded.keys())
    for r in st[drop].itertuples():
        n_f = int((fi.pmid == r.pmid).sum())
        note('studies', f'PMID {r.pmid} ({r.first_author} {r.year})', 'excluded publication',
             f'included genetic study ({r.design}; {n_f} findings)', f'left out: {excluded[r.pmid]}',
             basis[pubs.set_index('pmid').rule[r.pmid]])
    st, fi, me = st[~drop].copy(), fi[~fi.pmid.isin(excluded)].copy(), me[~me.pmid.isin(excluded)].copy()
    st['design'] = st.design.map(us_english)       # the extraction's design labels, in American spelling
    pi = pubs.set_index('pmid')
    eas = set(C['east_asian'])
    classes = cur['design_classes']
    rows = []
    for s in st.itertuples():
        inst = pi.institution_country.get(s.pmid, '') if s.pmid in pi.index else ''
        stated = '' if s.country in ('', 'not stated') else C['aliases'].get(s.country, s.country)
        place, src = (inst, "first author's institution") if inst else ((stated, 'stated in the abstract') if stated
                                                                        else ('', ''))
        if inst and stated and stated != inst:
            note('studies', f'PMID {s.pmid} ({s.first_author} {s.year})', 'place',
                 f'population stated: {stated}', f'place: {inst} (first author\'s institution)',
                 'The place of a study is the country of its first author\'s institution; the population the '
                 'abstract states is kept beside it.')
        c = design_class(s, classes)
        jp = place == 'Japan' or s.japanese == 'True'
        ea = jp or place in eas or s.east_asian == 'True'
        if jp != (s.japanese == 'True'):
            note('studies', f'PMID {s.pmid} ({s.first_author} {s.year})', 'Japanese flag', s.japanese, str(jp),
                 'Japanese: done in Japan (first author\'s institution) or a Japanese population stated.')
        rows.append({'place': place, 'place_source': src, 'institution_country': inst,
                     'affiliation': pi.affiliation.get(s.pmid, '') if s.pmid in pi.index else '',
                     'japanese': jp, 'east_asian': ea, 'design_class': c['key'], 'design_class_label': c['label'],
                     'publication_flag': flagged.get(s.pmid, '')})
    st = pd.concat([st.reset_index(drop=True), pd.DataFrame(rows, columns=[
        'place', 'place_source', 'institution_country', 'affiliation', 'japanese', 'east_asian', 'design_class',
        'design_class_label', 'publication_flag'])], axis=1)
    st = st.loc[:, ~st.columns.duplicated(keep='last')]
    for c in st.design_class_label.unique():
        note('studies', c, 'design class', '', f'{int((st.design_class_label == c).sum())} studies',
             'The first class of config/curation.yaml (design_classes) whose rule the study meets.')

    # ── findings ─────────────────────────────────────────────────────────────
    F = cur['findings']
    kinds = [(k, F[k]) for k in ('within_disease', 'other_trait', 'functional', 'background')]
    used = set()
    about, about_note, sets, members = [], [], [], []
    for f in fi.itertuples():
        lab, nt = ABOUT, ''
        for k, spec in kinds:
            for j, it in enumerate(spec['items']):
                if matches(f, it):
                    lab, nt = spec['label'], it.get('note', '')
                    used.add((k, j))
                    note('findings', f'PMID {f.pmid} {f.symbol or f.gene or "(no gene)"}', 'about',
                         f'{f.result}: about {short}', f'{f.result}: about {lab} ({nt})', spec['basis'])
        gs = ''
        for j, it in enumerate(F['gene_sets']):
            if str(f.pmid) == str(it['pmid']) and it['match'].lower() in (f.quote + ' ' + f.gene).lower():
                gs = it['label']
                used.add(('gene_sets', j))
        about.append(lab)
        about_note.append(nt)
        sets.append(gs)
        members.append('; '.join(x for x in re.findall(r'\b[A-Z][A-Z0-9]{2,}\b', f.gene)) if gs else '')
    unused = [(k, it) for k, spec in kinds for j, it in enumerate(spec['items']) if (k, j) not in used] + \
        [('gene_sets', it) for j, it in enumerate(F['gene_sets']) if ('gene_sets', j) not in used]
    assert not unused, f'curation items that match no finding (stale config): {unused}'
    fi['association_with'], fi['association_note'], fi['gene_set'] = about, about_note, sets
    fi['gene_set_members'] = members
    names = [variant_name(f, cur['variants']) for f in fi.itertuples()]
    fi['variant_name'] = [n for n, _r, _h in names]
    fi['rsid'] = [r for _n, r, _h in names]
    for f, (n, r, how) in zip(fi.itertuples(), names):
        if how and n != f.variant:
            note('findings', f'PMID {f.pmid} {f.symbol}', 'variant name', f.variant or '(none written)',
                 n + (f' ({r})' if r else ''), f'Matched in the {how}; config/curation.yaml (variants).')
    words = F['review_words']
    open_ = fi[(fi.association_with == ABOUT) & fi.result.isin(['associated', 'enriched', 'other'])]
    for f in open_.itertuples():
        hit = [w for w in words if w.lower() in f.quote.lower()]
        if hit:
            note('findings', f'PMID {f.pmid} {f.symbol or f.gene}', 'to review', f'about {short} ({f.result})',
                 f'words of a finding within {short}: ' + ', '.join(hit),
                 f'Read the source sentence; list the finding in the curation file of {short} if it is not about '
                 f'{short} itself.')
    fi['counts_for_ladder'] = (fi.association_with == ABOUT) & fi.symbol.ne('')

    # ── screening ────────────────────────────────────────────────────────────
    sc['publication_status'] = sc.pmid.map(lambda p: excluded.get(p) or flagged.get(p) or '')
    hit = sc.pmid.isin(excluded)
    sc.loc[hit, 'include'] = 'False'
    sc.loc[hit, 'stage'] = 'publication status'
    sc.loc[hit, 'reason'] = 'left out: ' + sc.loc[hit, 'pmid'].map(excluded)
    assert set(st.pmid) == set(sc[sc.include == 'True'].pmid), 'screening and studies disagree'

    # ── text-mined literature ────────────────────────────────────────────────
    wrong = {}
    for it in cur['text_mining'] + (cur.get('shared_text_mining') or []):
        for w in it['words']:
            wrong.setdefault(it['gene'], {})[w] = it['meaning']
    listed = {(it['gene'], w) for it in cur['text_mining'] for w in it['words']}   # stale-checked
    dwords = disease_words(load_json(Path(a.raw) / 'disease.json'))
    ev = load_json(Path(a.raw) / 'evidence.json.gz')
    men = load_json(Path(a.raw) / 'mendelian.json.gz') if (Path(a.raw) / 'mendelian.json.gz').exists() else \
        {'mendelian_targets': {}, 'mendelian_evidence': []}
    # a gene with curated evidence for the disease keeps a word that names the disease: the word is then the old
    # name of the gene's locus (PPH1, the BMPR2 locus of pulmonary hypertension), not a mis-mapping; a gene whose
    # symbol is itself such a word is not one of them (below)
    own_curated = {t_['symbol'] for ts in men['mendelian_targets'].values() for t_ in ts
                   if t_['symbol'] not in dwords}
    lit, review = [], []
    for e in mining_links(ev):                                   # one link per gene and publication
        sym = e['target']['approvedSymbol']
        _t, words_all, _n, _b = sentences(e, a.sentences)
        for w in words_all:
            if w in dwords and w not in wrong.get(sym, {}) and sym not in own_curated:
                wrong.setdefault(sym, {})[w] = dwords[w]
        title, words, n, best = sentences(e, a.sentences, set(wrong.get(sym, {})))
        gone = [w for w in words_all if w not in words]
        for pmid in e['literature'] or ['']:
            if pmid in excluded:
                decision, why = 'dropped', f'publication left out: {excluded[pmid]}'
            elif words_all and not words:
                decision, why = 'dropped', 'mapped to the wrong gene: ' + '; '.join(
                    f'{w} = {wrong[sym][w]}' for w in gone)
            else:
                decision = 'kept'
                why = ('words removed: ' + '; '.join(f'{w} = {wrong[sym][w]}' for w in gone)) if gone else ''
            review.append({'ensembl_id': e['target']['id'], 'symbol': sym, 'pmid': pmid,
                           'matched_as': ' | '.join(words_all), 'kept_words': ' | '.join(words),
                           'decision': decision, 'reason': why, 'publication_flag': flagged.get(pmid, '')})
            if decision == 'kept':
                lit.append({'ensembl_id': e['target']['id'], 'symbol': sym, 'pmid': pmid,
                            'year': e.get('publicationYear'), 'title': title, 'matched_as': ' | '.join(words),
                            'n_sentences': n, 'sentences': ' || '.join(best), 'evidence_score': e['score'],
                            'publication_flag': flagged.get(pmid, '')})
    # explicit columns: a disease can have no text-mined publication at all (PVOD/PCH in Open Targets 26.09)
    rv = pd.DataFrame(review, columns=['ensembl_id', 'symbol', 'pmid', 'matched_as', 'kept_words', 'decision',
                                       'reason', 'publication_flag'])
    used_words = {(r.symbol, w) for r in rv.itertuples() for w in r.matched_as.split(' | ')}
    stale = sorted((g, w) for g, w in listed if (g, w) not in used_words)
    assert not stale, f'text-mining corrections that match no link (stale config): {stale}'
    for (sym, why), g in rv[rv.decision == 'dropped'].groupby(['symbol', 'reason']):
        note('text mining', sym, 'link dropped', f'{len(g)} publication{"s" if len(g) > 1 else ""}', why,
             'The curation rules (publications; text mining: the words listed for the disease, and words that name '
             'the disease or an ontology neighbor).')
    for (sym, why), g in rv[(rv.decision == 'kept') & rv.reason.ne('')].groupby(['symbol', 'reason']):
        note('text mining', sym, 'words removed', f'{len(g)} publication{"s" if len(g) > 1 else ""}', why,
             'The curation rules (text mining).')
    lit = pd.DataFrame(lit, columns=['ensembl_id', 'symbol', 'pmid', 'year', 'title', 'matched_as', 'n_sentences',
                                     'sentences', 'evidence_score', 'publication_flag'])
    lit = lit.sort_values(['symbol', 'year', 'pmid'], ascending=[True, False, True])
    assert not lit.duplicated(['ensembl_id', 'pmid']).any()

    # ── curated gene evidence of the disease itself ─────────────────────────
    core_dropped = sorted({t_['symbol'] for ts in men['mendelian_targets'].values() for t_ in ts
                           if t_['symbol'] in dwords})
    for s in core_dropped:
        note('curated genes', s, 'gene symbol names a disease', f'curated evidence for {short}',
             f'left out: "{s}" is {dwords[s]}',
             'A gene whose symbol is the abbreviation of the disease or of an ontology neighbor is mapped to it by '
             'the abbreviation, not by evidence.')

    # ── the readers' disputed items: the genetic studies, and the clinical relations when read from PubMed ─────
    mr = read_tsv(X / 'review.tsv') if (X / 'review.tsv').exists() else pd.DataFrame()
    C = Path(a.relations).parent
    if (C / 'relations_review.tsv').exists():
        rr = read_tsv(C / 'relations_review.tsv').assign(**{'pass': 'relations'})
        mr = pd.concat([mr, rr], ignore_index=True).fillna('') if len(mr) else rr
    reading = json.loads((C / 'relations_reading.json').read_text()) if (C / 'relations_reading.json').exists() else None
    # an excluded publication is counted, never listed: its disputed items leave the review table
    mr = mr[~mr.pmid.isin(excluded)] if len(mr) else mr

    # ── write ────────────────────────────────────────────────────────────────
    write_tsv(pubs, out / 'publications.tsv')
    write_tsv(st, out / 'studies.tsv')
    write_tsv(fi, out / 'findings.tsv')
    write_tsv(me, out / 'methods.tsv')
    write_tsv(sc, out / 'screening.tsv')
    write_tsv(lit, out / 'literature.tsv')
    write_tsv(rv, out / 'textmining_review.tsv')
    write_tsv(rel, out / 'relations.tsv')
    write_tsv(mr, out / 'model_review.tsv', columns=['pmid', 'pass', 'item', 'final', 'decision'])
    L = pd.DataFrame(logs, columns=LOG_COLS)
    write_tsv(L, out / 'curation_log.tsv')
    stats = {'checked_at': P['checked_at'], 'publications': len(pubs),
             'by_status': pubs.status.value_counts().to_dict(),
             # counted only: the excluded publications are listed in publications.tsv (status) and the curation
             # log of this folder, the record of the curation, and in no product
             'excluded': int((pubs.status == 'excluded').sum()),
             'studies': len(st), 'studies_excluded': int(drop.sum()), 'findings': len(fi),
             'findings_by_subject': fi.association_with.value_counts().to_dict(),
             'variants_named': int(sum(1 for _n, _r, h in names if h)),
             'design_classes': st.design_class_label.value_counts().to_dict(),
             'places': st.place.replace('', 'not known').value_counts().to_dict(),
             'japanese': int(st.japanese.sum()), 'east_asian': int(st.east_asian.sum()),
             'textmining_links': len(rv), 'textmining_kept': int((rv.decision == 'kept').sum()),
             'textmining_dropped': int((rv.decision == 'dropped').sum()),
             'textmining_words_removed': int(((rv.decision == 'kept') & rv.reason.ne('')).sum()),
             'to_review': int((L.rule == 'to review').sum()), 'log_rows': len(L),
             'core_curated_dropped': core_dropped, 'relations_excluded': int(len(bad)),
             'model_review': len(mr), 'relations_reading': reading}
    (out / 'curation.json').write_text(json.dumps(stats, indent=1, ensure_ascii=False))
    log(f'publications {len(pubs)} ({stats["by_status"]}); studies {len(st)} (-{stats["studies_excluded"]}); '
        f'findings {len(fi)} {stats["findings_by_subject"]}; text-mined links kept {stats["textmining_kept"]}, '
        f'dropped {stats["textmining_dropped"]}; {stats["to_review"]} findings to review')


if __name__ == '__main__':
    main()
