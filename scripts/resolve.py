#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The conditions of the knowledge graph and the evidence that relates
#           each to the disease. A condition enters by at least one of:
#             clinical       the seed: the curated list of the profile
#                            (related_conditions.yaml), or the relations read from
#                            PubMed by relations.py; a guideline, review or study
#                            whose PubMed abstract states the relation; every quote
#                            is checked against the abstract
#             ontology       the disease ontology of Open Targets: the parents of
#                            the disease, their other children, its children, and
#                            the children of the diseases the profile names; a
#                            parent with more children than --max-children is too
#                            broad a class and is left out (logged)
#             colocalisation a GWAS credible set of the condition colocalises with a
#                            credible set of a GWAS region of the disease (eCAVIAR
#                            CLPP >= --clpp-min, or coloc H4 >= --h4-min)
#             lead-variant overlap  for a published region without a credible set
#                            of the disease in Open Targets: a GWAS credible set of
#                            the condition contains the region's lead variant with
#                            PIP >= --clpp-min (the CLPP the pair would reach if the
#                            lead were the causal variant)
#           Shared-signal relations are graded by the number of the disease's
#           regions shared, by either route; regions are named as in
#           ot_common.disease_regions (the published gene, else the top
#           locus-to-gene gene, else the position).
#           Every seed identifier is checked against Open Targets. Each condition
#           gets a genetics mode: 'full' (its GWAS and curated / Mendelian
#           evidence are fetched), 'coloc_only' (only the colocalising credible
#           sets: shared-signal traits, conditions with more than --max-gwas
#           GWAS studies, broad classes, entries with more than
#           --max-descendants descendants in the disease ontology, whose own
#           genetics are those of their subtypes, and conditions related only
#           by evidence of co-occurrence (a category the profile marks
#           cooccurrence) from fewer than --comorbidity-min-pubs publications)
#           or 'none' (no Open Targets record). A broad class is marked (broad_class): related only by shared
#           signals, it leaves the graph (tidy_graph.py).
# Inputs  : the profile, the seed and its PubMed records, the disease's GWAS with
#           colocalisations, its published regions as looked up
# Outputs : conditions.tsv (one row per condition), relations.tsv (one row per
#           piece of evidence), colocalisation.tsv (one row per colocalising pair),
#           overlap.tsv, shared_funnel.json (what each shared-signal filter left out, by
#           region), fetch_list.tsv (the conditions whose genetics are fetched),
#           seed_check.tsv, ontology_log.tsv, not_related.tsv (seed conditions whose
#           clinical evidence all found no association); for a seed read from PubMed (the
#           folder of relations.py), its disputed items and counts, passed on to
#           the curation: relations_review.tsv, relations_reading.json
# ---------------------------------------------------------------------------
import argparse
import json
import re
import shutil
import sys
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from llm import in_text  # noqa: E402
from ot_common import (SHARED_SIGNAL, disease_regions, gql, load_json, load_profile, log, relates,  # noqa: E402
                       about_relation, relation_design, shared_relations,
                       write_tsv)

DISEASES = '''query($ids: [String!]!) { diseases(efoIds: $ids) { id name description parents { id name }
  children { id name } } }'''
N_STUDIES = '''query($ids: [String!]) { studies(diseaseIds: $ids, page: {index: 0, size: 1}) { count } }'''
CO_COLS = ['disease_locus', 'disease_gene', 'l2g_gene', 'disease_study_locus_id', 'node_id', 'condition', 'study_id',
           'study_type', 'trait', 'first_author', 'pub_date', 'pmid', 'project', 'n_samples', 'ancestry',
           'other_study_locus_id', 'other_lead', 'p_mantissa', 'p_exponent', 'method', 'clpp', 'h4',
           'n_coloc_variants', 'beta_ratio_sign']
OV_COLS = ['disease_locus', 'disease_gene', 'variant_id', 'node_id', 'condition', 'study_id', 'trait', 'first_author',
           'pub_date', 'pmid', 'project', 'n_samples', 'ancestry', 'other_study_locus_id', 'other_lead', 'p_mantissa',
           'p_exponent', 'finemapping', 'pip']


COND_COLS = ['node_id', 'name', 'category', 'ot_ids', 'genetics', 'genetics_note', 'mendelian', 'genes', 'seed_key',
             'sources', 'support', 'n_descendants', 'broad_class', 'direction', 'direction_basis']
DESCENDANTS = '''query($ids: [String!]!) { diseases(efoIds: $ids) { id descendants } }'''
REL_COLS = ['node_id', 'evidence_type', 'relation', 'detail', 'quote', 'design', 'design_class', 'population',
            'east_asian', 'pmid', 'citation', 'strength', 'loci', 'n_loci', 'n_studies', 'agreement']


def records(api, ids):
    out = {}
    ids = sorted(set(ids))
    for i in range(0, len(ids), 100):
        for d in gql(api, DISEASES, {'ids': ids[i:i + 100]})['diseases']:
            out[d['id']] = d
    return out


def ancestry(st):
    return '; '.join(f'{p["ldPopulation"]} {100 * p["relativeSampleSize"]:.0f}%'
                     for p in st.get('ldPopulationStructure') or [])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--profile', required=True, help='profile.json of the disease')
    ap.add_argument('--seed', default='', help='the seed: a YAML file, or the folder of relations.py (its '
                    'extracted_relations.yaml); default: the curated list of the profile')
    ap.add_argument('--seed-pubmed', required=True, help="the seed PMIDs' PubMed records")
    ap.add_argument('--core-gwas', required=True, help='gwas.json.gz of the disease fetch (with colocalisations)')
    ap.add_argument('--published', default='', help='published_loci.json.gz of the disease fetch, if any')
    ap.add_argument('--api', required=True)
    ap.add_argument('--clpp-min', type=float, default=0.01)
    ap.add_argument('--h4-min', type=float, default=0.8)
    ap.add_argument('--region-window-kb', type=int, default=500, help='lead variants this close are one region')
    ap.add_argument('--genome-wide-p', type=float, default=5e-8, help='P of a GWAS region of the disease')
    ap.add_argument('--region-min-maf', type=float, default=0.01,
                    help='a lead of a region from Open Targets rarer than this needs P below --rare-variant-p')
    ap.add_argument('--region-min-mac', type=float, default=25,
                    help='expected copies of the minor allele of such a lead in the smaller group of its study')
    ap.add_argument('--rare-variant-p', type=float, default=5e-9, help='P of a rare lead of such a region')
    ap.add_argument('--l2g-min', type=float, default=0.5)
    ap.add_argument('--max-gwas', type=int, default=100, help='a condition with more GWAS studies is coloc_only')
    ap.add_argument('--overlap-pip-min', type=float, default=0.1,
                    help='PIP of a published lead variant in a credible set of a condition (lead-variant overlap)')
    ap.add_argument('--duplicate-study-share', type=float, default=0.8,
                    help='traits related only by shared signals whose GWAS studies overlap this much are one trait')
    ap.add_argument('--max-descendants', type=int, default=100,
                    help='an entry with more descendants in the disease ontology is a broad class (coloc_only)')
    ap.add_argument('--max-children', type=int, default=30, help='a parent with more children is too broad')
    ap.add_argument('--comorbidity-min-pubs', type=int, default=2,
                    help='a condition related only by co-occurrence needs this many publications for its genetics')
    ap.add_argument('--out-dir', default='.')
    a = ap.parse_args()
    P = load_profile(a.profile)
    CORE, short = P.id, P.short
    # the disease itself: its ID and the entries whose evidence it takes (evidence_also) are never related conditions
    OWN = {CORE, *(P.d.get('evidence_also') or [])}
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    seed_path = Path(a.seed) if a.seed else None
    if seed_path and seed_path.is_dir():                              # read from PubMed: pass on the readers' disputes
        for src, dst in (('relations_review.tsv', 'relations_review.tsv'), ('relations.json', 'relations_reading.json')):
            if (seed_path / src).exists():
                shutil.copyfile(seed_path / src, out / dst)
        seed_path = seed_path / 'extracted_relations.yaml'
    seed = yaml.safe_load(seed_path.read_text()) if seed_path else (P.seed or {'conditions': []})
    # the words of the disease's own names: a title that names only the disease does not name a related condition
    own_words = {w for t in (P.name, P.short) for w in re.split(r'[^a-z0-9]+', str(t).lower()) if len(w) >= 3}
    extracted = seed.get('source') == 'extracted'
    pub = {r['pmid']: r for r in load_json(a.seed_pubmed)['records']}

    # ── clinical evidence: every quote must be in its abstract ───────────────
    conds, rel, check, alias = {}, [], [], {}
    for c in seed.get('conditions') or []:
        ids = c.get('ot_ids') or []
        node = ids[0] if ids else f'custom:{c["key"]}'
        alias.update({i: node for i in ids})
        conds[node] = {'node_id': node, 'name': c['name'], 'category': c['category'], 'ot_ids': ' '.join(ids),
                       'genetics': c.get('genetics', 'full' if ids else 'none'), 'genetics_note': '',
                       'mendelian': bool(c.get('mendelian', bool(ids))), 'genes': ' '.join(c.get('genes') or []),
                       'seed_key': c['key'], 'sources': {'clinical'},
                       'support': c.get('support') or ('curated' if not extracted else ''),
                       # what its evidence says of the risk of the disease (relations.py; a curated list: below)
                       'direction': c.get('direction', ''), 'direction_basis': c.get('direction_basis', '')}
        for e in (c.get('evidence') or []) + [dict(x, east_asian=True) for x in (c.get('east_asian') or [])]:
            r = pub.get(str(e['pmid']))
            ok = r is not None and in_text(e['quote'], r['title'] + ' ' + r['abstract'])
            check.append({'node_id': node, 'pmid': str(e['pmid']), 'found': r is not None, 'quote_in_abstract': ok})
            if not ok:
                continue
            rel.append({'node_id': node, 'evidence_type': 'clinical', 'relation': e.get('direction', ''),
                        'detail': e['statement'], 'quote': e['quote'], 'design': e.get('design', ''),
                        'design_class': relation_design(r, e.get('design', ''),
                                                        about_relation(c['name'], r.get('title', ''), own_words)),
                        'population': e.get('population', ''), 'east_asian': bool(e.get('east_asian')),
                        'pmid': str(e['pmid']), 'citation': f'{r["first_author"]} {r["year"]}, {r["journal"]}',
                        'strength': '', 'loci': '', 'agreement': e.get('agreement', 'curated')})
    chk = pd.DataFrame(check, columns=['node_id', 'pmid', 'found', 'quote_in_abstract'])
    bad = chk[~chk.quote_in_abstract.astype(bool)]
    assert not len(bad), f'seed quotes not found in their PubMed abstracts:\n{bad.to_string()}'
    no_ev = [n for n in conds if not any(r['node_id'] == n for r in rel)]
    assert not no_ev, f'seed conditions without verified evidence: {no_ev}'
    # a condition whose clinical evidence all found no association is not related to the disease: it leaves the seed
    # (it can still enter by the ontology or a shared signal below, without that evidence) and is listed
    not_related = [{'node_id': n, 'name': conds[n]['name'],
                    'items': sum(r['node_id'] == n for r in rel), 'reason': 'every clinical item found no association'}
                   for n in conds if not any(r['node_id'] == n and relates('clinical', r['relation']) for r in rel)]
    gone = {x['node_id'] for x in not_related}
    conds = {n: c for n, c in conds.items() if n not in gone}
    # a condition of a curated list: the direction its rows that relate it state (both: inconsistent, not weighed)
    for n, c in conds.items():
        if not c['direction']:
            d = {r['relation'] for r in rel if r['node_id'] == n and relates('clinical', r['relation'])}
            up, down = 'increased risk' in d, 'decreased risk' in d
            c['direction'] = ('inconsistent' if up and down else 'increased risk' if up else 'decreased risk' if down
                              else 'associated')
            if up and down:
                c['direction_basis'] = 'both directions in the curated list; not weighed'
    rel = [r for r in rel if r['node_id'] not in gone]

    # ── ontology ─────────────────────────────────────────────────────────────
    O = P.ontology
    recs = records(a.api, [CORE] + [i for c in seed.get('conditions') or [] for i in (c.get('ot_ids') or [])])
    missing = [i for c in seed.get('conditions') or [] for i in (c.get('ot_ids') or []) if i not in recs]
    assert not missing, f'seed IDs not in Open Targets: {missing}'
    core = recs[CORE]
    onto, olog = [], []
    for p in core['parents'] if O.get('parents', True) else []:
        kids = records(a.api, [p['id']])[p['id']]['children']
        if len(kids) > a.max_children:
            olog.append({'node_id': p['id'], 'name': p['name'], 'rule': f'parent of {short}',
                         'decision': f'left out: {len(kids)} children, more than {a.max_children} (too broad a class)'})
            continue
        onto.append((p['id'], f'parent of {short}'))
        for s in kids if O.get('siblings', True) else []:
            if s['id'] not in OWN:
                onto.append((s['id'], f'sibling of {short} (child of {p["name"]})'))
    for ch in core['children'] if O.get('children', True) else []:
        onto.append((ch['id'], f'child of {short}'))
    for pid in O.get('children_of') or []:
        r = records(a.api, [pid])[pid]
        for ch in r['children']:
            if ch['id'] not in OWN:
                onto.append((ch['id'], f'child of {r["name"]}'))
    seen = set()
    onto = [(i, w) for i, w in onto if i not in OWN and not (i in seen or seen.add(i))]
    orec = records(a.api, [i for i, _ in onto]) if onto else {}
    for i0, why in onto:
        i = alias.get(i0, i0)
        if i not in conds:
            conds[i] = {'node_id': i, 'name': orec[i0]['name'], 'category': P.ontology_category, 'ot_ids': i0,
                        'genetics': 'full', 'genetics_note': '', 'mendelian': True, 'genes': '', 'seed_key': '',
                        'sources': set(), 'support': ''}
        conds[i]['sources'].add('ontology')
        rel.append({'node_id': i, 'evidence_type': 'ontology', 'relation': why, 'detail': f'Open Targets disease '
                    f'ontology: {orec[i0]["name"]} is the {why}', 'quote': '', 'design': '', 'population': '',
                    'east_asian': False, 'pmid': '', 'citation': 'EFO / MONDO via Open Targets', 'strength': '',
                    'loci': '', 'agreement': ''})

    # ── colocalisation with the credible sets of the disease's GWAS regions ──
    g = load_json(a.core_gwas)
    published = load_json(a.published) if a.published else None
    regions = disease_regions(g, published, P, a.l2g_min, a.region_window_kb, a.genome_wide_p, a.region_min_maf,
                              None, a.region_min_mac, a.rare_variant_p)
    region_of = {cs: r for r in regions for cs in r['study_locus_ids']}
    own_studies = {s['id'] for s in g['studies']}       # the disease's own GWAS (with its descendants): not a relation
    # what each filter leaves out, region by region (shared_funnel.json; tidy_graph.py adds its own steps): pairs of
    # credible sets as Open Targets lists them, then one row per pair and Open Targets entry its study is mapped to
    funnel = {r['gene']: {'method': 'colocalization' if r['study_locus_ids'] else 'lead-variant overlap' if
                          r['source'] == 'published' and not r['disease_study_locus_id'] else '',
                          'listed': 0, 'below': 0, 'not_other_gwas': 0, 'no_entry': 0, 'itself': 0, 'rows': 0,
                          'joint': 0, 'once': 0} for r in regions}
    co = []
    for s in g['studies']:
        for c in s['credible_sets']:
            r = region_of.get(c['studyLocusId'])
            if r is None:          # a credible set below genome-wide significance names no region
                continue
            f = funnel[r['gene']]
            locus = (c['variant'] or {}).get('rsIds', [''])[0] or c['variant']['id']
            l2g = ((c.get('l2GPredictions') or {}).get('rows') or [{}])[0].get('target', {}).get('approvedSymbol', '')
            for x in c.get('colocalisation') or []:
                f['listed'] += 1
                if not ((x.get('clpp') or 0) >= a.clpp_min or (x.get('h4') or 0) >= a.h4_min):
                    f['below'] += 1
                    continue
                o = x['otherStudyLocus']
                st = o['study']
                # another study's GWAS: not the disease's own, not a molecular QTL study (no disease to relate)
                if st['id'] in own_studies or (x.get('rightStudyType') or 'gwas') != 'gwas':
                    f['not_other_gwas'] += 1
                    continue
                if not st.get('diseases'):
                    f['no_entry'] += 1
                    continue
                for d in st['diseases']:
                    if d['id'] in OWN:
                        f['itself'] += 1
                        continue
                    f['rows'] += 1
                    co.append({'disease_locus': locus, 'disease_gene': r['gene'], 'l2g_gene': l2g,
                               'disease_study_locus_id': c['studyLocusId'],
                               'node_id': d['id'], 'condition': d['name'], 'study_id': st['id'],
                               'study_type': x['rightStudyType'], 'trait': st.get('traitFromSource') or '',
                               'first_author': st.get('publicationFirstAuthor') or '',
                               'pub_date': st.get('publicationDate') or '', 'pmid': st.get('pubmedId') or '',
                               'project': st.get('projectId') or '', 'n_samples': st.get('nSamples'),
                               'ancestry': ancestry(st), 'other_study_locus_id': o['studyLocusId'],
                               'other_lead': ((o.get('variant') or {}).get('rsIds') or [''])[0]
                               or (o.get('variant') or {}).get('id', ''),
                               'p_mantissa': o.get('pValueMantissa'), 'p_exponent': o.get('pValueExponent'),
                               'method': x['colocalisationMethod'], 'clpp': x.get('clpp'), 'h4': x.get('h4'),
                               'n_coloc_variants': x.get('numberColocalisingVariants'),
                               'beta_ratio_sign': x.get('betaRatioSignAverage')})
    co = pd.DataFrame(co, columns=CO_COLS)

    # ── lead-variant overlap, for the published regions without a credible set of the disease ──
    ov = []
    for r in regions:
        if r['source'] != 'published' or r['disease_study_locus_id'] or r['study_locus_ids']:
            continue
        f = funnel[r['gene']]
        for x in r.get('overlaps') or []:
            f['listed'] += 1
            if (x.get('pip') or 0) < a.overlap_pip_min:
                f['below'] += 1
                continue
            st = x['study']
            if st['id'] in own_studies:
                f['not_other_gwas'] += 1
                continue
            if not st.get('diseases'):
                f['no_entry'] += 1
                continue
            for d in st['diseases']:
                if d['id'] in OWN:
                    f['itself'] += 1
                    continue
                f['rows'] += 1
                ov.append({'disease_locus': r['rsid'], 'disease_gene': r['gene'], 'variant_id': r['variant_id'],
                           'node_id': d['id'], 'condition': d['name'], 'study_id': st['id'],
                           'trait': st.get('traitFromSource') or '', 'first_author': st.get('publicationFirstAuthor') or '',
                           'pub_date': st.get('publicationDate') or '', 'pmid': st.get('pubmedId') or '',
                           'project': st.get('projectId') or '', 'n_samples': st.get('nSamples'),
                           'ancestry': ancestry(st), 'other_study_locus_id': x['studyLocusId'],
                           'other_lead': ((x.get('variant') or {}).get('rsIds') or [''])[0]
                           or (x.get('variant') or {}).get('id', ''),
                           'p_mantissa': x.get('pValueMantissa'), 'p_exponent': x.get('pValueExponent'),
                           'finemapping': x.get('finemappingMethod') or '', 'pip': x.get('pip')})
    ov = pd.DataFrame(ov, columns=OV_COLS)

    # ── what a GWAS study can tell: one trait's own signal, counted once ──────
    # a study that tests several traits jointly (a multi-trait or pleiotropy analysis) gives no trait's own signal
    joint = [re.compile(x, re.I) for x in (P.curation.get('joint_trait_patterns') or [])]
    is_joint = lambda t: any(rx.search(str(t)) for rx in joint)  # noqa: E731
    # boolean masks of explicit dtype: an empty table's mask would otherwise select columns, not rows
    jc, jo = co.trait.map(is_joint).astype(bool), ov.trait.map(is_joint).astype(bool)
    n_joint = int(jc.sum() + jo.sum())
    for gname, n in pd.concat([co[jc].disease_gene, ov[jo].disease_gene]).value_counts().items():
        funnel[gname]['joint'] += int(n)
    co, ov = co[~jc].copy(), ov[~jo].copy()
    # a study counts once: one that Open Targets maps to a condition related by clinical or ontology evidence and
    # to other conditions supports that condition only (a PE study also mapped to 'pulmonary infarction')
    related = set(conds)

    def once(df):
        if not len(df):
            return df, 0
        node = df.node_id.map(lambda n: alias.get(n, n))
        by_study = node.groupby(df.study_id).agg(set)
        keep = pd.Series([n in related or not (by_study[s] & related) for s, n in zip(df.study_id, node)],
                         index=df.index, dtype=bool)
        for gname, n in df[~keep].disease_gene.value_counts().items():
            funnel[gname]['once'] += int(n)
        return df[keep].copy(), int((~keep).sum())
    co, n_once_co = once(co)
    ov, n_once_ov = once(ov)
    # traits related only by shared signals whose studies nearly coincide (one trait under two Open Targets entries)
    # are one trait, named by the one with the most studies
    sets = {}
    for df in (co, ov):
        for s_, n_ in zip(df.study_id, df.node_id):
            n2 = alias.get(n_, n_)
            if n2 and n2 not in related:
                sets.setdefault(n2, set()).add(s_)
    order = sorted(sets, key=lambda n: (-len(sets[n]), n))
    dup = {}
    for k, x in enumerate(order):
        if x in dup:
            continue
        for y in order[k + 1:]:
            if y not in dup and len(sets[x] & sets[y]) / len(sets[x] | sets[y]) >= a.duplicate_study_share:
                dup[y] = x
    alias.update(dup)
    log(f'shared signals: {n_joint} rows of joint-trait studies left out; {n_once_co + n_once_ov} rows of studies '
        f'counted for a condition related otherwise; {len(dup)} traits merged as near duplicates')

    # ── shared-signal relations, graded by the disease's regions shared ──────
    shared_rows, regions_all, shared_names = shared_relations(co, ov, alias, short)
    for x in shared_rows:
        i = x['node_id']
        if i not in conds:
            conds[i] = {'node_id': i, 'name': shared_names[i], 'category': SHARED_SIGNAL,
                        'ot_ids': ' '.join([i] + sorted(y for y, x_ in dup.items() if x_ == i)),
                        'genetics': 'coloc_only', 'genetics_note': 'related only by a shared GWAS signal',
                        'mendelian': False, 'genes': '', 'seed_key': '', 'sources': set(), 'support': ''}
        conds[i]['sources'].add(x['evidence_type'])
        rel.append(x)

    # ── broad classes: an entry with more descendants in the disease ontology than --max-descendants (heart
    # disorder, protein measurement) holds, as its own genetics, those of its many subtypes ──
    first = sorted({c['ot_ids'].split()[0] for c in conds.values() if c['ot_ids']})
    n_desc = {}
    for i in range(0, len(first), 50):
        for d in gql(a.api, DESCENDANTS, {'ids': first[i:i + 50]})['diseases']:
            n_desc[d['id']] = len(d.get('descendants') or [])
    for c in conds.values():
        c['n_descendants'] = n_desc.get(c['ot_ids'].split()[0], 0) if c['ot_ids'] else 0
        c['broad_class'] = c['n_descendants'] > a.max_descendants
        if c['broad_class'] and c['genetics'] == 'full':
            c['genetics'], c['genetics_note'] = 'coloc_only', (
                f'a broad class ({c["n_descendants"]} subtypes in the disease ontology, more than '
                f'{a.max_descendants}): its own genetics are those of its subtypes')

    # ── a condition related only by co-occurrence (a comorbidity), from too few publications, keeps its shared
    # signals only: its own genetics say nothing of the disease's causes until the relation is established ──
    cooccurrence = {c['name'] for c in P.categories if c.get('cooccurrence')}
    for n, c in conds.items():
        if c['genetics'] == 'full' and c['category'] in cooccurrence and 'ontology' not in c['sources']:
            pubs = {r['pmid'] for r in rel if r['node_id'] == n and r['evidence_type'] == 'clinical'
                    and relates('clinical', r['relation'])}
            if len(pubs) < a.comorbidity_min_pubs:
                c['genetics'], c['genetics_note'] = 'coloc_only', (
                    f'{c["category"].lower()} from {len(pubs)} publication{"s" if len(pubs) != 1 else ""}, fewer '
                    f'than {a.comorbidity_min_pubs}: its own genetics are not fetched')

    # ── depth of genetic retrieval: a condition with too many GWAS studies keeps its shared signals only ──
    for c in conds.values():
        if c['genetics'] == 'full' and c['ot_ids']:
            n = gql(a.api, N_STUDIES, {'ids': [c['ot_ids'].split()[0]]})['studies']['count']
            if n > a.max_gwas:
                c['genetics'], c['genetics_note'] = 'coloc_only', f'{n} GWAS studies, more than {a.max_gwas}'

    # explicit columns: a disease can have no related condition at all
    C = pd.DataFrame(conds.values(), columns=COND_COLS)
    C[['direction', 'direction_basis']] = C[['direction', 'direction_basis']].fillna('')   # ontology, shared signal
    C['sources'] = C.sources.map(lambda s: ' '.join(sorted(s)))
    R = pd.DataFrame(rel, columns=REL_COLS)
    R['design_class'] = R.design_class.fillna('')
    # clinical items that relate the condition (contrary findings, kept beside them, are not counted)
    rel_clin = R[(R.evidence_type == 'clinical')
                 & pd.Series([relates(t, r) for t, r in zip(R.evidence_type, R.relation)], index=R.index, dtype=bool)]
    C['n_clinical'] = C.node_id.map(rel_clin.groupby('node_id').size()).fillna(0).astype(int)
    C['coloc_loci'] = C.node_id.map(R[R.evidence_type == 'colocalisation'].set_index('node_id').n_loci).fillna(0)
    C['coloc_loci'] = C.coloc_loci.astype(int)
    C['overlap_loci'] = C.node_id.map(R[R.evidence_type == 'lead-variant overlap'].set_index('node_id').n_loci)
    C['overlap_loci'] = C.overlap_loci.fillna(0).astype(int)
    C['shared_regions'] = C.node_id.map(lambda n: ' '.join(sorted(regions_all.get(n, ()))))
    C['n_shared_regions'] = C.node_id.map(lambda n: len(regions_all.get(n, ())))
    C = C.sort_values(['category', 'name'])
    fetch = C[C.genetics.eq('full') & C.ot_ids.ne('')][['node_id', 'mendelian']]
    write_tsv(C, out / 'conditions.tsv')
    write_tsv(R, out / 'relations.tsv')
    write_tsv(co, out / 'colocalisation.tsv')
    write_tsv(ov, out / 'overlap.tsv')
    (out / 'shared_funnel.json').write_text(json.dumps(funnel, indent=1))
    write_tsv(chk, out / 'seed_check.tsv')
    write_tsv(pd.DataFrame(not_related, columns=['node_id', 'name', 'items', 'reason']), out / 'not_related.tsv')
    write_tsv(pd.DataFrame(olog, columns=['node_id', 'name', 'rule', 'decision']), out / 'ontology_log.tsv')
    fetch.assign(mendelian=fetch.mendelian.map(lambda x: 'true' if x else 'false')).to_csv(
        out / 'fetch_list.tsv', sep='\t', index=False, header=False)
    log(f'{len(C)} conditions ({C.sources.str.contains("clinical").sum()} clinical, '
        f'{C.sources.str.contains("ontology").sum()} ontology, {C.sources.str.contains("colocalisation").sum()} '
        f'colocalisation, {C.sources.str.contains("overlap").sum()} lead-variant overlap); {len(R)} evidence rows; '
        f'{len(fetch)} genetics fetches; {len(regions)} GWAS regions; {len(olog)} ontology classes left out')


if __name__ == '__main__':
    main()
