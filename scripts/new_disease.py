#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : A new disease from its Open Targets disease ID (EFO_, MONDO_, ...;
#           find it at https://platform.opentargets.org): writes
#           diseases/<key>/profile.yaml, every value with the rule that chose
#           it, ready to run (or to edit first). The ID is required: a disease
#           is never chosen from a name, which can name several.
#             name      the Open Targets name
#             short     --short, else the shortest acronym among its exact
#                       synonyms that spells its name, else the name itself (at
#                       most 11 characters), else its initials: made up, marked
#                       for checking and never searched
#             key       the short name in lower case (the folder name)
#             terms     the PubMed block: the name (always; 'A, B' also as 'B A'),
#                       the exact synonyms, the short name, and the names of the
#                       entries whose evidence it takes, in title or abstract,
#                       without a term Open Targets maps to a gene or to another
#                       disease, acronyms of two characters or fewer, and
#                       database forms of a name; an acronym only beside a word
#                       of the name
#           Everything else follows the defaults of the pipeline (README,
#           "Add a disease"). Open Targets is read; nothing else.
# Usage   : python3 scripts/new_disease.py MONDO_0017148 [--key hpah] [--short HPAH] [--force]
# ---------------------------------------------------------------------------
import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import gql, release  # noqa: E402
from profile import MADE_UP, evidence_also, search_terms, short_from  # noqa: E402

API = 'https://api.platform.opentargets.org/api/v4/graphql'
DISEASE = '''query($id: String!) { disease(efoId: $id) { id name synonyms { relation terms }
  parents { id name } children { id name } } }'''
ALSO = '''query($id: String!) { disease(efoId: $id) { id name synonyms { relation terms } } }'''
ID = re.compile(r'^(EFO|MONDO|HP|Orphanet|OTAR|DOID|GO|OBA|NCIT)_\d+$')


def find(api, disease_id):
    """The Open Targets record of a disease ID."""
    if not ID.match(disease_id):
        raise SystemExit(f'"{disease_id}" is not an Open Targets disease ID (EFO_, MONDO_, HP_, Orphanet_, ...): find the '
                         'ID of the disease at https://platform.opentargets.org and give it here')
    rec = gql(api, DISEASE, {'id': disease_id})['disease']
    if not rec:
        raise SystemExit(f'{disease_id}: not in Open Targets')
    return rec


def q(s):
    """A YAML string in single quotes."""
    return "'" + str(s).replace("'", "''") + "'"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('disease', help='the Open Targets disease ID (EFO_, MONDO_, HP_, Orphanet_, ...)')
    ap.add_argument('--key', default='', help='the folder name (default: the short name in lower case)')
    ap.add_argument('--short', default='', help='the short name (at most 11 characters), when the rule finds none or '
                    'makes one up')
    ap.add_argument('--diseases-dir', default=str(Path(__file__).resolve().parent.parent / 'diseases'))
    ap.add_argument('--api', default=API)
    ap.add_argument('--force', action='store_true', help='replace an existing profile.yaml')
    a = ap.parse_args()
    rec = find(a.api, a.disease)
    syn = [t for s in rec['synonyms'] or [] if s['relation'] == 'hasExactSynonym' for t in s['terms']]
    short, short_rule = (a.short, 'given (--short)') if a.short else short_from(rec['name'], syn)
    if not short or len(short) > 11 or "'" in short:
        raise SystemExit(f'{rec["id"]} {rec["name"]}: {short_rule if not short else "the short name is longer than 11 "
                                                      "characters"}; give --short, the abbreviation the literature uses '
                         '(at most 11 characters, no apostrophe)')
    also, warn = evidence_also(a.api, rec)
    also_recs = [gql(a.api, ALSO, {'id': x['id']})['disease'] for x in also]
    terms, notes = search_terms(a.api, rec['name'], syn, short, rec['id'], short_rule, also_recs)
    key = a.key or re.sub(r'[^a-z0-9]+', '_', short.lower()).strip('_')
    out = Path(a.diseases_dir) / key / 'profile.yaml'
    if out.exists() and not a.force:
        raise SystemExit(f'{out} exists; give another --key, or --force to replace it')
    rel = release(a.api)['data_version']
    lines = [f'# The disease profile of {rec["name"]}, written by scripts/new_disease.py from Open Targets {rel}.',
             '# Every value below follows the rule given beside it; edit any of them, then run the pipeline. Only `id`',
             '# is required; the related conditions, GWAS regions, gene themes and curation follow the defaults of the',
             '# pipeline (README, "Add a disease").',
             f'id: {rec["id"]}            # the Open Targets disease ID (required)',
             f'name: {q(rec["name"])}     # its Open Targets name',
             f'short: {q(short)}          # {short_rule}', '',
             *(['# CHECK: the short name above is made up; replace it with the abbreviation the literature uses', '']
               if MADE_UP in short_rule else []),
             'pubmed:',
             '  # the disease in title or abstract: its name, exact synonyms and short name'
             + (', without the terms that would find other things:' if any('left out' in n for n in notes) else ''),
             *[f'  #   {n}' for n in notes],
             f'  terms: {q(terms)}',
             '  # which studies count as studies of the disease, given to the language-model readers (optional):',
             '  # scope: <one or two sentences>', '',
             '# The neighbors of the disease in the ontology (all on by default); a parent with too many children is',
             '# left out by the pipeline (MaxOntologyChildren).',
             f'# parents:  {", ".join(p["name"] for p in rec["parents"] or []) or "none"}',
             f'# children: {", ".join(c["name"] for c in rec["children"] or []) or "none"}',
             'ontology:',
             '  category: Disease group   # the category of the ontology neighbors in the products', '',
             '# Other Open Targets entries whose evidence (targets, text-mined literature, GWAS, curated genes) is',
             '# fetched with the disease\'s own; the products name them. Proposed: an entry that a part of a grouping',
             '# name ("A and/or B") names, other than the disease and its descendants.',
             *[f'#   {w}' for w in warn],
             *([f'evidence_also: [{", ".join(x["id"] for x in also)}]'] +
               [f'#   {x["id"]} {x["name"]}: the name names "{x["part"]}"' for x in also] if also else
               ['# evidence_also: []'])]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text('\n'.join(lines) + '\n')
    print(f'{out}\n  {rec["id"]} {rec["name"]} ({short})\n  PubMed: {terms}' + ''.join(f'\n  {n}' for n in notes)
          + ''.join(f'\n  NOTE: {w}' for w in warn)
          + (f'\n  NOTE: the short name {short} is made up from the initials; replace it in the profile (or give '
             '--short) with the abbreviation the literature uses' if MADE_UP in short_rule else '')
          + ''.join(f'\n  evidence also from {x["id"]} {x["name"]} (the name names "{x["part"]}")' for x in also)
          + f'\nNext: nextflow run main.nf --diseases {key} --FetchDate $(date +%F) -c conf/local.config')


if __name__ == '__main__':
    main()
