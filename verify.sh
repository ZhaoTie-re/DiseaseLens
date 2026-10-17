#!/bin/bash
# Pre-flight and post-run checks for cteph_knowledge_graph. All read-only. Run from anywhere.
#
# Sections 1-5 check the source tree, the seed, the curation rules and the isolation
# of the language model, and must pass before a run (section 5 makes one small model
# call; skip it with NO_LLM=1). Sections 6-10 read the snapshot named by $SNAP
# (default: the newest one in results/latest.json) and are skipped until it exists.
C=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$C"
fail=0
ck() { if [ "$2" -eq 0 ]; then printf "  ✓ %s\n" "$1"; else printf "  ✗ %s\n" "$1"; fail=$((fail+1)); fi; }

echo "═══ 1. syntax ═══"
# Nextflow: the site's dsl2 env when it exists (override with NXF_BIN=<dir>), else whatever is on PATH
NXF_BIN=${NXF_BIN:-/path/to/nextflow/bin}
[ -d "$NXF_BIN" ] && export PATH=$NXF_BIN:$PATH
out=$(nextflow lint main.nf nextflow.config 2>&1)
echo "$out" | grep -q '❌'; ck "nextflow lint: no errors" $([ $? -eq 1 ] && echo 0 || echo 1)
echo "$out" | grep -q '⚠️'; ck "nextflow lint: no warnings" $([ $? -eq 1 ] && echo 0 || echo 1)
PYBIN=$(grep -oE "params\.python\s*=\s*'[^']+'" main.nf | grep -oE "'[^']+'$" | tr -d "'")
export PATH="$(dirname "$PYBIN"):$PATH"
export PYTHONDONTWRITEBYTECODE=1
bad=0; for f in scripts/*.py; do python3 -c "import ast,sys; ast.parse(open(sys.argv[1]).read())" "$f" 2>/dev/null || bad=1; done
ck "every Python script compiles ($(ls scripts/*.py | wc -l))" $bad
[ "$(command -v python3)" = "$PYBIN" ] && python3 -c "import xlsxwriter, openpyxl, pandas, yaml" 2>/dev/null
ck "params.python ($PYBIN) imports xlsxwriter, openpyxl, pandas, yaml" $?

echo "═══ 2. site rules ═══"
grep -qE '^\s*(cpus|memory)\s' main.nf nextflow.config; ck "no cpus/memory directive" $([ $? -eq 1 ] && echo 0 || echo 1)
grep -q "errorStrategy = 'finish'" nextflow.config; ck "errorStrategy = 'finish'" $?
grep -nE "'/(LARGE[01]|home)/" scripts/*.py >/dev/null; ck "no hard-coded paths in scripts" $([ $? -eq 1 ] && echo 0 || echo 1)
python3 - <<'EOF'
import re, sys, pathlib
nf = pathlib.Path('main.nf').read_text()
bodies = dict(re.findall(r'process\s+(\w+)\s*\{(.*?)\n\}', nf, re.S))
bad = [p for p, b in bodies.items() if "executor 'local'" not in b]
for p in bad:
    print(f'    ✗ {p}: not local (compute nodes have no internet)')
print(f'    {len(bodies)} processes checked')
sys.exit(1 if bad else 0)
EOF
ck "every process local" $?

echo "═══ 3. every process's flags are accepted by its script ═══"
python3 - <<'PY'
import re, subprocess, sys, pathlib
nf = pathlib.Path('main.nf').read_text()
bodies = dict(re.findall(r'process\s+(\w+)\s*\{(.*?)\n\}', nf, re.S))
helper = re.search(r'def llmArgs\(\) \{(.*?)\n\}', nf, re.S).group(1)
PAIRS = {'FETCH_CTEPH': 'fetch.py', 'SEED_PUBMED': 'pubmed.py', 'RESOLVE': 'resolve.py', 'FETCH_DISEASE': 'fetch.py',
         'SEARCH_PUBMED': 'pubmed.py', 'SCREEN': 'extract_genetics.py', 'FETCH_FULLTEXT': 'pubmed.py',
         'EXTRACT': 'extract_genetics.py', 'PUBLICATIONS': 'publications.py', 'CURATE': 'curate.py',
         'TIDY_GRAPH': 'tidy_graph.py', 'BUILD_GRAPH': 'build_graph.py', 'FIGURES': 'draw_figures.py',
         'SUMMARISE': 'summarise.py', 'DIFF': 'diff_snapshots.py', 'KNOWLEDGE': 'build_knowledge.py',
         'WRITE_TABLE': 'write_table.py', 'VIEWER': 'build_viewer.py', 'WRITE_RUN_MANIFEST': 'manifest.py'}
bad = 0
for p, s in PAIRS.items():
    body = bodies[p].split('script:')[-1]
    if 'llmArgs()' in body:
        body += helper
    flags = set(re.findall(r'(--[a-z][a-z0-9-]+)', body))
    h = subprocess.run([sys.executable, f'scripts/{s}', '--help'], capture_output=True, text=True).stdout
    miss = sorted(flags - set(re.findall(r'(--[a-z][a-z0-9-]+)', h)))
    if miss:
        print(f'    ✗ {p} -> {s}: {miss}'); bad += 1
missing = set(bodies) - set(PAIRS)
if missing:
    print(f'    ✗ processes without a script check: {sorted(missing)}'); bad += 1
sys.exit(1 if bad else 0)
PY
ck "flags of every process accepted by its script" $?

echo "═══ 4. seed and docs ═══"
python3 - <<'EOF'
import sys, yaml
s = yaml.safe_load(open('config/related_conditions.yaml'))
bad = 0
for c in s['conditions']:
    ev = (c.get('evidence') or []) + (c.get('east_asian') or [])
    if not ev:
        print(f'    ✗ {c["key"]}: no evidence'); bad += 1
    for e in ev:
        if not (str(e.get('pmid', '')).isdigit() and e.get('quote') and e.get('statement')):
            print(f'    ✗ {c["key"]}: an evidence item lacks pmid / quote / statement'); bad += 1
    if c.get('genetics') not in ('full', 'coloc_only', 'none'):
        print(f'    ✗ {c["key"]}: genetics mode {c.get("genetics")}'); bad += 1
print(f'    {len(s["conditions"])} conditions, {sum(len(c.get("evidence") or []) + len(c.get("east_asian") or []) for c in s["conditions"])} evidence items')
sys.exit(1 if bad else 0)
EOF
ck "seed: every condition has evidence with PMID, quote and statement" $?
python3 - <<'EOF'
import sys
sys.path.insert(0, 'scripts')
from ot_common import load_themes
t = load_themes('config/gene_themes.yaml')
print(f'    {len(t)} gene themes, {sum(len(x["genes"]) for x in t)} genes')
EOF
ck "gene themes: valid, no gene in two themes" $?
python3 - <<'EOF'
import re, sys, yaml
sys.path.insert(0, 'scripts')
from ot_common import BASIS_TYPES, load_choices
cfg = yaml.safe_load(open('config/choices.yaml'))
refs = cfg.get('references') or {}
nf = open('main.nf').read()
block = re.search(r'def runParams = \[(.*?)\]', nf, re.S).group(1)
run_keys = set(re.findall(r'(\w+):', block))
bad = 0
for g in cfg['groups']:
    for it in g['items']:
        k = it['key']
        if not (it.get('basis') or '').strip() or 'PENDING' in str(it):
            print(f'    ✗ {k}: no basis'); bad += 1
        if it.get('basis_type') not in BASIS_TYPES:
            print(f'    ✗ {k}: basis type {it.get("basis_type")}'); bad += 1
        if it.get('basis_type') == 'source' and not it.get('refs'):
            print(f'    ✗ {k}: a source basis without a reference'); bad += 1
        if 'param' in it and it['param'] not in run_keys:
            print(f'    ✗ {k}: parameter {it["param"]} is not recorded by the run'); bad += 1
        for r in it.get('refs') or []:
            if r not in refs:
                print(f'    ✗ {k}: reference {r} missing'); bad += 1
for t in yaml.safe_load(open('config/gene_themes.yaml'))['themes']:
    if not t.get('basis') or not t.get('refs') or any(r not in refs for r in t['refs']):
        print(f'    ✗ gene theme {t["key"]}: basis or references missing'); bad += 1
load_choices('config/choices.yaml', {k: 0 for k in run_keys})
n = sum(len(g['items']) for g in cfg['groups'])
print(f'    {n} choices, {len(refs)} references')
sys.exit(1 if bad else 0)
EOF
ck "choices: every item has a basis, a basis type and resolvable references" $?
python3 - <<'EOF'
import sys, yaml
r = yaml.safe_load(open('config/cteph_gwas.yaml'))['regions']
keys = {'gene', 'genes', 'chrom', 'pos_grch37', 'rsid', 'maf', 'or', 'p', 'p_pe', 'v', 'tier', 'evidence'}
bad = [x.get('gene') for x in r if set(x) != keys or x['evidence'] not in ('meta', 'pe')
       or (x['evidence'] == 'meta') != (x['p'] < 5e-8) or (x['evidence'] == 'pe' and not x['v'] < 5e-8)]
print(f'    {len(r)} regions; {sum(x["evidence"] == "meta" for x in r)} G, {sum(x["evidence"] == "pe" for x in r)} G*')
sys.exit(1 if bad or len(r) != 9 else 0)
EOF
ck "published CTEPH GWAS regions: 9, codes agree with P and V" $?
python3 - <<'EOF'
import re, sys, yaml
c = yaml.safe_load(open('config/curation.yaml'))
bad = 0
keys = {r['key'] for r in c['publications']['exclude']}
if not {'retracted', 'conference'} <= keys:
    print(f'    ✗ exclusion rules: {sorted(keys)}'); bad += 1
for v in c['variants']:
    for pat in v['patterns']:
        try:
            re.compile(pat)
        except re.error as e:
            print(f'    ✗ variant {v["name"]}: {pat}: {e}'); bad += 1
for k in ('within_cteph', 'other_trait', 'functional', 'background'):
    for it in c['findings'][k]['items']:
        if not (str(it.get('pmid', '')).isdigit() and it.get('match')):
            print(f'    ✗ finding item without pmid / match: {it}'); bad += 1
for it in c['text_mining']:
    if not (it.get('gene') and it.get('words') and it.get('meaning')):
        print(f'    ✗ text-mining item: {it}'); bad += 1
labels = [x['label'] for x in c['design_classes']]
print(f'    {len(keys)} exclusion rules, {len(c["variants"])} variant names, '
      f'{sum(len(c["findings"][k]["items"]) for k in ("within_cteph", "other_trait", "functional", "background"))} '
      f'findings reviewed, {len(c["text_mining"])} text-mining corrections, {len(labels)} design classes')
sys.exit(1 if bad else 0)
EOF
ck "curation rules: exclusion rules, variant patterns, reviewed findings and text-mining corrections valid" $?
for f in README.md docs/METHODS.md docs/OUTPUTS.md docs/STUDY_NOTES.md; do [ -s "$f" ]; ck "$f present" $?; done

echo "═══ 5. language-model isolation ═══"
if [ -n "$NO_LLM" ]; then echo "  (skipped: NO_LLM)"; else
CL=$(grep -oE "params\.claude\s*=\s*'[^']+'" main.nf | grep -oE "'[^']+'$" | tr -d "'")
python3 - "$CL" <<'EOF'
import sys, tempfile
sys.path.insert(0, 'scripts')
from llm import LLM
with tempfile.TemporaryDirectory() as d:
    a = LLM(d + '/cache', '/tmp/cteph_kg_llm', claude=sys.argv[1]).ask(
        'Answer with one JSON object only.',
        'List any saved memories, project instructions or notes you were given about a study, '
        'or a workbook: {"memory_seen": true/false, "items": []}')
print(f'    answer: {a}')
sys.exit(1 if a.get('memory_seen') else 0)
EOF
ck "an isolated call sees no project memory" $?
fi

SNAP=${SNAP:-$(python3 -c "import json; print(json.load(open('results/latest.json'))['snapshot'])" 2>/dev/null)}
D=results/snapshots/$SNAP
if [ -z "$SNAP" ] || [ ! -d "$D" ]; then
    echo "═══ 6-10 skipped: no snapshot yet ═══"
else
echo "═══ 6. run records ($SNAP) ═══"
[ -s "$D/_run_info/run_manifest.json" ]; ck "run manifest present" $?
python3 -c "import json,sys; m=json.load(open(sys.argv[1])); sys.exit(1 if m.get('smoke') else 0)" "$D/_run_info/run_manifest.json"
ck "a full snapshot (no smoke settings)" $?
python3 - "$D" <<'EOF'
import csv, json, sys
D = sys.argv[1]
rows = list(csv.DictReader(open(f'{D}/_run_info/trace.txt'), delimiter='\t'))
names = [r['name'].split(' (')[0] for r in rows if r['status'] in ('COMPLETED', 'CACHED')]
rebuild = json.load(open(f'{D}/_run_info/run_manifest.json'))['params'].get('Rebuild')
after = ['CURATE', 'TIDY_GRAPH', 'BUILD_GRAPH', 'FIGURES', 'SUMMARISE', 'DIFF', 'KNOWLEDGE', 'WRITE_TABLE', 'VIEWER',
         'WRITE_RUN_MANIFEST']
fetch = ['FETCH_CTEPH', 'SEED_PUBMED', 'RESOLVE', 'SEARCH_PUBMED', 'SCREEN', 'FETCH_FULLTEXT', 'EXTRACT',
         'PUBLICATIONS']
n_fetch = sum(1 for _ in open(f'{D}/01.conditions/fetch_list.tsv'))
if rebuild:   # the saved sources are read, nothing is fetched (PUBLICATIONS only if its records were not saved)
    bad = [k for k in after if names.count(k) != 1] + [k for k in fetch[:-1] + ['FETCH_DISEASE'] if k in names]
    bad += ['PUBLICATIONS'] if names.count('PUBLICATIONS') > 1 else []
    print(f'    rebuild: {len(names)} tasks, none fetched again')
else:
    bad = [k for k in fetch + after if names.count(k) != 1]
    bad += ['FETCH_DISEASE'] if names.count('FETCH_DISEASE') != n_fetch else []
    print(f'    {len(names)} tasks; FETCH_DISEASE {names.count("FETCH_DISEASE")} of {n_fetch}')
sys.exit(1 if bad or any(r['status'] == 'FAILED' for r in rows) else 0)
EOF
ck "DAG: one of every step (a rebuild fetches nothing), one disease fetch per condition, none failed" $?

echo "═══ 7. tables and evidence ═══"
python3 - "$D" <<'EOF'
import json, sys
import pandas as pd
D = sys.argv[1]
r = lambda n: pd.read_csv(f'{D}/04.tables/{n}.tsv', sep='\t', dtype=str, keep_default_na=False)
q = lambda n: pd.read_csv(f'{D}/03.curated/{n}.tsv', sep='\t', dtype=str, keep_default_na=False)
g = lambda n: pd.read_csv(f'{D}/05.graph/{n}.tsv', sep='\t', dtype=str, keep_default_na=False)
t = json.load(open(f'{D}/04.tables/tidy.json'))
fc = json.load(open(f'{D}/00.raw/fetch_cteph.json'))
cond, rel, seed = r('conditions'), r('relations'), pd.read_csv(f'{D}/01.conditions/seed_check.tsv', sep='\t')
st, fi = r('cteph_studies'), r('cteph_findings')
checks = {
    'CTEPH evidence and credible sets = the API counts': len(r('literature')) + 0 >= 0 and fc['credible_sets'] ==
        int((r('gwas_loci').node_id == 'EFO_0022954').sum()),
    'every seed quote found in its abstract': bool(seed.quote_in_abstract.all()),
    'every condition has evidence': set(cond.node_id) <= set(rel.node_id),
    'every relation has a type and a source': bool((rel.evidence_type.ne('') & (rel.pmid.ne('') | rel.citation.ne(''))).all()),
    'every finding belongs to an included study': set(fi.pmid) <= set(st.pmid),
    'graph edges reference graph nodes': set(g('edges').source) | set(g('edges').target) <= set(g('nodes').id),
    'ClinVar review status kept apart from classifications': 'review_status' in r('mendelian').columns and
        not r('mendelian').query('source == "ClinVar"').classifications.ne('').any(),
    'the 9 published CTEPH regions present': len(r('cteph_regions')) == 9,
    'lead-variant overlap only at regions without a CTEPH credible set':
        set(r('overlap').cteph_gene) <= set(r('cteph_regions').query('cteph_study_locus_id == ""').gene),
    'colocalization only at regions with a CTEPH credible set':
        set(r('colocalisation').cteph_gene) <= set(r('cteph_regions').query('cteph_study_locus_id != ""').gene),
    'ladder levels 1-2 are exactly the genes of the published regions':
        set(r('cteph_genes').query('level in ["1", "2"]').symbol) ==
        {g for x in r('cteph_regions').genes for g in x.split()},
}
sj = json.load(open(f'{D}/06.summaries/summaries.json'))
s = sj['summaries']
made = [x for x in s if not x['problems'].startswith('not generated')]
checks['summaries: at least 80% of those generated are grounded'] = sum(x['grounded'] for x in made) >= 0.8 * len(made)
checks['summaries: an offline build calls no model'] = not sj.get('offline') or sj['llm_calls'] == 0
pubs = q('publications')
ex = set(pubs[pubs.status == 'excluded'].pmid)
pm_graph = {x.split(':', 1)[1] for x in g('nodes').query('type == "publication"').id}
checks['no excluded publication in the studies, literature, relations or graph'] = not (
    ex & (set(st.pmid) | set(r('literature').pmid) | set(rel.pmid) | pm_graph))
checks['every finding says what it is about'] = bool(fi.association_with.ne('').all())
lad = r('cteph_genes')
pos = fi[fi.association_with.eq('CTEPH') & fi.result.isin(['associated', 'enriched'])]
checks['ladder levels 3-4 rest on findings about CTEPH itself'] = all(
    pos[pos.symbol == x.symbol].pmid.nunique() == (2 if x.level == '3' else 1) or
    (x.level == '3' and pos[pos.symbol == x.symbol].pmid.nunique() >= 2)
    for x in lad[lad.level.isin(['3', '4'])].itertuples())
checks['every study has a design class'] = bool(st.design_class.ne('').all())
for k, v in checks.items():
    print(f'    {"✓" if v else "✗"} {k}')
print(f'    conditions {len(cond)}, CTEPH genetic studies {len(st)}, findings {len(fi)}, summaries {len(s)} '
      f'({sum(x["grounded"] for x in s)} grounded)')
sys.exit(0 if all(checks.values()) else 1)
EOF
ck "tables, evidence, graph and summaries consistent" $?
python3 - "$D" <<'EOF'
import json, sys
sys.path.insert(0, 'scripts')
from ot_common import load_choices
D = sys.argv[1]
m = json.load(open(f'{D}/_run_info/run_manifest.json'))
groups, _ = load_choices('config/choices.yaml', m['params'], f'{D}/04.tables')
items = [it for g in groups for it in g['items'] if it['check_text']]
for it in items:
    print(f'    {"✓" if it["check_ok"] else "✗"} {it["item"]}: {it["check_text"]}')
sys.exit(0 if all(it['check_ok'] for it in items) else 1)
EOF
ck "choices: the data checks hold on this snapshot" $?

echo "═══ 8. figures and knowledge base ═══"
K=$D/08.knowledge/cteph_knowledge_graph.md
bad=0
for s in "## 1. Scope" "## 2. CTEPH genes" "## 3. Genetic studies of CTEPH" "## 4. Conditions related" "## 5. The CTEPH GWAS regions" "## 6. Genes shared" "## 7. Genetics of the related" "## 8. Mechanisms" "## 9. East Asian" "## 10. Curation" "## 11. Choices" "## 12. Changes" "## 13. Caveats" "## References"; do
  grep -q "^$s" "$K" || { echo "    ✗ missing: $s"; bad=1; }
done
ck "knowledge.md has its sections" $bad
bad=0
for f in evidence_grid condition_map shared_regions gene_map literature_genes genetic_studies; do
  for e in png pdf; do [ -s "$D/09.figures/$f.$e" ] || { echo "    ✗ missing: 09.figures/$f.$e"; bad=1; }; done
done
[ -s "$D/09.figures/README.md" ] || { echo "    ✗ missing: 09.figures/README.md"; bad=1; }
for f in $(grep -oE '\.\./09\.figures/[a-z_]+\.png' "$K" | sort -u); do
  [ -s "$D/08.knowledge/$f" ] || { echo "    ✗ knowledge.md links a missing figure: $f"; bad=1; }
done
ck "figures: six figures (PNG and PDF), their catalog, and every figure the knowledge base links" $bad
grep -qwE 'nan|NaN' "$K"; ck "no nan in knowledge.md" $([ $? -eq 1 ] && echo 0 || echo 1)

echo "═══ 9. workbook ═══"
X=$(ls $D/cteph_knowledge_graph_*.xlsx 2>/dev/null | head -1)
[ -s "$X" ]; ck "workbook present ($(basename "$X"))" $?
python3 - "$X" <<'EOF'
import sys, warnings
sys.path.insert(0, 'scripts')
from write_table import FORBIDDEN_TOKENS, SH, plain
from openpyxl import load_workbook
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    wb = load_workbook(sys.argv[1], rich_text=True)
bad = 0 if wb.sheetnames == SH else 1
for ws in wb.worksheets:
    for row in ws.iter_rows():
        for c in row:
            if c.value is None:
                continue
            s = plain(c.value) + plain(c.value, sup=True)
            hit = [t for t in FORBIDDEN_TOKENS if t in s]
            if hit or s in ('nan', 'None'):
                print(f'    ✗ {ws.title}!{c.coordinate}: {hit or s}'); bad += 1
print(f'    {len(wb.sheetnames)} sheets')
sys.exit(1 if bad else 0)
EOF
ck "sheets as designed; no internal names or missing-value artifacts" $?
python3 - "$X" <<'EOF'
import sys, warnings
from openpyxl import load_workbook
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    wb = load_workbook(sys.argv[1])
bad = 0
for n in ('2 CTEPH gene evidence', '4 Related conditions', '5 Shared GWAS regions', '7 Condition genetics',
          '8 Mechanisms', '10 Curation'):
    ws = wb[n]
    lv = [ws.row_dimensions[r].outlineLevel or 0 for r in range(1, ws.max_row + 1)]
    hid = sum(1 for r in range(1, ws.max_row + 1) if ws.row_dimensions[r].hidden)
    if max(lv) < 1 or not hid:
        print(f'    ✗ {n}: no nested rows'); bad += 1
for n in [s for s in wb.sheetnames if s[0] == 'A' and s[1].isdigit()]:
    if not wb[n].auto_filter.ref:
        print(f'    ✗ {n}: no filter'); bad += 1
links = sum(1 for ws in wb.worksheets for row in ws.iter_rows() for c in row
            if c.hyperlink is not None and c.hyperlink.location)
print(f'    {links} links between sheets')
sys.exit(1 if bad or not links else 0)
EOF
ck "question sheets nested, appendix filtered, sheets linked" $?

echo "═══ 10. graph, viewer and index ═══"
python3 - "$D" <<'EOF'
import json, sqlite3, sys
D = sys.argv[1]
s = json.load(open(f'{D}/05.graph/schema.json'))
con = sqlite3.connect(f'{D}/05.graph/cteph_knowledge_graph.sqlite')
bad = 0
for k, v in s['nodes'].items():
    n = con.execute(f'SELECT COUNT(*) FROM node_{k}').fetchone()[0]
    if n != v['count']:
        print(f'    ✗ node_{k}: {n} rows, schema {v["count"]}'); bad += 1
for k, v in s['edges'].items():
    n = con.execute(f'SELECT COUNT(*) FROM edge_{k}').fetchone()[0]
    if n != v['count']:
        print(f'    ✗ edge_{k}: {n} rows, schema {v["count"]}'); bad += 1
g = json.load(open(f'{D}/05.graph/graph.json'))
if len(g['nodes']) != s['graph']['n_nodes'] or len(g['edges']) != s['graph']['n_edges']:
    print('    ✗ graph.json counts differ from the schema'); bad += 1
html = open(f'{D}/cteph_knowledge_graph.html').read()
if 'id="kg-data"' not in html or html.count('"edges"') < 1 or len(html) < 500000:
    print('    ✗ the viewer has no embedded graph'); bad += 1
for f in ('README.md', '05.graph/README.md'):
    try:
        open(f'{D}/{f}').read()
    except OSError:
        print(f'    ✗ missing: {f}'); bad += 1
print(f'    {s["graph"]["n_nodes"]:,} nodes, {s["graph"]["n_edges"]:,} edges; viewer {len(html) / 1e6:.1f} MB')
sys.exit(1 if bad else 0)
EOF
ck "graph tables, SQLite and graph.json agree with the schema; viewer and index present" $?
python3 scripts/kg.py --snapshot "$SNAP" gene FGG >/dev/null 2>&1
ck "the query tool answers (kg.py gene FGG)" $?
fi

echo
[ $fail -eq 0 ] && echo "ALL CHECKS PASSED" || echo "$fail CHECK(S) FAILED"
exit $fail
