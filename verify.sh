#!/bin/bash
# Pre-flight and post-run checks of DiseaseLens. All read-only. Run from anywhere:
#   ./verify.sh                      every disease of diseases/ and every snapshot under results/
#   DISEASES=cteph,hpah ./verify.sh  these diseases only
#   NO_LLM=1 ./verify.sh             without the one model call of section 6
#   NO_NET=1 ./verify.sh             without network (no profile resolution against Open Targets)
#   RESULTS=<dir> ./verify.sh        the snapshots of another results folder (default: results)
#   DISEASES_DIR=<dir> ./verify.sh   the disease profiles of another folder (default: diseases), as for a run
#                                    with --diseases_dir
# Tools: nextflow on PATH (or NXF_BIN=<dir>); PYTHON=<python3 with requirements.txt> (default: python3).
# Sections 1-6 check the source tree, the configuration of every disease, privacy and the isolation of the
# language model, and must pass before a run. Sections 7-11 read the newest snapshot of each disease (SNAP=<date> for
# another) and are skipped for a disease without one.
C=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$C"
fail=0
ck() { if [ "$2" -eq 0 ]; then printf "  ✓ %s\n" "$1"; else printf "  ✗ %s\n" "$1"; fail=$((fail+1)); fi; }
[ -n "$NXF_BIN" ] && export PATH=$NXF_BIN:$PATH
PY=${PYTHON:-python3}
RES=${RESULTS:-results}
export PYTHONDONTWRITEBYTECODE=1
DDIR=${DISEASES_DIR:-diseases}
KEYS=${DISEASES:-$(for d in "$DDIR"/*/profile.yaml; do basename "$(dirname "$d")"; done | paste -sd, -)}
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

echo "═══ 1. syntax ═══"
out=$(nextflow lint main.nf nextflow.config 2>&1)
echo "$out" | grep -q '❌'; ck "nextflow lint: no errors" $([ $? -eq 1 ] && echo 0 || echo 1)
echo "$out" | grep -q '⚠️'; ck "nextflow lint: no warnings" $([ $? -eq 1 ] && echo 0 || echo 1)
bad=0; for f in scripts/*.py; do $PY -c "import ast,sys; ast.parse(open(sys.argv[1]).read())" "$f" 2>/dev/null || bad=1; done
ck "every Python script compiles ($(ls scripts/*.py | wc -l))" $bad
$PY -c "import xlsxwriter, openpyxl, pandas, yaml, matplotlib" 2>/dev/null
ck "$PY imports the packages of requirements.txt" $?

echo "═══ 2. pipeline rules ═══"
grep -qE '^\s*(cpus|memory)\s' main.nf nextflow.config; ck "no cpus/memory directive" $([ $? -eq 1 ] && echo 0 || echo 1)
grep -q "errorStrategy = 'finish'" nextflow.config; ck "errorStrategy = 'finish'" $?
$PY - <<'EOF'
import re, sys, pathlib
nf = pathlib.Path('main.nf').read_text()
bodies = dict(re.findall(r'process\s+(\w+)\s*\{(.*?)\n\}', nf, re.S))
bad = [p for p, b in bodies.items() if "executor 'local'" not in b]
for p in bad:
    print(f'    ✗ {p}: not local (compute nodes may have no internet)')
print(f'    {len(bodies)} processes checked')
sys.exit(1 if bad else 0)
EOF
ck "every process local" $?

echo "═══ 3. every process's flags are accepted by its script ═══"
$PY - <<'PY'
import re, subprocess, sys, pathlib
nf = pathlib.Path('main.nf').read_text()
bodies = dict(re.findall(r'process\s+(\w+)\s*\{(.*?)\n\}', nf, re.S))
PAIRS = {'PROFILE': 'profile.py', 'FETCH_CORE': 'fetch.py', 'SEED_PUBMED': 'pubmed.py', 'SEARCH_PUBMED': 'pubmed.py',
         'RELATIONS': 'relations.py', 'RESOLVE': 'resolve.py', 'FETCH_RELATED': 'fetch.py',
         'SCREEN': 'extract_genetics.py', 'FETCH_FULLTEXT': 'pubmed.py', 'EXTRACT': 'extract_genetics.py',
         'PUBLICATIONS': 'publications.py', 'CURATE': 'curate.py', 'TIDY_GRAPH': 'tidy_graph.py',
         'BUILD_GRAPH': 'build_graph.py', 'FIGURES': 'draw_figures.py', 'SUMMARISE': 'summarise.py',
         'DIFF': 'diff_snapshots.py', 'KNOWLEDGE': 'build_knowledge.py', 'WRITE_TABLE': 'write_table.py',
         'VIEWER': 'build_viewer.py', 'WRITE_RUN_MANIFEST': 'manifest.py', 'PROJECT_INDEX': 'manifest.py'}
bad = 0
for p, s in PAIRS.items():
    body = bodies[p].split('script:')[-1]
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

echo "═══ 4. shared rules and the configuration of each disease ═══"
$PY - <<'EOF'
import re, sys, yaml
bad = 0
s = yaml.safe_load(open('config/search.yaml'))
c = yaml.safe_load(open('config/curation.yaml'))
pats = s['prescreen']['genetic_terms'] + c['result_check']['negative'] + c['result_check']['positive'] + \
    [p for v in c['variants'] for p in v['patterns']]
for p in pats:
    try:
        re.compile(p)
    except re.error as e:
        print(f'    ✗ pattern {p}: {e}'); bad += 1
if not {'retracted', 'conference'} <= {r['key'] for r in c['publications']['exclude']}:
    print('    ✗ the exclusion rules lack retracted / conference'); bad += 1
print(f'    {len(pats)} patterns, {len(c["variants"])} variant names, {len(c["design_classes"])} design classes')
sys.exit(1 if bad else 0)
EOF
ck "shared rules: search blocks, pre-screen, variant and result-check patterns valid" $?
for k in ${KEYS//,/ }; do
  API=""; [ -z "$NO_NET" ] && API="--api https://api.platform.opentargets.org/api/v4/graphql"
  $PY scripts/profile.py --disease-dir "$DDIR/$k" --config-dir config $API --out $TMP/$k.json >/dev/null 2>$TMP/$k.err
  ck "$k: the profile resolves" $?; [ -s $TMP/$k.err ] && tail -2 $TMP/$k.err
  $PY - $TMP/$k.json "$k" <<'EOF'
import re, sys, json
sys.path.insert(0, 'scripts')
from ot_common import BASIS_TYPES, Profile, check_themes, load_choices
P = Profile(json.load(open(sys.argv[1])))
nf = open('main.nf').read()
block = re.search(r'def runParams = \[(.*?)\]', nf, re.S).group(1)
run_keys = set(re.findall(r'(\w+):', block)) | {'FetchDate', 'Rebuild', 'DiseaseCore', 'PubmedQuery'}
refs = P.choices.get('references') or {}
bad = 0
for g in P.choices['groups']:
    for it in g['items']:
        k = it['key']
        if not (it.get('basis') or '').strip() or '{short}' in str(it) or '{name}' in str(it):
            print(f'    ✗ {k}: no basis, or a name not filled in'); bad += 1
        if it.get('basis_type') not in BASIS_TYPES:
            print(f'    ✗ {k}: basis type {it.get("basis_type")}'); bad += 1
        if it.get('basis_type') == 'source' and not it.get('refs'):
            print(f'    ✗ {k}: a source basis without a reference'); bad += 1
        if 'param' in it and it['param'] not in run_keys:
            print(f'    ✗ {k}: parameter {it["param"]} is not recorded by the run'); bad += 1
        for r in it.get('refs') or []:
            if r not in refs:
                print(f'    ✗ {k}: reference {r} missing'); bad += 1
load_choices(P, {k: 0 for k in run_keys})
if P.themes:
    check_themes(P.themes)
    for t in P.themes:
        if not t.get('basis') or not t.get('refs') or any(r not in refs for r in t['refs']):
            print(f'    ✗ gene theme {t["key"]}: basis or references missing'); bad += 1
for c in (P.seed or {}).get('conditions') or []:
    ev = (c.get('evidence') or []) + (c.get('east_asian') or [])
    if not ev or any(not (str(e.get('pmid', '')).isdigit() and e.get('quote') and e.get('statement')) for e in ev):
        print(f'    ✗ {c["key"]}: an evidence item lacks pmid / quote / statement'); bad += 1
    if c.get('genetics') not in ('full', 'coloc_only', 'none'):
        print(f'    ✗ {c["key"]}: genetics mode {c.get("genetics")}'); bad += 1
codes = (P.gwas or {}).get('codes') or {}
for r in (P.gwas or {}).get('regions') or []:
    if r.get('evidence') not in codes or not r.get('rsid') or not r.get('genes'):
        print(f'    ✗ published region {r.get("gene")}: evidence code, rsid or genes missing'); bad += 1
want = P.checks.get('published_regions')
if want is not None and len((P.gwas or {}).get('regions') or []) != want:
    print(f'    ✗ {want} published regions expected'); bad += 1
cur = P.curation
for k in ('within_disease', 'other_trait', 'functional', 'background'):
    for it in cur['findings'][k]['items']:
        if not (str(it.get('pmid', '')).isdigit() and it.get('match')):
            print(f'    ✗ finding item without pmid / match: {it}'); bad += 1
for it in cur['text_mining']:
    if not (it.get('gene') and it.get('words') and it.get('meaning')):
        print(f'    ✗ text-mining item: {it}'); bad += 1
print(f'    {P.id} {P.name} ({P.short}): relations {P.relations["mode"]}, '
      f'{len((P.seed or {}).get("conditions") or [])} curated conditions, '
      f'{len((P.gwas or {}).get("regions") or [])} published regions, themes {"curated" if P.themes else "Reactome"}, '
      f'{sum(len(g["items"]) for g in P.choices["groups"])} choices, {len(refs)} references')
sys.exit(1 if bad else 0)
EOF
  ck "$k: choices with basis and references; seed, regions, themes and curation entries valid" $?
done
for f in README.md docs/METHODS.md docs/USAGE.md conf/local.config.example; do
  [ -s "$f" ]; ck "$f present" $?
done

echo "═══ 5. privacy of the repository ═══"
if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  me=$(whoami); grp=$(id -gn); parent=$(dirname "$C")
  hits=$(git ls-files -z | xargs -0 grep -IlE "(/home/|/Users/|$parent|\b$me\b|\b$grp\b|sk-ant-[A-Za-z0-9]|sk-[A-Za-z0-9]{32,}|(api_key|apikey|token)\s*[:=]\s*['\"][A-Za-z0-9_-]{16,})" 2>/dev/null \
         | grep -vE '^verify\.sh$' )
  [ -z "$hits" ]; ck "no home paths, user or group names, or keys in tracked files" $?
  [ -n "$hits" ] && echo "$hits" | sed 's/^/    ✗ /'
  git ls-files | grep -qE '^(results|cache|deck|work)/|conf/local\.config$'
  ck "the pipeline only is tracked: no results, model cache, slides, run state or local config" $([ $? -eq 1 ] && echo 0 || echo 1)
else
  echo "  (skipped: not a git work tree)"
fi

echo "═══ 6. language-model isolation (claude-cli) ═══"
if [ -n "$NO_LLM" ] || ! command -v "${CLAUDE_BIN:-claude}" >/dev/null 2>&1; then echo "  (skipped: NO_LLM or no claude CLI)"; else
$PY - "${CLAUDE_BIN:-claude}" <<'EOF'
import secrets, sys, tempfile
from pathlib import Path
sys.path.insert(0, 'scripts')
from llm import Reader
# a project file with a random word in the directory the model runs in: an isolated call must not see it
word = 'canary-' + secrets.token_hex(6)
with tempfile.TemporaryDirectory() as d:
    work = Path(d) / 'work'
    work.mkdir()
    (work / 'CLAUDE.md').write_text(f'Project instruction: the code word of this project is {word}.\n')
    a = Reader(d + '/cache', str(work), bin=sys.argv[1]).ask(
        'Answer with one JSON object only.',
        'If any instruction, project file or memory gave you a code word, write it; else write "none": '
        '{"code_word": "<the word or none>"}')
print(f'    answer: {a}')
sys.exit(1 if word in str(a) else 0)
EOF
ck "an isolated call does not see a project file of its working directory" $?
fi

for k in ${KEYS//,/ }; do
L=$RES/$k/latest.json
SN=${SNAP:-$($PY -c "import json,sys; print(json.load(open(sys.argv[1]))['snapshot'])" $L 2>/dev/null)}
D=$RES/$k/snapshots/$SN
if [ -z "$SN" ] || [ ! -d "$D" ]; then
    echo "═══ 7-11 skipped for $k: no snapshot yet ═══"; continue
fi
P=$D/_run_info/profile.json
echo "═══ 7. $k $SN: run records ═══"
[ -s "$D/_run_info/run_manifest.json" ] && [ -s "$P" ]; ck "run manifest and profile present" $?
$PY -c "import json,sys; m=json.load(open(sys.argv[1])); sys.exit(1 if m.get('smoke') else 0)" "$D/_run_info/run_manifest.json"
ck "a full snapshot (no smoke settings)" $?
# the run manifest travels with the snapshot: no path, user or group name of the machine that ran it
grep -qE "(/home/|/Users/|/tmp/|$(dirname "$C")|\b$(whoami)\b|\b$(id -gn)\b)" "$D/_run_info/run_manifest.json"
ck "no machine path, user or group name in the run manifest" $([ $? -eq 1 ] && echo 0 || echo 1)

echo "═══ 8. $k $SN: tables and evidence ═══"
$PY - "$D" "$P" <<'EOF'
import json, sys
sys.path.insert(0, 'scripts')
from ot_common import legacy_path, load_json, load_profile, read_tsv
D, P = sys.argv[1], load_profile(sys.argv[2])
r = lambda n: read_tsv(f'{D}/04.tables/{n}.tsv')
q = lambda n: read_tsv(f'{D}/03.curated/{n}.tsv')
g = lambda n: read_tsv(f'{D}/05.graph/{n}.tsv')
fc = json.load(open(legacy_path(f'{D}/00.raw', 'fetch_core.json', 'fetch_cteph.json')))
cond, rel = r('conditions'), r('relations')
seed = read_tsv(f'{D}/01.conditions/seed_check.tsv')
st, fi, rg, lad, men = r('disease_studies'), r('disease_findings'), r('disease_regions'), r('disease_genes'), r('mendelian')
loci = r('gwas_loci')
cs_of = {x for v in rg.disease_study_locus_id for x in v.split() if x}
checks = {
    # the API count, less the credible sets of GWAS studies whose publication the curation left out
    'credible sets of the disease = the API count': fc['credible_sets'] == int((loci.node_id == P.id).sum()) + sum(
        len(s_['credible_sets']) for s_ in load_json(f'{D}/00.raw/gwas.json.gz')['studies']
        if str(s_.get('pubmedId') or '') in set(q('publications').query('status == "excluded"').pmid)),
    'every seed quote found in its abstract': bool(seed.quote_in_abstract.isin(['True']).all()),
    'every condition has evidence': set(cond.node_id) <= set(rel.node_id),
    'every relation has a type and a source': bool((rel.evidence_type.ne('') & (rel.pmid.ne('') | rel.citation.ne(''))).all()),
    'every finding belongs to an included study': set(fi.pmid) <= set(st.pmid),
    'graph edges reference graph nodes': set(g('edges').source) | set(g('edges').target) <= set(g('nodes').id),
    'ClinVar review status kept apart from classifications': 'review_status' in men.columns and
        not men.query('source == "ClinVar"').classifications.ne('').any(),
    'lead-variant overlap only at regions without a credible set of the disease':
        set(r('overlap').disease_gene) <= set(rg[rg.disease_study_locus_id.eq('')].gene),
    'colocalization only at regions with a credible set of the disease':
        set(r('colocalisation').disease_gene) <= set(rg[rg.disease_study_locus_id.ne('')].gene),
    'code G and G* genes are exactly the genes of the regions':
        set(lad[lad.code.isin(['G', 'G*'])].symbol) <= {x for v in rg.genes for x in v.split()},
    'code C genes have curated evidence for the disease':
        set(lad[lad.code.eq('C')].symbol) <= set(men[men.node_id == P.id].symbol),
    'every included study and finding has its readers\' agreement': bool(
        st.agreement.isin(['both', 'adjudicated', 'one reader']).all() and
        fi.agreement.isin(['both', 'adjudicated', 'one reader']).all()) if 'agreement' in st else False,
}
if P.checks.get('published_regions') is not None:
    checks['the published regions present'] = int(rg.source.eq('published').sum()) == P.checks['published_regions']
sj = json.load(open(f'{D}/06.summaries/summaries.json'))
s = sj['summaries']
made = [x for x in s if not x['problems'].startswith('not generated')]
# a summary the checks did not pass is withheld, with its reason (the share withheld is reported, not judged)
checks['summaries: every summary withheld says why'] = all(x['grounded'] or str(x['problems']).strip() for x in s)
checks['summaries: an offline build calls no model'] = not sj.get('offline') or sj['llm_calls'] == 0
pubs = q('publications')
ex = set(pubs[pubs.status == 'excluded'].pmid)
pm_graph = {x.split(':', 1)[1] for x in g('nodes').query('type == "publication"').id}
checks['no excluded publication in the studies, literature, relations or graph'] = not (
    ex & (set(st.pmid) | set(r('literature').pmid) | set(rel.pmid) | pm_graph))
checks['every finding says what it is about'] = bool(fi.association_with.ne('').all())
pos = fi[fi.association_with.eq(P.about) & fi.result.isin(['associated', 'enriched'])]
checks[f'ladder levels 3-4 rest on findings about {P.short} itself'] = all(
    (pos[pos.symbol == x.symbol].pmid.nunique() >= 2) if x.level == '3' else pos[pos.symbol == x.symbol].pmid.nunique() == 1
    for x in lad[lad.level.isin(['3', '4'])].itertuples())
checks['every study has a design class'] = bool(st.design_class.ne('').all())
for k, v in checks.items():
    print(f'    {"✓" if v else "✗"} {k}')
print(f'    conditions {len(cond)}, regions {len(rg)}, genes on the ladder {len(lad)}, genetic studies {len(st)}, '
      f'findings {len(fi)}, summaries {len(s)} ({sum(x["grounded"] for x in s)} grounded, '
      f'{sum(not x["grounded"] for x in made)} withheld after the checks)')
sys.exit(0 if all(checks.values()) else 1)
EOF
ck "tables, evidence, graph and summaries consistent" $?
$PY - "$D" "$P" <<'EOF'
import json, sys
sys.path.insert(0, 'scripts')
from ot_common import load_choices, load_profile
D = sys.argv[1]
m = json.load(open(f'{D}/_run_info/run_manifest.json'))
groups, _ = load_choices(load_profile(sys.argv[2]), m['params'], f'{D}/04.tables')
items = [it for g in groups for it in g['items'] if it['check_text']]
for it in items:
    print(f'    {"✓" if it["check_ok"] else "✗"} {it["item"]}: {it["check_text"]}')
sys.exit(0 if all(it['check_ok'] for it in items) else 1)
EOF
ck "choices: the data checks hold on this snapshot" $?

echo "═══ 9. $k $SN: figures and knowledge base ═══"
K=$D/08.knowledge/${k}_knowledge_graph.md
bad=0
for s in "## 1. " "## 2. " "## 3. " "## 4. " "## 5. " "## 10. " "## 11. " "## 12. " "## 13. " "## References"; do
  grep -q "^$s" "$K" || { echo "    ✗ missing: $s"; bad=1; }
done
ck "the knowledge base has its sections" $bad
bad=0
for f in $($PY -c "import json,sys; print(' '.join(json.load(open(sys.argv[1])).get('figures', {}).keys()))" $D/09.figures/figures.json 2>/dev/null); do
  for e in png pdf; do [ -s "$D/09.figures/$f.$e" ] || { echo "    ✗ missing: 09.figures/$f.$e"; bad=1; }; done
done
[ -s "$D/09.figures/README.md" ] || { echo "    ✗ missing: 09.figures/README.md"; bad=1; }
for f in $(grep -oE '\.\./09\.figures/[a-z_]+\.png' "$K" | sort -u); do
  [ -s "$D/08.knowledge/$f" ] || { echo "    ✗ the knowledge base links a missing figure: $f"; bad=1; }
done
ck "figures (PNG and PDF), their catalog, and every figure the knowledge base links" $bad
grep -qwE 'nan|NaN' "$K"; ck "no nan in the knowledge base" $([ $? -eq 1 ] && echo 0 || echo 1)

echo "═══ 10. $k $SN: workbooks ═══"
for X in $D/${k}_knowledge_graph_$SN.xlsx $D/${k}_knowledge_graph_${SN}_summary.xlsx; do
[ -s "$X" ]; ck "$(basename "$X") present" $?
[ -s "$X" ] || continue
$PY - "$X" <<'EOF'
import sys, warnings
sys.path.insert(0, 'scripts')
from write_table import FORBIDDEN_TOKENS, plain
from openpyxl import load_workbook
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    wb = load_workbook(sys.argv[1], rich_text=True)
bad = 0 if wb.sheetnames[0] == 'Guide' else 1
names = set(wb.sheetnames)
for ws in wb.worksheets:
    for row in ws.iter_rows():
        for c in row:
            if c.value is None:
                continue
            s = plain(c.value) + plain(c.value, sup=True)
            hit = [t for t in FORBIDDEN_TOKENS if t in s]
            if hit or s in ('nan', 'None'):
                print(f'    ✗ {ws.title}!{c.coordinate}: {hit or s}'); bad += 1
            loc = c.hyperlink.location if c.hyperlink is not None else ''
            if loc and loc.split('!')[0].strip("'") not in names:
                print(f'    ✗ {ws.title}!{c.coordinate}: a link to a sheet that is not there ({loc})'); bad += 1
links = sum(1 for ws in wb.worksheets for row in ws.iter_rows() for c in row
            if c.hyperlink is not None and c.hyperlink.location)
print(f'    {len(wb.sheetnames)} sheets, {links} links between sheets')
sys.exit(1 if bad or not links else 0)
EOF
ck "$(basename "$X"): sheets linked, every link to a sheet that exists; no internal names or missing values" $?
done

echo "═══ 11. $k $SN: graph, viewer and index ═══"
$PY - "$D" "$k" <<'EOF'
import glob, json, sqlite3, sys
D, k = sys.argv[1], sys.argv[2]
s = json.load(open(f'{D}/05.graph/schema.json'))
con = sqlite3.connect(glob.glob(f'{D}/05.graph/*.sqlite')[0])
bad = 0
for kind, prefix in (('nodes', 'node'), ('edges', 'edge')):
    for t, v in s[kind].items():
        n = con.execute(f'SELECT COUNT(*) FROM {prefix}_{t}').fetchone()[0]
        if n != v['count']:
            print(f'    ✗ {prefix}_{t}: {n} rows, schema {v["count"]}'); bad += 1
g = json.load(open(f'{D}/05.graph/graph.json'))
if len(g['nodes']) != s['graph']['n_nodes'] or len(g['edges']) != s['graph']['n_edges']:
    print('    ✗ graph.json counts differ from the schema'); bad += 1
html = open(f'{D}/{k}_knowledge_graph.html').read()
if 'id="kg-data"' not in html:
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
EX=$($PY -c "import json,sys; d=json.load(open(sys.argv[1])); print(d.get('checks',{}).get('example_gene',''))" $P)
[ -n "$EX" ] && { $PY scripts/kg.py --root "$RES" --disease $k gene "$EX" >/dev/null 2>&1; ck "the query tool answers (kg.py --disease $k gene $EX)" $?; }
done
[ -s $RES/index.json ] && [ -s $RES/index.html ]; ck "the project index ($RES/index.json, index.html) present" $?

echo
[ $fail -eq 0 ] && echo "ALL CHECKS PASSED" || echo "$fail CHECK(S) FAILED"
exit $fail
