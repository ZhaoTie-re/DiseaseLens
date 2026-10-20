nextflow.enable.dsl = 2

// ═════════════════════════════════════════════════════════════════════════════
// DiseaseLens — for each disease of a run: the disease, the conditions related
// to it, and the genetics of both, as a dated, traceable knowledge graph (one
// snapshot per disease and fetch date, each compared with the snapshot before
// it). A disease is a folder diseases/<key>/ whose profile.yaml names its Open
// Targets ID; everything else has a default (README, "Add a disease").
//
// External knowledge only: Open Targets (GraphQL), PubMed (E-utilities) and
// Europe PMC. One line of stages per disease:
//   profile    the profile resolved (profile.py): names, queries, merged rules
//   sources    00  Open Targets (the disease; each related condition), PubMed (the
//                  genetic-study search; the seed or the relation search), Europe
//                  PMC Methods texts, and the PubMed record of every publication
//   evidence   01  related conditions and the evidence for each relation (a curated
//                  list, or relations read from PubMed by two model readers)
//              02  the genetic studies of the disease, read by two independent
//                  model readers, disagreements adjudicated, every item checked
//                  against its source text
//   curation   03  publications left out or flagged; variants named; findings
//                  classified; text-mined links mapped to the wrong gene removed
//   tables     04  the long tables, the GWAS regions and the evidence ladder
//   graph      05  typed nodes and edges, graph.json, a SQLite database, the schema
//   products   06  summaries · 07 changes · 08 knowledge base · 09 figures ·
//              the workbook · the viewer · the index; results/index.json
//
//   PROFILE ─► FETCH_CORE ─┬─► RESOLVE ─► FETCH_RELATED ────────────────┐
//              SEED_PUBMED ┤  (curated list)                            │
//   SEARCH_RELATIONS ─► RELATIONS ┘  (no curated list)                  │
//   SEARCH_PUBMED ─► SCREEN ─► FETCH_FULLTEXT ─► EXTRACT ───────────────┤
//                                   PUBLICATIONS ◄──────────────────────┤
//                                        └──► CURATE ─► TIDY_GRAPH ─► BUILD_GRAPH
//      BUILD_GRAPH ─┬─► FIGURES, SUMMARISE, DIFF ─► KNOWLEDGE ─► WRITE_TABLE ─┐
//                   └─► VIEWER ◄────────────────────────────────────────────┤
//                                  WRITE_RUN_MANIFEST ─► PROJECT_INDEX ◄──┘
//
// A disease marked `rebuild` re-derives its snapshot of `fetch_date` from the
// saved sources (00.raw, 00.pubmed, 01.conditions): nothing is fetched again
// except PubMed publication records not yet saved and Open Targets gene symbols
// (in the snapshot's own release); the readers answer from the cache, and only
// answers not yet cached are asked (none with --LlmOffline true).
// Every process is local: compute nodes may have no internet.
// Methodology: docs/METHODS.md. Outputs: docs/OUTPUTS.md.
// ═════════════════════════════════════════════════════════════════════════════

// ---- Helpers
def script_file(String name) {
    return file("${params.script_dir}/${name}", checkIfExists: true)
}

def toJson(x) {
    return groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(x))
}

def flag(v) {
    return v instanceof Boolean ? v : ['true', '1', 'yes'].contains(v.toString().toLowerCase())
}

// the diseases of the run: `diseases` as a list (params file: keys, or maps with key, fetch_date, rebuild) or as a
// comma-separated string (command line); `FetchDate` is the date of every disease that names none
def diseaseSpecs() {
    def raw = params.diseases instanceof List ? params.diseases : params.diseases.toString().tokenize(',')
    def rebuildAll = params.Rebuild instanceof Boolean ? params.Rebuild : false
    def rebuildKeys = params.Rebuild instanceof Boolean ? [] : params.Rebuild.toString().tokenize(',').collect { k -> k.trim() }
    def specs = raw.collect { x ->
        def m = x instanceof Map ? x : [key: x.toString().trim()]
        def key = diseaseKey(m.key.toString())
        def date = (m.fetch_date ?: params.FetchDate ?: '').toString()
        if (!(date ==~ /\d{4}-\d{2}-\d{2}(-[a-z0-9]+)?/)) {
            error "DiseaseLens | ${key}: a fetch date YYYY-MM-DD is required (--FetchDate or fetch_date), got '${date}'"
        }
        def dir = file("${params.diseases_dir}/${key}")
        if (!file("${dir}/profile.yaml").exists()) {
            error "DiseaseLens | ${key}: no ${dir}/profile.yaml (README, 'Add a disease')"
        }
        def root = "${params.results_root}/${key}/snapshots"
        def earlier = file(root).exists()
            ? file(root).listFiles().findAll { d -> d.isDirectory() && d.name != date && file("${d}/_run_info/run_manifest.json").exists() }
                .collect { d -> d.name }.sort()
            : []
        // the profile's files in the task hash: a change to any of them re-runs the disease
        def own = []
        dir.eachFileRecurse(groovy.io.FileType.FILES) { f -> own << f }
        own = own.sort { f -> f.toString() }.collect { f -> f.text.md5() }
        def shared = ['curation.yaml', 'choices.yaml', 'search.yaml'].collect { n -> file("${params.config_dir}/${n}").text.md5() }
        return [key: key, fetch_date: date, rebuild: flag(m.rebuild ?: false) || rebuildAll || rebuildKeys.contains(key),
                curated: file("${dir}/related_conditions.yaml").exists(), disease_dir: dir.toString(),
                out_dir: "${root}/${date}".toString(), snapshots_root: root.toString(),
                earlier: earlier.join(','), workbook: "${key}_knowledge_graph_${date}.xlsx".toString(),
                profile_md5: (own + shared).join('')]
    }
    if (!specs) {
        error "DiseaseLens | no disease: give --diseases <key>[,<key>] or a params file with `diseases`"
    }
    return specs
}

// a disease named by its folder (diseases/<key>), or by its Open Targets ID: the folder whose profile has that ID
def diseaseKey(String given) {
    if (!(given ==~ /^[A-Za-z]+_\d+$/)) {
        return given
    }
    def hits = file(params.diseases_dir).listFiles().findAll { d -> file("${d}/profile.yaml").exists() }
        .findAll { d -> file("${d}/profile.yaml").readLines().any { line -> line.startsWith('id:') && line.substring(3).tokenize(' #')[0] == given } }
    if (hits.size() != 1) {
        error "DiseaseLens | ${given}: ${hits ? 'several disease folders have this ID' : 'no disease folder has this ID yet'}; " +
            "write one with: python3 scripts/new_disease.py ${given}"
    }
    return hits[0].name
}

// the model configuration written for the steps that read with models (no key: keys stay in the environment)
def llmJson() {
    def L = params.llm
    def common = [bin: L.bin ?: 'claude']
    return toJson([readers: L.readers.collect { r -> common + r }, adjudicator: L.adjudicator ? common + L.adjudicator : null,
                   cache_dir: params.llm_cache_dir, workdir: params.llm_workdir, parallel: L.parallel ?: 3,
                   timeout: L.timeout ?: 600, max_tokens: L.max_tokens ?: 16000])
}

def readerText(r) {
    return "${r.provider ?: 'claude-cli'} ${r.model ?: '(default model)'}${r.effort ? ', effort ' + r.effort : ''}"
}

def readersText() {
    def L = params.llm
    def txt = L.readers.withIndex().collect { r, i -> "${r.name ?: ['A', 'B'][i]}: ${readerText(r)}" }.join('; ')
    return L.adjudicator && L.readers.size() == 2 ? "${txt}; adjudicator: ${readerText(L.adjudicator)}" : txt
}

// a channel of [spec, ...] keyed by the disease key, for join()
def keyed(ch) {
    return ch.map { row -> [row[0].key] + row.drop(1) }
}


// ═════════════════════════════════════════════════════════════════════════════
// PROCESSES — each takes the spec of its disease (d) first
// ═════════════════════════════════════════════════════════════════════════════

/* PROFILE — the disease's profile resolved and checked (profile.py). */
process PROFILE {
    executor 'local'
    tag "${d.key}"
    publishDir path: { "${d.out_dir}/_run_info" }, mode: 'copy'

    input:
    tuple val(d), path(disease_dir), path(config_dir)
    path script
    path lib

    output:
    tuple val(d), path('profile.json'), emit: json

    script:
    def saved = d.rebuild ? "--saved-disease '${d.out_dir}/00.raw/disease.json' --saved-search '${d.out_dir}/00.pubmed'" :
        "--api '${params.ot_api}'"
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    # profile files: ${d.profile_md5}
    ${params.python} ${script} --disease-dir ${disease_dir} --config-dir ${config_dir} ${saved} --out profile.json
    """
}

/* FETCH_CORE — the disease: targets, text-mined and GWAS evidence, target pathways, its GWAS credible sets with their
 * colocalizations, its curated gene evidence, and its published GWAS regions (if the profile has some). */
process FETCH_CORE {
    executor 'local'
    tag "${d.key}"
    publishDir path: { "${d.out_dir}" }, mode: 'copy'

    input:
    tuple val(d), path(profile)
    path script
    path lib

    output:
    tuple val(d), path('00.raw'), emit: raw

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --mode core --api '${params.ot_api}' --profile ${profile} --fetch-date ${d.fetch_date} \\
        --threshold ${params.MendelianMin} --page-size ${params.PageSize} --evidence-page ${params.EvidencePage} \\
        --out-dir 00.raw
    """
}

/* SEED_PUBMED — the PubMed records of every PMID the curated list cites. */
process SEED_PUBMED {
    executor 'local'
    tag "${d.key}"
    publishDir path: { "${d.out_dir}/00.pubmed" }, mode: 'copy'

    input:
    tuple val(d), path(profile)
    path script
    path lib

    output:
    tuple val(d), path('seed_pubmed.json.gz'), emit: records

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    # fetch date: ${d.fetch_date}
    ${params.python} ${script} --mode seed --profile ${profile} --out seed_pubmed.json.gz
    """
}

/* SEARCH_PUBMED — the PubMed records a query of the profile finds, up to the fetch date: the genetic studies
 * (`pubmed_query`) or the clinical relations (`relations`). */
process SEARCH_PUBMED {
    executor 'local'
    tag "${d.key} ${what}"
    publishDir path: { "${d.out_dir}/00.pubmed" }, mode: 'copy'

    input:
    tuple val(d), path(profile), val(what)
    path script
    path lib

    output:
    tuple val(d), path("${what == 'relations' ? 'relations_search' : 'pubmed_search'}.json.gz"), emit: records

    script:
    def out = what == 'relations' ? 'relations_search' : 'pubmed_search'
    def key = what == 'relations' ? "['relations']['query']" : "['pubmed_query']"
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    Q=\$(${params.python} -c "import json,sys; print(json.load(open(sys.argv[1]))${key})" ${profile})
    ${params.python} ${script} --mode search --query "\$Q" --max-date ${d.fetch_date.take(10)} --out ${out}.json.gz
    """
}

/* RELATIONS — the clinical relations of a disease without a curated list, read from PubMed by the model readers. */
process RELATIONS {
    executor 'local'
    tag "${d.key}"
    publishDir path: { "${d.out_dir}" }, mode: 'copy', pattern: '01.relations'
    publishDir path: { "${d.out_dir}/00.pubmed" }, mode: 'copy', pattern: 'seed_pubmed.json.gz'

    input:
    tuple val(d), path(profile), path(records)
    val llm
    val offline
    path script
    path lib
    path llmlib

    output:
    tuple val(d), path('01.relations'), path('seed_pubmed.json.gz'), emit: seed
    path '01.relations'

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    cat > llm.json <<'JSON'
${llm}
JSON
    ${params.python} ${script} --records ${records} --profile ${profile} --llm llm.json --offline ${offline} \\
        --batch ${params.LlmBatch} --limit ${params.LlmLimit} --api '${params.ot_api}' --out-dir 01.relations
    cp 01.relations/relations_pubmed.json.gz seed_pubmed.json.gz
    """
}

/* RESOLVE — the related conditions and their evidence (clinical quotes checked, ontology, colocalization). */
process RESOLVE {
    executor 'local'
    tag "${d.key}"
    publishDir path: { "${d.out_dir}" }, mode: 'copy'

    input:
    tuple val(d), path(profile), path(seed), path(seed_pubmed), path(raw)
    path script
    path lib
    path llmlib

    output:
    tuple val(d), path('01.conditions'), emit: dir
    tuple val(d), path('01.conditions/fetch_list.tsv'), emit: fetch_list

    script:
    def seedArg = d.curated ? '' : "--seed ${seed}"
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    PUB=""; [ -e ${raw}/published_loci.json.gz ] && PUB="--published ${raw}/published_loci.json.gz"
    ${params.python} ${script} --profile ${profile} ${seedArg} --seed-pubmed ${seed_pubmed} \\
        --core-gwas ${raw}/gwas.json.gz \$PUB --api '${params.ot_api}' --clpp-min ${params.ClppMin} \\
        --h4-min ${params.H4Min} --l2g-min ${params.L2gMin} --max-gwas ${params.MaxGwasStudies} \\
        --max-children ${params.MaxOntologyChildren} --region-window-kb ${params.RegionWindowKb} \\
        --genome-wide-p ${params.GenomeWideP} --region-min-maf ${params.RegionMinMaf} \\
        --region-min-mac ${params.RegionMinMac} --rare-variant-p ${params.RareVariantP} --out-dir 01.conditions
    """
}

/* FETCH_RELATED — every related condition's GWAS and curated / Mendelian gene evidence. */
process FETCH_RELATED {
    executor 'local'
    tag "${d.key}"
    publishDir path: { "${d.out_dir}/00.raw" }, mode: 'copy'

    input:
    tuple val(d), path(fetch_list)
    path script
    path lib

    output:
    tuple val(d), path('diseases'), emit: dir

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    # fetch date: ${d.fetch_date}
    mkdir -p diseases
    ${params.python} ${script} --mode related --api '${params.ot_api}' --list ${fetch_list} \\
        --threshold ${params.MendelianMin} --page-size ${params.PageSize} --parallel 3 --out-dir diseases
    """
}

/* SCREEN — pass A of the extraction (abstracts, both readers, adjudicated); lists the included records with a
 * PMCID for the full texts. */
process SCREEN {
    executor 'local'
    tag "${d.key}"

    input:
    tuple val(d), path(profile), path(records)
    val llm
    val offline
    path script
    path lib
    path llmlib

    output:
    tuple val(d), path('fulltext_ids.txt'), emit: ids

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    cat > llm.json <<'JSON'
${llm}
JSON
    ${params.python} ${script} --records ${records} --profile ${profile} --llm llm.json --offline ${offline} \\
        --api '${params.ot_api}' --batch ${params.LlmBatch} --limit ${params.LlmLimit} --prescreen ${params.Prescreen} \\
        --fulltext-list fulltext_ids.txt --out-dir screen
    """
}

/* FETCH_FULLTEXT — the Methods sections of the included open-access full texts (Europe PMC). */
process FETCH_FULLTEXT {
    executor 'local'
    tag "${d.key}"
    publishDir path: { "${d.out_dir}/00.pubmed" }, mode: 'copy'

    input:
    tuple val(d), path(ids)
    path script
    path lib

    output:
    tuple val(d), path('fulltext.json.gz'), emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --mode fulltext --ids ${ids} --out fulltext.json.gz
    """
}

/* EXTRACT — passes A (cached) and B: the genetic studies of the disease, their findings and methods. */
process EXTRACT {
    executor 'local'
    tag "${d.key}"
    publishDir path: { "${d.out_dir}" }, mode: 'copy'

    input:
    tuple val(d), path(profile), path(records), path(fulltext), val(release), path(raw)
    val llm
    val offline
    path script
    path lib
    path llmlib

    output:
    tuple val(d), path('02.literature'), emit: dir

    script:
    def expect = d.rebuild ? "--expect-release '${release}'" : ''
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    cat > llm.json <<'JSON'
${llm}
JSON
    ${params.python} ${script} --records ${records} --fulltext ${fulltext} --profile ${profile} --llm llm.json \\
        --offline ${offline} --api '${params.ot_api}' ${expect} --raw ${raw} --batch ${params.LlmBatch} --limit ${params.LlmLimit} \\
        --prescreen ${params.Prescreen} --out-dir 02.literature
    """
}

/* PUBLICATIONS — the PubMed record of every publication the graph uses (status, issue, affiliation). */
process PUBLICATIONS {
    executor 'local'
    tag "${d.key}"
    publishDir path: { "${d.out_dir}/00.pubmed" }, mode: 'copy'

    input:
    tuple val(d), path(search), path(seed_records), path(raw), path(diseases), path(conditions)
    path script
    path lib

    output:
    tuple val(d), path('publications.json.gz'), emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --pubmed ${search} ${seed_records} --raw ${raw} --diseases \$(ls ${diseases}/disease.*.json.gz 2>/dev/null || true) \\
        --relations ${conditions}/relations.tsv --out publications.json.gz
    """
}

/* CURATE — what the graph holds: exclusions, variant names, what each finding is about, text-mining review. */
process CURATE {
    executor 'local'
    tag "${d.key}"
    publishDir path: { "${d.out_dir}" }, mode: 'copy'

    input:
    tuple val(d), path(profile), path(publications), path(literature), path(raw), path(conditions)
    path script
    path lib

    output:
    tuple val(d), path('03.curated'), emit: dir

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --publications ${publications} --literature ${literature} --raw ${raw} \\
        --relations ${conditions}/relations.tsv --profile ${profile} --sentences ${params.SentencesPerPaper} \\
        --out-dir 03.curated
    """
}

/* TIDY_GRAPH — the long tables, the GWAS regions and the evidence ladder. */
process TIDY_GRAPH {
    executor 'local'
    tag "${d.key}"
    publishDir path: { "${d.out_dir}" }, mode: 'copy'

    input:
    tuple val(d), path(profile), path(raw), path(conditions), path(diseases), path(literature), path(curated)
    path script
    path lib

    output:
    tuple val(d), path('04.tables'), emit: tables

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --profile ${profile} --raw ${raw} --conditions ${conditions} \\
        --diseases \$(ls ${diseases}/disease.*.json.gz 2>/dev/null || true) --extract ${literature} --curated ${curated} \\
        --eas-min ${params.EasMin} --l2g-min ${params.L2gMin} --region-window-kb ${params.RegionWindowKb} \\
        --genome-wide-p ${params.GenomeWideP} --region-min-maf ${params.RegionMinMaf} \\
        --region-min-mac ${params.RegionMinMac} --rare-variant-p ${params.RareVariantP} --out-dir 04.tables
    """
}

/* BUILD_GRAPH — the typed graph: node and edge tables, graph.json, the SQLite database, the schema. */
process BUILD_GRAPH {
    executor 'local'
    tag "${d.key}"
    publishDir path: { "${d.out_dir}" }, mode: 'copy'

    input:
    tuple val(d), path(profile), path(curated), path(tables)
    path script
    path lib

    output:
    tuple val(d), path('05.graph'), emit: dir

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --profile ${profile} --curated ${curated} --tables ${tables} --out-dir 05.graph
    """
}

/* FIGURES — the snapshot's figures in the shared figure style (scripts/plot_style.py), with one catalog README. */
process FIGURES {
    executor 'local'
    tag "${d.key}"
    publishDir path: { "${d.out_dir}" }, mode: 'copy'

    input:
    tuple val(d), path(profile), path(tables)
    path style
    path script
    path lib

    output:
    tuple val(d), path('09.figures'), emit: dir

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=\$PWD/.mpl
    ${params.python} ${script} --profile ${profile} --tables ${tables} --plot-style ${style} \\
        --l2g-min ${params.L2gMin} --curated-min ${params.MendelianMin} --clpp-min ${params.ClppMin} \\
        --h4-min ${params.H4Min} --out-dir 09.figures
    """
}

/* SUMMARISE — grounded summaries: mechanism themes, methods and findings of the genetic studies, each checked by the
 * second reader (offline: cache only). */
process SUMMARISE {
    executor 'local'
    tag "${d.key}"
    publishDir path: { "${d.out_dir}/06.summaries" }, mode: 'copy'

    input:
    tuple val(d), path(profile), path(tables)
    val llm
    val offline
    path script
    path lib
    path llmlib

    output:
    tuple val(d), path('summaries.tsv'), path('summaries.json'), emit: out

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    cat > llm.json <<'JSON'
${llm}
JSON
    ${params.python} ${script} --tables ${tables} --profile ${profile} --llm llm.json \\
        --min-targets ${params.MinThemeTargets} --top-targets ${params.ThemeTargets} --max-themes ${params.MaxThemes} \\
        --offline ${offline} --out summaries.tsv --out-json summaries.json
    """
}

/* DIFF — the changes since the newest earlier snapshot of the disease. */
process DIFF {
    executor 'local'
    tag "${d.key}"
    publishDir path: { "${d.out_dir}/07.changes" }, mode: 'copy'

    input:
    tuple val(d), path(tables)
    path script
    path lib

    output:
    tuple val(d), path('changes.tsv'), path('changes.json'), emit: out

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --tables ${tables} --current ${d.fetch_date} --snapshots-root '${d.snapshots_root}' \\
        --earlier '${d.earlier}' --out changes.tsv --out-json changes.json
    """
}

/* KNOWLEDGE — the knowledge base of the snapshot. */
process KNOWLEDGE {
    executor 'local'
    tag "${d.key}"
    publishDir path: { "${d.out_dir}/08.knowledge" }, mode: 'copy'

    input:
    tuple val(d), path(profile), path(tables), path(curated), path(graph), path(summaries), path(summaries_json), path(changes), path(changes_json), path(seed_pubmed)
    val run_params
    path script
    path lib

    output:
    tuple val(d), path("${d.key}_knowledge_graph.md"), path("${d.key}_knowledge_graph.json"), emit: out

    script:
    def rp = toJson(run_params + [FetchDate: d.fetch_date, Rebuild: d.rebuild])
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    cat > run_params.json <<'JSON'
${rp}
JSON
    ${params.python} ${script} --profile ${profile} --tables ${tables} --curated ${curated} \\
        --graph-schema ${graph}/schema.json --summaries ${summaries} --changes ${changes} \\
        --changes-json ${changes_json} --seed-pubmed ${seed_pubmed} --run-params run_params.json \\
        --l2g-min ${params.L2gMin} --out-md ${d.key}_knowledge_graph.md --out-json ${d.key}_knowledge_graph.json
    """
}

/* WRITE_TABLE — the workbook, re-opened and checked before it is published. */
process WRITE_TABLE {
    executor 'local'
    tag "${d.key}"
    publishDir path: { "${d.out_dir}" }, mode: 'copy'

    input:
    tuple val(d), path(profile), path(tables), path(curated), path(graph), path(summaries), path(summaries_json), path(changes), path(changes_json), path(knowledge_md), path(knowledge_json)
    val run_params
    path script
    path lib

    output:
    tuple val(d), path("${d.workbook}"), path('workbook.json'), emit: out
    path "${d.workbook.replace('.xlsx', '_summary.xlsx')}"

    script:
    def rp = toJson(run_params + [FetchDate: d.fetch_date, Rebuild: d.rebuild])
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    cat > run_params.json <<'JSON'
${rp}
JSON
    ${params.python} ${script} --profile ${profile} --tables ${tables} --curated ${curated} \\
        --graph-schema ${graph}/schema.json --summaries ${summaries} --changes ${changes} \\
        --changes-json ${changes_json} --knowledge-json ${knowledge_json} --run-params run_params.json \\
        --l2g-min ${params.L2gMin} --eas-min ${params.EasMin} --out ${d.workbook} \\
        --out-summary ${d.workbook.replace('.xlsx', '_summary.xlsx')} --out-stats workbook.json
    """
}

/* VIEWER — the snapshot in one HTML file (offline; every edge with its source link). */
process VIEWER {
    executor 'local'
    tag "${d.key}"
    publishDir path: { "${d.out_dir}" }, mode: 'copy'

    input:
    tuple val(d), path(profile), path(graph), path(curated), path(tables), path(summaries), path(summaries_json), path(figures)
    val run_params
    path script
    path lib

    output:
    tuple val(d), path("${d.key}_knowledge_graph.html"), emit: html

    script:
    def rp = toJson(run_params + [FetchDate: d.fetch_date, Rebuild: d.rebuild])
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    cat > run_params.json <<'JSON'
${rp}
JSON
    ${params.python} ${script} --profile ${profile} --graph ${graph} --curated ${curated} --tables ${tables} \\
        --summaries ${summaries} --figures ${figures} --run-params run_params.json --out ${d.key}_knowledge_graph.html
    """
}

/* WRITE_RUN_MANIFEST — provenance, the index of the snapshot, and the pointer to the disease's newest snapshot. */
process WRITE_RUN_MANIFEST {
    executor 'local'
    tag "${d.key}"
    publishDir path: { "${d.out_dir}/_run_info" }, mode: 'copy', pattern: 'run_manifest.json'
    publishDir path: { "${params.results_root}/${d.key}" }, mode: 'copy', pattern: 'latest.json'
    publishDir path: { "${d.out_dir}" }, mode: 'copy', pattern: 'README.md'

    input:
    tuple val(d), path(profile), path(tables), path(curated), path(graph), path(viewer), path(workbook), path(workbook_stats), path(changes), path(changes_json)
    val run_json
    path script
    path lib

    output:
    tuple val(d), path('run_manifest.json'), path('latest.json'), emit: out
    path 'README.md'

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    cat > run.json <<'JSON'
${run_json}
JSON
    ${params.python} ${script} --profile ${profile} --run-json run.json --tidy ${tables}/tidy.json \\
        --workbook-stats ${workbook_stats} --changes-json ${changes_json} --snapshot ${d.fetch_date} \\
        --earlier '${d.earlier}' --snapshots-root '${d.snapshots_root}' --workbook ${d.workbook} \\
        --curation ${curated}/curation.json --graph-schema ${graph}/schema.json \\
        --out run_manifest.json --out-latest latest.json --out-readme README.md
    """
}

/* PROJECT_INDEX — results/index.json and results/README.md: every disease and its newest snapshot. */
process PROJECT_INDEX {
    executor 'local'
    publishDir "${params.results_root}", mode: 'copy'

    input:
    val keys
    path latest, stageAs: 'latest/*/latest.json'
    path script
    path lib

    output:
    path 'index.json'
    path 'index.html'
    path 'README.md'

    script:
    def given = keys.withIndex().collect { k, i -> "${k}=latest/${i + 1}/latest.json" }.join(' ')
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --project '${params.results_root}' --latest ${given} --out-index index.json \\
        --out-html index.html --out-readme README.md
    """
}


// ═════════════════════════════════════════════════════════════════════════════
// WORKFLOW
// ═════════════════════════════════════════════════════════════════════════════
workflow {
    def specs = diseaseSpecs()
    def spec = specs.collectEntries { s -> [(s.key): s] }
    def offline = flag(params.LlmOffline)
    def llm = llmJson()
    def smoke = (params.LlmLimit as int) > 0 || (params.MaxThemes as int) > 0
    specs.each { s ->
        log.info "DiseaseLens | ${s.key} | snapshot ${s.fetch_date}${s.rebuild ? ' (rebuild)' : ''} | relations ${s.curated ? 'curated' : 'read from PubMed'} | earlier: ${s.earlier ?: 'none'} | out ${s.out_dir}"
    }
    log.info "DiseaseLens | readers ${readersText()}${offline ? ' (cache only)' : ''}"
    if (smoke) {
        log.warn "DiseaseLens | smoke settings (LlmLimit ${params.LlmLimit}, MaxThemes ${params.MaxThemes}): these snapshots are partial"
    }
    def lib = script_file('ot_common.py')
    def llmlib = script_file('llm.py')
    def pubmed = script_file('pubmed.py')
    def fetch = script_file('fetch.py')

    PROFILE(channel.fromList(specs).map { s -> tuple(s, file(s.disease_dir), file(params.config_dir)) },
            script_file('profile.py'), lib)
    def prof = PROFILE.out.json
    def fetching = prof.filter { s, _p -> !s.rebuild }
    def saved = prof.filter { s, _p -> s.rebuild }

    // ---- sources and evidence: fetched ----
    FETCH_CORE(fetching, fetch, lib)
    SEED_PUBMED(fetching.filter { s, _p -> s.curated }, pubmed, lib)
    SEARCH_PUBMED(fetching.map { s, p -> tuple(s, p, 'genetics') }
                      .mix(fetching.filter { s, _p -> !s.curated }.map { s, p -> tuple(s, p, 'relations') }),
                  pubmed, lib)
    def searches = SEARCH_PUBMED.out.records.branch { _s, f ->
        relations: f.name.startsWith('relations')
        genetics: true
    }
    RELATIONS(keyed(prof).join(keyed(searches.relations)).map { k, p, r -> tuple(spec[k], p, r) }, llm, offline,
              script_file('relations.py'), lib, llmlib)
    // the seed of each disease: its curated list (in the profile) or the relations read from PubMed
    def seedFetched = keyed(SEED_PUBMED.out.records).map { k, r -> [k, file("${spec[k].disease_dir}/related_conditions.yaml"), r] }
        .mix(keyed(RELATIONS.out.seed))
    RESOLVE(keyed(prof).join(seedFetched).join(keyed(FETCH_CORE.out.raw))
                .map { k, p, sd, sp, raw -> tuple(spec[k], p, sd, sp, raw) },
            script_file('resolve.py'), lib, llmlib)
    FETCH_RELATED(RESOLVE.out.fetch_list, fetch, lib)
    SCREEN(keyed(prof).join(keyed(searches.genetics)).map { k, p, r -> tuple(spec[k], p, r) }, llm, offline,
           script_file('extract_genetics.py'), lib, llmlib)
    FETCH_FULLTEXT(SCREEN.out.ids, pubmed, lib)

    // ---- sources and evidence: saved (rebuild) ----
    def savedKeyed = keyed(saved)
    def raw = keyed(FETCH_CORE.out.raw).mix(savedKeyed.map { k, _p -> [k, file("${spec[k].out_dir}/00.raw", checkIfExists: true)] })
    def search = keyed(searches.genetics).mix(savedKeyed.map { k, _p -> [k, file("${spec[k].out_dir}/00.pubmed/pubmed_search.json.gz", checkIfExists: true)] })
    def fulltext = keyed(FETCH_FULLTEXT.out.json).mix(savedKeyed.map { k, _p -> [k, file("${spec[k].out_dir}/00.pubmed/fulltext.json.gz", checkIfExists: true)] })
    def seedPubmed = seedFetched.map { k, _sd, sp -> [k, sp] }
        .mix(savedKeyed.map { k, _p -> [k, file("${spec[k].out_dir}/00.pubmed/seed_pubmed.json.gz", checkIfExists: true)] })
    def conditions = keyed(RESOLVE.out.dir).mix(savedKeyed.map { k, _p -> [k, file("${spec[k].out_dir}/01.conditions", checkIfExists: true)] })
    def related = keyed(FETCH_RELATED.out.dir).mix(savedKeyed.map { k, _p -> [k, file("${spec[k].out_dir}/00.raw/diseases", checkIfExists: true)] })
    // the release of a saved snapshot: the gene symbols of a rebuild are looked up in it
    def release = raw.map { k, r -> [k, new groovy.json.JsonSlurper().parse(file("${r}/meta.json")).data_version ?: ''] }

    EXTRACT(keyed(prof).join(search).join(fulltext).join(release).join(raw)
                .map { k, p, s, f, rel, r -> tuple(spec[k], p, s, f, rel, r) },
            llm, offline, script_file('extract_genetics.py'), lib, llmlib)
    def literature = keyed(EXTRACT.out.dir)

    // the publication records: saved by a rebuild when they exist, else read from PubMed
    def savedPubs = savedKeyed.filter { k, _p -> file("${spec[k].out_dir}/00.pubmed/publications.json.gz").exists() }
        .map { k, _p -> [k, file("${spec[k].out_dir}/00.pubmed/publications.json.gz")] }
    def needPubs = conditions.filter { k, _c -> !(spec[k].rebuild && file("${spec[k].out_dir}/00.pubmed/publications.json.gz").exists()) }
    PUBLICATIONS(needPubs.join(search).join(seedPubmed).join(raw).join(related)
                     .map { k, c, s, sp, r, dz -> tuple(spec[k], s, sp, r, dz, c) },
                 script_file('publications.py'), lib)
    def publications = keyed(PUBLICATIONS.out.json).mix(savedPubs)

    // ---- curation, tables, graph, products ----
    CURATE(keyed(prof).join(publications).join(literature).join(raw).join(conditions)
               .map { k, p, pu, li, r, c -> tuple(spec[k], p, pu, li, r, c) },
           script_file('curate.py'), lib)
    def curated = keyed(CURATE.out.dir)
    TIDY_GRAPH(keyed(prof).join(raw).join(conditions).join(related).join(literature).join(curated)
                   .map { k, p, r, c, dz, li, cu -> tuple(spec[k], p, r, c, dz, li, cu) },
               script_file('tidy_graph.py'), lib)
    def tables = keyed(TIDY_GRAPH.out.tables)
    BUILD_GRAPH(keyed(prof).join(curated).join(tables).map { k, p, cu, t -> tuple(spec[k], p, cu, t) },
                script_file('build_graph.py'), lib)
    def graph = keyed(BUILD_GRAPH.out.dir)
    FIGURES(keyed(prof).join(tables).map { k, p, t -> tuple(spec[k], p, t) }, file(params.plot_style, checkIfExists: true),
            script_file('draw_figures.py'), lib)
    SUMMARISE(keyed(prof).join(tables).map { k, p, t -> tuple(spec[k], p, t) }, llm, offline,
              script_file('summarise.py'), lib, llmlib)
    def summaries = keyed(SUMMARISE.out.out)
    DIFF(TIDY_GRAPH.out.tables, script_file('diff_snapshots.py'), lib)
    def changes = keyed(DIFF.out.out)

    // the parameters that shape the results: recorded in the run manifest, and shown with their basis (the choices
    // registry) in the knowledge base and the workbook; the disease's own values (its ID, its queries) are added from
    // its profile
    def runParams = [ClppMin: params.ClppMin, H4Min: params.H4Min, MendelianMin: params.MendelianMin,
                     L2gMin: params.L2gMin, EasMin: params.EasMin, LlmReaders: readersText(), LlmBatch: params.LlmBatch,
                     LlmLimit: params.LlmLimit, MaxThemes: params.MaxThemes, Prescreen: params.Prescreen,
                     SentencesPerPaper: params.SentencesPerPaper, ThemeTargets: params.ThemeTargets,
                     MinThemeTargets: params.MinThemeTargets, MaxGwasStudies: params.MaxGwasStudies,
                     MaxOntologyChildren: params.MaxOntologyChildren, RegionWindowKb: params.RegionWindowKb,
                     GenomeWideP: params.GenomeWideP, RegionMinMaf: params.RegionMinMaf,
                     RegionMinMac: params.RegionMinMac, RareVariantP: params.RareVariantP, LlmOffline: offline]
    KNOWLEDGE(keyed(prof).join(tables).join(curated).join(graph).join(summaries).join(changes).join(seedPubmed)
                  .map { k, p, t, cu, g, su, sj, ch, cj, sp -> tuple(spec[k], p, t, cu, g, su, sj, ch, cj, sp) },
              runParams, script_file('build_knowledge.py'), lib)
    def knowledge = keyed(KNOWLEDGE.out.out)
    WRITE_TABLE(keyed(prof).join(tables).join(curated).join(graph).join(summaries).join(changes).join(knowledge)
                    .map { k, p, t, cu, g, su, sj, ch, cj, km, kj -> tuple(spec[k], p, t, cu, g, su, sj, ch, cj, km, kj) },
                runParams, script_file('write_table.py'), lib)
    VIEWER(keyed(prof).join(graph).join(curated).join(tables).join(summaries).join(keyed(FIGURES.out.dir))
               .map { k, p, g, cu, t, su, sj, fi -> tuple(spec[k], p, g, cu, t, su, sj, fi) },
           runParams, script_file('build_viewer.py'), lib)

    def scripts = ['ot_common.py', 'llm.py', 'profile.py', 'fetch.py', 'pubmed.py', 'relations.py', 'resolve.py',
                   'extract_genetics.py', 'publications.py', 'curate.py', 'tidy_graph.py', 'build_graph.py',
                   'draw_figures.py', 'summarise.py', 'diff_snapshots.py', 'build_knowledge.py', 'write_table.py',
                   'build_viewer.py', 'manifest.py']
    def run = [
        launched  : workflow.start.toString(),
        session   : workflow.sessionId.toString(),
        command   : workflow.commandLine.replace(projectDir.toString(), '<repository>').replace(System.getProperty('user.home'), '~'),
        nextflow  : workflow.nextflow.version.toString(),
        params    : runParams,
        readers   : readersText(),
        smoke     : smoke,
        script_md5: scripts.collectEntries { s -> [(s): script_file(s).text.md5()] },
    ]
    WRITE_RUN_MANIFEST(keyed(prof).join(tables).join(curated).join(graph).join(keyed(VIEWER.out.html))
                           .join(keyed(WRITE_TABLE.out.out)).join(changes)
                           .map { k, p, t, cu, g, v, wb, ws, ch, cj -> tuple(spec[k], p, t, cu, g, v, wb, ws, ch, cj) },
                       toJson(run), script_file('manifest.py'), lib)
    // the newest pointer of each disease of this run, passed in (its published copy may not be there yet)
    def pointers = WRITE_RUN_MANIFEST.out.out.map { s, _m, l -> [s.key, l] }.toSortedList { x, y -> x[0] <=> y[0] }
    PROJECT_INDEX(pointers.map { xs -> xs.collect { x -> x[0] } }, pointers.map { xs -> xs.collect { x -> x[1] } },
                  script_file('manifest.py'), lib)
}
