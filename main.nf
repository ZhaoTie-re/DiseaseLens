nextflow.enable.dsl = 2

// ═════════════════════════════════════════════════════════════════════════════
// DiseaseLens — a disease, the conditions related to it, and the genetics of
// both, as a dated, traceable knowledge graph that can be repeated (one run per
// --FetchDate, each compared with the snapshot before it). Configured for its
// first application, CTEPH (config/, --DiseaseCore, --PubmedQuery).
//
// External knowledge only: Open Targets (GraphQL), PubMed (E-utilities) and
// Europe PMC; nothing of our own cohort's results. One line of stages:
//   sources    00  Open Targets (CTEPH; each related condition), PubMed (the
//                  genetic-study search, the seed), Europe PMC Methods texts, and
//                  the PubMed record of every publication the graph uses
//   evidence   01  related conditions and the evidence for each relation
//              02  the genetic studies of CTEPH, extracted by an isolated language
//                  model (every item checked against its source text)
//   curation   03  publications left out (retracted, conference material) or
//                  flagged; variants named; findings classified by what they are
//                  about; text-mined links mapped to the wrong gene removed; the
//                  place and design class of each study (deterministic rules,
//                  config/curation.yaml, every rule applied logged)
//   tables     04  the long tables and the evidence ladder of the CTEPH genes
//   graph      05  typed nodes and edges, graph.json, a SQLite database, the schema
//   products   06  grounded summaries · 07 changes · 08 knowledge base ·
//              09 figures · the workbook · the viewer (one HTML file) · the index
//
//   FETCH_CTEPH ─┬──► RESOLVE ─► FETCH_DISEASE × condition ──────────────┐
//   SEED_PUBMED ─┘                                                        │
//   SEARCH_PUBMED ─► SCREEN ─► FETCH_FULLTEXT ─► EXTRACT ─────────────────┤
//                                     PUBLICATIONS ◄──────────────────────┤
//                                          └──► CURATE ─► TIDY_GRAPH ─► BUILD_GRAPH
//        BUILD_GRAPH ─┬─► FIGURES, SUMMARISE, DIFF ─► KNOWLEDGE ─► WRITE_TABLE ─┐
//                     └─► VIEWER ◄────────────────────────────────────────────┤
//                                                       WRITE_RUN_MANIFEST ◄──┘
//
// --Rebuild true re-derives the snapshot of --FetchDate from its saved sources
// (00.raw, 00.pubmed, 01.conditions, 02.literature): nothing is queried again
// except PubMed for the publication records when they are not yet saved, and no
// language model is called (--LlmOffline is implied: summaries come from the cache
// or are reported as not generated). Every process is local: compute nodes have
// no internet. Methodology: docs/METHODS.md. Outputs: docs/OUTPUTS.md.
// ═════════════════════════════════════════════════════════════════════════════

// =============================================================================
// SITE CONFIGURATION — override on the command line (`--key value`).
// FetchDate, out_dir and run_info_dir live in nextflow.config.
// =============================================================================
// ---- Paths
params.base_dir      = "${projectDir}"                                    // this repository
params.script_dir    = "${params.base_dir}/scripts"
params.seed          = "${params.base_dir}/config/related_conditions.yaml"  // curated, PubMed-verified relations
params.gene_themes   = "${params.base_dir}/config/gene_themes.yaml"        // gene themes of the CTEPH literature
params.plot_style    = "${params.base_dir}/scripts/plot_style.py"         // the figure style (vendored)
params.choices       = "${params.base_dir}/config/choices.yaml"            // every choice, its value source and basis
params.cteph_gwas    = "${params.base_dir}/config/cteph_gwas.yaml"         // the published CTEPH GWAS regions
params.curation      = "${params.base_dir}/config/curation.yaml"           // the curation rules
params.llm_cache_dir = "${params.base_dir}/cache/llm"                       // model answers, kept across snapshots
params.llm_workdir   = '/tmp/cteph_kg_llm'                                  // the model runs here, away from the project

// ---- Open Targets
params.ot_api        = 'https://api.platform.opentargets.org/api/v4/graphql'
params.DiseaseCore   = 'EFO_0022954'   // chronic thromboembolic pulmonary hypertension
params.ClppMin       = 0.01            // colocalization with a CTEPH locus: eCAVIAR CLPP ...
params.H4Min         = 0.8             // ... or coloc H4
params.MendelianMin  = 0.5             // smallest curated / Mendelian source score kept
params.L2gMin        = 0.5             // top locus-to-gene score for a gene-condition link
params.EasMin        = 0.5             // East Asian share of a GWAS sample to flag it
params.PageSize      = 100
params.EvidencePage  = 500

// ---- PubMed: the genetic studies of CTEPH
params.PubmedQuery   = '("chronic thromboembolic pulmonary hypertension"[tiab] OR CTEPH[tiab]) AND (gene*[tiab] OR genetic*[tiab] OR polymorphism*[tiab] OR variant*[tiab] OR mutation*[tiab] OR genome-wide[tiab] OR GWAS[tiab] OR HLA[tiab] OR sequencing[tiab] OR genotyp*[tiab] OR allele*[tiab] OR exome[tiab])'

// ---- Language model (isolated `claude -p`)
params.claude        = 'claude'
params.LlmModel      = ''      // '' = the CLI's default model
params.LlmEffort     = ''
params.LlmParallel   = 3
params.LlmBatch      = 5       // abstracts per screening call
params.Prescreen     = true    // drop reviews / records without a genetic-study term before the model reads them
params.LlmLimit      = 0       // smoke run: the first N PubMed records only (0 = all)
params.MaxThemes     = 0       // smoke run: the N largest mechanism themes only, no methods / findings summaries
params.SentencesPerPaper = 3
params.ThemeTargets  = 12      // targets whose sentences feed a theme summary
params.MinThemeTargets = 3
params.LlmOffline    = false   // true: cached model answers only, never call the model

// ---- Rebuild: re-derive the snapshot of FetchDate from its saved sources (00-02), no model call
params.Rebuild       = false

// ---- Tools (absolute: an activated env is not first on PATH at this site)
params.python = 'python3'  // pandas, xlsxwriter, openpyxl, yaml

// ---- Helpers
def script_file(String name) {
    return file("${params.script_dir}/${name}", checkIfExists: true)
}

def toJson(x) {
    return groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(x))
}

def llmArgs() {
    return "--cache-dir '${params.llm_cache_dir}' --workdir '${params.llm_workdir}' --model '${params.LlmModel}' " +
        "--effort '${params.LlmEffort}' --parallel ${params.LlmParallel} --claude '${params.claude}'"
}


/* FETCH_CTEPH — the core disease: targets, text-mined and GWAS evidence, target pathways,
 * its GWAS credible sets with their colocalizations, and the published GWAS regions
 * (config/cteph_gwas.yaml) looked up in Open Targets. The fetch date is in the task hash. */
process FETCH_CTEPH {
    executor 'local'
    tag "${params.DiseaseCore}"
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    val fetch_date
    path published
    path script
    path lib

    output:
    path '00.raw', emit: raw

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --mode cteph --api '${params.ot_api}' --disease ${params.DiseaseCore} \\
        --fetch-date ${fetch_date} --published ${published} --page-size ${params.PageSize} \\
        --evidence-page ${params.EvidencePage} --out-dir 00.raw
    """
}

/* SEED_PUBMED — the PubMed records of every PMID the curated seed cites. */
process SEED_PUBMED {
    executor 'local'
    tag 'seed'
    publishDir "${params.out_dir}/00.pubmed", mode: 'copy'

    input:
    path seed
    val fetch_date
    path script
    path lib

    output:
    path 'seed_pubmed.json.gz', emit: records

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    # fetch date: ${fetch_date}
    ${params.python} ${script} --mode seed --seed ${seed} --out seed_pubmed.json.gz
    """
}

/* RESOLVE — the related conditions and their evidence (clinical quotes checked, ontology, colocalization). */
process RESOLVE {
    executor 'local'
    tag 'conditions'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path raw
    path seed
    path seed_pubmed
    path script
    path lib
    path llm

    output:
    path '01.conditions', emit: dir
    path '01.conditions/fetch_list.tsv', emit: fetch_list

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --seed ${seed} --seed-pubmed ${seed_pubmed} --cteph-gwas ${raw}/gwas.json.gz \\
        --published ${raw}/published_loci.json.gz --api '${params.ot_api}' --clpp-min ${params.ClppMin} \\
        --h4-min ${params.H4Min} --out-dir 01.conditions
    """
}

/* FETCH_DISEASE — one related condition's GWAS and curated / Mendelian gene evidence. */
process FETCH_DISEASE {
    executor 'local'
    maxForks 3
    tag "${disease}"
    publishDir "${params.out_dir}/00.raw/diseases", mode: 'copy'

    input:
    tuple val(disease), val(mendelian)
    val fetch_date
    path script
    path lib

    output:
    path "disease.${disease}.json.gz", emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    # fetch date: ${fetch_date}
    ${params.python} ${script} --mode disease --api '${params.ot_api}' --disease ${disease} --mendelian ${mendelian} \\
        --threshold ${params.MendelianMin} --page-size ${params.PageSize} --out-dir .
    """
}

/* SEARCH_PUBMED — the PubMed records the CTEPH-genetics query finds, up to the fetch date. */
process SEARCH_PUBMED {
    executor 'local'
    tag 'search'
    publishDir "${params.out_dir}/00.pubmed", mode: 'copy'

    input:
    val fetch_date
    path script
    path lib

    output:
    path 'pubmed_search.json.gz', emit: records

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    cat > query.txt <<'Q'
${params.PubmedQuery}
Q
    ${params.python} ${script} --mode search --query "\$(cat query.txt)" --max-date ${fetch_date} --out pubmed_search.json.gz
    """
}

/* SCREEN — pass A of the extraction (abstracts); lists the included records with a PMCID. */
process SCREEN {
    executor 'local'
    tag 'screen'

    input:
    path records
    path script
    path lib
    path llm

    output:
    path 'fulltext_ids.txt', emit: ids

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --records ${records} ${llmArgs()} --api '${params.ot_api}' --batch ${params.LlmBatch} \\
        --limit ${params.LlmLimit} --prescreen ${params.Prescreen} --fulltext-list fulltext_ids.txt --out-dir screen
    """
}

/* FETCH_FULLTEXT — the Methods sections of the included open-access full texts (Europe PMC). */
process FETCH_FULLTEXT {
    executor 'local'
    tag 'fulltext'
    publishDir "${params.out_dir}/00.pubmed", mode: 'copy'

    input:
    path ids
    path script
    path lib

    output:
    path 'fulltext.json.gz', emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --mode fulltext --ids ${ids} --out fulltext.json.gz
    """
}

/* EXTRACT — passes A (cached) and B: the genetic studies of CTEPH, their findings and methods. */
process EXTRACT {
    executor 'local'
    tag 'extract'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path records
    path fulltext
    path script
    path lib
    path llm

    output:
    path '02.literature', emit: dir

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --records ${records} --fulltext ${fulltext} ${llmArgs()} --api '${params.ot_api}' \\
        --batch ${params.LlmBatch} --limit ${params.LlmLimit} --prescreen ${params.Prescreen} --out-dir 02.literature
    """
}

/* PUBLICATIONS — the PubMed record of every publication the graph uses (status, issue, affiliation). */
process PUBLICATIONS {
    executor 'local'
    tag 'publications'
    publishDir "${params.out_dir}/00.pubmed", mode: 'copy'

    input:
    path search
    path seed_records
    path raw
    path diseases
    path conditions
    path script
    path lib

    output:
    path 'publications.json.gz', emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --pubmed ${search} ${seed_records} --raw ${raw} --diseases ${diseases} \\
        --relations ${conditions}/relations.tsv --out publications.json.gz
    """
}

/* CURATE — what the graph holds: exclusions, variant names, what each finding is about, text-mining review. */
process CURATE {
    executor 'local'
    tag 'curation'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path publications
    path literature
    path raw
    path conditions
    path curation
    path script
    path lib

    output:
    path '03.curated', emit: dir

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --publications ${publications} --literature ${literature} --raw ${raw} \\
        --relations ${conditions}/relations.tsv --curation ${curation} --sentences ${params.SentencesPerPaper} \\
        --out-dir 03.curated
    """
}

/* TIDY_GRAPH — the long tables and the evidence ladder, from the curated literature. */
process TIDY_GRAPH {
    executor 'local'
    tag 'tables'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path raw
    path conditions
    path diseases
    path literature
    path curated
    path themes
    path script
    path lib

    output:
    path '04.tables', emit: tables

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --cteph-raw ${raw} --conditions ${conditions} --diseases ${diseases} \\
        --extract ${literature} --curated ${curated} --eas-min ${params.EasMin} \\
        --themes ${themes} --l2g-min ${params.L2gMin} --out-dir 04.tables
    """
}

/* BUILD_GRAPH — the typed graph: node and edge tables, graph.json, the SQLite database, the schema. */
process BUILD_GRAPH {
    executor 'local'
    tag 'graph'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path curated
    path tables
    path themes
    path script
    path lib

    output:
    path '05.graph', emit: dir
    path '05.graph/schema.json', emit: schema

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --curated ${curated} --tables ${tables} --themes ${themes} --out-dir 05.graph
    """
}

/* VIEWER — the snapshot in one HTML file (offline; every edge with its source link). */
process VIEWER {
    executor 'local'
    tag 'viewer'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path graph
    path curated
    path tables
    path summaries
    path script
    path lib

    output:
    path 'cteph_knowledge_graph.html', emit: html

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --graph ${graph} --curated ${curated} --tables ${tables} --summaries ${summaries} \\
        --out cteph_knowledge_graph.html
    """
}

/* FIGURES — the snapshot's figures (evidence grid of the CTEPH genes, condition map, shared CTEPH GWAS regions,
 * gene x condition map, literature by gene) in the shared figure style (scripts/plot_style.py), with one catalog
 * README. */
process FIGURES {
    executor 'local'
    tag 'figures'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path tables
    path themes
    path style
    path script
    path lib

    output:
    path '09.figures', emit: dir

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=\$PWD/.mpl
    ${params.python} ${script} --tables ${tables} --themes ${themes} --plot-style ${style} \\
        --l2g-min ${params.L2gMin} --curated-min ${params.MendelianMin} --clpp-min ${params.ClppMin} \\
        --out-dir 09.figures
    """
}

/* SUMMARISE — grounded summaries: mechanism themes, methods and findings of the genetic studies (offline: cache only). */
process SUMMARISE {
    executor 'local'
    tag 'summaries'
    publishDir "${params.out_dir}/06.summaries", mode: 'copy'

    input:
    path tables
    val offline
    path script
    path lib
    path llm

    output:
    path 'summaries.tsv', emit: tsv
    path 'summaries.json', emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --tables ${tables} ${llmArgs()} --min-targets ${params.MinThemeTargets} \\
        --top-targets ${params.ThemeTargets} --max-themes ${params.MaxThemes} --offline ${offline} \\
        --out summaries.tsv --out-json summaries.json
    """
}

/* DIFF — the changes since the newest earlier snapshot. */
process DIFF {
    executor 'local'
    tag 'changes'
    publishDir "${params.out_dir}/07.changes", mode: 'copy'

    input:
    path tables
    val current
    val earlier
    path script
    path lib

    output:
    path 'changes.tsv', emit: tsv
    path 'changes.json', emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --tables ${tables} --current ${current} --snapshots-root '${params.snapshots_root}' \\
        --earlier '${earlier}' --out changes.tsv --out-json changes.json
    """
}

/* KNOWLEDGE — the knowledge base of the snapshot. */
process KNOWLEDGE {
    executor 'local'
    tag 'knowledge'
    publishDir "${params.out_dir}/08.knowledge", mode: 'copy'

    input:
    path tables
    path curated
    path schema
    path summaries
    path changes
    path changes_json
    path seed
    path seed_pubmed
    path themes
    path choices
    val run_params
    path script
    path lib

    output:
    path 'cteph_knowledge_graph.md', emit: md
    path 'cteph_knowledge_graph.json', emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    cat > run_params.json <<'JSON'
${run_params}
JSON
    ${params.python} ${script} --tables ${tables} --curated ${curated} --graph-schema ${schema} \\
        --summaries ${summaries} --changes ${changes} \\
        --changes-json ${changes_json} --seed ${seed} --seed-pubmed ${seed_pubmed} --themes ${themes} \\
        --choices ${choices} --run-params run_params.json --l2g-min ${params.L2gMin} \\
        --out-md cteph_knowledge_graph.md --out-json cteph_knowledge_graph.json
    """
}

/* WRITE_TABLE — the workbook, re-opened and checked before it is published. */
process WRITE_TABLE {
    executor 'local'
    tag 'workbook'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path tables
    path curated
    path schema
    path summaries
    path changes
    path changes_json
    path knowledge_json
    path themes
    path choices
    val run_params
    val workbook
    path script
    path lib

    output:
    path "${workbook}", emit: xlsx
    path 'workbook.json', emit: stats

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    cat > run_params.json <<'JSON'
${run_params}
JSON
    ${params.python} ${script} --tables ${tables} --curated ${curated} --graph-schema ${schema} \\
        --summaries ${summaries} --changes ${changes} \\
        --changes-json ${changes_json} --knowledge-json ${knowledge_json} --themes ${themes} \\
        --choices ${choices} --run-params run_params.json \\
        --l2g-min ${params.L2gMin} --eas-min ${params.EasMin} --out ${workbook} --out-stats workbook.json
    """
}

/* WRITE_RUN_MANIFEST — provenance, the index of the snapshot, and the pointer to the newest snapshot. */
process WRITE_RUN_MANIFEST {
    executor 'local'
    tag 'manifest'
    publishDir "${params.run_info_dir}", mode: 'copy', pattern: 'run_manifest.json'
    publishDir "${params.results_root}", mode: 'copy', pattern: 'latest.json'
    publishDir "${params.out_dir}", mode: 'copy', pattern: 'README.md'

    input:
    path tables
    path curated
    path schema
    path viewer
    path workbook_stats
    path changes_json
    val run_json
    val current
    val earlier
    val workbook
    path script
    path lib

    output:
    path 'run_manifest.json'
    path 'latest.json'
    path 'README.md'

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    cat > run.json <<'JSON'
${run_json}
JSON
    ${params.python} ${script} --run-json run.json --tidy ${tables}/tidy.json --workbook-stats ${workbook_stats} \\
        --changes-json ${changes_json} --snapshot ${current} --earlier '${earlier}' \\
        --snapshots-root '${params.snapshots_root}' --workbook ${workbook} --curation ${curated}/curation.json \\
        --graph-schema ${schema} --out run_manifest.json --out-latest latest.json --out-readme README.md
    """
}


// ═════════════════════════════════════════════════════════════════════════════
// WORKFLOW
// ═════════════════════════════════════════════════════════════════════════════
workflow {
    def fetchDate = params.FetchDate ? params.FetchDate.toString() : ''
    if (!(fetchDate ==~ /\d{4}-\d{2}-\d{2}/)) {
        error "DiseaseLens | --FetchDate YYYY-MM-DD is required (the snapshot name), got '${fetchDate}'"
    }
    if (!params.out_dir.toString().endsWith("/snapshots/${fetchDate}")) {
        error "DiseaseLens | out_dir ${params.out_dir} is not the snapshot of ${fetchDate}"
    }
    def smoke = (params.LlmLimit as int) > 0 || (params.MaxThemes as int) > 0
    if (smoke && !fetchDate.endsWith('-smoke') && !params.out_dir.toString().contains('smoke')) {
        log.warn "DiseaseLens | smoke settings (LlmLimit ${params.LlmLimit}, MaxThemes ${params.MaxThemes}): this snapshot is partial"
    }
    def root = file(params.snapshots_root)
    def earlier = root.exists()
        ? root.listFiles().findAll { d -> d.isDirectory() && d.name != fetchDate && file("${d}/_run_info/run_manifest.json").exists() }
            .collect { d -> d.name }.sort()
        : []
    def earlierCsv = earlier.join(',')
    def workbook = "cteph_knowledge_graph_${fetchDate}.xlsx"
    log.info "DiseaseLens | snapshot ${fetchDate} | earlier snapshots: ${earlier ? earlierCsv : 'none'} | out ${params.out_dir}"

    def lib = script_file('ot_common.py')
    def llm = script_file('llm.py')
    def seed = file(params.seed, checkIfExists: true)
    def themes = file(params.gene_themes, checkIfExists: true)
    def choices = file(params.choices, checkIfExists: true)
    def curation = file(params.curation, checkIfExists: true)
    def published = file(params.cteph_gwas, checkIfExists: true)
    def rebuild = params.Rebuild as boolean
    def offline = rebuild || (params.LlmOffline as boolean)
    def raw = null
    def diseases = null
    def searchRecs = null
    def seedRecs = null
    def conditions = null
    def literature = null
    def publications = null
    if (rebuild) {
        // the saved sources of the snapshot: nothing is fetched or extracted again
        def d = params.out_dir
        log.info "DiseaseLens | rebuild of ${fetchDate} from its saved sources; no model call"
        raw        = channel.value(file("${d}/00.raw", checkIfExists: true))
        diseases   = channel.fromPath("${d}/00.raw/diseases/disease.*.json.gz", checkIfExists: true).collect()
        searchRecs = channel.value(file("${d}/00.pubmed/pubmed_search.json.gz", checkIfExists: true))
        seedRecs   = channel.value(file("${d}/00.pubmed/seed_pubmed.json.gz", checkIfExists: true))
        conditions = channel.value(file("${d}/01.conditions", checkIfExists: true))
        literature = channel.value(file("${d}/02.literature", checkIfExists: true))
        def saved = file("${d}/00.pubmed/publications.json.gz")
        if (saved.exists()) {
            publications = channel.value(saved)
        } else {
            PUBLICATIONS(searchRecs, seedRecs, raw, diseases, conditions, script_file('publications.py'), lib)
            publications = PUBLICATIONS.out.json
        }
    } else {
        FETCH_CTEPH(fetchDate, published, script_file('fetch.py'), lib)
        SEED_PUBMED(seed, fetchDate, script_file('pubmed.py'), lib)
        RESOLVE(FETCH_CTEPH.out.raw, seed, SEED_PUBMED.out.records, script_file('resolve.py'), lib, llm)
        def fetchList = RESOLVE.out.fetch_list.splitCsv(sep: '\t').map { r -> [r[0], r[1]] }
        FETCH_DISEASE(fetchList, fetchDate, script_file('fetch.py'), lib)
        SEARCH_PUBMED(fetchDate, script_file('pubmed.py'), lib)
        SCREEN(SEARCH_PUBMED.out.records, script_file('extract_genetics.py'), lib, llm)
        FETCH_FULLTEXT(SCREEN.out.ids, script_file('pubmed.py'), lib)
        EXTRACT(SEARCH_PUBMED.out.records, FETCH_FULLTEXT.out.json, script_file('extract_genetics.py'), lib, llm)
        raw = FETCH_CTEPH.out.raw
        diseases = FETCH_DISEASE.out.json.collect()
        searchRecs = SEARCH_PUBMED.out.records
        seedRecs = SEED_PUBMED.out.records
        conditions = RESOLVE.out.dir
        literature = EXTRACT.out.dir
        PUBLICATIONS(searchRecs, seedRecs, raw, diseases, conditions, script_file('publications.py'), lib)
        publications = PUBLICATIONS.out.json
    }
    CURATE(publications, literature, raw, conditions, curation, script_file('curate.py'), lib)
    TIDY_GRAPH(raw, conditions, diseases, literature, CURATE.out.dir, themes, script_file('tidy_graph.py'), lib)
    BUILD_GRAPH(CURATE.out.dir, TIDY_GRAPH.out.tables, themes, script_file('build_graph.py'), lib)
    FIGURES(TIDY_GRAPH.out.tables, themes, file(params.plot_style, checkIfExists: true),
            script_file('draw_figures.py'), lib)
    SUMMARISE(TIDY_GRAPH.out.tables, offline, script_file('summarise.py'), lib, llm)
    DIFF(TIDY_GRAPH.out.tables, fetchDate, earlierCsv, script_file('diff_snapshots.py'), lib)
    // the parameters that shape the results: recorded in the run manifest, and shown with their basis
    // (config/choices.yaml) in the knowledge base and the workbook
    def runParams = [FetchDate: fetchDate, DiseaseCore: params.DiseaseCore, ClppMin: params.ClppMin,
                     H4Min: params.H4Min, MendelianMin: params.MendelianMin, L2gMin: params.L2gMin,
                     EasMin: params.EasMin, PubmedQuery: params.PubmedQuery, LlmModel: params.LlmModel,
                     LlmEffort: params.LlmEffort, LlmBatch: params.LlmBatch, LlmLimit: params.LlmLimit,
                     MaxThemes: params.MaxThemes, Prescreen: params.Prescreen,
                     SentencesPerPaper: params.SentencesPerPaper, ThemeTargets: params.ThemeTargets,
                     MinThemeTargets: params.MinThemeTargets, Rebuild: rebuild, LlmOffline: offline]
    KNOWLEDGE(TIDY_GRAPH.out.tables, CURATE.out.dir, BUILD_GRAPH.out.schema, SUMMARISE.out.tsv, DIFF.out.tsv,
              DIFF.out.json, seed, seedRecs, themes, choices, toJson(runParams), script_file('build_knowledge.py'),
              lib)
    WRITE_TABLE(TIDY_GRAPH.out.tables, CURATE.out.dir, BUILD_GRAPH.out.schema, SUMMARISE.out.tsv, DIFF.out.tsv,
                DIFF.out.json, KNOWLEDGE.out.json, themes, choices, toJson(runParams), workbook,
                script_file('write_table.py'), lib)
    VIEWER(BUILD_GRAPH.out.dir, CURATE.out.dir, TIDY_GRAPH.out.tables, SUMMARISE.out.tsv,
           script_file('build_viewer.py'), lib)

    def scripts = ['ot_common.py', 'llm.py', 'fetch.py', 'pubmed.py', 'resolve.py', 'extract_genetics.py',
                   'publications.py', 'curate.py', 'tidy_graph.py', 'build_graph.py', 'draw_figures.py',
                   'summarise.py', 'diff_snapshots.py', 'build_knowledge.py', 'write_table.py', 'build_viewer.py',
                   'manifest.py']
    def run = [
        launched  : workflow.start.toString(),
        session   : workflow.sessionId.toString(),
        command   : workflow.commandLine,
        nextflow  : workflow.nextflow.version.toString(),
        params    : runParams + [seed_md5: seed.text.md5(), gene_themes_md5: themes.text.md5(),
                                 choices_md5: choices.text.md5(), curation_md5: curation.text.md5()],
        smoke     : smoke,
        earlier   : earlier,
        script_md5: scripts.collectEntries { s -> [(s): script_file(s).text.md5()] },
    ]
    WRITE_RUN_MANIFEST(TIDY_GRAPH.out.tables, CURATE.out.dir, BUILD_GRAPH.out.schema, VIEWER.out.html,
                       WRITE_TABLE.out.stats, DIFF.out.json, toJson(run), fetchDate, earlierCsv, workbook,
                       script_file('manifest.py'), lib)
}
