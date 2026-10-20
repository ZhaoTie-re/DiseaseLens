# DiseaseLens

**Dated, traceable knowledge graphs of the genetics of a disease and of the conditions related to it, for any disease
in Open Targets.**

Give DiseaseLens the Open Targets ID of a disease. It collects what public sources say about the disease's genetics,
reads the genetic literature with language models, checks every item against its source text, and writes one dated
snapshot per disease: an offline viewer, a workbook organized by question, a knowledge base with figures, and the
graph as tables, JSON and SQLite. Several diseases can be built in one run.

## What the graph holds

- **Genes of the disease**, each on an evidence ladder (strongest level wins):

  | Level | Code | Rule |
  |---|---|---|
  | 1 | G / C | a gene of a genome-wide significant GWAS region of the disease (G), or curated gene–disease evidence for it (C) |
  | 2 | G\* | a gene of a region significant only in a secondary analysis of a published GWAS |
  | 3–7 | | from the genetic studies in PubMed: associated in two or more studies; in one; a variant identified without a test; tested, not associated; other result |

  Levels 3–6 count only findings about the disease itself; a finding about an outcome among patients is kept and
  marked, and does not raise a level.
- **Genetic studies of the disease**, with design, place, sample sizes, findings and methods, each finding with its
  source sentence.
- **Related conditions**, each with its evidence: a verbatim quote of a PubMed abstract (clinical), the disease
  ontology, or a GWAS signal shared with a region of the disease. A study that found no association is shown as such
  and never relates a condition.
- **Shared genetics**: the regions and genes the disease shares with each condition, and the condition's own GWAS
  and curated genes.
- **Context**: mechanism themes of the text-mined literature, and the evidence from East Asian populations.

The graph has 7 node types (condition, gene, region, publication, GWAS study, theme, pathway) and 12 edge types;
every edge keeps its evidence and its source (PMID or database record).

## Sources and method, in brief

| Source | What it gives |
|---|---|
| [Open Targets Platform](https://platform.opentargets.org) (GraphQL API) | the disease and its ontology; targets and text-mined literature; GWAS studies, credible sets, locus-to-gene scores and colocalizations; curated gene–disease evidence (ClinGen, Gene2Phenotype, PanelApp, Orphanet, ClinVar, UniProt, gene burden); the same for each related condition |
| PubMed (NCBI E-utilities) | the genetic studies of the disease; the clinical relations (guidelines, reviews, cohort and case-control studies); the record of every publication the graph uses |
| Europe PMC | the Methods sections of open-access full texts; affiliations PubMed lacks |
| Ensembl REST | the position of a published lead variant that Open Targets does not index |

1. **Profile.** The disease is its Open Targets ID; its names and PubMed terms are derived from the record, each with
   its rule.
2. **Sources.** Open Targets is read in its current release; PubMed is searched for the records it held on the fetch
   date. Everything fetched is saved with the snapshot.
3. **Related conditions.** Clinical relations are taken from a curated list when the disease has one, else read from
   PubMed by the model readers; the ontology neighbors and the colocalizations with the disease's GWAS regions are
   added by rule.
4. **Literature.** A rule pre-screen; then every record is read by **two language models that never see each other's
   answer**; a **third model adjudicates** where they differ; every gene, variant, number and quote must occur in the
   source text, or the item is dropped and logged.
5. **Curation** by rule: retracted articles and conference material are left out; each study gets its place and
   design class; each variant one name; each finding what it is about; text-mined links to the wrong gene are removed.
6. **Tables and ladder**, **graph**, then the **products**, with summaries written by one model, checked claim by claim
   by the other, and ruled on by the adjudicator.

Every threshold and list is a parameter or a configuration file, recorded in each snapshot with its value and its
basis (the choices registry). The full methods are in [docs/METHODS.md](docs/METHODS.md), every output file in
[docs/OUTPUTS.md](docs/OUTPUTS.md), and the design decisions with their history in
[docs/STUDY_NOTES.md](docs/STUDY_NOTES.md).

## Build the knowledge graph of your disease

### 1. Install

```bash
git clone https://github.com/ZhaoTie-re/DiseaseLens.git
cd DiseaseLens
conda env create -f environment.yml && conda activate diseaselens    # Python, packages, Java 17, Nextflow
```

Without conda: Python 3.10+ with `pip install -r requirements.txt`, Java 17 and [Nextflow](https://www.nextflow.io).
The machine that launches the pipeline needs internet access (Open Targets, NCBI, Europe PMC, Ensembl, the model
backend); every process runs on it.

### 2. Configure the language models

```bash
cp conf/local.config.example conf/local.config      # git ignores it
```

Choose one block in `conf/local.config`:

| You have | Block | Then |
|---|---|---|
| the Claude Code command-line tool, logged in | (1) `claude-cli` | nothing (`claude update` if it does not know a model) |
| an Anthropic API key | (2) `anthropic` | `pip install anthropic`; `export ANTHROPIC_API_KEY=...` |
| an OpenAI-compatible endpoint (OpenAI, DeepSeek, vLLM, Ollama, ...) | (3) `openai`, with `base_url` | `export OPENAI_API_KEY=...` (a local server may need none) |

- The two readers and the adjudicator are **three different models, each named**: a backend's default changes with
  its version, and the adjudicator must be neither reader (the pipeline refuses it otherwise).
- **Keys** are read from the environment only; the configuration names the variable, never the key.
- **Isolation.** The command-line tool runs with every customization off, in a scratch directory: no file, memory or
  instruction of yours reaches a model. Logged in to a claude.ai account, the tool adds the account's e-mail address
  to what the model sees; any e-mail address not in the source text is removed from every answer before it is
  cached. The API backends send the prompt alone.
- Set `python` in `conf/local.config` if the Python with the packages is not the first `python3` on your `PATH`.

Then check the setup:

```bash
NO_LLM=1 ./verify.sh      # sections 1-6 must pass; without NO_LLM it also tests the isolation with one small call
```

### 3. Write the profile of the disease

Find the disease's ID at [platform.opentargets.org](https://platform.opentargets.org) (EFO_, MONDO_, HP_,
Orphanet_ ...; never a name, which can name several diseases), then:

```bash
python3 scripts/new_disease.py MONDO_0017148              # writes diseases/hpah/profile.yaml
python3 scripts/new_disease.py <ID> --short ABC --key abc  # your own short name or folder name
```

The profile is written from the Open Targets record, every value with the rule that chose it:

| Value | Rule |
|---|---|
| `name` | the Open Targets name |
| `short` (and the folder name) | the shortest acronym among the exact synonyms that spells the name; else the name itself if it has at most 11 characters; else the initials, **made up and marked `# CHECK`**: a label only, never searched; replace it with the abbreviation the literature uses (or give `--short`) |
| `pubmed: terms` | the name (always), the exact synonyms, the short name and the names of `evidence_also` entries, in title or abstract; left out: a one-word term that is also a gene (PAH, phenylalanine hydroxylase), a term that names another disease, acronyms of two characters or fewer, database forms ("A, B" is also searched as "B A"); an acronym only beside a word of the name. Each term left out is listed with its reason |
| `evidence_also` | for a grouping name ("A and/or B"), the Open Targets entries of its parts whose evidence is fetched with the disease's own (PVOD/PCH: the literature of pulmonary veno-occlusive disease is held under its own entry) |

Read the profile before the first run. You can:
- add a former name of the disease to `terms` (HPAH: "familial primary pulmonary hypertension"), so that older
  studies are found;
- add a `scope`: one or two sentences on which studies count as studies of the disease. Useful for a subtype, whose
  patients are often studied together with others;
- keep everything else at its default, or add curated files to the disease folder, each replacing a default:

| File in `diseases/<key>/` | Replaces |
|---|---|
| `related_conditions.yaml` | the clinical relations read from PubMed: a curated list, each statement with a PMID and its verbatim quote |
| `gwas_regions.yaml` | the regions from Open Targets: the regions of a published GWAS table |
| `gene_themes.yaml` | Reactome top-level pathways as gene functions |
| `curation.yaml` | (adds) findings that are not about the disease; text-mining corrections |
| `choices.yaml` | (adds) the bases of your curated files in the choices registry |
| `text/<product>.yaml` | (adds) prose about the disease in the products |

`diseases/cteph/` is a fully curated example; `diseases/hpah/`, `ipah/` and `pvod_pch/` are profiles written from
the ID.

### 4. A small test run

```bash
nextflow run main.nf -c conf/local.config --diseases hpah --FetchDate $(date +%F) \
    --LlmLimit 10 --MaxThemes 2 --results_root results/_smoke
```

The first 10 PubMed records only, into `results/_smoke`: a few minutes and a few dozen model calls. A disease can be
named by its folder or by its ID (`--diseases MONDO_0017148`).

### 5. The full run

```bash
nextflow run main.nf -c conf/local.config --diseases hpah --FetchDate $(date +%F)
nextflow run main.nf -c conf/local.config --diseases hpah,ipah --FetchDate $(date +%F)     # several diseases
nextflow run main.nf -c conf/local.config -params-file projects/pulmonary_hypertension.yaml
```

| Option | Effect |
|---|---|
| `--diseases a,b` | the diseases, by folder or ID; in a params file a list, each entry `{key, fetch_date, rebuild}` if they differ |
| `--FetchDate YYYY-MM-DD` | the snapshot date: PubMed records that entered PubMed up to it; Open Targets as served on the day of the run |
| `-resume` | after a stop (an error, a network failure, a model account's usage limit), continue where it stopped |
| `--LlmOffline true` | cached model answers only; no model is called (a record without a cached answer is reported as not read) |
| `--Rebuild true` (or `--Rebuild key`) | re-derive the snapshot of that date from its saved sources (step 7) |
| `--results_root <dir>` | write elsewhere than `results/` |

**Time and cost** grow with the genetic literature of the disease. With three Claude models through the
command-line tool, first builds took: IPAH, 911 genetic-study records (336 read after the pre-screen) and 587
relation records, about US$60 and four hours; PVOD/PCH, about 190 records, about US$20. Every answer is cached as
soon as it arrives, so a stopped run resumes without paying twice.

### 6. Check and open the results

```bash
./verify.sh                    # every disease: configuration, privacy, and the checks of each snapshot
DISEASES=hpah ./verify.sh      # one disease
```

`results/index.html` lists every disease with its newest snapshot (`results/<key>/latest.json` points to it). In
`results/<key>/snapshots/<date>/`:

| To | Open |
|---|---|
| browse the graph: search, a page per gene, condition, region and study, every source linked; "How it is built" explains every rule as applied to this disease | `<key>_knowledge_graph.html`, in any browser, offline |
| read the evidence by question | `<key>_knowledge_graph_<date>.xlsx`; `..._summary.xlsx` without the appendix, to send |
| read the knowledge base, with figures | `08.knowledge/<key>_knowledge_graph.md`, `09.figures/` |
| see what the two models disagreed on, and how it was decided | `03.curated/model_review.tsv`; sheet 10 of the workbook |
| see every rule applied by the curation | `03.curated/curation_log.tsv` |
| query from the shell or Python | `05.graph/<key>_knowledge_graph.sqlite` with `scripts/kg.py` |

```bash
python3 scripts/kg.py --disease hpah gene BMPR2           # level, function, every finding with its sentence
python3 scripts/kg.py --disease cteph condition "pulmonary embolism"
python3 scripts/kg.py --disease cteph studies --place Japan
python3 scripts/kg.py --disease ipah sql "SELECT label, disease_level, disease_code FROM node_gene WHERE disease_level <= 2"
python3 scripts/kg.py --disease cteph --json gene FGG      # JSON instead of Markdown
```

In Python, with `scripts/` on the path: `from kg import KG; KG(disease='hpah').gene('BMPR2')`. The query tool needs
only the standard library.

### 7. Update, rebuild, change a rule

- **A later snapshot.** Run again with a new `--FetchDate`. Sources are read again; only new or changed records go to
  the models (answers are cached per record); the snapshot is compared with the previous one (`07.changes`, sheet 12).
- **A rule changed.** Every threshold is a parameter of [nextflow.config](nextflow.config) (set it on the command
  line or in a params file); the shared lists are in `config/` (`search.yaml`: the PubMed blocks and the pre-screen;
  `curation.yaml`: the curation rules; `choices.yaml`: the registry with the basis of each choice). A disease's own
  entries go into its folder. Every product shows the values of the run.
- **Rebuild.** `--Rebuild true --FetchDate <the snapshot's date>` re-derives a snapshot from its saved sources after a
  change of rule, without reading the sources again. It does not re-derive the related conditions (their resolution
  queries Open Targets); a change to the rules of the relations takes a new fetch.

### 8. Share and cite

Send the files: the viewer is one HTML file that works offline; the summary workbook opens in Excel and WPS and
prints on A4; zipping `results/` keeps `index.html` and every link working from disk. Every viewer page carries the
snapshot, its sources and how to cite it.

Cite a snapshot as: DiseaseLens, `<disease>` knowledge graph, snapshot `<date>`, Open Targets `<release>`,
https://github.com/ZhaoTie-re/DiseaseLens. Open Targets data are open (CC0 1.0); abstracts and full texts remain
under their publishers' copyright or their articles' licenses, and the graph keeps short source sentences and
quotes as evidence, each linked to its PubMed record.

### Troubleshooting

| Symptom | Cause, and what to do |
|---|---|
| the run stops at a model account's usage limit | wait for the limit to reset, then the same command with `-resume` |
| `no short name: give short` | the name is long and has no acronym among its synonyms: `new_disease.py <ID> --short ABC`, or set `short` |
| `WARNING: the query finds no record` | the terms find nothing in PubMed: add the names the literature uses to `pubmed: terms` |
| `PubMed reports an error for the query` | a malformed term in `pubmed: terms` or `query` |
| `curation items that match no finding (stale config)` | an entry of the disease's `curation.yaml` no longer matches the literature: read it against the new snapshot, then edit or remove it |
| `two readers need an adjudicator` | name a third model in `llm.adjudicator` |
| a model name the command-line tool does not know | `claude update` |
| `Open Targets serves <release>, the snapshot is <release>` | a rebuild needs the release of its snapshot; build a new snapshot instead |

## Repository layout

| Path | Content |
|---|---|
| `main.nf`, `nextflow.config` | the pipeline and its parameters |
| `diseases/<key>/` | one folder per disease: its profile and optional curated files |
| `config/` | the shared rules: PubMed blocks and pre-screen, curation rules, the choices registry |
| `conf/local.config.example` | your settings (Python, language models), to copy to `conf/local.config` |
| `projects/` | an example params file: several diseases in one run |
| `scripts/` | the step scripts, `new_disease.py`, the query tool `kg.py`, `plot_style.py` (the figure style) |
| `verify.sh` | pre-flight and post-run checks (read-only) |
| `docs/` | [METHODS.md](docs/METHODS.md), [OUTPUTS.md](docs/OUTPUTS.md), [STUDY_NOTES.md](docs/STUDY_NOTES.md) |

The repository holds the pipeline only. What a run makes stays on the machine that made it, and git ignores it:
`results/` (the snapshots), `cache/` (the model answers and a local usage log), `work/` and `.nextflow/` (the run
state), and `conf/local.config`.

## Requirements

- Nextflow (tested with 26.04) and Java 17.
- Python 3 (tested with 3.12) with [requirements.txt](requirements.txt); `anthropic` only for the Anthropic backend.
- Internet access from the launching machine.

## Author

Zhao Tie ([@ZhaoTie-re](https://github.com/ZhaoTie-re)).
