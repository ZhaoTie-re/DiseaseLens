# CTEPH knowledge graph

A dated, reproducible knowledge graph of chronic thromboembolic pulmonary hypertension (CTEPH, [`EFO_0022954`](https://platform.opentargets.org/disease/EFO_0022954)). It covers the genes with genetic evidence for CTEPH, the genetic studies behind them, the conditions related to CTEPH, and the genetics those conditions share with it.

The graph is built from external knowledge only: the Open Targets Platform, the published CTEPH GWAS (Liley et al., 2024), PubMed and Europe PMC. It was built as background for interpreting the signals of a CTEPH whole-genome sequencing study, and none of that study's results enter the graph.

Each run writes one snapshot, named by its fetch date. A snapshot holds:
- the sources as fetched;
- the curated evidence;
- the typed graph (tables, `graph.json`, SQLite);
- the products: an offline HTML viewer, an Excel workbook, a knowledge base with figures, and an index.

Every choice (threshold, list, grouping, rule) is recorded with its value and its basis, and every statement links to its publication or database record.

**Current snapshot: [2026-10-16](results/snapshots/2026-10-16/README.md)** (Open Targets 26.09; PubMed searched up to 2026-10-16; every publication checked in PubMed on 2026-10-17).

## At a glance (snapshot 2026-10-16)

**CTEPH genes.** 67 genes sit on a seven-level evidence ladder (levels 1–7: 5, 5, 3, 9, 31, 11, 3).
- The nine published GWAS regions:
  - genome-wide significant (G): *ABO*, *FGG*, *F11*, *MYH7B*, *HLA-DRA*;
  - significant only after leveraging pulmonary embolism (G\*): *F2*, *TSPAN15*, *SLC44A2*, *F5*/*NME7*.
- Associated in two or more genetic studies (level 3): *FGA*, *HLA-B*, *HLA-DPB1*.

**Genetic studies of CTEPH.** 49 studies with 138 findings, from 471 PubMed records; 15 were done in Japan and 11 in China.

**Curation.**
- 2,151 publications were checked in PubMed; retracted articles and conference material are not part of the graph.
- 15 findings concern something other than CTEPH itself, for example an outcome within CTEPH. They are kept and marked, and do not raise a gene's level.
- 101 text-mined links mapped to the wrong gene are dropped.

**Related conditions.** 136 conditions, entering through four kinds of evidence (a condition can have several):

| Evidence | Conditions |
|---|---:|
| Clinical | 27 |
| Ontology | 11 |
| Colocalization | 67 |
| Lead-variant overlap | 70 |

Venous thromboembolism, deep vein thrombosis, pulmonary embolism and VWF level each share six of the seven CTEPH regions that Open Targets can assess.

**Graph.** 3,077 nodes and 5,816 edges (7 node types, 12 edge types).

## What it answers

The workbook, the knowledge base and the viewer follow the same twelve questions, in four parts.

**A. CTEPH genetics**
1. Which genes have genetic evidence for CTEPH, how strong is it, and what do the genes do?
2. What does each study report for each gene, in its own words?
3. Which genetic studies of CTEPH exist (design, place, size), and how were they done?

**B. Related conditions**

4. Which conditions are related to CTEPH, and on what evidence?
5. Which CTEPH GWAS regions do they share?
6. Which genes do they share?
7. What are their own genetics?

**C. Context**

8. Which mechanisms does the CTEPH literature describe?
9. What evidence comes from East Asian and Japanese populations?

**D. Provenance**

10. What did the curation leave out or correct, and why?
11. On what basis was each choice made?
12. What changed since the previous snapshot?

**The evidence ladder.** Each gene takes the strongest level that applies:

| Level | Code | Meaning |
|---|---|---|
| 1 | G | genome-wide significant in the CTEPH GWAS |
| 2 | G\* | significant only after leveraging pulmonary embolism |
| 3 | | associated with CTEPH in two or more studies |
| 4 | | associated with CTEPH in one study |
| 5 | | variant identified, no association test |
| 6 | | tested, not associated |
| 7 | | other result only |

Levels 1–2 come from the nine regions in Table 1 of Liley et al. (2024). Levels 3–7 come from the genetic studies of CTEPH in PubMed. Levels 3–6 count only findings about CTEPH itself.

## Open the current snapshot

| To | Open |
|---|---|
| browse the graph: search, a page per gene, condition, region and study, every source linked | [`cteph_knowledge_graph.html`](results/snapshots/2026-10-16/cteph_knowledge_graph.html). Download it, or clone the repository, and open it in any browser; it works offline. GitHub shows only its source. |
| read the evidence by question | [`cteph_knowledge_graph_2026-10-16.xlsx`](results/snapshots/2026-10-16/cteph_knowledge_graph_2026-10-16.xlsx) |
| read the knowledge base, with figures | [`08.knowledge/cteph_knowledge_graph.md`](results/snapshots/2026-10-16/08.knowledge/cteph_knowledge_graph.md) |
| see the figures | [`09.figures/`](results/snapshots/2026-10-16/09.figures/README.md) |
| query from code | [`05.graph/`](results/snapshots/2026-10-16/05.graph/README.md): node and edge tables, `graph.json`, and `cteph_knowledge_graph.sqlite` |
| see a summary of the build and the results | [`deck/cteph_knowledge_graph.pptx`](deck/cteph_knowledge_graph.pptx) (18 slides; [PDF preview](deck/cteph_knowledge_graph.preview.pdf)) |

`results/latest.json` points to the newest snapshot and to each of its products. Its paths are relative to `results/`.

`scripts/kg.py` queries the SQLite file of a snapshot, from the shell or from Python. It needs only the Python standard library.

```bash
python3 scripts/kg.py gene FGG                 # level, function, every finding with its study and sentence, conditions
python3 scripts/kg.py condition "Pulmonary embolism"
python3 scripts/kg.py region ABO
python3 scripts/kg.py publication 20677013
python3 scripts/kg.py studies --place Japan
python3 scripts/kg.py search fibrinogen
python3 scripts/kg.py neighbors HLA-B --edge reports
python3 scripts/kg.py sql "SELECT label, cteph_level FROM node_gene WHERE cteph_level <= 3"
python3 scripts/kg.py --json gene JAK2          # JSON instead of Markdown
```

In Python, with `scripts/` on the path: `from kg import KG; KG().gene('FGG')`.

## How it is built

The pipeline is one line of stages. Each stage writes one numbered folder of the snapshot.

| Stage | Folder | Steps |
|---|---|---|
| sources | `00.raw`, `00.pubmed` | Open Targets (CTEPH; each related condition); PubMed (the genetic-study search, the seed); Europe PMC Methods texts; the PubMed record of every publication the graph uses |
| evidence | `01.conditions`, `02.literature` | the related conditions and their evidence; the genetic studies of CTEPH, extracted by an isolated language model, with every item checked against its source text |
| curation | `03.curated` | publications left out (retracted, conference material) or flagged; the place and design class of each study; one name per variant; what each finding is about; text-mined links mapped to the wrong gene removed. Deterministic rules ([`config/curation.yaml`](config/curation.yaml)); every rule applied is logged |
| tables | `04.tables` | the long tables and the evidence ladder |
| graph | `05.graph` | typed nodes and edges, `graph.json`, a SQLite database, the schema |
| products | `06`–`09`, top level | grounded summaries, changes since the previous snapshot, the knowledge base, figures; the workbook; the viewer; the index (`README.md`) |

Language-model answers (screening, extraction and summaries) are cached in `cache/llm/`. A later snapshot sends only new PubMed records to the model.

Every choice has an entry in [`config/choices.yaml`](config/choices.yaml) with:
- its value, read from the run;
- its basis type: source, data check, design or user;
- its basis, with references.

Sheet 11 of the workbook and section 11 of the knowledge base print the choices with each snapshot, and `verify.sh` checks them. The methods are in [docs/METHODS.md](docs/METHODS.md).

## Requirements

- Nextflow (DSL2, strict syntax; tested with 26.04.4).
- Python 3 (tested with 3.12) with the packages in [`requirements.txt`](requirements.txt) (pandas, numpy, matplotlib, PyYAML, XlsxWriter, openpyxl).
- Internet access from the launch host for:
  - the Open Targets GraphQL API;
  - NCBI E-utilities;
  - the Europe PMC REST API.

  Every process runs locally on the launch host.
- The `claude` CLI (Claude Code; tested with 2.1.220). It is used only for records not yet in `cache/llm/`. `--LlmOffline true` never calls it.
- For the deck only: Node.js with `pptxgenjs`, and LibreOffice for the PDF preview.

The site paths of the tools are the defaults in the site configuration at the top of `cteph_knowledge_graph.nf`. Override them on the command line:
- `--python <python3 with the packages>`;
- `--claude <claude CLI>`.

## Run

```bash
git clone https://github.com/ZhaoTie-re/cteph_knowledge_graph.git && cd cteph_knowledge_graph
./verify.sh                                     # pre-flight, sections 1-5 (NO_LLM=1 skips the one model probe)
nextflow run cteph_knowledge_graph.nf --FetchDate $(date +%F) -ansi-log false
./verify.sh                                     # post-run, sections 6-10, on the newest snapshot
```

- **`--FetchDate`** names the snapshot. A new date queries Open Targets and PubMed again, writes `results/snapshots/<date>/`, compares it with the previous snapshot (`07.changes`), and moves `results/latest.json` to it.
- **`--Rebuild true --FetchDate <date>`** re-derives an existing snapshot from its saved sources (00–02) after a change to the curation, tables or products. Nothing is fetched again, except the PubMed publication records once if they are missing, and no model is called.
- **`--LlmOffline true`** uses cached model answers only.
- **A smoke run** adds `--LlmLimit 10 --MaxThemes 2 --results_root <scratch>`.
- **`nextflow lint cteph_knowledge_graph.nf nextflow.config`** is the safe syntax check. `nextflow run --help` would launch the pipeline.
- **On a cluster whose compute nodes have no internet access** (as on the one where this was built), launch from the login node. `verify.sh` puts the site's Nextflow environment on `PATH` when it exists; set `NXF_BIN=<dir>` to use another one.

## Update the configuration

| File | What it holds | Notes |
|---|---|---|
| [`config/curation.yaml`](config/curation.yaml) | exclusion rules, variant names, findings not about CTEPH, text-mining corrections, design classes, country names | An entry that matches nothing stops the build. A new finding with outcome words is logged for review in `03.curated/curation_log.tsv`. |
| [`config/related_conditions.yaml`](config/related_conditions.yaml) | the curated clinical relations | Each item needs a PMID, a verbatim quote from its title or abstract, and a one-line statement. The build stops if a quote is not in the abstract, or if it cites a publication the curation leaves out. |
| [`config/gene_themes.yaml`](config/gene_themes.yaml) | the functions (gene themes) | Each theme needs a basis and references. |
| [`config/cteph_gwas.yaml`](config/cteph_gwas.yaml) | the published CTEPH regions | Copied from Table 1 of the publication. |
| [`config/choices.yaml`](config/choices.yaml) | every choice and its basis | A new parameter or rule needs an entry, or `verify.sh` fails. |

## Repository layout

| Path | Content |
|---|---|
| `cteph_knowledge_graph.nf`, `nextflow.config` | the pipeline and its execution profile |
| `verify.sh` | pre-flight and post-run checks (read-only) |
| `config/` | the curated inputs and rules (above) |
| `scripts/` | the step scripts (below), the query tool `kg.py`, and `plot_style.py` (the figure style, vendored) |
| `docs/` | [METHODS.md](docs/METHODS.md) (how each item is derived), [OUTPUTS.md](docs/OUTPUTS.md) (every file of a snapshot), [STUDY_NOTES.md](docs/STUDY_NOTES.md) (design decisions and their history) |
| `deck/` | the slides on the build and the current snapshot, generated from the snapshot ([deck/README.md](deck/README.md)) |
| `results/` | `latest.json` and `snapshots/<FetchDate>/` |
| `cache/llm/` | the cached language-model answers, kept across snapshots |

The scripts, by stage:

| Stage | Scripts |
|---|---|
| sources | `fetch.py` (Open Targets; Ensembl for lead variants), `pubmed.py` (PubMed, Europe PMC), `publications.py` (the PubMed record of every publication) |
| evidence | `resolve.py` (conditions and evidence); `extract_genetics.py` and `llm.py` (literature screening and extraction) |
| curation | `curate.py` |
| tables and graph | `tidy_graph.py`, `build_graph.py` |
| products | `summarise.py`, `diff_snapshots.py`, `draw_figures.py`, `build_knowledge.py`, `write_table.py`, `build_viewer.py`, `manifest.py` |
| shared | `ot_common.py` |

## Data sources and reuse

- **Open Targets Platform**, release 26.09. The data are open (CC0 1.0).
- **PubMed** (NCBI E-utilities) and **Europe PMC**. Abstracts and full texts remain under their publishers' copyright or their articles' licenses. The graph keeps only short source sentences and verbatim quotes as evidence, each linked to its PubMed record.
- **The published CTEPH GWAS regions** are from Liley, J., Newnham, M., Bleda, M., et al. (2024). Shared and distinct genomics of chronic thromboembolic pulmonary hypertension and pulmonary embolism. *American Journal of Respiratory and Critical Care Medicine, 209*(12), 1477–1485.

## Author

Zhao Tie ([@ZhaoTie-re](https://github.com/ZhaoTie-re)).
