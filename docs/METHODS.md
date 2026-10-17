# Methods

The methods of DiseaseLens, described for its first application, CTEPH. The steps are the same for any disease configured as in the README ("Configure a disease"); where CTEPH appears below, read the configured disease.

Every threshold, list and grouping named here is also an entry of `config/choices.yaml`, the choices registry, with its value in the run and its basis (section 8). The knowledge base (section 11) and the workbook (sheet 11) print the registry with each snapshot.

The pipeline runs in one line of stages, each writing one numbered folder of the snapshot: sources (00), related conditions (01), the extracted genetic studies (02), curation (03), long tables (04), the graph (05), then the products: summaries (06), changes (07), knowledge base (08), figures (09), the workbook, the viewer and the index.

## 1. Sources
- Open Targets Platform. The GraphQL API, `https://api.platform.opentargets.org/api/v4/graphql`. The release (data and API version) is recorded in every snapshot.
- PubMed. NCBI E-utilities (esearch, efetch), without a key; requests are spaced 0.4 s apart.
- Europe PMC. REST `fullTextXML`, for the Methods sections of open-access full texts; REST search, for the first-author affiliation of the older records that PubMed holds without one.
- PubMed records of every publication the graph uses (`scripts/publications.py`): publication types, journal issue and pagination, whether PubMed holds an abstract, the comments and corrections linked to the record (RetractionIn, ErratumIn, ...), and the first author's affiliation. Read once per snapshot and saved (`00.pubmed/publications.json.gz`); the curation (section 6) rests on them.
- Ensembl REST (GRCh38). Only to place the published CTEPH lead variants that the Open Targets search does not index (section 2).
- The curated seed. `config/related_conditions.yaml` holds the clinical relations, each with a verbatim quote from a PubMed abstract: 27 conditions with 76 evidence items, and 6 facts about CTEPH genetics, from 56 PMIDs. Every quote is checked word for word against the efetch record at each build.
- The published CTEPH GWAS regions. `config/cteph_gwas.yaml` copies Table 1 of Liley et al. 2024 (PMID 38470220; PMC11208965), with quotes for the tier definitions and for V.
- Where it runs. Every query runs from the launch host, because compute nodes have no internet. Requests are retried with exponential backoff on server errors, HTTP 429 and network failures.

## 2. The genetics of CTEPH: published regions, codes and the evidence ladder
### The GWAS regions (`fetch.py --published`, `resolve.py`, `tidy_graph.py`)
The only CTEPH GWAS (Liley 2024) reports 9 genome-wide significant regions. Open Targets holds no summary statistics for it and fine-maps only the curated top hits it can, so the regions are taken from the publication, not from Open Targets:
- G, five regions significant on the CTEPH meta-analysis P value (P < 5 × 10⁻⁸; tiers 1–3): ABO, FGG, F11, MYH7B, HLA-DRA.
- G*, four regions significant only in the co-analysis with pulmonary embolism (conditional false discovery rate V < 5 × 10⁻⁸, tier 3): F2, TSPAN15, SLC44A2, F5/NME7. They are marked apart because part of their signal is that of PE (the user, 2026-10-17).

Each region keeps the gene the publication names. The locus-to-gene genes of Open Targets are shown beside it, not used to name it: at the FGG region the top L2G gene is FGA (0.44), below the score this project accepts, while the publication names FGG. (The first build named the region FGA from L2G; the user found the error against Table 1, 2026-10-17.)

For each lead variant, the pipeline looks for:
1. its Open Targets variant ID: the Open Targets search by rsID, else the Ensembl GRCh38 alleles checked against Open Targets;
2. a CTEPH credible set of the Liley study containing it;
3. the credible sets of other GWAS that contain it, with its PIP in each.

On release 26.09, Open Targets holds 7 of the 9 lead variants and has a CTEPH credible set at 3 regions (ABO, FGG, F11); the lead variants of HLA-DRA (rs17202899) and F5/NME7 (rs796548658) are not in Open Targets.

### The evidence ladder (`ot_common.gene_ladder`)
Every gene with genetic evidence for CTEPH takes its strongest level:
1. G: genome-wide significant in the CTEPH GWAS;
2. G*: genome-wide significant only after leveraging PE;
3. associated in two or more genetic studies of CTEPH;
4. associated in one study;
5. a variant identified, without an association test (case reports, sequencing series);
6. tested, not associated;
7. other result only.

Levels 1–2 follow the GWAS. Level 3 ranks above level 4 because replication is what makes an association credible rather than a chance finding (Kraft 2009). Studies are counted, not independent samples: some share cases.

Levels 3–6 count only findings about CTEPH itself (section 6): an association with an outcome or phenotype within CTEPH (survival, clinical worsening, exercise capacity, the thrombosis history of patients) modifies the disease but is not evidence that the gene raises its risk. Such findings are kept and shown, and a gene with nothing else stays at level 7 (GNB3; ACE, tested and not associated with CTEPH itself, is at level 6).

### Function (`config/gene_themes.yaml`)
The second axis of the ladder is the function of the gene: eight gene themes, each with its basis and references (fibrinogen; ABO, VWF, FVIII, FXI; thrombophilia and hemostasis; HLA and immunity; clonal hematopoiesis, MPN and platelets; PAH genes; endothelium and vascular tone; pharmacogenetics). For example, the PAH genes are those ClinGen classified as definitive or moderate for PAH (Welch 2023), and the clonal genes those most frequently mutated in clonal hematopoiesis (Jaiswal 2014; Genovese 2014). A gene belongs to the first theme that lists it. A finding without a gene symbol joins the first theme whose genes or words occur in its gene text. The rest are listed under "Other genes" (logged by the build, to be placed in the config) or "No gene named".

## 3. Which conditions are related, and why (`scripts/resolve.py`)
A condition enters the graph through at least one kind of evidence. Every relation is kept as an evidence row with its source.

1. Clinical. The seed lists guidelines, reviews, meta-analyses and cohort or case-control studies whose abstracts state a relation to CTEPH. The categories follow the 2022 ESC/ERS guidelines and the ERS statement on CTEPH:
   - precursor / spectrum;
   - risk factor;
   - comorbidity;
   - intermediate trait;
   - PH group.

   Studies that found no association are kept, labeled with their direction. At every build, each quote must occur in the abstract fetched that day (normalized Unicode, dashes and spaces), or the build stops. Evidence from East Asian populations is flagged.
2. Ontology. From the Open Targets disease ontology (EFO / MONDO):
   - the parent of CTEPH (pulmonary hypertension);
   - its siblings;
   - the children of pulmonary arterial hypertension.
3. Colocalization, computed by Open Targets and selected here. The pipeline runs no colocalization itself.
   - **What Open Targets provides.** Open Targets fine-maps GWAS credible sets and, for each pair of overlapping credible sets, reports the coloc posterior H4 and the eCAVIAR CLPP. In release 26.09:
     - the three CTEPH credible sets (*ABO*, *FGG*, *F11*) are fine-mapped with PICS;
     - every colocalization carries the method label `COLOC_PIP_ECAVIAR`, with both values.
   - **What the pipeline does.** `fetch.py` downloads every colocalization of the CTEPH credible sets. `resolve.py` keeps a condition when one of its GWAS credible sets reaches CLPP ≥ 0.01 or H4 ≥ 0.8 with a CTEPH credible set. The thresholds and their basis are in `config/choices.yaml`.
   - **In this snapshot,** 307 of the 308 rows of `01.conditions/colocalisation.tsv` pass both thresholds, and 1 passes H4 only.
4. Lead-variant overlap (weaker), computed by this pipeline from the Open Targets credible sets. At a region whose lead variant Open Targets holds without a CTEPH credible set, a GWAS credible set of the condition contains that lead variant with PIP ≥ 0.01 (the same value as the CLPP threshold). Basis: if the CTEPH lead variant were the causal variant, the CLPP of the pair would equal its PIP in the other trait's credible set; without a CTEPH credible set the CTEPH side is not fine-mapped, hence weaker. Overlap is never used at a region with a credible set (the user chose this weaker class, 2026-10-17).

GWAS signal counts. A GWAS signal is one credible set of one GWAS study of the condition that colocalizes with the CTEPH credible set (3) or contains the CTEPH lead variant (4). A study has one such credible set per region (2026-10-16: one exception, at SLC44A2), and every count is a number of distinct GWAS studies (`study_id`).
- **Per condition and region:** the dot size of the regions figure, and the hover text of sheet 5.
- **Per condition, across regions:** each study counted once; this is the count in the evidence of the condition (sheet 4, knowledge base section 4).
- **Per region:** the per-condition counts added up (knowledge base section 5). A study that Open Targets maps to several conditions counts once for each.

Studies of one trait from different consortia, biobanks or releases each count and can share participants, so a count measures consistency across studies, not independent replication.

A relation of the shared-signal kinds is graded by the number of CTEPH regions it shares, by either route. Multi-region relations are stronger: ABO is the region most often shared alone, and a trait sharing only ABO mostly reflects the pleiotropy of ABO.

Genetics mode. Each condition gets one:
- `full`: its GWAS and curated gene evidence are fetched;
- `coloc_only`: only the colocalizing or overlapping credible sets are kept (the shared-signal traits, and platelet count, whose 125 GWAS are too many);
- `none`: no Open Targets record (procedures such as splenectomy, and CTEPD), or too broad a phenotype (cancer).

## 4. Genetics of the related conditions (`scripts/fetch.py`)
GWAS. `studies(diseaseIds)` gives every GWAS study of the condition: its sample size, cases, controls, cohorts, ancestry (`ldPopulationStructure`) and whether summary statistics are available. For each study, every credible set is fetched:
- the lead variant;
- its own P (mantissa and exponent, which avoids float underflow below 10⁻³⁰⁸) and β;
- the fine-mapping method;
- the top three locus-to-gene (L2G) genes. A gene is linked to a condition when it is the top L2G gene with score ≥ 0.5 (Mountjoy 2021).

A study is East Asian when its EAS share is ≥ 0.5.

Curated / Mendelian evidence. The sources are ClinGen, Gene2Phenotype, Genomics England PanelApp, Orphanet, ClinVar (eva), UniProt and gene burden. For each source:
- the targets whose datasource score is ≥ 0.5 for the condition or its ontology descendants are listed (one page ordered by that source's score). The datasource score is Open Targets' harmonic sum of the source's evidence scores, normalized to 0–1 (dividing by about 1.645); a single evidence item therefore passes only with a score of at least about 0.82;
- every evidence row for those targets is then read, with its classification, clinical significance, specific disease and PMIDs. For ClinVar the record's confidence is its review status, so it is kept in its own column (`review_status`) and not mixed with the classifications of the other sources.

The workbook and the knowledge base call this curated gene evidence: it is mostly rare-variant, Mendelian evidence, but gene burden and UniProt are not curated gene–disease validity in the ClinGen sense. The tables keep the name `mendelian`.

CTEPH itself:
- every associated target;
- its text-mined (Europe PMC) evidence;
- the targets' Reactome top-level pathways;
- its GWAS credible sets with all the colocalizations Open Targets computed for them, and the published regions (section 2).

Completeness. Every list is checked against the count the API reports, and its IDs must be unique. Index/size pages are not stable where the sort key ties (2026-10-16: the 615 CTEPH targets paged by 100 differed from one full page in 10 targets). A list longer than its first page is therefore re-read as one page of its full length.

## 5. The genetic studies of CTEPH in the literature (`scripts/pubmed.py`, `extract_genetics.py`)
1. Search. PubMed, with the query in `params.PubmedQuery`: CTEPH in title or abstract AND a genetics term. Publication date up to the fetch date.
2. Pre-screen (deterministic, before any model call).
   - Records typed Review, Systematic Review, (Practice) Guideline, Editorial, Comment, News or Consensus Development Conference are excluded.
   - So are records whose title and abstract contain no genetic-study term: polymorphism, genotype, allele, mutation, sequencing, exome, genome-wide, GWAS, HLA, SNP, rsID, germline, somatic, clonal hematopoiesis / CHIP, heritability, Mendelian randomization, copy number, pathogenic / rare variant, missense, haplotype, linkage.
   - On 2026-10-16 this kept 98 of 471 records. Benchmark: the 15 original genetic studies the seed cites are all found, kept and included (re-checked at every build; workbook sheet 10).
   - Excluded records are listed in `screening.tsv` with the rule that excluded them (`stage = rule`). The user chose this over reading all records (2026-10-16): Open Targets has no literature summaries, and its text mining misses and mis-maps genetic studies.
3. Screening and extraction (pass A).
   - An isolated language model reads the remaining abstracts, five per call. It decides whether each record is an original genetic study (or a case report) of CTEPH.
   - For the included ones it extracts:
     - design;
     - population (country, ancestry, Japanese / East Asian);
     - cases and controls, and the control type;
     - genes (as written, plus the HGNC symbol), variants;
     - findings (result, effect, P, and the source sentence);
     - replication;
     - the methods the abstract states.
4. Methods (pass B). For every included study, the methods come from the Methods section of its open-access full text (Europe PMC). When no full text can be fetched, they come from its PubMed title and abstract instead (the user, 2026-10-16). Each field records its source. The fields follow the reporting items of STREGA where they apply (STREGA does not cover sequencing or HLA typing, which are added):
   - cohort and recruitment, diagnosis, controls, DNA source;
   - genotyping platform, sequencing, HLA typing, variant calling and quality control, imputation;
   - statistical test, covariates, multiple testing, replication, power, software.

   Each field comes with its source passage. Europe PMC answers some non-open-access PMCIDs with a persistent HTTP 500; that is treated as "no full text".
5. Grounding. Every gene (as written or its symbol), variant, finding quote and sample size must occur in the abstract. Every methods value needs a quote that occurs in the text it came from. "Japanese" or "East Asian" must be stated in the text. Items that fail are dropped and logged (`extract_log.tsv`). Official symbols are kept only when Open Targets knows them by exact name.
6. Isolation. Each call is `claude -p --safe-mode --tools "" --strict-mcp-config --no-session-persistence`, with auto-memory disabled, run from `/tmp/cteph_kg_llm`. No project memory, file or instruction reaches the model; `verify.sh` section 5 tests this.
7. Caching. Answers are cached by hash of (model, effort, system prompt, prompt) in `cache/llm/`, so later snapshots call the model only for new records.

## 6. Curation (`scripts/publications.py`, `scripts/curate.py`, `config/curation.yaml`)
Between what was fetched and extracted and what the graph holds, a deterministic step applies the rules of `config/curation.yaml` and logs every rule it applies, row by row (`03.curated/curation_log.tsv`). No rule calls a language model.
1. Publications left out. Every publication the graph uses (the genetic studies, the seed, the text-mined literature, the GWAS publications, the references of the curated gene evidence) is looked up in PubMed. Left out of every table:
   - retracted articles and retraction notices (publication type Retracted Publication or Retraction of Publication, or a RetractionIn link);
   - conference material (publication type Congress or Conference Proceedings, or a journal-supplement record without an abstract).

   The user decided both (2026-10-17), and that the publications left out are not listed in any product (they remain in `03.curated/publications.tsv` with their reason). A publication left out leaves the genetic studies and their findings, the text-mined links, the GWAS studies and the reference lists of the curated gene evidence; the seed may not cite one (the build stops). Preprints, errata and expressions of concern are kept, flagged.
2. Place of a study. The country of the first author's institution, from the PubMed affiliation (Europe PMC when PubMed has none), mapped by country names, aliases, US states and a few cities; the country the abstract states is kept beside it. The place sets the Japanese and East Asian flags, as does a Japanese or East Asian population stated in the abstract. A multinational study takes the place of its lead institution.
3. Design class. One class per study, the first that applies: case report; genome-wide (GWAS, whole-genome sequencing, Mendelian randomization); clonal hematopoiesis; exome or gene panel; HLA typing; candidate gene; other.
4. Variant names. The abstracts write one variant several ways (factor V Leiden, FVL, 1691G>A, "factor V leidan"; HLA-B*5201, B(*)5201). Patterns per gene give each one name and, for the classic variants, its dbSNP rsID; the name as written is kept. A variant written only as an rsID takes it as its name.
5. What a finding is about. CTEPH itself (the default); a phenotype or outcome within CTEPH; another trait (Mendelian randomization with CTEPH as the exposure); a functional experiment; a background statement (a sentence restating earlier knowledge). Each finding of the non-default kinds was read against its source sentence and is listed with that sentence fragment; an entry that matches no finding stops the build. A later finding with words of an outcome (survival, worsening, predictor, ...) is logged for review. Findings that name a group of genes are marked as gene sets.
6. Text-mined links. Europe PMC text mining maps words to genes by dictionary; a word that names something else in the CTEPH literature makes a false link (sGC, soluble guanylate cyclase, mapped to SGCB; PTPA and BPA, pulmonary angioplasty; aPL, antiphospholipid antibodies; TF, tissue factor, mapped to transferrin; CHIP, clonal hematopoiesis). Each matched word was read against its sentences (2026-10-17). A link whose every word is wrong is dropped, with its sentences; a wrong word is removed from a link that keeps another; a correction that matches no link stops the build.

## 7. Mechanisms and summaries (`scripts/summarise.py`)
- Mechanism themes. The text-mined CTEPH targets are grouped by Reactome top-level pathway. For each theme with at least 3 targets, a summary of at most 150 words is written. Its input is up to two sentences from the best-scored publications of each of the theme's 12 best literature targets.
- Methods and findings. Summaries of the methods and of the findings of the CTEPH genetic studies (at most 350 words each) are written from the extracted tables and their tallies.
- Inputs. The curated tables: a target enters a theme only with a curated text-mined publication.
- Checks. Every [PMID] a summary cites must be one of its inputs, and every gene it names must occur in its source text. A summary that fails is kept but withheld from the knowledge base and the workbook.
- Offline (`--LlmOffline`, implied by `--Rebuild`). The model is never called: a summary whose inputs are unchanged is read from the cache; one whose inputs changed (for example, because the curation removed a publication or a mis-mapped target) is reported as not generated, never invented.

## 8. The choices registry (`config/choices.yaml`)
Each item of the registry is one choice: what it decides, its value, its basis type and its basis, and where its result is used. Basis types:
- source: a publication or the documentation of a data source, with references; every statement was checked against the full text or documentation (2026-10-17);
- data check: a check on the snapshot's data, recomputed at every build (for example, that ABO is the region most often shared alone, and the literature benchmark of section 5);
- design: a decision of this project, with its reason;
- user: decided with the user, with the date.

Values that are parameters are read from the run itself (`run_params.json`, also in the run manifest), never typed into the registry. `verify.sh` checks that every item has a basis and resolvable references, that every parameter named is recorded by the run, and, after a run, that every data check holds.

## 9. Tables, graph, changes, figures, knowledge base, workbook, viewer
- `tidy_graph.py`. Writes the long tables (docs/OUTPUTS.md) from the curated literature, with the evidence ladder. Gene–condition genetic evidence: GWAS (for CTEPH the published region genes, G or G*; for the related conditions the top L2G gene ≥ 0.5), curated, and literature (a positive finding about CTEPH itself).
- `build_graph.py`. The typed graph (`05.graph`, schema in its README): nodes condition, gene, region, publication, gwas_study, theme, pathway; edges related_to, supports, has_genetic_evidence, reports, studies, region_of, names_gene, shares_region, in_theme, in_pathway, mentions, gwas_of. Each edge type is checked against the node types it joins. Written as node and edge tables (one per type), node-link `graph.json`, and a SQLite database that also holds every curated and long table (numeric columns typed).
- `build_viewer.py`. One HTML file, offline and without libraries: overview (counts, the evidence ladder, regions × conditions), lists of genes, conditions, regions and studies, the curation, and a page per node with its neighbors and every edge with its source link (PubMed, Open Targets, dbSNP).
- `manifest.py`. The run manifest, `results/latest.json` (with the path of each product, relative to `results/`) and the index of the snapshot (`README.md`).
- `scripts/kg.py`. Not a step: a query tool over the SQLite file (gene, condition, region, publication, studies, search, neighbors, sql), from the shell or from Python.
- `diff_snapshots.py`. Compares stable keys with the newest earlier snapshot; a table the earlier snapshot lacks is listed as not compared.
- Populations. A literature study's population is the country its abstract states, else its ancestry, else "not stated"; its place is the country of its first author's institution (section 6). It is Japanese when done in Japan or stating a Japanese population, and East Asian when it is Japanese, states an East Asian population, or was done in China, Korea, Taiwan, Hong Kong or Mongolia.
- `draw_figures.py`. The figures of the snapshot, in the figure style of the CTEPH project (`scripts/plot_style.py`, a vendored copy), with one catalog README (docs/OUTPUTS.md).
- `build_knowledge.py`. Writes the Markdown/JSON knowledge base in the order of the workbook (thirteen sections and the references, in four parts: A CTEPH genetics, B related conditions, C context, D provenance), with the figures embedded and generated summaries marked.
- `write_table.py`. Writes the workbook in the house style, organized by question in the same four parts (docs/OUTPUTS.md); gene symbols link to Open Targets, condition identifiers to their Open Targets page, findings and evidence to PubMed:
  - each question sheet opens with its answer in brief (counts computed from the tables, then the generated summary where there is one), and states the basis of the values it uses;
  - then comes a nested table, group › record › source, built with Excel row outlines: groups open, sources closed;
  - source passages sit in cell notes, and links join the same gene or condition across sheets;
  - the appendix holds the flat tables, with filters.

  The workbook is re-opened and checked before publication: sheets, rows, footnote letters, the reason for every "—", outline levels, every link between sheets, and no internal names in cells or notes.

## 10. Rebuild from the saved sources (`--Rebuild true`)
A snapshot can be re-derived from its own saved sources (00.raw, 00.pubmed, 01.conditions, 02.literature) after a change to the curation, the tables or the products, so that it stays the snapshot of its fetch date. Nothing is fetched or extracted again; the publication records are read from PubMed only if the snapshot has not saved them yet. No language model is called (`LlmOffline` is implied). The run manifest records the rebuild (`params.Rebuild`), and `verify.sh` checks that a rebuild fetched nothing.
