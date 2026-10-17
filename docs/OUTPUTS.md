# Outputs

`results/latest.json` points to the newest complete snapshot and to each of its products (`snapshot`, `dir`, `index`, `viewer`, `knowledge`, `graph`, `database`, `workbook`, `release`, all `snapshots`); the paths are relative to `results/`.

`results/snapshots/<FetchDate>/` holds one snapshot. Its folders are numbered in the order the pipeline writes them; the products a reader opens are at the top level.

| Path | Content |
|---|---|
| `README.md` | the index of the snapshot: where to start, the headline counts, what each folder holds |
| `cteph_knowledge_graph.html` | the viewer: one HTML file, offline, with the whole graph embedded |
| `cteph_knowledge_graph_<FetchDate>.xlsx`, `workbook.json` | the workbook (below); its layout and rows per sheet |
| `00.raw/` | Open Targets responses for CTEPH: `meta.json`, `disease.json`, `associations.json.gz`, `evidence.json.gz` (text-mined and GWAS), `targets.json.gz` (Reactome themes), `gwas.json.gz` (credible sets with colocalizations), `published_loci.json.gz` (the 9 published CTEPH regions), `fetch_cteph.json` (counts); `diseases/disease.<id>.json.gz`, one related condition each |
| `00.pubmed/` | `seed_pubmed.json.gz` (the records the seed cites), `pubmed_search.json.gz` (the genetic-study search), `fulltext.json.gz` (Methods texts), `publications.json.gz` (the PubMed record of every publication the graph uses: publication types, issue, pagination, abstract, comments and corrections, first-author affiliation; with `checked_at`) |
| `01.conditions/` | `conditions.tsv`, `relations.tsv`, `colocalisation.tsv`, `overlap.tsv`, `seed_check.tsv`, `fetch_list.tsv` |
| `02.literature/` | the extraction as the model returned it after grounding: `screening.tsv`, `studies.tsv`, `findings.tsv`, `methods.tsv`, `extract_log.tsv`, `extract.json` |
| `03.curated/` | the curation (below) |
| `04.tables/` | the long tables (below), `tidy.json` (counts) |
| `05.graph/` | the graph (below) |
| `06.summaries/` | `summaries.tsv` / `.json`: kind, theme, text, grounded, problems ("not generated" in an offline build whose inputs changed), inputs; model calls, cached and not generated |
| `07.changes/` | `changes.tsv` (category, change, key, item, detail), `changes.json` (with the tables not compared) |
| `08.knowledge/` | `cteph_knowledge_graph.md`, `.json`: scope; A CTEPH genes, genetic studies; B related conditions, shared regions, shared genes, condition genetics; C mechanisms, East Asian evidence; D curation, choices, changes; caveats; references. The figures are embedded from `../09.figures/` |
| `09.figures/` | the figures (PNG 300 dpi and PDF) and their catalog `README.md`, `figures.json` (below) |
| `_run_info/` | trace, report, timeline, DAG; `run_manifest.json` (parameters incl. `Rebuild` and `LlmOffline`, the md5 of the seed, gene themes, curation rules and choices, script md5s, release, counts, curation summary, graph counts, changes) |

## 03.curated

| Table | One row per |
|---|---|
| `publications.tsv` | publication the graph uses: status (kept, flagged, excluded, not checked), rule and reason, where it is used, the PubMed record fields, first-author affiliation and its country |
| `studies.tsv` | genetic study of CTEPH kept: the extracted fields, plus place (country of the first author's institution) and its source, Japanese / East Asian flags, design class, publication flag |
| `findings.tsv` | finding of a kept study: the extracted fields, plus the variant's name and rsID, what the finding is about (`association_with`, with a note), gene set, whether it counts for the ladder |
| `methods.tsv`, `screening.tsv` | as extracted, without the studies left out; a record left out by its status is marked in the screening (`stage = publication status`) |
| `literature.tsv` | text-mined CTEPH target × publication kept, with the words it was matched as (wrong words removed) and its best sentences |
| `textmining_review.tsv` | every text-mined link: matched words, words kept, decision (kept, dropped) and reason |
| `curation_log.tsv` | every rule applied: step, item, rule, before, after, basis |
| `curation.json` | the counts |

## 04.tables

| Table | One row per |
|---|---|
| `cteph_regions.tsv` | published CTEPH GWAS region: gene as published, lead variant, GRCh37 position, MAF, OR, P, P in PE, V, tier, code (meta = G, pe = G*), Open Targets variant ID and how it was found, CTEPH credible set, credible sets of other traits containing the lead variant |
| `cteph_genes.tsv` | gene with genetic evidence for CTEPH: evidence level (1–7) and code, function (gene theme), studies by result (findings about CTEPH itself), studies with a finding about an outcome within CTEPH, East Asian studies, GWAS region and its statistics, PMIDs |
| `conditions.tsv` | condition: category, identifiers, evidence sources, genetics mode, clinical items, CTEPH regions shared, GWAS studies (East Asian), credible sets, curated genes |
| `relations.tsv` | piece of relation evidence: type, relation, statement, verbatim quote, design, population, East Asian, PMID, citation, strength, CTEPH regions |
| `colocalisation.tsv`, `overlap.tsv` | CTEPH region × GWAS credible set of a condition sharing it (CLPP / H4; lead-variant PIP) |
| `gwas_studies.tsv`, `gwas_loci.tsv` | GWAS study of CTEPH or a condition (publications left out removed); credible set with its lead variant, P, fine-mapping and top three L2G genes |
| `mendelian.tsv`, `mendelian_records.tsv` | curated gene evidence: condition × gene × source; every evidence row (references left out removed) |
| `cteph_targets.tsv`, `literature.tsv` | CTEPH target (score, curated text-mined publications and those dropped, Reactome themes); target × publication with sentences |
| `cteph_studies.tsv`, `cteph_findings.tsv`, `cteph_methods.tsv`, `cteph_screening.tsv` | the curated literature tables (as in 03.curated) |
| `gene_condition.tsv` | gene × condition × evidence kind (GWAS, GWAS (PE-adjusted), Mendelian, Literature): records, details, East Asian |

## 05.graph

| File | Content |
|---|---|
| `README.md` | the schema: node and edge types, identifiers, meaning, counts, columns; the evidence ladder; example SQL; the tables of the SQLite file |
| `nodes_<type>.tsv`, `edges_<type>.tsv` | one table per node type (condition, gene, region, publication, gwas_study, theme, pathway) and edge type (related_to, supports, has_genetic_evidence, reports, studies, region_of, names_gene, shares_region, in_theme, in_pathway, mentions, gwas_of) |
| `nodes.tsv`, `edges.tsv` | every node (id, type, label, url) and edge (source, target, type) |
| `graph.json` | node-link JSON (`networkx.node_link_graph(json.load(f), edges="edges")`), without the source quotes |
| `schema.json` | the types, counts and the ladder |
| `cteph_knowledge_graph.sqlite` | the node and edge tables, `nodes` and `edges`, every curated (`curated_<name>`) and long table (`table_<name>`), and `meta`; numeric columns typed |

## 09.figures

| Figure | Shows |
|---|---|
| `evidence_grid` | the CTEPH genes: evidence level (rows, strongest first) × function (gene theme, columns) |
| `genetic_studies` | every genetic study of CTEPH by year and design class, colored by place, sized by cases; the genome-wide studies labeled |
| `literature_genes` | genes reported by two or more genetic studies of CTEPH: studies by result (including associations with an outcome within CTEPH), grouped by gene theme |
| `condition_map` | CTEPH at the center, the related conditions by category; line style = kind of evidence, node size = amount of evidence, red outline = two or more CTEPH regions shared |
| `shared_regions` | related conditions × the 9 published CTEPH regions: filled dot = colocalization (color = CLPP), open dot = lead-variant overlap, size = GWAS signals |
| `gene_map` | every gene studied genetically in CTEPH by evidence level × CTEPH and the related conditions with gene evidence (G / G* / C / L; – studied without a positive result) |

## The workbook

Organized by question in four parts. Each question sheet opens with its answer in brief, then a nested table (Excel row outlines): open a row with + at the left, or a whole level with the buttons 1 2 3. Every row links to its source: findings and evidence to PubMed, gene symbols and condition identifiers to Open Targets, genes and conditions to the same item on other sheets; cells with a note hold the source passage, or the variant as written.

| Part | Sheet | Question |
|---|---|---|
| | Guide | the questions by part with their answers in brief, how to read, the snapshot and its sources (versions, dates, links), caveats |
| A CTEPH genetics | 1 CTEPH genes at a glance | Which genes have genetic evidence for CTEPH, how strong is it, and what do they do? |
| | 2 CTEPH gene evidence | What does each study report for each CTEPH gene? (level › gene › finding: place, design, variant, result, what it is about, source sentence) |
| | 3 Genetic studies | Which genetic studies of CTEPH are there, and how were they done? (design class › study, its methods by step) |
| B Related conditions | 4 Related conditions | Which conditions are related to CTEPH, and on what evidence? |
| | 5 Shared GWAS regions | Which CTEPH GWAS regions do the related conditions share? |
| | 6 Shared genes | Which genes do CTEPH and its related conditions share? |
| | 7 Condition genetics | What are the genetics of each related condition? |
| C Context | 8 Mechanisms | Which mechanisms does the CTEPH literature describe? |
| | 9 East Asian | What evidence comes from East Asian and Japanese populations? |
| D Provenance | 10 Curation | What was left out or corrected, and why? |
| | 11 Choices and basis | On what basis was each item and value chosen? |
| | 12 Changes | What changed since the previous snapshot? |
| Appendix | A1 GWAS studies … A6 Screening | flat tables with filters: GWAS studies, GWAS loci, shared signals, curated genes, CTEPH studies (with place, design class, flag), the screening of every PubMed record |

## The viewer

`cteph_knowledge_graph.html` opens in any browser without a network. Tabs: Overview (counts; the evidence ladder; the CTEPH GWAS regions × the conditions sharing them), Genes, Related conditions, GWAS regions, Genetic studies (filter by place and design), Curation. Each node has a page (`#/n/<id>`, bookmarkable): its facts, a diagram of its neighbors, and one sortable, filterable table per edge type with every attribute and its source link (PubMed, Open Targets, dbSNP).
