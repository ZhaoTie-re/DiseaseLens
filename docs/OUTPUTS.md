# Outputs

Every disease has its own folder of results, `results/<key>/`, where `<key>` is the name of its folder in
`diseases/` (`cteph`, `hpah`). Results are written on the machine that runs the pipeline and are not part of the
repository (git ignores `results/` and `cache/`). Below, `<short>` is the disease's short name (CTEPH, HPAH) and `<date>` the fetch date
of a snapshot.

| Path | Content |
|---|---|
| `results/index.html`, `results/index.json`, `results/README.md` | every disease of the results folder: its newest snapshot, release and headline counts, with links to its viewer, workbooks and knowledge base (paths relative to `results/`) |
| `results/<key>/latest.json` | the disease's newest complete snapshot and each of its products (`snapshot`, `dir`, `index`, `viewer`, `knowledge`, `graph`, `database`, `workbook`, `release`, all `snapshots`, `disease`, `counts`); paths relative to `results/<key>/` |
| `results/<key>/snapshots/<date>/` | one snapshot (below) |
| `results/_runs/` | the trace, report and timeline of the last run (they hold paths of the machine that ran it) |

## A snapshot

Its folders are numbered in the order the pipeline writes them; the products a reader opens are at the top level.

| Path | Content |
|---|---|
| `README.md` | the index of the snapshot: where to start, the headline counts, what each folder holds |
| `<key>_knowledge_graph.html` | the viewer: one HTML file, offline, with the whole graph embedded |
| `<key>_knowledge_graph_<date>.xlsx`, `workbook.json` | the full workbook (below); its layout and rows per sheet |
| `<key>_knowledge_graph_<date>_summary.xlsx` | the summary workbook: the cover and the question sheets, without the appendix |
| `00.raw/` | Open Targets responses for the disease: `meta.json` (release), `disease.json` (record, synonyms, parents, children), `associations.json.gz`, `evidence.json.gz` (text-mined and GWAS), `targets.json.gz` (Reactome pathways), `gwas.json.gz` (credible sets with colocalizations), `mendelian.json.gz` (curated gene evidence), `published_loci.json.gz` (the published regions, when the profile has some), `fetch_core.json` (counts); `diseases/disease.<id>.json.gz`, one related condition each |
| `00.pubmed/` | `pubmed_search.json.gz` (the genetic-study search), `seed_pubmed.json.gz` (the records the clinical relations cite), `relations_search.json.gz` (the relation search, for a disease without a curated list), `fulltext.json.gz` (Methods texts), `publications.json.gz` (the PubMed record of every publication the graph uses, with `checked_at`) |
| `01.relations/` | for a disease without a curated list: `extracted_relations.yaml` (the relations read from PubMed, in the schema of a curated list), `relations_review.tsv` (disputed between the readers), `relations_log.tsv` (dropped by grounding), `relations.json` (counts and the readers) |
| `01.conditions/` | `conditions.tsv`, `relations.tsv`, `colocalisation.tsv`, `overlap.tsv`, `seed_check.tsv` (every quote checked), `ontology_log.tsv` (parents left out as too broad), `fetch_list.tsv` |
| `02.literature/` | the extraction after grounding and adjudication: `screening.tsv`, `studies.tsv`, `findings.tsv`, `methods.tsv` (each with `agreement`), `review.tsv` (disputed items and result checks), `gene_pathways.tsv`, `extract_log.tsv` (dropped items, by reader), `extract.json` (counts, agreement, the readers and the models that answered) |
| `03.curated/` | the curation (below) |
| `04.tables/` | the long tables (below), `ladder.json` (the evidence ladder of the snapshot), `themes.json` (the gene themes used), `themes_auto.yaml` (automatic themes, a draft to curate), `tidy.json` (counts) |
| `05.graph/` | the graph (below) |
| `06.summaries/` | `summaries.tsv` / `.json`: kind, theme, text, grounded, problems, the second reader's check (`check`), inputs; model calls, cached and not generated |
| `07.changes/` | `changes.tsv` (category, change, key, item, detail), `changes.json` (with the tables not compared and why) |
| `08.knowledge/` | `<key>_knowledge_graph.md`, `.json`: scope; A the genetics of the disease (genes, genetic studies); B related conditions, shared regions, shared genes, condition genetics; C mechanisms, East Asian evidence; D curation (with the readers' agreement), choices, changes; caveats; references. The figures are embedded from `../09.figures/` |
| `09.figures/` | the figures (PNG 300 dpi and PDF), their catalog `README.md` and `figures.json` (below) |
| `_run_info/` | `profile.json` (the resolved profile of the disease), `run_manifest.json` (parameters, the readers, script md5s, release, counts, curation summary, graph counts, changes) |

## 03.curated

| Table | One row per |
|---|---|
| `publications.tsv` | publication the graph uses: status (kept, flagged, excluded, not checked), rule and reason, where it is used, the PubMed record fields, first-author affiliation and its country |
| `studies.tsv` | genetic study kept: the extracted fields and the readers' agreement, plus place (country of the first author's institution) and its source, Japanese / East Asian flags, design class, publication flag |
| `findings.tsv` | finding of a kept study: the extracted fields, agreement and result check, plus the variant's name and rsID, what the finding is about (`association_with`, with a note), gene set, whether it counts for the ladder |
| `methods.tsv`, `screening.tsv` | as extracted, without the studies left out; a record left out by its status is marked in the screening (`stage = publication status`) |
| `relations.tsv` | the clinical relations, without those citing a publication left out |
| `literature.tsv` | text-mined target × publication kept, with the words it was matched as (wrong words removed) and its best sentences |
| `textmining_review.tsv` | every text-mined link: matched words, words kept, decision (kept, dropped) and reason |
| `model_review.tsv` | every item the two readers disagreed on: pass (abstract, full text, or relations for a disease whose relations are read from PubMed), item, each reader's answer, the final value, the adjudicator's reason (its "Reading 1" and "Reading 2" given as reader A and reader B); and the findings whose result disagrees with their sentence |
| `curation_log.tsv` | every rule applied: step, item, rule, before, after, basis |
| `curation.json` | the counts (the excluded publications counted only), the curated genes left out for their name, the disputed items |

## 04.tables

| Table | One row per |
|---|---|
| `disease_regions.tsv` | GWAS region of the disease: name, genes, lead variant, position (GRCh37 as published, or GRCh38), code (G, G*), source (published, Open Targets), the published statistics (MAF, OR, P, secondary P, V, tier) where given, P as mantissa and exponent (Open Targets), Open Targets variant ID and how it was found, the disease's credible set, locus-to-gene gene and score, the credible sets of other GWAS holding the lead variant of a region without a credible set of the disease (as in `overlap.tsv`) |
| `disease_genes.tsv` | gene with genetic evidence for the disease: evidence level (1–7) and code (G, C, G*), function (gene theme), studies by result (findings about the disease itself), studies with a finding about an outcome within it, East Asian studies, GWAS region and its statistics, curated sources, PMIDs |
| `conditions.tsv` | condition (`custom:<name>` without an Open Targets entry): category, identifiers, evidence sources, genetics mode and why (`genetics_note`), clinical items that relate it, regions shared, GWAS studies (East Asian), credible sets, curated genes (blank when its genetics are not fetched) |
| `relations.tsv` | piece of relation evidence: type, relation, statement, verbatim quote, design, population, East Asian (the evidence of this relation: for clinical evidence, a study done in East Asia or of an East Asian population; for a shared signal, both GWAS mostly East Asian; the one flag the map, viewer, workbook and knowledge base read), PMID, citation, strength, regions, agreement |
| `colocalisation.tsv`, `overlap.tsv` | region of the disease × GWAS credible set of a condition sharing it (CLPP / H4; lead-variant PIP) |
| `gwas_studies.tsv`, `gwas_loci.tsv` | GWAS study of the disease or a condition (publications left out removed); credible set with its lead variant, P, fine-mapping and top three L2G genes |
| `mendelian.tsv`, `mendelian_records.tsv` | curated gene evidence of the disease and of each condition: condition × gene × source; every evidence row (references left out removed) |
| `disease_targets.tsv`, `literature.tsv` | target of the disease (score, curated text-mined publications and those dropped, Reactome themes); target × publication with sentences |
| `disease_studies.tsv`, `disease_findings.tsv`, `disease_methods.tsv`, `disease_screening.tsv` | the curated literature tables (as in 03.curated) |
| `disease_signals_not_regions.tsv` | genome-wide significant credible set of the disease's GWAS in Open Targets that is no region, its lead not credible and not replicated: lead variant, study, P, MAF, expected copies of the minor allele in the smaller group, the reason |
| `tidy.json` | the counts of the snapshot: the GWAS of the related conditions (a study of two conditions once) apart from the disease's own, the signals that are no region, the conditions not related, the excluded publications (their IDs only, for the comparison with a later snapshot) |
| `gene_condition.tsv` | gene × condition × evidence kind (GWAS, the secondary label of G*, Mendelian, Literature): records, details, East Asian |

## 05.graph

| File | Content |
|---|---|
| `README.md` | the schema: node and edge types, identifiers, meaning, counts, columns; the evidence ladder; example SQL; the tables of the SQLite file |
| `nodes_<type>.tsv`, `edges_<type>.tsv` | one table per node type (condition, gene, region, publication, gwas_study, theme, pathway) and edge type (related_to, supports, has_genetic_evidence, reports, studies, region_of, names_gene, shares_region, in_theme, in_pathway, mentions, gwas_of) |
| `nodes.tsv`, `edges.tsv` | every node (id, type, label, url) and edge (source, target, type) |
| `graph.json` | node-link JSON (`networkx.node_link_graph(json.load(f), edges="edges")`), without the source quotes |
| `schema.json` | the types, counts, the ladder and its codes, the region labels, the disease (`graph`) |
| `<key>_knowledge_graph.sqlite` | the node and edge tables, `nodes` and `edges`, every curated (`curated_<name>`) and long table (`table_<name>`), and `meta`; numeric columns typed |

## 09.figures

A figure without data is not drawn (for example `shared_regions` for a disease without GWAS regions), and the
catalog says why.

| Figure | Shows |
|---|---|
| `evidence_grid` | the disease's genes: evidence level (rows, strongest first) × function (gene theme, columns) |
| `genetic_studies` | every genetic study of the disease by year and design class, colored by place, sized by cases; the genome-wide studies labeled |
| `literature_genes` | genes reported by two or more genetic studies: studies by result (including associations with an outcome within the disease), grouped by gene theme |
| `condition_map` | the disease at the center, the related conditions by category; line style = kind of evidence, node size = amount of evidence, red outline = two or more regions shared |
| `shared_regions` | related conditions × the disease's GWAS regions: filled dot = colocalization (color = CLPP), open dot = lead-variant overlap, size = GWAS signals |
| `gene_map` | every gene studied genetically in the disease by evidence level × the disease and the related conditions with gene evidence (G / G* / C / L; – studied without a positive result) |

## The workbook

The workbook is organized by question, in four parts. It opens in Microsoft Excel, WPS and LibreOffice and prints on
A4 landscape.

Each question sheet opens with its answer in brief, then a nested table (Excel row outlines): open a row with + at
the left, or a whole level with the buttons 1 2 3. Every row links to its source: findings and evidence to PubMed,
gene symbols and condition identifiers to Open Targets, and genes and conditions to the same item on other sheets.
Cells with a note hold the source passage, or the variant as written.

| Part | Sheet | Question |
|---|---|---|
| | Guide | a cover: the disease, the answers in brief with the key numbers, how to read, the sources with their versions and dates, the readers and their agreement, how to cite; then the contents |
| A Genetics of the disease | 1 `<short>` genes at a glance | Which genes have genetic evidence for the disease, how strong is it, and what do they do? |
| | 2 `<short>` gene evidence | What does each study report for each gene? (level › gene › finding: place, design, variant, result, what it is about, readers, source sentence) |
| | 3 Genetic studies | Which genetic studies are there, and how were they done? (design class › study, its methods by step) |
| B Related conditions | 4 Related conditions | Which conditions are related to the disease, and on what evidence? |
| | 5 Shared GWAS regions | Which GWAS regions of the disease do the related conditions share? |
| | 6 Shared genes | Which genes do the disease and its related conditions share? |
| | 7 Condition genetics | What are the genetics of each related condition? |
| C Context | 8 Mechanisms | Which mechanisms does the literature describe? |
| | 9 East Asian | What evidence comes from East Asian and Japanese populations? |
| D Provenance | 10 Curation | What was left out or corrected, and why; what did the readers disagree on, and how was it decided? |
| | 11 Choices and basis | On what basis was each item and value chosen? |
| | 12 Changes | What changed since the previous snapshot? |
| Appendix (full workbook only) | A1 GWAS studies … A6 Screening | flat tables with filters: GWAS studies, GWAS loci, shared signals, curated genes, the disease's studies, the screening of every PubMed record |

## The viewer

`<key>_knowledge_graph.html` opens in any browser without a network, and can be sent as a file.
- Tabs: Overview (counts; the evidence ladder; the regions × the conditions sharing them), Genes, Related conditions,
  GWAS regions, Genetic studies (filter by place and design), Curation, How it is built.
- How it is built, in three layers; the general rules, each followed by a line giving the disease in this snapshot:
  - at a glance: sources → methods → knowledge graph (each source with its release, date and counts; the rules, the
    two readers and the adjudicator, the source-text check; the nodes by type, the edges, the products);
  - step by step: the eight stages of the build (the disease, sources, related conditions, genetic studies,
    curation, evidence and tables, graph, summaries), each with the sources and methods it uses and its counts;
  - in detail (folded): the sources; the associations (every edge type, where it comes from and by which rule, with
    its count); categories and levels; the language models (the three models, a diagram of one reading from record
    to graph, where models are used and how each answer is checked, how often the two readers agreed, what no model
    decides, every disputed item with both readings, the final reading and the decision); every rule of the choices
    registry with its value and basis.
- Related conditions: the condition map, interactive (the figure `09.figures/condition_map`, drawn from its data
  `condition_map.json`): hover over a condition for its evidence and the East Asian evidence of its relation, click
  it to open its page,
  select a category (a button or its sector) to list its conditions in the table below; the selection is part of the
  address (`#/conditions/<category>`). Under the map, what it shows and what it leaves out, of every related
  condition; the rules that bring a condition into the graph and classify it are on the How it is built page.
- Each node has a page (`#/n/<id>`, bookmarkable): its facts, a diagram of its neighbors, and one sortable,
  filterable table per edge type, with every attribute and its source link (PubMed, Open Targets, dbSNP).
- Every identifier in a text is a link, the text unchanged: a GWAS study (GCST…, FINNGEN_…, UKB_PPP_…) to its page in
  the graph, else in Open Targets; an rsID to dbSNP; a PMID to its page in the graph, else PubMed. (A workbook cell
  holds one link, so the workbook links each row's source instead.)
- Any page prints, or saves as PDF, without the navigation.
