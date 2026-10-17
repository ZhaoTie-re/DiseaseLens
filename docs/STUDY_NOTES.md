# Study notes

## Design decisions (2026-10-16, with the user)

**History.** The first build, `cteph_opentargets` (a component of the CTEPH sequencing project), took everything Open Targets holds on CTEPH: targets, drugs, mouse models, publications. The user redirected it. The graph is to center on:
- CTEPH and the conditions related to it, every relation with its evidence;
- the genetic studies of both, with their methods;
- other important content.

The first snapshot of that build was deleted on 2026-10-16 at the user's request.

**Relatedness evidence.**
- The sources are Open Targets data plus PubMed-verified literature, kept apart by type: clinical, ontology and colocalization.
- Open Targets literature "co-occurrence" counts are unusable: `literatureOcurrences(additionalIds)` returns the sum of the two entities' counts, not their overlap.

**Categories.** All four, plus the data-driven shared-signal traits:
- precursor / spectrum;
- risk factors and comorbidities;
- intermediate traits;
- PH group.

**Genetic evidence.**
- GWAS credible sets with L2G;
- curated Mendelian sources;
- the CTEPH genetic literature (PubMed), since Open Targets has a single CTEPH GWAS (Liley 2024).

**Methods.** The user asked that the methods of the genetic studies be summarized too, so they are extracted per study, from the abstract and the open-access Methods sections.

**Other content.**
- East Asian / Japanese evidence is flagged throughout and gathered in one section and one sheet.
- Text-mined CTEPH targets are kept with their source sentences and summarized by Reactome theme.
- Drugs and mouse models are dropped.

**Language model.** Summaries were requested, so the screening, extraction and summaries use an isolated, cached model. Everything else is deterministic.

**External only.** Nothing of the CTEPH sequencing project's own results is read or compared.

## What the seed research showed (PubMed, 2026-10-16)
- **Control choice.** Most classic risk factors were found against non-thromboembolic PH controls (Bonderman 2005, 2009; Jaïs 2005; Lang 2013). Against PE survivors (Dodson 2022), malignancy, thyroid replacement, intravascular devices, IBD and osteomyelitis were not enriched; only splenectomy and non-O blood group trended. The seed keeps these contrary results.
- **IBD and osteomyelitis** rest mainly on one Vienna case-control study (OR 67, 95% CI 7.9–8832).
- **Platelets.** The evidence conflicts: a review lists high platelet counts, while measured cohorts found lower counts with more activation.
- **FXI** has genetic evidence only (the F11 locus).
- **Lp(a)** rests on one small 1998 study.
- **Hereditary spherocytosis** has no supporting abstract and is not in the seed.
- **The 2022 ESC/ERS guidelines** have no PubMed abstract, so a focused review of them (Rosenkranz 2023) is quoted instead.
- **HLA naming.** Older abstracts write HLA-B*5201 / DPB1*0202 (now B*52:01 / DPB1*02:02). PubMed renders some as "B(*)5201", so the grounding check ignores parentheses.

## What the build showed
- **Colocalization.** It is eCAVIAR only (PICS fine-mapping; no coloc H4). The ABO locus colocalizes with many unrelated traits (pancreatic cancer, type 2 diabetes, protein levels), so relations are graded by the number of CTEPH loci shared.
- **PubMed records.** The genetics query returns book records (GeneReviews) as `PubmedBookArticle`; they are parsed too.
- **Extraction smoke test.**
  - 6 records.
  - **Correct:** HLA-B*5201 / DPB1*0202, NFKBIL1, MUC6 and the CHIP genes were extracted, with their sample sizes and effects.
  - **Dropped by grounding:** an HLA-DRB1 the abstract does not name, and a "Mutations of factor II" variant that the abstract does not write.
  - **From the full text:** WES / MuTect2 / Holm / HLA-HD methods were read from the open-access Methods section of a study.
- **Summary grounding.** A gene check against the theme's target list was too strict: sentences name other genes. The check is now that a gene named must occur in the source text.

## History
- 2026-10-16: redesigned from `cteph_opentargets`; seed researched and verified; component built. Then:
  - At the user's choice, a deterministic pre-screen was added: 471 → 98 records, all 16 known genetic studies kept. Open Targets has no literature summaries, and its text mining misses and mis-maps the genetic studies.
  - The full-text fetch failed on PMC1729855, an HTTP 500 from Europe PMC. Such failures now count as "no full text".
  - At the user's request, the methods of studies without a full text now come from their abstract.
  - The user found the 16-sheet workbook "still cluttered, without logic, not intuitive" and chose a question-led workbook whose rows open level by level, with the CTEPH literature grouped by gene / pathway theme. It now has a guide, seven question sheets (answer in brief, then group › record › source as Excel row outlines; quotes in cell notes; links between sheets), the changes, and six flat appendix tables with filters. The knowledge base follows the same order.
  - In the same change, "Mendelian" evidence is shown as **curated gene evidence** (it includes gene burden and UniProt), and the ClinVar review status, which Open Targets returns as the record's confidence, is kept apart from the classifications of the other sources.
  - The first build's snapshot (`cteph_opentargets`) was deleted at the user's request.
- 2026-10-17: the user asked for a more professional and intuitive knowledge graph and deck, and that every item and value show its basis. Then:
  - The user found that the CTEPH GWAS regions were wrong: the pipeline named them by the top locus-to-gene gene of the three Open Targets credible sets (FGA at the FGG region) and missed the other six regions of Liley 2024, Table 1. The regions now come from the publication (`config/cteph_gwas.yaml`), named as there, with their Open Targets status beside them. Open Targets holds 7 of the 9 lead variants and a CTEPH credible set at 3.
  - At the user's choice, a weaker kind of shared signal was added for the regions without a CTEPH credible set: lead-variant overlap (a credible set of the condition contains the CTEPH lead variant, PIP ≥ ClppMin).
  - At the user's request, the four regions significant only after leveraging PE are marked G*, apart from G.
  - The CTEPH genes are now organized by evidence strength (a ladder of seven levels) and, as the second axis, function (gene themes, each with its basis and references). The user asked to see both; strength is the primary axis.
  - Every choice now has an entry in `config/choices.yaml` with its value (read from the run), basis type and basis; statements attributed to a source were checked against the full text or documentation. Corrections found in that check: Open Targets' documentation states no colocalization thresholds (the Open Targets preprint does); the curated threshold applies to the harmonic-sum datasource score; STREGA does not cover sequencing or HLA typing; AQP1 is limited and KDR definitive for PAH (Welch 2023).
  - The workbook now follows one line of argument (CTEPH genes, related conditions, shared regions, shared genes, context, choices); the knowledge base follows the same order and embeds the pipeline's figures (`07.figures`).
  - The user kept the gene × condition map as a heatmap, with how each cell is derived printed beside it, and asked for fewer asterisks: the East Asian star markers were dropped from the figures (East Asian evidence has its own sheet and section) and generated summaries lose their bold markup in the knowledge base.
  - `diff_snapshots.py` looked for the earlier snapshot's tables in `01.tables` (they are in `03.tables`), so no earlier snapshot was ever found; fixed.
  - The user asked why the gene map was incomplete. It had three filters that were not stated: CTEPH genes with a positive result only (dropping the 12 genes at levels 6–7, among them F8, which has GWAS evidence in VTE, DVT and PE), genes shared with a related condition only (dropping the 19 CTEPH-specific genes, among them HLA-DRA and MYH7B of the GWAS and HLA-B / HLA-DPB1 of the Japanese studies), and, in the deck, levels 1–4 only. The map (figure, sheet 5, knowledge base) now lists all 68 genes studied in CTEPH by level, with the number of related conditions per gene; a CTEPH-specific gene is an empty row, a finding in itself. The traits related by a shared signal only have no gene columns by design, and the figure now says so.

- 2026-10-17 (later): the user deleted the deck and asked for the knowledge graph to be curated and rebuilt, without new language-model calls: conference abstracts and retracted articles left out, the identified content organized, the structure and its presentation made clear, intuitive, professional and traceable. Then:
  - A curation step (`curate.py`, rules in `config/curation.yaml`) now sits between the extraction and the tables, fed by the PubMed record of every publication the graph uses (`publications.py`). It left out 4 of 2,151 publications (two retracted articles, two conference records; two of them had been included as genetic studies). At the user's request they are not listed in any product, and the products state the current snapshot without comparing it with earlier builds.
  - Reading the 51 extracted studies showed what needed organizing. The place of 28 studies was blank or "not stated" and 8 Japanese studies (Chen 2010, Momoi 2024, ...) were not flagged; the place is now the country of the first author's institution (PubMed affiliation; Europe PMC for older records): 47 of 49 placed, 15 in Japan. One variant was written up to five ways (factor V Leiden, FVL, 1691G>A, "factor V leidan", rs6025); 45 findings now carry one name and, for the classic variants, the dbSNP rsID.
  - 11 findings are associations with an outcome or phenotype within CTEPH (survival, clinical worsening, exercise capacity, DVT history), not with CTEPH itself. They counted toward the ladder: ACE and GNB3 were at level 3 on outcome associations alone. Levels 3–6 now count only findings about CTEPH itself (ACE 3 → 6, GNB3 3 → 7, BMPR2 4 → 5); the outcome findings stay, marked. A background sentence (BMPR2 in heritable PAH), a functional experiment (miR-759) and a Mendelian randomization with CTEPH as the exposure are marked too.
  - Europe PMC text mining mapped words to the wrong gene in the CTEPH literature (30 gene–word pairs); the largest was sGC (soluble guanylate cyclase, the riociguat target) mapped to SGCB, sarcoglycan β, in 51 publications. 103 text-mined links were dropped; 382 targets keep their literature (409 before).
  - The snapshot was re-derived from its saved sources with a new `--Rebuild` mode (nothing fetched again, no model call): 16 of the 29 generated summaries were reused from the cache, 13 are reported as not generated because the curation changed their inputs.
  - Structure: the folders follow the pipeline (00 sources … 03 curated, 04 tables, 05 graph, 06–09 products), with an index, the workbook and a new viewer at the top level. The graph is now typed (7 node and 12 edge types, schema in `05.graph/README.md`) and also written as a SQLite database, queried with `scripts/kg.py`. The workbook and the knowledge base follow four parts (A CTEPH genetics, B related conditions, C context, D provenance); the study methods moved next to the gene evidence (sheet 3), and a curation sheet (10) lists every exclusion and correction. For traceability, gene symbols and condition identifiers link to Open Targets, every finding and quote to PubMed, and the guide lists each source with its version and date.
  - A deck on the build and the current snapshot (`deck/`, 18 slides) is generated from the snapshot; it states the current state only and lists no excluded publication. All text of the repository (code, configuration, documentation, products, deck) is in American English; verbatim quotes, publication titles, ontology labels and the extracted methods text stay as written.

- 2026-10-17 (last): the component left the author's CTEPH sequencing project for its own repository, with its code, configuration, model cache and snapshot unchanged and nothing re-run. To run from any location, the pipeline now locates itself (`projectDir`), the figure style is vendored (`scripts/plot_style.py`), `results/latest.json` holds paths relative to `results/`, and `verify.sh` runs from its own folder. The Nextflow work directory was removed: the snapshot is complete, and `--Rebuild` re-derives it from its saved sources (00–02).
