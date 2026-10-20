# Study notes

## Design decisions (2026-10-16, with the user)

**History.** The first build took everything Open Targets holds on CTEPH: targets, drugs, mouse models, publications. The user redirected it. The graph is to center on:
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

**External only.** Nothing of the user's own studies is read or compared.

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
- 2026-10-16: redesigned from the first build; seed researched and verified; component built. Then:
  - At the user's choice, a deterministic pre-screen was added: 471 → 98 records, all 16 known genetic studies kept. Open Targets has no literature summaries, and its text mining misses and mis-maps the genetic studies.
  - The full-text fetch failed on PMC1729855, an HTTP 500 from Europe PMC. Such failures now count as "no full text".
  - At the user's request, the methods of studies without a full text now come from their abstract.
  - The user found the 16-sheet workbook "still cluttered, without logic, not intuitive" and chose a question-led workbook whose rows open level by level, with the CTEPH literature grouped by gene / pathway theme. It now has a guide, seven question sheets (answer in brief, then group › record › source as Excel row outlines; quotes in cell notes; links between sheets), the changes, and six flat appendix tables with filters. The knowledge base follows the same order.
  - In the same change, "Mendelian" evidence is shown as **curated gene evidence** (it includes gene burden and UniProt), and the ClinVar review status, which Open Targets returns as the record's confidence, is kept apart from the classifications of the other sources.
  - The first build's snapshot was deleted at the user's request.
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

- 2026-10-17 (last): the component moved to its own repository, with its code, configuration, model cache and snapshot unchanged and nothing re-run. To run from any location, the pipeline now locates itself (`projectDir`), the figure style is vendored (`scripts/plot_style.py`), `results/latest.json` holds paths relative to `results/`, and `verify.sh` runs from its own folder. The Nextflow work directory was removed: the snapshot is complete, and `--Rebuild` re-derives it from its saved sources (00–02).

- 2026-10-18: renamed **DiseaseLens** (folder `DiseaseLens`; GitHub `ZhaoTie-re/DiseaseLens`, formerly `cteph_knowledge_graph`). The user repositioned it as a pipeline that reviews the genetics of any disease, with CTEPH as its first application. The pipeline script became `main.nf`; the disease is still configured for CTEPH (`--DiseaseCore`, `--PubmedQuery`, `config/`). Labels, output names and model prompts still name CTEPH; a disease profile that sets them is on the roadmap. The snapshot 2026-10-16 is unchanged apart from the "Generated by" line of its README, knowledge base and graph README.

- 2026-10-19: generalized to any disease, several per run, at the user's request ("for others to build the knowledge
  graph of their own disease"). Then:
  - **A disease is its Open Targets ID.** The user required the ID, never a name. `scripts/new_disease.py <ID>`
    writes `diseases/<key>/profile.yaml`, each value with the rule that chose it. Every rule starts from the ID and
    can be configured: the profile, optional curated files, the shared rules of `config/`, and the run's parameters.
    The search rules came from testing several diseases:
    - an acronym that also names a gene is not searched (PAH, phenylalanine hydroxylase; PPH1, an alias of BMPR2);
    - an exact synonym that names another disease is not searched ("pulmonary hypertension" for PAH);
    - an acronym is searched only beside a word of the name: HPAH alone also found human phenylalanine
      hydroxylase and halogenated polycyclic aromatic hydrocarbons, 57 of 331 records.
  - **Defaults where CTEPH had curated files.**
    - GWAS regions from Open Targets.
    - Curated gene–disease evidence of the disease itself (code C, at level 1 with G; CTEPH has none in Open Targets).
    - Reactome themes.
    - Clinical relations read from PubMed guidelines, meta-analyses, cohort and case-control studies.
    - Generic text-mining rules: a word that names the disease or an ontology neighbor is a false link, and a curated
      gene named like such a word is left out (PAH for HPAH).
    - An ontology parent with more than 30 children is too broad to give neighbors ("hereditary disease", 1,070
      children, for HPAH).
  - **Two independent model readers, at the user's request.**
    - Every reading is done twice, and a third model adjudicates where the readers disagree. Only agreed or
      adjudicated items count, and every dispute is listed.
    - The three models are named, not left to the CLI's default, after a CLI update changed the default from Opus 5 to
      Opus 5.5: reader A Opus 5, the model of the cached CTEPH answers; reader B Sonnet 5.5; adjudicator Opus 5.5.
    - The adjudicator reads the readings in an order that alternates by record.
    - A first HPAH test showed the adjudicator, then the model of reader A, siding with A on complications that the
      prompt excludes. Hence the third model, and a clearer rule that consequences of the disease are not related
      conditions.
  - **Privacy.**
    - Keys come from the environment only, and no path of the machine that ran is kept.
    - The command-line backend adds the account's e-mail address to what the model sees, so any e-mail address not in
      the source text is removed from every answer before it is cached.
    - Isolation is tested with a project file holding a random word that a model must not see.
  - **Several diseases per run** (`results/<key>/`, `results/index.html`). Products are made for sharing as files: a
    cover sheet, A4 printing, a summary workbook, an information bar and print style in the viewer.
  - **The generalized code was checked on the 2026-10-16 CTEPH snapshot:** every table, the graph, the figures and the
    workbook numbers were the same. At the user's choice that snapshot was then dropped, and CTEPH was built anew with
    HPAH on 2026-10-19.
  - **The viewer** (user, 2026-10-19): the information bar was removed (the header keeps the snapshot, release and
    counts; the workbook's cover keeps the readers and how to cite). The Related conditions page holds the condition
    map, interactive, above the table, which follows the selection. The rules are explained once, on a How it is built
    page of every viewer (sources → methods → graph at a glance, the build step by step, then the details), as
    general rules each followed by how they apply to the disease; the map keeps only its key and what it shows.
  - **East Asian evidence of a relation** (user, 2026-10-20): the flag of a relation, and the diamond of the map, say
    where the evidence that relates a condition to the disease comes from: clinical evidence from an East Asian
    population (named in the abstract, when the readers read it), or a shared signal whose two GWAS, the disease's and
    the condition's, have mostly East Asian samples. A condition's own East Asian GWAS do not make its relation East
    Asian (they are listed with its genetics); the ontology and a lead-variant overlap never do. An earlier draft
    counted the condition's own GWAS, and the shared-signal flag had never been set. Clinical evidence from a study
    done in East Asia counts as for the genetic studies (user, 2026-10-20: a Japanese cohort, PMID 28890480, whose
    abstract names no country had been flagged not East Asian, while the map, on its older rule, showed a diamond).
  - **Rules tightened after the first two-disease build** (user, 2026-10-19), for every disease:
    - Curated gene evidence (code C, and the curated genes of the related conditions): a ClinVar link counts only
      with a pathogenic or likely pathogenic record of two stars or more, and a gene needs an approved symbol. The
      curated-gene check of HPAH had found three level-1 genes resting on single-submitter or unreviewed ClinVar
      records (CAPNS1, NOP58, and a gene without a symbol).
    - Relations read from PubMed: not a sign, symptom, imaging or laboratory finding of the disease, a molecular or
      cellular process, an animal or cell experiment, or a comparison between groups of patients; the readers are
      compared on the category too. The first HPAH reading had kept CT findings of TBX4 carriers and mouse and
      cell experiments as related conditions.
    - Relations read from PubMed (user, 2026-10-20, while building PVOD/PCH): a result for a broader group counts
      only when the disease's patients are most of it or are given apart (a pediatric PH cohort with 9 suspected
      PVOD had lent its list of syndromes to PVOD/PCH and HPAH); a demographic attribute (sex, age, race or
      ethnicity) is not a related condition.
    - A finding of no association is not an association (user, 2026-10-20: IPAH's lipoprotein(a) level, inherited
      thrombophilia and liver disease severity had entered as related conditions on studies that found none).
      The direction stays a decisive item of the two readers and the adjudicator.
    - Conditions without an Open Targets entry are merged when their names differ only by use / exposure or an
      alias in parentheses (user, 2026-10-20: IPAH had anorexigen and fenfluramine in four and two names).
    - Open Targets evidence of a disease (user, 2026-10-20): the disease with its descendants, for targets,
      text-mined literature, GWAS and curated genes alike; and other entries a profile names (evidence_also),
      proposed by a new profile for a grouping name: PVOD/PCH (MONDO_0018554) has no literature or GWAS evidence of
      its own, while the literature of pulmonary veno-occlusive disease is mapped to MONDO_0009937.
    - Summaries that the second reader and the adjudicator find unsupported are withheld with the reason; the share
      withheld is reported, not held to a threshold.
    - GWAS regions from Open Targets need a credible lead variant (2026-10-20, from a review of the IPAH build): all
      fourteen IPAH regions came from one PheCode GWAS of a biobank (GCST90480139, 2,353 cases) with leads of
      frequency 0 to 0.5 percent, and eight of them gave genes the code G. A first fix left out every lead rarer than
      1 percent; the user corrected it ("rare variants can be evidence too, as long as they are credible"). A lead is
      now credible with at least 25 expected copies of its minor allele in the smaller group of its study (Neale lab
      2017), and, when rare, P < 5 × 10⁻⁹ (Pulit 2017); or when another GWAS of the disease replicates it. Leads that
      are not credible are listed with the reason. The IPAH signals had 0 to 23 expected copies in the cases.
  - **A review of the whole build for any disease** (user, 2026-10-20: "rethink the logic of the knowledge graph,
    check for errors and misunderstandings; it must work for any disease in Open Targets"). Four read-only reviews
    (sources and profile; model reading; conditions, curation and tables; graph and products) were checked against
    the four snapshots and the code; every confirmed finding was fixed as a general rule:
    - **Sources.** The disease's own name is always searched (a database form "A, B" also as "B A"); a short name
      made up from initials is never searched and is marked for checking; a synonym that names one of the
      disease's descendants or an entry whose evidence it takes is kept (PVOD lost "pulmonary capillary
      hemangiomatosis"); the names of those entries are searched (PVOD/PCH missed 30 records of "pulmonary
      venoocclusive disease"). PubMed is searched by entry date, so that a date finds the same records later.
      Open Targets lists longer than 3,000 rows (PAH has 4,549 targets) are read in pages until complete; an rsID
      must match exactly (the search also returns longer rsIDs). A first author without an affiliation is looked
      up in Europe PMC, not given another author's institution (16 records had been).
    - **Model reading.** A record's answer is cached on its own, so a record is read once whatever batch it falls
      in. Short gene names and population words are matched as whole words: "Chin" had matched "which in", so the
      check of a reader's East Asian flag passed for a reason that would pass a European study too (the study it
      was found in, from Peking University, is East Asian, by its place and by the city its abstract names; East
      Asian cities now count as a statement of place); "Asian" alone no longer counts. A methods field stated twice keeps both
      quotes (74 had been dropped unlogged). Findings without a gene are compared by variant; inclusion agreement
      counts only records disputed on inclusion; a malformed answer is no answer, never a value; a summary not
      checked by the second reader, or without a ruling on each disputed claim, is withheld. Abbreviations no
      longer merge conditions (PE: pulmonary embolism, pleural effusion); the relation prompt names no example
      consequences of pulmonary hypertension.
    - **Conditions and tables.** GWAS regions from Open Targets need a credible lead (above); a credible set near a
      published lead belongs to that region; two regions never share a name. The disease's own GWAS and the entries
      whose evidence it takes are never related conditions. A publication the curation leaves out removes its GWAS
      from the regions and shared signals, and each condition's sources and category are taken from the curated
      relations. The literature links of non-Japanese East Asian studies are flagged East Asian (28 of 67 IPAH
      links had not been). A gene with curated evidence keeps text-mined links whose word names the disease (PPH1
      for BMPR2 in IPAH).
    - **Products.** Excluded publications are counted and listed nowhere (they had remained in the viewer's
      review table, the SQLite file and the run manifest). Counts are labeled by what they count (the readers'
      answers, or the items kept after curation; the related conditions' GWAS without the disease's own). The
      gene map, the regions figure and the regions page state every condition they leave out and why ("No
      related condition shares a region" had been printed for IPAH, where one trait shares one). A condition
      whose genetics are not fetched shows no counts instead of zeros.
    - **The recompute itself showed one more error.** A relation's quote had to contain its condition's name as
      written, so a quote that used an abbreviation or another word order dropped the relation even when the
      adjudicator kept it: PVOD/PCH lost mitomycin-C therapy ("MMC therapy is a potent inducer of PVOD"), IPAH lost
      noncirrhotic portal fibrosis ("NCPF") and autoimmune thyroid disease ("AITD"), HPAH lost fenfluramine exposure
      ("exposure to fenfluramine"). A quote now names a condition also by its abbreviation in parentheses or by every
      word of its name but the words of use or exposure.
    - Not changed: the genetics prompt (its byte identity keeps the cached answers of every disease); a profile's
      `scope` states which patients count for a subtype. A rebuild does not re-read the relations (METHODS §11).
  - **The repository holds the pipeline only** (user, 2026-10-19). Snapshots, the model cache and the slides stay on the
    machine that made them and are shared as files; they were also removed from the repository's history, with the
    run reports and the paths of the machine that ran them.
