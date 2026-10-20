# Methods

How DiseaseLens derives every item of a disease's knowledge graph. The steps are the same for every disease; a
disease enters only through its profile (section 1). CTEPH (`EFO_0022954`) and HPAH (`MONDO_0017148`) are used as
examples.

Every threshold, list and grouping named here is also an entry of the choices registry, with its value in the run
and its basis (section 10). The registry is `config/choices.yaml`, merged with the disease's own
`diseases/<key>/choices.yaml`. The knowledge base and the workbook print the merged registry with each snapshot.
Every rule is either a parameter of the run (`nextflow.config`) or a file of `config/` or `diseases/<key>/`; no
threshold, list or pattern is hidden in the code.

The stages of a disease each write one numbered folder of its snapshot:
- sources (00);
- related conditions (01);
- the genetic studies, as extracted (02);
- curation (03);
- long tables (04);
- the graph (05);
- the products: summaries (06), changes (07), knowledge base (08), figures (09), the workbook, the viewer and the
  index.

## 1. The disease profile (`scripts/new_disease.py`, `scripts/profile.py`)

A disease is identified by its Open Targets disease ID (EFO, MONDO, HP, Orphanet), never by a name, which can name
several diseases. Its folder `diseases/<key>/` holds `profile.yaml`, in which only `id` is required.

**From the ID.** `new_disease.py <ID>` writes the profile from the Open Targets record of the disease, each value
with the rule that chose it:
- **name**: the Open Targets name.
- **short name**: the shortest acronym among the exact synonyms that spells the name: each word of the name gives
  its first letter, and may give further letters of itself, in order (CTEPH: Chronic ThromboEmbolic Pulmonary
  Hypertension; HPAH). Without one, the name itself when it has at most 11 characters (Asthma); else the initials
  of the name, stop words left out: a short name made up here, marked in the profile for checking, a label of the
  products that is never searched (`--short` gives one). At most 11 characters, because workbook sheet names are
  limited to 31.
- **PubMed terms**: the name, the exact synonyms, the short name (not one made up), and the names and exact synonyms
  of the entries whose evidence the disease takes (`evidence_also`), in title or abstract. The name is always kept,
  whatever its form; a name in database form "A, B" is also searched as "B A". Left out are:
  - a term Open Targets maps to a gene (symbol or alias: PAH is phenylalanine hydroxylase; PPH1 is an alias of
    BMPR2);
  - a term Open Targets maps to another disease, other than the disease, its descendants and its `evidence_also`
    entries ("pulmonary hypertension" is an exact synonym of pulmonary arterial hypertension, and the name of another
    disease); the parts of a grouping name ("A and/or B") are exempt;
  - an acronym of two characters or fewer;
  - database forms of a name (with parentheses, brackets, " - ", NOS, a comma that does not list, or in capitals
    with spaces).

  An acronym that is kept is searched only together with a word of the name (HPAH AND heritable, pulmonary,
  arterial or hypertension). Alone, an acronym names other things too: HPAH is also human phenylalanine hydroxylase,
  and halogenated polycyclic aromatic hydrocarbons. This rule removed 57 of 331 records from the HPAH search.

  The terms left out are listed in the profile with the reason.

**What a profile can add.** A profile can also set a scope (which studies count as studies of the disease), the
ontology rule, the categories of related conditions, short labels and abbreviations, and expectations for the
checks. A disease folder can also hold curated files, each replacing a default:

| File | Replaces | CTEPH |
|---|---|---|
| `related_conditions.yaml` | clinical relations read from PubMed (section 4) | 27 conditions, each statement with a PMID and a verbatim quote |
| `gwas_regions.yaml` | GWAS regions from Open Targets (section 3) | the nine regions of Liley et al. (2024), Table 1, with G and G* |
| `gene_themes.yaml` | Reactome top-level pathways as gene functions | eight themes, each with its basis and references |
| `curation.yaml` | (adds to the shared rules) findings not about the disease; text-mining corrections | 18 findings, 30 corrections |
| `choices.yaml` | (adds to, or replaces items of, the shared registry) | the bases of its curated files |
| `text/<product>.yaml` | (adds) disease prose in the products | the Liley table, the PE co-analysis |

**Resolution.** `profile.py` merges the profile and the disease's files with the shared rules of `config/`
(`search.yaml`, `curation.yaml`, `choices.yaml`) into one `profile.json`. That file is saved in the snapshot
(`_run_info/profile.json`), and every later step reads it. A change to any file of the disease folder, its subfolders
included, re-runs the disease. A rebuild records the queries that were run (saved with the snapshot's searches), never
new ones.

## 2. Sources

- **Open Targets Platform.** The GraphQL API, `https://api.platform.opentargets.org/api/v4/graphql`; the release is
  recorded in every snapshot.
- **PubMed.** NCBI E-utilities (esearch, efetch), without a key; requests are spaced 0.4 s apart. A search finds the
  records that PubMed held on the fetch date (Entrez date, the date a record entered PubMed), so that the same date
  finds the same records later. A search that finds more than 9,999 records, the most PubMed returns at once, is
  split by date until each part fits. An error PubMed reports for a query stops the run; a query that finds nothing
  is logged as a warning.
- **Europe PMC.** REST `fullTextXML` for the Methods sections of open-access full texts; REST search for the
  first-author affiliation of records that PubMed holds without one. Only a collective first author (a consortium)
  takes the first affiliation of the author list; no other first author is given another author's institution.
- **The PubMed record of every publication the graph uses** (`scripts/publications.py`): publication types, journal
  issue and pagination, whether PubMed holds an abstract, the comments and corrections linked to the record
  (RetractionIn, ErratumIn, ...), and the first author's affiliation. Read once per snapshot and saved
  (`00.pubmed/publications.json.gz`); the curation (section 7) rests on them.
- **Ensembl REST (GRCh38).** Only to place published lead variants that the Open Targets search does not index.

For the disease itself, Open Targets is read for:
- every associated target, with its datasource scores;
- its text-mined (Europe PMC) and GWAS evidence;
- the targets' Reactome top-level pathways;
- its GWAS studies and credible sets, with every colocalization computed for them;
- its curated gene–disease evidence (section 3);
- the published regions of its profile, if any.

Every list is checked against the count the API reports, and its IDs must be unique. Paged lists are not stable
where the sort key ties, so a list longer than its first page is re-read as one page of its full length; a list
longer than the largest page the API serves (3,000) is read in such pages, again, until the union of its rows (by ID)
holds them all.

Every query runs from the machine that launches the pipeline. Requests are retried with exponential backoff on server
errors, HTTP 429 and network failures.

## 3. The genetics of the disease: regions, curated genes and the evidence ladder

### GWAS regions (`ot_common.disease_regions`)

- **Published regions.** When the profile gives a published table (`gwas_regions.yaml`), its regions are the
  disease's regions, each named by the gene the publication names. For each lead variant the pipeline looks up:
  1. its Open Targets variant ID: the hits of the Open Targets search whose rsIDs include it (the search also
     returns longer rsIDs that begin with it), else the Ensembl GRCh38 alleles checked against Open Targets; at a
     site with several alleles, the first that holds credible sets;
  2. a credible set of the disease at it;
  3. the credible sets of other GWAS that contain it, with its PIP in each.

  The table's `evidence` labels map to codes. CTEPH: `meta` → G, the CTEPH meta-analysis P < 5 × 10⁻⁸; `pe` → G\*,
  genome-wide significant only after leveraging pulmonary embolism (Liley et al. 2024).
- **Regions from Open Targets.** The genome-wide significant credible sets of the disease's own GWAS (lead P <
  `GenomeWideP`, 5 × 10⁻⁸) that no published region holds, whose lead variant is credible. Rare variants are
  evidence when they are credible (user, 2026-10-20); the frequency of a lead is that of the study's largest ancestry
  group, as Open Targets gives it from gnomAD (a lead absent from gnomAD has frequency 0). A lead is credible when:
  - its study carries enough copies of its minor allele: at least `RegionMinMac` (25) expected in the smaller group
    of the study (cases, or all samples for a study without cases and controls), 2 × n × MAF. A case-control test
    needs about 25 copies to reach genome-wide significance, and fewer give miscalibrated tests (Neale lab, 2017);
  - and, when it is rarer than `RegionMinMaf` (0.01), its P is below `RareVariantP` (5 × 10⁻⁹), the threshold of a
    study that tests rare variants too (Pulit et al., 2017).

  A lead that is not credible still makes a region when another GWAS study of the disease has a genome-wide
  significant credible set within `RegionWindowKb`. Otherwise it is no region, and it is listed with its reason
  (`04.tables/disease_signals_not_regions.tsv`, the regions page of the viewer, sheet 5 of the workbook). IPAH: the
  fourteen genome-wide significant credible sets of one PheCode GWAS of a biobank (2,353 cases) had 0 to 23 expected
  copies of their minor allele in the cases, and no other GWAS replicated them. Credible sets whose lead
  variants lie within `RegionWindowKb` (500 kb) are one region, represented by its credible set with the smallest P;
  a credible set within that distance of a published lead belongs to the published region. The region is named
  by its top locus-to-gene gene when the score reaches `L2gMin` (0.5), and that gene gets the code G; below it, the
  region is named by the position of its lead variant and gives no gene a code. A name an earlier region has is
  followed by the position. P is written from Open Targets' mantissa and exponent (a number below 10⁻³⁰⁸ has no
  floating-point value).
- **A disease without either** (HPAH: no GWAS in Open Targets 26.09) has no region. Its products say so.

### Curated gene evidence of the disease (code C)

The curated / Mendelian sources of Open Targets are ClinGen, Gene2Phenotype, Genomics England PanelApp, Orphanet,
ClinVar (eva), UniProt and gene burden.
- A target is kept when the datasource score of a source reaches `MendelianMin` (0.5) for the disease or its ontology
  descendants. The score is Open Targets' harmonic sum of the source's evidence scores scaled to 0–1, so a single
  evidence item passes only with a score of about 0.82.
- Every evidence row of those targets is then read, with its classification, clinical significance, specific
  disease and PMIDs. For ClinVar, the record's confidence is its review status, kept in its own column.
- **A ClinVar link counts only with a reviewed pathogenic record** (`ot_common.curated_links`): a pathogenic or likely
  pathogenic record whose review status has two stars or more (criteria provided by multiple submitters without
  conflicts, an expert panel, or a practice guideline). A single submitter, or a record without assertion criteria,
  does not make a gene a curated gene of a disease; the gene keeps any other source that passes. The same rule holds
  for the related conditions. (HPAH: CAPNS1 and NOP58 had only such ClinVar records.)
- **A gene without an approved symbol is left out** (its Ensembl ID stands in place of a symbol; HPAH:
  ENSG00000289490). The links left out are counted in `04.tables/tidy.json`.
- **A gene whose symbol is a word that names the disease or an ontology neighbor is left out** (section 7). Such a
  gene is linked by its abbreviation, not by evidence (HPAH: PAH, phenylalanine hydroxylase).

### The evidence ladder (`ot_common.gene_ladder`)

Every gene with genetic evidence for the disease takes its strongest level:
1. **G**: a gene of a genome-wide significant GWAS region of the disease; or **C**: curated gene–disease evidence for
   the disease (code C, when no G applies);
2. **G\***: a gene of a region significant only in a secondary analysis of a published GWAS (CTEPH: after leveraging
   pulmonary embolism);
3. associated in two or more genetic studies of the disease;
4. associated in one study;
5. a variant identified, without an association test (case reports, sequencing series);
6. tested, not associated;
7. other result only.

Level 3 ranks above level 4 because replication is what makes an association credible rather than a chance finding
(Kraft 2009). Studies are counted, not independent samples: some share cases.

Levels 3–6 count only findings about the disease itself (section 7). An association with an outcome or phenotype
within the disease (survival, clinical worsening, exercise capacity) modifies the disease, but is not evidence that
the gene raises its risk. Such findings are kept and shown; a gene with nothing else stays at level 7.

The ladder a snapshot uses, with its labels, is written to `04.tables/ladder.json`. Level 1 names C only when a gene
of the snapshot has curated evidence; CTEPH has none in Open Targets, so its level 1 is G only.

### Function (gene themes)

The second axis of the ladder is the function of a gene. It comes from:
- **Curated themes** (`gene_themes.yaml`), each with its basis and references. CTEPH has eight: fibrinogen; ABO,
  VWF, FVIII, FXI; thrombophilia and hemostasis; HLA and immunity; clonal hematopoiesis, MPN and platelets; PAH genes
  (ClinGen, Welch 2023); endothelium and vascular tone; warfarin pharmacogenetics. A gene belongs to the first theme
  that lists it.
- **Otherwise, automatically:** each gene joins its Reactome top-level pathway that holds the most of the disease's
  genes, ties broken alphabetically. The top level "Disease", a mirror of disease variants of normal pathways, is not
  used. The automatic themes are also written as `04.tables/themes_auto.yaml`, a draft to curate.

Genes without a theme are listed as "Other genes", and findings without a gene as "No gene named".

## 4. Which conditions are related, and why (`scripts/relations.py`, `scripts/resolve.py`)

A condition enters the graph through at least one kind of evidence, and every relation is kept as an evidence row
with its source.

1. **Clinical.** A guideline, review, meta-analysis, or cohort or case-control study states the relation, and its
   sentence is quoted verbatim.
   - **With a curated list** (CTEPH), the list is the evidence. Its categories follow the 2022 ESC/ERS guidelines and
     the ERS statement on CTEPH: precursor / spectrum, risk factor, comorbidity, intermediate trait.
   - **Without one** (HPAH), the relations are read from PubMed:
     - The search is the disease's terms, AND guidelines, consensus statements, meta-analyses, systematic reviews,
       cohort and case-control studies, AND relation words (risk, associat\*, comorbid\*, complicat\*,
       predispos\*, preval\*, inciden\*), up to the fetch date (`config/search.yaml`).
     - Two language-model readers read each abstract independently (section 6). For each relation they give the
       condition as written, its usual name, its category, its direction, the design, the population, a one-line
       statement and a verbatim quote.
     - A relation is stated for the disease itself: for a group of patients only when the disease's patients are
       most of it, or when the record gives their result apart (not a broader group in which they are a minority).
     - What is not a relation: a consequence of the disease or of its treatment, an outcome, a prognostic marker
       or a measurement among patients, a sign, symptom, imaging, histological or laboratory finding of the disease
       itself, a molecular or cellular process (such as inflammation), a treatment, a gene, a demographic attribute
       (sex, age, race or ethnicity), the disease or one of its subtypes, and a comparison between groups of
       patients (such as carriers and non-carriers of a variant).
       Only evidence from human patients or populations counts, not an animal model or a cell or tissue experiment.
     - A relation is kept when its quote is in the abstract and names its condition, and when both readers report
       it (condition, direction and category) or the adjudicator confirms it. The quote names the condition by its
       full name, its name without parentheses, an abbreviation the name gives in parentheses ("MMC therapy" for
       mitomycin-C (MMC) therapy), or every word of the name but those that say how it was met ("exposure to
       fenfluramine" for fenfluramine exposure); a different term is not the condition.
     - Each condition is mapped to the Open Targets disease whose name or exact synonym matches (case, spaces,
       hyphens and spelling aside), else it has no genetics; names mapped to one entry are one condition. Names
       without an entry are merged when they share a form: without the words that say how a condition was met
       (use, exposure, intake, therapy, treatment, administration), each alias in parentheses a form of its own
       ("Anorexigen use", "Appetite suppressant (anorexigen) exposure" are one condition), except an alias that is
       an abbreviation (one word of at most five letters or digits, or in capitals: PE is pulmonary embolism and
       pleural effusion); the merge is logged. The disease itself, and the entries whose evidence it takes, are
       left out, and a condition supported by one publication is marked as such.

   In both cases, every quote is re-checked against the abstract fetched for the snapshot, and the build stops if
   it is not found.
   - **A finding of no association is not an association.** The direction of each relation (increased risk,
     decreased risk, associated, no association found) is judged by the two readers apart and the adjudicator
     where they differ, or given by the curated list. A condition enters by clinical evidence only through a
     relation that found an association; one whose clinical evidence all found none is not a related condition
     (`01.conditions/not_related.tsv`, named on the How it is built page). A study that found no association is
     kept, labeled, beside the relations of a condition that is related otherwise, but it is not counted as a
     clinical item and does not draw the map's clinical line or its East Asian diamond.
   - **East Asian evidence of a relation** is the population of the evidence that relates the condition to the
     disease: clinical evidence from a study done in East Asia (the country of the first author's institution, as
     for the genetic studies) or of an East Asian population (the final reading's population and the abstract both
     name an East Asian country, people or major city as a whole word, and neither names a country outside East
     Asia: a multinational study is not East Asian), or a shared signal whose two
     GWAS, the disease's and the condition's, have mostly East Asian samples
     (`EasMin`). A condition's own East Asian GWAS do not make its relation East Asian; the ontology and a
     lead-variant overlap never do.
2. **Ontology.** From the Open Targets disease ontology (EFO / MONDO): the parents of the disease, their other
   children, the disease's own children, and the children of any disease the profile names (CTEPH: the children of
   pulmonary arterial hypertension).
   - A parent with more than `MaxOntologyChildren` (30) children is too broad a class to give neighbors: it and its
     children are left out, and `01.conditions/ontology_log.tsv` names it. HPAH: "hereditary disease", with 1,070
     children.
3. **Colocalization,** computed by Open Targets and selected here; the pipeline runs no colocalization itself.
   - Open Targets fine-maps GWAS credible sets and reports, for each pair of overlapping credible sets, the coloc
     posterior H4 and the eCAVIAR CLPP.
   - A condition is kept when one of its credible sets reaches CLPP ≥ `ClppMin` (0.01) or H4 ≥ `H4Min` (0.8) with a
     credible set of a region of the disease. The disease's own GWAS (with its descendants') and the entries whose
     evidence it takes are not related conditions.
4. **Lead-variant overlap** (weaker), computed by this pipeline, and only at a published region whose lead variant
   Open Targets holds without any credible set of the disease's GWAS. A credible set of the condition must contain that lead
   variant with PIP ≥ `ClppMin`. If the lead were the causal variant, the CLPP of the pair would equal its PIP;
   without a credible set of the disease that side is not fine-mapped, hence weaker.

**GWAS signal counts.** A GWAS signal is one credible set of one GWAS study of the condition that colocalizes with,
or contains the lead variant of, a region of the disease. Counts are numbers of distinct GWAS studies, per condition
and region, per condition, and per region. Studies of one trait from different consortia or releases each count and
can share participants, so a count measures consistency across studies, not independent replication. A relation of
the shared-signal kinds is graded by the number of regions it shares: two or more is stronger than one.

**After the curation.** The sources, category and number of clinical items of each condition are taken again from the
relations the curation keeps (section 7); a condition left without evidence that relates it leaves the graph, and a
GWAS study the curation leaves out defines no region and no shared signal.

**Depth of the genetics of a condition.**
- `full`: its GWAS and curated gene evidence are fetched;
- `coloc_only`: only its shared signals are kept. This applies to traits related only by a shared signal, and to a
  condition with more than `MaxGwasStudies` (100) GWAS studies;
- `none`: no Open Targets record.

A condition whose genetics are not fetched has no counts of its own (blank, not zero), and `genetics_note` says why.

## 5. Genetics of the related conditions (`scripts/fetch.py`)

- **GWAS.** Every GWAS study of the condition, with its sample size, cases, controls, cohorts, ancestry and summary
  statistics. For each study, every credible set: its lead variant, P (mantissa and exponent) and β, the fine-mapping
  method, and the top three locus-to-gene genes. A gene is linked to a condition when it is the top L2G gene with a
  score ≥ `L2gMin` (Mountjoy 2021). A study is East Asian when its EAS share is ≥ `EasMin` (0.5).
- **Curated evidence,** as for the disease itself (section 3).

The conditions are fetched three at a time.

## 6. The genetic studies of the disease, read by two independent models (`scripts/extract_genetics.py`, `scripts/llm.py`)

1. **Search.** PubMed: the disease's terms in title or abstract AND the genetic-study block of `config/search.yaml`,
   up to the fetch date.
2. **Pre-screen** (deterministic, before any model reads a record). A record is excluded when:
   - it is typed Review, Systematic Review, (Practice) Guideline, Editorial, Comment, News or Consensus Development
     Conference; or
   - its title and abstract name no genetic-study term (the patterns of `config/search.yaml`).

   CTEPH, 2026-10-16: 98 of 471 records were kept, including all 15 original genetic studies the curated seed cites
   (the benchmark, re-checked at every build).
3. **Reading the abstracts (pass A).**
   - Each record is read by two language models, five records per call. The models never see each other's answer.
   - Each decides whether the record is an original genetic study (or case report) of the disease, within the
     profile's scope when one is given.
   - For an included record, each extracts:
     - design;
     - population (country, ancestry, Japanese / East Asian);
     - cases, controls and the control type;
     - genes (as written, plus the HGNC symbol) and variants;
     - findings (result, effect, P, and the source sentence);
     - replication;
     - the methods the abstract states.
4. **Methods (pass B).** For every included study, both readers extract the methods from the Methods section of its
   open-access full text (Europe PMC), or from the abstract when there is none. The fields follow STREGA where it
   applies: cohort, diagnosis, controls, DNA source, platform, sequencing, HLA typing, variant calling and quality
   control, imputation, statistical test, covariates, multiple testing, replication, power, software. Each field
   comes with its source passage.
5. **Grounding.** Every gene, variant, finding quote and sample size must occur in the source text (a gene or variant
   name of at most four characters as a whole word in its own case: "F5" is not in "one of 5 patients"), and every
   methods value needs a quote that occurs in its text; a field stated twice (two cohorts, two platforms) keeps both.
   A Japanese or East Asian population must be named as a whole word: a country, people or major city (Japan,
   Japanese, China, Chinese, Korea, Korean, Taiwan, Hong Kong, Mongolia, East Asian; Beijing or Peking, Shanghai,
   Tokyo, Seoul, ...); "Asian" alone, which also names South Asia, is not enough. Independently of the readers, a
   study done in East Asia (the first author's institution, section 7) is East Asian. Items
   that fail are dropped and logged (`extract_log.tsv`), for each reader. An answer of the wrong shape counts as no
   answer for its item, never as a value. Gene symbols are kept only when Open Targets knows them by exact name.
6. **Comparison** (deterministic, per record). The decisive items are compared:
   - inclusion;
   - sample sizes;
   - the Japanese and East Asian flags;
   - each finding, by its gene (symbol, else as written; a finding without a gene, by its variant) and result,
     "enriched" counting as "associated";
   - each methods field: two quotes agree when one contains the other, or when they share at least 60% of the words
     of the shorter (a field stated twice agrees when any quote of one reading agrees with one of the other).

   Descriptive fields (reason, design wording, country, control type) come from the first reader.
7. **Adjudication.**
   - Where the readers differ, or one gave no answer, a third call reads the record with both readings and gives
     the final reading. Its answer is grounded in the same way.
   - Only agreed or adjudicated items enter the graph, each with `agreement` (both / adjudicated). With one reader,
     its reading is kept (agreement: one reader); two readers need an adjudicator (the configuration is refused
     otherwise). The agreement on inclusion counts a record as adjudicated only when the readers differed on its
     inclusion. A model whose answer is still no JSON after three calls gives no answer for that question, which is
     logged; a backend that keeps failing (network, command-line tool) stops the run, with every answer so far
     cached.
   - Every disputed item is listed, with both readings, the final value and the adjudicator's reason, in
     `02.literature/review.tsv` and `03.curated/model_review.tsv`.
8. **Result check.** A finding whose result disagrees with the wording of its own sentence is flagged for review:
   "not associated" without a negative or non-significant phrase, or "associated" with only a negative one. The
   patterns are in `config/curation.yaml`.
9. **Isolation.**
   - **Command-line backend:** `claude -p --safe-mode --tools "" --strict-mcp-config --no-session-persistence`, with
     auto-memory disabled, run from a scratch directory. No memory, file or instruction of the user's reaches the
     model; `verify.sh` section 6 tests this.
   - **API backends:** only the prompt is sent.
10. **Caching and reproducibility.** Answers are cached by hash of (model, effort, system prompt, prompt) in
    `cache/llm/`, with the model that answered (an OpenAI-compatible model also by its endpoint). The answer of a
    batch is also cached record by record, so a record is read once whatever batch it falls in: a later snapshot
    whose search adds or removes records asks only about the new or changed ones, and a rebuild with the same cache
    gives the same result. The prompts of CTEPH are those of its first build, character for character, so its cached
    answers are reused.

## 7. Curation (`scripts/publications.py`, `scripts/curate.py`)

Between what was fetched and extracted and what the graph holds, a deterministic step applies the curation rules of
the profile (`config/curation.yaml` merged with `diseases/<key>/curation.yaml`). Every rule it applies is logged,
row by row (`03.curated/curation_log.tsv`). No rule calls a language model.

1. **Publications left out.** Every publication the graph uses is looked up in PubMed. Left out of every table:
   - retracted articles and retraction notices (publication type Retracted Publication or Retraction of
     Publication, or a RetractionIn link);
   - conference material (publication type Congress or Conference Proceedings, or a journal-supplement record
     without an abstract).

   Left-out publications are counted and not listed in any product (the viewer, the workbook, the knowledge base,
   the graph and its SQLite file, the run manifest, the changes); they remain in `03.curated/publications.tsv` with
   their reason, the record of the curation. A GWAS study they report defines no region and no shared signal, and
   their disputed items leave the review table. A curated list of relations may not cite one (the build stops); a relation read from PubMed that does is
   left out. Preprints, errata and expressions of concern are kept, flagged.
2. **Place of a study.** The country of the first author's institution, from the PubMed affiliation (Europe PMC when
   PubMed has none), mapped by country names, aliases, US states and cities. The country the abstract states is kept
   beside it. The place sets the Japanese and East Asian flags, as does a population stated in the abstract.
3. **Design class.** One class per study, the first that applies: case report; genome-wide; clonal hematopoiesis;
   exome or gene panel; HLA typing; candidate gene; other.
4. **Variant names.** One name per variant and, for the classic variants, its dbSNP rsID, from patterns per gene. The
   name as written is kept beside it.
5. **What a finding is about.** The disease itself (the default); a phenotype or outcome within it; another trait
   (Mendelian randomization with the disease as the exposure); a functional experiment; a background statement.
   - The findings of each kind are listed in the disease's curation file, each read against its source sentence; an
     entry that matches no finding stops the build.
   - A later finding with words of an outcome (survival, worsening, predictor, ...) is logged for review.
6. **Text-mined links.** Europe PMC text mining maps words to genes by dictionary.
   - A matched word that is the name or an exact synonym of the disease, of a parent, of a child or of an entry
     whose evidence it takes (Open Targets) is a false link (HPAH: "PAH", mapped to phenylalanine hydroxylase),
     except for a gene with curated evidence for the disease, for which the word is the old name of its locus (IPAH:
     "PPH1", the BMPR2 locus of familial primary pulmonary hypertension).
   - Other false words are read against their sentences and listed for the disease (CTEPH: sGC, soluble guanylate
     cyclase, mapped to SGCB; PTPA and BPA, pulmonary angioplasty; TF, tissue factor, mapped to transferrin).
   - A link whose every word is wrong is dropped, with its sentences; a wrong word is removed from a link that keeps
     another. A listed correction that matches no link stops the build.
7. **Curated genes named like a disease.** A gene of the curated evidence of the disease whose symbol is such a word
   is left out (section 3).

## 8. Mechanisms and summaries (`scripts/summarise.py`)

- **Mechanism themes.** The text-mined targets of the disease are grouped by Reactome top-level pathway. For each
  theme with at least `MinThemeTargets` (3) targets, a summary of at most 150 words is written. Its input is up to two
  sentences from the best-scored publications of each of the theme's `ThemeTargets` (12) best literature targets.
- **Methods and findings.** Summaries of at most 350 words each, written from the extracted tables and their tallies.
- **Checks.**
  - Every [PMID] a summary cites must be one of its inputs, and every gene symbol of its material (three characters
    or more, not a name of the disease: PAH is also a gene) it names must occur in its source text.
  - The second reader reads each summary with its material and lists the claims it does not support, and the
    adjudicator gives a ruling, true or false, on each.
  - A summary that fails, that the second reader did not check, or whose disputed claims have no ruling each, is
    kept but withheld from the knowledge base and the workbook, with the reason.
- **Offline** (`--LlmOffline true`). No model is called. A summary whose inputs are unchanged is read from the cache;
  one whose inputs changed is reported as not generated, never invented.

## 9. Tables, graph, changes, figures, knowledge base, workbook, viewer

- `tidy_graph.py` writes the long tables (docs/OUTPUTS.md), the regions, the ladder and the themes. The
  gene–condition genetic evidence is:
  - GWAS: for the disease, the genes of its regions (G, or G\*); for a related condition, its top L2G gene ≥
    `L2gMin`;
  - curated (C);
  - literature (L): a positive finding about the disease itself.
- `build_graph.py` writes the typed graph (`05.graph`, schema in its README). Node types: condition, gene, region,
  publication, gwas_study, theme, pathway. Edge types: related_to, supports, has_genetic_evidence, reports, studies,
  region_of, names_gene, shares_region, in_theme, in_pathway, mentions, gwas_of. Each edge type is checked against
  the node types it joins. The graph is written as node and edge tables, `graph.json` and a SQLite database that also
  holds every curated and long table.
- `build_viewer.py` writes one HTML file, offline and without libraries: an overview, lists of genes, conditions,
  regions and studies, the curation, and a page per node with its neighbors and every edge with its source link.
- `diff_snapshots.py` compares stable keys with the newest earlier snapshot of the disease. A table the earlier
  snapshot lacks is listed as not compared; rows of publications the curation now leaves out are counted, not
  listed.
- **Counts.** The GWAS of the related conditions are counted apart from the disease's own, a study of two conditions
  once. A product that leaves conditions out of a map or a table says how many, and why (the gene map: traits
  related only by a shared signal, conditions with more GWAS than are fetched, conditions without an Open Targets
  entry, conditions with no gene of the map).
- **Populations.** A literature study's population is the country its abstract states, else its ancestry, else "not
  stated"; its place is the country of its first author's institution (section 7).
- `draw_figures.py` writes the figures, in the shared figure style (`scripts/plot_style.py`), with one catalog
  README. A figure without data (the regions of a disease without regions) is left out and the catalog says so.
- `build_knowledge.py` writes the Markdown/JSON knowledge base, in the order of the workbook.
- `write_table.py` writes the workbook, organized by question in four parts: A, the genetics of the disease; B,
  related conditions; C, context; D, provenance. Each question sheet:
  - opens with its answer in brief and states the basis of the values it uses;
  - then gives a nested table, group › record › source, with Excel row outlines;
  - keeps source passages in cell notes and links the same gene or condition across sheets.

  The appendix holds the flat tables, with filters. The workbook is re-opened and checked before publication.
- `manifest.py` writes the run manifest, the disease's `latest.json`, the index of the snapshot and, for the whole
  results folder, `results/index.json` and `results/README.md`.
- `scripts/kg.py` is not a step: it is a query tool over the SQLite file of a disease (gene, condition, region,
  publication, studies, search, neighbors, sql), from the shell or from Python.

## 10. The choices registry

Each item of the registry is one choice: what it decides, its value, its basis type and its basis, and where its
result is used. Basis types:
- **source:** a publication or the documentation of a data source, with references;
- **data check:** a check on the snapshot's data, recomputed at every build (the regions of the disease in Open
  Targets; the region most often shared alone; the curated genes found by the extraction; the readers' agreement; a
  benchmark of known studies);
- **design:** a decision of this project, with its reason;
- **user:** decided with the user, with the date.

Values that are parameters are read from the run itself (`run_params.json`, also in the run manifest), never typed
into the registry; the disease's ID and query come from its profile. `verify.sh` checks:
- that every item has a basis and resolvable references;
- that every parameter named is recorded by the run;
- after a run, that every data check holds.

## 11. Rebuild from the saved sources (`--Rebuild true`)

A snapshot can be re-derived from its own saved sources (00.raw, 00.pubmed, 01.conditions) after a change to a rule,
so that it stays the snapshot of its fetch date. Nothing is fetched again, with two exceptions:
- the publication records, when the snapshot has not saved them yet;
- Open Targets gene symbols, looked up only after the API is shown to serve the snapshot's own release.

The readers answer from the cache, and only questions not yet cached go to a model (none with `--LlmOffline true`).
The run manifest records the rebuild.

What a rebuild does not re-derive: the related conditions (`01.conditions`) and, for a disease without a curated
list, the relations read from PubMed (`01.relations`). Their resolution queries Open Targets, and a new condition
would need genetics that the snapshot did not fetch. A change to the rules of the relations therefore takes a new
fetch of the disease.
