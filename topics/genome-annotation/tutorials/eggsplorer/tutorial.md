---
layout: tutorial_hands_on

title: "Exploring the functional profile of a bacterial genome with Eggsplorer"
zenodo_link: ""
questions:
  - How can the functional annotation of a whole bacterial proteome be summarised and explored?
  - Which metabolic, respiratory, transport and secretion capabilities can be read from eggNOG-mapper annotations?
  - How far can we trust a function that is predicted only from orthology?
objectives:
  - Download the annotated proteome of a reference bacterium from NCBI directly into Galaxy
  - Annotate all proteins with eggNOG-mapper
  - Generate an interactive Eggsplorer report from the eggNOG-mapper annotations
  - Interpret the COG and GO profiles, the cell overview, the metabolic function scores and the element cycle diagrams
  - Recognise the limits of annotation transfer by orthology
time_estimation: 1H
level: Introductory
key_points:
  - NCBI Datasets brings a reference proteome into Galaxy by its assembly accession, without manual downloads.
  - eggNOG-mapper transfers COG, GO, EC, KEGG, CAZy, TC and Pfam annotations from orthologous groups to every protein.
  - Eggsplorer turns the eggNOG-mapper table into one interactive report that can be explored offline and shared as a single file.
  - The cell overview and the element cycle diagrams show what the genome can potentially do, not what the cell does under a given condition.
  - Completeness scores depend on KO assignments; steps that rest on the same proteins, or on non-catalytic subunits, need careful interpretation.
requirements:
  - type: "internal"
    topic_name: introduction
  - type: "internal"
    topic_name: genome-annotation
    tutorials:
      - functional
contributions:
  authorship:
    - vebaev
tags:
  - microgalaxy
  - prokaryote
  - functional-annotation
subtopic: prokaryote
edam_ontology:
  - topic_0085 # Functional genomics
  - topic_3301 # Microbiology
  - topic_0092 # Data visualisation
---

When a bacterial genome is sequenced and annotated, we obtain a list of several thousand genes and proteins. **Functional annotation** tools such as **eggNOG-mapper** {% cite Cantalapiedra2021 %} then assign each protein to orthologous groups of the eggNOG database {% cite HuertaCepas2019 %} and transfer their functional descriptions: COG categories {% cite Tatusov2000 %}, Gene Ontology terms, EC numbers, KEGG orthologs (KO) {% cite Kanehisa2000 %}, transporter (TC) families, carbohydrate-active enzyme (CAZy) families and Pfam domains.

The result is a large table with one row per protein and more than twenty columns. It is complete, but it is hard to answer simple biological questions from it: *Does this bacterium respire nitrate? Which secretion systems does it have? Is it motile? Can it fix CO~2~?*

**Eggsplorer** {% cite Baev2026 %} is a Galaxy tool that turns the eggNOG-mapper annotation table into a single, self-contained interactive HTML report. Among others, the report contains:

- an **overview** of how many proteins were annotated and at which taxonomic level the annotations were transferred,
- the **COG** profile as a radar chart and the most frequent informative **GO** terms,
- a **cell overview** that places the respiratory chain, ATP synthase, uptake and secretion systems, transporter families, motility and groups of metabolic functions on a fixed drawing of a bacterial cell,
- the **completeness of 178 metabolic functions**, scored with the definitions of KEGG-Decoder {% cite Graham2018 %},
- **nitrogen, sulfur and carbon cycle diagrams**,
- interactive **Sankey diagrams** for transporters and enzymes, **CAZy** families grouped by substrate {% cite Zheng2023 %}, **Pfam** domains and a searchable gene table.

In this tutorial we follow the full path from a genome accession to an interpreted report, using the reference strain ***Pseudomonas aeruginosa* PAO1** {% cite Stover2000 %} (RefSeq assembly `GCF_000006765.1`). *P. aeruginosa* is an opportunistic pathogen and a model Gram-negative bacterium with a versatile, branched respiratory chain {% cite Arai2011 %}, several secretion systems and a large regulatory repertoire. Because its biology is very well described, it is a good genome on which to check what the report shows.

> <comment-title>Where does this tutorial fit?</comment-title>
>
> This tutorial continues the [Functional annotation of protein sequences]({% link topics/genome-annotation/tutorials/functional/tutorial.md %}) tutorial. There, eggNOG-mapper is run on a few proteins and its output table is inspected by eye. Here we run it on a complete bacterial proteome and use Eggsplorer to explore the result.
{: .comment}

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get the proteome of *P. aeruginosa* PAO1

NCBI distributes genome assemblies together with their annotation: genomic sequence, GFF3, GenBank file, and the predicted protein sequences. The **NCBI Datasets Genomes** tool {% cite OLeary2024 %} downloads any of these files into Galaxy from an assembly accession.

> <hands-on-title>Download the PAO1 proteome</hands-on-title>
>
> 1. Create a new history for this tutorial and give it a name, for example `Eggsplorer PAO1`
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. {% tool [NCBI Datasets Genomes](toolshed.g2.bx.psu.edu/repos/iuc/ncbi_datasets/datasets_download_genome/18.39.0+galaxy0) %} with the following parameters:
>    - In *"Query"*:
>        - *"Choose how to find genomes to download"*: `By NCBI assembly or BioProject accession`
>            - *"Enter accession or read from file ?"*: `Enter accessions`
>                - *"Enter comma separated list of accessions"*: `GCF_000006765.1`
>    - In *"Output options"*:
>        - *"Include"*: select only `amino acid sequences (protein)` (deselect `genomic sequence (genome)`)
>        - *"Decompress FASTA"*: `Yes`
>
>    > <comment-title>Outputs</comment-title>
>    > The tool creates a tabular **Data Report** with the assembly metadata and a **collection** `NCBI Genome Datasets: protein fasta` with one element, `GCF_000006765.1`, which holds the protein sequences.
>    {: .comment}
>
> 3. Inspect the **Data Report** and the protein FASTA file
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Which organism and assembly name does the Data Report show?
> 2. How many protein sequences does the FASTA file contain?
>
> > <solution-title></solution-title>
> >
> > 1. *Pseudomonas aeruginosa* PAO1, assembly `ASM676v1`.
> > > 2. 5479 proteins, one per `>` header line. Click on the dataset name to expand it: Galaxy shows the number of sequences in the dataset preview.
> >
> {: .solution}
>
{: .question}

> <tip-title>Working with a collection</tip-title>
>
> The protein FASTA arrives as a collection with a single element. You can run the next tools directly on the collection: Galaxy runs the tool once for each element and keeps the results in a collection. The element name, `GCF_000006765.1`, is passed along to every result.
{: .tip}

# Annotate the proteins with eggNOG-mapper

**eggNOG-mapper** first searches every protein against the eggNOG database (with DIAMOND by default) to find its best match, the *seed ortholog*. It then places the protein into orthologous groups at several taxonomic levels and transfers the annotations of the narrowest suitable group to the protein.

> <hands-on-title>Functional annotation with eggNOG-mapper</hands-on-title>
>
> 1. {% tool [eggNOG Mapper](toolshed.g2.bx.psu.edu/repos/galaxyp/eggnog_mapper/eggnog_mapper/2.1.13+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Fasta sequences to annotate"*: `NCBI Genome Datasets: protein fasta` (output of **NCBI Datasets Genomes** {% icon tool %})
>    - *"Version of eggNOG Database"*: select the latest version available
>    - In *"Output Options"*:
>        - *"Exclude header lines and stats from output files"*: `Yes`
>
>    > <comment-title>Runtime</comment-title>
>    > Annotating about 5500 proteins takes roughly 30 to 60 minutes on a public Galaxy server, depending on the load. This is a good moment to read the next section and look at the example figures of the report.
>    {: .comment}
>
{: .hands_on}

eggNOG-mapper produces several outputs. We need the **annotations** table (`*.emapper.annotations`): one row per protein, with the seed ortholog and its e-value and score, the orthologous groups, the taxonomic level of annotation transfer (`max_annot_lvl`), the COG category, a description and preferred gene name, and the GO, EC, KEGG, TC, CAZy and Pfam columns.

> <question-title></question-title>
>
> 1. Open the annotations table. What is in the `max_annot_lvl` column for most proteins, and what does it mean?
> 2. Find the protein with the preferred name `nosZ`. What is its function, and which KO is assigned to it?
>
> > <solution-title></solution-title>
> >
> > 1. 5300 of the 5479 proteins carry `1236|Gammaproteobacteria`. This is the broadest taxonomic level from which the annotation was taken. The orthologous groups themselves are much narrower: the last group in the `eggNOG_OGs` column is at the level of the *Pseudomonas aeruginosa* group for 4745 proteins. The closer the annotated relatives of a protein, the more reliable the transferred function.
> > 2. `nosZ` encodes nitrous oxide reductase (EC 1.7.2.4, KO `K00376`), the last enzyme of denitrification, which reduces N~2~O to N~2~.
> >
> {: .solution}
>
{: .question}

The table already contains everything we need, but it is not easy to see the biology of the bacterium in five thousand rows. This is what Eggsplorer is for.

# Generate the Eggsplorer report

> <hands-on-title>Create the interactive report</hands-on-title>
>
> 1. {% tool [Eggsplorer](toolshed.g2.bx.psu.edu/repos/galaxyp/eggsplorer/eggsplorer/1.1.0+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"eggNOG-mapper annotations"*: the `annotations` output of **eggNOG Mapper** {% icon tool %}
>    - *"Sample name"*: `Pseudomonas aeruginosa PAO1` (the title of the report; if left empty, the dataset name is used)
>    - *"Add KofamScan KO assignments as a second source?"*: `No`
>    - *"Number of items in 'top' charts"*: `25`
>    - *"Report language"*: `English`
>    - *"Also output summary tables"*: `Yes`
>
> 2. Open the report with the {% icon galaxy-eye %} icon
>
{: .hands_on}

> <tip-title>The report shows no charts</tip-title>
>
> Galaxy removes JavaScript from HTML datasets unless the administrator has allowed it for the tool. If the report opens without charts, download it with the {% icon galaxy-save %} icon and open the file in your web browser. The report is a single file and works offline: the plotting library and the fonts are embedded and no database is queried.
{: .tip}

The tool also produces a collection of tabular **summary tables** (COG counts, top GO terms, KEGG pathways and modules, EC numbers, TC families, CAZy families and substrates, Pfam domains and metabolic function scores). They hold the same numbers as the report and are useful for further analysis in Galaxy or for comparing several genomes.

# Explore the report

The report is organised from general to specific. We will walk through its main sections.

## Overview and level of annotation transfer

The first section gives the number of proteins, the share that received a functional annotation, the number of KEGG orthologs (KOs) and the coverage of every eggNOG-mapper column. The next section shows the taxonomic levels at which the annotations were transferred.

> <question-title></question-title>
>
> 1. What share of the PAO1 proteins received a functional annotation? How many unique KOs were found?
> 2. Which annotation fields have the lowest coverage? Why would you expect fewer CAZy and TC assignments than COG assignments?
>
> > <solution-title></solution-title>
> >
> > 1. 5163 of 5479 proteins (about 94 %) are annotated, and 3455 proteins carry 2790 unique KOs. This is typical for a well-studied model bacterium; for a genome from a poorly studied lineage the share is lower.
> > 2. CAZy (1 % of the proteins), TC (15.5 %), GO (24 %) and EC (27 %) have the lowest coverage, whereas COG categories (94 %) and Pfam domains (93 %) cover almost all proteins. Nearly every protein can be placed into a COG category, if necessary into category S (function unknown), but only carbohydrate-active enzymes receive a CAZy family, only transporters a TC family and only enzymes an EC number. The low GO coverage reflects that eggNOG-mapper transfers GO terms only from orthologous groups with experimental or curated GO annotation.
> >
> {: .solution}
>
{: .question}

## COG categories and Gene Ontology

The COG radar chart groups the proteins into the 26 functional categories of the COG database, arranged in the four super-categories *Information storage and processing*, *Cellular processes and signalling*, *Metabolism* and *Poorly characterised*. Hover over a point to see the category name and the number of proteins.

> <question-title></question-title>
>
> 1. Which COG category contains the most proteins?
> 2. Which are the two largest well-characterised categories? What do they tell you about the lifestyle of *P. aeruginosa*?
>
> > <solution-title></solution-title>
> >
> > 1. Category S, *Function unknown*, with 1176 proteins. Even in a model bacterium, many proteins belong to orthologous groups without a known function.
> > 2. *Amino acid transport and metabolism* (E, 552 proteins) and *Transcription* (K, 523 proteins). Category K contains mostly transcriptional regulators; the most frequent Pfam domains of the genome (`HTH_1`, `LysR_substrate`, `Response_reg`) tell the same story. *P. aeruginosa* has one of the largest regulatory repertoires among bacteria and a broad catabolic range, which fits a versatile organism that lives in soil and water and colonises plants and animals.
> >
> {: .solution}
>
{: .question}

The GO section lists the most frequent informative terms of each ontology (biological process, molecular function, cellular component). The annotations are first mapped to current GO terms and propagated to their parents; very general and redundant terms are hidden.

> <comment-title>Frequency is not enrichment</comment-title>
>
> The GO section shows how often each term occurs in the genome. It is a profile, not an enrichment test: it does not compare the genome with a background and gives no p-values.
{: .comment}

## Cell overview

The **cell overview** draws the components that were found in the genome where they act in the cell, on a fixed template. The respiratory chain and ATP synthase are at the top, uptake systems on the left, secretion and export on the right, and the largest transporter families and motility at the bottom. The tiles in the centre summarise the groups of metabolic functions. The shade gives the completeness of each component: dark for complete (1.00), light for 0.50 to 0.99. Components below 0.50 are not drawn but listed under the drawing.

![Cell overview of P. aeruginosa PAO1](./images/pao1_cell_overview.png "Cell overview of <i>P. aeruginosa</i> PAO1 drawn by Eggsplorer 1.1.0. The respiratory chain and ATP synthase are on top, uptake systems on the left, secretion and export on the right, and the largest transporter families, the flagellum and chemotaxis at the bottom. The central tiles summarise the groups of metabolic functions.")

> <question-title></question-title>
>
> 1. Is the cell envelope drawn as Gram-negative (diderm) or Gram-positive (monoderm)? Which genes support this?
> 2. Which terminal oxidases are drawn? How does this compare with what is known about *P. aeruginosa*?
> 3. Which protein secretion systems does PAO1 have according to the report? Which of them is best known for its role in virulence?
> 4. Some components are drawn in a lighter shade, for example T1SS (0.67). What does this mean?
>
> > <solution-title></solution-title>
> >
> > 1. Diderm. Eggsplorer draws an outer membrane when it finds BamA together with lipid A synthesis (LpxA, LpxC, LpxD) or LptD. In PAO1 the outer-membrane markers `lpxC`, `lptD` and `bamA` are present.
> > 2. The aa~3~, cbb~3~, bo~3~ and bd-type oxidases, all complete, on top of a complete chain from complex I (NDH-1), NDH-2 and Na^+^-NQR through the bc~1~ complex. This agrees with the five terminal oxidases described for *P. aeruginosa*: aa~3~, two cbb~3~ isoforms, bo~3~ and the cyanide-insensitive bd-type oxidase {% cite Arai2011 %}. The two cbb~3~ isoforms share KOs and therefore appear as one component.
> > 3. Besides the general Sec and Tat pathways: type I, II, III and VI secretion systems. The type III secretion system injects effector toxins (ExoS, ExoT, ExoU, ExoY) directly into host cells and is a major virulence factor of *P. aeruginosa*. The type VI system delivers toxins into competing bacteria.
> > 4. Only part of the genes (KOs) defined for this component were found. This can mean that a subunit is truly missing, that it is too divergent to be assigned the expected KO, or that the system is built from a different set of proteins than the definition assumes. Click the tile in the report, or the function in the next section, to see which genes were found and which are missing.
> >
> {: .solution}
>
{: .question}

## Metabolic function completeness

This section lists every metabolic function with at least one step present. The functions and their scoring rules come from KEGG-Decoder {% cite Graham2018 %}: a score of 1.00 means that all required KOs, or all subunits of a complex, were found. Eggsplorer also reports a **corrected** score that fixes known issues of the original definitions. Click on a function to see the supporting genes, their copy numbers, the KOs found and the missing steps that would raise the score.

> <question-title></question-title>
>
> 1. Look up *dissim nitrate reduction*, *nitrite reduction*, *nitric oxide reduction* and *nitrous-oxide reduction*. What can you conclude?
> 2. What is the score of *RuBisCo*? Is PAO1 an autotroph?
>
> > <solution-title></solution-title>
> >
> > 1. All four steps of denitrification score 1.00 (`narGH`/`napAB`, `nirS`, `norBC`, `nosZ`). PAO1 can use nitrate as terminal electron acceptor and grow anaerobically by denitrification, as is well known for *P. aeruginosa* {% cite Arai2011 %}.
> > 2. 0. Without RuBisCO there is no Calvin–Benson–Bassham cycle: PAO1 is a heterotroph. A few enzymes of the cycle may still be found because they also act in the pentose phosphate pathway and gluconeogenesis.
> >
> {: .solution}
>
{: .question}

## Element cycles

The nitrogen, sulfur and carbon cycle diagrams show the transformations the genome can potentially carry out. The colour and width of each arrow give the completeness of the genes for that step, labels name the genes and the number of proteins, and grey dashed arrows are steps without genes. A step whose catalytic subunit is missing is drawn pale, and steps that rest on the same proteins are drawn as green dashed arrows and listed under the diagrams.

![Nitrogen and sulfur cycles of P. aeruginosa PAO1](./images/pao1_cycles.png "Nitrogen (left) and sulfur (right) cycle diagrams of <i>P. aeruginosa</i> PAO1 from the Eggsplorer report.")

> <question-title></question-title>
>
> 1. Follow the nitrogen cycle from nitrate (NO~3~^-^) to dinitrogen (N~2~). Which genes carry each step?
> 2. Is nitrogen fixation possible? Compare with the root-associated *Pseudomonas stutzeri* A1501 {% cite Yan2008 %}, which fixes N~2~.
> 3. The arrow from nitrite back to nitrate (nitrite oxidation, `nxrAB`) is drawn dashed, although *P. aeruginosa* is not a nitrifier. Why?
> 4. In the sulfur cycle, only the inner arrows from sulfate to sulfide are coloured. What do they represent, and why is the outer, dissimilatory pathway absent?
>
> > <solution-title></solution-title>
> >
> > 1. Nitrate reduction (`narGH`/`napAB`, 4 genes), nitrite reduction to NO (`nirS`), NO reduction (`norBC`) and N~2~O reduction (`nosZ`), all complete: the full denitrification pathway. The arrow from nitrite to ammonium rests on `nirBD`, the assimilatory nitrite reductase that provides nitrogen for biosynthesis.
> > 2. No. The nitrogenase genes `nifHDK` are absent. *P. stutzeri* A1501 has a whole *nif* island and combines nitrogen fixation with denitrification {% cite Yan2008 %}; running the same workflow on its proteome (`GCF_000013785.1`) is a good exercise.
> > 3. The nitrite oxidoreductase of nitrifiers (NxrAB) and the respiratory nitrate reductase (NarGH) are related enzymes and share KEGG orthologs. The nitrite oxidation step is therefore supported by exactly the same proteins as nitrate reduction, and Eggsplorer draws it dashed and lists it under the diagram. A KO tells us which family a protein belongs to, not always which reaction it catalyses in the cell.
> > 4. The inner arrows are **assimilatory** sulfate reduction (`sat` or `cysDN`, `cysH`, `cysJI` or `sir`; KEGG module M00176): the cell reduces sulfate to sulfide to build cysteine and methionine. The outer, **dissimilatory** pathway (`sat`, `aprAB`, `dsrAB`) uses sulfate as electron acceptor for anaerobic respiration and is found in sulfate-reducing bacteria, not in pseudomonads. Two further steps are present: sulfide oxidation by sulfide:quinone oxidoreductase (`sqr`) and thiosulfate oxidation to tetrathionate by TsdA.
> >
> {: .solution}
>
{: .question}

## Transporters, enzymes, CAZymes and domains

The last sections are interactive Sankey diagrams from TC class to subclass and family, and from EC class down to individual enzymes with their ENZYME names. CAZy families are shown by class and grouped by the glycans they act on, and the most frequent Pfam domains and domain combinations are listed. The searchable gene table at the end links every summary back to individual proteins.

> <question-title></question-title>
>
> 1. In the transporter Sankey diagram, which TC class contains the most proteins in PAO1? Which family dominates it?
> 2. Type `pvd` into the search field of the gene table. What do the genes you find make, and why would a pathogen need it?
>
> > <solution-title></solution-title>
> >
> > 1. Class 3, *Primary active transporters* (470 proteins), dominated by the ATP-binding cassette (ABC) superfamily 3.A.1 with 278 proteins. It is followed by class 2, *Electrochemical potential-driven transporters* (214), and class 1, *Channels and pores* (135). Class 3 also contains the type II, III and VI secretion systems, which are classified as transporters in TCDB.
> > 2. The `pvd` genes (`pvdA`, `pvdE`, `pvdH`, `pvdI`, `pvdL`, `pvdS` and others), together with the receptor `fpvA`, make, export and take up pyoverdine, the main siderophore of *P. aeruginosa*: a fluorescent molecule that binds iron(III) with very high affinity. Inside a host, iron is kept bound to proteins such as transferrin; siderophores let the bacterium take it back, so pyoverdine production is linked to virulence. Note that several of these genes have generic descriptions (e.g. *Non-ribosomal peptide synthetase*), so a search by description alone would miss them.
> >
> {: .solution}
>
{: .question}

# Conclusion

In this tutorial we downloaded the proteome of *Pseudomonas aeruginosa* PAO1 from NCBI with a single accession, annotated all its proteins with eggNOG-mapper and summarised the annotation with Eggsplorer. The report recovered the known biology of this model bacterium: a Gram-negative envelope, a branched respiratory chain with five terminal oxidases, complete denitrification, assimilatory sulfate reduction, no carbon fixation, a flagellum and type I, II, III and VI secretion systems.

Keep in mind what this kind of analysis can and cannot tell you. The report shows the **potential** of the genome, inferred by orthology; it does not show which genes are expressed under a given condition. KO assignments can blur related reactions, as we saw for nitrate reductase and nitrite oxidoreductase. For higher confidence in KO assignments, Eggsplorer can combine eggNOG-mapper with KofamScan, which assigns KOs with profile HMMs and adaptive thresholds, and keep only the assignments both tools agree on (*consensus* mode).

The same workflow can be run on any bacterial proteome: one from NCBI, as here, or the proteins predicted by Bakta or Prokka for your own genome assembly. Because the report is a single file that works offline, it is easy to share with colleagues or to compare genomes side by side.

> <comment-title>Citing Eggsplorer</comment-title>
>
> If you use Eggsplorer in your work, please cite the Eggsplorer paper {% cite Baev2026 %} (submitted; the full reference will be added here once it is published), together with eggNOG-mapper {% cite Cantalapiedra2021 %}. If you use the metabolic function scores, please also cite KEGG-Decoder {% cite Graham2018 %} and KEGG {% cite Kanehisa2000 %}; for the CAZy-by-substrate chart, dbCAN3 {% cite Zheng2023 %}.
{: .comment}
