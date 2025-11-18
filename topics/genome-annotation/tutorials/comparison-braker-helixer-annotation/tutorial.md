---
layout: tutorial_hands_on

title: Comparison of two annotation tools - Helixer and Braker3
zenodo_link: https://zenodo.org/records/14770765
tags:
  - eukaryota
  - Braker3
  - Helixer
  - jbrowse1
  - biodiversity

edam_ontology:
- topic_0196 # Sequence Assembly
- topic_0622 # Genomics
- topic_0362 # Genome annotation


questions:
    - How to annotate an eukaryotic genome with Braker3?
    - How to annotate an eukaryotic genome with Helixer?
    - How to evaluate and visualize annotated genomic features?

objectives:
    - Load genome into Galaxy
    - Annotate genome with Braker3 and Helixer
    - Evaluate annotations quality with BUSCO
    - Compare the annotation of two annotation tools
    - View annotations in JBrowse

time_estimation: 8h
level: Intermediate
key_points:
    - Braker3 and Helixer allow to perform structural annotation of an eukaryotic genome
    - BUSCO allows to inspect the quality of an annotation
contributions:
  authorship:
    - rlibouba

  funding:
    - eurosciencegateway
abbreviations:
    GPU: Graphics Processing Unit
    UTR: Untranslated region

requirements:
 - type: internal
   topic_name: genome-annotation
   tutorials:
     - repeatmasker
     - helixer
     - braker3

subtopic: eukaryote
---

Annotating the eukaryotic genome represents a somewhat more complex challenge than that of prokaryotes, mainly due to the generally larger size of eukaryotic genomes and their greater number of genes. This annotation can be carried out at different levels of precision, ranging from simple identification of coding and non-coding parts to detailed structural labeling, including for example the precise location of exons, introns and other regulatory elements.

In this tutorial, we will use two software tools, Helixer and Braker3, to annotate the genome sequence of a small eukaryote: [*Mucor mucedo*](https://en.wikipedia.org/wiki/Mucor_mucedo) (a plant pathogenic fungus).
The idea of this tutorial is to compare annotation with two annotation tools that differ in their approach and operation.

[Helixer](https://github.com/weberlab-hhu/Helixer) is an annotation software with a new and different approach: it performs evidence-free predictions (no need for RNA-seq data or sequence alignments), using {GPU}, with a much faster execution time. The annotation is based on the development and use of a cross-species deep learning model. The software is used to configure and train models for *ab initio* prediction of gene structure. In other words, it identifies the base pairs in a genome that belong to the UTR/CDS/Intron genes.
There's a tutorial describing the steps involved in [annotating a genome with Helixer]({% link topics/genome-annotation/tutorials/helixer/tutorial.md %}).

[Braker3](https://github.com/Gaius-Augustus/BRAKER) is an automated bioinformatics tool that uses RNA-seq and protein sequences to annotate genome. It integrates GeneMark-ETP and AUGUSTUS software to predict genes with a high degree of precision. By combining the results of these two tools, Braker3 generates a final file containing gene annotations with strong extrinsic support (i.e. based on external experimental data). Braker3 facilitates genome annotation by leveraging transcriptomic and protein sequences to produce more reliable and robust gene predictions.
There's a tutorial describing the steps involved in [annotating a genome with Braker3]({% link topics/genome-annotation/tutorials/braker3/tutorial.md %}).

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data upload

To annotate our genome both using **Helixer** and **Braker3**, we will use the same genome sequence, in fasta format.

Unlike Helixer, **Braker3** needs two other files:
- Some **RNAseq data** in bam format. These reads are already aligned on the genome, and Braker3 will use it as evidences to annotate genes.
- A set of **protein sequences**, like UniProt/SwissProt. It is important to have good quality, curated sequences here, and the UniProt/SwissProt databank fits very well.

Galaxy uses Job Cache to avoid recalculating identical analyses. If a tool has already been run with the same parameters and inputs, previous results are reused, reducing waiting time and resource consumption.

> <hands-on-title>Using Job cache</hands-on-title>
>
> {% snippet faqs/galaxy/re_use_equivalent_jobs.md %}
>
{: .hands_on}

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    https://zenodo.org/records/15051665/files/genome_masked.fasta
>    https://zenodo.org/records/14770765/files/RNAseq_braker.bam
>    https://zenodo.org/records/14770765/files/swissprot_subset.fasta
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
{: .hands_on}

# Preparing RNAseq data for Braker3 annotation

In this tutorial, the alignments from RNA-seq are already prepared.
This data was obtained using the **RNA STAR** tool, as explained in detail in the [Braker3]({% link topics/genome-annotation/tutorials/braker3/tutorial.md %}) tutorial.

Note that when using Braker3 with RNA STAR to align RNA-seq data, the **--outSAMstrandField intronMotif** parameter must be added.
This parameter adds specific intron information to the alignment file (BAM).
This information is necessary for Braker3 to correctly understand and use the alignments to annotate genes.
Without this parameter, Braker3 may not function correctly or may produce incomplete results.

Here's how to run RNA STAR with this specific parameter (**"Read alignement tags to include in the BAM output"** option), before annotating with Braker3:

> <hands-on-title></hands-on-title>
>
> 1. {% tool [RNA STAR](toolshed.g2.bx.psu.edu/repos/iuc/rgrnastar/rna_star/2.7.11a+galaxy1) %} with the following parameters:
>    - *"Single-end or paired-end reads"*: `Paired-end (as individual datasets)`
>        - {% icon param-file %} *"RNA-Seq FASTQ/FASTA file, forward reads"*: `rnaseq_R1.fq.gz` (Input dataset)
>        - {% icon param-file %} *"RNA-Seq FASTQ/FASTA file, reverse reads"*: `rnaseq_R2.fq.gz` (Input dataset)
>    - *"Custom or built-in reference genome"*: `Use reference genome from history and create temporary index`
>        - {% icon param-file %} *"Select a reference genome"*: `genome_masked.fasta` (Input dataset)
>        - *"Length of the SA pre-indexing string"*: `11`
>    - In *"BAM output format specification"*:
>        - In *"Read alignement tags to include in the BAM output"*:
>             - Select `XS (strand flag, see parameter help below)`
>
{: .hands_on}

# Structural annotation
## By Helixer

We can run [**Helixer**](https://github.com/weberlab-hhu/Helixer) to perform the structural annotation of the genome.

Helixer requires only the genome sequence.
We also need to choose between 4 different lineages: *invertebrate*, *vertebrate*, *land plant* or *fungi*. Select the one that fits the best to the species you are studying: *fungi* in our case. Helixer is shipped with these 4 models that were trained specifically to annotate genes from each of these lineages. Advanced users can upload their own lineage model in .h5 format with the *"Lineage model"* option.

As an option, we can also enter a species name.


> <hands-on-title>Genome annotation with Helixer</hands-on-title>
>
> 1. {% tool [Helixer](toolshed.g2.bx.psu.edu/repos/genouest/helixer/helixer/0.3.3+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Genomic sequence"*: `genome_masked.fasta` (Input dataset)
>    - In *"Available lineages"*: "*select*" `fungi`
>    - In *"Species name"*: `Mucor mucedo`
>
{: .hands_on}

> <comment-title>Advanced parameters</comment-title>
>
> Depending on the lineage, the parameters *"Subsequence length"*, *"Overlap offset"* and *"Overlap corelength"* are adjusted to corresponding default values (listed in the help of each option).
>
> This is due in particular to the size of the genomes. Indeed, it is recommended to increase the value of *"Subsequence length"* for genomes containing large genes. This is particularly important for vertebrates and invertebrates.
>
> The default values used by Galaxy are the ones recommended by Helixer authors. If you wish to modify these default values, you can do so by entering your values in the *"Subsequence length"*, *"Overlap offset"* and *"Overlap corelength"* parameters.
>
{: .comment}

> <comment-title>Don't wait</comment-title>
>
> This step can take a bit of time to run: although Helixer runs much faster than many other annotation tools (typically <20min for this tutorial), it requires specific hardware (GPU) that is often available in limited quantity on computing systems. It means your job can be placed in queue for a longer time than a more standard Galaxy job.
>
> While it runs, we can already schedule the following steps. Galaxy will run them automatically as soon as the Helixer annotation is ready.
{: .comment}

Helixer produces a single output dataset: a GFF3 file. The GFF3 format is a standard bioinformatics format for storing genome annotations. Each row describes a genomic entity, with columns detailing its identifier, location, score and other attributes.

## By Braker3

We have the following required inputs:

- A masked genome in FASTA format
- RNAseq sequence alignments in BAM format
- Protein annotation in FASTA format

With these files, we can run [**Braker3**](https://github.com/Gaius-Augustus/BRAKER) to perform the structural annotation of the genome.

> <hands-on-title>Genome annotation with Braker3</hands-on-title>
>
> 1. {% tool [Braker3](toolshed.g2.bx.psu.edu/repos/iuc/braker3/braker3/3.0.8+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Assembly to annotate"*: `genome_masked.fasta` (Input dataset)
>    - *"Species name"*: `Mucor mucedo`
>    - {% icon param-file %} *"RNA-seq mapped to genome to train Augustus/GeneMark"*: `RNASeq_braker.bam` (Input dataset)
>    - {% icon param-file %} *"Proteins to map to genome"*: `SwissProt_subset.fasta` (Input dataset)
>    - {% icon param-file %} *"Fungal genome"*: select `Yes`
>    - *"Output format"*: `GFF3`
>
{: .hands_on}

> <comment-title>On parameters</comment-title>
>
> - If you are working with fungi, we recommend that you activate the **--fungus** parameter. This enables Braker3 to run the GeneMark-ETP algorithm with a branch point model.
{: .comment}

> <comment-title>Don't wait</comment-title>
>
> - This step will take a long time to complete: it could take several hours for this tutorial. While it's running, we can already schedule the following steps. Galaxy will automatically execute them as soon as the structural annotation is ready.
> - There are several reasons for this: the size of the genome, the quality and quantity of the RNA-Seq data, the supervised learning model for prediction, and more.
{: .comment}

Braker3 generates an output file in GTF format by default. But it is possible to generate the output file in GFF3 output. This file contains informations on gene locations, exons, introns, etc.
GFF3 files are generally preferred for complex and varied annotations, while GTF is often used for transcript annotations in RNA-seq studies.

If you can't wait to have the Braker3 annotation and move on to the comparison step between our two annotation tools, we invite you to import the Braker3 annotation:

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    https://zenodo.org/records/14770765/files/braker3.gff3
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
{: .hands_on}

# Quality evaluation

## Evaluation with **Busco**

[BUSCO](http://busco.ezlab.org/) (Benchmarking Universal Single-Copy Orthologs) is a widely used tool to evaluate the quality of a genome assembly and annotation. By comparing genomes from various related and distantly related species, the authors determined sets of ortholog genes that are present in single copy in (almost) all the species of a clade (Bacteria, Fungi, Plants, Insects, Mammalians, …).
Most of these genes are essential for the organism to live, and are expected to be found in any newly sequenced and annotated genome from the corresponding clade. Using this data, BUSCO evaluates the "completeness" of genome annotation by assessing the proportion of these essential genes (also named BUSCOs) found in a set of (predicted) transcript or protein sequences.

### Braker3 annotation

We want to run BUSCO on the protein sequences predicted from gene sequences of the Braker3 annotation.

So first generate these sequences:

> <hands-on-title>Extract protein sequences</hands-on-title>
>
> 1. {% tool [GFFread](toolshed.g2.bx.psu.edu/repos/devteam/gffread/gffread/2.2.1.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input GFF3 or GTF feature file"*: output of {% tool [Braker3](toolshed.g2.bx.psu.edu/repos/iuc/braker3/braker3/3.0.8+galaxy0) %}
>    - In *"Reference Genome"* select: `From your history` (Input dataset)
>    - *"Genome Reference Fasta"*: `genome_masked.fasta` (Input dataset)
>    - In *"Select fasta outputs"* select: `fasta file with spliced exons for each GFF transcript (-y)`
>    - *"full GFF attribute preservation (all attributes are shown)"*: `Yes`
>    - *"decode url encoded characters within attributes"*: `Yes`
>    - *"warn about duplicate transcript IDs and other potential problems with the given GFF/GTF records"*: `Yes`
>
{: .hands_on}

The parameters for running BUSCO on these protein sequences:

> <hands-on-title>BUSCO in proteome mode</hands-on-title>
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `gffread: pep.fa`
>    - *"Mode"*: `annotated gene sets (protein)`
>    - *"Auto-detect or select lineage?"*: `Select lineage`
>        - *"Lineage"*: `Mucorales`
>    - In *"Advanced Options"*:
>        - *"Which outputs should be generated"*: select all outputs
{: .hands_on}

Several output files are generated:

- short summary: statistical summary of the quality of genomic assembly or annotation, including total number of genes evaluated, percentage of complete genes, percentage of partial genes, etc.
- full table: list of universal orthologs found in the assembled or annotated genome, with information on their completeness, location in the genome, quality score, etc.
- missing buscos: list of orthologs not found in the genome, which may indicate gaps in assembly or annotation.
- summary image: graphics and visualizations to visually represent the results of the evaluation, such as bar charts showing the proportion of complete, partial and missing genes.
- GFF: contain information on gene locations, exons, introns, etc.

This gives information about the completeness of the Braker3 annotation. A good idea is to compare this first result with the one you get on the initial genome sequence, and see if the annotation tool found all the genes that BUSCO finds in the raw genome sequence. To do this, you can restart BUSCO on the masked genome with the same parameters, but selecting the **Genome assemblies (DNA)** mode.

> <comment-title>What can we deduce from these results?</comment-title>
>
> - 94.4% of genes are complete, so the annotation is of high quality in terms of genomic completeness.
> - The duplication rate is low, with 7.6% of genes duplicated which may indicate either real biological duplications or annotation errors.
> - So the Braker3 annotation looks like a good one, with high completeness and low duplication.
{: .comment}

### Helixer annotation

We are going to do the same thing again, but with the Helixer annotation. The first step is to generate the gene sequences predicted by the annotation:

> <hands-on-title>Extract protein sequences</hands-on-title>
>
> 1. {% tool [GFFread](toolshed.g2.bx.psu.edu/repos/devteam/gffread/gffread/2.2.1.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input GFF3 or GTF feature file"*: output of {% tool [Helixer](toolshed.g2.bx.psu.edu/repos/genouest/helixer/helixer/0.3.3+galaxy1) %}
>    - In *"Reference Genome"* select: `From your history` (Input dataset)
>    - *"Genome Reference Fasta"*: `genome_masked.fasta` (Input dataset)
>    - In *"Select fasta outputs"* select: `protein fasta file with the translation of CDS for each record (-y)`
>    - *"full GFF attribute preservation (all attributes are shown)"*: `Yes`
>    - *"decode url encoded characters within attributes"*: `Yes`
>    - *"warn about duplicate transcript IDs and other potential problems with the given GFF/GTF records"*: `Yes`
>
{: .hands_on}

To run BUSCO on these protein sequences predicted by Helixer annotation:

> <hands-on-title>BUSCO in proteome mode</hands-on-title>
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `gffread: pep.fa`
>    - *"Mode"*: `annotated gene sets (protein)`
>    - *"Auto-detect or select lineage?"*: `Select lineage`
>        - *"Lineage"*: `Mucorales`
>    - In *"Advanced Options"*:
>        - *"Which outputs should be generated"*: select all outputs
>
{: .hands_on}

> <comment-title>What can we deduce from these results?</comment-title>
>
> - 94.6% of genes are complete, so the annotation is of high quality in terms of genomic completeness.
> - The duplication rate is low, with 1.3% of genes duplicated.
> - So the Helixer annotation looks like a good one, with high completeness and low duplication.
>
{: .comment}

### Compare the two annotations with BUSCO

Braker3 identifies 94.4% complete genes, while Helixer detects 94.6%. Helixer therefore finds slightly more complete genes.

As for duplications, Braker3 detects 7.6%, compared with only 1.3% for Helixer. Helixer's annotation is therefore more precise, with fewer redundant copies.

These differences can be explained by the methods used by each tool to annotate genes.

- BRAKER3 uses an approach based on biological evidence, combining RNA-Seq data and protein alignments. This enables it to identify genes with good accuracy, but can also lead to duplications if several isoforms are poorly distinguished.
- Helixer, on the other hand, relies on a deep learning model that directly analyzes the genomic sequence. This approach delivers a more consistent prediction and reduces duplications, but may be limited if the model has not been well trained on similar genomes.

## Evaluation with **OMArk**

[OMArk](https://github.com/DessimozLab/OMArk) is proteome quality assessment software. It provides measures of proteome completeness, characterises the consistency of all protein-coding genes with their homologues and identifies the presence of contamination by other species. OMArk is based on the OMA orthology database, from which it exploits orthology relationships, and on the OMAmer software for rapid placement of all proteins in gene families.

OMArk's analysis is based on HOGs (Hierarchical Orthologous Groups), which play a central role in its assessment of the completeness and coherence of gene sets. HOGs make it possible to compare the genes of a given species with groups of orthologous genes conserved across a taxonomic clade.

### Braker3 annotation

> <hands-on-title>OMArk on extracted protein sequences</hands-on-title>
>
> 1. {% tool [OMArk](toolshed.g2.bx.psu.edu/repos/iuc/omark/omark/0.3.0+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Protein sequences"*: `gffread: pep.fa` (Braker3 annotation)
>    - *"OMAmer database*: select `LUCA-v2.0.0`
>    - In *"Which outputs should be generated"*: select `Detailed summary`
>
{: .hands_on}

> <comment-title>What can we deduce from these results?</comment-title>
>
> - Number of conserved HOGs: OMArk has identified a set of 5622 HOGs which are thought to be conserved in the majority of species in the Mucorineae clade.
> - 78.66% of genes are complete, so the annotation is of good quality in terms of genomic completeness.
> - Number of proteins in the whole proteome: 14 209. Of which 77.59% are present and 15.69% of the proteome does not share sufficient similarities with known gene families.
> - No contamination detected.
> - The OMArk analysis is based on the Mucorineae lineage, a more recent and specific clade than that used in the BUSCO assessment, which selected the Mucorales as the reference group.
{: .comment}

### Helixer annotation

> <hands-on-title>OMArk on extracted protein sequences</hands-on-title>
>
> 1. {% tool [OMArk](toolshed.g2.bx.psu.edu/repos/iuc/omark/omark/0.3.0+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Protein sequences"*: `gffread: pep.fa` (Helixer annotation)
>    - *"OMAmer database*: select `LUCA-v2.0.0`
>    - In *"Which outputs should be generated"*: select `Detailed summary`
>
{: .hands_on}

> <comment-title>What can we deduce from these results?</comment-title>
>
> - Number of conserved HOGs: OMArk has identified a set of 5622 HOGs which are thought to be conserved in the majority of species in the Mucorineae clade.
> - 85.52% of genes are complete, so the annotation is of good quality in terms of genomic completeness.
> - Number of proteins in the whole proteome: 19 299. Of which 62.83% are present and 30.94% of the proteome does not share sufficient similarities with known gene families.
> - No contamination detected.
> - The OMArk analysis is based on the Mucorineae lineage, a more recent and specific clade than that used in the BUSCO assessment, which selected the Mucorales as the reference group.
{: .comment}

### Compare the two annotations with OMArk

Both annotations identify the same set of 5622 HOGs conserved in the majority of species in the Mucorineae clade. Furthermore, no contamination is detected in either case, meaning that no protein has been incorrectly assigned to another lineage.

The annotation obtained with Braker3 displays a percentage of 78.66% complete genes, indicating good annotation quality in terms of genomic completeness. Helixer, on the other hand, shows higher completeness, with 85.52% complete genes, suggesting better genomic coverage.

With regard to the total proteome, Braker3 identified 14,209 proteins, 77.59% of which are consistent with known gene families, while 17.35% of proteins do not share sufficient similarities with these families. Helixer, on the other hand, generates a higher number of proteins (19,299), but only 62.83% of proteins are consistent with known gene families, and 30.94% do not share sufficient similarities. This suggests that Helixer generates a greater number of proteins, but a higher proportion have no clear correspondence with existing gene families.

With our example dataset, Braker3 seems to be a more suitable tool if you're looking for accurate and consistent annotation, with more reliable proteins. Helixer, on the other hand, offers a more complete coverage of the genome, but with a higher number of inconsistent proteins, which may require refinement of the results to improve annotation quality.

# Visualization with a genome browser

You can visualize the annotations generated using a genomic browser like [JBrowse](https://jbrowse.org/jbrowse1.html). This browser enables you to navigate along the chromosomes of the genome, to view the structure of each predicted gene, as well as to compare the genes predicted by each method.

> <hands-on-title>JBrowse visualisation</hands-on-title>
>
> 1. {% tool [JBrowse](toolshed.g2.bx.psu.edu/repos/iuc/jbrowse/jbrowse/1.16.11+galaxy1) %} with the following parameters:
>    - *"Reference genome to display"*: `Use a genome from history`
>        - {% icon param-file %} *"Select the reference genome"*: `genome_masked.fasta` (Input dataset)
>    - In *"Track Group"*:
>        - {% icon param-repeat %} *"Insert Track Group"*
>            - *"Track Category"*: `Braker3 nnotation`
>            - In *"Annotation Track"*:
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `GFF/GFF3/BED Features`
>                        - {% icon param-file %} *"GFF/GFF3/BED Track Data"*: `gff3` (output of **Braker3** {% icon tool %})
>    - In *"Track Group"*:
>        - {% icon param-repeat %} *"Insert Track Group"*
>            - *"Track Category"*: `Helixer nnotation`
>            - In *"Annotation Track"*:
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `GFF/GFF3/BED Features`
>                        - {% icon param-file %} *"GFF/GFF3/BED Track Data"*: `gff3` (output of **Helixer** {% icon tool %})
{: .hands_on}

Click on the newly created dataset's eye to display it. You will see a JBrowse genome browser. You can have a look at the [JBrowse tutorial]({% link topics/visualisation/tutorials/jbrowse/tutorial.md %}) for a more in-depth description of JBrowse.

We have embedded a copy of the resulting JBrowse here, if something went wrong during one of the steps you can always just check this output:

{% capture jbtracks %}DNA%2Caa6aca08d6fd82b94c7a86712a632760_0{% endcapture %}
{% snippet topics/visualisation/faqs/visualizations_jbrowse.html datadir="data" tracks=jbtracks %}

Check that **Braker3 annotation** and **Helixer annotation** are selected in the top left-hand corner. If not, select it to view the annotation.

As confirmed by OMArk, the Braker3 annotation shows a higher duplication rate (7.8%). On JBrowse, these duplications appear as repeated genes in close proximity to each other. If you haven't seen them yet, we recommend you explore scaffold 1 at **scaffold_1:536081..544070**.

Annotation with Helixer seems to generate longer genes than Braker3. This may be explained by its annotation method based on deep neural networks. These longer genes could correspond to extra exons or poorly predicted introns. If you'd like to see an example of this, please visit scaffold 2 at **scaffold_2:494676..508755**.

Helixer also detects a larger number of genes. OMArk analysis showed a higher quantity of proteins in the Helixer annotation (19,299 proteins) compared to Braker3 (14,209 proteins). If you wish to observe this difference, we recommend you explore scaffold 1 at **scaffold_1:585671..591035** or scaffold 13 at **scaffold_13:100856..107419**.

# Conclusion

Congratulations on reaching the end of this tutorial! You now know how to perform a structural annotation of a new eukaryotic genome, using Helixer and Braker3. And you've learned how to evaluate its quality and how to visualize it using JBrowse.

Through this comparative analysis, you've seen that both annotation tools have their strengths with this example dataset:

* Helixer offers a high level of gene completeness and avoids redundant duplications thanks to its deep learning-based approach. It's particularly effective when a consistent, duplication-free prediction is preferred.
* Braker3, leveraging RNA-Seq and protein evidence, delivers more biologically grounded annotations, with a better match to known gene families, making it a solid choice when biological accuracy and interpretability are key. As Braker3 relies on evidence data, the quality of the results depends on the quality and quantity of available data.

When annotating a new genome sequence, it is a good idea to use different methods (Helixer, Braker3, but possibly others like Funannotate, Maker, Braker2, ...), and to perform this kind of comparison to determine the best result.
