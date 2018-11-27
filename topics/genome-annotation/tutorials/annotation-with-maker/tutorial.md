---
layout: tutorial_hands_on

title: Genome annotation with Maker
zenodo_link: https://doi.org/10.5281/zenodo.1488687
tags:
  - eukaryote
questions:
  - How to annotate an eukaryotic genome?
  - How to evaluate and visualize annotated genomic features?
objectives:
  - Load genome into Galaxy
  - Annotate genome with Maker
  - Evaluate annotation quality with BUSCO
  - View annotations in JBrowse
time_estimation: 4h
key_points:
  - Maker allows to annotate a eukaryotic genome.
  - BUSCO and JBrowse allow to inspect the quality of an annotation.
contributors:
  - abretaud
---


# Introduction
{:.no_toc}

Genome annotation of eukaryotes is a little more complicated than for prokaryotes: eukaryotic genomes are usually larger than prokaryotes, with more genes. The sequences determining the beginning and the end of a gene are generally less conserved than the prokaryotic ones. Many genes also contain introns, and the limits of these introns (acceptor and donor sites) are not highly conserved.

In this tutorial we will use a software tool called Maker to annotate the genome sequence of a small eukaryote: Schizosaccharomyces pombe (a yeast).

Maker is able to annotate both prokaryotes and eukaryotes. It works by aligning as many evidences as possible along the genome sequence, and then reconciliating all these signals to determine probable gene structures.

The evidences can be transcript or protein sequences from the same (or closely related) organism. These sequences can come from public databases (like NR or GenBank) or from your own experimental data (transcriptome assembly from an RNASeq experiment for example). Maker is also able to take into account repeated elements.

Maker uses ab-initio predictors (like [Augustus](http://bioinf.uni-greifswald.de/augustus/) or [SNAP](https://github.com/KorfLab/SNAP)) to improve its predictions: these software tools are able to make gene structure predictions by analysing only the genome sequence with a statistical model.

In this tutorial you will learn how to perform a genome annotation, and how to evaluate its quality. You will see how training ab-initio predictors is an important step to produce good results. Finally, you will learn how to use the [JBrowse](http://jbrowse.org/) genome browser to visualise the results.

More information about Maker can be found [here](http://www.yandell-lab.org/software/maker.html).

This tutorial was inspired by the MAKER Tutorial for [WGS Assembly and Annotation Winter School 2018 ](http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018), don't hesitate to consult it for more information on Maker, and on how to run it with command line.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data upload

To annotate a genome using Maker, you need the following files:

- The genome sequence in fasta format
- A set of transcripts or EST sequences, preferably from the same organism.
- A set of protein sequences, usually from closely related species or from a curated sequence database like UniProt/SwissProt.

 Maker will align the transcript and protein sequences on the genome sequence to determine gene positions.

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create and name a new history for this tutorial.
> 2. Import the following files from [Zenodo](https://doi.org/10.5281/zenodo.1488687) or from the shared data library
>
>    ```
>    https://zenodo.org/api/files/4385871d-9632-4fae-9aaf-f8ed692163d1/augustus_training_1.tar.gz
>    https://zenodo.org/api/files/4385871d-9632-4fae-9aaf-f8ed692163d1/augustus_training_2.tar.gz
>    https://zenodo.org/api/files/4385871d-9632-4fae-9aaf-f8ed692163d1/S_pombe_chrIII.fasta
>    https://zenodo.org/api/files/4385871d-9632-4fae-9aaf-f8ed692163d1/S_pombe_genome.fasta
>    https://zenodo.org/api/files/4385871d-9632-4fae-9aaf-f8ed692163d1/S_pombe_trinity_assembly.fasta
>    https://zenodo.org/api/files/4385871d-9632-4fae-9aaf-f8ed692163d1/Swissprot_no_S_pombe.fasta
>    ```
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype for `augustus_training_1.tar.gz` an `augustus_training_2.tar.gz` are set to `augustus`
>
>    {% include snippets/change_datatype.md datatype="augustus" %}
>
{: .hands_on}

You have four main datasets:

- `S_pombe_trinity_assembly.fasta` contains the EST sequences
- `Swissprot_no_S_pombe.fasta` contanis the protein sequences from SwissProt
- `S_pombe_genome.fasta` contains the full genome sequence
- `S_pombe_chrIII.fasta` contains only a fraction of the full genome sequence, ie the chromosome III

For the rest of this tutorial, you need to choose between `S_pombe_chrIII.fasta` and `S_pombe_genome.fasta`. If you have time, use the full genome (`S_pombe_genome.fasta`), it will take more computing time, but the results will be closer to real-life data. If you want to get results faster, use the chromosome III (`S_pombe_chrIII.fasta`). In the rest of this tutorial, we will refer to the file you choose as `the genome`.

The two other datasets (`augustus_training_1.tar.gz` an `augustus_training_2.tar.gz`) will be used later in the tutorial.

# Genome quality evaluation

The quality of a genome annotation is highly dependent on the quality of the genome sequences. It is impossible to obtain a good quality annotation with a poorly assembled genome sequence. Annotation tools will have trouble finding genes if the genome sequence is highly fragmented, if it contains chimeric sequences, or if there are a lot of sequencing errors.

Before running the full annotation process, you need first to evaluate the quality of the sequence. It will give you a good idea of what you can expect from it at the end of the annotation.

> ### {% icon hands_on %} Hands-on: Get genome sequence statistics
>
> 1. **Fasta Statistics** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: select the genome sequence from your history
>
{: .hands_on}

Have a look at the statistics:

- `num_seq`: the number of contigs (or scaffold or chromosomes), compare it to expected chromosome numbers
- `len_min`, `len_max`, `len_N50`, `len_mean`, `len_median`: the distribution of contig sizes
- `num_bp_not_N`: the number of bases that are not N, it should be as close as possible to the total number of bases (`num_bp`)

These statistics are useful to detect obvious problems in the genome assembly, but it gives no information about the quality of the sequence content. We want to evaluate if the genome sequence contains all the genes we expect to find in the considered species, and if their sequence are correct.

[BUSCO](http://busco.ezlab.org/) (Benchmarking Universal Single-Copy Orthologs) is a tool allowing to answer this question: by comparing genomes from various more or less related species, the authors determined sets of ortholog genes that are present in single copy in (almost) all the species of a clade (Bacteria, Fungi, Plants, Insects, Mammalians, ...). Most of these genes are essential for the organism to live, and are expected to be found in any newly sequenced genome from the corresponding clade. Using this data, BUSCO is able to evaluate the proportion of these essential genes (also named BUSCOs) found in a genome sequence or a set of (predicted) transcript or protein sequences. This is a good evluation of the "completeness" of the genome or annotation.

We will first run this tool on the genome sequence to evaluate its quality.


> ### {% icon hands_on %} Hands-on: Run Busco on the genome
>
> 1. **Busco** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: select the genome sequence from your history
>    - *"Mode"*: `Genome`
>    - *"Lineage"*: `fungi_odb9`
>
>    > ### {% icon comment %} Comment
>    >
>    > We select fungi_odb9 as the genome we will annotate is a fungi. It is usually better to select the most specific lineage for the species you study. Large lineages (like Metazoa) will consist of fewer genes, but with a strong support. More specific lineages (like Hymenoptera) will have more genes, but with a weaker support (has they are found in fewer genomes).
>    {: .comment}
>
{: .hands_on}

BUSCO produces three output datasets

- A short summary: summarizes the results of BUSCO (see below)
- A full table: lists all the BUSCOs that were searched for, with the corresponding status (was it found in the genome? how many times? where?)
- A table of missing BUSCOs: this is the list of all genes that were not found in the genome

>    ![BUSCO genome summary](../../images/busco_genome_summary.png "BUSCO short summary. Here BUSCO searched for 290 genes in a genome, and found 282 of them complete (269 in single-copy, and 13 duplicated). 5 where found but are fragmented, and 3 were not found at all.")

> ### {% icon question %} Questions
>
> 1. Do you think the genome quality is good enough for performing the annotation?
>
> > ### {% icon solution %} Solution
> >
> > 1. The genome consists of the exepected 4 chromosomes sequences, with very few N, which is the ideal case. Most of the BUSCO genes are found as complete single copy, and very few are fragmented, which means that our genome have a good quality as it contains most of the expected content. That's a very good material to perform an annotation.
> >
> {: .solution}
>
{: .question}

> ### {% icon comment %} Comment
>
> If you chose to use only the chromosome III sequence (`S_pombe_chrIII.fasta`), the statistics will be different. The genome size will be lower, with only 1 chromosome. The BUSCO result will also show a lot of missing genes: it is expected as all the BUSCO genes that are not on the chromosome III cannot be found by the tool. Keep it in mind when comparing these results with the other BUSCO results later.
{: .comment}

# First Maker annotation round

## Maker

For this first round, we configure Maker to construct gene models only by aligning ESTs and proteins to the genome. This will produce a first draft annotation that we will improve in the next steps.

> ### {% icon hands_on %} Hands-on: First draft annotation with Maker
>
> 1. **Maker** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select the genome sequence from your history
>    - *"Re-annotate using an existing Maker annotation"*: `No`
>    - *"Organism type"*: `Eukaryotic`
>    - In *"EST evidences (for best results provide at least one of these)"*:
>        - *"Infer gene predictions directly from all ESTs"*: `Yes`
>        - {% icon param-file %} *"ESTs or assembled cDNA"*: `S_pombe_trinity_assembly.fasta`
>    - In *"Protein evidences (for best results provide at least one of these)"*:
>        - *"Infer gene predictions directly from all protein alignments"*: `Yes`
>        - {% icon param-file %} *"Protein sequences"*: `Swissprot_no_S_pombe.fasta`
>    - In *"Ab-initio gene prediction"*:
>        - *"Prediction with Augustus"*: `Don't use Augustus to predict genes`
>    - In *"Repeat masking"*:
>        - *"Repeat library source"*: `Disable repeat masking (not recommended)`
>
>    > ### {% icon comment %} Comment
>    >
>    > For this tutorial repeat masking is disabled, which is not the recommended setting. When doing a real-life annotation, you should either provide your own repeats library, or use [RepBase](https://www.girinst.org/repbase/). As RepBase is not free to use, you will need to download manually the RepBaseRepeatMaskerEdition file and upload the enclosed EMBL format file (`RMRBSeqs.embl`) to use with Maker.
>    {: .comment}
>
{: .hands_on}

Maker produces three GFF3 datasets:

- The final annotation: the final consensus gene models produced by Maker
- The evidences: the alignements of all the data Maker used to construct the final annotation (ESTs and proteins that we used)
- A GFF3 file containing both the final annotation and the evidences

## Annotation statistics

We need now to evaluate this first draft annotation produced by Maker.

First, use the `Genome annotation statistics` that will compute some general statistics on the annotation.

> ### {% icon hands_on %} Hands-on: Get annotation statistics
>
> 1. **Genome annotation statistics** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotation to analyse"*: `final annotation` (output of **Maker** {% icon tool %})
>    - *"Reference genome"*: `Use a genome from history`
>        - {% icon param-file %} *"Corresponding genome sequence"*: select the genome sequence from your history
>
>
{: .hands_on}


## Busco

Just as we did for the genome at the beginning, we can use BUSCO to check the quality of this first Maker annotation. Instead of looking for known genes in the genome sequence, BUSCO will inspect the transcript sequences of the genes predicted by Maker.

First we need to compute all the transcript sequences from the Maker annotation.

> ### {% icon hands_on %} Hands-on: Extract transcript sequences
>
> 1. **gffread** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input GFF3 or GTF feature file"*: `final annotation` (output of **Maker** {% icon tool %})
>    - *"Reference Genome"*: `select the genome sequence from your history`
>        - *"Select fasta outputs"*:
>           - `fasta file with spliced exons for each GFF transcript (-w exons.fa)`
>    - *"full GFF attribute preservation (all attributes are shown)"*: `Yes`
>    - *"decode url encoded characters within attributes"*: `Yes`
>    - *"warn about duplicate transcript IDs and other potential problems with the given GFF/GTF records"*: `Yes`
>
{: .hands_on}

Now run BUSCO with the predicted transcript sequences:

> ### {% icon hands_on %} Hands-on: Run BUSCO
>
> 1. **Busco** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `exons` (output of **gffread** {% icon tool %})
>    - *"Mode"*: `Transcriptome`
>    - *"Lineage"*: `fungi_odb9`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. How many genes and transcripts where predicted by Maker?
> 2. What is the medium length of these genes?
> 3. How do the BUSCO statistics compare to the ones at the genome level?
>
> > ### {% icon solution %} Solution
> >
> > These numbers are for the full genome only (not the chromosome III).
> > 1. 713 genes, 752 transcripts
> > 2. 993 bp
> > 3. 40 complete single-copy, 6 duplicated, 9 fragmented, 235 missing. This is far from what BUSCO found in the genome sequence, which means the quality of this first draft is not very good.
> >
> {: .solution}
>
{: .question}


The statistics are not really satisfactory at this stage, but it's normal: Maker only used the EST and protein evidences to guess the gene positions. Let's see now how this first draft can be improved.

# Ab-initio predictors first training

Maker can use several ab-initio predictors to annotate a genome. "Ab-initio" means that these predictors are able to predict the structure of genes in a genome based only on its sequence and on a species-specific statistical model. They don't use any evidence (e.g. EST or proteins) to predict genes, but they need to be trained with a set of already predicted genes.

Maker is able to use the EST and protein evidences, and to combine them with the result of several ab-initio predictors to predict consensus gene models. It allows to detect genes in regions where no EST or protein align, and also to refine gene structures in regions where there are EST and/or protein evidences and ab-initio predictions.

We will use two of the most widely used ab-initio predictors: SNAP and Augustus. Before using it within Maker, we need to train them with the first annotation draft we produced in the previous steps. We know the quality of this draft is not perfect, but only the best scoring genes (ie the ones having the strongest evidences) will be retained to train the predictors.

> ### {% icon hands_on %} Hands-on: Train SNAP
>
> 1. **Train SNAP** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select the genome sequence from your history
>    - {% icon param-file %} *"Maker annotation to use for training"*: `final annotation` (output of **Maker** {% icon tool %})
>
{: .hands_on}


Augustus is trained in a very similar way.

> ### {% icon hands_on %} Hands-on: Train Augustus
>
> 1. **Train Augustus** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select the genome sequence from your history
>    - {% icon param-file %} *"Annotation to use for training"*: `final annotation` (output of **Maker** {% icon tool %})
>
{: .hands_on}

Both SNAP and Augustus produce a statistical model representing the observed general structure of genes in the analysed genome. Maker will use these models to create a new annotation for our genome.

The Augustus training usually take around 2 hours to complete, to continue this tutorial without waiting for the result, you can use the file `augustus_training_1.tar.gz` imported from Zenodo.

# Second Maker annotation round

## Maker

We need now to run a new round of Maker. As the evidences were already aligned on the genome on the first run, we can reuse these alignments as is.
But this time, enable ab-initio gene prediction, and input the output of **Train SNAP** {% icon tool %} and **Train Augustus** {% icon tool %} tools.
We also disable infering gene predictions directly from all ESTs and proteins: now we want Maker to infer gene predictions by reconciliating evidence alignments *and* ab-initio gene predictions.

> ### {% icon hands_on %} Hands-on: Second draft annotation with Maker
>
> 1. **Maker** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select the genome sequence from your history
>    - *"Organism type"*: `Eukaryotic`
>    - *"Re-annotate using an existing Maker annotation"*: `Yes`
>        - {% icon param-file %} *"Previous Maker annotation"*: `evidences` (output of the previous **Maker** {% icon tool %} run)
>        - *"Re-use ESTs"*: `Yes`
>        - *"Re-use alternate organism ESTs"*: `Yes`
>        - *"Re-use protein alignments"*: `Yes`
>        - *"Re-use repeats"*: `Yes`
>    - In *"EST evidences (for best results provide at least one of these)"*:
>        - *"Infer gene predictions directly from all ESTs"*: `No`
>    - In *"Protein evidences (for best results provide at least one of these)"*:
>        - *"Infer gene predictions directly from all protein alignments"*: `No`
>    - In *"Ab-initio gene prediction"*:
>        - *"SNAP model"*: `snap model` (output of **Train SNAP** {% icon tool %})
>        - *"Prediction with Augustus"*: `Run Augustus with a custom prediction model`
>            - {% icon param-file %} *"Augustus model"*: `augustus model` (output of **Train Augustus** {% icon tool %})
>    - In *"Repeat masking"*:
>        - *"Repeat library source"*: `Disable repeat masking (not recommended)`
>
{: .hands_on}


## Annotation statistics

Do we get a better result from Maker after this second run? Let's run the same tools as after the first run, and compare the results.

> ### {% icon hands_on %} Hands-on: Get annotation statistics
>
> 1. **Genome annotation statistics** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotation to analyse"*: `final annotation` (output of **Maker** {% icon tool %} second run)
>    - *"Reference genome"*: `Use a genome from history`
>        - {% icon param-file %} *"Corresponding genome sequence"*: select the genome sequence from your history
>
>
{: .hands_on}

## Busco

> ### {% icon hands_on %} Hands-on: Extract transcript sequences
>
> 1. **gffread** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input GFF3 or GTF feature file"*: `final annotation` (output of **Maker** {% icon tool %} second run)
>    - *"Reference Genome"*: `select the genome sequence from your history`
>        - *"Select fasta outputs"*:
>           - `fasta file with spliced exons for each GFF transcript (-w exons.fa)`
>    - *"full GFF attribute preservation (all attributes are shown)"*: `Yes`
>    - *"decode url encoded characters within attributes"*: `Yes`
>    - *"warn about duplicate transcript IDs and other potential problems with the given GFF/GTF records"*: `Yes`
>
{: .hands_on}

Now run BUSCO with the predicted transcript sequences:

> ### {% icon hands_on %} Hands-on: Run BUSCO
>
> 1. **Busco** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `exons` (output of **gffread** {% icon tool %})
>    - *"Mode"*: `Transcriptome`
>    - *"Lineage"*: `fungi_odb9`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. How does the second annotation compare to the previous one? Did the ab-initio predictors training improve the results?
> 2. How do you explain these changes?
>
> > ### {% icon solution %} Solution
> >
> > 1. The annotation looks much better: 253 complete single-copy instead of 40, 5,036 genes instead of 713.
> > 2. Using ab-initio predictors allowed to find much more genes in regions where EST or protein alignments were not sufficient to predict genes.
> >
> {: .solution}
>
{: .question}


# Ab-initio predictors second training

To get better results, we are going to perform a second training of SNAP and Augustus, and then run Maker for a third (final) time.

> ### {% icon hands_on %} Hands-on: Train SNAP and Augustus
>
> 1. **Train SNAP** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select the genome sequence from your history
>    - {% icon param-file %} *"Maker annotation to use for training"*: `final annotation` (output of **Maker** {% icon tool %}, second run)
>    - *"Number of gene model to use for training"*: `"1000"`
>
> 2. **Train Augustus** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select the genome sequence from your history
>    - {% icon param-file %} *"Annotation to use for training"*: `final annotation` (output of **Maker** {% icon tool %}, second run)
>
{: .hands_on}

The Augustus training usually take around 2 hours to complete, to continue this tutorial without waiting for the result, you can use the file `augustus_training_2.tar.gz` imported from Zenodo.

# Third (last) Maker annotation round

## Maker

Let's run the final round of Maker, in the same way as we did for the second run.

> ### {% icon hands_on %} Hands-on: Final annotation with Maker
>
> 1. **Maker** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select the genome sequence from your history
>    - *"Organism type"*: `Eukaryotic`
>    - *"Re-annotate using an existing Maker annotation"*: `Yes`
>        - {% icon param-file %} *"Previous Maker annotation"*: `evidences` (output of the previous second **Maker** {% icon tool %} run)
>        - *"Re-use ESTs"*: `Yes`
>        - *"Re-use alternate organism ESTs"*: `Yes`
>        - *"Re-use protein alignments"*: `Yes`
>        - *"Re-use repeats"*: `Yes`
>    - In *"EST evidences (for best results provide at least one of these)"*:
>        - *"Infer gene predictions directly from all ESTs"*: `No`
>    - In *"Protein evidences (for best results provide at least one of these)"*:
>        - *"Infer gene predictions directly from all protein alignments"*: `No`
>    - In *"Ab-initio gene prediction"*:
>        - {% icon param-file %} *"SNAP model"*: `snap model` (output of **Train SNAP** {% icon tool %})
>        - *"Prediction with Augustus"*: `Run Augustus with a custom prediction model`
>        - {% icon param-file %} *"Augustus model"*: `augustus model` (output of **Train Augustus** {% icon tool %})
>    - In *"Repeat masking"*:
>        - *"Repeat library source"*: `Disable repeat masking (not recommended)`
>
{: .hands_on}


## Annotation statistics

Do we get a better result from Maker after this third run? Let's run the same tools as after the first and second run, and compare the results.

> ### {% icon hands_on %} Hands-on: Get annotation statistics
>
> 1. **Genome annotation statistics** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotation to analyse"*: `final annotation` (output of **Maker** {% icon tool %} third run)
>    - *"Reference genome"*: `Use a genome from history`
>        - {% icon param-file %} *"Corresponding genome sequence"*: select the genome sequence from your history
>
>
{: .hands_on}

## Busco

> ### {% icon hands_on %} Hands-on: Extract transcript sequences
>
> 1. **gffread** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input GFF3 or GTF feature file"*: `final annotation` (output of **Maker** {% icon tool %} third run)
>    - *"Reference Genome"*: `select the genome sequence from your history`
>        - *"Select fasta outputs"*:
>           - `fasta file with spliced exons for each GFF transcript (-w exons.fa)`
>           - `fasta file with spliced CDS for each GFF transcript (-x cds.fa)`
>           - `protein fasta file with the translation of CDS for each record (-y pep.fa)`
>    - *"full GFF attribute preservation (all attributes are shown)"*: `Yes`
>    - *"decode url encoded characters within attributes"*: `Yes`
>    - *"warn about duplicate transcript IDs and other potential problems with the given GFF/GTF records"*: `Yes`
>
{: .hands_on}

Now run BUSCO with the predicted transcript sequences:

> ### {% icon hands_on %} Hands-on: Run BUSCO
>
> 1. **Busco** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `exons` (output of **gffread** {% icon tool %})
>    - *"Mode"*: `Transcriptome`
>    - *"Lineage"*: `fungi_odb9`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. How do you explain the evolution of the gene number?
> 2. How do the third annotation compare to the previous ones? Did the second ab-initio predictors training improve the results?
>
> > ### {% icon solution %} Solution
> >
> > 1. We now get 4,767 genes instead of 5,036 after the second round. But you'll notice that the average gene length have increased from 1,562 to 1,658. It means that in this third round, Maker was able to merge some genes that were considered separate beforehand.
> > 2. The annotation looks better: BUSCO found 256 complete single-copy instead of 253, and genes are longer, and have more exons. It means Maker was able to predict genes with more complex structure.
> >
> {: .solution}
>
{: .question}

Usually no more than two rounds of training is needed to get the best results from the ab-initio predictors. You can try to retrain Augustus and SNAP, but you will probably notice very few changes. We will keep the final annotation we obtained for the rest of this tutorial.

# Improving gene naming

If you look at the content of the `final annotation` dataset, you will notice that the gene names are long, complicated, and not very readable. That's because Maker assign them automatic names based on the way it computed each gene model. We are now going to automatically assign more readable names.

> ### {% icon hands_on %} Hands-on: Change gene names
>
> 1. **Map annotation ids** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Maker annotation where to change ids"*: `final annotation` (output of **Maker** {% icon tool %})
>    - *"Prefix for ids"*: `"TEST_"`
>    - *"Justify numeric ids to this length"*: `"6"`
>
>    > ### {% icon comment %} Comment
>    >
>    > Genes will be renamed to look like: `TEST_001234`. You can replace `TEST_` by anything you like, usually an uppercase short prefix.
>    {: .comment}
>
{: .hands_on}

Look at the generated dataset, it should be much more readable, and ready for an official release.


# Visualising the results

With Galaxy, you can visualize the annotation you have generated using JBrowse. This allows you to navigate along the chromosomes of the genome and see the structure of each predicted gene.

> ### {% icon hands_on %} Hands-on: Visualize annotations in JBrowse
>
> 1. **JBrowse** {% icon tool %} with the following parameters:
>    - *"Reference genome to display"*: `Use a genome from history`
>        - {% icon param-file %} *"Select the reference genome"*: select the genome sequence from your history
>    - *"JBrowse-in-Galaxy Action"*: `New JBrowse Instance`
>    - In *"Track Group"*:
>        - Click on *"Insert Track Group"*:
>        - In *"1: Track Group"*:
>            - *"Track Category"*: `Maker annotation`
>            - In *"Annotation Track"*:
>                - Click on *"Insert Annotation Track"*:
>                - In *"1: Annotation Track"*:
>                    - *"Track Type"*: `GFF/GFF3/BED/GBK Features`
>                    - {% icon param-files %} *"GFF/GFF3/BED Track Data"*: select the final annotation of each **Maker** {% icon tool %} run
>                    - *"This is match/match_part data"*: `No`
>
{: .hands_on}

Enable the three different tracks on the left side of JBrowse, then navigate along the genome and compare the three different annotations. You should see how Maker progressively produced more complex gene models.

> ### {% icon question %} Questions
>
> Navigate to the sequence named `NC_003421.2` (NCBI identifier for Chromosome III), between positions `41800` and `45200` (or between `143850` and `148763` if you used the full genome sequence sequence).
> ![JBrowse navigation](../../images/jbrowse_navigate.png "Navigating to the given sequence and positions.")
> 1. How did the annotation improved in this region after each Maker round?
>
> > ### {% icon solution %} Solution
> >
> > 1. At the end of the first round, no gene model was predicted by Maker in this region, because not enough EST or protein could be aligned there.
> > After the second round, Maker was able to predict a first gene model with no intron. Notice the name of the model beginning with `snap_masked`: it means that Maker used mainly a gene prediction from SNAP to construct this gene model.
> > After the third round, there are two gene models, more complex as they contain introns. Training Augustus and SNAP allowed to refine the gene structures and to detect a new gene.
> >
> {: .solution}
>
{: .question}

## More visualisation

You might want to understand how a specific gene model was predicted by Maker. You can easily visualise the evidences used by Maker (EST ailgnements, protein aligments, ab-initio predictions, ...) by using JBrowse too.

> ### {% icon hands_on %} Hands-on: Visualize evidences in JBrowse
>
> 1. **JBrowse** {% icon tool %} with the following parameters:
>    - *"Reference genome to display"*: `Use a genome from history`
>        - {% icon param-file %} *"Select the reference genome"*: select the genome sequence from your history
>    - *"JBrowse-in-Galaxy Action"*: `Update existing JBrowse Instance`
>    - *"Previous JBrowse Instance"*: select the result from the previous **JBrowse** {% icon tool %} run
>    - In *"Track Group"*:
>        - Click on *"Insert Track Group"*:
>        - In *"1: Track Group"*:
>            - *"Track Category"*: `Maker evidences`
>            - In *"Annotation Track"*:
>                - Click on *"Insert Annotation Track"*:
>                - In *"1: Annotation Track"*:
>                    - *"Track Type"*: `GFF/GFF3/BED/GBK Features`
>                    - {% icon param-files %} *"GFF/GFF3/BED Track Data"*: select the "evidences" output of each **Maker** {% icon tool %} run
>                    - *"This is match/match_part data"*: `Yes`
>                        - *"Match Part Feature Type"*: Leave empty
>
{: .hands_on}

You will now see three new tracks displaying all the evidences used by Maker to generate consensus gene models.

# Conclusion
{:.no_toc}

Congratulations, you finished this tutorial! You learned how to annotate an eukaryotic genome using Maker, how to evaluate the quality of the annotation, and how to visualize it using the JBrowse genome browser.

# What's next?

After generating your annotation, you will probably want to automatically assign functional annotation to each predicted gene model. You can do it by using Blast, InterProScan, or Blast2GO for example.

An automatic annotation of an eukaryotic genome is rarely perfect. If you inspect some predicted genes, you will probably find some mistakes made by Maker, e.g. wrong exon/intron limits, splitted genes, or merged genes. Setting up a manual curation project using [Apollo](http://genomearchitect.org/) helps a lot to manually fix these errors.
