---
layout: tutorial_hands_on
topic_name: genome-annotation
tutorial_name: annotation-with-maker
---

# Introduction
{:.no_toc}

<!-- This is a comment. -->

General introduction about the topic and then an introduction of the
tutorial (the questions and the objectives). It is nice also to have a
scheme to sum up the pipeline used during the tutorial. The idea is to
give to trainees insight into the content of the tutorial and the (theoretical
and technical) key concepts they will learn.

**Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data upload

Give some background about what the trainees will be doing in the section.

Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Import the following files from [Zenodo]() or from a data
>    library named `TODO` if available (ask your instructor)
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    > ### {% icon tip %} Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field
>    > * Press **Start**
>    >
>    > By default, Galaxy uses the url as the name, so please rename them to something more pleasing.
>    {: .tip}
>
>    > ### {% icon tip %} Tip: Importing data from a data library
>    >
>    > * Go into "Shared data" (top panel) then "Data libraries"
>    > * Click on "Training data" and then "Genome annotation with Maker"
>    > * Select interesting file
>    > * Click on "Import selected datasets into history"
>    > * Import in a new history
>    {: .tip}
>
{: .hands_on}


# Genome quality evaluation

The quality of a genome annotation is highly dependent on the quality of the genome sequences. It is impossible to obtain a good quality annotation with a poorly assembled genome sequence. Annotation tools will have trouble finding genes if the genome sequence is highly fragmented, if it contains chimeric sequences, or if there are a lot of sequencing errors.

Before running the full annotation process, it is a good idea first to evaluate the quality of the sequence. It will give you a good idea of what you can expect from it at the end of the annotation.

First have a look at the basic statistics: scaffold number compared to expected chromosome numbers, N50, ...

> ### {% icon hands_on %} Hands-on: Get genome sequence statistics
>
> 1. **Fasta Statistics** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: select the genome sequence from your history
>
{: .hands_on}

These statistics are useful to detect obvious problems in the genome assembly, but it gives no information about the quality of the sequence content. We want to know if the genome sequence contains all the genes we expect to find in the considered species, and if their sequence are correct.

[BUSCO](http://busco.ezlab.org/) (Benchmarking Universal Single-Copy Orthologs) is a tool allowing to answer this question: by comparing genomes from various more or less related species, the authors were able to determine sets of genes that are present in single copy in (almost) all the species of a clade (Bacteria, Fungi, Plants, Insects, Mammalians, ...). Most of these genes are essential for the organism to live, and are expected to be found in any newly sequences genome from the corresponding clade. Using this data, BUSCO is able to evaluate the number of genes found in a genome sequence or a set of (predicted) transcript or protein sequences.

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

> ### {% icon question %} Questions
>
> 1. Are the BUSCO results good for our genome?
>
> > ### {% icon solution %} Solution
> >
> > 1. Most of the BUSCO genes are found as complete single copy, which means that our genome have a good quality as it contains most of the expected content.
> >
> {: .solution}
>
{: .question}


# First Maker annotation round

## Maker

Maker blablabla

For this first round, we configure Maker to construct gene models only by aligning ESTs and proteins to the genome. This will produce a first draft annotation that we will improve in the steps.

> ### {% icon hands_on %} Hands-on: First draft annotation with Maker
>
> 1. **Maker** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select the genome sequence from your history
>    - *"Organism type"*: `Eukaryotic`
>    - In *"EST evidences (for best results provide at least one of these)"*:
>        - *"Infer gene predictions directly from all ESTs"*: `Yes`
>    - In *"Protein evidences (for best results provide at least one of these)"*:
>        - *"Infer gene predictions directly from all protein alignments"*: `Yes`
>    - In *"Repeat masking"*:
>        - *"Enable repeat masking with RepeatMasker"*: `No`
>
>    ***TODO***: *Check repeatmasking option*
>
{: .hands_on}

Maker produces three gff datasets:

- the final annotation: the gene models produced by Maker
- the evidences: the alignement of all the data Maker used to construct the final annotation (ESTs and proteins that we used)
- a gff containing both the final annotation and the evidences

## Annotation statistics

We need now to evaluate this first annotation produced by Maker.

First, use the `Genome annotation statistics` that will compute some general statistics on the annotation.

> ### {% icon hands_on %} Hands-on: Get annotation statistics
>
> 1. **Genome annotation statistics** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotation to analyse"*: `final annotation` (output of **Maker** {% icon tool %})
>
>
{: .hands_on}


## Busco

Just as we did for the genome at the beginning, we can use BUSCO to check the quality of this first Maker annotation. Instead of looking for known genes in the genome sequence, BUSCO will look for known transcript sequences of the genes predicted by Maker.

First we need to compute all the transcript sequences from the Maker annotation.

> ### {% icon hands_on %} Hands-on: Task description
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

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Busco** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `output_exons` (output of **gffread** {% icon tool %})
>    - *"Mode"*: `Transcriptome`
>    - *"Lineage"*: `fungi_odb9`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. How many genes where predicted by Maker?
> 2. What is the medium length of these genes?
> 2. How do the BUSCO statistics compare to the ones at the genome level?
>
> > ### {% icon solution %} Solution
> >
> > 1. I don't know, I should check. TODO
> > 2. I don't know, I should check. TODO
> > 3. I don't know, I should check. TODO
> >
> {: .solution}
>
{: .question}

Let's see now how this first draft can be improved.

# Ab-initio predictors first training

Maker can use several ab-initio predictors to annotate a genome. "Ab-initio" means that these predictors are able to predict the structure of genes in a genome based only on its sequence and on a species-specific statistical model. They don't use any evidence (e.g. EST or proteins) to predict genes, but they need to be trained with a set of already predicted genes.

Today we will use two of the most widely used ab-initio predictors SNAP and Augustus. Before using it within Maker, we need to train them with the first annotation draft we produced in the previous steps. We know the quality of this draft is not perfect, but only the best scoring genes (ie the ones having the strongest evidences) will be retained to train the predictors.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Train SNAP** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select the genome sequence from your history
>    - {% icon param-file %} *"Maker annotation to use for training"*: `final annotation` (output of **Maker** {% icon tool %})
>    - *"Number of gene model to use for training"*: `"1000"` ***TODO should we reduce it with our small genome?***
>
>    > ### {% icon comment %} Comment
>    >
>    > The parameter "Number of gene model to use for training" set to `"1000"` means that SNAP will be trained on the 1000 best-scoring genes.
>    {: .comment}
>
{: .hands_on}


Augustus is trained in a very similar way.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Train Augustus** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select the genome sequence from your history
>    - {% icon param-file %} *"Annotation to use for training"*: `final annotation` (output of **Maker** {% icon tool %})
>
{: .hands_on}

The Augustus training usually take around 2 hours to complete, to continue this tutorial without waiting for the result, you can get it from Zenodo. ***TODO***

Both SNAP and Augustus produce a statistical model representing the observed general structure of genes in the analysed genome. Maker will use these models to create a new annotation for our genome.

# Second Maker annotation round

## Maker

We need now to run a new round of Maker. As the evidences were already aligned on the genome on the first run, we can reuse these alignments.
This time, enable ab-initio gene prediction, and input the output of **Train SNAP** {% icon tool %} and **Train Augustus** {% icon tool %} tools.
We also disable infering gene predictions directly from all ESTs and proteins: now we want Maker to infer gene predictions by reconciliating evidence alignments *and* ab-initio gene predictions.

> ### {% icon hands_on %} Hands-on: Second draft annotation with Maker
>
> 1. **Maker** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select the genome sequence from your history
>    - *"Organism type"*: `Eukaryotic`
>    - *"Re-annotate using an existing Maker annotation"*: `Yes`
>        - *"Previous Maker annotation"*: `evidences` (output of the previous **Maker** {% icon tool %} run)
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
>        - *"Augustus model"*: `augustus model` (output of **Train Augustus** {% icon tool %})
>    - In *"Repeat masking"*:
>        - *"Enable repeat masking with RepeatMasker"*: `No`
>
{: .hands_on}


## Annotation statistics

Do we get a better result from Maker after this second run? Let's run the same tools as after the first run, and compare the results.

> ### {% icon hands_on %} Hands-on: Get annotation statistics
>
> 1. **Genome annotation statistics** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotation to analyse"*: `final annotation` (output of **Maker** {% icon tool %} second run)
>
>
{: .hands_on}

## Busco

> ### {% icon hands_on %} Hands-on: Task description
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

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Busco** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `output_exons` (output of **gffread** {% icon tool %})
>    - *"Mode"*: `Transcriptome`
>    - *"Lineage"*: `fungi_odb9`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. How do the second annotation compare to the previous one? Did the ab-initio predictors training improve the results?
>
> > ### {% icon solution %} Solution
> >
> > 1. Oh yes, better, because of x and y ***TODO***
> >
> {: .solution}
>
{: .question}


# Ab-initio predictors second training

To get better results, we are going to perform a second training of SNAP and Augustus, and then run Maker for a third (final) time. Experience shows that no more than two rounds of training is needed to get the best results from the ab-initio predictors. ***TODO: ref + put it at the end?***

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Train SNAP** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select the genome sequence from your history
>    - {% icon param-file %} *"Maker annotation to use for training"*: `final annotation` (output of **Maker** {% icon tool %}, second run)
>    - *"Number of gene model to use for training"*: `"1000"`
>
{: .hands_on}

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Train Augustus** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select the genome sequence from your history
>    - {% icon param-file %} *"Annotation to use for training"*: `final annotation` (output of **Maker** {% icon tool %}, second run)
>
{: .hands_on}

# Third (last) Maker annotation round

## Maker

Let's run the final round of Maker, in the same way as we did for the second run.

> ### {% icon hands_on %} Hands-on: Final annotation with Maker
>
> 1. **Maker** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select the genome sequence from your history
>    - *"Organism type"*: `Eukaryotic`
>    - *"Re-annotate using an existing Maker annotation"*: `Yes`
>        - *"Previous Maker annotation"*: `evidences` (output of the previous second **Maker** {% icon tool %} run)
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
>        - *"Augustus model"*: `augustus model` (output of **Train Augustus** {% icon tool %})
>    - In *"Repeat masking"*:
>        - *"Enable repeat masking with RepeatMasker"*: `No`
>
{: .hands_on}


## Annotation statistics

Do we get a better result from Maker after this third run? Let's run the same tools as after the first and second run, and compare the results.

> ### {% icon hands_on %} Hands-on: Get annotation statistics
>
> 1. **Genome annotation statistics** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotation to analyse"*: `final annotation` (output of **Maker** {% icon tool %} second run)
>
>
{: .hands_on}

## Busco

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **gffread** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input GFF3 or GTF feature file"*: `final annotation` (output of **Maker** {% icon tool %} second run)
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

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Busco** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `output_exons` (output of **gffread** {% icon tool %})
>    - *"Mode"*: `Transcriptome`
>    - *"Lineage"*: `fungi_odb9`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. How do the second annotation compare to the previous ones? Did the second ab-initio predictors training improve the results?
>
> > ### {% icon solution %} Solution
> >
> > 1. Oh yes, a little better, because of x and y ***TODO***
> >
> {: .solution}
>
{: .question}


## Improving gene naming

We want to keep the results of the third Maker run. If you look at the content of the `final annotation` dataset, you will notice that the gene names are very long, complicated, and not very readable. That's because Maker assign them automatic names based on the way it computed each gene model. We are now going to automatically assign more readable names.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Map annotation ids** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Maker annotation where to change ids"*: `final annotation` (output of **Maker** {% icon tool %})
>    - *"Prefix for ids"*: `"TEST"`
>    - *"Justify numeric ids to this length"*: `"6"`
>
>    > ### {% icon comment %} Comment
>    >
>    > You can replace `TEST` by anything you like, usually an uppercase short prefix.
>    {: .comment}
>
{: .hands_on}

Look at the generated dataset, it should be much more readable.


## Visualising the results

With Galaxy, you can visualize the annotation you have generated using JBrowse. This allows you to navigate along the chromosomes of the genome and see the structure of each predicted gene.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **JBrowse** {% icon tool %} with the following parameters:
>    - *"Reference genome to display"*: `Use a genome from history`
>    - *"Select the reference genome"*: select the genome sequence from your history
>    - *"Produce Standalone Instance"*: `Yes`
>    - *"Genetic Code"*: ``
>    - Insert a *"Track Group"*:
>        - In *"1: Track Group"*:
>            - *"Track Category"*: `Maker annotation`
>            - Insert an *"Annotation Track"*:
>                - In *"1: Annotation Track"*:
>                    - *"Track Type"*: `GFF/GFF3/BED/GBK Features`
>                    - *"GFF/GFF3/BED Track Data"*: select the final annotation of each Maker run
>                    - *"This is match/match_part data"*: `No`
>                    - *"Index this track"*: `No`
>    - *""*: `""`
>
{: .hands_on}

***TODO Go to region XXX for example***

> ### {% icon question %} Questions
>
> 1. How do the three annotations compare in the genome browser? Is it consistent with the annotation statistics and BUSCO results you obtained before?
>
> > ### {% icon solution %} Solution
> >
> > 1. TODO
> >
> {: .solution}
>
{: .question}


# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
