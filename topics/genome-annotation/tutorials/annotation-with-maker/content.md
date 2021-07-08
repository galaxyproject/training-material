{% if include.short %}
  {% assign genome_name = "S_pombe_chrIII.fasta" %}
{% else %}
  {% assign genome_name = "genome.fasta"%}
{% endif %}

# Introduction
{:.no_toc}

Genome annotation of eukaryotes is a little more complicated than for prokaryotes: eukaryotic genomes are usually larger than prokaryotes, with more genes. The sequences determining the beginning and the end of a gene are generally less conserved than the prokaryotic ones. Many genes also contain introns, and the limits of these introns (acceptor and donor sites) are not highly conserved.

In this tutorial we will use a software tool called Maker {% cite Campbell2014 %} to annotate the genome sequence of a small eukaryote: *Schizosaccharomyces pombe* (a yeast).

Maker is able to annotate both prokaryotes and eukaryotes. It works by aligning as many evidences as possible along the genome sequence, and then reconciliating all these signals to determine probable gene structures.

The evidences can be transcript or protein sequences from the same (or closely related) organism. These sequences can come from public databases (like NR or GenBank) or from your own experimental data (transcriptome assembly from an RNASeq experiment for example). Maker is also able to take into account repeated elements.

Maker uses ab-initio predictors (like [Augustus](http://bioinf.uni-greifswald.de/augustus/) or [SNAP](https://github.com/KorfLab/SNAP)) to improve its predictions: these software tools are able to make gene structure predictions by analysing only the genome sequence with a statistical model.

In this tutorial you will learn how to perform a genome annotation, and how to evaluate its quality. {% unless include.short %}You will see how training ab-initio predictors is an important step to produce good results. {% endunless %}Finally, you will learn how to use the [JBrowse](http://jbrowse.org/) genome browser to visualise the results.

More information about Maker can be found [here](http://www.yandell-lab.org/software/maker.html).

This tutorial was inspired by the MAKER Tutorial for [WGS Assembly and Annotation Winter School 2018](http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018), don't hesitate to consult it for more information on Maker, and on how to run it with command line.


> ### {% icon comment %} Note: Two versions of this tutorial
>
> Because this tutorial consists of many steps, we have made two versions of it, one long and one short.
>
> {% if include.short %}
> This is the **shortened version**. We will skip the training of ab-initio predictors
> and use pre-trained data instead. We will also annotate only the third chromosome of the genome. If you would like
> to learn how to perform the training steps, please see the [longer version of tutorial]({% link topics/genome-annotation/tutorials/annotation-with-maker/tutorial.md %})
> {% else %}
> This is the **extended version**. We will perform the complete training of ab-initio predictors and discuss the results in detail.
> If you would like to run through the tutorial a bit quicker and focus on the main
> analysis steps, please see the [shorter version of this tutorial]({% link topics/genome-annotation/tutorials/annotation-with-maker-short/tutorial.md %})
> {% endif %}
{: .comment}

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
- A set of transcripts or [EST sequences](https://en.wikipedia.org/wiki/Expressed_sequence_tag), preferably from the same organism.
- A set of protein sequences, usually from closely related species or from a curated sequence database like UniProt/SwissProt.

 Maker will align the transcript and protein sequences on the genome sequence to determine gene positions.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create and name a new history for this tutorial.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the following files from [Zenodo](https://doi.org/10.5281/zenodo.4406623) or from the shared data library
>
>    {% if include.short %}
>    ```
>    https://zenodo.org/api/files/647ad552-19a8-46d9-aad8-f81f56860582/S_pombe_chrIII.fasta
>    https://zenodo.org/api/files/647ad552-19a8-46d9-aad8-f81f56860582/S_pombe_trinity_assembly.fasta
>    https://zenodo.org/api/files/647ad552-19a8-46d9-aad8-f81f56860582/Swissprot_no_S_pombe.fasta
>    https://zenodo.org/api/files/647ad552-19a8-46d9-aad8-f81f56860582/augustus_training_2.tar.gz
>    https://zenodo.org/api/files/647ad552-19a8-46d9-aad8-f81f56860582/snap_training_2.snaphmm
>    ```
>    {% else %}
>    ```
>    https://zenodo.org/api/files/647ad552-19a8-46d9-aad8-f81f56860582/S_pombe_chrIII.fasta
>    https://zenodo.org/api/files/647ad552-19a8-46d9-aad8-f81f56860582/S_pombe_genome.fasta
>    https://zenodo.org/api/files/647ad552-19a8-46d9-aad8-f81f56860582/S_pombe_trinity_assembly.fasta
>    https://zenodo.org/api/files/647ad552-19a8-46d9-aad8-f81f56860582/Swissprot_no_S_pombe.fasta
>    https://zenodo.org/api/files/647ad552-19a8-46d9-aad8-f81f56860582/augustus_training_1.tar.gz
>    https://zenodo.org/api/files/647ad552-19a8-46d9-aad8-f81f56860582/augustus_training_2.tar.gz
>    ```
>    {% endif %}
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype for {% unless include.short %}`augustus_training_1.tar.gz` and{% endunless %} `augustus_training_2.tar.gz` is set to `augustus`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="augustus" %}
>
{: .hands_on}

You have the following main datasets:

- `S_pombe_trinity_assembly.fasta` contains EST sequences from *S. pombe*, assembled from RNASeq data with Trinity
- `Swissprot_no_S_pombe.fasta` contains a subset of the SwissProt protein sequence database (public sequences from *S. pombe* were removed to stay as close as possible to real-life analysis)
- `S_pombe_chrIII.fasta` contains only the third chromosome from the full genome of *S. pombe*
{% unless include.short %}- `S_pombe_genome.fasta` contains the full genome sequence of *S. pombe*

> ### {% icon hands_on %} Hands-on: Choose your Genome
>
> 1. You need to choose between `S_pombe_chrIII.fasta` and `S_pombe_genome.fasta`:
>
>    - If you have time: use the full genome (`S_pombe_genome.fasta`), it will take more computing time, but the results will be closer to real-life data.
>    - If you want to get results faster: use the chromosome III (`S_pombe_chrIII.fasta`).
>
> 2. Rename the file you will use to `genome.fasta`. E.g. if you are using `S_pombe_chrIII.fasta`, rename it to `genome.fa`
>
>    {% snippet faqs/galaxy/datasets_rename.md name="genome.fa" %}
>
{: .hands_on}
{% endunless %}

The other datasets will be used later in the tutorial.

# Genome quality evaluation

The quality of a genome annotation is highly dependent on the quality of the genome sequences. It is impossible to obtain a good quality annotation with a poorly assembled genome sequence. Annotation tools will have trouble finding genes if the genome sequence is highly fragmented, if it contains chimeric sequences, or if there are a lot of sequencing errors.

Before running the full annotation process, you need first to evaluate the quality of the sequence. It will give you a good idea of what you can expect from it at the end of the annotation.

> ### {% icon hands_on %} Hands-on: Get genome sequence statistics
>
> 1. {% tool [Fasta Statistics](toolshed.g2.bx.psu.edu/repos/iuc/fasta_stats/fasta-stats/1.0.1) %} with the following parameters:
>    - {% icon param-file %} *"fasta or multifasta file"*: select `{{ genome_name }}` from your history
>
{: .hands_on}

Have a look at the statistics:

- `num_seq`: the number of contigs (or scaffold or chromosomes), compare it to expected chromosome numbers
- `len_min`, `len_max`, `len_N50`, `len_mean`, `len_median`: the distribution of contig sizes
- `num_bp_not_N`: the number of bases that are not N, it should be as close as possible to the total number of bases (`num_bp`)

These statistics are useful to detect obvious problems in the genome assembly, but it gives no information about the quality of the sequence content. We want to evaluate if the genome sequence contains all the genes we expect to find in the considered species, and if their sequence are correct.

> ### {% icon comment %} Comment
>
> {% if include.short %}
> Keep in mind that we are running this tutorial only on the chromosome III instead of the whole genome.
> {% else %}
> If you chose to use only the chromosome III sequence (`S_pombe_chrIII.fasta`), the statistics will be different. You will only see 1 contig.
> {% endif %}
{: .comment}

[BUSCO](http://busco.ezlab.org/) (Benchmarking Universal Single-Copy Orthologs) is a tool allowing to answer this question: by comparing genomes from various more or less related species, the authors determined sets of ortholog genes that are present in single copy in (almost) all the species of a clade (Bacteria, Fungi, Plants, Insects, Mammalians, ...). Most of these genes are essential for the organism to live, and are expected to be found in any newly sequenced genome from the corresponding clade. Using this data, BUSCO is able to evaluate the proportion of these essential genes (also named BUSCOs) found in a genome sequence or a set of (predicted) transcript or protein sequences. This is a good evaluation of the "completeness" of the genome or annotation.

We will first run this tool on the genome sequence to evaluate its quality.


> ### {% icon hands_on %} Hands-on: Run Busco on the genome
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/4.1.4) %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: select `{{ genome_name }}` from your history
>    - *"Mode"*: `Genome`
>    - *"Lineage"*: `Fungi`
>
>    > ### {% icon comment %} Comment
>    >
>    > We select `Fungi` as we will annotate the genome of *Schizosaccharomyces pombe* which belongs to the Fungi kingdom. It is usually better to select the most specific lineage for the species you study. Large lineages (like Metazoa) will consist of fewer genes, but with a strong support. More specific lineages (like Hymenoptera) will have more genes, but with a weaker support (has they are found in fewer genomes).
>    {: .comment}
>
{: .hands_on}

BUSCO produces three output datasets

- A short summary: summarizes the results of BUSCO (see below)
- A full table: lists all the BUSCOs that were searched for, with the corresponding status (was it found in the genome? how many times? where?)
- A table of missing BUSCOs: this is the list of all genes that were not found in the genome

![BUSCO genome summary](../../images/busco_genome_summary.png "Example of BUSCO short summary (not from this tutorial). Here BUSCO searched for 290 genes in a genome, and found 282 of them complete (269 in single-copy, and 13 duplicated). 5 where found but are fragmented, and 3 were not found at all.")

> ### {% icon question %} Questions
>
> Do you think the genome quality is good enough for performing the annotation?
>
> > ### {% icon solution %} Solution
> >
> > The genome consists of the expected number of chromosome sequences ({% if include.short %}1{% else %}4, or 1 if you chose to annotate chromosome III only{% endif %}), with very few N, which is the ideal case.
> > {% if include.short %}As we only analysed chromosome III, many BUSCO genes are missing, but still ~100 are found as complete single copy, and very few are found fragmented, which means that our genome have a good quality, at least on this single chromosome. That's a very good material to perform an annotation.{% else %}
> > If you used the full genome, most of the BUSCO genes are found as complete single copy, and very few are fragmented, which means that our genome have a good quality as it contains most of the expected content. That's a very good material to perform an annotation. If you only analysed chromosome III, many BUSCO genes are missing, but still ~100 are found as complete single copy, and very few are found fragmented, which means that our genome have a good quality, at least on this single chromosome.{% endif %}
> >
> {: .solution}
>
{: .question}

> ### {% icon comment %} Comment
>
> {% if include.short %}
> Keep in mind that we are running this tutorial only on the chromosome III instead of the whole genome. The BUSCO result will also show a lot of missing genes: it is expected as all the BUSCO genes that are not on the chromosome III cannot be found by the tool.
> {% else %}
> If you chose to use only the chromosome III sequence (`S_pombe_chrIII.fasta`), the statistics will be different. The genome size will be lower, with only 1 chromosome. The BUSCO result will also show a lot of missing genes: it is expected as all the BUSCO genes that are not on the chromosome III cannot be found by the tool. Keep it in mind when comparing these results with the other BUSCO results later.
> {% endif %}
{: .comment}

{% if include.short %}
# Maker

Let's run Maker to predict gene models! Maker will use align [ESTs](https://en.wikipedia.org/wiki/Expressed_sequence_tag) and proteins to the genome, and it will run ab initio predictors (SNAP and Augustus) using pre-trained models for this organism (have a look at the [longer version of tutorial]({% link topics/genome-annotation/tutorials/annotation-with-maker/tutorial.md %}) to understand how they were trained).

> ### {% icon hands_on %} Hands-on: Annotation with Maker
>
> 1. {% tool [Maker](toolshed.g2.bx.psu.edu/repos/iuc/maker/maker/2.31.11) %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select `{{ genome_name }}` from your history
>    - *"Organism type"*: `Eukaryotic`
>    - *"Re-annotate using an existing Maker annotation"*: `No`
>    - In *"EST evidences (for best results provide at least one of these)"*:
>        - {% icon param-file %} *"ESTs or assembled cDNA"*: `S_pombe_trinity_assembly.fasta`
>    - In *"Protein evidences (for best results provide at least one of these)"*:
>        - {% icon param-file %} *"Protein sequences"*: `Swissprot_no_S_pombe.fasta`
>    - In *"Ab-initio gene prediction"*:
>        - *"SNAP model"*: `snap_training_2.snaphmm`
>        - *"Prediction with Augustus"*: `Run Augustus with a custom prediction model`
>            - {% icon param-file %} *"Augustus model"*: `augustus_training_2.tar.gz`
>    - In *"Repeat masking"*:
>        - *"Repeat library source"*: `Disable repeat masking (not recommended)`
>
>    > ### {% icon comment %} Comment
>    >
>    > For this tutorial repeat masking is disabled, which is not the recommended setting. When doing a real-life annotation, you should either use [Dfam](https://www.dfam.org) or provide your own repeats library.
>    {: .comment}
>
{: .hands_on}

{% else %}
# First Maker annotation round

## Maker

For this first round, we configure Maker to construct gene models only by aligning [ESTs](https://en.wikipedia.org/wiki/Expressed_sequence_tag) and proteins to the genome. This will produce a first draft annotation that we will improve in the next steps.

> ### {% icon hands_on %} Hands-on: First draft annotation with Maker
>
> 1. {% tool [Maker](toolshed.g2.bx.psu.edu/repos/iuc/maker/maker/2.31.11) %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select `{{ genome_name }}` from your history
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
>    > For this tutorial repeat masking is disabled, which is not the recommended setting. When doing a real-life annotation, you should either use [Dfam](https://www.dfam.org) or provide your own repeats library.
>    {: .comment}
>
{: .hands_on}

{% endif %}

Maker produces three GFF3 datasets:

- The final annotation: the final consensus gene models produced by Maker
- The evidences: the alignments of all the data Maker used to construct the final annotation (ESTs and proteins that we used)
- A GFF3 file containing both the final annotation and the evidences

{% if include.short %}
# Annotation statistics
{% else %}
## Annotation statistics
{% endif %}

We need now to evaluate this {% unless include.short %}first draft {% endunless %}annotation produced by Maker.

First, use the `Genome annotation statistics` that will compute some general statistics on the annotation.

> ### {% icon hands_on %} Hands-on: Get annotation statistics
>
> 1. {% tool [Genome annotation statistics](toolshed.g2.bx.psu.edu/repos/iuc/jcvi_gff_stats/jcvi_gff_stats/0.8.4) %} with the following parameters:
>    - {% icon param-file %} *"Annotation to analyse"*: `final annotation` (output of {% tool [Maker](toolshed.g2.bx.psu.edu/repos/iuc/maker/maker/2.31.11) %})
>    - *"Reference genome"*: `Use a genome from history`
>        - {% icon param-file %} *"Corresponding genome sequence"*: select `{{ genome_name }}` from your history
>
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. How many genes where predicted by Maker?
> 2. What is the mean gene locus size of these genes?
>
> > ### {% icon solution %} Solution
> >
> > {% if include.short %}
> > 1. 864 genes
> > 2. 1793 bp
> > {% else %}
> > 1. 712 genes (503 for the chromosome III only)
> > 2. 1570 bp (1688bp for the chromosome III only)
> > {% endif %}
> >
> {: .solution}
>
{: .question}


{% if include.short %}
# Busco
{% else %}
## Busco
{% endif %}

Just as we did for the genome at the beginning, we can use BUSCO to check the quality of this {% unless include.short %}first {% endunless %}Maker annotation. Instead of looking for known genes in the genome sequence, BUSCO will inspect the transcript sequences of the genes predicted by Maker. This will allow us to see if Maker was able to properly identify the set of genes that Busco found in the genome sequence at the beginning of this tutorial.

First we need to compute all the transcript sequences from the Maker annotation, using {% tool [GFFread](toolshed.g2.bx.psu.edu/repos/devteam/gffread/gffread/2.2.1.1) %}. This tool will compute the sequence of each transcript that was predicted by {% tool [Maker](toolshed.g2.bx.psu.edu/repos/iuc/maker/maker/2.31.11) %} and write them all in a FASTA file.

> ### {% icon hands_on %} Hands-on: Extract transcript sequences
>
> 1. {% tool [GFFread](toolshed.g2.bx.psu.edu/repos/devteam/gffread/gffread/2.2.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Input GFF3 or GTF feature file"*: `final annotation` (output of {% tool [Maker](toolshed.g2.bx.psu.edu/repos/iuc/maker/maker/2.31.11) %})
>    - *"Reference Genome"*: `select `{{ genome_name }}` from your history`
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
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/4.1.4) %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `exons` (output of {% tool [GFFread](toolshed.g2.bx.psu.edu/repos/devteam/gffread/gffread/2.2.1.1) %})
>    - *"Mode"*: `Transcriptome`
>    - *"Lineage"*: `Fungi`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> How do the BUSCO statistics compare to the ones at the genome level?
>
> > ### {% icon solution %} Solution
> >
> > {% if include.short %}
> > 128 complete single-copy, 0 duplicated, 10 fragmented, 620 missing. This is in fact better than what BUSCO found in the genome sequence. That means the quality of this annotation is very good (by default BUSCO in genome mode can miss some genes, the advanced options can improve this at the cost of computing time). (Results can be very slightly different in your own history, it's normal).
> > {% else %}
> > Around 100 complete single-copy, and 650 missing. As the quality of this first draft is yey not very good, you should see better results after next rounds of Maker.
> > {% endif %}
> >
> {: .solution}
>
{: .question}


{% unless include.short %}
The statistics are not really satisfactory at this stage, but it's normal: Maker only used the EST and protein evidences to guess the gene positions. Let's see now how this first draft can be improved.

# Ab-initio predictors first training

Maker can use several ab-initio predictors to annotate a genome. "Ab-initio" means that these predictors are able to predict the structure of genes in a genome based only on its sequence and on a species-specific statistical model. They don't use any evidence (e.g. EST or proteins) to predict genes, but they need to be trained with a set of already predicted genes.

Maker is able to use the EST and protein evidences, and to combine them with the result of several ab-initio predictors to predict consensus gene models. It allows to detect genes in regions where no EST or protein align, and also to refine gene structures in regions where there are EST and/or protein evidences and ab-initio predictions.

We will use two of the most widely used ab-initio predictors: SNAP and Augustus. Before using it within Maker, we need to train them with the first annotation draft we produced in the previous steps. We know the quality of this draft is not perfect, but only the best scoring genes (ie the ones having the strongest evidences) will be retained to train the predictors.

> ### {% icon hands_on %} Hands-on: Train SNAP
>
> 1. {% tool [Train SNAP](toolshed.g2.bx.psu.edu/repos/iuc/snap_training/snap_training/2013_11_29+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select `{{ genome_name }}` from your history
>    - {% icon param-file %} *"Maker annotation to use for training"*: `final annotation` (output of {% tool [Maker](toolshed.g2.bx.psu.edu/repos/iuc/maker/maker/2.31.11) %})
>
{: .hands_on}


Augustus is trained in a very similar way.

> ### {% icon hands_on %} Hands-on: Train Augustus
>
> 1. {% tool [Train Augustus](toolshed.g2.bx.psu.edu/repos/bgruening/augustus_training/augustus_training/3.3.3) %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select `{{ genome_name }}` from your history
>    - {% icon param-file %} *"Annotation to use for training"*: `final annotation` (output of {% tool [Maker](toolshed.g2.bx.psu.edu/repos/iuc/maker/maker/2.31.11) %})
>
{: .hands_on}

Both SNAP and Augustus produce a statistical model representing the observed general structure of genes in the analysed genome. Maker will use these models to create a new annotation for our genome.

The Augustus training usually takes around 2 hours to complete, to continue this tutorial without waiting for the result, you can use the file `augustus_training_1.tar.gz` imported from Zenodo.

# Second Maker annotation round

## Maker

We need now to run a new round of Maker. As the evidences were already aligned on the genome on the first run, we can reuse these alignments as is.
But this time, enable ab-initio gene prediction, and input the output of {% tool [Train SNAP](toolshed.g2.bx.psu.edu/repos/iuc/snap_training/snap_training/2013_11_29+galaxy1) %} and {% tool [Train Augustus](toolshed.g2.bx.psu.edu/repos/bgruening/augustus_training/augustus_training/3.3.3) %} tools.
We also disable inferring gene predictions directly from all ESTs and proteins: now we want Maker to infer gene predictions by reconciliating evidence alignments *and* ab-initio gene predictions.

> ### {% icon hands_on %} Hands-on: Second draft annotation with Maker
>
> 1. {% tool [Maker](toolshed.g2.bx.psu.edu/repos/iuc/maker/maker/2.31.11) %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select `{{ genome_name }}` from your history
>    - *"Organism type"*: `Eukaryotic`
>    - *"Re-annotate using an existing Maker annotation"*: `Yes`
>        - {% icon param-file %} *"Previous Maker annotation"*: `evidences` (output of the previous {% tool [Maker](toolshed.g2.bx.psu.edu/repos/iuc/maker/maker/2.31.11) %} run)
>        - *"Re-use ESTs"*: `Yes`
>        - *"Re-use alternate organism ESTs"*: `Yes`
>        - *"Re-use protein alignments"*: `Yes`
>        - *"Re-use repeats"*: `Yes`
>    - In *"EST evidences (for best results provide at least one of these)"*:
>        - *"Infer gene predictions directly from all ESTs"*: `No`
>    - In *"Protein evidences (for best results provide at least one of these)"*:
>        - *"Infer gene predictions directly from all protein alignments"*: `No`
>    - In *"Ab-initio gene prediction"*:
>        - *"SNAP model"*: `snap model` (output of {% tool [Train SNAP](toolshed.g2.bx.psu.edu/repos/iuc/snap_training/snap_training/2013_11_29+galaxy1) %})
>        - *"Prediction with Augustus"*: `Run Augustus with a custom prediction model`
>            - {% icon param-file %} *"Augustus model"*: `augustus model` (output of {% tool [Train Augustus](toolshed.g2.bx.psu.edu/repos/bgruening/augustus_training/augustus_training/3.3.3) %})
>    - In *"Repeat masking"*:
>        - *"Repeat library source"*: `Disable repeat masking (not recommended)`
>
{: .hands_on}


## Annotation statistics

Do we get a better result from Maker after this second run? Let's run the same tools as after the first run, and compare the results.

> ### {% icon hands_on %} Hands-on: Get annotation statistics
>
> 1. {% tool [Genome annotation statistics](toolshed.g2.bx.psu.edu/repos/iuc/jcvi_gff_stats/jcvi_gff_stats/0.8.4) %} with the following parameters:
>    - {% icon param-file %} *"Annotation to analyse"*: `final annotation` (output of {% tool [Maker](toolshed.g2.bx.psu.edu/repos/iuc/maker/maker/2.31.11) %} second run)
>    - *"Reference genome"*: `Use a genome from history`
>        - {% icon param-file %} *"Corresponding genome sequence"*: select `{{ genome_name }}` from your history
>
>
{: .hands_on}

## Busco

> ### {% icon hands_on %} Hands-on: Extract transcript sequences
>
> 1. {% tool [GFFread](toolshed.g2.bx.psu.edu/repos/devteam/gffread/gffread/2.2.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Input GFF3 or GTF feature file"*: `final annotation` (output of {% tool [Maker](toolshed.g2.bx.psu.edu/repos/iuc/maker/maker/2.31.11) %} second run)
>    - *"Reference Genome"*: `select `{{ genome_name }}` from your history`
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
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/4.1.4) %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `exons` (output of {% tool [GFFread](toolshed.g2.bx.psu.edu/repos/devteam/gffread/gffread/2.2.1.1) %})
>    - *"Mode"*: `Transcriptome`
>    - *"Lineage"*: `Fungi`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. How does the second annotation compare to the previous one? Did the ab-initio predictors training improve the results?
> 2. How do you explain these changes?
>
> > ### {% icon solution %} Solution
> >
> > 1. The annotation looks much better: more BUSCO complete single-copy, and more genes
> > 2. Using ab-initio predictors allowed to find much more genes in regions where EST or protein alignments were not sufficient to predict genes.
> >
> {: .solution}
>
{: .question}


# Ab-initio predictors second training

To get better results, we are going to perform a second training of SNAP and Augustus, and then run Maker for a third (final) time.

> ### {% icon hands_on %} Hands-on: Train SNAP and Augustus
>
> 1. {% tool [Train SNAP](toolshed.g2.bx.psu.edu/repos/iuc/snap_training/snap_training/2013_11_29+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select `{{ genome_name }}` from your history
>    - {% icon param-file %} *"Maker annotation to use for training"*: `final annotation` (output of {% tool [Maker](toolshed.g2.bx.psu.edu/repos/iuc/maker/maker/2.31.11) %}, second run)
>    - *"Number of gene model to use for training"*: `"1000"`
>
> 2. {% tool [Train Augustus](toolshed.g2.bx.psu.edu/repos/bgruening/augustus_training/augustus_training/3.3.3) %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select `{{ genome_name }}` from your history
>    - {% icon param-file %} *"Annotation to use for training"*: `final annotation` (output of {% tool [Maker](toolshed.g2.bx.psu.edu/repos/iuc/maker/maker/2.31.11) %}, second run)
>
{: .hands_on}

The Augustus training usually take around 2 hours to complete, to continue this tutorial without waiting for the result, you can use the file `augustus_training_2.tar.gz` imported from Zenodo.

# Third (last) Maker annotation round

## Maker

Let's run the final round of Maker, in the same way as we did for the second run.

> ### {% icon hands_on %} Hands-on: Final annotation with Maker
>
> 1. {% tool [Maker](toolshed.g2.bx.psu.edu/repos/iuc/maker/maker/2.31.11) %} with the following parameters:
>    - {% icon param-file %} *"Genome to annotate"*: select `{{ genome_name }}` from your history
>    - *"Organism type"*: `Eukaryotic`
>    - *"Re-annotate using an existing Maker annotation"*: `Yes`
>        - {% icon param-file %} *"Previous Maker annotation"*: `evidences` (output of the previous second {% tool [Maker](toolshed.g2.bx.psu.edu/repos/iuc/maker/maker/2.31.11) %} run)
>        - *"Re-use ESTs"*: `Yes`
>        - *"Re-use alternate organism ESTs"*: `Yes`
>        - *"Re-use protein alignments"*: `Yes`
>        - *"Re-use repeats"*: `Yes`
>    - In *"EST evidences (for best results provide at least one of these)"*:
>        - *"Infer gene predictions directly from all ESTs"*: `No`
>    - In *"Protein evidences (for best results provide at least one of these)"*:
>        - *"Infer gene predictions directly from all protein alignments"*: `No`
>    - In *"Ab-initio gene prediction"*:
>        - {% icon param-file %} *"SNAP model"*: `snap model` (output of {% tool [Train SNAP](toolshed.g2.bx.psu.edu/repos/iuc/snap_training/snap_training/2013_11_29+galaxy1) %})
>        - *"Prediction with Augustus"*: `Run Augustus with a custom prediction model`
>        - {% icon param-file %} *"Augustus model"*: `augustus model` (output of {% tool [Train Augustus](toolshed.g2.bx.psu.edu/repos/bgruening/augustus_training/augustus_training/3.3.3) %})
>    - In *"Repeat masking"*:
>        - *"Repeat library source"*: `Disable repeat masking (not recommended)`
>
{: .hands_on}


## Annotation statistics

Do we get a better result from Maker after this third run? Let's run the same tools as after the first and second run, and compare the results.

> ### {% icon hands_on %} Hands-on: Get annotation statistics
>
> 1. {% tool [Genome annotation statistics](toolshed.g2.bx.psu.edu/repos/iuc/jcvi_gff_stats/jcvi_gff_stats/0.8.4) %} with the following parameters:
>    - {% icon param-file %} *"Annotation to analyse"*: `final annotation` (output of {% tool [Maker](toolshed.g2.bx.psu.edu/repos/iuc/maker/maker/2.31.11) %} third run)
>    - *"Reference genome"*: `Use a genome from history`
>        - {% icon param-file %} *"Corresponding genome sequence"*: select `{{ genome_name }}` from your history
>
>
{: .hands_on}

## Busco

> ### {% icon hands_on %} Hands-on: Extract transcript sequences
>
> 1. {% tool [GFFread](toolshed.g2.bx.psu.edu/repos/devteam/gffread/gffread/2.2.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Input GFF3 or GTF feature file"*: `final annotation` (output of {% tool [Maker](toolshed.g2.bx.psu.edu/repos/iuc/maker/maker/2.31.11) %} third run)
>    - *"Reference Genome"*: `select `{{ genome_name }}` from your history`
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
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/4.1.4) %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `exons` (output of {% tool [GFFread](toolshed.g2.bx.psu.edu/repos/devteam/gffread/gffread/2.2.1.1) %})
>    - *"Mode"*: `Transcriptome`
>    - *"Lineage"*: `Fungi`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> How do the third annotation compare to the previous ones? Did the second ab-initio predictors training improve the results?
>
> > ### {% icon solution %} Solution
> >
> > Depending on wether you annotated the full genome or only chromosome III, you should see nearly the same, or even less genes than in the previous Maker round. But you'll notice that the number of multi-exon genes have increased. It means that in this third round, Maker was able to predict more complex genes, for example merging some genes that were considered separate beforehand.
> >
> {: .solution}
>
{: .question}

Usually no more than two rounds of training is needed to get the best results from the ab-initio predictors. You can try to retrain Augustus and SNAP, but you will probably notice very few changes. We will keep the final annotation we obtained for the rest of this tutorial.

{% endunless %}

# Improving gene naming

If you look at the content of the `final annotation` dataset, you will notice that the gene names are long, complicated, and not very readable. That's because Maker assign them automatic names based on the way it computed each gene model. We are now going to automatically assign more readable names.

> ### {% icon hands_on %} Hands-on: Change gene names
>
> 1. {% tool [Map annotation ids](toolshed.g2.bx.psu.edu/repos/iuc/maker_map_ids/maker_map_ids/2.31.11) %} with the following parameters:
>    - {% icon param-file %} *"Maker annotation where to change ids"*: `final annotation` (output of {% tool [Maker](toolshed.g2.bx.psu.edu/repos/iuc/maker/maker/2.31.11) %})
>    - *"Prefix for ids"*: `TEST_`
>    - *"Justify numeric ids to this length"*: `6`
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
> 1. {% tool [JBrowse](toolshed.g2.bx.psu.edu/repos/iuc/jbrowse/jbrowse/1.16.10+galaxy0) %} with the following parameters:
>    - *"Reference genome to display"*: `Use a genome from history`
>        - {% icon param-file %} *"Select the reference genome"*: select `{{ genome_name }}` from your history
>    - *"JBrowse-in-Galaxy Action"*: `New JBrowse Instance`
>    - In *"Track Group"*:
>        - Click on *"Insert Track Group"*:
>        - In *"1: Track Group"*:
>            - *"Track Category"*: `Maker annotation`
>            - In *"Annotation Track"*:
>                - Click on *"Insert Annotation Track"*:
>                - In *"1: Annotation Track"*:
>                    - *"Track Type"*: `GFF/GFF3/BED Features`
> {% if include.short %}
>                    - {% icon param-files %} *"GFF/GFF3/BED Track Data"*: select the output of {% tool [Map annotation ids](toolshed.g2.bx.psu.edu/repos/iuc/maker_map_ids/maker_map_ids/2.31.11) %}
> {% else %}
>                    - {% icon param-files %} *"GFF/GFF3/BED Track Data"*: select the final annotation of each {% tool [Maker](toolshed.g2.bx.psu.edu/repos/iuc/maker/maker/2.31.11) %} run
> {% endif %}
>
{: .hands_on}

{% if include.short %}
Enable the track on the left side of JBrowse, then navigate along the genome and look at the genes that were predicted by Maker.
{% else %}
Enable the three different tracks on the left side of JBrowse, then navigate along the genome and compare the three different annotations. You should see how Maker progressively produced more complex gene models.

> ### {% icon question %} Questions
>
> Navigate to the position `NC_003421.2:143850..148763` (meaning: on the sequence named `NC_003421.2` (NCBI identifier for Chromosome III), between positions `143850` and `148763`).
> ![JBrowse navigation](../../images/jbrowse_navigate.png "Navigating to the given sequence and positions.")
> 1. How did the annotation improved in this region after each Maker round?
>
> > ### {% icon solution %} Solution
> >
> > 1. At the end of the first round, a first short gene model was predicted by Maker in this region.
> > After the second round, Maker was able to predict a second gene model in this region. Notice the name of the model beginning with `snap_masked`: it means that Maker used mainly a gene prediction from SNAP to construct this gene model.
> > After the third round, the two genes were merged into a single one. Training Augustus and SNAP allowed to refine the gene structures and to refine the gene structures found in this region.
> >
> {: .solution}
>
{: .question}
{% endif %}

{% unless include.short %}
## More visualisation

You might want to understand how a specific gene model was predicted by Maker. You can easily visualise the evidences used by Maker (EST alignments, protein alignments, ab-initio predictions, ...) by using JBrowse too.

> ### {% icon hands_on %} Hands-on: Visualize evidences in JBrowse
>
> 1. {% tool [JBrowse](toolshed.g2.bx.psu.edu/repos/iuc/jbrowse/jbrowse/1.16.10+galaxy0) %} with the following parameters:
>    - *"Reference genome to display"*: `Use a genome from history`
>        - {% icon param-file %} *"Select the reference genome"*: select `{{ genome_name }}` from your history
>    - *"JBrowse-in-Galaxy Action"*: `Update existing JBrowse Instance`
>    - *"Previous JBrowse Instance"*: select the result from the previous {% tool [JBrowse](toolshed.g2.bx.psu.edu/repos/iuc/jbrowse/jbrowse/1.16.10+galaxy0) %} run
>    - In *"Track Group"*:
>        - Click on *"Insert Track Group"*:
>        - In *"1: Track Group"*:
>            - *"Track Category"*: `Maker evidences`
>            - In *"Annotation Track"*:
>                - Click on *"Insert Annotation Track"*:
>                - In *"1: Annotation Track"*:
>                    - *"Track Type"*: `GFF/GFF3/BED Features`
>                    - {% icon param-files %} *"GFF/GFF3/BED Track Data"*: select the "evidences" output of each {% tool [Maker](toolshed.g2.bx.psu.edu/repos/iuc/maker/maker/2.31.11) %} run
>                    - *"This is match/match_part data"*: `Yes`
>
{: .hands_on}

You will now see new tracks displaying all the evidences used by Maker to generate consensus gene models.
{% endunless %}

# Conclusion
{:.no_toc}

Congratulations, you finished this tutorial! You learned how to annotate an eukaryotic genome using Maker, how to evaluate the quality of the annotation, and how to visualize it using the JBrowse genome browser.

# What's next?

After generating your annotation, you will probably want to automatically assign functional annotation to each predicted gene model. You can do it by using Blast, InterProScan, or Blast2GO for example.

An automatic annotation of an eukaryotic genome is rarely perfect. If you inspect some predicted genes, you will probably find some mistakes made by Maker, e.g. wrong exon/intron limits, splitted genes, or merged genes. Setting up a manual curation project using [Apollo](http://genomearchitect.org/) helps a lot to manually fix these errors. Check out the [Apollo tutorial]({% link topics/genome-annotation/tutorials/apollo/tutorial.md %}) for more details.
