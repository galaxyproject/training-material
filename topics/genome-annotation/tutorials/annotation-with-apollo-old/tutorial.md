---
layout: tutorial_hands_on

title: Genome annotation with Apollo
zenodo_link: https://doi.org/10.5281/zenodo.3270822
tags:
  - eukaryote
questions:
  - How to visualize your genome after automated annotations have been performed?
  - How to manually annotate genome after automated annotations have been performed?
  - How to evaluate and visualize annotated genomic features?
objectives:
  - Load genome into Galaxy
  - View annotations in JBrowse
  - Manually Annotate genome with Apollo
  - Export genomes
time_estimation: 2h
key_points:
  - Apollo allows a group to manually annotate a eukaryotic genome.
  - Use Apollo allow to inspect the quality of an annotation.
  - Use Apollo to edit annotations within your group.
  - Export manual annotations as GFF3.
contributors:
  - nathandunn
---


# Introduction
{:.no_toc}

After automatically editing annotations using [Prokker](https://training.galaxyproject.org/training-material/topics/genome-annotation/tutorials/annotation-with-prokka/tutorial.html) or [Maker](https://training.galaxyproject.org/training-material/topics/genome-annotation/tutorials/annotation-with-maker/tutorial.html), its important to visualize and then manually refine any additional data. 

This process is most often done as part of a group.  

This demo is inspired by the [Apollo User's Guide](http://genomearchitect.github.io/users-guide/), which provides additional guidance. 


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
>    {% include snippets/create_new_history.md %}
>
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
- `S_pombe_chrIII.fasta` contains only the third chromosome from the full genome


> ### {% icon hands_on %} Hands-on: Choose your Genome
>
> 1. You need to choose between `S_pombe_chrIII.fasta` and `S_pombe_genome.fasta`:
>
>    - If you have time: use the full genome (`S_pombe_genome.fasta`), it will take more computing time, but the results will be closer to real-life data.
>    - If you want to get results faster: use the chromosome III (`S_pombe_chrIII.fasta`).
>
> 2. Rename the file you will use to `genome.fasta`. E.g. if you are using `S_pombe_chrIII.fasta`, rename it to `genome.fa`
>
>    {% include snippets/rename_dataset.md name="genome.fa" %}
>
{: .hands_on}

The two other datasets (`augustus_training_1.tar.gz` an `augustus_training_2.tar.gz`) will be used later in the tutorial.


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

You might want to understand how a specific gene model was predicted by Maker. You can easily visualise the evidences used by Maker (EST alignments, protein alignments, ab-initio predictions, ...) by using JBrowse too.

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

## Create an Organism from Apollo sequences 

(can we use the maker data or should we use the other data library?)

## Create refinements

View data at a particular position

### Create structural edits

Drag exons, create isoforms, etc. create introns, split transcripts.

### Create structural variations

View data in many ways

### Edit functional data. 


### Edit names, etc. 



# Conclusion
{:.no_toc}

Congratulations, you finished this tutorial! You learned how to manually refine predicted eukaryotic genomes using Apollo and export them to other forms.

# What's next?

After generating your refined genome, you'll want to merge it back into the official gene sets.  

If a de novo set, you can export it as GFF3 and load it into a tool like [Tripal](http://tripal.info) to provide visualization.

