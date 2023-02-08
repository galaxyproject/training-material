---
layout: tutorial_hands_on

title: Preparing genomic data for phylogeny reconstruction
zenodo_link: https://zenodo.org/record/6610704
questions:
- How do I find a set of common proteins (orthologs) across related species or strains?
- How do I organize a set of orthologs to infer evolutionary relations between species or strains (phylogenetic reconstruction)?
objectives:
- Mask repetitive elements from a genome
- Annotate (predict protein-coding genes) the genomes of the samples to compare
- Find a set of common proteins across the samples (orthologs)
- Align orthologs across samples
time_estimation: 3H
requirements:
  -
    type: "internal"
    topic_name: assembly
    tutorials:
      - general-introduction

key_points:
- You now are able to
- Predict proteins in a nucleotide sequence *de-novo* using **funannotate_predict**
- Find orthologs across different samples with **orthofinder**
- Align orthologs with **ClustalW** in preparation for phylogeny reconstruction <!-- link to phylogeny reconstruction training. -->
tags:
  - phylogeny
  - data handling
  - functional annotation
contributors:
- roncoronimiguel
- brigidagallone

---

A robust and well-resolved phylogenetic classification is essential to understand genetic relationships within and between species and the evolution of their phenotypic diversity. In the last decade the genomic revolution has represented a drastic change in the amount of data used for phylogenetic inference. The single-gene approach using universal phylogenetic markers for the different lineages across the tree of life, is now being replaced by the assembly of taxon-rich and genome-scale data matrices, the so called phylogenomic approach.

Molecular sequence data can be used to construct a phylogeny by comparing differences between nucleotide or amino acid sequences across species or strains, a technique called phylogenomics. {% cite Young2019 %} have written a comprehensive review on the topic of phylogenomics.

In this tutorial we prepare genetic sequence data for phylogenetic reconstruction, using sequences from chromosome 5 of five strains of the yeast *Saccharomyces cerevisiae*.

To prepare the data, we will:
1. Predict protein coding genes from the genome using Funannotate
2. Find the proteins that are present in more than one genome, called orthologs, using Proteinortho and extract a subset with orthologs that are present in all samples.
4. Align each set of orthologs using ClustalW.

The resulting dataset is ready to be used for phylogenetic reconstruction.


**If you are starting from sequence reads, please follow
[An Introduction to Genome Assembly]({% link topics/assembly/tutorials/general-introduction/tutorial.md %}), and the appropriate genome assembly training for your sequencing technology from GTN's [Assembly]({% link topics/assembly/index.md %}) section**

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get data
For this training we will use a subset of the genome (chromosome 5) from four strains of *S. cerevisiae*. The GenBank annotated sequenced were produced using 'funannotate predict annotation' (Galaxy Version 1.8.9+galaxy2) on the nucleotide sequences sequences.

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{page.topic_name}}`
>     -> `{{page.title}}`):
>
>    ```
>    https://zenodo.org/record/6610704/files/BK006939.2.genbank
>    https://zenodo.org/record/6610704/files/BK006939.2.genome.fasta
>    https://zenodo.org/record/6610704/files/CM000925.1.genbank
>    https://zenodo.org/record/6610704/files/CM000925.1.genome.fasta
>    https://zenodo.org/record/6610704/files/CM005043.2.genbank
>    https://zenodo.org/record/6610704/files/CM005043.2.genome.fasta
>    https://zenodo.org/record/6610704/files/CM005299.1.genbank
>    https://zenodo.org/record/6610704/files/CM005299.1.genome.fasta
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>
> 3. Optional: Rename each dataset to its accession number followed by '.genome.fasta' or '.genbank'.
> 4. Group the datasets into [collections]({% link topics/galaxy-interface/tutorials/collections/tutorial.md %}). These will ease data handling and help minimize the clutter in your history. Make a collection of nucleotide sequences and another of protein sequences.
>
>    {% snippet faqs/galaxy/collections_build_list.md %}
>
>
{: .hands_on}

# Protein coding gene prediction

## Mask repetitive sequences

Before we can annotate the genome, we will prepare the data by masking repetitive sequences in the genome.
Repeat-rich regions can interfere with genome annotation tools. In this step we find and soft-mask repetitive regions in the genome. The annotation tool can then take this information into account ([see RepeatMasker tutorial for more details]({% link topics/genome-annotation/tutorials/repeatmasker/tutorial.md %})).
We use **RepeatMasker**
{% cite RepeatMasker %}, a program that screens DNA sequences for interspersed repeats and low complexity DNA sequences.
This program has historically made use of [RepBase](https://www.girinst.org/repbase/update/index.html) ({%cite Kohany2006-ks %}), a service of the Genetic Information Research Institute, but this database in no longer open access. Instead, we will use [Dfam](https://www.dfam.org/home) ({%cite Storer2021 %}) an open collection of Transposable Element DNA sequence alignments,  HMMs derived from Repbase sequences and consensus sequences. For this reason, the annotation of repetitive sequences might be incomplete.

RepeatMasker will only accept compact fasta headers. Before we can mask repetitive regions with RepeatMasker we must trim the NCBI long header (`BK006939.2 TPA_inf: Saccharomyces cerevisiae S288C chromosome V, complete sequence`) to leave only the accession number (`>BK006939.2`) by using a regular expression.

> <hands-on-title> Mask repetitive sequences </hands-on-title>
>
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_line/1.1.2) %} with the following parameters:
>    - {% icon param-collection %} *"File to process"*: `output` (Input dataset collection)
>    - In *"Replacement"*:
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `(>[^ ]+).+`
>            - *"Replace with:"*: `\1`
>
>               > <comment-title></comment-title>
>               >
>               >`\1` replaces the text in the header with the first matched text in the header, the accession number in this case.
>               {: .comment}
>
> 2. {% tool [RepeatMasker](toolshed.g2.bx.psu.edu/repos/bgruening/repeat_masker/repeatmasker_wrapper/4.1.2-p1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Genomic DNA"*: `outfile` (output of **Replace Text** {% icon tool %})
>    - *"Repeat library source"*: `DFam (curated only, bundled with RepeatMasker)`
>        - *"Select species name from a list?"*: `No`
>            - *"Repeat source species"*: `"Saccharomyces cerevisiae"`
>
>    > <comment-title></comment-title>
>    >
>    > If you don't select the species from the list, you must provide a species name ("Saccharomyces cerevisiae" in this example). In principal, all unique clade names occurring in [NCBI taxonomy database](https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html) can be used for species. Capitalization is ignored, multiple words need to bound by apostrophes. Not all "common" English names occur in the taxonomy database. Using Latin names is always safest.
>    {: .comment}
>
>3. Inspect the 'RepeatMasker masked sequence on data' output file. Scroll down and you will find stretches of 'N' on the location of repetitive sequences. The file is now ready for annotation with Funannotate.
{: .hands_on}

## Annotate with Funannotate

We will predict protein-coding genes from genomic sequences using [Funannotate](https://funannotate.readthedocs.io/) ({% cite Young2019 %}), which collects evidence from different ab-initio gene predictors as well as from RNA-seq or ESTs data. Funannotate has been developed for Fungi but it works with any Eukaryotic genome. The output of Funannotate is a list of ORFs and their translation in GenBank format.

> <comment-title></comment-title>
>
>If you would like to learn about genome annotation in more depth, the GTN has a [section]({% link topics/genome-annotation/index.md %}) dedicated to training on genome annotation, including a hands-on tutorial on [Funannotate]({% link topics/genome-annotation/tutorials/funannotate/tutorial.md %}).
{: .comment}

> <warning-title> Slow Step Ahead! </warning-title>
> Even for a small dataset, Funannotate can take a very long time to run. You can skip this step and use the Genbank files downloaded from Zenodo for the following step. These were generated using Funannotate as described in the hands-on below.
{: .warning}

> <hands-on-title> Annotate genome </hands-on-title>
>
> 1. {% tool [Funannotate predict annotation](toolshed.g2.bx.psu.edu/repos/iuc/funannotate_predict/funannotate_predict/1.8.9+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Assembly to annotate"*: `output_masked_genome` (output of **RepeatMasker** {% icon tool %})
>    - In *"Organism"*:
>        - *"Name of the species to annotate"*: `Saccharomyces cerevisiae`
>        - *"Is it a fungus species?"*: `Yes`
>        - *"Ploidy of assembly"*: `1`
>    - In *"Evidences"*:
>        - *"Select protein evidences"*: `Use UniProtKb/SwissProt (from selected Funannotate database)`
>
>    > <tip-title></tip-title>
>    >
>    > If available, include mRNA and/or ESTs in this evidence section to increase sensitivity of predictions.
>    {: .tip}
>
>    - In *"Busco"*:
>        - *"BUSCO models to align"*: `saccharomycetes`
>        - *"Initial Augustus species training set for BUSCO alignment"*: `saccharomyces`
>    - In *"Augustus settings (advanced)"*:
>        - *"Minimum number of models to train Augustus"*: `15`
>
>
>
>    > <comment-title></comment-title>
>    >
>    > When annotating full genomes, increase the number of *'Minimum number of models to train Augustus'* to an appropriate value. For the small sample dataset used here, values larger than 15 will result in failure.
>    {: .comment}
> 2. Inspect the output GenBank file. The FEATURES section contains the genome annotation of protein predictions, their location and their translation. Each predicted protein is given a unique ID, which will become the FASTA header in the next step.
{: .hands_on}

## Extract ORFs into FASTA files

In this step, we extract the protein sequences from the GenBank files into multi-FASTA files. Additionally, we modify the headers to include the accession number of the sample. This creates an unique accession-proteinID header for each predicted protein, and will allow us to retrieve them after we combine all predicted proteins into one multi-FASTA file.

So far, we have kept sequence files in a collection. Now we will combine all predicted proteins from all samples into one big multi-fasta file. Later, we will retrieve ortholog sequences from this file.

> <hands-on-title> Extract ORFs </hands-on-title>
>
> 1. {% tool [Extract ORF](toolshed.g2.bx.psu.edu/repos/bgruening/glimmer_gbk_to_orf/glimmer_gbk_to_orf/3.02) %} with the following parameters:
>    - {% icon param-file %} *"gene bank file"*: `annot_gbk` (output of **Funannotate predict annotation** {% icon tool %} - or the genbank files downloaded from Zenodo)
>
> 2. {% tool [Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regex1/1.0.1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `aa_output` (output of **Extract ORF** {% icon tool %})
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `>([^ ]+).+`
>            - *"Replacement"*: `>#{input_name}_\1`
>
> 1. {% tool [Collapse Collection](toolshed.g2.bx.psu.edu/repos/nml/collapse_collections/collapse_dataset/4.2) %} with the following parameters:
>    - {% icon param-file %} *"Collection of files to collapse into single dataset"*: `out_file1` (output of **Regex Find And Replace** {% icon tool %})
>    - *"Prepend File name"*: `Yes`
>
{: .hands_on}


> <question-title></question-title>
>
> 1. Is the number of predicted ORFs the same across the samples?
>
> > <solution-title></solution-title>
> >
> > 1. No, the number of predicted ORFs ranges from 193 to 199.
> >
> {: .solution}
>
{: .question}

# Find orthologs

Orthologs are genes in different species evolved from a common ancestral gene by a speciation (lineage-splitting) event and contain the information needed for building phylogenies. The input for this step is the collection of multi-FASTA extracted from the GenBank file and modified to have unique headers.
The result file 'orthology-groups' contains one row per orthogroup and one column per sample.

Next, we select orthogroups where all the species are represented by only one protein, 1:1 single copy orthologs (a total of 4 proteins per orthogroup for this dataset).

> <hands-on-title> Find and filter orthologs </hands-on-title>
>
> 1. {% tool [Proteinortho](toolshed.g2.bx.psu.edu/repos/iuc/proteinortho/proteinortho/6.0.14+galaxy2.9.1) %} with the following parameters:
>    - {% icon param-file %} *"Select the input fasta files (>2)"*: `out_file1` (output of **Regex Find And Replace** {% icon tool %})
>    - *"Activate synteny feature (POFF)"*: `no`
>
> 2. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `proteinortho` (output of **Proteinortho** {% icon tool %})
>    - *"With following condition"*: `c1==4 and c2==4`
>
> 3. Inspect the 'orthology-groups' tabular file.
{: .hands_on}


> <question-title></question-title>
>
> 1. Inspect the the 'orthology-groups' tabular file.
>
>    1. Do any samples contain more than one gene for any given orthogroup?
>    2. What is the name given to these homologous genes within the same genome?
>
>    > <solution-title></solution-title>
>    >
>    > 1. Yes, *CM000925.1* contains two genes on the first orthogroup.
>    > 2. Paralogs. These are gene copies created by a duplication event within the same genome.
>    >
>    {: .solution}
>
> 2. Look at the filtered orthogroups.
>     1. How many orthogroups are represented once only in all four samples?
>
>     > <solution-title></solution-title>
>     >
>     > 1. 173
>     >
>     {: .solution}
>
{: .question}




## Quality control with **Busco**

Busco is a dataset of nearl-universal, single-copy orthologs.
Here we use it for quality control of the assembly and annotation produced above.
It outputs a proxy of completeness, duplication and fragmentation of the annotation (Busco can also be used to assess the completeness of a genome assembly).

> <hands-on-title> QC with Busco </hands-on-title>
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/4.1.4) %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `out_file1` (output of **Regex Find And Replace** {% icon tool %})
>    - *"Mode"*: `Proteome`
>    - *"Lineage"*: `Saccharomycetes`
>    - In *"Advanced Options"*:
>        - *"Augustus species model"*: `Use the default species for selected lineage`
>
>    > <tip-title></tip-title>
>    > Make sure you select the 'lineage' closest to the species(s) you are analyzing.
>    >
>    {: .tip}
>
{: .hands_on}


> <question-title></question-title>
>
> 1. The [Complete ('C') metric](https://busco.ezlab.org/busco_userguide.html#complete) stands for complete Busco genes identified. Look at the Busco output file 'Short summary' for sample BK006939.2. What is the 'C' number?
> 2. Why is it so low?
>
> > <solution-title></solution-title>
> >
> > 1. 4.0%
> > 2. Our dataset contains only chromosome 5 of the yeast genome.
> >
> {: .solution}
>
{: .question}

## Extract proteins with Proteinortho

Next we extract 1:1 single copy orthologs and generate one multi fasta file per ortholog.

> <hands-on-title> Extract protein sequences </hands-on-title>
>
> 1. {% tool [Proteinortho grab proteins](toolshed.g2.bx.psu.edu/repos/iuc/proteinortho_grab_proteins/proteinortho_grab_proteins/6.0.14+galaxy2.9.1) %} with the following parameters:
>    - {% icon param-file %} *"Select the input fasta files"*: `output` (output of **Collapse Collection** {% icon tool %})
>    - *"Query type"*: `orthology-groups output file`
>        - {% icon param-file %} *"A orthology-groups file"*: `out_file1` (output of **Filter** {% icon tool %})
>
>
{: .hands_on}

The output is a collection of multi-fasta ortholog files. All species are represented in each file and are ready to be aligned.

# Align ortholog sequences

Alignment of sequences allows cross-species (or other taxonomic level, strain, taxa) comparison. An alignment is a hypothesis of positional homology between bases/amino acids of different sequences. A correct alignment is crucial for phylogenetic inference. Here we use ClustaW, a well-established pairwise sequence aligner.

First we modify the headers of the multi-fasta file, such that only the sample name is retained. This is important for future conatenation of alignment into a supermatrix (see Maximum likelihood GTN training)

> <hands-on-title> Align orthologs with ClustalW </hands-on-title>
>
> 1. {% tool [Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regex1/1.0.1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `listproteinorthograbproteins` (output of **Proteinortho grab proteins** {% icon tool %})
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `(>[^_]+).+`
>            - *"Replacement"*: `\1`
>
> 2. {% tool [ClustalW](toolshed.g2.bx.psu.edu/repos/devteam/clustalw/clustalw/2.1) %} with the following parameters:
>    - {% icon param-file %} *"FASTA file"*: `out_file1` (output of **Regex Find And Replace** {% icon tool %})
>    - *"Data type"*: `Protein sequences`
>    - *"Output alignment format"*: `FASTA format`
>    - *"Output complete alignment (or specify part to output)"*: `Complete alignment`
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Open the ClustalW output 'queryOrthoGroup121.fasta' and its corresponding multifasta input. You can compare them side by side activating the 'Scratchbook' on the top panel. What is the main difference between the sequences in the unaligned multifasta (input) and the ClustalW output multifasta?
>
> > <solution-title></solution-title>
> >
> > 1. Two of the sequences from this orthogroup are truncated (early stop codon). The alignment program inserts '-' to represent indels in the alignment.
> >
> {: .solution}
>
{: .question}



# Conclusion
In this tutorial, you have prepared genome sequence data for phylogenetic analysis. First, you have extracted information from these in the form of predicted protein sequences. You then grouped the predicted proteins into orthogroups and aligned them. These aligned sequences can now be used for reconstructing a phylogeny and building a phylogenetic tree.

