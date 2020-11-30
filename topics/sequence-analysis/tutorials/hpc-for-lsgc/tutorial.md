---
layout: tutorial_hands_on

title: "Large-scale genome comparison"

zenodo_link: ""

questions:
- How can we run pairwise genome comparisons using Galaxy?
- How can we extract further information from sequence comparisons in Galaxy?
objectives:
- Learn the basics of pairwise sequence comparison
- Learn how to run different tools in Galaxy to perform sequence comparison at fine and coarse-grained levels
- Learn how to post-process your sequence comparisons
time_estimation: 2H
key_points:
- We learnt that sequence comparison is a demanding problem and that there are several ways to approach it
- We learnt how to run sequence comparisons in Galaxy with different levels of precision
- Fine-grained and coarse-grained sequence comparison using GECKO and CHROMEISTER for smaller and larger sequences, respectively
- We learnt how to post-process our comparison by extracting alignments, performing realigments, etc.
requirements:
  -
    type: "internal"
    topic_name: introduction
    tutorials:
      - galaxy-intro-101

contributors:
- estebanpw

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

Sequence comparison is a core problem in bioinformatics. It is used widely in evolutionary studies, structural and functional analyses, assembly, metagenomics, etc. Despite its regular presence in everyday Life-sciences pipelines, it is still not a trivial step that can be overlooked. Therefore, understanding how sequence comparison works is key to developing efficient workflows that are central to so many other disciplines.

In the following tutorial, we will learn how to compare both small and large sequences using both seed-based alignment methods and alignment-free methods, and how to post-process our comparisons interactively to refine our results. Besides, we will also learn about the theoretical background of sequence comparison, including why some tools are suitable for some jobs and others are not.

This tutorial is divided into two large sections:

 - Fine-grained interactive sequence comparison: In this part of the tutorial we will use `GECKO` and `GECKO-MGV` to perform sequence alignment between small sequences. We will also identify, extract and re-align regions of interest.
 - Coarse-grained sequence comparison: In this part of the tutorial, we will tackle on how to compare massive sequences using `CHROMEISTER`, and alignment-free sequence comparison tool. We will generate visualization plots for the comparison of large plant genomes.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Fine-grained interactive sequence comparison

Imagine you are working on an evolutionary study regarding the species `mycoplasma hyopneumoniae`. In particular, you are interested in the strains `232` and `7422` and wish to compare their DNA sequence to know more about the evolutionary changes that took place between both. The workflow you will follow starts with (1) acquiring the data,(2)  getting it ready for Galaxy, (3) running the comparison and (4) inspecting and working with the resulting alignments. Let's go! 

## Preparing the data

First we will be uploading the data to Galaxy so that we can run our tools on it.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial and give it a descriptive name (e.g. "Sequence comparison hands-on"
>
>    {% include snippets/create_new_history.md %}
>    {% include snippets/rename_history.md %}
>
> 2. Import `mycoplasma-hyopneumoniae-232` and `mycoplasma-hyopneumoniae-7422` from [Zenodo](zenodoFolderLink). You can also download these two sequences from the NCBI from [here](https://www.ncbi.nlm.nih.gov/nuccore/NC_006360.1?report=fasta) and [here](https://www.ncbi.nlm.nih.gov/nuccore/NC_021831.1?report=fasta).
>
>    ```
>    https://zenodo.org/linkTo232
>    https://zenodo.org/linkTo7422
>    ```
>
>    {% include snippets/import_via_link.md %}
>
>    As default, Galaxy takes the link as name, so rename them.
>
> 3. Rename the files to `232.fasta` and `7422.fasta` and change the datatype to `fasta`.
>
>    {% include snippets/rename_dataset.md %}
>    {% include snippets/change_datatype.md %}
>
{: .hands_on}

If you were successful, both sequences should now be available as `.fasta` datasets in your history.


> ### {% icon question %} Questions
>
> 1. What do you think about the size of the sequences in regards to the difficulty of comparing them?
>
> > ### {% icon solution %} Solution
> > 1. It always depends on the focus of our study. For instance, if we were looking for optimal alignments, two 1 MB sequences are indeed large enough to make most approaches either fail or take a decent amount of time and resources. On the other hand, if we were looking for seed-based local alignments (e.g. `GECKO` ({% cite GECKO %}) or `BLAST` ({% cite BLAST %}) ), the comparison would require merely seconds (check the slides for more information).
> >
> {: .solution }
>
{: .question}

As we discussed in the previous section, running optimal aligning tools on such "small" data can be difficult (in fact, tools such as the well-known `EMBOSS needle` will require cuadratic amounts of time and memory, whereas tools such as `EMBOSS strecher` will require quadratic time ({% cite myers1988optimal %})). These limitations can make it impractical in many situations. Therefore, we will now learn how to overcome these limitations by employing seed-based methods, particularly `GECKO`.

## Running the comparison

We will now run a comparison between `mycoplasma hyopneumoniae 232` and `mycoplasma hyopneumoniae 7422` in Galaxy using `GECKO`.

> ### {% icon hands_on %} Hands-on: Comparing two mycoplasmas with GECKO
> 1. **GECKO** {% icon tool %} with the following parameters
>    - {% icon param-file %} *"Query sequence"*: `232.fasta`
>    - {% icon param-file %} *"Reference sequence"*: `7422.fasta`
>    - *"K-mer seed size"*: `16`
>    - *"Minimum length"*: `50`
>    - *"Minimum similarity"*: `60`
>    - *"Generate alignments file?"*: `Extract alignments (CSV and alignments file)`
> 2. Check out the files that have been generated, i.e. the `CSV` and the `Alignments` file. 
>    > ### {% icon question %} Questions
>    >
>    > 1. What information is provided in the `CSV` file?
>    > 2. And in the `Alignments` file?
>    > 2. What happens if we re-run the experiment with other parameters (e.g. change `Minimum length` to `5000` and `Minimum similarity` to `95`)?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. The `CSV` file contains a summary of the detected High-Scoring Segment Pairs (HSPs) or alignments. The file is divided in a few rows of metadata (e.g. containing the sequence names) and one row per alignment detected. Check out the `GECKO` help (bottom of the tool {% icon tool %} page to know what each column does!
>    > > 2. The `Alignments` file contains the actual regions of the query and reference sequence for each alignment detected. With this file, you can investigate individual alignments, find mutations and differences between the sequences, extract the aligned part, etc.
>    > > 3. Changing the parameters affect the number of alignments GECKO will detect. In fact, if we use parameters that are too restrictive, then we might not get any alignments at all! On the other hand, if we use parameters that are too permissive, then we might get a lot of noise in the output. Thus, it is very important to understand what the parameters do. Never leave your parameters default, always know what they do! Check out the help section to get information about the parameters.
>    >  {: .solution }
>    {: .question}
>
{: .hands_on}

## Interactive post-processing

todo

# Coarse-grained sequence comparison

todo

> ### {% icon question %} Questions
>
> 1. We already discussed the mycoplasma sequences in regards to their size. What do you think about the size of the two chromosomes that we will be comparing now?
>
> > ### {% icon solution %} Solution
> > 1. Chromosome-like sequences are arguably some of the largest DNA sequences you can find. Regarding the comparison, optimal chromosome comparison typically requires either large clusters to be run or special hardware accelerators (such as GPUs). On the other hand, for seed-based local alignment we can still use common approaches, but they will take a considerable amount of time. Finally, if we use alignment-free methods such as CHROMEISTER ({% cite perez2019ultra %}), we can perform a comparison in less than 5 minutes (check the slides for more information).
> >
> {: .solution }
>
{: .question}



# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
