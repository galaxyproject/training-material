---
layout: tutorial_hands_on

title: "De Bruijn Graph Assembly"
zenodo_link: "https://doi.org/10.5281/zenodo.582600"
questions:
  - "What are the factors that affect genome assembly?"
  - "How does Genome assembly work?"
objectives:
  - "Perform an optimised Velvet assembly with the Velvet Optimiser"
  - "Compare this assembly with those we did in the basic tutorial"
  - "Perform an assembly using the SPAdes assembler."
time_estimation: "2h"
level: Introductory
key_points:
  - "We learned about how the choice of k-mer size will affect assembly outcomes"
  - "We learned about the strategies that assemblers use to make reference genomes"
  - "We performed a number of assemblies with Velvet and SPAdes."
  - "You should use SPAdes or another more modern assembler than Velvet for actual assemblies now."
contributors:
  - slugger70
  - hexylena
  - shiltemann
---

# Optimised de Bruijn Graph assemblies using the Velvet Optimiser and SPAdes

In this activity, we will perform *de novo* assemblies of a short read set using the Velvet Optimiser and the SPAdes assemblers. We are using the Velvet Optimiser for illustrative purposes. For real assembly work, a more suitable assembler should be chosen - such as SPAdes.

The Velvet Optimiser is a script written by Simon Gladman to optimise the k-mer size and coverage cutoff parameters for Velvet. More information can be found [here](https://github.com/slugger70/VelvetOptimiser)

SPAdes is a de novo genome assembler written by Pavel Pevzner's group in St. Petersburg. More details on it can be found [here](http://cab.spbu.ru/software/spades/)



> <agenda-title></agenda-title>
>
> In this tutorial, we will deal with:
>
> 1. [Get the data](#get-the-data)
> 2. [Assemble with the Velvet Optimiser](#assembly-with-the-velvet-optimiser)
> 3. [Assemble with SPAdes](#assemble-with-spades)
{: .agenda}

# Get the data

We will be using the same data that we used in the introductory tutorial, so if you have already completed that and have the data, skip this section.

> <hands-on-title>Getting the data</hands-on-title>
>
> 1. Create and name a new history for this tutorial.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the sequence read raw data (\*.fastq) from [Zenodo](https://zenodo.org/record/582600)
>
>    ```
>    https://zenodo.org/record/582600/files/mutant_R1.fastq
>    https://zenodo.org/record/582600/files/mutant_R2.fastq
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Rename the files {% icon galaxy-pencil %}
>    - The name of the files are the full URL, let's make the names a little clearer
>    - Change the names to just the last part, `Mutant_R1.fastq`, `Mutant_R2.fastq`  respectively
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
>    > <question-title></question-title>
>    >
>    > 1. What are four key features of a FASTQ file?
>    > 2. What is the main difference between a FASTQ and a FASTA file?
>    {: .question}
>
>
{: .hands_on}

# Assembly with the Velvet Optimiser

We will perform an assembly with the Velvet Optimiser, which automatically runs and optimises the output of the Velvet assembler ({% cite Velvet2008 %}). It will automatically choose a suitable value for the k-mer size (**k**). It will then go on to optimise the coverage cutoff (**cov_cutoff**) which corrects for read errors. It will use the "*n50*" metric for optimising the k-mer size and the "*total number of bases in contigs*" for optimising the coverage cutoff.

> <hands-on-title>Assemble with the Velvet Optimiser</hands-on-title>
>
>  1. **Velvet Optimiser** {% icon tool %}: Optimise your assembly with the following parameters:
>    - *"Start k-mer size"*: `45`
>    - *"End k-mer size"*: `73`
>    - *"Input file type"*: `Fastq`
>    - *"Single or paired end reads"*: `Paired`
>    - {% icon param-file %} *"Select first set of reads"*: `mutant_R1.fastq`
>    - {% icon param-file %} *"Select second set of reads"*: `mutant_R2.fastq`
>
{: .hands_on}

Your history will now contain a number of new files:

* Velvet optimiser contigs
  * A fasta file of the final assembled contigs
* Velvet optimiser contig stats
  * A table of the lengths (in k-mer length) and coverages (k-mer coverages) for the final contigs.

Have a look at each file.


> <hands-on-title>Get contig statistics for Velvet Optimiser contigs</hands-on-title>
>
> 1. **Fasta Statistics** {% icon tool %}: Produce a summary of the velvet optimiser contigs:
>    - {% icon param-file %} *"fasta or multifasta file"*: Select your velvet optimiser contigs file
>
> 2. View the output
>
>    > <question-title></question-title>
>    >
>    > Compare the output we got here with the output of the simple assemblies obtained in the introductory tutorial.
>    > 1. What are the main differences between them?
>    > 2. Which has a higher "n50"? What does this mean?
>    {: .question}
>
{: .hands_on}

Tables of results from **(a)** Simple assembly and **(b)** optimised assembly.

**(a)** ![The results of the contigs from Simple assembly.](../../images/image12.png)

**(b)** ![The results of the contigs from Optimised assembly. In contrast to simple assembly produced much higher n_50, while num_seq is lower.](../../images/optstats.png)

> <details-title>Further reading on assembly with Velvet</details-title>
> - Heuristic Resolution of Repeats and Scaffolding in the Velvet Short-Read de Novo Assembler ({% cite Zerbino2009 %})
>
{: .details}

## Visualisation of the Assembly

Now that we've assembled the genomes, let's visualise this assembly using [Bandage](https://rrwick.github.io/Bandage/) ({% cite Wick2015 %}). This tool will let us better understand how the assembly graph really looks, and can give us a feeling for if the genome was well assembled or not.

Currently VelvetOptimiser does not include the LastGraph output, so we will manually run `velveth` and `velvetg` with the optimised parameters.

> <hands-on-title>Manually running velvetg/h</hands-on-title>
>
> 1. Locate the output called "VelvetOptimiser: Contigs" in your history
>
> 2. Click the (i) information icon
>
> 3. Check the tool `stderr` in the information page for the optimised k-mer value
{: .hands_on}

> <question-title></question-title>
> What was the optimal k-mer value? (referred to as *"hash"* in the stderr log)
> > <solution-title></solution-title>
> > 55
> {: .solution}
{: .question}

With this information in hand, let's run velvet:

> <hands-on-title>Manually running velvetg/h</hands-on-title>
>
> 1. **velveth** {% icon tool %}: Prepare a dataset for the Velvet velvetg Assembler
>    - *"Hash length"*: `55`
>    - *"Insert Input Files"*:
>      - 1: Input Files
>        - *"file format"*: `fastq`
>        - *"read type"*: `shortPaired reads`
>        - *"Dataset"*: `mutant_R1.fastq`
>    - *"Insert Input Files"*:
>      - 2: Input Files
>        - *"file format"*: `fastq`
>        - *"read type"*: `shortPaired reads`
>        - *"Dataset"*: `mutant_R2.fastq`
>
> 2. **velvetg** {% icon tool %}: Velvet sequence assembler for very short reads
>    - *"Velvet dataset"*: output from **velveth** {% icon tool %}
>    - *"Generate velvet LastGraph file"*: `Yes`
>    - *"Coverage cutoff"*: `Specify Cutoff Value`
>      - *"Remove nodes with coverage below"*: `1.44`
>    - *"Using Paired Reads"*: `Yes`
>
{: .hands_on}

The LastGraph contains a detailed representation of the De Bruijn graph, which can give us an idea how velvet has assembled the genome and potentially resolved any conflicts.

> <hands-on-title>Bandage</hands-on-title>
>
> 1. **Bandage Image** {% icon tool %}: visualize de novo assembly graphs
>    - *"Graphical Fragment Assembly"*: The "LastGraph" output of **velvetg** {% icon tool %}
>    - *"Produce jpg, png or svg file?"*: `.svg`
>
> 2. Execute
> 3. View the output file
{: .hands_on}

And now you should be able to see the graph that velvet produced:

![velvet graph](../../images/bandage-velvet.svg)

## Interpreting Bandage Graphs

k-mer size has a [significant effect](https://github.com/rrwick/Bandage/wiki/Effect-of-kmer-size) on the assembly. You can play around with various k-mers to see this effect in practice.

k-mer | graph
----- | -----
21    | [![21](../../images/bandage-velvet-21.svg)](../../images/bandage-velvet-21.svg)
33    | [![33](../../images/bandage-velvet-33.svg)](../../images/bandage-velvet-33.svg)
53    | [![53](../../images/bandage-velvet-53.svg)](../../images/bandage-velvet-53.svg)
77    | [![77](../../images/bandage-velvet-77.svg)](../../images/bandage-velvet-77.svg)

The next thing to be aware of is that there can be multiple valid interpretations of a graph, all equally valid in absence of other data. The following is taken verbatim [from Bandage's wiki](https://github.com/rrwick/Bandage/wiki/Simple-example):

> For a simple case, imagine a bacterial genome that contains a single repeated element in two separate places in the chromosome:
>
> ![Simple example 1](https://camo.githubusercontent.com/03628b49f50ccf7a9c565d7712bfc70c7764cbeb/687474703a2f2f72727769636b2e6769746875622e696f2f42616e646167652f696d616765732f77696b692f73696d706c655f6578616d706c655f312e706e67)
>
> A researcher (who does not yet know the structure of the genome) sequences it, and the resulting 100 bp reads are assembled with a de novo assembler:
>
> ![Simple example 2](https://camo.githubusercontent.com/a51f384b83fbb97590ce86b8ec14d4ebd1bb60d1/687474703a2f2f72727769636b2e6769746875622e696f2f42616e646167652f696d616765732f77696b692f73696d706c655f6578616d706c655f322e706e67)
>
> Because the repeated element is longer than the sequencing reads, the assembler was not able to reproduce the original genome as a single contig. Rather, three contigs are produced: one for the repeated sequence (even though it occurs twice) and one for each sequence between the repeated elements.
>
> Given only the contigs, the relationship between these sequences is not clear. However, the assembly graph contains additional information which is made apparent in Bandage:
>
> ![Simple example 3](https://camo.githubusercontent.com/406648509cf478ac0b2ab9f2447aec4e7575b7dd/687474703a2f2f72727769636b2e6769746875622e696f2f42616e646167652f696d616765732f77696b692f73696d706c655f6578616d706c655f332e706e67)
>
> There are two principal underlying sequences compatible with this graph: two separate circular sequences that share a region in common, or a single larger circular sequence with an element that occurs twice:
>
> ![Simple example 4](https://camo.githubusercontent.com/58d0aa7eff4cfd3d36c9210e9f6a2f0265396715/687474703a2f2f72727769636b2e6769746875622e696f2f42616e646167652f696d616765732f77696b692f73696d706c655f6578616d706c655f342e706e67)
>
> Additional knowledge, such as information on the approximate size of the bacterial chromosome, can help the researcher to rule out the first alternative. In this way, Bandage has assisted in turning a fragmented assembly of three contigs into a completed genome of one sequence.
{: .quote}

# Assemble with SPAdes

We will now perform an assembly with the much more modern SPAdes assembler ({% cite Bankevich2012 %}). It goes through a similar process to Velvet in the fact that it uses and simplifies de Bruijn graphs but it uses multiple values for k-mer size and combines the resultant graphs. This combination produces very good assemblies. When using SPAdes it is typical to choose at least 3 k-mer sizes. One low, one medium and one high. We will use 33, 55 and 91.

> <hands-on-title>Assemble with SPAdes</hands-on-title>
>
> 1. **SPAdes** {% icon tool %}: Assemble the reads:
>
>    - *"Run only assembly"*: `yes`
>    - *"K-mers to use separated by commas"*: `33,55,91` [note: no spaces!]
>    - *"Coverage cutoff"*: `auto`
>    - {% icon param-file %} *"Files -> forward reads"*: `mutant_R1.fastq`
>    - {% icon param-file %} *"Files -> reverse reads"*: `mutant_R2.fastq`
>    - *"Output final assembly graph with scaffolds?"*: `Yes`
>
{: .hands_on}

You will now have 5 new files in your history:

* two Fasta files, one for contigs and one for scaffolds
* two statistics files, one for contigs and one for scaffolds
* the SPAdes log file.

Examine each file, especially the stats files.

![Contig stats file with NODE_5 being the shortest contig with the highest coverage and NODE_1 being the opposite.](../../images/contig_stats.png)

> <question-title></question-title>
>
> 1. Why would one of the contigs have much higher coverage than the others?
> 2. What could this represent?
>
{: .question}


> <hands-on-title>Visualize assembly with Bandage</hands-on-title>
>
> 1. **Bandage** {% icon tool %} with the following parameters:
>    - *"Graphical Fragment Assembly"*: `assembly graph with scaffolds` output from **SPAdes** {% icon tool %}
>
> 2. Examine the output image {% icon galaxy-eye %}
>
{: .hands_on}

The visualized assembly should look something like this:

![bandage spades](../../images/bandage_spades.svg)


> <question-title></question-title>
>
> Which assembly looks better to you? Why?
>
{: .question}



> <hands-on-title>Get contig statistics for SPAdes contigs</hands-on-title>
>
> 1. **Fasta Statistics** {% icon tool %}: Produce a summary of the SPAdes contigs:
>    - {% icon param-file %} *"fasta or multifasta file"*: Select your velvet optimiser contigs file
>
> 2. Look at the output file.
>
>    > <question-title></question-title>
>    >
>    > Compare the output we got here with the output of the simple assemblies obtained in the introductory tutorial.
>    > 1. What are the main differences between them?
>    > 2. Did SPAdes produce a better assembly than the Velvet Optimiser?
>    {: .question}
>
{: .hands_on}
