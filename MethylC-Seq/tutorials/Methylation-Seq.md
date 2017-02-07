---
layout: tutorial_hands_on
topic_name: DNA methylation
tutorial_name: dna_methylation
---

# Galaxy workshop on DNA Methylation data analysis

> ### Agenda
>
> In this tutorial we will:
>
> 1. [Load data and quality control](#Load\ data\ and\ quality\ control)
> 2. [Align the data](#alignment)
> 3. [Methylation bias and metric extraction](#Methylation\ bias\ and\ metric\ extraction)
> 4. [Visualize the mapped data](#visualization)
> 5. [Metilene](#Metilene)
> 
> 
> {: .agenda}


# Load data and quality control
> ### :pencil2: Hands-on: Get the data and look at the quality
> 
> We load now one example data set which will be used for the tutorial. 
>
> 1. Load the dataset from: XXX
>
> 2. **FastQC**
> 
>   > ### :bulb: Tip: Search for tools
>    >
>    > * Clink into the search field on the left
>    > * Type **fastqc**
>    > * Select **FastQC**
>    > * Select the uploaded dataset as the fastq file
>    > * Choose as a reference genome XXX
>    {: .tip}
> The computation will take a while and we continue with the theory. 
>
> 3. Go to the webpage result page and have a closer look at 'Per base sequence content'
>
>    > ### :question: Questions
>    >
>    > - Note the GC distribution and percentage of "T" and "C". Why is this so weird?
>    > - Is everything as expected?
>    > 
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>The attentive audience of the theory part knows: Every C-meth stays a C and every normal C becomes a T during the bisulfite conversion. </li>
>    >    <li>Yes it is. Always be careful and have the specific characteristics of your data in mind during the interpretation of FastQC results.</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
{: .hands_on}

# Alignment

> ### :pencil2: Hands-on: Mapping with bwameth
> 
> We will map now the imported dataset against a reference genome.
> 
> 1. **Galaxy** :wrench:: Search for the tool 'bwameth'
> 2. **bwameth** :wrench:: Chose 
>
>    > ### :question: Questions
>    >
>    > -  Why we need other alignment tools for bisulfite sequencing data?
>    > 
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>You may have noticed that all the C's are C-meth's and a T can be a T or a C. A mapper for methylation data needs to find out what is what.</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
{: .hands_on}

# Methylation bias and metric extraction

> ### :pencil2: Hands-on: Methylation extraction with PileOMeth / MethylDackel
> 
> We will extract the methylation on the resulting BAM file of the alignment step.
> 
> 1. **Galaxy** :wrench:: Search for the tool 'PileOMeth'
> 2. **PileOMeth** :wrench:: Chose 
>
>    > ### :question: Questions
>    >
>    > - Is the methylation bias as expected? 
>    > - 
>    > 
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>Some answer.</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
> 
> You now have a visualization of if/where there's methylation bias and modified bedGraph files with methylation metrics in them. This could be used for downstream statistical analysis (typically in an R package).
>
{: .hands_on}

# Visualization 

> ### :pencil2: Hands-on: Methylation extraction with PileOMeth / MethylDackel
> 
> We will extract the methylation on the resulting BAM file of the alignment step.
> 
> 1. **Galaxy** :wrench:: Search for the tool 'PileOMeth'
> 2. **PileOMeth** :wrench:: Chose 
>
>    > ### :question: Questions
>    >
>    > - Is the methylation bias as expected? 
>    > - 
>    > 
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>Some answer.</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
> 
> You now have a visualization of if/where there's methylation bias and modified bedGraph files with methylation metrics in them. This could be used for downstream statistical analysis (typically in an R package).
>
{: .hands_on}

# Metilene 

> ### :pencil2: Hands-on: Metilene
> 
> We will extract the methylation on the resulting BAM file of the alignment step.
> 
> 1. **Galaxy** :wrench:: Search for the tool 'PileOMeth'
> 2. **PileOMeth** :wrench:: Chose 
>
>    > ### :question: Questions
>    >
>    > - Is the methylation bias as expected? 
>    > - 
>    > 
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>Some answer.</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
> 
> You now have a visualization of if/where there's methylation bias and modified bedGraph files with methylation metrics in them. This could be used for downstream statistical analysis (typically in an R package).
>
{: .hands_on}
