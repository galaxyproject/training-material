---
layout: tutorial_hands_on
topic_name: your_topic
tutorial_name: your_tutorial_name
---

# Introduction
{:.no_toc}

<!-- This is a comment. -->

Here we will demonstrate genome analyses strategies for understanding structural differences between a newly assembled genome and a set of published, annotated genomes. 

> ### Agenda
>
> In this tutorial we begin with a new genome assembly just produced in the [Unicycler tutorial]({{site.baseurl}}/topics/assembly/tutorials/unicycler-assembly/tutorial.html). This is an assembly of *E. coli* C, which we will be comparing to assemblies of all other complete genes of *E. coli*.
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Finding and loading all complete *E. coli* genomes

*E. coli* is one of the most studied organisms. Naturally, there are hundreds of complete genomes. Here we will shows how to uploaded all (!) complete *E. coli* genomes as once. 

## Preparing the data

Our initial objective is to compare our assembly against all complete *E. coli* genomes to identify the most related ones and to identify any interesting genome alterations. In order to do this we need to align our assembly against all other genomes. And in order to do that we need to first obtain all these other genomes. 

[NCBI](https://www.ncbi.nlm.nih.gov/) is the resource that would store all complete *E. coli* genomes. Specifically, they can be found [here](https://www.ncbi.nlm.nih.gov/genome/genomes/167). As we will see this list contains over 500 genomes and so uploading them by hand will likely result in carpal tunnel syndrome, which we want to prevent. Galaxy has several features that are specifically designed for uploading and managing large sets of similar types of data. The following two **Hands-on** section show how they can be used to import all completed *E. coli* genomes into Galaxy. 

<!--
{% icon hands_on %} will render the hands_on icon as specified in
_config.yml in the root of this repository.
-->

> ### {% icon hands_on %} Hands-on: Preparing a list of all complete *E. coli* genomes
>
>Open [the NCBI list of of *E. coli* genomes](https://www.ncbi.nlm.nih.gov/genome/genomes/167) in a new window and position two browser windows (one the tutorial and the one you just opened) side by side. Then follow the steps in the following video. 
>
>---------------------
>
><div style="padding:56.25% 0 0 0;position:relative;"><iframe src="https://player.vimeo.com/video/271328293?title=0&byline=0&portrait=0" style="position:absolute;top:0;left:0;width:100%;height:100%;" frameborder="0" webkitallowfullscreen mozallowfullscreen allowfullscreen></iframe></div><script src="https://player.vimeo.com/api/player.js"></script>
>
{: .hands_on}

## Getting complete *E. coli* genomes into Galaxy

Now that the list is formatted as a table in a spreadsheet it is time to upload it into Galaxy. There is a problem though &uarr; the URLs (web addresses) in the list do not actually point to sequence files that we would need to perform alignments. Instead they point to directories. For example, this URL:

```
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/008/865/GCA_000008865.1_ASM886v1					
```

points to a directory (rather than a file) containing the following files:

![GenBank assembly files for an E. coli strain](../../images/genbank_dir.png "A list of files for an E. coli assembly. For further analyses we need datasets ending with '_genomic.fna.gz'.")

So to download sequence files we need to edit URLs by adding filenames to them. For example, in the case of URL shown above we need to add `/GCA_000008865.1_ASM886v1_genomic.fna.gz` to the end to get this:

```
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/008/865/GCA_000008865.1_ASM886v1/GCA_000008865.1_ASM886v1_genomic.fna.gz
```

this can be done as a two step process where we first copy the end part of the existing URL (`/GCA_000008865.1_ASM886v1_genomic.fna.gz`) and then add a fixed string `_genomic.fna.gz` to the end of it. Doing by hand is crazy and trying it in a spreadsheet is complicated. Fortunately Galaxy's new rule-based unloader helps with that as shown in the next **Hands-on** section.

Short introduction about this subpart.

> ### {% icon hands_on %} Hands-on: Data upload
>
>Here we copy data from the spreadsheet described in the previous section into Galaxy's rule-based uploader to download several hunder complete genomes into a Collection. Follow the steps in the video below.
>
>----------------------
>
><div style="padding:56.25% 0 0 0;position:relative;"><iframe src="https://player.vimeo.com/video/271336444?title=0&byline=0&portrait=0" style="position:absolute;top:0;left:0;width:100%;height:100%;" frameborder="0" webkitallowfullscreen mozallowfullscreen allowfullscreen></iframe></div><script src="https://player.vimeo.com/api/player.js"></script>
{: .hands_on}

## Preparing assembly

Before doing any analyses we need to upload assembly produced in [Unicycler tutorial]({{site.baseurl}}/topics/assembly/tutorials/unicycler-assembly/tutorial.html) from Zenodo. 

 > ### {% icon hands_on %} Uploading *E. coli* assembly into Galaxy
 >
 > 1. Open upload tool (Upload icon on the top of the left pane)
 > 2. Click **Paste/Fetch data** button (Bottom of the interface box)
 > 3. Paste `https://zenodo.org/record/1251125/files/Ecoli_C_assembly.fna` into the box.
 > 4. Set **Type** to `fasta`
 > 5. Click **Start**
{: .hands_on}

The assembly we just uploaded has two issues that need to be addressed before proceeding with our analysis:

 1. It contains two sequences: the one of *E. coli* C genome (the one we really needed) and another representing phage phiX174 (a by product of Illumina sequencing were it is used a spike in DNA). 
 2. Sequences have unwieldy names like `>1 length=4576293 depth=1.00x circular=true`. We need to rename it to something more meaningful.

 Let's fix these two problems.

 > ### {% icon tip %} Tip: Finding tools mentioned in this tutorial
 >Galaxy instances contain hundreds of tools. As a result it is sometimes hard to find tools mentioned in tutorials such as this one. 
 >
 >To help with this challenge Galaxy has a search box at the top of the left panel. Use this box to find the tools mentioned here.
 >![](../../images/tool_search.png "Use search box to find tools!")
 {: .tip}

> ### {% icon hands_on %} Hands-on: Fixing assembly
>
> 1. First we will use **Filter sequences by length** {% icon tool %} tool to remove phiX 174 genome. Set parameters as follows:
>   - **Fasta file** to the dataset you've just uploaded. (It will have a name something like `https://zenodo.org/record/1251125/files/Ecoli_C_assembly.fna` ).
>   - **Minimal length** to `10000` (this was phiX174, which is around 5,000 bp, will be filtered out)
>   - **Maxumum length** you do not need to change.
> 2. Second we will use **Text transformation with sed** {% icon tool %} tool to convert name of the assembly from `>1 length=4576293 depth=1.00x circular=true` to `>Ecoli_C`. Set parameters as follows:
>   - **File to process** should be set to the output of the previous step
>   - Inside **SED program** box enter the following expression (so called [Regular Expression](https://en.wikipedia.org/wiki/Regular_expression): `s/^\>1.*$/\>Ecoli_C/`
>
>	> ### {% icon tip %} Highlight: SED editor and Regular Expressions
>	>The expression `s/^>1.*$/>Ecoli_C/` contains several pieces that you need to understand. Let's write it top-to-bottom and explain:
>	>
>	> - `s` - tell SED to *Substitute*
>	> - `/` - opens a section of the commands telling SED *what* to substitute. 
>	> - `^` - tell SED to start looking at *the beginning* of each line
>	> - `>` - is the first character we want to match. Remember that name of the sequence in FASTA files starts with `>`
>	> - `1` - is the number present is our old name (`>1 length=4576293 depth=1.00x circular=true` to `>Ecoli_C`)
>	> - `.` - dot has a special meaning. It signifies *any* character. 
>	> - `*` - is a *quantifier*. From [Wikipedia](https://en.wikipedia.org/wiki/Regular_expression): "The asterisk indicates zero or more occurrences of the preceding element. For example, ab*c matches `ac`, `abc`, `abbc`, `abbbc`, and so on."
>	> - `$` - signifies *the end* of a line
> 	> - `/` - is *the end* of the *what to substitute* section. It also serves as the beginning of *what to substitute WITH* section
> 	> - `>` - is the required element of the FASTA sequence name
>	> - `Ecoli_C` is the *name* we want the sequence to have
>	> - `/` - is the end of SED command
>	>
>	>So in short we are replacing `>1 length=4576293 depth=1.00x circular=true` with `>Ecoli_C`. The *Regular expression* `^\>1.*$` is used here to represent `>1 length=4576293 depth=1.00x circular=true`.<br>
>	>Detailed description of regular expressions is outside of the scope of this tutorial, but there are other great resources. Start with [Software Carpentry Regular Expressions tutorial](http://v4.software-carpentry.org/regexp/index.html)!
>	{: .tip}
>
>    > ### {% icon question %} Questions
>    >
>    > 1. What is the meaning of `^` character is SED expression?
>    > 2. Where do you go to learn more about regular expressions?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>It tells SED to start matching from the beginning of the string.</li>
>    >    <li>Software carpentry at <a href="https://software-carpentry.org">https://software-carpentry.org</a></li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
{: .hands_on}

# Generating alignments

Now everything is loaded and ready to go. We will now align our assembly against each of the *E. coli* genomes we have uploaded into the collection. To do this we will use [LASTZ](https://lastz.github.io/lastz/) an aligner designed to align long sequences. 

> ### {% icon hands_on %} Hands-on: Running LASTZ
> 1. Open **LASTZ** interface 
> 2. Change **Select TARGET sequnce(s) to align against** to `from your history`
> 3. In **Select a reference dataset** click on the folder icon (![](../../images/folder-o.png)) and the collection containing all *E. coli* genomes we uploaded below. 
> 4. In **Select QUERY sequence(s)** choose our assembly which was prepared in the previous step.
> 5. Find section of LASTZ interface called **Chaining** and expand it.
> 6. Set **Perform chaining of HSPs with no penalties** to `Yes`
> 7. Run LASTZ by clicking **Execute** button
{: .hands_on}

> ### {% icon warning %} It will take a while!
> Please understand that alignment is not an instantaneous process: allow several hours for these jobs to clear.
{: .warning-box}

# Finding closely related assemblies


Short introduction about this subpart.

> ### {% icon comment %} Comment
>
> Do you want to learn more about the principles behind mapping? Follow our [training](../../NGS-mapping)
{: .comment}

# Conclusion
{:.no_toc}

Conclusion about the technical key points. And then relation between the technics and the biological question to end with a global view.
