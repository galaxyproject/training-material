---
layout: tutorial_hands_on
topic_name: introduction
tutorial_name: galaxy-intro-strands
---

# Introduction to Galaxy 
{:.no_toc}

This practical aims to familiarize you with the Galaxy user interface. It will teach you how to perform basic tasks such as importing data, running tools, working with histories, creating workflows, and sharing your work.

This practical teaches the same basic content as [Galaxy 101](/topics/introduction/tutorials/galaxy-intro-101/tutorial.html), but requires less knowledge of biology to understand the question.

> ###  {% icon comment %} Audience
> This practical is for those who are new to Galaxy, genomics, and bioinformatics.  If you aren't new to bioinformatics you can just do the items listed in the Hands-On boxes ({% icon hands_on %}), or you can try one of the other introductory tutorials.
{: .comment}

> ### Agenda
>
> In this tutorial, we will:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Pretreatments

> ### {% icon requirements %} Requirements
>
> To run this practical you will need
>
> 1. An internet-connected computer.  Galaxy can run on your laptop without an internet connection, but this practical requires access to resources on the web. 
> 1. A web browser. Firefox and Google Chrome work well, as does Safari.  Internet Explorer is known to have issues with Galaxy so avoid using that.
> 1. Access to a Galaxy instance.  Galaxy is available in many ways. If you are doing this practical as part of a workshop, the instructor will tell you which instance to use. If you are doing this on your own, you can use [usegalaxy.org](https://usegalaxy.org).
{: .comment}


> ### {% icon question %} Our Motivating Question
> *I wonder if genes on opposite strands ever overlap with each other, and if so, how common is that?*
{: .question}

To explore this question we need a basic understanding of *genomes, chromosomes, strands, genes,* and *exons*.

> ### {% icon comment %} Definitions 1
>
> * **Genome**
>> The genome is the collection of all DNA native to an organism. For humans, the genome is all of a person's chromosomes. 
>
> * **Chromosome**
>> The largest unit of DNA organization in an organism.  Humans have two copies of 23 chromosomes.  Chromosomes are *linear* in humans, and all animals and plants.  (Bacteria have *circular* chromosomes.)
>
> * **Strand**
>> Chromosomes are *double-stranded*.  One is the forward strand, is typically drawn on top, and moves from left to right. The other, reverse strand, is typically drawn on the bottom and moves from right to left.  Genes can occur on either strand.  A single gene will have parts on only one stand.
>
> * **Gene**
>> "What is a gene?" is actually a hotly debated question.  For our purposes, a gene is a section of DNA on chromosome strand that creates a molecule used by an organism.
>
> * **Exon**
>> In humans (and in all plants and animals) the molecules that are built from genes are often only built from a part of the DNA in a gene. The sections of DNA that can produce the molecules are called exons.
{: .comment}

[ADD IMAGE HERE]

Now lets refine our question slightly

> ### {% icon question %} Our Revised Motivating Question
> *I wonder if **exons** on opposite strands ever overlap with each other, and if so, how common is that?*
{: .question}



Conceptually our question looks like

[IMAGE SHOWING CHROM STRANDS WITH GENES ON BOTH STRANDS. once without overlapping genes, and once with.]

The first case is common.  We want to know if the second case happens.

## Get human gene definitions

To answer this question we need to know where genes start and stop on human chromosomes.  That seems like a simple question, but if you are new to bioinformatics it's actually a hard question to answer.  Web searches will land you at any number of useful places on the web, but without a lot of background knowledge it's hard to know what you want:  *What's the difference between sequence and annotation?  What are FASTA, BED, GTF, GFF3, and VCF?  What are GRCh37, GRCh38, hg19, and hg38 (and what happened to hg20 through hg37 - are they okay)?*

It turns out that for this particular question (and for many others), most **Galaxy** instances can help us find this information.

In your web browser, go to [your Galaxy instance](#-requirements) and *log in or register*.

The Galaxy interface consists of three main parts. The available tools are listed on the left, your analysis history is recorded on the right, and the middle panel will show the home page, tool forms, and dataset content.

![Galaxy interface](../../images/galaxy_interface.png)

> ### {% icon hands_on %} Hands-on: Create history
>
> 1. Make sure you start from an empty analysis history.
>
>    > ### {% icon tip %} Creating a new history
>    >
>    > * Click the **gear icon** at the top of the history panel
>    > * Select the option **Create New** from the menu
>    {: .tip}
>
{: .hands_on}

## Get data into Galaxy

There are [many ways to get data into a Galaxy instance](/topics/introduction/tutorials/galaxy-intro-get-data/slides.html). We are going to use the **Get Data** toolbox in the **Tools** panel on the left.

> ### {% icon hands_on %} Open **Get Data** toolbox
>
>  ![The Get Data toolbox](../../images/101_01.png)  *Click* on the **Get Data** toolbox to expand it.
{: .hands_on}

The **Get Data** toolbox contains a list of data sources that are available on this server.  **Upload file** is quite useful for getting data from your computer or from the web (see the [Getting data into Galaxy slides](topics/introduction/tutorials/galaxy-intro-get-data/slides.html#4).  Today we are going to use the **UCSC Main table browser**.

### Get exons

> ### {% icon hands_on %} Go to UCSC
>
> ![Click on UCSC Main table browser](../../images/101_01.png) *Click* on **UCSC Main table browser** to go to UCSC.
{: .hands_on}

This will take you to the UCSC Table Browser:

 ![UCSC table browser tool, first screen for genes](../../images/101_02.png)

The [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables) provides access to all the data that is shown in the [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgGateway) (see [box](#-ucsc-genome-browser) below). If you are working on a species that UCSC supports (like human) then the Table Browser is a great place to get genomic data.

The Table Browser has a daunting number of options. Fortunately,they are all set to commonly used defaults, greatly simplifying things, and most of the options are already set to what we want:

* **clade:** `Mammal`
* **genome:** `Human`
* **assembly:** `Dec. 2013 (GRCh38/hg38)`
* **group:** `Genes and Gene Predictions`
* **track:** `GENCODE v24`

**clade** and **genome** seem pretty clear.  **assembly** asks which version/definition of the human genome we want.  (Any will do for our question, but UCSC is suggesting `hg38`, which is also the most recent.)  **group** is set to `Genes and Gene Predictions` which sounds like what we want. So far so good.

**track** has a bewildering list of options. UCSC suggests `GENCODE v24`.  A web search leads us to the [GENCODE web site](https://www.gencodegenes.org/) which prominently states:

> The GENCODE project produces high quality reference gene annotation and ...

Time for a few more definitions.

> ### {% icon comment %} Definitions
>
> * **Reference genome**
>> A reference genome is the *genome of a single individual* that has been thoroughly studied, to the point that we know exactly what most of that individual's DNA is.  In practice a reference genome is used as shared map by researchers working on that organism. Reference genomes are updated periodically as techniques improve.
>
> * **Sequence**
>> A genome's sequence describes the DNA in that genome, down to the A, C, T, and G (single nucleotide) level including the exact location where each is.  Given a reference genome, you can ask questions like, "What's the DNA on chromosome 2 between positions 1,678,901 and 1.688,322?"
>
> * **Genome/Gene annotation**
>> The sequence tells us what DNA is where, but it doesn't tell us anything about the function of that DNA.  *Annotation* is additional information about particular regions of the genome like where genes, repeats, promotors, and centromeres are, or how active a particular gene is. 
{: .comment}

The **track** option asks us which set of annotations do we want to get?  There are so many choices because annotation is the result of analysis and interpretation, and there are many ways to do this. (And in this case, many of the options aren't even genes or gene predictions.) 

GENCODE is "high-quality" and  "gene annotation." That sounds like a good thing to use.  Lets stay with the default: `GENCODE V24`.

So far we haven't changed *anything* from the defaults.  Lets change something.  The default  **region** is the whole genome, which can be done, but it's a lot of information. For this exercise lets use just one (small) chromosome.

> ### {% icon hands_on %} Limit the region and get the data.
>
> * Say that we just want chromosome 22 
>   * For **region** select `position`.
>   * In the text box next to `position` enter `chr22` (case matters).
>
> * *Click* the **get output** button.
>   * And, that doesn't actually get us the output.  It sends us to a second UCSC page that asks us exactly what we want. Which is good, because we don't want the default.
>
> * Under **Create one BED record per** *select* **Coding exons**.
>   * This will get only the parts of genes that actually produce active molecules.
>
> * *Click* the **Send query to Galaxy** button at the bottom of the form.
{: .hands_on}

This returns us to Galaxy, first displaying a big green box (that's good!) and then returning us to the view we started with.  Excepnt that we now have an item in our history, the dataset from UCSC.

> ### {% icon comment %} UCSC Genome Browser
>
> If you aren't yet familiar with the [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgGateway), then it's worth spending some time learning how to use it.  *Genome browsers* are software for viewing genomic information graphically.  The UCSC Genome Browser (and most genome browsers) typically display different types of *annotation* about a region of a genome.  This is displayed as a stack of *tracks* and each track contains a different type of information. 
>
> Genome browsers are useful for seeing information in context and for seeing (and discovering) correlations between different types of information.
>
> The UCSC Genome Browser has information on over 100 animals, and their [Archaeal Genome Browser](http://archaea.ucsc.edu/cgi-bin/hgGateway?db=pyrFur2) has genomic information on well over 100 microbial species.
{: .comment}

### History Item Status

Watch your new history item.  It will go through three statuses before it's done.

* **Grey**: Item is waiting to start (wating for data transfer to start)
* **Yellow**: Item is running (data is actively being transferred).
* **Green**: Item has finished successfully (data transfer complete).

Occassionally you will also see a 4th status

* **Red**: Item did not finish successfully.

See the *Galaxy History Item Status* practical for more. [TODO]

## Examine the data

> ### {% icon hands_on %} Look at the data.
> Once the dataset is green, *click* on the dataset name (something like **UCSC Main on Human...**)
{: .hands_on}

This expands the dataset and shows you information about it, and a preview of its contents.

[IMAGE of PREVIEW]

1. The preview tells us  several things:
1. The dataset has over 4000 regions, meaning that there are over 4000 genes on chromosome 22.
1. The dataset is in **bed** format.  BED is one of several standard formats for representing genome annotation.  BED is a tabular format that we'll expand on below.  We got BED format because BED was preselected as the output format in the UCSC table browser.
1. The dataset's "database" is **hg38**.  This says which revision of the reference genome this data maps too.  hg38 is the latest human reference genome.  hg38 was also selected by default in UCSC.
1. Finally, it shows us the first 5 rows in the dataset.

The dataset preview is informative, but you can't see much of the actual dataset.  Lets use one of the dataset icons to see the whole dataset

> ### {% icon hands_on %} Look at all the data.
> *Click* on the **eye icon** to view the contents of the dataset.
{: .hands_on}

It will look something like this:
![Contents of the `UCSC Main on Human: knownGene` dataset](../../images/101_exons.png)

### BED Format

BED is one of several well-established tabular formats for genomic data.  Other formats include GFF3 and GTF.  For the type of analysis we are doing today, BED format is easiest to work with.

BED was created to power the UCSC genome browser.  BED files contains between 3 and 15 columns.  Our example BED file describes genes and contains 12 columns.

We care about columns 1, 2, 3, and 6:

| # | Column Name | Meaning |
| ---- | ---- | ---- |
| 1 | Chromosome | The name of the chromosome this gene is on. |
| 2 | Start | Where on the chromosome the gene starts. |
| 3 | End | Where on the chromosome the gene ends. |
| 6 | Strand | Which strand the gene is on.  "+" means forward (top, left to right), "-" means reverse (bottom, right to left) |

See the [BED format description at UCSC](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) for a full description of all the columns.

## Naming

Galaxy allow you to name your analyses (your histories) and your datasets.  We only have one history ("Unnamed history") and one dateset ("UCSC Main on Human:...") so far, but it's a good idea to

1. Always name your histories
2. Name your input and final output datasets, and any significant intermediate datasets.

You don't have to do this.  Galaxy is quite happy for you to have an infinite number of "Unnamed history" histories, and to have all your datasets be obscurely named.  However, once you've run your first 5 unnamed analyses, all with obscurely named datasets, you'll might wish you would have named everything.

> ### {% icon hands_on %} Name your stuff
> 
> 1. **Name your history** to be meaningful and easy to find.
>    - *Click* on the title of the history and enter something like **Intro - Strands** as the name.  Hit the `enter` key on your keyboard to save it.
>   ![Rename the history](../../../../shared/images/rename_history.png)
> 1. **Rename your dataset** 
>    - *Click* on the **pencil icon** to edit the dataset attributes.
>    - In the next screen change the name of the dataset to something like `Genes` or `Genes chr22`.
>    - *Click* the **Save** button at the bottom of the screen.
>
>    Your history should now look something like this:
>
>    ![Rename dataset](../../images/101_rename.png)
{: .hands_on}

## We've got the data - what's our plan for answering the question?

You have to know what's possible, before you can build a plan.  If you don't have experience with data analysis then you might not have any idea how you would answer our question.  Before we dive in using a particular solution, think about how you might solve this.  If you don't have any experience with tools, then think about how you might solve it manually, using pencil and paper (it will help to assume you have an infinite supply of helpers to do the actual work).

Here's a high level description of how we'll answer this question.

1. Split the genes dataset in two: one for genes on the forward strand, and one for genes on the reverse strand.
1. Compare the two datasets to see which ones, if any, overlap.
1. Figure out how many (or what percentage) or our genes overlap with another gene.

It turns out that all of these steps are easy in Galaxy.

### Split the genes into forward and reverse datasets

How might we do this?  Column 6 contains the strand information.  Can we split genes into two datasets based on the value of Column 6.  How?  Lets take a look at our available tools.  And *whoa! There are over 40 toolboxes, and several hundred tools.* How are we going to find a tool that can do the split?

First, we can try the *tools search box.*  Think of terms that might describe what we want to do and type them in the search box.  Do you see anything promising?  Explore a little.

If you haven't already searched with it, enter `split` in the search box.  Near the top of the results is

> **Filter** data on any column using simple expressions.

That might work.

> ### {% icon hands_on %} Open the Filter tool
> * *Click* on **Filter** to open the Filter tool in the middle panel.
> * Take a look at the **Syntax** and **Example** sections to understand what the tool does.
{: .hands_on}

It doesn't say anything about Filter being able to split a file into multiple files.  It does look like we can use Filter to get only genes on the forward strand, or only genes on the reverse strand.  We would have to run Filter twice, once for forward strand genes, and once for reverse strand genes. Let's do that.

The tool form for Filter looks like

[IMAGE]

(You may have noticed during your search for tools that *all* tools have a similar look and feel.)

> ### {% icon hands_on %} Run the Filter tool to get genes on the forwrd strand.
> 
> The filter tool has 3 fields:
>
> 1. **Dataset**: This pulldown will list any dataset from your history that this tool can work on.  In your case that's probably only one dataset.  Make sure this is set to your `Genes` dataset.
> 1. **Condition**: this free text field is where we specify which records we want in the output dataset.  *Enter* `c6 == "+"` in the text box.
>   * This specifies that column 6 (the strand) must be equal to (`==` is Python for *is equal to*) a plus sign.
> 1. **Header lines to skip**: Leave this as `0`. Our dataset does not have any header lines.
>
> Finally, *click* the **Execute** button.
{: .hands_on}

This adds another dataset to your history.  This one should contain only genes on the forward strand.  Once the dataset is green, *click* the eye icon to confirm this.  We also recommend that you rename this dataset to something like `Genes, forward strand` (remember how?).

Now we want to get the genes on the reverse strand.  There are actually many ways to get this.  Here are two of them.

> ### {% icon hands_on %} Get genes on the reverse strand
> 
> **Method 1**
>
> 1. Open the dataset preview by *clicking* on the name of the `Genes, forward strand` dataset.  This show a different set of icons than the uploaded `Genes` dataset did.
> 1. *Click* the **looping arrow** ("Run this job again") icon.  This won't actually run the job again.  What it will do is bring up the Filter tool form with *the exact same settings that were used to produce this dataset.*
> 1. Rather than run Filter again with the same settings, *change* **Condition** to `c6 == "-"`
> 1. Finally, *click* the **Execute** button.
>
> **Method 2**
> 1. *Click* on **Filter** in the tool panel to open the Filter tool in the middle panel.
> 1. *Fill* the form as before, *except*:
>    * Make sure the **Dataset** pulldown is set to the `Genes` dataset.
>    * *Set* **Condition** to `c6 == "-"`.
> 1. Finally, *click* the **Execute** button.
{: .hands_on}

The rerun button can be a huge help as you run more complex tools.

> ### {% icon tip %} Empty result?
>
> If you used Method 2 and didn't explicitly set the dataset, then you ran Filter on the `Genes, forward strand` dataset. None of the genes in the forward strand dataset have "-" in column 6 so all of them were filtered out from the resutlt.
>
> Try again and set the dataset to your `Genes` dataset.
{: .tip}

> ### {% icon hands_on %}
> * *Rename* your new dataset to something like `Genes, reverse strand` 
{: .hands_on}

Your history should now have (at least) 3 datasets in it, with names something like:

* `Genes`
* `Genes, forward strand`
* `Genes, reverse strand`

The number of genes in the `forward` plus `reverse` datasets should be the same as in the `Genes` dataset.  If they aren't can you figure out why?

### Check for overlaps

Genes are an example of a *genomic interval*.

> ### {% icon comment %} Definitions 3
>
> * **Genomic interval**
>> In Galaxy, a *genomic interval* is a something that spans part of a chromosome (or some other frame of reference).  Genes and exons are common examples of genomic intervals.  Even a chromosome is a genetic interval, albeit a very long one.
{: .comment}

Galaxy excels at answering questions about genomic intervals and different sets of genomic intervals relate to each other.  Lets take a look.

> ### {% icon hands_on %} Genomic Interval Tools
>
> * In the tool panel, *open* the **Operate on Genomic Intervals** toolbox.  It's typically past the **NGS** toolboxes.
> * *Explore* the tools in this toolbox, looking for something that we can use to see which genes on opposite strands overlap.
{: .hands_on}

Of the tools in the **Operate on Genomic Intervals** toolbox, **Join** and particularly **Intersect** have the most promise.  Let's try **Intersect**.

> ### {% icon hands_on %} Genomic Interval Tools
>
> * Run the **Intersect** tool
>   * In the tool panel, *click* **Intersect** in the **Operate on Genomic Intervals** toolbox.
>   * *Set* **Return** to `Overlapping Intervals`.
>     > This looks like it might return whole genes, while `Overlapping pieces` may return only the parts that overlap.  We suspect that whole genes might be more useful.
>   * *Set* **of** (the first dataset) to `Genes, forward strand`
>   * *Set* **that intersect** (the second dataset) to `Genes, reverse strand`
>   * *Set* **for at least** to `1`
>     > This will return genes with even just one position overlapping.
>   * *Click* **Execute**.
>
> * Now repeat the intersect, but make the first dataset be the reverse genes, and the second be the forward genes.
>
> * Finally give both of the new datasets meaningful names, like `Overlapping forward genes` and `Overlapping reverse genes`
{: .hands_on}

At this point, we have answered our question: Compare the number of genes in the `Overlapping forward` and `Overlapping reverse` datasets with the number of genes in the full `Genes` dataset.  Are there genes that overlap?  If so, how common is it?

## Try something

Run again with the entire genome.

Exons!  Eek!  What if we want to ask this question, where we only look for genes with exons that overlap?  How can we do that?  (Hint, start by going back to UCSC, and selecting exons.)

Are there genes on the same strand that overlap with each other?

## Final thoughts

### Why not use Excel for this?

You could use Excel or another spreadsheet program to do these calculations.  Here, we learned how to use Galaxy by answering a question.  You could just as easily learn Excel by answering the same question, and if the goal is to learn how to do something, then either would be great. But what if you are working on a question where your analysis matters?  Maybe you are working with human clinical data trying to diagnose a set of symptoms, or you are working on research that will eventually be published and maybe earn you a Nobel Prize?

In these cases your analysis, *and the ability to reproduce it exactly*, is vitally important.  Excel won't help you with this. It doesn't track changes and it offers very little insight to others on how you got from your initial data to your conclusions.

Galaxy, on the other hand, *automatically records every step of your analysis.*  And when you are done, you can share your analysis with anyone.  You can even include a link to it in a paper (or your acceptance speech).  In addition, you can create a reusable recipe (a "workflow" in Galaxy) from your analysis that others (or yourself) can use on other datasets.

Another challenge with spreadsheet programs is that they don't scale to support next generation sequencing datasets, which often reach gigabytes or even terabytes in size.

# Conclusion
{: .no_toc}

:tada: Well done! :clap: You have just performed your first analysis in Galaxy. 
