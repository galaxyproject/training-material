---
layout: tutorial_hands_on

title: "Using data from NCBI's Sequence Read Archive (SRA): A
SARS-CoV-2 variant analysis example"
zenodo_link: ''
questions:
- Which biological questions are addressed by the tutorial?
- Learn how to get and use data from the Sequence Read Archive in Galaxy.
objectives:
- Understand how Galaxy and the Short Read Archive interact.
- Be able to go from Galaxy to the Short Reach Archive, query SRA, use the SRA Run Selector to send selected metadata to Galaxy, and then import sequence data from SRA into Galaxy.
time_estimation: ''
key_points:
- Sequence data in the SRA can be directly imported into Galaxy.
contributors:
- mvdbeek
- tnabtaf
- blankenberg

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

The aim of this tutorial is to learn how to Galaxy and NCBI's Short Read Archive interact (SRA) with each other.  This tutorial uses a COVID-19 variant calling example, but it isn't about variant calling per se.

SRA, like many data sources is Galaxy aware.  It has support for sending information directly to Galaxy, and it tracks which Galaxy instance it was invoked from.  Getting sequence data from SRA is a multi-step process.  This tutorial explains each step in the process and then demonstrates a particular example of how to use SRA data in Galaxy.

At the completion of this tutorial you know

* how to go from Galaxy to SRA.
* a basic understanding of how to select data in SRA
  * see the [Search in SRA documentation](https://www.ncbi.nlm.nih.gov/sra/docs/srasearch/) from the SRA team for a fuller introduction to this topic
* how to send SRA *metadata* (such as accession numbers) to Galaxy
* how to use that SRA metadata to import *sequence* data from SRA into Galaxy.
* how to run a simple variant analysis in Galaxy using that data
  * See these [variant analysis tutorials](/training-material/topics/variant-analysis) for a more in-depth explanation of variant analysis.

{% raw %} `{% cite Batut2018 %}`{% endraw %}.
This will be rendered like this: {% cite Batut2018 %}, and links to a
[bibliography section](#bibliography) which will automatically be created at the end of the
tutorial.


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# The Sequence Read Archive

The [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) is the primary archive of *unassembled reads*  for the [US National Institutes of Health (NIH)](https://www.ncbi.nlm.nih.gov/).  SRA is a great place to get the sequencing data that underlie publications and studies.

This tutorial covers how to get sequence data from SRA into Galaxy using a direct connection between the two.

> ### {% icon comment %} Comment
>
> You will also hear SRA referred to as the *Short Read Archive*, its original name.
>
{: .comment}


## Accessing SRA

SRA can be reached either directly through it's website, or through the tool panel on Galaxy. 

> ### {% icon comment %} Comment
>
> Initially the tool panel option exists only on the [usegalaxy.org server](https://usegalaxy.org/).  Support for the direct connection to SRA will be included in the 20.05 release of Galaxy
{: .comment}


> ### {% icon hands_on %} Hands-on: Explore SRA Entrez
>
> 1. Go to [usegalaxy.org](https://usegalaxy.org/)
> 1. If your history is not already empty, than start a new history (see here for more on Galaxy histories)
> 1. **Click** `Get Data` at the top of the tool panel.
> 1. **Click** `SRA Server` in the list of tools shown under `Get Data`.
>    This takes you the [Sequence Read Archive home page](https://www.ncbi.nlm.nih.gov/sra).  A search box is shown at the top of the page.  Try searching for something you are interested in, such as `dolphin` or `kidney` or `dolphin kidney` and then **click** the  `Search` button.
>
>    This returns a list of *SRA Experiments* that match your search string.  SRA Experiments, also know as *SRX entries*, contain sequence data from a particular experiment, as well as an explanation of the experiment itself and any other related data. You can explore the returned experiments by clicking on their name.  See [Understanding SRA Search Results](https://www.ncbi.nlm.nih.gov/books/NBK56913/) in the [SRA Knowledge Base](https://www.ncbi.nlm.nih.gov/books/n/helpsrakb/) for more.
>
>    When you enter text in the SRA search box, you are using [SRA's Entrez search interface](https://www.ncbi.nlm.nih.gov/sra/docs/srasearch/).  Entrez supports both simple text searches, and very precise searches that check specific metadata and use arbitrarily complex logical expressions.  Entrez allows you to scale up your searches from basic to advanced as you narrow your searches.  The syntax of advanced searches can seem daunting, but SRA provides a graphical [Advanced Search Builder](https://www.ncbi.nlm.nih.gov/sra/advanced/) to generate the specific syntax.  And, as we shall see below, the SRA Run Selector provides an even friendlier user interface for narrowing our selected data.
>
>    Play around with the SRA Entrez interface, including the advanced query builder, to see if you can identify a set of SRA experiments that are relevant to one of your research areas.
{: .hands_on}


> ### {% icon hands_on %} Hands-on: Generate list of matching experiments using Entrez
> 
> Now that you have a basic familiarity with SRA Entrez, let's find the sequences used in this tutorial.
>
> 1. If you aren't already there, **navigate** back to the [Sequence Read Archive search page](https://www.ncbi.nlm.nih.gov/sra)
> 1. **Clear** any search text from the search box.
> 1. ***TODO***:**Type** `our excellent first search` in the search box and **click** `Search`.
>
>    This returns a longish list of SRA experiments that match our search, and that list is far too long to use in a tutorial exercise.  At this point we could use the advanced Entrez query builder we learned about above.
>
>    But we won't.  Instead lets send the *too long for a tutorial* list results we have to the SRA Run Selector, and use its friendlier interface to narrow our results.
{: .hands_on}


> ### {% icon hands_on %} Hands-on: Go from Entrez to SRA Run Selector
>
> This text appears in a box at the top of the search results
> 
> ***TODO*** View results as an expanded interactive table using the RunSelector.  <u>Send results to Run selector</u>
>
> > ### {% icon comment %} What if you don't see the Run Selector Link?
> >
> > You may have noticed this text earlier when you were exploring Entrez search.  This text only appears some of the time, when the number of search results falls within a fairly broad window.  You won't see it if you only have a few results, and you won't see it if you have more results than the Run Selector can accept.
> >
> > *You need to get to Run Selector to send your results to Galaxy.* What if you don't have enough results to trigger this link being shown?  In that case you call get to the Run Selector by **clicking** on the `Send to` pulldown menu at the top right of the results panel.  To get to Run Selector, **select** `Run Selector` and then **click** the `Go` button. 
> {: .comment}
>
>
> 1. **Click** `Send results to Run selector` at the top of the search results panel. (If you don't see this link, then see the comment directly above.)
{: .hands_on}

## SRA Run Selector

We learned earlier how to narrow our search results by using Entrez's advanced syntax.  However, we didn't take advantage of that power when we were in Entrez.  Instead we used a simple search and then sent all the results to the Run Selector.  We don't yet have the (short) list of results we want to run analysis on. *What are we doing?*

Well, we are doing something fairly common, and using Entrez and the Run Selector how they are designed to be used:

 * Use the Entrez interface to narrow your results down to a size that the Run Selector can consume.
 * Send those Entrez results to the SRA Run Selector
 * Use the Run Selector's much friendlier interface to
    1. More easily understand the data we have
    1. Narrow those results using that knowledge.

 
> ### {% icon comment %} Run Selector is both more and less than Entrez
>
> Run Selector can do most, but not all of what Entrez search syntax can do.  Run selector uses *faceted search* technology which is easy to use, and powerful, but which has inherent limits.  Specifically, Entrez will work better when searching on attributes that have tens, hundreds, or thousands of different values.  Run Selector will work better searching attributes with fewer than 20 different values.  Fortunately, that describes most searches. 
{: .comment}


The Run Selector window is divided into several panels:

**`Filters List`**: In the upper left hand corner.  This is where we will refine our search.
**`Select`**: A summary of what was initially passed to Run Selector, and how much of that we have selected so far.  (And so far, we haven't selected any of it.)  Also note the tantalizing, but still grayed out, `Galaxy` button.
**`Found x Items`** Initially, this is the list of items sent to Run Selector from Entrez.  This list will shrink as we apply filters to it.


> ### {% icon comment %} Why did the number of found items *go up?*
>
> Recall that the Entrez interface lists SRA experiments (SRX entries).  Run Selector lists *runs* &mdash; sequencing datasets &mdash; and there are *one or more* runs per experiment. We have the same data as before, we are now just seeing it in finer detail.
{: .comment}




## Get data: The SRA data source

### Search SRA

[Search in SRA documentation](https://www.ncbi.nlm.nih.gov/sra/docs/srasearch/)

Search for `sars-cov-2`.  And we get back over 24,000 results (as of June 2020).  Refine that search by adding `"library layout paired"[Properties]` ([from SRA search documentation](https://www.ncbi.nlm.nih.gov/sra/docs/srasearch/)). That reduces our result size by almost 75% and has the added advantage of having all data in the same paired end format.

Wait, wait, wait.  This is viral sequence.  Does it matter if we use paired end or single end?

"platform oxford_nanopore"[Properties] 

Need someway to narrow the results that make sense, and makes downstream analysis simpler too.  We could narrow results here, or in run selector or both.

Maybe demo syntax in Entrez, give up and say there is a simpler way in the run selector.  Pass everything there.


### Switch to the Run Selector

`View results as an expanded interactive table using the RunSelector. Send results to Run selector`

Note the filters list box in the upper left.

Good time to talk about DATASTORE filetype?

Host and Host_scientific name




### Send to Galaxy

## It's not sequence, it's *metadata!*

That section title won't work  Doesn't translate.

explain how to see what we got, and then what we got, and that we still need the sequence

Search for SRA in tools.  Faster looks good.

Wait, it says we need just the accession number column.


### Get just the accession ids

Two operations
Keep column 1 only
Remove leading lines



## Now, get the sequence data

Back to Faster


# Cool! Now what can we do with it?

TBD


> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]() or from the shared data library
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}

## 


# Title of the section usually corresponding to a big step in the analysis

It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:


## Sub-step with **Faster Download and Extract Reads in FASTQ**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Faster Download and Extract Reads in FASTQ** {% icon tool %} with the following parameters:
>    - *"select input type"*: `List of SRA accession, one per line`
>        - {% icon param-file %} *"sra accession list"*: `output` (Input dataset)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **fastp**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **fastp** {% icon tool %} with the following parameters:
>    - *"Single-end or paired reads"*: `Paired Collection`
>        - {% icon param-file %} *"Select paired collection(s)"*: `list_paired` (output of **Faster Download and Extract Reads in FASTQ** {% icon tool %})
>    - In *"Output Options"*:
>        - *"Output JSON report"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Map with BWA-MEM**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Map with BWA-MEM** {% icon tool %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `output` (Input dataset)
>    - *"Single or Paired-end reads"*: `Paired Collection`
>        - {% icon param-file %} *"Select a paired collection"*: `output_paired_coll` (output of **fastp** {% icon tool %})
>    - *"Set read groups information?"*: `Do not set`
>    - *"Select analysis mode"*: `1.Simple Illumina mode`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **MarkDuplicates**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **MarkDuplicates** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select SAM/BAM dataset or dataset collection"*: `bam_output` (output of **Map with BWA-MEM** {% icon tool %})
>    - *"If true do not write duplicates to the output file instead of writing them with appropriate flags set"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Realign reads**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Realign reads** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Reads to realign"*: `outFile` (output of **MarkDuplicates** {% icon tool %})
>    - *"Choose the source for the reference genome"*: `History`
>        - {% icon param-file %} *"Reference"*: `output` (Input dataset)
>    - In *"Advanced options"*:
>        - *"How to handle base qualities of 2?"*: `Keep unchanged`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Samtools stats**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Samtools stats** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"BAM file"*: `outFile` (output of **MarkDuplicates** {% icon tool %})
>    - *"Set coverage distribution"*: `No`
>    - *"Output"*: `One single summary file`
>    - *"Filter by SAM flags"*: `Do not filter`
>    - *"Use a reference sequence"*: `No`
>    - *"Filter by regions"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Insert indel qualities**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Insert indel qualities** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Reads"*: `realigned` (output of **Realign reads** {% icon tool %})
>    - *"Indel calculation approach"*: `Dindel`
>        - *"Choose the source for the reference genome"*: `History`
>            - {% icon param-file %} *"Reference"*: `output` (Input dataset)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Call variants**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Call variants** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input reads in BAM format"*: `output` (output of **Insert indel qualities** {% icon tool %})
>    - *"Choose the source for the reference genome"*: `History`
>        - {% icon param-file %} *"Reference"*: `output` (Input dataset)
>    - *"Call variants across"*: `Whole reference`
>    - *"Types of variants to call"*: `SNVs and indels`
>    - *"Variant calling parameters"*: `Configure settings`
>        - In *"Coverage"*:
>            - *"Minimal coverage"*: `50`
>        - In *"Base-calling quality"*:
>            - *"Minimum baseQ"*: `30`
>            - *"Minimum baseQ for alternate bases"*: `30`
>        - In *"Mapping quality"*:
>            - *"Minimum mapping quality"*: `20`
>    - *"Variant filter parameters"*: `Preset filtering on QUAL score + coverage + strand bias (lofreq call default)`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SnpEff eff:**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **SnpEff eff:** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Sequence changes (SNPs, MNPs, InDels)"*: `variants` (output of **Call variants** {% icon tool %})
>    - *"Output format"*: `VCF (only if input is VCF)`
>    - *"Create CSV report, useful for downstream analysis (-csvStats)"*: `Yes`
>    - *"Annotation options"*: ``
>    - *"Filter output"*: ``
>    - *"Filter out specific Effects"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SnpSift Extract Fields**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **SnpSift Extract Fields** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Variant input file in VCF format"*: `snpeff_output` (output of **SnpEff eff:** {% icon tool %})
>    - *"Fields to extract"*: `CHROM POS REF ALT QUAL DP AF SB DP4 EFF[*].IMPACT EFF[*].FUNCLASS EFF[*].EFFECT EFF[*].GENE EFF[*].CODON`
>    - *"multiple field separator"*: `,`
>    - *"empty field text"*: `.`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **MultiQC**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **MultiQC** {% icon tool %} with the following parameters:
>    - In *"Results"*:
>        - {% icon param-repeat %} *"Insert Results"*
>            - *"Which tool was used generate logs?"*: `fastp`
>                - {% icon param-file %} *"Output of fastp"*: `report_json` (output of **fastp** {% icon tool %})
>        - {% icon param-repeat %} *"Insert Results"*
>            - *"Which tool was used generate logs?"*: `Samtools`
>                - In *"Samtools output"*:
>                    - {% icon param-repeat %} *"Insert Samtools output"*
>                        - *"Type of Samtools output?"*: `stats`
>                            - {% icon param-file %} *"Samtools stats output"*: `output` (output of **Samtools stats** {% icon tool %})
>        - {% icon param-repeat %} *"Insert Results"*
>            - *"Which tool was used generate logs?"*: `Picard`
>                - In *"Picard output"*:
>                    - {% icon param-repeat %} *"Insert Picard output"*
>                        - *"Type of Picard output?"*: `Markdups`
>                        - {% icon param-file %} *"Picard output"*: `metrics_file` (output of **MarkDuplicates** {% icon tool %})
>        - {% icon param-repeat %} *"Insert Results"*
>            - *"Which tool was used generate logs?"*: `SnpEff`
>                - {% icon param-file %} *"Output of SnpEff"*: `csvFile` (output of **SnpEff eff:** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
