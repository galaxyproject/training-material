---
layout: tutorial_hands_on

title: Preparing genomic data for phylogeny reconstruction
zenodo_link: https://zenodo.org/record/6524847#.YnUycVxByV4
questions:
- Which biological questions are addressed by the tutorial?
- Which bioinformatics techniques are important to know for this type of data?
objectives:
- The learning objectives are the goals of the tutorial
- They will be informed by your audience and will communicate to them and to yourself
  what you should focus on during the course
- They are single sentences describing what a learner should be able to do once they
  have completed the tutorial
- You can use Bloom's Taxonomy to write effective learning objectives
time_estimation: 3H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- contributor1
- contributor2

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

General introduction about the topic and then an introduction of the
tutorial (the questions and the objectives). It is nice also to have a
scheme to sum up the pipeline used during the tutorial. The idea is to
give to trainees insight into the content of the tutorial and the (theoretical
and technical) key concepts they will learn.

You may want to cite some publications; this can be done by adding citations to the
bibliography file (`tutorial.bib` file next to your `tutorial.md` file). These citations
must be in bibtex format. If you have the DOI for the paper you wish to cite, you can
get the corresponding bibtex entry using [doi2bib.org](https://doi2bib.org).

With the example you will find in the `tutorial.bib` file, you can add a citation to
this article here in your tutorial like this:
{% raw %} `{% cite Batut2018 %}`{% endraw %}.
This will be rendered like this: {% cite Batut2018 %}, and links to a
[bibliography section](#bibliography) which will automatically be created at the end of the
tutorial.


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

# Title for your first section

Give some background about what the trainees will be doing in the section.
Remember that many people reading your materials will likely be novices,
so make sure to explain all the relevant concepts.

## Title for a subsection
Section and subsection titles will be displayed in the tutorial index on the left side of
the page, so try to make them informative and concise!

# Hands-on Sections
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
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    https://zenodo.org/api/files/8e32cfe7-7f9f-4443-9c65-68242f601cc2/BK006939.2.genome.fasta
>    https://zenodo.org/api/files/8e32cfe7-7f9f-4443-9c65-68242f601cc2/BK006939.2.prot.fasta
>    https://zenodo.org/api/files/8e32cfe7-7f9f-4443-9c65-68242f601cc2/CM000925.1.genome.fasta
>    https://zenodo.org/api/files/8e32cfe7-7f9f-4443-9c65-68242f601cc2/CM000925.1.prot.fasta
>    https://zenodo.org/api/files/8e32cfe7-7f9f-4443-9c65-68242f601cc2/CM005043.2.genome.fasta
>    https://zenodo.org/api/files/8e32cfe7-7f9f-4443-9c65-68242f601cc2/CM005043.2.prot.fasta
>    https://zenodo.org/api/files/8e32cfe7-7f9f-4443-9c65-68242f601cc2/CM005299.1.genome.fasta
>    https://zenodo.org/api/files/8e32cfe7-7f9f-4443-9c65-68242f601cc2/CM005299.1.prot.fasta
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

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


## Sub-step with **Replace Text**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_line/1.1.2) %} with the following parameters:
>    - {% icon param-collection %} *"File to process"*: `output` (Input dataset collection)
>    - In *"Replacement"*:
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `(>[^ ]+).+`
>            - *"Replace with:"*: `\1`
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

## Sub-step with **RepeatMasker**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [RepeatMasker](toolshed.g2.bx.psu.edu/repos/bgruening/repeat_masker/repeatmasker_wrapper/4.1.2-p1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Genomic DNA"*: `outfile` (output of **Replace Text** {% icon tool %})
>    - *"Repeat library source"*: `DFam (curated only, bundled with RepeatMasker)`
>        - *"Select species name from a list?"*: `No`
>            - *"Repeat source species"*: `Saccharomyces cerevisiae`
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

## Sub-step with **Funannotate predict annotation**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Funannotate predict annotation](toolshed.g2.bx.psu.edu/repos/iuc/funannotate_predict/funannotate_predict/1.8.9+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Assembly to annotate"*: `output_masked_genome` (output of **RepeatMasker** {% icon tool %})
>    - In *"Organism"*:
>        - *"Name of the species to annotate"*: `Saccharomyces cerevisiae`
>        - *"Is it a fungus species?"*: `Yes`
>    - In *"Evidences"*:
>        - *"Select protein evidences"*: `Use UniProtKb/SwissProt (from selected Funannotate database)`
>    - In *"Busco"*:
>        - *"BUSCO models to align"*: `saccharomycetes`
>        - *"Initial Augustus species training set for BUSCO alignment"*: `saccharomyces`
>    - In *"Augustus settings (advanced)"*:
>        - *"Minimum number of models to train Augustus"*: `15`
>    - In *"EVM settings (advanced)"*:
>        - *"Split contigs into partitions for EVM processing?"*: `Yes`
>    - *"Which outputs should be generated"*: ``
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

## Sub-step with **Extract ORF**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Extract ORF](toolshed.g2.bx.psu.edu/repos/bgruening/glimmer_gbk_to_orf/glimmer_gbk_to_orf/3.02) %} with the following parameters:
>    - {% icon param-file %} *"gene bank file"*: `annot_gbk` (output of **Funannotate predict annotation** {% icon tool %})
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

## Sub-step with **Regex Find And Replace**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regex1/1.0.1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `aa_output` (output of **Extract ORF** {% icon tool %})
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `>([^ ]+).+`
>            - *"Replacement"*: `>#{input_name}_\1`
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

## Sub-step with **Collapse Collection**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Collapse Collection](toolshed.g2.bx.psu.edu/repos/nml/collapse_collections/collapse_dataset/4.2) %} with the following parameters:
>    - {% icon param-file %} *"Collection of files to collapse into single dataset"*: `out_file1` (output of **Regex Find And Replace** {% icon tool %})
>    - *"Prepend File name"*: `Yes`
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

## Sub-step with **Proteinortho**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Proteinortho](toolshed.g2.bx.psu.edu/repos/iuc/proteinortho/proteinortho/6.0.14+galaxy2.9.1) %} with the following parameters:
>    - {% icon param-file %} *"Select the input fasta files (>2)"*: `out_file1` (output of **Regex Find And Replace** {% icon tool %})
>    - *"Activate synteny feature (POFF)"*: `no`
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

## Sub-step with **Busco**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/4.1.4) %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `out_file1` (output of **Regex Find And Replace** {% icon tool %})
>    - *"Mode"*: `Proteome`
>    - *"Lineage"*: `Saccharomycetes`
>    - In *"Advanced Options"*:
>        - *"Augustus species model"*: `Use the default species for selected lineage`
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

## Sub-step with **Filter**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `proteinortho` (output of **Proteinortho** {% icon tool %})
>    - *"With following condition"*: `c1==4 and c2==4`
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

## Sub-step with **Proteinortho grab proteins**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Proteinortho grab proteins](toolshed.g2.bx.psu.edu/repos/iuc/proteinortho_grab_proteins/proteinortho_grab_proteins/6.0.14+galaxy2.9.1) %} with the following parameters:
>    - {% icon param-file %} *"Select the input fasta files"*: `output` (output of **Collapse Collection** {% icon tool %})
>    - *"Query type"*: `orthology-groups output file`
>        - {% icon param-file %} *"A orthology-groups file"*: `out_file1` (output of **Filter** {% icon tool %})
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

## Sub-step with **Regex Find And Replace**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regex1/1.0.1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `listproteinorthograbproteins` (output of **Proteinortho grab proteins** {% icon tool %})
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `(>[^_]+).+`
>            - *"Replacement"*: `\1`
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

## Sub-step with **ClustalW**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [ClustalW](toolshed.g2.bx.psu.edu/repos/devteam/clustalw/clustalw/2.1) %} with the following parameters:
>    - {% icon param-file %} *"FASTA file"*: `out_file1` (output of **Regex Find And Replace** {% icon tool %})
>    - *"Data type"*: `Protein sequences`
>    - *"Output alignment format"*: `FASTA format`
>    - *"Output complete alignment (or specify part to output)"*: `Complete alignment`
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

## Sub-step with **ClipKIT. Alignment trimming software for phylogenetics.**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [ClipKIT. Alignment trimming software for phylogenetics.](toolshed.g2.bx.psu.edu/repos/padge/clipkit/clipkit/0.1.0) %} with the following parameters:
>    - {% icon param-file %} *"Alignment file"*: `output` (output of **ClustalW** {% icon tool %})
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

## Sub-step with **PhyKit - Alignment-based functions**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [PhyKit - Alignment-based functions](toolshed.g2.bx.psu.edu/repos/padge/phykit/phykit_alignment_based/0.1.0) %} with the following parameters:
>    - *"Select tool for processing the alignment(s)"*: `Concatenate alignments.`
>        - {% icon param-file %} *"alignment list file. File should contain a single column list of alignment sequence files to concatenate into a single matrix."*: `trimmed_output` (output of **ClipKIT. Alignment trimming software for phylogenetics.** {% icon tool %})
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