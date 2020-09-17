---
layout: tutorial_hands_on

title: Title of the tutorial
zenodo_link: ''
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


## Sub-step with **msconvert**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [msconvert](toolshed.g2.bx.psu.edu/repos/galaxyp/msconvert/msconvert/3.0.19052.0) %} with the following parameters:
>    - {% icon param-collection %} *"Input unrefined MS data"*: `output` (Input dataset collection)
>    - *"Do you agree to the vendor licenses?"*: `Yes`
>    - *"Output Type"*: `mgf`
>    - In *"Data Processing Filters"*:
>        - *"Apply peak picking?"*: `Yes`
>        - *"Apply m/z refinement with identification data?"*: `Yes`
>        - *"(Re-)calculate charge states?"*: `no`
>        - *"Filter m/z Window"*: `Yes`
>        - *"Filter out ETD precursor peaks?"*: `Yes`
>        - *"De-noise MS2 with moving window filter"*: `Yes`
>    - In *"Scan Inclusion/Exclusion Filters"*:
>        - *"Filter MS Levels"*: `Yes`
>    - In *"General Options"*:
>        - *"Sum adjacent scans"*: `Yes`
>        - *"Output multiple runs per file"*: `Yes`
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

## Sub-step with **Search GUI**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Search GUI](toolshed.g2.bx.psu.edu/repos/galaxyp/peptideshaker/search_gui/3.3.10.1) %} with the following parameters:
>    - {% icon param-file %} *"Protein Database"*: `output` (Input dataset)
>    - {% icon param-file %} *"Input Peak Lists (mgf)"*: `output` (output of **msconvert** {% icon tool %})
>    - In *"Search Engine Options"*:
>        - *"DB-Search Engines"*: ``
>    - In *"Protein Digestion Options"*:
>        - *"Digestion"*: `Select Enzymes`
>    - In *"Precursor Options"*:
>        - *"Fragment Tolerance"*: `0.2`
>        - *"Maximum Charge"*: `6`
>    - In *"Protein Modification Options"*:
>        - *"Fixed Modifications"*: ``
>        - *"Variable Modifications"*: ``
>    - In *"Andvanced Options"*:
>        - *"SearchGUI Options"*: `Default`
>        - *"X!Tandem Options"*: `Advanced`
>            - *"X!Tandem: Quick Acetyl"*: `Yes`
>            - *"X!Tandem: Quick Pyrolidone"*: `Yes`
>            - *"X!Tandem: Maximum Valid Expectation Value"*: `100.0`
>            - *"X!Tandem peptide model refinement"*: `Don't refine`
>        - *"OMSSA Options"*: `Default`
>        - *"MSGF Options"*: `Default`
>        - *"MS Amanda Options"*: `Default`
>        - *"TIDE Options"*: `Default`
>        - *"MyriMatch Options"*: `Default`
>        - *"Comet Options"*: `Default`
>        - *"DirectTag Options"*: `Default`
>        - *"Novor Options"*: `Default`
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

## Sub-step with **Peptide Shaker**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Peptide Shaker](toolshed.g2.bx.psu.edu/repos/galaxyp/peptideshaker/peptide_shaker/1.16.36.3) %} with the following parameters:
>    - {% icon param-file %} *"Compressed SearchGUI results"*: `searchgui_results` (output of **Search GUI** {% icon tool %})
>    - *"Specify Advanced PeptideShaker Processing Options"*: `Advanced Processing Options`
>        - *"The PTM probabilistic score to use for PTM localization"*: `A-score`
>    - *"Specify Advanced Filtering Options"*: `Advanced Filtering Options`
>        - *"Maximum Peptide Length"*: `60`
>    - *"Specify Contact Information for mzIdendML"*: `GalaxyP Project contact (Not suitable for PRIDE submission)`
>    - In *"Exporting options"*:
>        - *"Creates a mzIdentML file"*: `Yes`
>        - *"Compress results into a single zip file"*: `Yes`
>        - *"Reports to be generated"*: ``
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

## Sub-step with **Select**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Select](Grep1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `output_psm` (output of **Peptide Shaker** {% icon tool %})
>    - *"that"*: `NOT Matching`
>    - *"the pattern"*: `con_`
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

## Sub-step with **Select**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Select](Grep1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `output_peptides` (output of **Peptide Shaker** {% icon tool %})
>    - *"that"*: `NOT Matching`
>    - *"the pattern"*: `con_`
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

## Sub-step with **Replace Text**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_column/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `out_file1` (output of **Select** {% icon tool %})
>    - In *"Replacement"*:
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"in column"*: `c10`
>            - *"Find pattern"*: `.raw.mgf`
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

## Sub-step with **Cut**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c6`
>    - {% icon param-file %} *"From"*: `out_file1` (output of **Select** {% icon tool %})
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

## Sub-step with **FlashLFQ**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [FlashLFQ](toolshed.g2.bx.psu.edu/repos/galaxyp/flashlfq/flashlfq/1.0.3.0) %} with the following parameters:
>    - {% icon param-file %} *"identification file"*: `outfile` (output of **Replace Text** {% icon tool %})
>    - {% icon param-file %} *"spectrum files"*: `output` (output of **msconvert** {% icon tool %})
>    - *"match between runs"*: `Yes`
>    - *"Use experimnetal design for normalization or protein fold-change analysis"*: `Yes`
>        - {% icon param-file %} *"ExperimentalDesign.tsv"*: `output` (Input dataset)
>        - *"Perform Bayesian protein fold-change analysis"*: `Yes`
>            - *"control condition for Bayesian protein fold-change analysis"*: ``
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
>    - {% icon param-file %} *"Filter"*: `out_file1` (output of **Cut** {% icon tool %})
>    - *"With following condition"*: `len(c1)<=50`
>    - *"Number of header lines to skip"*: `1`
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
> 1. {% tool [Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regex1/1.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `quantifiedPeptides` (output of **FlashLFQ** {% icon tool %})
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Base Sequence`
>            - *"Replacement"*: `peptide`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Intensity_`
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

## Sub-step with **Unipept**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Unipept](toolshed.g2.bx.psu.edu/repos/galaxyp/unipept/unipept/4.3.0) %} with the following parameters:
>    - *"Unipept application"*: `peptinfo: Tryptic peptides and associated EC and GO terms and lowest common ancestor taxonomy`
>        - *"Equate isoleucine and leucine"*: `Yes`
>        - *"retrieve extra information"*: `Yes`
>        - *"group responses by GO namespace (biological process, molecular function, cellular component)"*: `Yes`
>        - *"allfields"*: `Yes`
>    - *"Peptides input format"*: `tabular`
>        - {% icon param-file %} *"Tabular Input Containing Peptide column"*: `out_file1` (output of **Filter** {% icon tool %})
>        - *"Select column with peptides"*: `c1`
>    - *"Choose outputs"*: ``
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

## Sub-step with **Unipept**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Unipept](toolshed.g2.bx.psu.edu/repos/galaxyp/unipept/unipept/4.0.0) %} with the following parameters:
>    - *"Unipept application"*: `pept2lca: lowest common ancestor`
>        - *"Equate isoleucine and leucine"*: `Yes`
>        - *"allfields"*: `Yes`
>    - *"Peptides input format"*: `tabular`
>        - {% icon param-file %} *"Tabular Input Containing Peptide column"*: `out_file1` (output of **Filter** {% icon tool %})
>        - *"Select column with peptides"*: `c1`
>    - *"Choose outputs"*: ``
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

## Sub-step with **Cut**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1,c3`
>    - {% icon param-file %} *"From"*: `output_ec_tsv` (output of **Unipept** {% icon tool %})
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

## Sub-step with **Query Tabular**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.0.0) %} with the following parameters:
>    - In *"Database Table"*:
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `output_go_tsv` (output of **Unipept** {% icon tool %})
>            - In *"Table Options"*:
>                - *"Specify Name for Table"*: `Goterm`
>                - *"Specify Column Names (comma-separated list)"*: `peptide,total_protein_count,go_term,protein_count,go_name,go_funct`
>    - *"SQL Query to generate tabular output"*: `SELECT Goterm.*
FROM Goterm 
WHERE ((1.0*Goterm.protein_count)/(1.0*Goterm.total_protein_count)) >= 0.05


`
>    - *"include query result column headers"*: `Yes`
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

## Sub-step with **Replace Text**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_line/1.1.2) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `output_tsv` (output of **Unipept** {% icon tool %})
>    - In *"Replacement"*:
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `#peptide`
>            - *"Replace with:"*: `peptide`
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

## Sub-step with **Query Tabular**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.0.0) %} with the following parameters:
>    - In *"Database Table"*:
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `out_file1` (output of **Cut** {% icon tool %})
>            - In *"Table Options"*:
>                - *"Specify Name for Table"*: `ec`
>                - *"Specify Column Names (comma-separated list)"*: `peptide,go_ec`
>    - *"SQL Query to generate tabular output"*: `SELECT *
FROM ec`
>    - *"include query result column headers"*: `Yes`
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
>    - {% icon param-file %} *"Filter"*: `output` (output of **Query Tabular** {% icon tool %})
>    - *"With following condition"*: `c5=='biological process'`
>    - *"Number of header lines to skip"*: `1`
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
>    - {% icon param-file %} *"Filter"*: `output` (output of **Query Tabular** {% icon tool %})
>    - *"With following condition"*: `c5=='cellular component'`
>    - *"Number of header lines to skip"*: `1`
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
>    - {% icon param-file %} *"Filter"*: `output` (output of **Query Tabular** {% icon tool %})
>    - *"With following condition"*: `c5=='molecular function'`
>    - *"Number of header lines to skip"*: `1`
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