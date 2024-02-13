---
layout: tutorial_hands_on

title: Secondary metabolite discovery
zenodo_link: https://zenodo.org/records/10652998
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

<!-- This is a comment. -->

**Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

> <agenda-title></agenda-title>
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

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    https://zenodo.org/api/records/10652998/files/MIBiG_compounds_3.0.sdf/content
>    https://zenodo.org/api/records/10652998/files/gbk2features.ipynb/content
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

> <details-title> More details about the theory </details-title>
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:


## Sub-step with **NCBI Accession Download**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NCBI Accession Download](toolshed.g2.bx.psu.edu/repos/iuc/ncbi_acc_download/ncbi_acc_download/0.2.8+galaxy0) %} with the following parameters:
>    - *"Select source for IDs"*: `Direct Entry`
>        - *"ID List"*: `{'id': 0, 'output_name': 'output'}`
>    - *"Molecule Type"*: `Nucleotide`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Molecule to fingerprint**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Molecule to fingerprint](toolshed.g2.bx.psu.edu/repos/bgruening/chemfp/ctb_chemfp_mol2fps/1.5) %} with the following parameters:
>    - {% icon param-file %} *"Molecule file"*: `output` (Input dataset)
>    - *"Type of fingerprint"*: `Open Babel FP2 fingerprints`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Antismash**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Antismash](toolshed.g2.bx.psu.edu/repos/bgruening/antismash/antismash/6.1.1+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Sequence file in GenBank,EMBL or FASTA format"*: `output` (output of **NCBI Accession Download** {% icon tool %})
>    - *"Taxonomic classification of input sequence"*: `Bacteria`
>    - *"Outputs"*: ``
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Collapse Collection**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Collapse Collection](toolshed.g2.bx.psu.edu/repos/nml/collapse_collections/collapse_dataset/5.1.0) %} with the following parameters:
>    - {% icon param-file %} *"Collection of files to collapse into single dataset"*: `genbank` (output of **Antismash** {% icon tool %})
>    - *"Prepend File name"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Interactive JupyTool and notebook**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Interactive JupyTool and notebook](interactive_tool_jupyter_notebook) %} with the following parameters:
>    - *"Do you already have a notebook?"*: `Load a previous notebook`
>        - {% icon param-file %} *"IPython Notebook"*: `output` (Input dataset)
>        - *"Execute notebook and return a new one."*: `Yes`
>    - In *"User inputs"*:
>        - {% icon param-repeat %} *"Insert User inputs"*
>            - *"Name for parameter"*: `dataset`
>            - *"Choose the input type"*: `Dataset`
>                - {% icon param-file %} *"Select value"*: `output` (output of **Collapse Collection** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Text reformatting**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Text reformatting](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_awk_tool/1.1.2) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `output_collection` (output of **Interactive JupyTool and notebook** {% icon tool %})
>    - *"AWK Program"*: `{print $8, $2-$5-$6}`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Remove beginning**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Remove beginning](Remove beginning1) %} with the following parameters:
>    - {% icon param-file %} *"from"*: `outfile` (output of **Text reformatting** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Remove duplicated molecules**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Remove duplicated molecules](toolshed.g2.bx.psu.edu/repos/bgruening/openbabel_remduplicates/openbabel_remDuplicates/3.1.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Molecular input file"*: `out_file1` (output of **Remove beginning** {% icon tool %})
>    - *"Select descriptor for molecule comparison"*: `Canonical SMILES`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Molecule to fingerprint**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Molecule to fingerprint](toolshed.g2.bx.psu.edu/repos/bgruening/chemfp/ctb_chemfp_mol2fps/1.5) %} with the following parameters:
>    - {% icon param-file %} *"Molecule file"*: `outfile` (output of **Remove duplicated molecules** {% icon tool %})
>    - *"Type of fingerprint"*: `Open Babel FP2 fingerprints`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Natural Product**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Natural Product](toolshed.g2.bx.psu.edu/repos/bgruening/natural_product_likeness/ctb_np-likeness-calculator/2.1) %} with the following parameters:
>    - {% icon param-file %} *"Molecule file"*: `outfile` (output of **Remove duplicated molecules** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Drug-likeness**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Drug-likeness](toolshed.g2.bx.psu.edu/repos/bgruening/qed/ctb_silicos_qed/2021.03.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Molecule data in SDF or SMILES format"*: `outfile` (output of **Remove duplicated molecules** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Similarity Search**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Similarity Search](toolshed.g2.bx.psu.edu/repos/bgruening/simsearch/ctb_simsearch/0.2) %} with the following parameters:
>    - *"Subject database/sequences"*: `Chemfp fingerprint file`
>        - *"Query Mode"*: `Query molecules are stores in a separate file`
>            - {% icon param-file %} *"Target molecules"*: `outfile` (output of **Molecule to fingerprint** {% icon tool %})
>        - *"select the k nearest neighbors"*: `1`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
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

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.