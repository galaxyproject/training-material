---
layout: tutorial_hands_on

title: EncyclopeDIA
zenodo_link: ''
questions:
- How to perform quantitative analysis using DIA data with the help of EncyclopeDIA?
- How to perform quantitation with or without Chromatogram Libraries?
objectives:
- Understand the difference between DDA and DIA methods
- Performing quantitative analysis using DIA data
- Understand the purpose and benefits of using a Chromatogram Library for detection of peptides
- You can use Bloom's Taxonomy to write effective learning objectives
time_estimation: 3H
key_points:
- With SearchToLib, Chromatogram Libraries can be created
- Learning conversion of file types using MSConvert
contributors:
- emmaleith
- subinamehta
- jj-umn
- pratikdjagtap
- timothygriffin


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
{% raw %} `{% cite FernandezCosta2020 %}`{% endraw %}.
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

# Conversion of file types

msconvert is the first tool in this EncyclopeDIA workflow and before analysis of such DIA data may begin, the data file must be converted to the correct file type (mzML) from the raw data. As msconvert exists on the Galaxy platform, conversion of files to the necessary type is straightforward, and can be incorporated into the workflow itself as opposed to a separate precursor. Both the Gas Phase Fractionated DIA data and the Experimental DIA data are run through msconvert to convert to the necessary file type (mzML) for the following steps, creation of the chromatogram library and analysis of DIA data through EncyclopeDIA.

![Alternative text](../../images/image_name "Legend of the image")


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
> 1. {% tool [msconvert](toolshed.g2.bx.psu.edu/repos/galaxyp/msconvert/msconvert/3.0.19052.1) %} with the following parameters:
>    - {% icon param-collection %} *"Input unrefined MS data"*: `output` (Input dataset collection)
>    - *"Do you agree to the vendor licenses?"*: `Yes`
>    - *"Output Type"*: `mzML`
>    - In *"Data Processing Filters"*:
>        - *"Apply peak picking?"*: `Yes`
>        - *"Apply m/z refinement with identification data?"*: `Yes`
>        - *"(Re-)calculate charge states?"*: `no`
>        - *"Filter m/z Window"*: `Yes`
>        - *"Filter out ETD precursor peaks?"*: `Yes`
>        - *"De-noise MS2 with moving window filter"*: `Yes`
>        - *"Demultiplex overlapping or MSX spectra"*: `Yes`
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
> 1. WHy is converting to mzml necessary?
> 2. Can you use any other tool for conversion?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **msconvert**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [msconvert](toolshed.g2.bx.psu.edu/repos/galaxyp/msconvert/msconvert/3.0.19052.1) %} with the following parameters:
>    - {% icon param-collection %} *"Input unrefined MS data"*: `output` (Input dataset collection)
>    - *"Do you agree to the vendor licenses?"*: `Yes`
>    - *"Output Type"*: `mzML`
>    - In *"Data Processing Filters"*:
>        - *"Apply peak picking?"*: `Yes`
>        - *"Apply m/z refinement with identification data?"*: `Yes`
>        - *"(Re-)calculate charge states?"*: `no`
>        - *"Filter m/z Window"*: `Yes`
>        - *"Filter out ETD precursor peaks?"*: `Yes`
>        - *"De-noise MS2 with moving window filter"*: `Yes`
>        - *"Demultiplex overlapping or MSX spectra"*: `Yes`
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

# Chromatogram Library Generation

Write about Chromatogram Libraries and add few images
![Alternative text](../../images/image_name "Legend of the image")
![Alternative text](../../images/image_name "Legend of the image")

SearchToLib is the tool responsible for the generation of the Chromatogram Library in this EncyclopeDIA workflow. A library is generated using the spectral files converted previously, a background protein database in FASTA format, as well as a .dlib spectral library. Outputs from this tool include the Chromatogram Library in an .elib format, as well as a text log file.

## **SearchToLib**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [SearchToLib](toolshed.g2.bx.psu.edu/repos/galaxyp/encyclopedia_searchtolib/encyclopedia_searchtolib/0.9.5.0) %} with the following parameters:
>    - {% icon param-file %} *"Spectrum files in  mzML format"*: `output` (output of **msconvert** {% icon tool %})
>    - {% icon param-file %} *"Library: Chromatagram .ELIB or Spectrum .DLIB"*: `output` (Input dataset)
>    - {% icon param-file %} *"Background proteome protein fasta database"*: `output` (Input dataset)
>    - In *"Parameter Settings"*:
>        - *"Set Acquisition Options"*: `No - use default options`
>        - *"Set Tolerance Options"*: `No - use default options`
>        - *"Set Percolator Options"*: `No - use default options`
>        - *"Set Peak Options"*: `No - use default options`
>        - *"Set Window Options"*: `No - use default options`
>        - *"Set Modifications Options"*: `No - use default options`
>        - *"Set Search Options"*: `No - use default options`
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
> 1. What happens if you already have a library?
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

# Without Chromatogram library

Write about Walnut if Chromatogram Libraries are absent and add few images
![Alternative text](../../images/image_name "Legend of the image")
![Alternative text](../../images/image_name "Legend of the image")

## Sub-step with **Walnut**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [EncyclopeDIA Quantify](toolshed.g2.bx.psu.edu/repos/galaxyp/encyclopedia_quantify/encyclopedia_quantify/0.9.5.0) %} with the following parameters:
>    - {% icon param-file %} *"Spectrum files in  mzML format"*: `output` (output of **msconvert** {% icon tool %})
>    - {% icon param-file %} *"Library: Chromatagram .ELIB or Spectrum .DLIB"*: `elib` (output of **SearchToLib** {% icon tool %})
>    - {% icon param-file %} *"Background proteome protein fasta database"*: `output` (Input dataset)
>    - In *"Parameter Settings"*:
>        - *"Set Acquisition Options"*: `No - use default options`
>        - *"Set Tolerance Options"*: `No - use default options`
>        - *"Set Percolator Options"*: `No - use default options`
>        - *"Set Peak Options"*: `No - use default options`
>        - *"Set Window Options"*: `No - use default options`
>        - *"Set Modifications Options"*: `No - use default options`
>        - *"Set Search Options"*: `No - use default options`
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

# Analysis of DIA data through EncyclopeDIA

Encyclopedia is the tool used for DIA data analysis through searching peptides against the generated Chromatogram Library. Utilizing the generated Chromatogram library, as well as the experimental DIA data (mzML format), and the background protein database used previously, EncyclopeDIA searches the experimental DIA data against these libraries. Generated are a log .txt file and two quantitation outputs, one for protein quantitation and one for peptide quantitation. 

![Alternative text](../../images/image_name "Legend of the image")

## Sub-step with **EncyclopeDIA Quantify**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [EncyclopeDIA Quantify](toolshed.g2.bx.psu.edu/repos/galaxyp/encyclopedia_quantify/encyclopedia_quantify/0.9.5.0) %} with the following parameters:
>    - {% icon param-file %} *"Spectrum files in  mzML format"*: `output` (output of **msconvert** {% icon tool %})
>    - {% icon param-file %} *"Library: Chromatagram .ELIB or Spectrum .DLIB"*: `elib` (output of **SearchToLib** {% icon tool %})
>    - {% icon param-file %} *"Background proteome protein fasta database"*: `output` (Input dataset)
>    - In *"Parameter Settings"*:
>        - *"Set Acquisition Options"*: `No - use default options`
>        - *"Set Tolerance Options"*: `No - use default options`
>        - *"Set Percolator Options"*: `No - use default options`
>        - *"Set Peak Options"*: `No - use default options`
>        - *"Set Window Options"*: `No - use default options`
>        - *"Set Modifications Options"*: `No - use default options`
>        - *"Set Search Options"*: `No - use default options`
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
