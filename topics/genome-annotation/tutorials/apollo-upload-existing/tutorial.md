---
layout: tutorial_hands_on
topic_name: genome-annotation
tutorial_name: apollo-upload-existing
title: “Upload Existing Genome into Apollo”
---

# Introduction
{:.no_toc}

<!-- These are notes for myself. This is ONLY in markdown. -->

GenBank allows previously uploaded genomes to be downloaded onto your local device. These genomes can then be uploaded into Apollo for annotation and editing. Here, the genome of phage P22 will be used as an example.

> ### Agenda
>
> In this tutorial, we will explore the path of a genome file:
>
> 1. From GenBank to Local Device
> 
> 2. From Local Device to Apollo
>
> 3. Preparing for Annotation in Apollo
>
{: .agenda}

# From GenBank to Local Device

The NIH genetic sequence database, GenBank, is an open access collection of annotated, publicly available nucleotide sequences and protein translations.

<!-- Here, it may be helpful to link the Genome Annotation tutorial within Galaxy GitHub. -—>

<!--
{% icon hands_on %} will render the hands_on icon as specified in
_config.yml in the root of this repository.
-->

<!-- Better separation of what belongs in steps vs. what belongs in tip/comment box. -->

1. Access [the NCBI website.](https://www.ncbi.nlm.nih.gov “NCBI Homepage”)

2. Using specific key words, search for your desired genome in the search box at the top of the page.

3. Look through the results for hits that appear to align with your search. Most likely, there are multiple hits for the same organism. It is important that the file contains the whole genome, and is the most recently updated version of that genome.

> ### {% icon tip %}Tip: Searching
> * Adjust the drop-down menu to the left of the search bar to ‘Nucleotide.’ This will narrow your results to genome, gene, and transcript sequence data from databases including GenBank.
> * If your search yields many seemingly unrelated results, consider using the ‘Advanced Search’ option to further refine results.
{: .tip}

> ### {% icon hands_on %} Hands-on: Genome download
> * When you’ve found the optimal genome file, click on the ‘Send to:’ drop down menu towards the top right of the page.
> * Once the download parameters have been set, click ‘Create file,’ and the genome file will be downloaded to your computer.
>  ![](../../images/apollo-upload-existing-tutorial-screenshots/1 Genbank File Download Settings.png)
>
>    > #### {% icon details %} More details on the ....
>    >
>    > Add more details in Markdown. By default the box is collapsed. And is expanded when clicked
>    >
>    {: .details}
{: .hands_on}

# From Local Device to Apollo

Now that the genome has been downloaded in the proper format, it is ready to be transferred to Apollo.

> ### {% icon hands_on %} Hands-on: Genome loading into Apollo
>
> 1. Open [Galaxy](https://cpt.tamu.edu/galaxy-pub “Galaxy”) and log in.
> 2. On the right side of the web page is the History column. In the top right corner of the History column are three symbols.
>    > * If you already have an active history open, click on the gear symbol. Click on ‘Create new history.’
>    > * If you do not have an active history, move on to step 3.
> ![](../../images/apollo-upload-existing-tutorial-screenshots/6 New History.png)
<!-- Weird formatting of pictures -->
> 3. Once you have an empty, new history ready, click on the ‘load your own data’ option in the blue box in the empty history. An empty, white box should appear on your screen. ![](../../images/apollo-upload-existing-tutorial-screenshots/2 New History Created.png)![](../../images/apollo-upload-existing-tutorial-screenshots/4 Uploading Genome Drag n Drop.png)
> 4. You can drag the downloaded genome file into the box, or search files on your local device by clicking ‘Choose local file’ and locating the downloaded genome on your local device.
> 5. Click ‘Start’ to begin loading the GenBank genome file into Apollo. Upon completion of this loading, there should be a green box now in your History column that says, “1: sequence.gb.”
>    > * If you click on the eye symbol of this step, it will open up the contents of that step. Because it is the full GenBank file, it should resemble what was on the NCBI webpage for that genome when you downloaded it.
{: .hands_on}

# Preparing Genome for Annotation in Apollo

You have successfully downloaded a genome from NCBI GenBank, and placed said genome in Apollo. The next step is loading the full GenBank file into Apollo so that the genome may be annotated/edited. To do this, you will need to import and run a workflow.

1. At the top of the Galaxy/Apollo page, click on the “Shared Data” drop-down menu. Select “workflows,” and you will be brought to a new page containing workflows that execute different functions within Apollo.

2. Find “Load GenBank (.gb, .gbk) into Apollo (v#)” where ‘#’ is the number indicating the most recent version of this workflow. Click on the drop-down menu on that workflow, and select “Import.”

3. Return to 

COMMIT TEST


> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Step1
> 2. **My Tool** {% icon tool %} with the following parameters
>   - *"param1"*: the file `myfile`
>   - *"param2"*: `42`
>   - *"param3"*: `Yes`
>
> 3. **My Tool** {% icon tool %} with the following parameters
>   - {% icon param-text %} *"My text parameter"*: `my value`
>   - {% icon param-file %} *"My input file"*: `my file`
>   - {% icon param-files %} *"My multiple file input or collection"*: `my collection`
>   - {% icon param-select %} *"My select menu"*: `my choice`
>   - {% icon param-check %} *"My check box"*: `yes`
>
>    > ### {% icon question %} Questions
>    >
>    > 1. Question1?
>    > 2. Question2?
>    >
>    >    > ### {% icon solution %} Solution
>    >    >
>    >    > 1. Answer for question1
>    >    > 2. Answer for question2
>    >    >
>    >    {: .solution}
>    >
>    {: .question}
>
> 3. Step3
{: .hands_on}

# Conclusion
{:.no_toc}

Conclusion about the technical key points. And then relation between the techniques and the biological question to end with a global view.
