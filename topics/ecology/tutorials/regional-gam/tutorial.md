---
layout: tutorial_hands_on
topic_name: ecology
tutorial_name: regional-gam
---

# Introduction
{:.no_toc}

<!-- This is a comment. -->

General introduction about the topic and then an introduction of the 
tutorial (the questions and the objectives). It is nice also to have a 
scheme to sum up the pipeline used during the tutorial. The idea is to 
give to trainees insight into the content of the tutorial and the (theoretical 
and technical) key concepts they will learn.

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

Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word `TODO`, there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Import the following files from [Zenodo](https://zenodo.org/record/1321885) or from a data
>    library named `TODO` if available (ask your instructor)
>
>    ```
>    TODO: add the files by the ones on Zenodo here (if not added)
>    TODO: remove the useless files (if added)
>    TODO: so that they can easily be copy-pasted into Galaxy's upload dialog
>    https://zenodo.org/api/files/51a1b5db-ff05-4cda-83d4-3b46682f921f/gatekeeper_CM%20.RData
>    https://zenodo.org/api/files/51a1b5db-ff05-4cda-83d4-3b46682f921f/regional%20GAM%20data.csv
>    ```
>
>    > ### {% icon tip %} Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field
>    > * Press **Start**
>    >
>    > By default, Galaxy uses the url as the name, so please rename them to something more pleasing.
>    {: .tip}
>
>    > ### {% icon tip %} Tip: Importing data from a data library
>    >
>    > * Go into "Shared data" (top panel) then "Data libraries"
>    > * Click on "Training data" and then "Analyses of metagenomics data"
>    > * Select interesting file
>    > * Click on "Import selected datasets into history"
>    > * Import in a new history
>    {: .tip}
>
{: .hands_on}


# Title of the section usually corresponding to a big step

Description of the step: some background and some theory. Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part. 

<-- Consider adding a detail box to expand the theory -->

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
> 
{: .details}

> ### {% icon hands_on %} Hands-on: TODO: task description
>
> 1. **CSV to Tabular** {% icon tool %} with the following parameters:
>   - {% icon param-file %} *"CSV file"*: `output` (output of **Input dataset** {% icon tool %})
>   - *"Separator"*: `","`
>   - *"Header in file"*: `Yes`
>
>   TODO: check parameter descriptions
>   TODO: some of these parameters may be the default values and can be removed
>         unless they have some didactic value.
>
>   <-- Consider adding a comment or tip box -->
>
>   > ### {% icon comment %}} Comment
>   >
>   > A comment about the tool or something else. This box can also be in the main text
>   {: .comment}
>
{: .hands_on}

<-- Consider adding a question to test the learners understanding of the previous exercise -->

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


# Title of the section usually corresponding to a big step

Description of the step: some background and some theory. Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part. 

<-- Consider adding a detail box to expand the theory -->

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
> 
{: .details}

> ### {% icon hands_on %} Hands-on: TODO: task description
>
> 1. **Compter** {% icon tool %} with the following parameters:
>   - {% icon param-file %} *"Sur le jeu de données"*: `output` (output of **CSV to Tabular** {% icon tool %})
>   - *"Compter les occurrences des valeurs présentes dans la(les) colonne(s)"*: `c["1"]`
>   - *"Délimité par"*: ``
>   - *"Comment les résultats doivent t'ils être triés ?"*: ``
>
>   TODO: check parameter descriptions
>   TODO: some of these parameters may be the default values and can be removed
>         unless they have some didactic value.
>
>   <-- Consider adding a comment or tip box -->
>
>   > ### {% icon comment %}} Comment
>   >
>   > A comment about the tool or something else. This box can also be in the main text
>   {: .comment}
>
{: .hands_on}

<-- Consider adding a question to test the learners understanding of the previous exercise -->

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


# Title of the section usually corresponding to a big step

Description of the step: some background and some theory. Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part. 

<-- Consider adding a detail box to expand the theory -->

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
> 
{: .details}

> ### {% icon hands_on %} Hands-on: TODO: task description
>
> 1. **Compter** {% icon tool %} with the following parameters:
>   - {% icon param-file %} *"Sur le jeu de données"*: `output` (output of **CSV to Tabular** {% icon tool %})
>   - *"Compter les occurrences des valeurs présentes dans la(les) colonne(s)"*: `c["3"]`
>   - *"Délimité par"*: ``
>   - *"Comment les résultats doivent t'ils être triés ?"*: ``
>
>   TODO: check parameter descriptions
>   TODO: some of these parameters may be the default values and can be removed
>         unless they have some didactic value.
>
>   <-- Consider adding a comment or tip box -->
>
>   > ### {% icon comment %}} Comment
>   >
>   > A comment about the tool or something else. This box can also be in the main text
>   {: .comment}
>
{: .hands_on}

<-- Consider adding a question to test the learners understanding of the previous exercise -->

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


# Title of the section usually corresponding to a big step

Description of the step: some background and some theory. Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part. 

<-- Consider adding a detail box to expand the theory -->

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
> 
{: .details}

> ### {% icon hands_on %} Hands-on: TODO: task description
>
> 1. **Compter** {% icon tool %} with the following parameters:
>   - {% icon param-file %} *"Sur le jeu de données"*: `output` (output of **CSV to Tabular** {% icon tool %})
>   - *"Compter les occurrences des valeurs présentes dans la(les) colonne(s)"*: `c["2"]`
>   - *"Délimité par"*: ``
>   - *"Comment les résultats doivent t'ils être triés ?"*: ``
>
>   TODO: check parameter descriptions
>   TODO: some of these parameters may be the default values and can be removed
>         unless they have some didactic value.
>
>   <-- Consider adding a comment or tip box -->
>
>   > ### {% icon comment %}} Comment
>   >
>   > A comment about the tool or something else. This box can also be in the main text
>   {: .comment}
>
{: .hands_on}

<-- Consider adding a question to test the learners understanding of the previous exercise -->

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


# Title of the section usually corresponding to a big step

Description of the step: some background and some theory. Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part. 

<-- Consider adding a detail box to expand the theory -->

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
> 
{: .details}

> ### {% icon hands_on %} Hands-on: TODO: task description
>
> 1. **Trouver et Remplacer des patterns dans des colonnes** {% icon tool %} with the following parameters:
>   - {% icon param-file %} *"Selectionner les cellules à partir de"*: `output` (output of **CSV to Tabular** {% icon tool %})
>   - *"la colonne"*: `c"2"`
>   - In *"Check"*:
>      - In *"1: Check"*:
>         - *"Trouver l'expression suivante"*: `(\.[0-9]+)`
>         - *"Remplacement"*: ``
>
>   TODO: check parameter descriptions
>   TODO: some of these parameters may be the default values and can be removed
>         unless they have some didactic value.
>
>   <-- Consider adding a comment or tip box -->
>
>   > ### {% icon comment %}} Comment
>   >
>   > A comment about the tool or something else. This box can also be in the main text
>   {: .comment}
>
{: .hands_on}

<-- Consider adding a question to test the learners understanding of the previous exercise -->

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


# Title of the section usually corresponding to a big step

Description of the step: some background and some theory. Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part. 

<-- Consider adding a detail box to expand the theory -->

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
> 
{: .details}

> ### {% icon hands_on %} Hands-on: TODO: task description
>
> 1. **Tabular to CSV** {% icon tool %} with the following parameters:
>   - {% icon param-file %} *"tabular file"*: `out_file1` (output of **Trouver et Remplacer des patterns dans des colonnes** {% icon tool %})
>   - *"output csv Separator"*: `","`
>   - *"Header in file"*: `Yes`
>
>   TODO: check parameter descriptions
>   TODO: some of these parameters may be the default values and can be removed
>         unless they have some didactic value.
>
>   <-- Consider adding a comment or tip box -->
>
>   > ### {% icon comment %}} Comment
>   >
>   > A comment about the tool or something else. This box can also be in the main text
>   {: .comment}
>
{: .hands_on}

<-- Consider adding a question to test the learners understanding of the previous exercise -->

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


# Title of the section usually corresponding to a big step

Description of the step: some background and some theory. Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part. 

<-- Consider adding a detail box to expand the theory -->

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
> 
{: .details}

> ### {% icon hands_on %} Hands-on: TODO: task description
>
> 1. **Flight curve** {% icon tool %} with the following parameters:
>   - {% icon param-file %} *"Fichier de comptage"*: `output` (output of **Tabular to CSV** {% icon tool %})
>
>   TODO: check parameter descriptions
>   TODO: some of these parameters may be the default values and can be removed
>         unless they have some didactic value.
>
>   <-- Consider adding a comment or tip box -->
>
>   > ### {% icon comment %}} Comment
>   >
>   > A comment about the tool or something else. This box can also be in the main text
>   {: .comment}
>
{: .hands_on}

<-- Consider adding a question to test the learners understanding of the previous exercise -->

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



# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.