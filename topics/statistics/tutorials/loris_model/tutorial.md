---
layout: tutorial_hands_on
level: Intermediate
title: Building the LORIS LLR6 Model Using Galaxy ML Tools
zenodo_link: https://zenodo.org/records/13885908
questions:
- How can I build the LORIS LLR6 model published by Chang et al., 2024?
- Which tools can I use in Galaxy to obtain a logistic regression model with my desired hyperparameters?
- How can I evaluate the model to confirm its performance?
objectives:
- Build the LORIS LLR6 model using machine learning tools available in Galaxy for reproducibility.
- Use the PyCaret tool to generate a new model using the LORIS dataset published by Chang et al., 2024.
- Evaluate the created models for reproducibility and robustness by comparing them with the original LORIS LLR6 model.
time_estimation: 2H
key_points:
- Use Galaxy tools to build the identical LORIS LLR6 model published by Chang et al., 2024.
- Use the PyCaret tool to train a new model and confirm the robustness of the study published by Chang et al., 2024.
contributors:
- paulocilasjr 
- qchiujunhao 
- jgoecks
tags:
- LORIS Score Model
- Machine Learning
- PyCaret
- Galaxy-ML tools
---

> <comment-title>Building the LORIS LLR6 Model</comment-title>
>
> The PyCaret Model Comparison tool described in this tutorial is only available at: 
> [Cancer-Galaxy](https://cancer.usegalaxy.org)
>
> Galaxy-ML tools > PyCaret Model Comparison 
>
{:  .comment}

Using the MNIST image dataset of handwritten digits as input, we will build an image recognition model with the Galaxy-Ludwig tool. The goal is to classify the handwritten digit in each image.

To accomplish this, three steps are needed: (i) upload Ludwig files and image files to Galaxy (ii) Set up and running the Ludwig experiment function on Galaxy, and (iii) Evaluate the image classification model. As a bonus step, we'll also explore (iv) improving the model's classification performance (Figure 1).

![schema of the whole process of training model and test.](../../images/galaxy-ludwig/explain_model_schema.png "Overview of the steps process to obtain the handwritten classification model and testing it.")

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Introduction

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
>    https://zenodo.org/api/records/13885908/files/Chowell_test_No_Response.tsv/content
>    https://zenodo.org/api/records/13885908/files/Chowell_train_Response.tsv/content
>    https://zenodo.org/api/records/13885908/files/Chowell_test_Response.tsv/content
>    https://zenodo.org/api/records/13885908/files/MSK1_Response.tsv/content
>    https://zenodo.org/api/records/13885908/files/MSK1_No_Response.tsv/content
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


## Sub-step with **Hyperparameter Search**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Hyperparameter Search](toolshed.g2.bx.psu.edu/repos/bgruening/sklearn_searchcv/sklearn_searchcv/1.0.11.0) %} with the following parameters:
>    - *"Select a hyperparameter search algorithm"*: `GridSearchCV - Exhaustive search over specified parameter values for an estimator `
>    - {% icon param-file %} *"Choose the dataset containing pipeline/estimator object"*: `outfile` (output of **Pipeline Builder** {% icon tool %})
>    - In *"Advanced Options for SearchCV"*:
>        - *"Select the primary metric (scoring):"*: `default with estimator`
>        - *"Select the cv splitter:"*: `default splitter`
>    - *"Select input type:"*: `tabular data`
>        - {% icon param-file %} *"Training samples dataset:"*: `output` (Input dataset)
>        - *"Does the dataset contain header:"*: `Yes`
>        - *"Choose how to select data by column:"*: `All columns EXCLUDING some by column index number(s)`
>            - *"Select target column(s):"*: `c['22']`
>        - {% icon param-file %} *"Dataset containing class labels or target values:"*: `output` (Input dataset)
>        - *"Does the dataset contain header:"*: `Yes`
>        - *"Choose how to select data by column:"*: `Select columns by column index number(s)`
>            - *"Select target column(s):"*: `c22`
>    - *"Whether to hold a portion of samples for test exclusively, nested CV?"*: `Nope`
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
>    - {% icon param-file %} *"from"*: `output` (Input dataset)
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
>    - {% icon param-file %} *"from"*: `output` (Input dataset)
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

## Sub-step with **Ensemble methods**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Ensemble methods](toolshed.g2.bx.psu.edu/repos/bgruening/sklearn_ensemble/sklearn_ensemble/1.0.11.0) %} with the following parameters:
>    - *"Select a Classification Task"*: `Load a model and predict`
>        - {% icon param-file %} *"Models"*: `outfile_object` (output of **Hyperparameter Search** {% icon tool %})
>        - {% icon param-file %} *"Data (tabular)"*: `output` (Input dataset)
>        - *"Does the dataset contain header:"*: `Yes`
>        - *"Select the type of prediction"*: `Predict class labels`
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

## Sub-step with **Ensemble methods**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Ensemble methods](toolshed.g2.bx.psu.edu/repos/bgruening/sklearn_ensemble/sklearn_ensemble/1.0.11.0) %} with the following parameters:
>    - *"Select a Classification Task"*: `Load a model and predict`
>        - {% icon param-file %} *"Models"*: `outfile_object` (output of **Hyperparameter Search** {% icon tool %})
>        - {% icon param-file %} *"Data (tabular)"*: `output` (Input dataset)
>        - *"Does the dataset contain header:"*: `Yes`
>        - *"Select the type of prediction"*: `Predict class labels`
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
>    - {% icon param-file %} *"from"*: `outfile_predict` (output of **Ensemble methods** {% icon tool %})
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
>    - {% icon param-file %} *"from"*: `outfile_predict` (output of **Ensemble methods** {% icon tool %})
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

## Sub-step with **Plot confusion matrix, precision, recall and ROC and AUC curves**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Plot confusion matrix, precision, recall and ROC and AUC curves](toolshed.g2.bx.psu.edu/repos/bgruening/plotly_ml_performance_plots/plotly_ml_performance_plots/0.4) %} with the following parameters:
>    - {% icon param-file %} *"Select input data file :"*: `out_file1` (output of **Remove beginning** {% icon tool %})
>    - {% icon param-file %} *"Select predicted data file :"*: `out_file1` (output of **Remove beginning** {% icon tool %})
>    - {% icon param-file %} *"Select trained model :"*: `outfile_object` (output of **Hyperparameter Search** {% icon tool %})
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

## Sub-step with **Plot confusion matrix, precision, recall and ROC and AUC curves**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Plot confusion matrix, precision, recall and ROC and AUC curves](toolshed.g2.bx.psu.edu/repos/bgruening/plotly_ml_performance_plots/plotly_ml_performance_plots/0.4) %} with the following parameters:
>    - {% icon param-file %} *"Select input data file :"*: `out_file1` (output of **Remove beginning** {% icon tool %})
>    - {% icon param-file %} *"Select predicted data file :"*: `out_file1` (output of **Remove beginning** {% icon tool %})
>    - {% icon param-file %} *"Select trained model :"*: `outfile_object` (output of **Hyperparameter Search** {% icon tool %})
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
