---
layout: tutorial_hands_on
title: Best Practices for Citing Galaxy
tags:
  - fair
  - how-to-cite-galaxy
questions:
  - How to cite Galaxy in Publications?
objectives:
  - Apply the Galaxy citation guidelines to your own work.
  - Prepare your methodology section based on your Galaxy Histories.
time_estimation: 20m
key_points:
  - The Galaxy Project primary publication should always be cited.
  - When analysing data using a public Galaxy server, include an acknowledgement for the Galaxy server(s) used.
  - Analysis methods can be made transparent and reproducible by publishing Galaxy Histories and Workflows.
  - Public data used in an analysis should always be cited.
  - A list of tool citations can be extracted directly from your Galaxy History.
contributions:
  authorship:
    - tflowers15
    - GarethPrice-Aus
    - mirelaminkova
    - MarisaJL
  funding:
    - unimelb
    - melbournebioinformatics
    - AustralianBioCommons
    - surf

---

# Overview

This is a short primer on how to cite Galaxy in your publication and also how to prepare your Methodology based on your Galaxy history or histories. Following these steps will make it easy to cite Galaxy, share your methods, and adhere to best practices for FAIR and open science.


> <agenda-title></agenda-title>
> In this tutorial we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Citing Galaxy
Citing Galaxy comes in two main actions: 

1. Citing the global **Galaxy Project** that provides the core code that all Galaxy instances operate on. 
2. Citing the **instance of Galaxy**, you did your work on.

## Citing Galaxy Project
Please cite Galaxy Project’s primary publication: 

- ***The Galaxy platform for accessible, reproducible, and collaborative data analyses: 2024 update***
  - Nucleic Acids Research, Volume 52, Issue W1, 5 July 2024, Pages W83–W94, [https://doi.org/10.1093/nar/gkae410](https://doi.org/10.1093/nar/gkae410)

This is also found on the Galaxy Project Hub with more specific examples of how to cite individual aspects of Galaxy Project - [https://galaxyproject.org/citing-galaxy/](https://galaxyproject.org/citing-galaxy/)

> <comment-title></comment-title>
> The main Galaxy Project publication will be included in the list if you follow the instructions below to extract tool citations from your History.
>
{: .comment}

## Citing a Galaxy Service

### Service Descriptions
Each public Galaxy service benefits from being able to track the amazing work done by researchers using its service. Provided below are the DOIs for the Galaxy Australia, France and Main `usegalaxy.*` services. The linked papers are available on Zenodo in the Galaxy Project Zenodo Community at [https://zenodo.org/communities/galaxy/](https://zenodo.org/communities/galaxy/):

- Galaxy Australia - [https://doi.org/10.5281/zenodo.17298823](https://doi.org/10.5281/zenodo.17298823)
- Galaxy France - [https://doi.org/10.5281/zenodo.17298858](https://doi.org/10.5281/zenodo.17298858)
- Galaxy Main - [https://doi.org/10.5281/zenodo.17298864](https://doi.org/10.5281/zenodo.17298864)

We recommend including the relevant service name and DOI link in your Methods.


### Acknowledging a Service
Each public Galaxy service is provided through a complex mixture of funding, governance and infrastructure providers. As such, each service provides a cut/paste **Acknowledgement** statement that supports their particular configuration and this can be found on each service. If your publication supports Acknowledgements, we strongly recommend you include the service-recommended **Acknowledgements**.

**Acknowledgement** statements for the four main `usegalaxy.*` services can be found using the following links:

- Galaxy Australia - [https://site.usegalaxy.org.au/about#acknowledgement-statement](https://site.usegalaxy.org.au/about#acknowledgement-statement)
- Galaxy Europe - [https://usegalaxy-eu.github.io/about](https://usegalaxy-eu.github.io/about)
- Galaxy France - [https://www.france-bioinformatique.fr/en/services/acknowledgement/](https://www.france-bioinformatique.fr/en/services/acknowledgement/)
- Galaxy Main - [https://usegalaxy.org/#about](https://usegalaxy.org/#about)

# Writing a Methods section based on your Galaxy work

## Making a History public

We provide the following recommendations and process for making your History persistent and your options for public archiving of your History. Sharing your History will make it easier for others to understand or replicate your analysis. You can also add **annotations** and **tags** to ensure your history aligns to the FAIR principles.

### Preparing a History

We recommend cleaning up your History when preparing to make your History public.

> <hands-on-title>Clean up Unneeded Datasets</hands-on-title>
>
> 1. `Delete` and `Delete (permanently)` (purge) any failed jobs.
> 2. `Delete` and `Delete (permanently)` (purge) any datasets/job outputs that were used for testing and are not required for the final, published analysis.
> 3. `Delete` and `Delete (permanently)` (purge) any datasets that you don't want to share publicly (e.g. restricted data, unpublished results) - but make sure you save them to a private history first.
> 
> {% snippet faqs/galaxy/datasets_deleting.md %}
>
{: .hands_on}

We also suggest providing a description/context of the History in the `History Annotation`.

> <hands-on-title>Annotating a History</hands-on-title>
>
> 1. Add an annotation to your History.
>
> {% snippet faqs/galaxy/histories_annotation.md %}
>
{: .hands_on}

Users can add tags to their Galaxy histories to organise and connect them to specific analyses. Tags improve the findability and traceability of your work, aligning with the FAIR principles. Applying relevant tags before publishing your history is recommended to enhance its clarity and usefulness for both yourself and others.

> <hands-on-title>Adding Tags to a History</hands-on-title>
>
> 1. Add Tags to your History
>
> {% snippet faqs/galaxy/histories_tagging.md %}
>
{: .hands_on}

### Sharing a History

Galaxy provides a number of options for sharing histories which allows others to import and access the datasets, parameters, and steps of your history.

> <hands-on-title>Sharing a History</hands-on-title>
>
> 1. Choose one or several ways to share your History.
>
> {% snippet faqs/galaxy/histories_sharing.md %}
>
{: .hands_on}

### Exporting a History to Zenodo

Integration of InvenioRDM and Zenodo into Galaxy now enables you to export your history directly to Zenodo. Publishing your history on Zenodo will generate a DOI for your history, which you can include in your publication and share with other researchers to make your work easily citable.

For instructions on connecting to Zenodo and exporting your history, see the following tutorial on the Galaxy Training Network:

- [Integrating InvenioRDM-compatible Repositories with Galaxy]({% link topics/fair/tutorials/inveniordm-integration/tutorial.md %})

## Making a Workflow public

We recommend publishing a workflow alongside your data and in your article. In addition to providing a complete description of your method, the published workflow will also make your approach reproducible and citeable.

If you have not created and used a workflow for your analysis, Galaxy can automatically create a workflow based on the analysis you have performed in a history.

> <hands-on-title>Extract a Workflow from a History</hands-on-title>
>
> 1. Extract a Workflow from a Galaxy History
>
> {% snippet faqs/galaxy/workflows_extract_from_history.md %}
>
{: .hands_on}

For instructions on preparing your workflow, see the following tutorial on the Galaxy Training Network:

- [Annotate, prepare tests and publish Galaxy workflows in workflow registries]({% link topics/galaxy-interface/tutorials/workflow-fairification/tutorial.md %})

Two repository options that we recommend for publishing your Workflows are:

- WorkflowHub ([https://workflowhub.eu/](https://workflowhub.eu/))
- Dockstore ([https://dockstore.org/](https://dockstore.org/)) 

## Citing Public Workflows

If you used a workflow that was published on Galaxy, in Galaxy's workflow library, [the IWC](https://iwc.galaxyproject.org/)], or a repository such as WorkflowHub, then we recommend you cite the original workflow. You can cite a workflow as you would any other software. If the information is available, we recommend including the creator’s name, the date of creation, the name of the repository, and the workflow’s DOI or URL.

## Citing Public Data

If you used a public data source in your work, we recommend you cite both the data source publication DOI and the dataset itself. Most citation styles will include options for citing datasets. You should include the repository name and accession code in your Methods section to make it easy to find.

## Extracting Tool Citations from a Galaxy History

Galaxy can all citations of the tools used in your history. You can directly extract them from your history. It is the quickest way to add all these tools to your references, and you might find this helpful when writing up your methods.

> <hands-on-title>Extracting Tool Citations</hands-on-title>
>
> 1. Export a list of tool citations from your History
>
> {% snippet faqs/galaxy/cite-tools-used-in-a-history.md %}
> 
> > <comment-title></comment-title>
> > 
> > You will also get the Galaxy Project publication in all History citation exports!
> > 
> {: .comment}
> 
{: .hands_on}

 
# Conclusion

In this tutorial, we have covered how to cite Galaxy Project and the specific Galaxy server that you used and how to prepare your methodology based on your Galaxy History.



