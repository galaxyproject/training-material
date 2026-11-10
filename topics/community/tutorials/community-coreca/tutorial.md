---
layout: tutorial_hands_on
title: Creation of the CoReCa (Community Research Catalogue) for your community
level: Introductory
subtopic: sig

questions:
- How can I create a Community Resource Catalogue (CoReCa) for a community?
objectives:
- Create CoReCa for communities
time_estimation: 1H
key_points:
- The CoReCa is a list of resources that enables communities to list specific tools, workflows and bioconda recipies.
tags:
- Community
- SIG
contributions:
  authorship:
    - scorreard
  editing:
    - scorreard
---

The rapid growth of research infrastructure has led to the development of numerous bioinformatics tools and databases, but information about them remains fragmented across disparate, often outdated cataloging efforts, hindering resource discovery and utilization.

To address this critical gap, the ELIXIR Microbiome Community collaborated with the **Research Software Ecosystem** to create MiCoReCa (Microbiome Community Resource Catalogue), a comprehensive, dynamic, open-access catalogue of microbiome-related bioinformatics resources.

This resource is now doable by other communities, and are grouped under the **CoReCa (Community Resource Catalogue)**.

The aim of this tutorial is to explain how to create the files necessary to create the CoReCa for your community.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# How to generate a CoReCa?

## Add your community to the CoReCa communities

You first need to check if your community is in the [coreca repo](https://github.com/research-software-ecosystem/coreca/tree/main/communities), a central resource for communities.
If the community is already there, you can move to the next step of this tutorial.

If your community is not already included, follow these steps :
> <hands-on-title>Add your community to the coreca repository</hands-on-title>
>
> You need to create a new folder in the coreca/communities folder within coreca code source.
> 1. If not already done, fork the [coreca repository](https://github.com/research-software-ecosystem/coreca)
> 2. Go to the `communities` folder
> 3. Click on **Add file** in the drop-down menu at the top
> 4. Select **Create a new file**
> 5. Fill in the `Name of your file` field with:  <name of your community> + `/metadata/keywords.yml`
>    This will create a folder for your community, a new folder for your community called metadata and a file called keywords.yml.
> 6. Click on **Commit changes**
> 7. Fill in the commit message with something like `Add <community-name> community`
> 8. Click on `Create a new branch for this commit and start a pull request`
> 9. Create the pull request by following the instructions
>
{: .hands_on}

## Extract relevant resources based on EDAM terms.

The generation of a CoReCa mostly rely on a list of keywords (EDAM terms and operations) that are used to filter the resources.

> <hands-on-title>Add the EDAM categories to the keywords file</hands-on-title>
> 1. Open or create the file named `keywords.yml` in your community metadata folder (`communities/<your community>/metadata/keywords.yml`)
> 2. Add the EDAM terms relevant to your community in the `keywords.yml` file you started above, with 1 term per row. 
>    For example:
>    ```
>    edam:
>        topics:
>            - Metagenomics
>            - Metatranscriptomics
>        operations:
>            - Read binning
>
>    keywords:
>        - metage.*
>        - metabar.*
>
>    acronyms:
>        - ITS
>        - OTU
>    ```
> 3. Click on **Commit changes** at the top
> 4. Fill in the commit message with something like `Add EDAM categories for community <community-name>`
> 3. Click on `Create a new branch for this commit and start a pull request`
> 4. Create the pull request by following the instructions
>
{: .hands_on}

The Pull Request will be reviewed. Make sure to respond to any feedback.

On the Sunday following the merging of the previous pull request, several folders will be generated in you community folder:
- The `bìoconda` folder with several json and tsv files
- The `workflowhub` folder with several json and tsv files
- The `rsec` folder with several folder within it

## Curate filtered resources

As of today, only the bioconda recipies and the workflows are possible to curate.
The files generated at the previous step contain all the conda recipies / workflows associated with the EDAM categories that you included in the keywords.yml file. However, not all of these tools might be interesting for your community.

CoReCa allows for an additional optional manual curation for conda recipies and workflows that can be selected by the community curator (maybe that is you!).

The curation information must be stored in the files called `*_status.tsv` 

### Curation of Bioconda recipies

> <hands-on-title>Review Bioconda recipies in your community folder</hands-on-title>
>
> 1. Download the `recipes_status.tsv` file located in `communities/<community-name>/bioconda/recipes_status.tsv`
> 2. Open `recipes_status.tsv` with a Spreadsheet Software.
> 3. Review each line corresponding to a bioconda recipie.
>
>    You can also just review some tools. 
> 1. Add `TRUE` to the `keep` column if the recipie should be kept, and `FALSE` if not.
> 3. Export the new table as TSV.
> 4. Submit the TSV as `recipes_status.tsv` in your community folder.
> 5. Make a Pull request to submit the curation information
> 6. Wait for the Pull Request to be merged
>
{: .hands_on}

On the Sunday following the merging of the previous pull request, `recipes_curated.json` and `recipes_curated.tsv` files will be updated in the `communities/<community-name>/bioconda` folder reflecting your curation effort.
You or other community members can review bioconda recipies in your community from this page to make step-by-step changes or updates to the `recipes_status.tsv` file as needed.

We recommend collaborative work be done with a cloud-based system such as a Google Sheets page.

### Curation of workflows

> <hands-on-title>Review workflows in your community folder</hands-on-title>
>
> 1. Download the `workflows_status.tsv` file located in `communities/<community-name>/workflowhub/workflows_status.tsv`
> 2. Open `workflows_status.tsv` with a Spreadsheet Software.
> 3. Review each line corresponding to a workflow.
>
>    You can also just review some tools. 
> 1. Add `TRUE` to the `To keep` column if the recipie should be kept, and `FALSE` if not.
> 3. Export the new table as TSV.
> 4. Submit the TSV as `workflows_status.tsv` in your community folder.
> 5. Make a Pull request to submit the curation information
> 6. Wait for the Pull Request to be merged
>
{: .hands_on}

On the Sunday following the merging of the previous pull request, `workflows_curated.json` and `workflows_curated.tsv` files will be updated in the `communities/<community-name>/workflowhub` folder reflecting your curation effort.
You or other community members can review workflows in your community from this page to make step-by-step changes or updates to the `workflows_status.tsv` file as needed.

We recommend collaborative work be done with a cloud-based system such as a Google Sheets page.

# Conclusion

You now have curated resources listing resources available for your community.
