---
layout: tutorial_hands_on

title: "InterMine integration with Galaxy"
zenodo_link: "https://zenodo.org/record/3407174"
questions:
    - How to export your query results from your InterMine of choice to Galaxy?
    - How to export a list of identifiers from Galaxy to your InterMine of choice?
objectives:
    - Learn how to import/export data from/to InterMine instances
    - Understand the InterMine Interchange Dataset
time_estimation: 1h
contributors:
    - danielabutano
    - yochannah
subtopic: analyse
---

# Introduction


InterMine ({% cite Smith2012 %}) is a well-establish platform to integrate and access life sciences data.
It provides the integrated data via a web interface and RESTful web services.

Other organizations download and deploy InterMine on their servers:
there are more than 30 instances over the world (registered at [registry.intermine.org](http://registry.intermine.org)), covering many organism,
including human data, model animals, plants and drug targets.

InterMine has been integrated with Galaxy: the InterMine tool server in Galaxy allows
to import the data returned by any InterMine search and viceversa, using the InterMine Interchange format
it's possible to export a list of identifiers from Galaxy into any InterMine instance of your choice.

Learn more in this tutorial.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Import data from InterMine

> <hands-on-title>Import</hands-on-title>
> Search Galaxy for `InterMine` (not case sensitive; `intermine` is fine too), and click on **InterMine Server** under **Get Data**.
>
> 1. {% tool [InterMine Server](intermine) %}
>
> 2. This will redirect you to the InterMine registry, which shows a full list of InterMines and the various organisms they support. Find an InterMine that has the organism type you’re working with, and click on it to redirect to that InterMine.
>
> 3. Once you arrive at your InterMine of choice, you can run a query as normal - this could be a search, a list results page, a template, or a query in the query builder. Eventually you’ll be presented with an InterMine results table.
>
> 4. Click on **Export** (top right). This will bring up a modal window.
> 5. Select **Send to Galaxy** and double-check the *"Galaxy Location"* is correct.
> 6. Click on the **Send to Galaxy** button on the bottom right of the pop-up window.
>
>    > <tip-title>Enable popups</tip-title>
>    >
>    > If you get an error when you click on the **Send to Galaxy** button, please make sure to allow popups and try again.
>    {: .tip}
>
{: .hands_on}
You have now exported your query results from InterMine to Galaxy.


# Export identifiers into InterMine

## Get data

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Import some fly data from [Zenodo](https://zenodo.org/record/3407174) or from the data library
>
>    ```
>    https://zenodo.org/record/3407174/files/GenesLocatedOnChromosome4.tsv
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 2. Rename the dataset to `GenesLocatedOnChromosome4`
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
> 3. Inspect the data
>
{: .hands_on}

The dataset contains the secondary identifier and the symbol of the *Drosophila melanogaster* genes and their location on the chromosome 4

> <question-title></question-title>
>
> Do the data contain the type, e.g `Protein` or `Gene`?
>
> > <solution-title></solution-title>
> >
> > No, they don't. So we have to specify it, when we create the InterMine Interchange file
> >
> {: .solution}
>
{: .question}

## Create InterMine Interchange dataset

We will use **Create InterMine Interchange Dataset** {% icon tool %} in order to generate an intermediate file which will be used to send the identifiers (e.g. gene identifiers) to InterMine. This file requires the identifier's type (e.g. `Gene`), the identifier (e.g `WBGene00007063`) and, optionally, the organims's name.

> <hands-on-title>Generate InterMine file</hands-on-title>
>
> 1. {% tool [Create InterMine Interchange dateset](toolshed.g2.bx.psu.edu/repos/iuc/intermine_galaxy_exchange/galaxy_intermine_exchange/0.0.1) %} with the following parameters:
>    - {% icon param-file %} *"Tabular file"*: select the `GenesLocatedOnChromosome4` dataset which contains some fly's genes
>    - *"Feature Type Column"*: `Column: 1`
>    - *"Feature Type"*: `Gene`
>    - *"Feature Identifier column"*: `Column: 2`
>
>    > <comment-title></comment-title>
>    > - In this example, because the `GenesLocatedOnChromosome4` dataset does not contain the type we have to specify it, in the *"Feature Type"*
>    > - *"Feature Type"*: this is type of the identifiers you are exporting to InterMine, in this example `Gene`. It must be a class in the InterMine data model.
>    > - *"Feature Identifier column"*: select a column from the input file which contains the identifier. We have selected Column 2, which contains the gene symbol.
>    > - *"Feature Identifier"*: This could be, as an example, a gene symbol like `GATA1` or another other identifier, e.g. `FBGN0000099` or perhaps a  protein accession. In our example we do not have to edit anything because the values for this field are contained in the `GenesLocatedOnChromosome4` dataset, in *Column 2*.
>    > - *"Organism Name column"*: select a column from the input file which contains the organism's name, if you have multiple organisms in the same dataset.
>    > - *"Organism Name"*: alternatively you can directly provide the organism's name. The organims' name is not mandatory, but is good to provide if it is known. It does not have to be precise
>    {: .comment}
>
> 2. Click on **Execute**
>
{: .hands_on}

## Send identifiers to InterMine

Once the generation of the interchange dataset has been completed, open the green box related to **Create InterMine Interchange on data**.

> <hands-on-title>Send data</hands-on-title>
>
> 1. Click on view intermine at **Registry** to be redirected to the InterMine registry, which shows a full list of InterMines and the various organisms they support.
> 2. Find an InterMine that has the organism type you’re working with, in our case FlyMine, and click on the **Send to** green button to export the identifiers to.
>    3. You are redirected to FlyMine, in the List Analysis page showing the identifiers you have just exported from Galaxy.
>
{: .hands_on}

# Conclusion

You have now exported your identifiers from Galaxy to InterMine.