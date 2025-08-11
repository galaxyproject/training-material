---
layout: tutorial_hands_on

title: Prepare data from CbioPortal for Flexynesis integration
zenodo_link: https://zenodo.org/records/16287482
questions:
- How to download data from cBioPortal?
- How to prepare omics data for Flexynesis integration.
objectives:
- Download Breast cancer data from Metabric through cBioportal using Flexynesis
- Clean and preprocess genomics data
- Format data for downstream integration analysis
time_estimation: 1H
key_points:
- cBioportal is a repository for accessible and interpretable cancer genomic data
- Flexynesis comprehensive tool can be used to make data ready for integration.
contributors:
- Nilchia
- plushz
- bgruening
answer_histories:
    - label: "usegalaxy.eu"
      history: https://usegalaxy.eu/u/nilchia/h/final-flexynesis-get-data-from-cbioportal-test-wf
      date: 2025-08-01
---

The cBioPortal is an open-access web platform that provides intuitive access to large-scale cancer genomics datasets {% cite Gao2013 %} {% cite cbioportal-website %}. Originally developed to make complex molecular profiling data more accessible to the broader research community, cBioPortal hosts data from thousands of cancer studies and tens of thousands of tumor samples.

In this tutorial, we will work with data from the METABRIC consortium {% cite metabric-website %}, one of the landmark breast cancer genomics studies available through cBioPortal. This dataset contains comprehensive molecular and clinical data from over 2,000 breast cancer patients, including gene expression profiles, copy number alterations, mutation data, and clinical outcomes.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Import data from CbioPortal

> <hands-on-title> Import data from cBioPortal </hands-on-title>
>
> 1. {% tool [Flexynesis cBioPortal import](toolshed.g2.bx.psu.edu/repos/bgruening/flexynesis_cbioportal_import/flexynesis_cbioportal_import/0.2.20+galaxy3) %} with the following parameters:
>    - *"I certify that I am not using this tool for commercial purposes."*: `Yes`
>    - *"cBioPortal study ID"*: `brca_metabric`
>
{: .hands_on}

> <question-title></question-title>
>
> 1. What modalities are imported to Galaxy?
>
> > <solution-title></solution-title>
> >
> > 1. * Clinical data
> >    * Copy number alteration
> >    * Methylation
> >    * Gene expression
> >    * Mutation
> >
> {: .solution}
>
{: .question}


# Data cleanup

Now we need to clean up our data. This means removing comment lines in the matrix, removing duplicate samples and so on.

## Clinical data

> <hands-on-title> Prepare clinical data </hands-on-title>
>
> Here we'll First extract the clinical data from the collection, then we remove the comment lines, and finally we remove the extra index column.
>
> 1. {% tool [Extract dataset](__EXTRACT_DATASET__) %} with the following parameters:
>    - {% icon param-file %} *"Input List"*: `datasets` (output of **Flexynesis cBioPortal import** {% icon tool %})
>    - *"How should a dataset be selected?"*: `Select by element identifier`
>        - *"Element identifier:"*: `data_clinical_patient`
>
>
> 2. {% tool [Table Compute](toolshed.g2.bx.psu.edu/repos/iuc/table_compute/table_compute/1.2.4+galaxy2) %} with the following parameters:
>    - *"Input Single or Multiple Tables"*: `Single Table`
>        - {% icon param-file %} *"Table"*: `data_clinical_patient` (output of **Extract dataset** {% icon tool %})
>        - In *"Advanced File Options "*:
>            - *"Header begins at line N"*: `4`
>        - *"Type of table operation"*: `No operation (just reformat on output)`
>
> 3. {% tool [Advanced Cut](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cut_tool/9.5+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"File to cut"*: `table` (output of **Table Compute** {% icon tool %})
>    - *"Operation"*: `Discard`
>    - *"Cut by"*: `fields`
>        - *"Is there a header for the data's columns ?"*: `Yes`
>            - *"List of Fields"*: `c1:`
>
> 4. Rename the data `clinical data - cleaned`
>
{: .hands_on}

> <question-title></question-title>
>
> 1. How many samples are there in the clinical data?
>
> > <solution-title></solution-title>
> >
> > 1. Click on the data, you will see the data has **2510** lines, so **2509** samples.
> >
> {: .solution}
>
{: .question}

## Omics data

Time to prepare our omics data. We are interested in mutation and gene expression.

> <hands-on-title> Prepare gene expression data </hands-on-title>
>
> Here we'll First extract the clinical data from the collection, then we remove duplicate genes from the matrix, and finally we remove the ENTREZ Ids.
>
> 1. {% tool [Extract dataset](__EXTRACT_DATASET__) %} with the following parameters:
>    - {% icon param-file %} *"Input List"*: `data` (output of **Flexynesis cBioPortal import** {% icon tool %})
>    - *"How should a dataset be selected?"*: `Select by element identifier`
>        - *"Element identifier:"*: `data_mrna_illumina_microarray`
>
> 2. {% tool [Sort](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_sort_header_tool/9.5+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Sort Query"*: `data_mrna_illumina_microarray` (output of **Extract dataset** {% icon tool %})
>    - *"Number of header lines"*: `1`
>    - In *"Column selections"*:
>        - {% icon param-repeat %} *"Insert Column selections"*
>            - *"on column"*: `c1`
>    - *"Output unique values"*: `Yes`
>
> 3. {% tool [Advanced Cut](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cut_tool/9.5+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"File to cut"*: `output` (output of **Sort** {% icon tool %})
>    - *"Operation"*: `Discard`
>    - *"Cut by"*: `fields`
>        - *"Is there a header for the data's columns ?"*: `Yes`
>            - *"List of Fields"*: `c2: Entrez_Gene_Id`
>
> 4. Rename the data `expression data - cleaned`
>
{: .hands_on}

> <question-title></question-title>
>
> 1. How many samples and genes are there in the gene expression data?
>
> > <solution-title></solution-title>
> >
> > 1. Click on the data, you will see the data has **20386** lines and 1981 columns, so **20385** genes and 1980 samples.
> >
> {: .solution}
>
{: .question}

> <hands-on-title> Prepare mutation data </hands-on-title>
>
> Here we'll First extract the clinical data from the collection, then we remove comment lines from the matrix, and remove the extra index column, and finally, we create a binarized matrix which indicates the number of mutations per genes.
>
> 1. {% tool [Extract dataset](__EXTRACT_DATASET__) %} with the following parameters:
>    - {% icon param-file %} *"Input List"*: `data` (output of **Flexynesis cBioPortal import** {% icon tool %})
>    - *"How should a dataset be selected?"*: `Select by element identifier`
>        - *"Element identifier:"*: `data_mutations`
>
> 2. {% tool [Table Compute](toolshed.g2.bx.psu.edu/repos/iuc/table_compute/table_compute/1.2.4+galaxy2) %} with the following parameters:
>    - *"Input Single or Multiple Tables"*: `Single Table`
>        - {% icon param-file %} *"Table"*: `data_mutations` (output of **Extract dataset** {% icon tool %})
>        - In *"Advanced File Options "*:
>            - *"Header begins at line N"*: `1`
>        - *"Type of table operation"*: `No operation (just reformat on output)`
>
> 3. {% tool [Advanced Cut](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cut_tool/9.5+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"File to cut"*: `table` (output of **Table Compute** {% icon tool %})
>    - *"Operation"*: `Discard`
>    - *"Cut by"*: `fields`
>        - *"Is there a header for the data's columns ?"*: `Yes`
>            - *"List of Fields"*: `c1:`
>
> 4. {% tool [Flexynesis utils](toolshed.g2.bx.psu.edu/repos/bgruening/flexynesis_utils/flexynesis_utils/0.2.20+galaxy3) %} with the following parameters:
>    - *"I certify that I am not using this tool for commercial purposes."*: `Yes`
>    - *"Flexynesis utils"*: `Binarize mutation data`
>        - {% icon param-file %} *"Mutation data"*: `table` (output of **Advanced Cut** {% icon tool %})
>        - *"Column in the mutation file with genes"*: `Column: 1`
>        - *"Column in the mutation file with samples"*: `Column: 17`
>
> 5. Rename the data `mutation data - cleaned`
>
{: .hands_on}

> <question-title></question-title>
>
> 1. How many samples and genes are there in the mutation data?
>
> > <solution-title></solution-title>
> >
> > 1. Click on the data, you will see the data has **174** lines and **2370** columns, so **174** genes and **2370** samples.
> >
> {: .solution}
>
{: .question}

# Split data to train and test

In the last step we split our clinical and omics data into train and test with ratio of 0.7 (70% as training and 30% as test)

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Flexynesis utils](toolshed.g2.bx.psu.edu/repos/bgruening/flexynesis_utils/flexynesis_utils/0.2.20+galaxy2) %} with the following parameters:
>    - *"I certify that I am not using this tool for commercial purposes."*: `Yes`
>    - *"Flexynesis utils"*: `Split data to train and test`
>        - {% icon param-file %} *"Clinical data"*: `clinical data - cleaned` (output of **Advanced Cut** {% icon tool %})
>        - {% icon param-files %} *"Omics data"*: `expression data - cleaned` (output of **Advanced Cut** {% icon tool %}), `mutation data - cleaned` (output of **Flexynesis utils** {% icon tool %})
>
{: .hands_on}

# Conclusion

In this tutorial, we showed how to download and prepare multi-modal cancer genomics data from cBioPortal for integration using Flexynesis. Working with the METABRIC breast cancer dataset, we covered the complete data preparation pipeline from raw data access to analysis-ready formats.