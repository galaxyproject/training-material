---
layout: tutorial_hands_on

title: Trio Analysis using Synthetic Datasets from RD-Connect GPAP
subtopic: human-genetics-cancer
zenodo_link: 'https://doi.org/10.5281/zenodo.6483454'
questions:
- How do you import data from the EGA?
- How to download files with HTSGET in Galaxy?
objectives:
- Requesting DAC access and importing data from the EGA.
time_estimation: 2H
key_points:
- Downloading whole datasets with HTSGET is safe and easy with Galaxy.
contributions:
  authorship:
  - JasperO98
  editing:
  - wm75
  - hexylena
  - shiltemann
tags:
- ega
- get-data
subtopic: upload
---

In this tutorial we will also make use of the HTSGET protocol, which is a program to download our data securely and savely. This protocol has been implemented in the {% tool [EGA Download Client](toolshed.g2.bx.psu.edu/repos/iuc/ega_download_client/pyega3/4.0.0+galaxy0) %} tool, so we don't have to leave Galaxy to retrieve our data.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% include _includes/cyoa-choices.html option1="DAC" option2="Test" default="Test"
       text="If you want to learn about getting DAC access to a datset, you can take this detour. Otherwise the tutorial will use the default test dataset available to everyone." %}

<div class="DAC" markdown="1">

## Getting DAC access

Our test data is stored in EGA, which can be easily accessed using the EGA Download Client. Our specific EGA dataset accession ID is: "EGAD00001008392". However, before you can access this data you need to request DAC access to this dataset. This can be requested by emailing to <helpdesk@ega-archive.org>, don’t forget to mention the dataset ID! When the EGA grants you access it will create an account for you, if you don't have it already. Next, you should link your account to your Galaxy account by going to the homepage on Galaxy and at the top bar click **User > Preferences > Manage Information**. Now add your email and password of your (new) EGA account under: **Your EGA (european Genome Archive) account**. After that you can check if you can log in and see if you have access to the dataset.

</div>

> <hands-on-title>Check log-in and authorized datasets</hands-on-title>
>
> 1. {% tool [EGA Download Client](toolshed.g2.bx.psu.edu/repos/iuc/ega_download_client/pyega3/4.0.0+galaxy0) %} with the following parameters:
>    - *"What would you like to do?"*: `List my authorized datasets`
>
> > <comment-title>Check if the dataset is listed.</comment-title>
> >
> > Check if your dataset is listed in the output of the tool. If not you can look at the header of the output to find out why it is not listed. When the header does not provide any information you could have a look at the error message by clicking on the **eye** {% icon galaxy-eye %} of the output dataset and then click on the icon **view details** {% icon details %}. The error should be listed at **Tool Standard Error**.
> {: .comment}
>
{: .hands_on}

## Download list of files

When you have access to the EGA dataset, you can download all the needed files. However, the EGA dataset contains many different filetypes and cases, but we are only interested in the VCFs from case 5 and, to reduce execution time, the variants on chromosome 17. To be able to donwload these files we first need to request the list of files from which we can download. Make sure to use **version 4+** of the {% tool [EGA Download Client](toolshed.g2.bx.psu.edu/repos/iuc/ega_download_client/pyega3/4.0.0+galaxy0) %}.

> <hands-on-title>Request list of files in the dataset</hands-on-title>
>
> 1. {% tool [EGA Download Client](toolshed.g2.bx.psu.edu/repos/iuc/ega_download_client/pyega3/4.0.0+galaxy0) %} **version 4+** with the following parameters:
>    - *"What would you like to do?"*: `List files in a datasets`
>    - *"EGA Dataset Accession ID?"*: `EGAD00001003338` (or `EGAD00001008392` if you have DAC access)
>
> {% snippet faqs/galaxy/tools_change_version.md %}
>
> {: .comment}
{: .hands_on}


## Filter list of files

Now that we have listed all the files, we need to filter out the files we actually need. We can do this by using a simple regular expression or regex. With regex it will be easy to find or replace patterns within a textfile.

> <hands-on-title>Filter out VCFs from list of files</hands-on-title>
>
> 1. {% tool [Search in textfiles](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_grep_tool/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `List of EGA datasets` (output of **EGA Download Client** {% icon tool %})
>    - *"Type of regex"*: `Extended (egrep)`
>    - *"Regular Expression"*: `^HG.*vcf.gz$` (This regular expression finds things beginning with `HG` and ending with `vcf.gz`
>    - *"Match type"*: `case sensitive`
>
{: .hands_on}

## Download files
After the filtering you should have a tabular file with 3 lines each containing the ID of a VCF file from case 5, as shown below.
```
EGAF00005007327›1›      850737› f3dee64b466efe334b2cac77f5c2f710›       HG01775.chrY.vcf.gz
EGAF00007243775›1›      23033›  51cfb69bf3b9416ff425381a58c18a2b›       HG00408.novoBreak__256r__4.100100-10100100__7.200100-9000100.vcf.gz
EGAF00007243779›1›      15340›  ebad4425191a89d3e970c02190a87175›       HG01890.HGSVC__145r__1.900100-10001000__18.2001000-90001000.vcf.gz
```

> <hands-on-title>Download listed VCFs</hands-on-title>
>
> 1. {% tool [EGA Download Client](toolshed.g2.bx.psu.edu/repos/iuc/ega_download_client/pyega3/4.0.0+galaxy0) %} with the following parameters:
>    - *"What would you like to do?"*: `Download multiple files (based on a file with IDs)`
>        - {% icon param-file %} *"Table with IDs to download"*: `Filtered list of files` (output of **Search in textfiles** {% icon tool %})
>        - *"Column containing the file IDs"*: `Column: 1`
>
> {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}
