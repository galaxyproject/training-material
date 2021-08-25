---
layout: tutorial_hands_on

title: Submitting raw sequencing reads to ENA
zenodo_link: https://doi.org/10.5281/zenodo.5163611
questions:
- How do you submit raw viral sequence reads to the European Nucleotide Archive?
objectives:
- Add ENA Webin credentials to your Galaxy user information
- Use Galaxy's 'ENA upload tool' to interactively generate metadata
- Use a metadata template to upload bulk metadata
- Submit raw sequencing reads and metadata to ENA's test server
time_estimation: "1h"
level: Intermediate
key_points:
- Use Galaxy's 'ENA Upload tool' to submit raw SARS-CoV-2 reads to ENA
- You need to include your ENA Webin credentials in Galaxy
- For small submission use 'ENA Upload tool' interactive metadata forms feature
- For bulk submissions use a spreadsheet metadata template


requirements:
  -
    type: "internal"
    topic_name: sequence-analysis
    tutorials:
      - human-reads-removal
  -
    type: "internal"
    topic_name: galaxy-interface
    tutorials:
      - collections
  -
    type: "internal"
    topic_name: galaxy-interface
    tutorials:
      - upload-data-to-ena
tags:
  - covid19

contributors:
- roncoronimiguel

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

Raw reads contain valuable information, such as coverage depth and quality scores, that is lost in a consensus sequence.
Submission of raw SARS-CoV-2 reads to public repositories allows reuse of data and reproducibility of analysis and enables discovery of minor allelic variants and [intrahost variation]{% cite Maier2021 %}.

The European Nucleotide Archive is an Open and FAIR repository of nucleotide data. As part of the International Nucleotide Sequence Database Collaboration (INSDC), ENA also indexes data from the NCBI and DDBJ {% cite Arita2020 %}. Data submitted to ENA must be accompanied by sufficient metadata. You can learn more from this [introductory slide deck](https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/upload-data-to-ena/slides.html) or directly from [ENA](https://ena-docs.readthedocs.io/en/latest/submit/general-guide/metadata.html).


In this tutorial we will show you how to use Galaxy's 'ENA Upload tool' to submit SARS-CoV-2 raw sequencing reads and its associated metadata to ENA {% cite Roncoroni2021 %}. You will learn to add your ENA Webin credentials to Galaxy, input metadata interactively or via a metadata template and submit the reads to ENA (test) server.
Specifically, we will use one ONT sequencing file to demonstrate interactive metadata input and two sets of PE Illumina reads to demonstrate how to use the metadata template. Data will be submitted to ENA's test server and will not be public.

> ### {% icon comment %} Nature of the input data
> We will use data derived from sequencing data of bronchoalveolar lavage fluid (BALF) samples obtained from early COVID-19 patients in China as our input data.
> Human traces have been removed in [Galaxy](https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/human-reads-removal/tutorial.html).
>
{: .comment}


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Adding ENA Webin credentials to your Galaxy user information

In order to submit data to ENA, you need to have a valid Webin account. If you don't have one already you can register for one [here](https://www.ebi.ac.uk/ena/submit/sra/#registration). Webin credentials need to be included in your Galaxy user information before you can use the {% tool [ENA Upload tool](toolshed.g2.bx.psu.edu/repos/iuc/ena_upload/ena_upload/0.3.2) %}.

> ### {% icon hands_on %} Hands-on: Add Webin credentials to your Galaxy user information
>
> 1. If you have not already done so, log in to usegalaxy.eu
> 2. Navigate to *"User"* > *"Preferences"* on the top menu
>   - Click on <i class="fa fa-user" aria-hidden="true"></i> **Manage Information**
>     - Scroll down to *"Your ENA Webin account details"* and fill in your ***ENA Webin ID*** and ***Password***
>![ENA Webin Account details in Galaxy](../../images/upload-data-to-ena/ENA-credentials.png "ENA Webin Account details")
{: .hands_on}


# Option 1: submitting to ENA using interactive metadata generator

In this first example, you will submit one ONT sequence file using {% tool [ENA Upload tool](toolshed.g2.bx.psu.edu/repos/iuc/ena_upload/ena_upload/0.3.2) %} interactive metadata forms. This method is only convenient for small submissions. For bulk submissions, we recommend you use the metadata template described below in Option 2.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Upload the ONT data from Zenodo via URLs
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    The URL for our example data is this:
>
>    ```
>    https://zenodo.org/record/5176347/files/SRR10902284_ONT.fq.gz
>    ```
>
>
{: .hands_on}

Once the data is uploaded, we fill the metadata using {% tool [ENA Upload tool](toolshed.g2.bx.psu.edu/repos/iuc/ena_upload/ena_upload/0.3.2) %}'. Interactive metadata forms are nested to fit [ENA's metadata model](https://ena-docs.readthedocs.io/en/latest/submit/general-guide/metadata.html). Briefly, you add Samples to a Study, Experiments to Samples and Runs to Experiments.

We recommend always submitting to the test server before submitting to the public one.
After you confirm that all the data and metadata looks ok, you can go ahead and submit to the public ENA server.




> ### {% icon hands_on %} Hands-on: add metadata interactively and submit a single sequence to ENA
>
> 1. {% tool [ENA Upload tool](toolshed.g2.bx.psu.edu/repos/iuc/ena_upload/ena_upload/0.3.2) %}:
>    - *"Submit to test ENA server?"*: `yes`
>    - *"Would you like to submit pregenerated table files or interactively define the input structures?"*: `Interactive generation of the study structure`
>    - *"Does your submission contains viral samples?"*: `yes`
> 2. Fill all metadata boxes and make sure that:
>    - *"Enter the species of the sample"*: `Severe acute respiratory syndrome coronavirus 2`
>    - *"Enter the taxonomic ID corresponding to the sample species"*: `2697049`
>    - *"Host subject id"*: avoid using ID that can be use to trace samples back to patients
>    - *"Host scientific name"*: `Homo sapiens`
>    - *"Library layout"*: `SINGLE`
>    - *"Select the sequencing platform used"*: `Oxford Nanopore`
>    - *"Instrument model"*: `minION`
>    - *"Runs executed within this experiment"*
>      - *`File(s) associated with this run`*: {% icon param-files %} select the uploaded ONT fastq.gz file
>    - *"Affiliation center"*: your institution
>
{: .hands_on}

> ### {% icon warning %} Submit to the test server first
> Make sure  *"Submit to test ENA server?"*: `yes`. Otherwise your data will be submitted to the public server.
{: .warning}

Four metadata tables (Study, Sample, Experiment and Run), and a metadata ticket with submission information are generated. You can confirm a successful submission at ENA [test server](https://wwwdev.ebi.ac.uk/ena/submit/sra) (or the [public server](https://www.ebi.ac.uk/ena/submit/sra), if you chose it).







# Option 2: submitting to ENA using the metadata template

For larger submissions, interactive metadata input can be tedious and not practical. In the second example, you will submit two sets of Illumina PE sequence files and input metadata using a [template spreadsheet](https://drive.google.com/file/d/1Gx78GKh58PmRjdmJ05DBbpObAL-3oUFX/view?usp=sharing). For this exercise, we provide you with a pre-filled template and encourage you to explore it.

> ### {% icon hands_on %} Hands-on: Upload and inspect data
>
> 1. Upload the ONT data from Zenodo via URLs:
>
>    ```
>    https://zenodo.org/record/5176347/files/GTN_tutorial_mock_metadata_template.xlsx
>    https://zenodo.org/record/5176347/files/SRR10903401_1.fastq.gz
>    https://zenodo.org/record/5176347/files/SRR10903401_2.fastq.gz
>    https://zenodo.org/record/5176347/files/SRR10903402_1.fastq.gz
>    https://zenodo.org/record/5176347/files/SRR10903402_2.fastq.gz
>    ```
> 2. Arrange the data into a paired dataset collection
>
>    {% snippet faqs/galaxy/collections_build_list_paired.md %}
>
>    For the example datasets this means:
>    - You need to tell Galaxy about the suffix for your forward and reverse reads, respectively:
>      - set the text of *unpaired forward* to: `_1.fastq.gz`
>      - set the text of *unpaired reverse* to: `_2.fastq.gz`
>      - click: `Auto-pair`
>
>      All datasets should now be moved to the *paired section* of the dialog, and the middle column there should show that only the sample accession numbers, *i.e.* `SRR10903401` and `SRR10903402`, will be used as the pair names.
>
>    - Make sure *Hide original elements* is checked to obtain a cleaned-up history after building the collection.
>    - Click *Create Collection*
>
> 3. Inspect the `GTN_tutorial_mock_metadata.xlsx` (filled-in template) file by clicking on the {% icon galaxy-eye %} (eye) icon
>
>    > ### {% icon question %} Questions
>    >
>    > 1. How many metadata sheets are there?
>    > 2. Which metadata section is different from the corresponding section in the interactive metadata input?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. There are four metadata sheets, one per metadata object (Study, Sample, Experiment, Run)
>    > > 2. The Sample section is more extensive in the template spreadsheet, because it contains `Mandatory`, `Recommended` and `Optional` fields, whereas the interactive metadata `Sample` form contains only `Mandatory` ones.
>    >  {: .solution }
>    {: .question}
>
{: .hands_on}


As before, the submission is done to the test server before submitting to the public one.

> ### {% icon hands_on %} Hands-on: add metadata interactively and submit a single sequence to ENA
>
> 1. {% tool [ENA Upload tool](toolshed.g2.bx.psu.edu/repos/iuc/ena_upload/ena_upload/0.3.2) %}:
>    - *"Submit to test ENA server?"*: `yes`
>    - *"Would you like to submit pregenerated table files or interactively define the input structures?"*: `User generated metadata tables based on Excel template`
>    - *"Does your submission contains viral samples?"*: `yes`
>    - *"Select Excel (xlsx) file based on templates"*: {% icon param-files %} select the uploaded .xlsx template file
>    - *"Select runs input format"*: `Input from a paired collection`
>      - *"List of paired-end runs files"*: : select the PE collection containing the PE sequencing reads
>    - *"Affiliation center"*: your institution
>
{: .hands_on}

> ### {% icon warning %} Submit to the test server first
> Make sure  *"Submit to test ENA server?"*: `yes`. Otherwise your data will be submitted to the public server.
{: .warning}

Four metadata tables (Study, Sample, Experiment and Run), and a metadata ticket with submission information are generated. You can confirm a successful submission at ENA [test server](https://wwwdev.ebi.ac.uk/ena/submit/sra) (or the [public server](https://www.ebi.ac.uk/ena/submit/sra), if you chose it).
