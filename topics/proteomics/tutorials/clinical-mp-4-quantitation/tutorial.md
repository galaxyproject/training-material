---
layout: tutorial_hands_on

title: "Clinical Metaproteomics 4: Quantitation"
zenodo_link: "https://doi.org/10.5281/zenodo.10105821"
questions:
- How to perform quantitation?
objectives:
- Perform quantitation using MaxQuant and extract microbial and human proteins and peptides.
time_estimation: 3H
key_points:
- Quantified Microbial and Human peptides/proteins can be analyzed separately so that the results are more comparative.
contributions:
  authorship:
    - subinamehta
    - katherine-d21
    - dechendb
  editing:
    - pratikdjagtap
    - timothygriffin
requirements:
  -
    type: "internal"
    topic_name: proteomics
subtopic: clinical-metaproteomics
follow_up_training:

    -
        type: "internal"
        topic_name: proteomics
        tutorials:
            - clinical-mp-5-data-interpretation
tags: [label-TMT11]
redirect_from:
- /topics/proteomics/tutorials/clinical-mp-quantitation/tutorial

recordings:
- captioners:
  - katherine-d21
  date: '2024-06-21'
  galaxy_version: '23.1'
  length: 8M
  youtube_id: _e0l1AtfZ4Y
  speakers:
  - katherine-d21
---


The next step of the clinical metaproteomics workflow is the quantification workflow. Running a quantification workflow in proteomics is essential for several critical purposes. It allows researchers to measure and compare the abundance of proteins or peptides in biological samples, offering valuable insights into biomarker discovery, comparative analysis, and differential expression studies. Quantitative proteomics helps reveal the functional roles of proteins, the stoichiometry of protein complexes, and the effects of drugs on protein expression in pharmacological studies. Additionally, it serves as a quality control measure, validating initial protein identifications, and providing data normalization for increased accuracy. Quantitative data are indispensable for hypothesis testing, systems biology, and their clinical relevance in areas such as disease diagnosis, prognosis, and therapeutic decision-making. In summary, the quantitation workflow in proteomics is a cornerstone for deciphering the complexities of protein expression and regulation, facilitating a wide array of biological and clinical applications.

In this current workflow, we perform Quantification using the MaxQuant tool and the output will be interpreted in our next module.

![Quantitation workflow]({% link topics/proteomics/images/clinical-mp/clinical-mp-quantification.JPG %})



> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


## Get data

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    https://zenodo.org/records/10105821/files/PTRC_Skubitz_Plex2_F10_9Aug19_Rage_Rep-19-06-08.raw
>    https://zenodo.org/records/10105821/files/PTRC_Skubitz_Plex2_F11_9Aug19_Rage_Rep-19-06-08.raw
>    https://zenodo.org/records/10105821/files/PTRC_Skubitz_Plex2_F13_9Aug19_Rage_Rep-19-06-08.raw
>    https://zenodo.org/records/10105821/files/PTRC_Skubitz_Plex2_F15_9Aug19_Rage_Rep-19-06-08.raw
>    https://zenodo.org/records/10105821/files/Experimental-Design_Discovery_MaxQuant.tabular
>    https://zenodo.org/records/10105821/files/Quantitation_Database_for_MaxQuant.fasta
>    ```
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
> 5. Add to each database a tag corresponding to input files.
> 6. Create a dataset of the RAW files.
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
{: .hands_on}


# Import Workflow

> <hands-on-title>Running the Workflow</hands-on-title>
>
> 1. **Import the workflow** into Galaxy:
>
>    {% snippet faqs/galaxy/workflows_run_trs.md path="topics/proteomics/tutorials/clinical-mp-4-quantitation/workflows/WF4_Quantitation_Workflow.ga" title="Quantitation Workflow" %}
>
> 2. Run **Workflow** {% icon workflow %} using the following parameters:
>    - *"Send results to a new history"*: `No`
>    - {% icon param-file %} *" Quantitation_Database-For-MaxQuant * "*: `Quantitation_Database_for_MaxQuant.fasta`
>    - {% icon param-file %} *" Experimental-Design Discovery MaxQuant"*: `Experimental-Design_Discovery_MaxQuant.tabular`
>    - {% icon param-file %} *" Input Raw-files"*: `RAW dataset collection`
>
>    {% snippet faqs/galaxy/workflows_run.md %}
>
{: .hands_on}


# Peptide quantification

In the [Discovery Module](https://github.com/subinamehta/training-material/blob/main/topics/proteomics/tutorials/clinical-mp-discovery/tutorial.md), we used MaxQuant to identify peptides for verification. Now, we will again use MaxQuant to further quantify the PepQuery-verified peptides, both microbial and human. More information about quantitation using MaxQuant is available, including [Label-free data analysis](https://gxy.io/GTN:T00218) and [MaxQuant and MSstats for the analysis of TMT data](https://gxy.io/GTN:T00220).

The outputs we are most interested in consist of the `MaxQuant Evidence file`, `MaxQuant Protein Groups`, and `MaxQuant Peptides`. The `MaxQuant Peptides` file will allow us to group them to generate a list of quantified microbial peptides.

> <hands-on-title> Quantify verified peptides (from PepQuery2) </hands-on-title>
>
> 1. {% tool [MaxQuant](toolshed.g2.bx.psu.edu/repos/galaxyp/maxquant/maxquant/1.6.17.0+galaxy4) %} with the following parameters:
>    - In *"Input Options"*:
>        - {% icon param-file %} *"FASTA files"*: `Quantitation Database for MaxQuant` (Input dataset)
>    - In *"Search Options"*:
>        - {% icon param-file %} *"Specify an experimental design template (if needed). For detailed                           instructions see the help text."*: `output` (Input dataset)
>        - *"minimum peptide length"*: `8`
>        - *"Match between runs"*: `Yes`
>        - *"Maximum peptide length for unspecific searches"*: `50`
>    - In *"Protein quantification"*:
>        - *"Use only unmodified peptides"*: `Yes`
>            - *"Modifications used in protein quantification"*: `Oxidation (M)`
>        - In *"LFQ Options"*:
>            - *"iBAQ (calculates absolute protein abundances by normalizing to copy number and not protein mass)"*: `No`
>    - In *"Parameter Group"*:
>        - {% icon param-repeat %} *"Insert Parameter Group"*
>            - {% icon param-collection %} *"Infiles"*: `output` (Input dataset collection)
>            - *"fixed modifications"*: `Carbamidomethyl (C)`
>            - *"variable modifications"*: `Oxidation (M)`
>            - *"enzyme"*: `Trypsin/P`
>            - *"Quantitation Methods"*: `reporter ion MS2`
>                - *"isobaric labeling"*: `TMT11plex`
>                - *"Filter by PIF"*: `Yes`
>    - In *"Output Options"*:
>        - *"Select the desired outputs."*: `Protein Groups` `mqpar.xml` `Peptides` `Evidence` `MSMS`
>
>
{: .hands_on}


> <question-title></question-title>
>
> 1. Why can we switch back to using RAW files for MaxQuant, instead of using MGF files?
>
> > <solution-title></solution-title>
> >
> > 1. MaxQuant prefers RAW format compared to MGF as it has more information compared to MGF.
> >
> {: .solution}
>
{: .question}
> <question-title></question-title>
>
> 1. Previously, we used MaxQuant in the Discovery workflow. Why are we using MaxQuant again, instead of Search GUI/PeptideShaker?
>
> > <solution-title></solution-title>
> >
> > 1. We are using MaxQuant for quantification purposes only. SearchGUI Peptide Shaker doesn't have the capability to perform quantification of peptides or proteins.
> >
> {: .solution}
>
{: .question}

## Using Text Manipulation Tools to Manage MaxQuant Outputs

> <hands-on-title> Select microbial protein groups from MaxQuant with Select </hands-on-title>
>
> 1. {% tool [Select](Grep1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `proteinGroups` (output of **MaxQuant** {% icon tool %})
>    - *"that"*: `NOT Matching`
>    - *"the pattern"*: `(_HUMAN)|(_REVERSED)|(REV_)|(CON)|(con)`
>
>
> 2. {% tool [Select](Grep1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `peptides` (output of **MaxQuant** {% icon tool %})
>    - *"that"*: `NOT Matching`
>    - *"the pattern"*: `(_HUMAN)|(_REVERSED)|(REV_)|(CON)|(con)`
>
>
> 3. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1`
>    - {% icon param-file %} *"From"*: `out_file1` (output of **Select** {% icon tool %})
>
>
>
> 4. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1`
>    - {% icon param-file %} *"From"*: `out_file1` (output of **Select** {% icon tool %})
>
>
{: .hands_on}

## Generating a list of quantified proteins and peptides

> <hands-on-title> Group quantified proteins </hands-on-title>
>
> 1. {% tool [Group](Grouping1) %} with the following parameters:
>    - {% icon param-file %} *"Select data"*: `out_file1` (output of **Cut** {% icon tool %})
>    - *"Group by column"*: `c1`
>
>
{: .hands_on}
> <hands-on-title> Group quantified peptides </hands-on-title>
>
> 1. {% tool [Group](Grouping1) %} with the following parameters:
>    - {% icon param-file %} *"Select data"*: `out_file1` (output of **Cut** {% icon tool %})
>    - *"Group by column"*: `c1`
>
>
{: .hands_on}


# Conclusion

In summary, the implementation of a quantitation workflow using MaxQuant represents a significant advancement in quantitative proteomic research. This approach enables precise measurement of protein and peptide abundances, enhancing our ability to unravel the complexities of biological systems. This workflow is instrumental in biomarker discovery, comparative analysis, and understanding differential protein expression by offering detailed insights into quantitative changes across different experimental conditions. Its capacity to generate accurate data supports a wide spectrum of applications, including disease research, drug development, and systems biology investigations. Furthermore, the MaxQuant-based quantitation workflow ensures data quality, enabling reliable and reproducible results. It serves as a vital step for quality control, allowing researchers to draw meaningful conclusions from proteomic experiments confidently.
