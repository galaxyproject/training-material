---
layout: tutorial_hands_on

title: "Genome Alignment"
questions:
objectives:
key_points:
zenodo_url: "https://zenodo.org/record/8010522"
time_estimation:
contributions:
  authorship: [mbourgey]
  editing: [shiltemann]
  funding: [bioinformatics-ca,erasmusplus]

subtopic: alignment-assembly
priority: 4
---

# Data Import

We will be working on a [CageKid](http://www.cng.fr/cagekid/) sample pair, patient `C0098`. The CageKid project is part of ICGC and is focused on renal cancer in many of it’s forms. The raw data can be found on EGA and calls, RNA and DNA, can be found on the [ICGC](https://dcc.icgc.org/) portal.

For practical reasons we subsampled the reads from the sample because running the whole dataset would take way too much time and resources.

> <hands-on-title>Obtaining our data</hands-on-title>
>
> 1. Make sure you have an empty analysis history. Give it a name.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. **Import Data.**
>    - Import the sample FASTQ files to your history, either from a shared data library (if available on your Galaxy),
>      or from Zenodo using the URLs listed in the box below (click {% icon param-repeat %} to expand):
>
>      > <solution-title>List of Zenodo URLs</solution-title>
>      > ```
>      > {{zenodo_url}}/files/run62DPDAAXX_8-normal.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/run62DPDAAXX_8-normal.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/run62DU0AAXX_8-tumor.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/run62DU0AAXX_8-tumor.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/run62DU6AAXX_8-tumor.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/run62DU6AAXX_8-tumor.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/run62DUUAAXX_8-tumor.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/run62DUUAAXX_8-tumor.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/run62DUYAAXX_7-tumor.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/run62DUYAAXX_7-tumor.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/run62DVGAAXX_1-normal.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/run62DVGAAXX_1-normal.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/run62DVMAAXX_4-tumor.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/run62DVMAAXX_4-tumor.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/run62DVMAAXX_5-tumor.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/run62DVMAAXX_5-tumor.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/run62DVMAAXX_6-tumor.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/run62DVMAAXX_6-tumor.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/run62DVMAAXX_7-tumor.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/run62DVMAAXX_7-tumor.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/run62DVMAAXX_8-tumor.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/run62DVMAAXX_8-tumor.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/run62JREAAXX_3-tumor.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/run62JREAAXX_3-tumor.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/run62JREAAXX_4-tumor.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/run62JREAAXX_4-tumor.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/run62JREAAXX_5-tumor.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/run62JREAAXX_5-tumor.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/run62JREAAXX_6-tumor.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/run62JREAAXX_6-tumor.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/run62JREAAXX_7-tumor.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/run62JREAAXX_7-tumor.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/run62JREAAXX_8-tumor.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/run62JREAAXX_8-tumor.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/run62MK3AAXX_5-normal.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/run62MK3AAXX_5-normal.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/runA81DF6ABXX_1-normal.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/runA81DF6ABXX_1-normal.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/runA81DF6ABXX_2-normal.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/runA81DF6ABXX_2-normal.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/runAC0756ACXX_4-tumor.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/runAC0756ACXX_4-tumor.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/runAC0756ACXX_5-tumor.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/runAC0756ACXX_5-tumor.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/runAD08C1ACXX_1-tumor.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/runAD08C1ACXX_1-tumor.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/runBC04D4ACXX_2-normal.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/runBC04D4ACXX_2-normal.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/runBC04D4ACXX_3-normal.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/runBC04D4ACXX_3-normal.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/runBD06UFACXX_4-normal.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/runBD06UFACXX_4-normal.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/runBD06UFACXX_5-normal.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/runBD06UFACXX_5-normal.64.pair2.fastq.gz
>      > {{zenodo_url}}/files/runBD08K8ACXX_1-tumor.64.pair1.fastq.gz
>      > {{zenodo_url}}/files/runBD08K8ACXX_1-tumor.64.pair2.fastq.gz
>      >
>      > ```
>      {: .details }
>
>      {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>      {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Import adapters file in the same way you uploaded the sample data
>
>    ```
>    {{zenodo_url}}/files/adapters.fa
>    ```
>
{: .hands_on}




