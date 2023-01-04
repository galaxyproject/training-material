---
layout: tutorial_hands_on

title: "NCBI BLAST+ against the MAdLand"
zenodo_link: 'https://zenodo.org/record/4710649'
questions:
- How can we perform Blast analysis on Galaxy ?
- Why can be use MadLand DB for sequence comparison ? 

objectives:
- Test objective
time_estimation: "20m"
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- Deepti Varshney

---


# Introduction

<!-- This is a comment. -->

MAdLand is a collection of fully sequenced plant and algal genomes with a focus on non-seed plants and streptophyte algae. It includes, for comparison, genomes of fungi, animals, SAR group, bacteria and archaea. It is developed and maintained by the [Rensing lab](http://plantco.de). The species are abbreviated by a 5 letter code, which consists of the first three letters of the genus and the first two of the species name, e.g. CHABR for Chara braunii. We add the gene ID to that and additional shortcuts, like whether it is plastome encoded (pt) or transcriptome-based (tr, in cases when no genome is available yet).


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Register your account on Galaxy 

If you are a first-time user of Galaxy, follow the [tutorial](https://training.galaxyproject.org/training-material/faqs/galaxy/account_create.html) to create an account and start using the Galaxy workspace.

# How to upload data

> <hands-on-title> Data Upload </hands-on-title>
> 1. Create a new history for this tutorial and give it a proper name
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>    {% snippet faqs/galaxy/histories_rename.md %}


# NCBI Blast+ on Galaxy ##

After successfully logging in to the Galaxy server, Go to the [NCBI-Blast+](https://usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/devteam/ncbi_blast_plus/ncbi_blastp_wrapper/2.10.1+galaxy2) tool. For more details for BLAST analysis, we recommand you to follow the [Similarity-searches-blast](https://training.galaxyproject.org/training-material/topics/genome-annotation/tutorials/genome-annotation/tutorial.html#similarity-searches-blast) tutorial.


## Sub-step with **NCBI BLAST+ blastp or blastx**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NCBI BLAST+ blastp](toolshed.g2.bx.psu.edu/repos/devteam/ncbi_blast_plus/ncbi_blastp_wrapper/2.10.1+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Protein query sequence(s)"*: `output` (Input dataset)
>    - *"Subject database/sequences"*: `Locally installed BLAST database`
>        - *"Protein BLAST database"*: `MAdLandDB`

