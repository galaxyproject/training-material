---
layout: tutorial_hands_on
title: Creating an Official Gene Set
zenodo_link: ""
questions:
- I have several genomes assemblies that are not annotated (or I do not trust annotations)
- I am interested to compare structure of a particular gene across these genome assemblies
- How do I do that?
objectives:
- Validate your genes and create an official gene set from them.
time_estimation: 30M
key_points:
- Genoest' OGS Tools makes it easy to validate and release a build of your genome and its annotations.
contributions:
  authorship:
  - abretaud
  funding:
  - gallantries

draft: true
---

Here's how we manage annotations at BIPAA:

- Automatic annotation with Maker
- Apollo server for manual curation, following our guidelines
- Each night, our Apollo report application checks that each annotated gene follows our guidelines. A report is available for each user with a list of errors and warnings that needs to be fixed, and a list of valid gene.
- Regularly we generate a new OGS by merging the automatic annotation with Apollo annotated genes (the valid ones only). We use the ogs_merge script to do that.
- To submit the OGS to ENA @ EBI in EMBL format, we (would like to) use gff2embl



> <agenda-title></agenda-title>
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}


