---
title: Adding Custom Builds to database
area: datasets
description: Some tools and functions require that the 'database' attribute that is assigned. A Custom Reference Genome is set up as a Custom Build prior to use. Examples are the tools Featurecounts, Extract Genomic DNA, certain Picard tools, and the functions under Visualization.
box_type: tip
layout: faq
contributors: [jennaj, Nurzhamalyrys]
---


Once Custom Build is created, it is added to the list Database/Build: on the dataset 'Edit Attributes' and 'Upload File' tool forms and is available for 'Visualizations'. These can be assigned or used just like any other reference genome.

1. Start with an existing FASTA Custom Reference Genome in your history. It is very important make sure the format is correct.
2. Go to the top "User" menu and select "Custom Builds".
3. Enter in the labels (no spaces and no special characters other than "_").
4. Select the fasta Custom Reference Genome.
5. Submit and wait for the build to finish loading before assigning to a dataset or using to start a new Visualization.

Note: It is fine to navigate away from this form and come back to it later to check for status. The larger the fasta file and busier the Galaxy instance is, the longer the processing will take.
