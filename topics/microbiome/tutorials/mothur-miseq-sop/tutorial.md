---
layout: tutorial_hands_on

title: "16S Microbial Analysis with mothur (extended)"
zenodo_link: "https://doi.org/10.5281/zenodo.800651"
priority: 1000
questions:
  - "What is the effect of normal variation in the gut microbiome on host health?"
objectives:
  - "Analyze of 16S rRNA sequencing data using the mothur toolsuite in Galaxy"
  - "Using a mock community to assess the error rate of your sequencing experiment"
  - "Visualize sample diversity using Krona and Phinch"
time_estimation: "6h"
key_points:
  - "16S rRNA gene sequencing analysis results depend on the many algorithms used and their settings"
  - "Quality control and cleaning of your data is a crucial step in order to obtain optimal results"
  - "Adding a mock community to serve as a control sample can help you asses the error rate of your experimental setup"
  - "We can explore alpha and beta diversities using Krona and Phinch for dynamic visualizations"
contributors:
  - shiltemann
  - bebatut
  - tnabtaf
subtopic: metabarcoding
tags:
  - metabarcoding
  - 16S
  - microgalaxy
redirect_from:
  - /topics/metagenomics/tutorials/mothur-miseq-sop/tutorial
edam_ontology:
- topic_3697 # Microbial ecology
- topic_0637 # Taxonomy
- topic_0080 # Sequence analysis
- topic_4038 # Metabarcoding
---

{% include topics/microbiome/tutorials/mothur-miseq-sop/content.md short=false %}
