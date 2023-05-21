---
layout: news
title: "New Galaxy training: Genome-wide alternative splicing analysis"
contributions:
  authorship: [gallardoalba]

cover: news/images/isoform_usage_post.jpg
coveralt: "Schematic of an isoform switch and detection pipeline. Data is annotated and a prediction is made for isoform switch consequences"
tags: [new tutorial]
tutorial: topics/transcriptomics/tutorials/differential-isoform-expression/tutorial.html
---

The GTN hosts a new training for analyzing alternative splicing at genome-wide scale!

# New Galaxy training: Genome-wide alternative splicing analysis

Isoform switching is a biological phenomenon that plays a crucial role in the regulation and diversity of gene expression. It refers to the process by which alternative splicing of pre-mRNA generates different mRNA isoforms that can produce distinct protein products under different conditions. Despite the fact that alternative splicing has been demonstrated to play a fundamental roles in a wide range of biological contexts, such as tissue-specific expression, disease development and evolutionary adaptation, isoform-based expression anaysis pipelines are still rare when compared with gene-level quantification studies.

This new Galaxy training explores the isoform switching concept and provides a step-by-step guide for analyzing this regulatory mechanism. It makes use of [IsoformSwitchAnalyzeR](https://bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html), a R package which integrates a multi-step analysis for isoform detection and functional prediction.

![Overview of computational pipeline for isoform detection and functional prediction.]({% link news/images/isoform_usage_post.jpg %} "Adapted from Vitting-Seerup et Sandelin, 2017.")

A second important aspect of this pipeline is that it allows to evaluate alternative splicing at two different levels: gene-specific and genome-wide scale. This dual-approach allows not only to indentify specific genes, but also prividing an analytical instrument for the study of genome-wide regulatory networks. Finally, a third important element of this training, directly derived from the operational capabilities of IsoformSwitchAnalyzeR, is that it allows to integrate multiple layers of information in the analysis. Concretely, in order to evaluate the functional consequences of the isoform switching events, IsoformSwitchAnalyzeR allows to integrate information bout coding-potential, protein-domains, signal peptides and instrinsically disordered regions.
