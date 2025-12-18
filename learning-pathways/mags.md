---
layout: learning-pathway  # (uncomment this line to activate it)
type: use
cover-image: assets/images/microgalaxy-logo.png
cover-image-alt: "microgalaxy logo"

editorial_board:
- bebatut
- paulzierep
funding:
- elixir-europe
- deNBI

title: Building and Annotating Metagenome-Assembled Genomes (MAGs) from Metagenomics Reads
description: |
  This learning path will guide you through the process of constructing and analyzing **Metagenome-Assembled Genomes (MAGs)** using the Galaxy platform. You will explore the key steps involved in transforming raw metagenomic data into high-quality MAGs, from preprocessing to functional annotation.

  By the end of this path, you will be able to:
  - **List and describe** the essential steps in MAGs construction, including quality control, assembly, binning, and refinement.
  - **Define** core concepts such as MAGs, binning, and functional annotation, and understand their significance in metagenomic analysis.
  - **Explain** the importance of preprocessing metagenomic reads, focusing on quality control and contamination removal.
  - **Compare** the quality of MAGs using metrics like completeness and contamination, and assess their suitability for downstream analysis.
  - **Evaluate** the reliability of taxonomic assignments and functional annotations by leveraging reference databases.
  - **Analyze** the relative abundance of microbial taxa in samples and infer ecological dynamics.
  - **Identify** genomic features annotated by tools like Bakta, including coding sequences (CDS), rRNA, and tRNA.
  - **Interpret** functional annotation results to uncover metabolic pathways, virulence factors, and other biological roles within microbial communities.

  This path is designed to equip you with both the theoretical knowledge and practical skills needed to confidently construct, evaluate, and analyze MAGs in your research.

tags: [microbiome]

pathway:
  - section: "Module 0: Introduction to Galaxy – Navigating the Platform and Performing Your First Analysis"
    description: |
        Before diving into metagenomics, it's essential to become comfortable with the tools you'll be using. This module is designed to introduce you to the **Galaxy platform**—a user-friendly, web-based environment for bioinformatics analysis.

        Through a combination of **video tutorials and hands-on exercises**, you will:
        - Familiarize yourself with the **Galaxy interface**, including its key features and navigation.
        - Learn how to **import data**, organize your workspace, and use basic tools.
        - Complete a **guided, simple analysis** to gain confidence in running workflows and interpreting results.

        By the end of this module, you'll be ready to tackle more advanced analyses in subsequent modules.
    tutorials:
      - name: galaxy-intro-short
        topic: introduction
      - name: galaxy-intro-101
        topic: introduction

  - section: "Module 1: Quality Control – Ensuring High-Quality Metagenomic Data"
    description: |
        High-quality data is the foundation of reliable metagenomic analysis. **Poor-quality reads**—whether due to low base-calling accuracy, adapter contamination, or insufficient length—can introduce errors, bias assemblies, and compromise your results.

        In this module, you will:
        - Understand the **importance of quality control** in metagenomic workflows and its impact on downstream analyses.
        - Learn how to **assess, trim, and filter raw sequencing data** to retain only high-quality reads.
        - Use Galaxy tools to **remove contaminants, trim adapters, and filter low-quality sequences**, ensuring your data is clean and ready for further analysis.

        By the end of this module, you'll be equipped to confidently prepare your metagenomic data for assembly and other advanced analyses.
    tutorials:
      - name: quality-control
        topic: sequence-analysis

  - section: "Module 2: Contamination and Host Reads Removal – Purifying Your Metagenomic Dataset"
    description: |
        Metagenomic datasets frequently contain **non-microbial sequences**, such as **host DNA** (e.g., from honey bees) or **external contaminants** (e.g., human DNA introduced during sample handling or sequencing). These sequences can distort downstream analyses, leading to **misassemblies, incorrect taxonomic assignments, and biased functional interpretations**.

        In this module, you will:
        - Recognize the **sources and impacts** of contamination in metagenomic datasets.
        - Learn how to **identify and filter out host and contaminant sequences** using Galaxy tools.
        - Ensure your dataset is **enriched for microbial reads**, improving the accuracy of MAG reconstruction and enabling more reliable biological insights.

        By the end of this module, you'll be able to confidently clean your metagenomic data, setting the stage for high-quality MAG assembly and analysis.
    tutorials:
      - name: host-removal
        topic: microbiome

  - section: "Module 3: Assembly – Reconstructing and Assessing Contigs from Metagenomic Reads"
    description: |
        The foundation of MAG reconstruction lies in **assembly**—the computational process of piecing together fragmented metagenomic reads into longer genomic sequences called **contigs**. Think of it as solving a complex jigsaw puzzle: your goal is to identify reads that "fit together" by detecting overlapping sequences.

        In this module, you will:
        - Understand the **principles and challenges** of metagenomic assembly.
        - Learn how to use Galaxy tools to **assemble high-quality contigs** from your preprocessed reads.
        - Explore strategies to **optimize assembly parameters** for improved accuracy and completeness.
        - **Assess the quality of your assembly** using metrics such as contig length distribution, N50, and coverage, ensuring your contigs are suitable for downstream analysis.

        By the end of this module, you'll be equipped to transform your cleaned metagenomic data into contiguous sequences and evaluate their quality, setting the stage for successful MAG reconstruction.
    tutorials:
      - name: metagenomics-assembly
        topic: microbiome

  - section: "Module 4: Binning – From Contigs to Refined Microbial Genomes"
    description:  |
        Metagenomic **binning** is the process of grouping assembled contigs into discrete **bins**, each representing a potential microbial genome. By analyzing **sequence composition, coverage, and similarity**, binning allows researchers to reconstruct individual genomes from complex microbial communities.

        However, initial bins often contain **fragmented, redundant, or contaminated sequences**, which can compromise downstream analyses. To address this, **bin refinement** and **de-replication** are essential steps to improve the quality, completeness, and non-redundancy of your Metagenome-Assembled Genomes (MAGs).

        In this module, you will:
        - Understand how **binning algorithms** classify contigs into bins based on genomic signatures.
        - Use Galaxy tools to **perform binning** and assign sequences to their likely microbial origins.
        - **Evaluate bin quality** using metrics such as completeness, contamination, and strain heterogeneity.
        - Learn techniques for **refining bins**, including merging, splitting, and contamination removal.
        - Explore **de-replication** to identify and retain only the highest-quality representative MAG from sets of similar genomes.

        By the end of this module, you'll be able to reconstruct, refine, and validate high-quality MAGs, ensuring they are ready for taxonomic and functional analysis.
    tutorials:
      - name: metagenomics-binning
        topic: microbiome

#  - section: "Module 5: Taxonomic Assignment of MAGs – Identifying and Classifying Microbial Genomes"
#    description: |
#        Taxonomic assignment is a **pivotal step** in metagenomic analysis, allowing researchers to **identify and classify** microbial genomes based on their genetic content. By assigning taxonomic labels to your Metagenome-Assembled Genomes (MAGs), you can uncover the **diversity of microbial communities**, explore their ecological roles, and compare populations across different environments.
#
#        In this module, you will:
#        - Learn how **taxonomic assignment tools** use reference databases to classify MAGs at various taxonomic levels (e.g., phylum, genus, species).
#        - Use Galaxy to **assign taxonomic labels** to your MAGs and interpret the results.
#        - Explore how taxonomic information can reveal **microbial community structure** and ecological dynamics.
#        - Assess the **reliability of taxonomic assignments** and understand the limitations of reference databases.
#
#        By the end of this module, you'll be able to confidently assign taxonomic identities to your MAGs, enabling deeper insights into microbial diversity and function.
#    tutorials:
#      - name: 
#        topic: microbiome

  - section: "Module 6: Functional Annotation of MAGs – Applying Genomic Approaches to Metagenome-Assembled Genomes"
    description: |
        Functional annotation is a **fundamental process** in genomic analysis, whether you're working with microbial isolates or Metagenome-Assembled Genomes (MAGs). By applying the same robust approaches used for isolates, you can **identify and characterize genes** in MAGs, revealing their roles in **metabolic pathways, environmental interactions, and ecological functions**.

        In this module, you will:
        - Learn how **functional annotation tools** (such as Bakta) predict and classify genomic features, including coding sequences (CDS), rRNA, tRNA, and more.
        - Use Galaxy to **annotate your MAGs**, uncovering their biological potential and functional capabilities.
        - Explore **antimicrobial resistance (AMR) gene detection** as part of functional annotation, identifying genes associated with resistance mechanisms.
        - Interpret the **ecological and functional roles** of your MAGs, including genes linked to **pathogenicity, nutrient cycling, symbiosis, and AMR**.
        - Assess the **reliability of annotations** and understand the importance of reference databases in ensuring accurate predictions.

        By the end of this module, you'll be able to analyze MAGs with the same confidence and precision as microbial isolates, gaining deeper insights into their ecological roles and functional potential—including their resistance profiles.
    tutorials:
      - name: bacterial-genome-annotation
        topic: genome-annotation
      - name: amr-gene-detection
        topic: genome-annotation


  - section: "Recommended follow-up tutorial"
    description: |
        While this learning pathway goes into detail about all steps required for MAGs generation, a full end-to-end training was developed, which executes all workflows for MAGs generation and annotation instead of running individual steps. This tutorial can be used as a guide to execute the workflows on your own data. It explains how you can perform quality control, remove host contamination, and run the main MAGs workflow to obtain quality-controlled, taxonomically and functionally annotated MAGs. 
    tutorials:
      - name: mags-building
        topic: microbiome
---


