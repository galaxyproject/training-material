---
title: "Average Nucleotide Identity (ANI): A Measure of Genomic Similarity"
box_type: details
layout: faq
contributors: [bebatut]
---

**Average Nucleotide Identity (ANI)** is a computational metric used to quantify the **genomic similarity** between two microbial genomes. It measures the **mean sequence identity** across all **orthologous regions**—regions of the genome that are shared and aligned between the two organisms. ANI is widely used in microbial genomics and metagenomics to:
- **Compare the similarity of bacterial or archaeal genomes.**
- **Define species boundaries** (e.g., a threshold of **95–96% ANI** is commonly used to delineate bacterial species).
- **Identify redundant genomes** in datasets, such as Metagenome-Assembled Genomes (MAGs), for de-replication.

**How ANI Is Calculated**
1. **Genome Alignment:** The genomes of two organisms are compared using **whole-genome alignment tools** (e.g., BLAST, MUMmer, or FastANI). These tools identify regions of the genomes that are **homologous** (shared due to common ancestry).
2. **Identity Calculation:** For each aligned region, the **percentage of identical nucleotides** is calculated. ANI is then computed as the **mean identity** across all aligned regions that meet a specified length threshold (e.g., regions \\( \geq \\) 1,000 base pairs).
3. **Normalization:** ANI accounts for **unaligned regions** (e.g., due to genomic rearrangements or horizontal gene transfer) by focusing only on the aligned portions of the genomes. This ensures the metric reflects **conserved genomic similarity** rather than absolute sequence coverage.

**Interpreting ANI Values**
ANI values range from **0% to 100%**, where:
- **ANI \\( \approx \\) 100%:** The genomes are nearly identical, likely representing **strains of the same species** or **clones**.
- **ANI \\( \geq \\) 95–96%:** The genomes likely belong to the **same species** (a widely accepted threshold for bacterial species delineation).
- **ANI \\( \approx \\) 80–95%:** The genomes belong to **closely related species** within the same genus.
- **ANI \\( < \\) 80%:** The genomes are **distantly related** and likely belong to different genera or higher taxonomic ranks.

**Applications of ANI**
1. **Species Delineation:** ANI is a **gold standard** for defining microbial species, replacing or complementing traditional methods like **DNA-DNA hybridization (DDH)**. For example, an ANI \\(\geq \\) 95% is often used to confirm that two genomes belong to the same species.
2. **De-Replication of MAGs:** In metagenomics, ANI is used to **identify and remove redundant MAGs** from datasets. MAGs with ANI values above a threshold (e.g., **99%**) are considered redundant, and only the highest-quality representative is retained for downstream analysis.
3. **Strain-Level Comparisons:** ANI can distinguish between **closely related strains** of the same species, helping researchers study microbial diversity at fine taxonomic resolutions.
4. **Taxonomic Classification:** ANI is used to **assign unknown genomes** to known taxonomic groups by comparing them to reference genomes in databases.

**Tools for Calculating ANI**

Several bioinformatics tools are available to compute ANI, including:
- {% tool [FastANI](toolshed.g2.bx.psu.edu/repos/iuc/fastani/fastani/1.3) %} ({% cite Jain_2018 %}): A fast, alignment-free tool for estimating ANI between genomes.
- **[PyANI](https://github.com/widdowquinn/pyani):** A Python-based tool that supports multiple ANI calculation methods (e.g., BLAST, MUMmer).

**Limitations of ANI**

While ANI is a powerful metric, it has some limitations:
- **Dependence on Alignment:** ANI requires **sufficiently aligned regions** to be accurate. Highly divergent or rearranged genomes may yield unreliable ANI values.
- **Threshold Variability:** The ANI threshold for species delineation (e.g., 95%) may vary depending on the microbial group or study context.
- **Computational Requirements:** Calculating ANI for large datasets (e.g., thousands of MAGs) can be computationally intensive.


Average Nucleotide Identity (ANI) is a **fundamental metric** in microbial genomics, providing a robust way to compare genomes, define species, and refine metagenomic datasets. By leveraging ANI, researchers can ensure the **accuracy and reliability** of their genomic analyses, from species classification to de-replication of MAGs.