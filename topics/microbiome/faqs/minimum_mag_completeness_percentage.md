---
title: "Minimum MAG completeness percentage"
box_type: details
layout: faq
contributors: [bebatut]
---

Choosing a **minimum genome completeness threshold** is a critical but complex decision in de-replication and bin refinement. There is a **trade-off between computational efficiency and genome quality**:
- **Lower completeness thresholds** allow more genomes to be included but reduce the accuracy of similarity comparisons.
- **Higher completeness thresholds** improve accuracy but may exclude valuable genomes.


**Impact of Genome Completeness on Aligned Fractions**

When genomes are incomplete, the **aligned fraction**—the proportion of the genome that can be compared—decreases. For example, if you randomly sample **20% of a genome** twice, the aligned fraction between these subsets will be low, even if they originate from the same genome.

This effect is illustrated below, where lower completeness thresholds result in a **wider range of aligned fractions**, reducing the reliability of similarity metrics like ANI.

![This scatter plot visualizes the relationship between alignment coverage (y-axis) and minimum aligned fraction (min_frac) (x-axis) for five distinct microbial genomes, each represented by a different color: blue (GCA_000390985.1_Entamoeba_faecalis), green (GCA_000492815.1_Klebsiella_oxytoca), red (GCA_000988385.1_ASM98838v1), purple (GCA_001472555.1_SMART_628), and yellow (GCA_000821625.1_18PV). Each dot corresponds to a comparison between genome subsets of varying completeness, showing how alignment coverage decreases as the minimum aligned fraction decreases. The plot highlights that lower completeness thresholds result in reduced alignment coverage, which can impact the accuracy of genome similarity metrics like ANI, especially in de-replication workflows. The clustering of points at higher min_frac values indicates more reliable comparisons, while scattered points at lower min_frac values reflect the challenges of comparing incomplete genomes.]({{site.baseurl}}/topics/microbiome/faqs/images/minimum_mag_completeness_percentage_figure_a.png "Five genomes subset to 10–100% completeness, showing how aligned fractions vary with completeness. Source: <a href='https://drep.readthedocs.io/en/latest/choosing_parameters.html#importance-of-genome-completeness'>dRep documentation</a>")


**Effect on Mash ANI**

Incomplete genomes also **artificially lower Mash ANI values**, even for identical genomes. As completeness decreases, the **reported Mash ANI drops**, even when comparing identical genomes.

![This scatter plot compares the Average Nucleotide Identity (ANI) (y-axis) against alignment coverage (x-axis) for three different ANI calculation algorithms: Mash (green dots), ANIm (red dots), and gANI (blue dots). The plot illustrates how ANI values vary with alignment coverage, showing that gANI and ANIm consistently achieve high ANI values (~1.00) even at lower alignment coverage, indicating robust performance across a range of genome completeness. In contrast, Mash exhibits a clear dependency on alignment coverage, with ANI values decreasing significantly as coverage drops below ~0.6, reflecting its sensitivity to genome incompleteness. This highlights the importance of choosing the appropriate ANI algorithm based on the completeness and quality of the genomes being compared.]({{site.baseurl}}/topics/microbiome/faqs/images/minimum_mag_completeness_percentage_figure_b.png "An identical E. coli genome is subset to fractions ranging from 10% - 100% and fractions are compared. When lower amounts of the genome align (due to incompleteness), Mash ANI is severely impacted. Source: <a href='https://drep.readthedocs.io/en/latest/choosing_parameters.html#importance-of-genome-completeness'>dRep documentation</a>")

This is problematic because **Mash is used for primary clustering** in tools like dRep. If identical genomes are split into different primary clusters due to low Mash ANI, they will **never be compared by the secondary algorithm** and thus won't be de-replicated.

**Practical Implications for De-Replication**

1. **Primary Clustering Thresholds:**
   
   If you set a **minimum completeness of 50%**, identical genomes subset to 50% may only achieve a **Mash ANI of \\( \approx \\) 96%**. To ensure these genomes are grouped in the same primary cluster, the **primary clustering threshold must be \\( \leq \\) 96%**. Otherwise, they may be incorrectly separated.

2. **Computational Trade-Offs:**
   
   **Lower primary thresholds** increase the size of primary clusters, leading to **more secondary comparisons** and **longer runtimes**. **Higher thresholds** improve speed but risk missing true matches.

3. **Unknown Completeness:**

   In practice, the **true completeness of genomes is often unknown**. Tools like **CheckM** estimate completeness using single-copy gene inventories, but these estimates are not perfect in particular for **phages and plasmids**, explaining why they are not supported in dRep. In general though, checkM is pretty good at accessing genome completeness:

   ![This scatter plot illustrates the relationship between true genome subset percentages (x-axis) and completeness as estimated by CheckM (y-axis) for five distinct microbial genomes, each represented by a different color: blue (GCA_000390985.1_Entamoeba_faecalis), green (GCA_000492815.1_Klebsiella_oxytoca), red (GCA_000988385.1_ASM98838v1), purple (GCA_001472555.1_SMART_628), and yellow (GCA_000821625.1_18PV). Each point corresponds to a subset of the genome, showing how CheckM's completeness estimates align closely with the true subset percentages. The plot demonstrates that CheckM accurately predicts genome completeness, with points clustering near the diagonal line, indicating a strong correlation between the true and estimated completeness values. This highlights CheckM's reliability for assessing genome completeness in metagenomic studies.]({{site.baseurl}}/topics/microbiome/faqs/images/minimum_mag_completeness_percentage_figure_c.png "Source: <a href='https://drep.readthedocs.io/en/latest/choosing_parameters.html#importance-of-genome-completeness'>dRep documentation</a>")


**Guidelines for Setting Completeness Thresholds**

- **Avoid thresholds below 50% completeness**: Genomes below this threshold are often too fragmented for reliable comparisons, and secondary algorithms may fail.
- **Adjust Mash ANI thresholds accordingly**: If you lower the **secondary ANI threshold**, also lower the **Mash ANI threshold** to ensure incomplete but similar genomes are grouped together.

Balancing **genome completeness** and **computational efficiency** is key to effective de-replication. While lower completeness thresholds include more genomes, they **reduce alignment accuracy** and **increase runtime**. Aim for a **minimum completeness of \\( \geq \\)50%** and adjust clustering thresholds to avoid splitting identical genomes.