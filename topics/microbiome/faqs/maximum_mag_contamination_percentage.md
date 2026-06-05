---
title: "Maximum MAG contamination percentage"
box_type: details
layout: faq
contributors: [bebatut]
---

**Contamination** in Metagenome-Assembled Genomes (MAGs) refers to the presence of sequences from organisms other than the target organism. During de-replication, setting a **maximum contamination threshold** ensures that only high-quality, representative genomes are retained for downstream analyses.

**Why Set a Maximum Contamination Threshold?**

1. **Data Quality:**
   
   High contamination can distort taxonomic classification, functional annotation, and comparative genomics. A contamination threshold ensures that only high-purity MAGs are included.

2. **Accurate De-Replication:**
   
   Contaminated MAGs may cluster incorrectly, leading to misrepresentation of microbial diversity. A contamination threshold helps ensure that only genuinely similar genomes are grouped together.

3. **Functional and Ecological Insights:**
   Low-contamination MAGs provide more reliable insights into microbial functions and ecological roles.


**Default and Common Contamination Thresholds**

The **default contamination threshold in dRep is 25%**, which balances inclusivity and quality for general metagenomic studies. However, the choice of threshold depends on the specific goals of your analysis:

| **Threshold** | **Use Case**                                                                                     |
|---------------|--------------------------------------------------------------------------------------------------|
| **< 5%**      | Ideal for high-confidence analyses, such as reference genomes or species-level comparisons.      |
| **< 10%**     | Suitable for most metagenomic studies, balancing purity and genome diversity.                     |
| **< 25%**     | Default in dRep, allowing for broader genome inclusion while maintaining reasonable purity.       |

- **< 5% Contamination:**
  Used for **high-quality MAGs** intended for reference databases or detailed functional analyses.

- **< 10% Contamination:**
  A widely adopted threshold for **general metagenomic studies**, balancing purity and genome retention.

- **< 25% Contamination:**
  The **default in dRep**, this threshold is more permissive, allowing for a broader range of MAGs while still maintaining reasonable quality.

**How to Choose the Right Threshold**

- **Study Goals:**
  For **high-quality reference databases**, use a **< 5% threshold**. For **general analyses**, **< 10%** is recommended. The **default 25%** in dRep is suitable for exploratory studies.

- **Tool Recommendations:**
  Tools like **dRep** and **CheckM** provide contamination estimates and allow you to set thresholds based on your needs.

- **Trade-Offs:**
  Stricter thresholds (e.g., < 5%) will **exclude more MAGs**, potentially reducing dataset diversity. More permissive thresholds (e.g., < 25%) will **retain more MAGs** but may include lower-quality genomes.

Setting a **maximum contamination threshold** is essential for ensuring the quality of de-replicated MAGs. The **default threshold in dRep is 25%**, but you can adjust it based on your study goals and the trade-offs between purity and genome retention. By carefully selecting this threshold, you can optimize your MAG dataset for accurate downstream analyses.