---
title: "CheckM2 vs CheckM"
box_type: details
layout: faq
contributors: [bebatut]
---

CheckM2 ({%cite Chklovski2023CheckM2%}) is the successor of CheckM, but CheckM is still widely used, since its marker-based logic can be more interpretable in a biological sense. E.g., to date (2025-11-21), NCBI still allows submitting MAGs to GenBank if either checkM or checkM2 has a completeness of > 90% (see the [NCBI WGS/MAG submission guidelines](https://www.ncbi.nlm.nih.gov/genbank/wgsfaq/#metagen)).

**Key differences compared to CheckM**:

* CheckM relies primarily on lineage-specific single-copy marker genes to estimate completeness and contamination of microbial genomes.
* CheckM2 uses a machine-learning (gradient boost / ML) approach trained on simulated and experimental genomes, and does *not* strictly require a well-represented lineage in its marker database. 
* CheckM2 is reported to be more accurate and faster for both bacterial and archaeal lineages, especially when dealing with novel or very reduced-genome lineages (e.g., candidate phyla, CPR/DPANN) where classical marker-gene methods may struggle. 
* The database of CheckM2 can be updated more rapidly with new high-quality reference genomes, which supports scalability and improved performance over time.

If you're working with MAGs from underrepresented taxa (novel lineages) or very small genomes (streamlined bacteria/archaea), CheckM2 tends to give more reliable estimates of completeness/contamination. For more "standard" microbial genomes from well-studied taxa, CheckM may still work well, but you may benefit from the improved performance with CheckM2.