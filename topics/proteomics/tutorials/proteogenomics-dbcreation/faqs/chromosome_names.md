---
title: Why do we change the chromosome names in the Ensembl GTF to match the UCSC genome reference?
box_type: question
layout: faq
contributors: [subinamehta]
---

UCSC chromosome names begin with the prefix `chr`, but Ensembl chromosome names do not. For example, chromosome 19 would be denoted as `chr19` in UCSC, and as `19` in Ensemble. Most tools would view those as different when looking for matches/overlaps. Therefore it is always a good idea to make sure these match before you perform any downstream analysis.

