---
title: Why do we do dimension reduction and then clustering? Why not just cluster on the actual data?
area: gene
box_type: tip
layout: faq
---

The actual data has tens of thousands of genes, and so tens of thousands of variables to consider. Even after selecting for the most variable genes and the most high quality genes, we can still be left with > 1000 genes. Performing clustering on a dataset with 1000s of variables is possible, but computationally expensive. It is therefore better to perform dimension reduction to reduce the number of variables to a latent representation of these variables. These latent variables are ideally more than 10 but less than 50 to capture the variability in the data to perform clustering upon.
