---
title: "New Tutorial: Single-cell ATAC-seq standard processing with SnapATAC2"
tags: [new tutorial, single-cell, epigenetics]
contributions:
  authorship: [timonschlegel]
cover: "topics/single-cell/images/scatac-standard-snapatac2/snapatac2-pipeline.png"
coveralt: "Analysis Pipeline with SnapATAC2"
tutorial: topics/single-cell/tutorials/scatac-standard-processing-snapatac2/tutorial.html

layout: news
---

We are proud to announce that a new training, explaining the analysis of single cell ATAC-seq data with **SnapATAC2** and **Scanpy**, is now available in the Galaxy Training Network. 

![SnapATAC2 pipeline]({% link topics/single-cell/images/scatac-standard-snapatac2/snapatac2-pipeline.png %} "SnapATAC2 standard processing pipeline. The single cell data is first preprocessed and a high-quality count matrix is produced. Next, the dimensionality of the data is reduced and clusters visualized. Finally, clusters are manually annotated with marker genes.")

The tutorial consists of three sections: Preprocessing, Dimensionality Reduction and Clustering, and Cluster Annotation. During the preprocessing, cell-by-feature count matrices are created, which are filtered to produce high-quality AnnData count matrices. After that, the dimensionality of the data is reduced through nonlinear spectral embedding. The data is then projected onto two-dimensional space and clusters are identified. To annotate the clusters, the gene activity for each cell is measured and the activity of marker genes is visualized with the **Scanpy** package. The activity profile of each cluster is then used to determine the correct cell type. 
