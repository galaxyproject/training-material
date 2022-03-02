---
title: Add genome and annotations to IGV from Galaxy
box_type: tip
area: igv
layout: faq
contributors: [delphine-l,bebatut,hexylena,shiltemann]
---

1. Upload a FASTA file with the reference genome and a GFF3 file with its annotation in the history (if not already there)
1. Install [IGV](https://software.broadinstitute.org/software/igv/download) (if not already installed)
2. Launch IGV on your computer
2. Expand the FASTA dataset with the genome in the history
3. Click on the `local` in `display with IGV` to load the genome into the IGV browser
4. Wait until all **Dataset status** are `ok`
5. Close the window

   An alert `ERROR Parameter "file" is required` may appear. Ignore it.

2. Expand the GFF3 dataset with the annotations of the genome in the history
3. Click on the `local` in `display with IGV` to load the annotation into the IGV browser
5. Switch to the IGV instance

   The annotation track should appear. Be careful that all files have the same genome ID
