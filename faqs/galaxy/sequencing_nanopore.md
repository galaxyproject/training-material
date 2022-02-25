---
title: "Nanopore sequencing"
area: sequencing
box_type: comment
layout: faq
contributors: [shiltemann]
---


Nanopore sequencing has several properties that make it well-suited for our purposes

1. Long-read sequencing technology offers **simplified** and less ambiguous genome **assembly**
2. Long-read sequencing gives the ability to **span repetitive genomic regions**
3. Long-read sequencing makes it possible to **identify large structural variations**

![How nanopore sequencing works]({% link topics/metagenomics/images/plasmid-metagenomics-nanopore/sequence_method.jpg %} "Using nanopore sequencing, a single molecule of DNA or RNA can be sequenced without the need for PCR amplification or chemical labeling of the sample. (Image from: <a href="https://nanoporetech.com/sites/default/files/s3/white-papers/WGS_Assembly_white_paper.pdf?submissionGuid=40a7546b-9e51-42e7-bde9-b5ddef3c3512">Nanopore sequencing: The advantages of long reads for genome assembly</a>)") <br><br>


When using Oxford Nanopore Technologies (ONT) sequencing, the change in
electrical current is measured over the membrane of a flow cell. When
nucleotides pass the pores in the flow cell the current change is translated
(basecalled) to nucleotides by a basecaller. A schematic overview is given in
the picture above.

When sequencing using a MinIT or MinION Mk1C, the basecalling software is
present on the devices. With basecalling the electrical signals are translated
to bases (A,T,G,C) with a quality score per base. The sequenced DNA strand will
be basecalled and this will form one read. Multiple reads will be stored in a
fastq file.

