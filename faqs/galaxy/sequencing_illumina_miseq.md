---
title: "Illumina MiSeq sequencing"
area: sequencing
box_type: comment
layout: faq
contributors: [shiltemann]
---


Illumina MiSeq sequencing is based on sequencing by synthesis. As the name
suggests, fluorescent labels are measured for every base that bind at a
specific moment at a specific place on a flow cell. These flow cells are
covered with oligos (small single strand DNA strands). In the library
preparation the DNA strands are cut into small DNA fragments (differs per
kit/device) and specific pieces of DNA (adapters) are added, which are
complementary to the oligos. Using bridge amplification large amounts of
clusters of these DNA fragments are made. The reverse string is washed away,
making the clusters single stranded. Fluorescent bases are added one by one,
which emit a specific light for different bases when added. This is happening
for whole clusters, so this light can be detected and this data is basecalled
(translation from light to a nucleotide) to a nucleotide sequence (Read). For
every base a quality score is determined and also saved per read. This
process is repeated for the reverse strand on the same place on the flow
cell, so the forward and reverse reads are from the same DNA strand. The
forward and reversed reads are linked together and should always be processed
together!

For more information watch this [video from Illumina](https://youtu.be/fCd6B5HRaZ8)

