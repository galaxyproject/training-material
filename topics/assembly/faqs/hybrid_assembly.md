---
title: Most tools seem to have options for assembly using long and short reads, what are the pros and cons of the different tools?
layout: faq
box_type: question
contributors: [annasyme,slugger70]
---

In our experience, when both long and short reads are allowed as input, the difference comes down to the order in which set is assembled first. For example, Unicycler assembles the short reads first (which can be good, because they are more accurate), and then scaffolds these into larger contigs using long reads. Other tools (or workflows) often assemble long reads first (which can also be good because these can span repeat regions), then correct this assembly with information from the more accurate short reads. There may also be other variations on long/short read assembly, and/or iterations of these types of steps (assemble, correct). My preference is to assemble long reads first, but that's because I'm really interested in covering repeat regions. If accuracy was the aim, rather than contig length, the short-reads-first approach may be better. For even more complexity ... I think some tools now allow input of "trusted contigs" - i.e. contigs assembled from other tools. Ryan Wick has a new tool called [Trycyler](https://github.com/rrwick/Trycycler/wiki ) that can take in multiple assemblies to make a consensus (bacterial genomes).
