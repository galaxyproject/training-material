---
title: Why does my assembly graph in Bandage look different to the one pictured in the tutorial?
layout: faq
box_type: question
contributors: [annasyme,slugger70]
---

The assembly process in Flye is heuristic, and the resulting assembly will not necessarily be exactly the same each time. This may happen even if running the same data with the same version of Flye. It can also happen with a different version of Flye.

To make things more complicated (stop reading now if you would like!)... the chloroplast genome has a structure that includes repeats (the inverted repeats), and, the small-single-copy region of the chloroplast exists in two orientations between these repeats. So, sometimes the assembly will be a perfect circle, sometimes the inverted repeats will be collapsed into one piece, and sometimes the small-single-copy region will be attached ambiguously.  To make things even more complicated...the chloroplast genome may even be a dynamic structure, due to flip flop recombination.

For more see [this article](https://academic.oup.com/gbe/article/11/12/3372/5637229)



