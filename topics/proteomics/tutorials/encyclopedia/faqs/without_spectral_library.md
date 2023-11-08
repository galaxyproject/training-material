---
title: Can EncyclopeDIA be run on a DIA-MS dataset without a spectral library?
box_type: question
layout: faq
contributors: [emmaleith]
---

Yes. In this GTN, the workflow presented is the Standard EncyclopeDIA workflow; however, there is a variation upon the Standard EncyclopeDIA workflow, named the WALNUT EncyclopeDIA workflow in which a spectral library is not required. Simply, the WALNUT variation of the workflow omits the DLIB spectral/PROSIT library input, hence requiring just the GPF DIA dataset collection, Experimental DIA dataset collection, and the FASTA Protein Database file. Therefore, the Chromatogram Library is generated using the GPF DIA dataset collection and the FASTA Protein Database alone. This method does generate fewer searches than if a spectral library is used. The Galaxy-P team tested the efficacy of the WALNUT workflow compared to the Standard EncyclopeDIA workflow, and more information on that comparison and those results can be found at this [link](http://galaxyp.org/wp-content/uploads/2020/10/BCC_2020-Pratik.pdf).




