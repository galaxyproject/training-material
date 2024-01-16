---
title: Can we polish the assembly with long reads too?
layout: faq
box_type: tip
contributors: [annasyme,slugger70]
---

Yes. In this tutorial, we only polish the assembly with the short reads. This may be enough for bacterial genomes. However, for an even better polish (usually), a common approach is to also polish the assembly with the long reads. A typical workflow for this would assemble with long reads, then polish with long reads (x 4 rounds, with Racon), polish with long reads again (x 1 round, with Medaka), then polish with short reads (x2 rounds with Pilon).



