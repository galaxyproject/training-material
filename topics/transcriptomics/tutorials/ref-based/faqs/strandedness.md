---
title:  In 'infer experiments' I get unequal numbers, but in the IGV it looks like it is unstranded. What does this mean?
area: estimation of strandedness
box_type: question
layout: faq
contributors: [bebatut]
---

It’s also often the case that elimination of the second strand is not perfect, and there are genuine cases of bidirectional transcription in the genome. 70 / 30 % as in your report is not a good result for a stranded library. You can treat this as a stranded library in your analysis, but for instance you couldn’t make the conclusion that a given gene is actually transcribed from the reverse strand. Likely that the library preparation didn’t work perfectly. This can depend on many factors, one is that you need to completely digest your DNA using a high quality DNase before doing the reverse transcription.
