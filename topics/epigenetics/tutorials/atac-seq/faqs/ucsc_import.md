---
title: "UCSC import: what should my file look like?"
box_type: question
layout: faq
contributors: [lldelisle,hrhotz,heylf]
---

~2020 lines, with the following header line:

```
bin    name    chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds        score   name2   cdsStartStat    cdsEndStat      exonFrames
```

Where:
- `txStart`: Transcript start site
- `cdsStart`: CodingSequence start site


**Note:** UCSC is updated frequently, you might get a slightly different number of lines. If you only get one row in this file, make sure you requested the entire chr22, not just one position.
