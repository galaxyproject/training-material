---
title:  Running more than one round of Pilon polishing
layout: faq
box_type: tip
contributors: [annasyme,slugger70]
---

Include the most recent polished assembly as input to the next round. You will also need to make a new bam file (here, we have `round1.bam` and `round2.bam`).

Round 1

```
assembly.fasta + illumina reads => BWA MEM => round1.bam
```
```
round1.bam + assembly.fasta => pilon => polished.fasta
```

Round 2

```
polished.fasta + illumina reads => BWA MEM => round2.bam
```
```
round2.bam + polished.fasta => pilon => polished2.fasta
```

**How to know when enough polishing iterations have run?**

There is no single answer, but a common way is to see when pilon stops making many polishing changes between rounds. So if round1 made 100 changes, and round2 made only 3, this seems like there would not be much more polishing to do.

**How can I see how many changes Pilon has made?**

There are two ways that I know of to see how many changes that Pilon made:

The first is to look at the tool *standard output* (`stdout`) from Pilon ([instructions]({% link faqs/galaxy/tools_logs.md %})).

Somewhere near the top of this log file will be a line that says how many corrections (changes) were made.

The second way is to count the number of lines in the changes file. To do this, use the tool called **Line/Word/Character count** {% icon tool %}, and select the line count option.
