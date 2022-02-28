---
title: Why don't the aligned read files have quality scores?
box_type: question
layout: faq
contributors: [jontrow,RareSeas]
---

Quality scores take up the majority of space in our compressed sequence files, so removing them makes the files much smaller (~80% or more). In addition, many uses don't require per-base quality scores to successfully complete their work (some pipelines even require fastq format but don't actually use the quality scores), so these files represent a faster route to completing many analyses. The full quality scores are still available in the original SRA Runs for anyone that requires them, using the SRA Tools available in Galaxy.

