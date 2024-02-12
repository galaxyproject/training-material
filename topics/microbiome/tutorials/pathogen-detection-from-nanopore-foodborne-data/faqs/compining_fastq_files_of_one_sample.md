---
title: After sequencing with MinKNOW software, we get many fastq files, do these files need to be combined into one file before uploading or is it possible to upload them all at once?
box_type: question
layout: faq
contributors: [EngyNasr]
redirect_from:
- /topics/metagenomics/tutorials/pathogen-detection-from-nanopore-foodborne-data/faqs/compining_fastq_files_of_one_sample
---

After sequencing with MinKNOW software, it is a good approach to combine the files from the same run before processing them. You could create a collection per run with all fastq files and then use the collection operation to concatenate all files in a collection.
