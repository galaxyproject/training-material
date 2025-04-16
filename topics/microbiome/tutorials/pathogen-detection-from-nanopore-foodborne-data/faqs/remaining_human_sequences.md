---
title: Isn't it awkward to find so many humans sequences there, since we filter for them before?
box_type: question
layout: faq
contributors: [EngyNasr]
redirect_from:
- /topics/metagenomics/tutorials/pathogen-detection-from-nanopore-foodborne-data/faqs/remaining_human_sequences
---

We see a lot that Kraken tends to assign many reads to human, despite they do not map to human genome. Due to resemblance between organisms and the limited species coverage of Kraken databases sometimes does happen that reads corresponding to higher organisms get mapped to humans. It was a very severe problem for the standard databases, because yeast genes were mis-assigned to human.
