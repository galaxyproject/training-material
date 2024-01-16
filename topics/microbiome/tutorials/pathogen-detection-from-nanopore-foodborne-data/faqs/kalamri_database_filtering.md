---
title: Is there a way to filter on the Kalimari database? 
box_type: question
layout: faq
contributors: [EngyNasr]
redirect_from:
- /topics/metagenomics/tutorials/pathogen-detection-from-nanopore-foodborne-data/faqs/kalamri_database_filtering
---

To filter the Kalamari database, e.g. leaving out milk bacteria only to detect spoilers or contaminants, but the Kalimera list contains a lot more than that, you can:
1. Look at a publication etc. to find a list of bacteria to remove. 
2. Change the regex `^.*Gallus|Homo|Bos.*$`  to `^.*Gallus|Homo|Bos|Bacterium1|Bacterium2...|BacteriumN.*$`

Milk pathogens are somewhat known, *Salmonella*, *Escherichia*... It might be easier to retain reads only mapping to pathogens instead

