---
title: Is there a way to filter on the Kalimari database? E.g. leaving out milk bacteria only to detect spoilers or contaminants, but the Kalimera list contains a lot more than that
box_type: question
layout: faq
contributors: [EngyNasr]
---

There are of 2 approaches:
1. You can look at a publication etc. to find a list of bacteria to remove. Then you would change the regex `^.*Gallus|Homo|Bos.*$`  to `^.*Gallus|Homo|Bos|Bacterium1|Bacterium2...|BacteriumN.*$`
2. Milk pathogens are somewhat known, salmonella, escherichia ...It might be easier to retain reads only mapping to pathogens instead

