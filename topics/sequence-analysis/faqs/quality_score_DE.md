---
title: "Qualitätswerte"
area: format
box_type: details
layout: faq
contributions:
  authorship: [bebatut, nakucher, hexylena]
  translation: [Tillsa, unode]
  funding: [biont]
---


Aber was bedeutet dieser Qualitätswert?

Der Qualitätsscore für jede Sequenz ist eine Zeichenkette, eine für jede Base der Nukleotidsequenz, um die Wahrscheinlichkeit einer falschen Identifizierung jeder Base zu charakterisieren. Der Score wird unter Verwendung der ASCII-Zeichentabelle (mit [einigen historischen Unterschieden](https://en.wikipedia.org/wiki/FASTQ_format#Encoding)) kodiert:

Um Platz zu sparen, zeichnet der Sequenzer ein [ASCII-Zeichen](http://drive5.com/usearch/manual/quality_score.html) auf, um die Punktzahlen 0-42 darzustellen. Zum Beispiel entspricht 10 einem "+" und 40 einem "I". FastQC weiß, wie dies zu übersetzen ist. Dies wird oft als "Phred"-Bewertung bezeichnet.

![Kodierung der Qualitätsbewertung mit ASCII-Zeichen für verschiedene Phred-Kodierungen. Oben ist die ASCII-Codefolge mit Symbolen für 33 bis 64, Großbuchstaben, weiteren Symbolen und dann Kleinbuchstaben dargestellt. Die Sanger-Kodierung reicht von 33 bis 73, während die Solexa-Kodierung verschoben ist und bei 59 beginnt und bis 104 geht. Illumina 1.3 beginnt bei 54 und geht bis 104, Illumina 1.5 ist um drei Ziffern nach rechts verschoben, endet aber immer noch bei 104. Illumina 1.8+ entspricht dem Sanger-Test, ist aber um einen einzigen Punkt breiter. Illumina]({% link topics/sequence-analysis/faqs/images/fastq-quality-encoding.png %})

Jedem Nukleotid ist also ein ASCII-Zeichen zugeordnet, das seinen [Phred-Qualitätsscore](https://en.wikipedia.org/wiki/Phred_quality_score), die Wahrscheinlichkeit eines falschen Basenrufs, darstellt:

| Phred Quality Score | Probability of incorrect base call | Base call accuracy |
| ------------------- | ---------------------------------- | ------------------ |
| 10                  | 1 in 10                            | 90%                |
| 20                  | 1 in 100                           | 99%                |
| 30                  | 1 in 1000                          | 99.9%              |
| 40                  | 1 in 10,000                        | 99.99%             |
| 50                  | 1 in 100,000                       | 99.999%            |
| 60                  | 1 in 1,000,000                     | 99.9999%           |


Was bedeutet 0-42? Wenn man diese Zahlen in eine Formel einsetzt, erhält man die Fehlerwahrscheinlichkeit für diese Basis. Die Formel lautet wie folgt: Q ist unser Qualitätswert (0-42) und P ist die Fehlerwahrscheinlichkeit:

```
Q = -10 log10(P)
```

Mit dieser Formel können wir berechnen, dass eine Qualitätsbewertung von 40 nur 0,00010 Wahrscheinlichkeit eines Fehlers bedeutet!

