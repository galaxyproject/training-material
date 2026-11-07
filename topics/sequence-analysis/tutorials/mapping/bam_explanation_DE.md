Eine BAM-Datei ([Binary Alignment Map](https://en.wikipedia.org/wiki/SAM_(file_format))) ist eine komprimierte Binärdatei, in der die Lesesequenzen gespeichert sind und in der angegeben ist, ob sie an eine Referenzsequenz (z. B. ein Chromosom) angeglichen wurden, und wenn ja, an welcher Position auf der Referenzsequenz sie angeglichen wurden.

> <hands-on-title>Inspektion einer BAM/SAM-Datei</hands-on-title>
> 
> 1. Untersuchen Sie die {% icon param-file %} Ausgabe von **{{ include.mapper }}** {% icon tool %}
> 
{: .hands_on}

Eine BAM-Datei (oder eine SAM-Datei, die nicht komprimierte Version) besteht aus:

- Ein Header-Abschnitt (die Zeilen, die mit `@` beginnen), der Metadaten enthält, insbesondere die Chromosomennamen und -längen (Zeilen, die mit dem Symbol `@SQ` beginnen)
- Ein Alignment-Abschnitt, bestehend aus einer Tabelle mit 11 Pflichtfeldern sowie einer variablen Anzahl von optionalen Feldern:

  | Col | Field | Type    | Brief Description                     |
  | --- | ----- | ------- | ------------------------------------- |
  | 1   | QNAME | String  | Query template NAME                   |
  | 2   | FLAG  | Integer | Bitwise FLAG                          |
  | 3   | RNAME | String  | References sequence NAME              |
  | 4   | POS   | Integer | 1- based leftmost mapping POSition    |
  | 5   | MAPQ  | Integer | MAPping Quality                       |
  | 6   | CIGAR | String  | CIGAR String                          |
  | 7   | RNEXT | String  | Ref. name of the mate/next read       |
  | 8   | PNEXT | Integer | Position of the mate/next read        |
  | 9   | TLEN  | Integer | Observed Template LENgth              |
  | 10  | SEQ   | String  | Segment SEQuence                      |
  | 11  | QUAL  | String  | ASCII of Phred-scaled base QUALity+33 |

> <question-title></question-title>
> 
> 1. Welche Informationen finden Sie in einer SAM/BAM-Datei?
> 2. Was sind die zusätzlichen Informationen im Vergleich zu einer FASTQ-Datei?
> 
> > <solution-title></solution-title>
> > 1. Sequenzen und Qualitätsinformationen, wie ein FASTQ
> > 2. Mapping-Informationen, Position des Read auf dem Chromosom, Mapping-Qualität, etc
> {: .solution }
{: .question}

