---
title: "Punteggi di qualità"
area: format
box_type: details
layout: faq
contributions:
  authorship:
    - bebatut
    - nakucher
    - hexylena
  translation:
    - lisanna
    - silviadg87
    - unode
  funding:
    - biont
---


Ma cosa significa questo punteggio di qualità?

Il punteggio di qualità per ogni sequenza è una stringa di caratteri, uno per ogni base della sequenza nucleotidica, utilizzata per caratterizzare la probabilità di errata identificazione di ogni base. Il punteggio è codificato utilizzando la tabella dei caratteri ASCII (con [alcune differenze storiche](https://en.wikipedia.org/wiki/FASTQ_format#Encoding)):

Per risparmiare spazio, il sequenziatore registra un [carattere ASCII](http://drive5.com/usearch/manual/quality_score.html) per rappresentare i punteggi da 0 a 42. Ad esempio, 10 corrisponde a "+". Ad esempio, 10 corrisponde a "+" e 40 a "I". FastQC sa come tradurlo. Questo viene spesso chiamato punteggio "Phred".

![Codifica del punteggio di qualità con caratteri ASCII per diverse codifiche Phred. La sequenza di codici ascii è mostrata in alto con i simboli da 33 a 64, le lettere maiuscole, altri simboli e poi le lettere minuscole. Sanger mappa da 33 a 73, mentre Solexa è spostato, partendo da 59 e arrivando a 104. Illumina 1.3 inizia a 54 e arriva a 104, Illumina 1.5 è spostato di tre posizioni a destra ma termina comunque a 104. Illumina 1.8+ torna al Sanger, tranne che per un singolo punteggio più ampio. Illumina]({{site.baseurl}}/topics/sequence-analysis/faqs/images/fastq-quality-encoding.png)

Quindi a ogni nucleotide è associato un carattere ASCII che rappresenta il suo [punteggio di qualità Phred](https://en.wikipedia.org/wiki/Phred_quality_score), la probabilità di una chiamata di base errata:

| Phred Quality Score | Probability of incorrect base call | Base call accuracy |
| ------------------- | ---------------------------------- | ------------------ |
| 10                  | 1 in 10                            | 90%                |
| 20                  | 1 in 100                           | 99%                |
| 30                  | 1 in 1000                          | 99.9%              |
| 40                  | 1 in 10,000                        | 99.99%             |
| 50                  | 1 in 100,000                       | 99.999%            |
| 60                  | 1 in 1,000,000                     | 99.9999%           |


Cosa rappresenta 0-42? Questi numeri, se inseriti in una formula, ci dicono la probabilità di errore per quella base. Questa è la formula, dove Q è il nostro punteggio di qualità (0-42) e P è la probabilità di errore:

```
Q = -10 log10(P)
```

Utilizzando questa formula, possiamo calcolare che un punteggio di qualità di 40 significa solo 0,00010 probabilità di errore!

