---
title: Puntuación de calidad
area: format
box_type: details
layout: faq
contributions:
  authorship:
    - bebatut
    - nakucher
    - hexylena
  translation:
    - hvelab
    - unode
  funding:
    - biont
---


¿Pero qué significa esta puntuación de calidad?

La puntuación de calidad de cada secuencia es una cadena de caracteres, uno por cada base de la secuencia de nucleótidos, que se utiliza para caracterizar la probabilidad de identificación errónea de cada base. La puntuación se codifica utilizando la tabla de caracteres ASCII (con [algunas diferencias históricas](https://en.wikipedia.org/wiki/FASTQ_format#Encoding)):

Para ahorrar espacio, el secuenciador registra un [carácter ASCII](http://drive5.com/usearch/manual/quality_score.html) para representar las puntuaciones 0-42. Por ejemplo, 10 corresponde a "+" y 40 a "I". FastQC sabe cómo traducir esto. A menudo se denomina puntuación "Phred".

![Codificación de la puntuación de calidad con caracteres ASCII para diferentes codificaciones Phred. La secuencia de códigos ascii se muestra en la parte superior con símbolos para 33 a 64, letras mayúsculas, más símbolos y luego letras minúsculas. Sanger mapea de 33 a 73 mientras que solexa está desplazado, empezando en 59 y llegando hasta 104. Illumina 1.3 comienza en 54 y llega hasta 104, Illumina 1.5 se desplaza tres puntuaciones a la derecha, pero aún así termina en 104. Illumina 1.8+ se remonta al Sanger excepto una sola puntuación más ancha. Illumina]({{site.baseurl}}/topics/sequence-analysis/faqs/images/fastq-quality-encoding.png)

Así que hay un carácter ASCII asociado a cada nucleótido, que representa su [puntuación de calidad Phred](https://en.wikipedia.org/wiki/Phred_quality_score), la probabilidad de una llamada de base incorrecta:

| Puntuación de calidad Phred | Probabilidad de una llamada de base incorrecta | Precisión de llamada de base |
| ------------------- | ---------------------------------- | ------------------ |
| 10                  | 1 in 10                            | 90%                |
| 20                  | 1 in 100                           | 99%                |
| 30                  | 1 in 1000                          | 99.9%              |
| 40                  | 1 in 10,000                        | 99.99%             |
| 50                  | 1 in 100,000                       | 99.999%            |
| 60                  | 1 in 1,000,000                     | 99.9999%           |


¿Qué representa 0-42? Estos números, cuando se introducen en una fórmula, nos indican la probabilidad de error para esa base. Esta es la fórmula, donde Q es nuestra puntuación de calidad (0-42) y P es la probabilidad de error:

```
Q = -10 log10(P)
```

Utilizando esta fórmula, podemos calcular que una puntuación de calidad de 40 significa sólo 0,00010 de probabilidad de error

