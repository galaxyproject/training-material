---
title: Añadir una etiqueta
description: Las etiquetas pueden ayudarte a organizar mejor tu historial y a hacer un seguimiento de los conjuntos de datos.
area: datasets
layout: faq
box_type: tip
contributions:
  authorship:
    - bebatut
    - wm75
    - hexylena
    - shiltemann
    - nekrut
  translation:
    - hvelab
    - unode
  funding:
    - biont
---


Los conjuntos de datos se pueden etiquetar. Esto simplifica el seguimiento de los conjuntos de datos a través de la interfaz de Galaxy. Las etiquetas pueden contener cualquier combinación de letras o números, pero no pueden contener espacios.

**Para etiquetar un conjunto de datos**:

1. Haga clic en el conjunto de datos para expandirlo
2. Haz click en **Añadir Etiquetas** {% icon galaxy-tags %}
3. Añade {% if include.tag %} una etiqueta llamada `{{include.tag}}`{% else %} tag text{% endif %}. Las etiquetas que empiecen por `#` se propagarán automáticamente a los resultados de las herramientas que utilicen este conjunto de datos (véase más abajo).
4. Pulse <kbd>Intro</kbd>
5. Compruebe que la etiqueta aparece debajo del nombre del conjunto de datos

**¡Las etiquetas que empiezan por `#` son especiales!**

Se llaman **etiquetas de nombre**. La característica única de estas etiquetas es que se *propagan*: si un conjunto de datos se etiqueta con una etiqueta de nombre, todos los derivados (hijos) de este conjunto de datos heredarán automáticamente esta etiqueta (véase más abajo). La figura siguiente explica por qué es tan útil. Considere el siguiente análisis (los números entre paréntesis corresponden a los números de los conjuntos de datos en la siguiente figura):

1. un conjunto de lecturas hacia adelante y hacia atrás (conjuntos de datos 1 y 2) se mapea contra una referencia usando {% tool Bowtie2 %} generando el conjunto de datos 3;
1. dataset 3 se utiliza para calcular la cobertura de lectura utilizando {% tool BedTools Genome Coverage %} *por separado* para las cadenas `+` y `-`. Esto genera dos conjuntos de datos (4 y 5 para más y menos, respectivamente);
1. los conjuntos de datos 4 y 5 se utilizan como entrada para los conjuntos de datos {% tool Macs2 broadCall %} que generan los conjuntos de datos 6 y 8;
1. los conjuntos de datos 6 y 8 son intersecados con las coordenadas de los genes (conjunto de datos 9) usando {% tool BedTools Intersect %} generando los conjuntos de datos 10 y 11.

![Una historia sin etiquetas de nombre versus una historia con etiquetas de nombre]({% link shared/images/histories_why_nametags.svg %})

Ahora considere que este análisis se realiza sin etiquetas de nombre. Esto se muestra en el lado izquierdo de la figura. Es difícil determinar qué conjuntos de datos contienen datos "positivos" y cuáles "negativos". Por ejemplo, ¿el conjunto de datos 10 contiene datos "más" o datos "menos"? Probablemente "menos", pero ¿está seguro? En el caso de un historial pequeño, como el que se muestra aquí, es posible rastrearlo manualmente, pero a medida que aumenta el tamaño del historial se convierte en todo un reto.

La parte derecha de la figura muestra exactamente el mismo análisis, pero utilizando etiquetas de nombre. Cuando se realizó el análisis, los conjuntos de datos 4 y 5 estaban etiquetados con `#plus` y `#minus`, respectivamente. Cuando se utilizaron como entradas para Macs2, los conjuntos de datos resultantes 6 y 8 las heredaron automáticamente, y así sucesivamente... Por lo tanto, es fácil rastrear las dos ramas (positiva y negativa) de este análisis.

Más información en un [tutorial #nametag dedicado]({% link topics/galaxy-interface/tutorials/name-tags/tutorial.md %}).


<!-- Image is here = https://docs.google.com/drawings/d/1iiNsau6ddiE2MV9qMyekUq2mrpDHHcc02bXtcFEAnhY/edit?usp=sharing -->


