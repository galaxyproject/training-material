---
layout: tutorial_hands_on
title: Análisis de datos RNA-Seq basados en referencias
subtopic: introduction
priority: 2
tags:
- bulk
- rna-seq
- collections
- drosophila
- QC
- cyoa
level: Introductory
zenodo_link: https://zenodo.org/record/6457007
questions:
- ¿Cuáles son los pasos para procesar datos de RNA-Seq?
- ¿Cómo identificar genes diferencialmente expresados en múltiples condiciones experimentales?
-  ¿Qué funciones biológicas se ven afectadas por la expresión diferencial de genes?
objectives:
- Revisar un informe de calidad de secuencias generado por Falco/MultiQC para datos de RNA-Seq
- Explicar el principio y la especificidad del mapeo de datos de RNA-Seq a un genoma de referencia eucariota
- Seleccionar y ejecutar una herramienta de mapeo de última generación para datos de RNA-Seq
- Evaluar la calidad de los resultados del mapeo
- Describir el proceso para estimar la direccionalidad de la librería (strandness)
- Estimar el número de lecturas por gen
- Explicar la normalización de conteos necesaria antes de comparar muestras
- Construir y ejecutar un análisis de expresión génica diferencial
- Analizar la salida de DESeq2 para identificar, anotar y visualizar genes diferencialmente expresados
- Realizar un análisis de enriquecimiento de ontologías génicas (GO)
- Realizar y visualizar un análisis de enriquecimiento de rutas KEGG
time_estimation: 8h
key_points:
- Se debe usar una herramienta de mapeo con conciencia de empalmes para datos RNA‑Seq eucariotas
- Numerosos factores deben tenerse en cuenta al ejecutar un análisis de expresión génica diferencial
follow_up_training:
- type: internal
  topic_name: transcriptomics
  tutorials:
  - rna-seq-viz-with-heatmap2
  - rna-seq-viz-with-volcanoplot
  - rna-seq-genes-to-pathways
contributions:
  authorship:
  - bebatut
  - malloryfreeberg
  - moheydarian
  - erxleben
  - pavanvidem
  - blankclemens
  - mblue9
  - nsoranzo
  - pvanheus
  - lldelisle
  editing:
  - hexylena
  - clsiguret
  translation:
  - hvelab
  - unode
  funding:
  - deNBI
  - biont
recordings:
- youtube_id: AeiW3IItO_c
  speakers:
  - lldelisle
  captioners:
  - lldelisle
  date: '2023-05-15'
  galaxy_version: '23.01'
  length: 2H50M
  cyoa: true
- captioners:
  - hexylena
  - shiltemann
  date: '2021-02-15'
  galaxy_version: '21.01'
  length: 2H30M
  youtube_id: j4onRSN650A
  speakers:
  - bebatut
lang: es
tags:
  - deutsch
  - english
  - italiano
translations:
  - de
  - en
  - it
---



En los últimos años, la secuenciación del ARN (abreviada RNA-Seq) se ha convertido en una tecnología muy utilizada para analizar el transcriptoma celular en continuo cambio, es decir, el conjunto de todas las moléculas de ARN de una célula o de una población de células. Uno de los objetivos más comunes de RNA-Seq es el perfilado de la expresión génica mediante la identificación de genes o rutas moleculares que se expresan de forma diferencial (DE) entre dos o más condiciones biológicas. Este tutorial muestra un flujo de trabajo computacional para la detección de genes y rutas de expresión diferencial a partir de datos de ARN-Seq, proporcionando un análisis completo de un experimento de ARN-Seq en células de *Drosophila* tras la eliminación de un gen regulador.

En el estudio de {% cite brooks2011conservation %}, los autores identificaron genes y vías reguladas por el gen *Pasilla* (el homólogo en *Drosophila* de las proteínas reguladoras del splicing en mamíferos Nova-1 y Nova-2) usando datos de RNA-Seq. Redujeron el gen *Pasilla* (*PS*) en *Drosophila melanogaster* mediante ARN de interferencia (ARNi). A continuación, se aisló el ARN total y se utilizó para preparar bibliotecas de ARN-Seq de extremo único y de extremo pareado para muestras tratadas (sin PS) y no tratadas. Estas bibliotecas se secuenciaron para obtener las lecturas de RNA-Seq de cada muestra. Los datos de RNA-Seq de las muestras tratadas y no tratadas pueden compararse para identificar los efectos de la depleción del gen *Pasilla* en la expresión génica.

En este tutorial, ilustramos paso a paso el análisis de los datos de expresión génica utilizando 7 de los conjuntos de datos originales:

- 4 muestras no tratadas: [GSM461176](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461176), [GSM461177](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461177), [GSM461178](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461178), [GSM461182](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461182)
- 3 muestras tratadas (gen *Pasilla* agotado por RNAi): [GSM461179](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461179), [GSM461180](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461180), [GSM461181](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461181)

Cada muestra constituye una réplica biológica separada de la condición correspondiente (tratada o no tratada). Además, dos de las muestras tratadas y dos de las no tratadas proceden de un ensayo de secuenciación de extremo pareado, mientras que las muestras restantes proceden de un experimento de secuenciación de extremo único.

> <comment-title>Datos completos</comment-title>
> 
> Los datos originales están disponibles en NCBI Gene Expression Omnibus (GEO) con el número de acceso [GSE18508](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18508). Las lecturas de ARN-Seq en bruto se han extraído de los archivos Sequence Read Archive (SRA) y se han convertido en archivos FASTQ.
> 
{: .comment}

> <agenda-title></agenda-title>
> 
> En este tutorial, nos ocuparemos de:
> 
> 1. TOC
> {:toc}
{: .agenda}

# Carga de datos

En la primera parte de este tutorial utilizaremos los archivos de 2 de las 7 muestras para demostrar cómo calcular el recuento de lecturas (una medida de la expresión génica) a partir de archivos FASTQ (control de calidad, mapeo, recuento de lecturas). Proporcionamos los archivos FASTQ para las otras 5 muestras por si desea reproducir todo el análisis más adelante.

En la segunda parte del tutorial, se utilizan los recuentos de lecturas de las 7 muestras para identificar y visualizar los genes, familias de genes y rutas moleculares DE debido a la depleción del gen *PS*.

> <hands-on-title>Carga de datos</hands-on-title>
> 
> 1. Crear un nuevo historial para este ejercicio de RNA-Seq
> 
>    {% snippet faqs/galaxy-es/histories_create_new.md %}
> 
> 2. Importar los pares de archivos FASTQ desde [Zenodo]({{ page.zenodo_link }}) o una biblioteca de datos:
>    - `GSM461177` (sin tratar): `GSM461177_1` y `GSM461177_2`
>    - `GSM461180` (tratado): `GSM461180_1` y `GSM461180_2`
> 
>    ```text
>    {{ page.zenodo_link }}/files/GSM461177_1.fastqsanger
>    {{ page.zenodo_link }}/files/GSM461177_2.fastqsanger
>    {{ page.zenodo_link }}/files/GSM461180_1.fastqsanger
>    {{ page.zenodo_link }}/files/GSM461180_2.fastqsanger
>    ```
> 
>    {% snippet faqs/galaxy-es/datasets_import_via_link.md %}
> 
>    {% snippet faqs/galaxy-es/datasets_import_from_data_library.md %}
> 
>    > <comment-title></comment-title>
>    > 
>    > Tenga en cuenta que estos son los archivos completos de las muestras y ~1.5Gb cada uno, por lo que puede tardar algunos minutos en importarse.
>    > 
>    > Para una rápida ejecución de los pasos FASTQ un pequeño subconjunto de cada archivo FASTQ (~ 5Mb) se puede encontrar aquí en [Zenodo]({{ page.zenodo_link }}):
>    > 
>    > ```text
>    > {{ page.zenodo_link }}/files/GSM461177_1_subsampled.fastqsanger
>    > {{ page.zenodo_link }}/files/GSM461177_2_subsampled.fastqsanger
>    > {{ page.zenodo_link }}/files/GSM461180_1_subsampled.fastqsanger
>    > {{ page.zenodo_link }}/files/GSM461180_2_subsampled.fastqsanger
>    > ```
>    > 
>    {: .comment}
> 
> 3. Compruebe que el tipo de datos es `fastqsanger` (por ejemplo, **no** `fastq`). Si no lo es, cambie el tipo de datos a `fastqsanger`.
> 
>    {% snippet faqs/galaxy-es/datasets_change_datatype.md datatype="fastqsanger" %}
> 
> 4. Cree una colección emparejada llamada `2 PE fastqs`, nombre sus pares con el nombre de la muestra seguido de los atributos: `GSM461177_untreat_paired` y `GSM461180_treat_paired`.
> 
>    {% snippet faqs/galaxy-es/collections_build_list_paired.md %}
> 
{: .hands_on}

{% include topics/sequence-analysis/tutorials/quality-control/fastq_question_ES.md %}

Las lecturas son datos brutos de la máquina de secuenciación sin ningún tratamiento previo. Es necesario evaluar su calidad.

# Control de calidad

Durante la secuenciación se introducen errores, como la llamada de nucleótidos incorrectos. Esto se debe a las limitaciones técnicas de cada plataforma de secuenciación. Los errores de secuenciación pueden sesgar el análisis y dar lugar a una interpretación errónea de los datos. También puede haber adaptadores si las lecturas son más largas que los fragmentos secuenciados, y recortarlos puede mejorar el número de lecturas mapeadas.

El control de calidad de la secuencia es, por tanto, un primer paso esencial en su análisis. Utilizaremos herramientas similares a las descritas en el tutorial ["Control de calidad"]({% link topics/sequence-analysis/tutorials/quality-control/tutorial.md %}):
- [Falco](https://falco.readthedocs.io/en/latest/), que es una reescritura optimizada para la eficiencia de [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), para crear un informe de la calidad de la secuencia
- [MultiQC](https://multiqc.info/) ({% cite ewels2016multiqc %}) para agregar los informes generados
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) ({% cite marcel2011cutadapt %}) para mejorar la calidad de las secuencias mediante el recorte y el filtrado.

Lamentablemente, la versión actual de MultiQC (la herramienta que utilizamos para combinar informes) no admite colecciones de listas de pares. Primero tendremos que transformar nuestra lista de pares en una lista simple.

> <details-title>¿Qué significa esto exactamente?</details-title>
> 
> La situación actual está en la parte superior y la herramienta **Aplanar colección** la transformará en la situación que se muestra en la parte inferior:
> ![Flatten](../../images/ref-based/flatten.png "Aplanar la lista de pares a lista")
{: .details}

> <hands-on-title>Control de calidad</hands-on-title>
> 
> 1. {% tool [Flatten collection](__FLATTEN__) %} con los siguientes parámetros convertir la lista de pares en una lista simple:
>     - *"Input Collection "*: `2 PE fastqs`
> 
> 2. {% tool [Falco](toolshed.g2.bx.psu.edu/repos/iuc/falco/falco/1.2.4+galaxy0) %} con los siguientes parámetros:
>    - {% icon param-collection %} *"Raw read data from your current history "*: Salida de **Flatten collection** {% icon tool %} seleccionada como **Dataset collection**
> 
>    {% snippet faqs/galaxy-es/tools_select_collection.md %}
> 
> 3. Inspeccionar la página web de salida de **Falco** {% icon tool %} para la muestra `GSM461177_untreat_paired` (forward y reverse)
> 
>    > <question-title></question-title>
>    > 
>    > ¿Cuál es la longitud de la lectura?
>    > 
>    > > <solution-title></solution-title>
>    > > 
>    > > La longitud de lectura de ambas parejas es de 37bp.
>    > > 
>    > {: .solution}
>    {: .question}
> 
>    Como es tedioso inspeccionar todos estos informes individualmente los combinaremos con {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.11+galaxy1) %}.
> 
> 4. {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.24.1+galaxy0) %} para agregar los informes Falco con los siguientes parámetros:
>    - En *"Results "*:
>        - *"Results "*
>            - *"Which tool was used generate logs?*: `FastQC` (**Falco** es un sustituto directo de *FastQC* y podemos pasar su salida a MultiQC como si hubiera sido generada por la herramienta original.)
>                - En *"FastQC output"*:
>                    - {% icon param-repeat %} *"Insert FastQC output"*
>                        - {% icon param-collection %} *"FastQC output"*: `Falco on collection N: RawData` (salida de **Falco** {% icon tool %})
> 
> 5. Inspeccionar la página web de salida de MultiQC para cada FASTQ
> 
>    > <question-title></question-title>
>    > 
>    > 1. ¿Qué opina de la calidad de las secuencias?
>    > 2. ¿Qué debemos hacer?
>    > 
>    > > <solution-title></solution-title>
>    > > 
>    > > 1. Todo parece correcto para 3 de los archivos. El `GSM461177_untreat_paired` tiene 10,6 millones de secuencias emparejadas y el `GSM461180_treat_paired` 12,3 millones de secuencias emparejadas. Pero en `GSM461180_treat_paired_reverse` (reverse reads de GSM461180) la calidad disminuye bastante al final de las secuencias.
>    > > 
>    > >    Todos los archivos excepto `GSM461180_treat_paired_reverse` tienen una alta proporción de lecturas duplicadas (esperable en datos RNA-Seq).
>    > > 
>    > >    ![Recuentos de secuencias](../../images/ref-based/fastqc_sequence_counts_plot.png "Recuentos de secuencias")
>    > > 
>    > >    La "calidad de la secuencia por base" es globalmente buena, con un ligero descenso al final de las secuencias. Para `GSM461180_treat_paired_reverse`, la disminución es bastante grande.
>    > > 
>    > >    ![Calidad de la secuencia](../../images/ref-based/fastqc_per_base_sequence_quality_plot.png "Calidad de la secuencia")
>    > > 
>    > >    La puntuación de calidad media de las lecturas es bastante alta, pero la distribución es ligeramente diferente para `GSM461180_treat_paired_reverse`.
>    > > 
>    > >    ![Puntuaciones de calidad por secuencia](../../images/ref-based/fastqc_per_sequence_quality_scores_plot.png "Puntuaciones de calidad por secuencia")
>    > > 
>    > >    Las lecturas no siguen realmente una distribución normal del contenido de GC, excepto para `GSM461180_treat_paired_reverse`.
>    > > 
>    > >    ![Contenido GC por secuencia](../../images/ref-based/fastqc_per_sequence_gc_content_plot.png "Contenido GC por secuencia")
>    > > 
>    > >    La proporción de N en las lecturas (bases que no pudieron ser llamadas) es baja.
>    > > 
>    > >    ![Contenido N por base](../../images/ref-based/fastqc_per_base_n_content_plot.png "Contenido N por base")
>    > > 
>    > >    Secuencias duplicadas: >10 a >500
>    > > 
>    > >    ![Niveles de duplicación de secuencias](../../images/ref-based/fastqc_sequence_duplication_levels_plot.png "Niveles de duplicación de secuencias")
>    > > 
>    > >    Casi no hay adaptadores conocidos ni secuencias sobrerrepresentadas.
>    > > 
>    > > 2. Si la calidad de las lecturas es mala, deberíamos:
>    > >    1. Compruebe qué está mal y piense en las posibles razones de la mala calidad de las lecturas: puede deberse al tipo de secuenciación o a lo que hemos secuenciado (gran cantidad de secuencias sobrerrepresentadas en datos transcriptómicos, porcentaje sesgado de bases en datos Hi-C)
>    > >    2. Pregunte al centro de secuenciación al respecto
>    > >    3. Realizar algún tratamiento de calidad (teniendo cuidado de no perder demasiada información) con algún recorte o eliminación de malas lecturas
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

Debemos recortar las lecturas para eliminar las bases que se secuenciaron con alta incertidumbre (es decir, bases de baja calidad) en los extremos de las lecturas, y también eliminar las lecturas de mala calidad general.

{% include topics/sequence-analysis/tutorials/quality-control/paired_end_question_ES.md forward="GSM461177_untreat_paired_forward" reverse="GSM461177_untreat_paired_reverse" %}

> <hands-on-title>Recorte de FASTQs</hands-on-title>
> 
> 1. {% tool [Cutadapt](toolshed.g2.bx.psu.edu/repos/lparsons/cutadapt/cutadapt/4.9+galaxy1) %} con los siguientes parámetros para recortar las secuencias de baja calidad:
>    - *"Single-end or Paired-end reads?*: `Paired-end Collection`
>       - {% icon param-collection %} *"Paired Collection "*: `2 PE fastqs`
>    - En *"Other Read Trimming Options"*
>       - *"Quality cutoff(s) (R1)"*: `20`
>    - En *"Read Filtering Options"*
>       - *"Minimum length (R1)"*: `20`
>    - En *"Additional outputs to generate"*
>       - Seleccionar: `Report: Cutadapt's per-adapter statistics. You can use this file with MultiQC.`
> 
>      {% snippet topics/sequence-analysis/tutorials/quality-control/trimming_question_ES.md %}
> 
> 2. {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.11+galaxy1) %} para agregar los informes Cutadapt con los siguientes parámetros:
>    - En *"Results"*:
>        - {% icon param-repeat %} *"Results"*
>            - *"Which tool was used generate logs?"*: `Cutadapt/Trim Galore!`
>               - {% icon param-collection %} *"Output of Cutadapt"*: `Cutadapt on collection N: Report` (salida de **Cutadapt** {% icon tool %}) seleccionada como **Dataset collection**
>
>    > <question-title></question-title>
>    > 
>    > 1. ¿Cuántos pares de secuencias se han eliminado porque al menos una lectura era más corta que el corte de longitud?
>    > 2. ¿Cuántos pares de bases se han eliminado de las lecturas directas debido a la mala calidad? ¿Y de las lecturas inversas?
>    > 
>    > > <solution-title></solution-title>
>    > > 
>    > > 1. 147.810 (1,4%) lecturas eran demasiado cortas para `GSM461177_untreat_paired` y 1.101.875 (9%) para `GSM461180_treat_paired`.
>    > >
>    > >    ![lecturas filtradas por Cutadapt](../../images/ref-based/cutadapt_filtered_reads_plot.png "Lecturas filtradas por Cutadapt")
>    > >
>    > > 2. La salida MultiQC sólo proporciona la proporción de pb recortados en total, no para cada lectura. Para obtener esta información, es necesario volver a los informes individuales. Para `GSM461177_untreat_paired`, se han recortado 5.072.810 pb de las lecturas directas (Lectura 1) y 8.648.619 pb de las lecturas inversas (Lectura 2) debido a la calidad. En el caso de `GSM461180_treat_paired`, se han recortado 10.224.537 pb de las lecturas anteriores y 51.746.850 pb de las inversas. Esto no es una sorpresa; vimos que al final de las lecturas la calidad disminuía más para las lecturas inversas que para las lecturas delanteras, especialmente para `GSM461180_treat_paired`.
>    > >
>    > {: .solution }
>    {: .question}
{: .hands_on}

# Mapeo

Para dar sentido a las lecturas, primero tenemos que averiguar de dónde proceden las secuencias en el genoma, para poder determinar a qué genes pertenecen. Cuando se dispone de un genoma de referencia del organismo, este proceso se conoce como alineación o "mapeo" de las lecturas con la referencia. Esto equivale a resolver un rompecabezas, pero, por desgracia, no todas las piezas son únicas.

> <comment-title></comment-title>
> 
> ¿Desea obtener más información sobre los principios en los que se basa el mapeo? Siga nuestra [formación]({% link topics/sequence-analysis/tutorials/mapping/tutorial.md %}).
> 
{: .comment}

En este estudio, los autores utilizaron células de *Drosophila melanogaster*. Por lo tanto, debemos mapear las secuencias de calidad controlada al genoma de referencia de *Drosophila melanogaster*.

{% include topics/sequence-analysis/tutorials/mapping/ref_genome_explanation_ES.md answer_3="El genoma de *Drosophila melanogaster* es conocido y ensamblado y puede utilizarse como genoma de referencia en este análisis. Tenga en cuenta que pueden publicarse nuevas versiones de genomas de referencia si el ensamblaje mejora, para este tutorial vamos a utilizar la versión 6 del ensamblaje del genoma de referencia de *Drosophila melanogaster* [(dm6)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4383921/)."%}

En los transcriptomas de eucariotas, la mayoría de las lecturas proceden de ARNm procesados que carecen de intrones:

![Types of RNA-Seq reads](../../images/ref-based/rna-seq-reads.png "Los tipos de lecturas RNA-seq (adaptación de la Figura 1a de {% cite kim2015hisat %}): lecturas que mapean completamente dentro de un exón (en rojo), lecturas que abarcan más de 2 exones (en azul), lecturas que abarcan más de 2 exones (en púrpura)")

Por lo tanto no pueden ser simplemente mapeadas de vuelta al genoma como hacemos normalmente con los datos de ADN. Se han desarrollado mapeadores empalmados para mapear eficientemente las lecturas derivadas de transcripciones contra un genoma de referencia:

![Splice-aware alignment](../../images/transcriptomics_images/splice_aware_alignment.png "Principio de los mapeadores empalmados: (1) identificación de las lecturas que abarcan un único exón, (2) identificación de las uniones de empalme en las lecturas no mapeadas")

> <details-title>Más detalles sobre los distintos mapeadores empalmados</details-title>
> 
> En los últimos años se han desarrollado varios mapeadores empalmados para procesar la explosión de datos de RNA-Seq.
> 
> [**TopHat**](https://ccb.jhu.edu/software/tophat/index.shtml) ({% cite trapnell2009tophat %}) fue una de las primeras herramientas diseñadas específicamente para abordar este problema. En **TopHat** las lecturas se mapean contra el genoma y se separan en dos categorías: (1) las que se mapean, y (2) las que están inicialmente sin mapear (IUM). las "pilas" de lecturas que representan posibles exones se extienden en busca de posibles sitios de empalme donante/aceptor y se reconstruyen las posibles uniones de empalme. A continuación, se mapean los IUM en estas uniones.
> 
> ![TopHat](../../images/transcriptomics_images/tophat.png "TopHat (Figura 1 de {% cite trapnell2009tophat %})")
> 
> **TopHat** se ha mejorado posteriormente con el desarrollo de **TopHat2** ({% cite kim2013tophat2 %}):
> 
> ![TopHat2](../../images/transcriptomics_images/13059_2012_Article_3053_Fig6_HTML.jpg "TopHat2 (Figura 6 de {% cite kim2013tophat2 %})")
> 
> Para optimizar y acelerar aún más la alineación de lecturas empalmadas, se desarrolló [**HISAT2**](https://ccb.jhu.edu/software/hisat2/index.shtml) ({% cite kim2019graph %}). Utiliza un índice de grafo jerárquico [FM](https://en.wikipedia.org/wiki/FM-index) (HGFM), que representa el genoma completo y las posibles variantes, junto con índices locales superpuestos (cada uno de los cuales abarca ~57 kb) que cubren colectivamente el genoma y sus variantes. Esto permite encontrar ubicaciones semilla iniciales para potenciales alineaciones de lecturas en el genoma utilizando el índice global y refinar rápidamente estas alineaciones utilizando un índice local correspondiente:
> 
> ![Índice FM del gráfico jerárquico en HISAT/HISAT2](../../images/transcriptomics_images/hisat.png "Índice FM del gráfico jerárquico en HISAT/HISAT2 (Figura S8 de {% cite kim2015hisat %})")
> 
> Una parte de la lectura (flecha azul) se mapea primero en el genoma utilizando el índice FM global. **HISAT2** intenta entonces extender el alineamiento directamente utilizando la secuencia del genoma (flecha violeta). En (**a**) tiene éxito y esta lectura se alinea ya que reside completamente dentro de un exón. En (**b**) la extensión encuentra un desajuste. Ahora **HISAT2** se aprovecha del índice FM local que se solapa con esta localización para encontrar el mapeo apropiado para el resto de esta lectura (flecha verde). La (**c**) muestra una combinación de estas dos estrategias: el comienzo de la lectura se mapea usando el índice FM global (flecha azul), se extiende hasta que alcanza el final del exón (flecha violeta), se mapea usando el índice FM local (flecha verde) y se extiende de nuevo (flecha violeta).
> 
> [**STAR** aligner](https://github.com/alexdobin/STAR) ({% cite dobin2013star %}) es una alternativa rápida para mapear lecturas RNA-Seq contra un genoma de referencia utilizando un [suffix array](https://en.wikipedia.org/wiki/Suffix_array) sin comprimir. Funciona en dos etapas. En la primera etapa realiza una búsqueda de semillas:
> 
> ![STAR's seed search](../../images/transcriptomics_images/star.png "STAR's seed search (Figure 1 from {% cite dobin2013star %})")
> 
> Aquí una lectura se divide entre dos exones consecutivos. **STAR** empieza a buscar un prefijo mapeable máximo (MMP) desde el principio de la lectura hasta que ya no puede coincidir continuamente. Después de este punto, empieza a buscar un MMP para la parte no coincidente de la lectura (**a**). En el caso de emparejamientos erróneos (**b**) y regiones no alineables (**c**), los MMP sirven como anclas a partir de las cuales extender los alineamientos.
> 
> En la segunda etapa **STAR** stitches MMPs to generate read-level alignments that (contrary to MMPs) can contain mismatches and indels. Se utiliza un esquema de puntuación para evaluar y priorizar las combinaciones de cosido y para evaluar las lecturas que mapean múltiples localizaciones. **STAR** es extremadamente rápido, pero requiere una cantidad considerable de RAM para funcionar con eficacia.
> 
{: .details}

## Mapeo

Mapearemos nuestras lecturas en el genoma de *Drosophila melanogaster* usando **STAR** ({% cite dobin2013star %}).

> <hands-on-title>Mapeo empalmado</hands-on-title>
> 
> 1. Importe la anotación de genes Ensembl para *Drosophila melanogaster* (`Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz`) de la biblioteca de Datos Compartidos si está disponible o de [Zenodo]({{ page.zenodo_link }}/files/Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz) a su actual historial Galaxy
> 
>    ```text
>    {{ page.zenodo_link }}/files/Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz
>    ```
> 
>    1. Cambie el nombre del conjunto de datos si es necesario
>    2. Compruebe que el tipo de datos es `gtf` y no `gff`, y que la base de datos es `dm6`
> 
>    > <comment-title>¿Cómo obtener el archivo de anotación?</comment-title>
>    > 
>    > Los archivos de anotación de los organismos modelo pueden estar disponibles en la biblioteca de Datos Compartidos (la ruta a ellos cambiará de un servidor Galaxy a otro). También puede recuperar el archivo de anotación de UCSC (utilizando la herramienta **UCSC Main**).
>    > 
>    > Para generar este archivo específico, el archivo de anotación se descargó de Ensembl, que proporciona una base de datos de transcritos más completa, y se adaptó para que funcionara con el genoma dm6, instalado en servidores Galaxy compatibles.
>    > 
>    {: .comment}
> 
> 2. {% tool [RNA STAR](toolshed.g2.bx.psu.edu/repos/iuc/rgrnastar/rna_star/2.7.11a+galaxy0) %} con los siguientes parámetros para mapear sus lecturas en el genoma de referencia:
>    - *"Single-end or paired-end reads "*: `Paired-end (as collection)`
>       - {% icon param-collection %} *"RNA-Seq FASTQ/FASTA paired reads "*: the `Cutadapt on collection N: Reads` (output of **Cutadapt** {% icon tool %})
>    - *"Custom or built-in reference genome" "*: `Use a built-in index`
>       - *"Reference genome with or without an annotation "*: `use genome reference without builtin gene-model but provide a gtf`
>           - *"Select genome reference"*: `Fly (Drosophila melanogaster): dm6 Full`
>           - {% icon param-file %} *"Gene model (gff3,gtf) file for splice junctions "*: el `Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz` importado
>           - *"Length of the genomic sequence around annotated junctions "*: `36`
> 
>             Este parámetro debe ser la longitud de las lecturas - 1
>                - *"Per gene/transcript output "*: `Per gene read counts (GeneCounts)`
>    - *"Compute coverage "*:
>       - `Yes in bedgraph format`
> 
> 3. {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.11+galaxy1) %} para agregar los logs STAR con los siguientes parámetros:
>    - En *"Results "*:
>        - *"Results "*
>            - *"Which tool was used generate logs?"*: `STAR`
>                - En *"STAR output "*:
>                    - {% icon param-repeat %} *"Insertar STAR output"*
>                        - *"Type of STAR output? "*: `Log`
>                            - {% icon param-collection %} *"STAR log output "*: `RNA STAR on collection N: log` (salida de **RNA STAR** {% icon tool %})
> 
>    > <question-title></question-title>
>    > 
>    > 1. ¿Qué porcentaje de lecturas se mapean exactamente una vez para ambas muestras?
>    > 2. ¿Cuáles son las otras estadísticas disponibles?
>    > 
>    > > <solution-title></solution-title>
>    > > 
>    > > 1. Más del 83% para `GSM461177_untreat_paired` y del 79% para `GSM461180_treat_paired`. Podemos continuar con el análisis, ya que sólo los porcentajes inferiores al 70% deben investigarse para detectar una posible contaminación.
>    > > 2. También tenemos acceso al número y porcentaje de lecturas que están mapeadas en varias localizaciones, mapeadas en demasiadas localizaciones diferentes, no mapeadas porque son demasiado cortas.
>    > > 
>    > >    ![Puntuaciones de alineación STAR](../../images/ref-based/star_alignment_plot.png "Puntuaciones de alineación")
>    > > 
>    > >    Podríamos haber sido probablemente más estrictos en la longitud mínima de lectura para evitar estas lecturas no mapeadas debido a su longitud.
>    > {: .solution}
>    {: .question}
> 
{: .hands_on}

Según el informe **MultiQC**, alrededor del 80% de las lecturas de ambas muestras se mapean exactamente una vez con el genoma de referencia. Podemos proseguir con el análisis, ya que sólo los porcentajes inferiores al 70% deben investigarse para detectar una posible contaminación. Ambas muestras tienen un porcentaje bajo (menos del 10%) de lecturas que se corresponden con múltiples ubicaciones en el genoma de referencia. Esto está dentro del rango normal para la secuenciación de lectura corta de Illumina, pero puede ser menor para los nuevos conjuntos de datos de secuenciación de lectura larga que pueden abarcar regiones repetidas más grandes en el genoma de referencia y será mayor para las bibliotecas de extremo 3'.

La salida principal de **STAR** es un archivo BAM.

{% include topics/sequence-analysis/tutorials/mapping/bam_explanation_ES.md mapper="RNA STAR" %}

## Inspección de los resultados del mapeo

El archivo BAM contiene información de todas nuestras lecturas, lo que dificulta su inspección y exploración en formato de texto. Una potente herramienta para visualizar el contenido de los archivos BAM es Integrative Genomics Viewer (**IGV**, {% cite robinson2011integrative %}).

> <hands-on-title>Inspección de los resultados del mapeo</hands-on-title>
> 
> 1. Instale [**IGV**](https://software.broadinstitute.org/software/igv/download) (si no está ya instalado)
> 2. Iniciar IGV localmente
> 3. Haga clic en la colección `RNA STAR on collection N: mapped.bam` (salida de **RNA STAR** {% icon tool %})
> 4. Expandir el {% icon param-file %} archivo `GSM461177_untreat_paired`.
> 5. Haga clic en el icono {% icon galaxy-barchart %} visualizar en el bloque de archivo `GSM461177`.
> 6. En el panel central haga clic en `local` en `display with IGV (local, D. melanogaster (dm6))` para cargar las lecturas en el navegador IGV
>    > <comment-title></comment-title>
>    > 
>    > Para que este paso funcione, necesitará tener IGV o [Java Web Start](https://www.java.com/en/download/faq/java_webstart.xml) instalado en su máquina. Sin embargo, las preguntas de esta sección también pueden responderse inspeccionando las capturas de pantalla de IGV que aparecen a continuación.
>    > 
>    > Consulte la [documentación de IGV](https://software.broadinstitute.org/software/igv/AlignmentData) para obtener más información.
>    > 
>    {: .comment}
> 
> 6. **IGV** {% icon tool %}: Zoom a `chr4:540,000-560,000` (Cromosoma 4 entre 540 kb y 560 kb)
> 
>    > <question-title></question-title>
>    > 
>    > ![Captura de pantalla de la vista IGV en el cromosoma 4](../../images/transcriptomics_images/junction_igv_screenshot.png "Captura de pantalla de IGV en el cromosoma 4")
>    > 
>    > 1. ¿Qué información aparece en la parte superior como picos grises?
>    > 2. ¿Qué indican las líneas de conexión entre algunas de las lecturas alineadas?
>    > 
>    > > <solution-title></solution-title>
>    > > 
>    > > 1. El gráfico de cobertura: la suma de lecturas mapeadas en cada posición
>    > > 2. Indican eventos de unión (o sitios de empalme), *es decir* lecturas que se mapean a través de un intrón
>    > {: .solution}
>    {: .question}
> 
> 7. **IGV** {% icon tool %}: Inspeccione las uniones de empalme utilizando un **trazado Sashimi**
> 
>    > <comment-title>Creación de un gráfico Sashimi</comment-title>
>    > 
>    > - Haga clic con el botón derecho en el archivo BAM (en IGV)
>    > - Seleccione **Sashimi Plot** del menú
>    {: .comment}
> 
> > 
> > <question-title></question-title>
> > 
> > ![Captura de pantalla de un gráfico Sashimi del cromosoma 4](../../images/transcriptomics_images/star_igv_sashimi.png "Captura de pantalla de un gráfico Sashimi del cromosoma 4")
> > 
> > 1. ¿Qué representa el gráfico de barras rojas verticales? ¿Y los arcos con números?
> > 2. ¿Qué significan los números de los arcos?
> > 3. ¿Por qué observamos diferentes grupos apilados de cajas azules enlazadas en la parte inferior?
> > 
> > > <solution-title></solution-title>
> > > 
> > > 1. La cobertura de cada pista de alineamiento se representa en un gráfico de barras rojas. Los arcos representan las uniones de empalme observadas, es decir, las lecturas que abarcan intrones.
> > > 2. Los números se refieren al número de lecturas de unión observadas.
> > > 3. Los diferentes grupos de cuadros enlazados en la parte inferior representan los diferentes transcritos de los genes en esta ubicación que están presentes en el archivo GTF.
> > > 
> > {: .solution}
> {: .question}
> 
> > <comment-title></comment-title>
> > 
> > Consulte la [documentación de IGV sobre gráficos Sashimi](https://software.broadinstitute.org/software/igv/Sashimi) para encontrar algunas pistas
> > 
> {: .comment}
{: .hands_on}

> <details-title>Comprobación adicional de la calidad de los datos</details-title>
> 
> La calidad de los datos y del mapeo puede comprobarse aún más, por ejemplo, inspeccionando el nivel de duplicación de lecturas, el número de lecturas mapeadas a cada cromosoma, la cobertura del cuerpo génico y la distribución de lecturas a través de las características.
> 
> *Estos pasos se han inspirado en los proporcionados en el [gran tutorial "RNA-Seq reads to counts"]({% link topics/transcriptomics/tutorials/rna-seq-reads-to-counts/tutorial.md %}) y se han adaptado a nuestros conjuntos de datos.*
> 
> #### Lecturas duplicadas
> 
> En el informe Falco/MultiQC, hemos visto que algunas lecturas están duplicadas:
> 
> ![Niveles de duplicación de secuencias](../../images/ref-based/fastqc_sequence_duplication_levels_plot.png "Niveles de duplicación de secuencias")
> 
> Las lecturas duplicadas pueden provenir de genes altamente expresados, por lo que normalmente se mantienen en el análisis de expresión diferencial RNA-Seq. Pero un alto porcentaje de duplicados puede indicar un problema, por ejemplo, una sobreamplificación durante la PCR de una biblioteca de baja complejidad.
> 
> **MarkDuplicates** de [Picard suite](http://broadinstitute.github.io/picard/) examina registros alineados de un archivo BAM para localizar lecturas duplicadas, es decir, lecturas que se mapean en la misma ubicación (basándose en la posición de inicio del mapeo).
> 
> > <hands-on-title>Comprobar lecturas duplicadas</hands-on-title>
> > 
> > 1. {% tool [MarkDuplicates](toolshed.g2.bx.psu.edu/repos/devteam/picard/picard_MarkDuplicates/2.18.2.4) %} con los siguientes parámetros:
> >    - {% icon param-collection %} *"Seleccionar conjunto de datos SAM/BAM o colección de conjuntos de datos "*: `RNA STAR on collection N: mapped.bam` (salida de **RNA STAR** {% icon tool %})
> > 
> > 2. {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.11+galaxy1) %} para agregar los logs de MarkDuplicates con los siguientes parámetros:
> >    - En *"Resultados "*:
> >        - *"Resultados "*
> >            - *"Which tool was used generate logs?"*: `Picard`
> >                - En *"salida Picard "*:
> >                    - {% icon param-repeat %} *"Insertar salida Picard "*
> >                        - *"¿Tipo de salida Picard? "*: `Markdups`
> >                        - {% icon param-collection %} *"Picard output "*: `MarkDuplicates on collection N: MarkDuplicate metrics` (salida de **MarkDuplicates** {% icon tool %})
> > 
> >    > <question-title></question-title>
> >    > 
> >    > ¿Cuáles son los porcentajes de lecturas duplicadas para cada muestra?
> >    > 
> >    > > <solution-title></solution-title>
> >    > > 
> >    > > La muestra `GSM461177_untreat_paired` tiene un 25,9% de lecturas duplicadas mientras que `GSM461180_treat_paired` tiene un 27,8%.
> >    > {: .solution}
> >    {: .question}
> {: .hands_on}
> 
> En general, obtener hasta un 50% de lecturas duplicadas se considera normal. Por lo tanto, nuestras dos muestras están bien.
> 
> #### Número de lecturas asignadas a cada cromosoma
> 
> Para evaluar la calidad de la muestra (por ejemplo, exceso de contaminación mitocondrial), podemos comprobar el sexo de las muestras, o para ver si algún cromosoma tiene genes altamente expresados, podemos comprobar el número de lecturas mapeadas a cada cromosoma utilizando **IdxStats** de la suite **Samtools**.
> 
> > <hands-on-title>Comprobar el número de lecturas asignadas a cada cromosoma</hands-on-title>
> > 
> > 1. {% tool [Samtools idxstats](toolshed.g2.bx.psu.edu/repos/devteam/samtools_idxstats/samtools_idxstats/2.0.4) %} con los siguientes parámetros:
> >    - {% icon param-collection %} *"Archivo BAM "*: `RNA STAR on collection N: mapped.bam` (salida de **RNA STAR** {% icon tool %})
> > 
> > 2. {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.11+galaxy1) %} para agregar los logs de idxstats con los siguientes parámetros:
> >    - En *"Resultados "*:
> >        - *"Resultados "*
> >            - *"Which tool was used generate logs?"*: `Samtools`
> >                - En *"Samtools output "*:
> >                    - {% icon param-repeat %} *"Insertar salida Samtools "*
> >                        - *"¿Tipo de salida Samtools? "*: `idxstats`
> >                            - {% icon param-collection %} *"Samtools idxstats output "*: `Samtools idxstats on collection N` (salida de **Samtools idxstats** {% icon tool %})
> > 
> >    > <question-title></question-title>
> >    > 
> >    > ![Samtools idxstats](../../images/ref-based/samtools-idxstats-mapped-reads-plot.png)
> >    > 
> >    > 1. ¿Cuántos cromosomas tiene el genoma de *Drosophila*?
> >    > 2. ¿Dónde se encuentran la mayoría de las lecturas?
> >    > 3. ¿Podemos determinar el sexo de las muestras?
> >    > 
> >    > > <solution-title></solution-title>
> >    > > 
> >    > > 1. El genoma de *Drosophila* tiene 4 pares de cromosomas: X/Y, 2, 3 y 4.
> >    > > 2. Las lecturas corresponden principalmente al cromosoma 2 (chr2L y chr2R), 3 (chr3L y chr3R) y X. Sólo unas pocas lecturas corresponden al cromosoma 4, lo que es de esperar dado que este cromosoma es muy pequeño.
> >    > > 3. A juzgar por el porcentaje de lecturas X+Y, la mayoría de las lecturas se asignan a X y sólo unas pocas a Y. Esto indica que probablemente no hay muchos genes en Y, por lo que las muestras son probablemente femeninas.
> >    > > 
> >    > >    ![Samtools idxstats](../../images/ref-based/samtools-idxstats-xy-plot.png)
> >    > {: .solution}
> >    {: .question}
> {: .hands_on}
> 
> #### Cobertura del cuerpo genético
> 
> Las diferentes regiones de un gen forman el cuerpo del gen. Es importante comprobar si la cobertura de lectura es uniforme en todo el cuerpo del gen. Por ejemplo, un sesgo hacia el extremo 5' de los genes podría indicar la degradación del ARN. Alternativamente, un sesgo hacia el extremo 3' podría indicar que los datos proceden de un ensayo 3'. Para evaluar esto, podemos utilizar la herramienta **Gene Body Coverage** del conjunto de herramientas RSeQC ({% cite wang2012rseqc %}). Esta herramienta escala todas las transcripciones a 100 nucleótidos (utilizando un archivo de anotación proporcionado) y calcula el número de lecturas que cubren cada posición de nucleótido (escalada). Como esta herramienta es muy lenta, calcularemos la cobertura sólo con 200.000 lecturas aleatorias.
> 
> > <hands-on-title>Comprobar la cobertura del cuerpo del gen</hands-on-title>
> > 
> > 1. {% tool [Samtools view](toolshed.g2.bx.psu.edu/repos/iuc/samtools_view/samtools_view/1.15.1+galaxy0) %} con los siguientes parámetros:
> >    - {% icon param-collection %} *"Conjunto de datos SAM/BAM/CRAM "*: `mapped_reads` (salida de **RNA STAR** {% icon tool %})
> >    - *"¿Qué te gustaría mirar?*: `A filtered/subsampled selection of reads`
> >        - En *"Configurar submuestreo "*:
> >            - *"Subsample alignment "*: `Specify a target # of reads`
> >                - *"Objetivo # de lecturas "*: `200000`
> >                - *"Semilla para el generador de números aleatorios "*: `1`
> >        - *"¿Qué le gustaría que se informara? "*: `All reads retained after filtering and subsampling`
> >            - *"Output format "*: `BAM (-b)`
> >    - *"Usar una secuencia de referencia "*: `No`
> > 
> > 2. {% tool [Convert GTF to BED12](toolshed.g2.bx.psu.edu/repos/iuc/gtftobed12/gtftobed12/357) %} para convertir el archivo GTF a BED:
> >    - {% icon param-file %} *"GTF File to convert "*: `Drosophila_melanogaster.BDGP6.32.109.gtf.gz`
> > 
> > 3. {% tool [Gene Body Coverage (BAM)](toolshed.g2.bx.psu.edu/repos/nilesh/rseqc/rseqc_geneBody_coverage/5.0.1+galaxy2) %} con los siguientes parámetros:
> >    - *"Ejecutar cada muestra por separado, o combinar varias muestras en un gráfico "*: `Run each sample separately`
> >        - {% icon param-collection %} *"Input .bam file "*: salida de la vista **Samtools** {% icon tool %}
> >    - {% icon param-file %} *"Reference gene model ": `Convert GTF to BED12 on data N: BED12` (salida de **Convert GTF to BED12** {% icon tool %})
> > 
> > 4. {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.11+galaxy1) %} para agregar los resultados RSeQC con los siguientes parámetros:
> >    - En *"Resultados "*:
> >        - *"Resultados "*
> >            - *"Which tool was used generate logs?"*: `RSeQC`
> >                - En *"RSeQC output "*:
> >                    - {% icon param-repeat %} *"Insertar salida RSeQC "*
> >                        - *"¿Tipo de salida RSeQC? "*: `gene_body_coverage`
> >                            - {% icon param-collection %} *"RSeQC gene_body_coverage output "*: `Gene Body Coverage (BAM) on collection N (text)` (salida de **Gene Body Coverage (BAM)** {% icon tool %})
> > 
> >    > <question-title></question-title>
> >    > 
> >    > ![Gene body coverage](../../images/ref-based/rseqc_gene_body_coverage_plot.png)
> >    > 
> >    > ¿Cómo es la cobertura a través de los cuerpos génicos? ¿Las muestras están sesgadas en 3' o 5'?
> >    > 
> >    > > <solution-title></solution-title>
> >    > > 
> >    > > En ambas muestras hay una cobertura bastante uniforme de los extremos 5' a 3' (a pesar de algo de ruido en el medio). Así que no hay sesgo obvio en ambas muestras.
> >    > {: .solution}
> >    {: .question}
> {: .hands_on}
> 
> #### Distribución de las lecturas por características
> 
> Con los datos de RNA-Seq, esperamos que la mayoría de las lecturas correspondan a exones en lugar de a intrones o regiones intergénicas. Antes de seguir adelante con el recuento y el análisis de expresión diferencial, puede ser interesante comprobar la distribución de las lecturas en las características conocidas de los genes (exones, CDS, 5' UTR, 3' UTR, intrones, regiones intergénicas). Por ejemplo, un número elevado de lecturas en regiones intergénicas puede indicar la presencia de contaminación por ADN.
> 
> Aquí utilizaremos la herramienta **Read Distribution** del conjunto de herramientas RSeQC ({% cite wang2012rseqc %}), que utiliza el archivo de anotación para identificar la posición de las diferentes características de los genes.
> 
> > <hands-on-title>Comprobar el número de lecturas asignadas a cada cromosoma</hands-on-title>
> > 
> > 1. {% tool [Read Distribution](toolshed.g2.bx.psu.edu/repos/nilesh/rseqc/rseqc_read_distribution/5.0.1+galaxy2) %} con los siguientes parámetros:
> >    - {% icon param-collection %} *"Archivo .bam/.sam de entrada "*: `RNA STAR on collection N: mapped.bam` (salida de **RNA STAR** {% icon tool %})
> >    - {% icon param-file %} *"Reference gene model ": Archivo BED12 (salida de **Convert GTF to BED12** {% icon tool %})
> > 
> > 2. {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.11+galaxy1) %} para agregar los resultados de la Distribución de Lecturas con los siguientes parámetros:
> >    - En *"Resultados "*:
> >        - *"Resultados "*
> >            - *"Which tool was used generate logs?"*: `RSeQC`
> >                - En *"RSeQC output "*:
> >                    - {% icon param-repeat %} *"Insertar salida RSeQC "*
> >                        - *"¿Tipo de salida RSeQC? "*: `read_distribution`
> >                            - {% icon param-collection %} *"RSeQC read_distribution output "*: `Read Distribution on collection N` (salida de **Distribución de lecturas** {% icon tool %})
> > 
> >    > <question-title></question-title>
> >    > 
> >    > ![Read Distribution](../../images/ref-based/rseqc_read_distribution_plot.png)
> >    > 
> >    > ¿Qué opina de la distribución de las lecturas?
> >    > 
> >    > > <solution-title></solution-title>
> >    > > 
> >    > > La mayoría de las lecturas están mapeadas a exones (>80%), sólo ~2% a intrones y ~5% a regiones intergénicas, que es lo que esperamos. Esto confirma que nuestros datos son datos RNA-Seq y que el mapeo fue exitoso.
> >    > {: .solution}
> >    {: .question}
> {: .hands_on}
> 
> Ahora que hemos comprobado los resultados del mapeo de lecturas, podemos pasar a la siguiente fase del análisis.
> 
{: .details}

Después del mapeo, ahora tenemos la información de dónde están localizadas las lecturas en el genoma de referencia y cómo de bien fueron mapeadas. El siguiente paso en el análisis de datos RNA-Seq es la cuantificación del número de lecturas mapeadas a características genómicas (genes, transcritos, exones, ...).

> <comment-title></comment-title>
> 
> La cuantificación depende tanto del genoma de referencia (el archivo FASTA) como de sus anotaciones asociadas (el archivo GTF). Es extremadamente importante utilizar un archivo de anotaciones que corresponda a la misma versión del genoma de referencia que utilizó para el mapeo (por ejemplo, `dm6` aquí), ya que las coordenadas cromosómicas de los genes suelen ser diferentes entre las distintas versiones del genoma de referencia.
{: .comment}

Aquí nos centraremos en los genes, ya que nos gustaría identificar los que se expresan diferencialmente debido al knockdown del gen Pasilla.

# Recuento del número de lecturas por gen anotado

Para comparar la expresión de genes individuales entre diferentes condiciones (*e.g.* con o sin depleción de PS), un primer paso esencial es cuantificar el número de lecturas por gen, o más específicamente el número de lecturas mapeadas a los exones de cada gen.

![Recuento del número de lecturas por gen anotado](../../images/transcriptomics_images/gene_counting.png "Recuento del número de lecturas por gen anotado")

> <question-title></question-title> 
> 
> En la imagen anterior,
> 
> 1. ¿Cuántas lecturas se encuentran para los diferentes exones?
> 2. ¿Cuántas lecturas se encuentran para los diferentes genes?
> 
> > <solution-title></solution-title>
> > 
> > 1. Número de lecturas por exón
> > 
> >    | Exon          | Number of reads |
> >    | ------------- | --------------- |
> >    | gene1 - exon1 | 3               |
> >    | gene1 - exon2 | 2               |
> >    | gene2 - exon1 | 3               |
> >    | gene2 - exon2 | 4               |
> >    | gene2 - exon3 | 3               |
> > 
> > 2. el gen1 tiene 4 lecturas, no 5, debido al empalme de la última lectura (gen1 - exón1 + gen1 - exón2). el gen2 tiene 6 lecturas, 3 de las cuales están empalmadas.
> > 
> {: .solution}
{: .question}

Existen dos herramientas principales para el recuento de lecturas: [**HTSeq-count**](http://htseq.readthedocs.io/en/release_0.9.1/count.html) ({% cite anders2015htseq %}) o **featureCounts** ({% cite liao2013featurecounts %}). Además, **STAR** permite contar lecturas mientras se mapea: sus resultados son idénticos a los de **HTSeq-count**. Mientras que esta salida es suficiente para la mayoría de los análisis, **featureCounts** ofrece más personalización sobre cómo contar lecturas (calidad mínima de mapeo, contar lecturas en lugar de fragmentos, contar transcripciones en lugar de genes, etc.).

En principio, el recuento de lecturas que se solapan con características genómicas es una tarea bastante sencilla. Pero es necesario determinar la estridencia de la biblioteca. De hecho, este es un parámetro de **featureCounts**. Por el contrario, **STAR** evalúa los recuentos en los tres posibles strandnesses pero usted todavía necesita esta información para extraer los recuentos que corresponden a su biblioteca.

## Estimación de la strandness

Los ARN a los que normalmente se dirigen los experimentos de ARN-Seq son monocatenarios (por ejemplo, ARNm) y, por lo tanto, tienen polaridad (extremos 5' y 3' funcionalmente distintos). Durante un experimento típico de RNA-Seq, la información sobre la cadena se pierde después de que ambas cadenas de cDNA se sinteticen, se seleccionen por tamaño y se conviertan en una biblioteca de secuenciación. Sin embargo, esta información puede ser muy útil para el paso de recuento de lecturas, especialmente para lecturas localizadas en el solapamiento de 2 genes que están en hebras diferentes.

![¿Por qué estrangulamiento?](../../images/ref-based/strandness_why.png "Si la información de estrangulamiento se perdió durante la preparación de la biblioteca, la lectura1 se asignará al gen1 localizado en la cadena directa, pero la lectura2 será 'ambigua', ya que puede asignarse al gen1 (cadena directa) o al gen2 (cadena inversa).")

Algunos protocolos de preparación de bibliotecas crean las denominadas bibliotecas de ARN-Seq *stranded* que conservan la información de la cadena ({% cite levin2010comprehensive %} ofrece una excelente descripción general). En la práctica, con los protocolos de RNA-Seq de Illumina es poco probable que se encuentre con todas las posibilidades descritas en este artículo. Lo más probable es que se encuentre con

- Datos de ARN-Seq sin cadena
- Datos de RNA-Seq trenzados generados por el uso de kits especializados de aislamiento de RNA durante la preparación de la muestra

> <details-title>Más detalles sobre el estrangulamiento</details-title>
> 
> ![Relación entre la orientación del ADN y el ARN](../../images/transcriptomics_images/dna_rna.png "Relación entre la orientación del ADN y el ARN")
> 
> La implicación de RNA-Seq stranded es que se puede distinguir si las lecturas derivan de transcritos codificados hacia delante o hacia atrás. En el siguiente ejemplo, los recuentos para el gen Mrpl43 sólo pueden estimarse eficientemente en una librería stranded ya que la mayor parte se solapa con el gen Peo1 en la orientación inversa:
> 
> ![Stranded RNA-Seq data look like this](../../images/ref-based/igv_stranded_screenshot.png "Non-stranded (top) vs. reverse strand-specific (bottom) RNA-Seq read alignment (using IGV, forward mapping reads are red and reverse mapping reads are blue )")
> 
> Dependiendo del enfoque, y de si se realiza secuenciación de extremo único o de extremo pareado, existen múltiples posibilidades sobre cómo interpretar los resultados del mapeo de estas lecturas con el genoma:
> 
> ![Efectos de los tipos de bibliotecas RNA-Seq](../../images/transcriptomics_images/rnaseq_library_type.png "Efectos de los tipos de bibliotecas RNA-Seq (Figura adaptada de la documentación de Sailfish)")
> 
{: .details}

Esta información debería estar incluida en los archivos FASTQ, ¡pregunte en su centro de secuenciación! Si no es así, intente encontrarla en el sitio donde descargó los datos o en la publicación correspondiente.

![¿Cómo estimar la hebra?](../../images/ref-based/strandness_cases.png "En una biblioteca de hebra directa, las lecturas se mapean principalmente en la misma hebra que los genes. En una biblioteca inversa trenzada, las lecturas se sitúan principalmente en la cadena opuesta. Con una biblioteca no trenzada, las lecturas se mapean en los genes en ambas cadenas independientemente de la orientación del gen (Ejemplo para una biblioteca de lectura de un solo extremo).")

Hay 4 formas de estimar la estridencia a partir de los resultados de **STAR** (elija la que prefiera)

1. Podemos hacer una inspección visual de las hebras de lectura en IGV (para conjuntos de datos Paired-end es menos fácil que con una sola lectura y cuando se tienen muchas muestras, esto puede ser doloroso).

   > <hands-on-title>Estimar el estrangulamiento con IGV para una biblioteca de extremos pareados</hands-on-title>
   > 
   > 1. Vuelva a su sesión IGV con el `GSM461177_untreat_paired` BAM abierto.
   > 
   >    > <tip-title>Si no lo tiene</tip-title>
   >    > 
   >    > No hay problema, sólo hay que rehacer los pasos anteriores:
   >    > 
   >    > 1. Iniciar IGV localmente
   >    > 2. Haga clic en la colección `RNA STAR on collection N: mapped.bam` (salida de **RNA STAR** {% icon tool %})
   >    > 3. Expandir el {% icon param-file %} archivo `GSM461177_untreat_paired`.
   >    > 4. Haga clic en `local` en `display with IGV local D. melanogaster (dm6)` para cargar las lecturas en el navegador IGV
   >    > 
   >    {: .tip}
   > 
   > 2. **IGV** {% icon tool %}
   >    1. Zoom a `chr3R:9,445,000-9,448,000` (Cromosoma 3 entre 9,445 kb a 9,448 kb), en la pista `mapped.bam`
   >    2. Haga clic con el botón derecho y seleccione `Color Aligments by` -> `first-in-pair strand`
   >    3. Haga clic con el botón derecho y seleccione `Squished`
   > 
   {: .hands_on}

   > <question-title> Pregunta </question-title>
   > 
   > ![Captura de pantalla de la vista IGV en ps](../../images/ref-based/group_strand_igv_screenshot.png "Captura de pantalla de IGV en ps")
   > 
   > 1. ¿Están las lecturas distribuidas uniformemente entre los 2 grupos (NEGATIVO y POSITIVO)?
   > 2. ¿Cuál es el tipo de strandness de la biblioteca?
   > 
   > > <solution-title></solution-title>
   > > 
   > > 1. Sí, vemos el mismo número de lecturas en ambos grupos.
   > > 2. Significa que la biblioteca no estaba trenzada.
   > > 
   > > > <comment-title>¿Cómo sería si la biblioteca estuviera trenzada?</comment-title>
   > > > 
   > > > ![Captura de pantalla del IGV para hebra frente a no hebra](../../images/ref-based/group_strand_igv_screenshot_RSvsUS.png "Captura de pantalla del IGV para no hebra (arriba) frente a hebra específica inversa (abajo)")
   > > > 
   > > > Observe que no hay ninguna lectura en el grupo POSITIVO para la cadena inversa específica.
   > > {: .comment}
   > {: .solution}
   {: .question}

2. Alternativamente, en lugar de utilizar el BAM se puede utilizar la cobertura de hebra generada por **STAR**. Usando **pyGenomeTracks** podremos visualizar la cobertura en cada hebra para cada muestra. Esta herramienta tiene un montón de parámetros para personalizar sus gráficos.

   > <hands-on-title>Estimar el strandness con pyGenometracks a partir de la cobertura STAR</hands-on-title>
   > 
   > 1. {% tool [pyGenomeTracks](toolshed.g2.bx.psu.edu/repos/iuc/pygenometracks/pygenomeTracks/3.8+galaxy2) %}:
   >    - *"Region of the genome to plot "*: `chr4:540,000-560,000`
   >    - En *""Include tracks in your plot "*:
   >        - {% icon param-repeat %} *"Insert Include tracks in your plot "*
   >            - *"Choose style of the track "*: `Bedgraph track`
   >                - *"Plot title "*: Deje este campo vacío para que el título del gráfico sea el nombre de la muestra.
   >                - {% icon param-collection %} *"Track file(s) bedgraph format "*: Seleccione `RNA STAR on collection N: Coverage Uniquely mapped strand 1`.
   >                - *"Color of track "*: Seleccione un color de su elección, por ejemplo azul
   >                - *"Minimum value "*: `0`
   >                - *"height "*: `3`
   >                - *"Show visualization of data range"*: `Yes`
   >        - {% icon param-repeat %} *"Insert Include tracks in your plot "*
   >            - *"Choose style of the track "*: `Bedgraph track`
   >                - *"Plot title "*: Deje este campo vacío para que el título del gráfico sea el nombre de la muestra.
   >                - {% icon param-collection %} *"Track file(s) bedgraph format "*: Seleccione `RNA STAR on collection N: Coverage Uniquely mapped strand 2`.
   >                - *"Color of track "*: Seleccione un color de su elección distinto del primero, por ejemplo rojo
   >                - *"Minimum value "*: `0`
   >                - *"height "*: `3`
   >                - *"Show visualization of data range "*: `Yes`
   >        - {% icon param-repeat %} *"Insert Include tracks in your plot "*
   >            - *"Choose style of the track "*: `Gene track / Bed track`
   >                - *"Plot title "*: `Genes`
   >                - {% icon param-file %} *"Track file(s) bed or gtf format "*: Seleccione `Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz`
   >                - *"height "*: `5`
   {: .hands_on}

   > <question-title></question-title>
   > 
   > ![pyGenomeTracks](../../images/ref-based/pyGenomeTracks.png "STAR coverage for strand 1 in blue and strand 2 in red")
   > 
   > 1. ¿Qué gen estamos buscando? ¿De qué cadena se trata?
   > 2. ¿Cuál es la cobertura media de cada cadena?
   > 3. ¿Cuál es el nivel de fragilidad de la biblioteca?
   > 
   > > <solution-title></solution-title>
   > > 
   > > 1. Vemos 3 transcritos llamados Thd1-RC, Thd1-RB y Thd1-RA del gen Thd1. El gen se encuentra en la cadena inversa.
   > > 2. La escala pasa a 1,5-2 en los 4 perfiles. La cobertura media debería estar en torno a 1,2-1,5
   > > 3. Deducimos que la biblioteca es no trenzada.
   > > 
   > > > <comment-title>¿Cómo sería si la biblioteca estuviera trenzada?</comment-title>
   > > > 
   > > > ![pyGenomeTracks USvsRS](../../images/ref-based/pyGenomeTracks_USvsRS.png "Cobertura STAR para el strand 1 en azul y el strand 2 en rojo para librerías unstranded y reverse stranded") Observe que la cobertura en la hebra 1 es muy baja para la muestra stranded_PE mientras que el gen es forward. Esto significa que la biblioteca de stranded_PE es de cadena inversa. Por el contrario, para unstranded_PE la escala es comparable para ambas cadenas.
   > > {: .comment}
   > {: .solution}
   {: .question}

3. Se puede utilizar la salida de **STAR** con los recuentos. En efecto, como se ha explicado antes, **STAR** evalúa el número de lecturas en los genes para los tres escenarios posibles: librería no trenzada, trenzado directo o trenzado inverso. La condición que atribuye más lecturas al gen debe ser la condición que corresponde a su biblioteca.

   > <hands-on-title>Estimar el estrangulamiento con recuentos STAR</hands-on-title>
   > 
   > 1. {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.11+galaxy1) %} para agregar los recuentos STAR con los siguientes parámetros:
   >    - En *"Results "*:
   >        - *"Results "*
   >            - *"Which tool was used generate logs?"*: `STAR`
   >                - En *"STAR output "*:
   >                    - {% icon param-repeat %} *"Insert STAR output "*
   >                        - *"Type of STAR output? "*: `Gene counts`
   >                            - {% icon param-collection %} *"STAR gene count output "*: `RNA STAR on collection N: reads per gene` (salida de **RNA STAR** {% icon tool %})
   > 
   {: .hands_on}

   > <question-title></question-title>
   > 
   > 1. ¿Qué porcentaje de lecturas se asignan a genes si la biblioteca es no trenzada/misma trenzada/trenzada inversa?
   > 2. ¿Cuál es el nivel de fragilidad de la biblioteca?
   > 
   > > <solution-title></solution-title>
   > > 
   > > ![STAR Gene counts unstranded](../../images/ref-based/star_gene_counts_unstranded.png "Gen cuenta sin cadena")
   > > ![STAR Gene counts same stranded](../../¡STAR Gene counts reverse stranded](../../images/ref-based/star_gene_counts_reverse.png "Recuento de genes trenzados inversamente")
   > > ![STAR Gene counts reverse stranded](../../images/ref-based/star_gene_counts_reverse.png "Gene counts reverse stranded")
   > > 
   > > 1. Alrededor del 75% de las lecturas se asignan a genes si la biblioteca es no trenzada, mientras que en los otros casos es alrededor del 40%.
   > > 2. Esto sugiere que la biblioteca no está trenzada.
   > > 
   > > > <comment-title>¿Cómo sería si la biblioteca estuviera trenzada?</comment-title>
   > > > 
   > > > ![STAR Gene counts unstranded USvsRS](../../images/ref-based/star_gene_counts_unstranded_USvsRS.png "Recuento de genes no trenzados para bibliotecas no trenzadas y trenzadas inversas")
   > > > ![STAR Gene counts same stranded USvsRS](../../images/ref-based/star_gene_counts_same_USvsRS.png "Gene counts same stranded for unstranded and reverse stranded library") 
   > > > ![STAR Gene counts reverse stranded USvsRS](../../images/ref-based/star_gene_counts_reverse_USvsRS.png "Gene counts reverse stranded for unstranded and reverse stranded library")
   > > > Observe que hay muy pocas lecturas atribuidas a genes para la misma cadena. Las cifras son comparables entre las bibliotecas no trenzadas y las trenzadas inversas, ya que muy pocos genes se solapan en las cadenas opuestas, pero aún así la cifra pasa del 63,6% (no trenzadas) al 65% (trenzadas inversas).
   > > {: .comment}
   > {: .solution}
   {: .question}

4. Otra opción es estimar estos parámetros con una herramienta llamada **Infer Experiment** del conjunto de herramientas RSeQC ({% cite wang2012rseqc %}).

   Esta herramienta toma los archivos BAM del mapeo, selecciona una submuestra de las lecturas y compara sus coordenadas genómicas y hebras con las del Reference gene model ": (de un archivo de anotación). Basándose en la hebra de los genes, puede determinar si la secuenciación es específica de la hebra y, en caso afirmativo, cómo son las hebras de las lecturas (directa o inversa).

   > <hands-on-title>Determinación de la cadena de la biblioteca mediante el Experimento Infer</hands-on-title>
   > 
   > 1. {% tool [Convert GTF to BED12](toolshed.g2.bx.psu.edu/repos/iuc/gtftobed12/gtftobed12/357) %} para convertir el archivo GTF a BED:
   >    - {% icon param-file %} *"GTF File to convert "*: `Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz`
   > 
   >    Es posible que ya haya convertido este archivo `BED12` a partir del conjunto de datos `Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz` anteriormente si realizó la parte detallada sobre las comprobaciones de calidad. En este caso, no es necesario volver a hacerlo una segunda vez
   > 
   > 2. {% tool [Infer Experiment](toolshed.g2.bx.psu.edu/repos/nilesh/rseqc/rseqc_infer_experiment/5.0.3+galaxy0) %} para determinar la robustez de la biblioteca con los siguientes parámetros:
   >    - {% icon param-collection %} *"Input .bam file "*: `RNA STAR on collection N: mapped.bam` (salida de **RNA STAR** {% icon tool %})
   >    - {% icon param-file %} *"Reference gene model "*: Archivo BED12 (salida de **Convert GTF to BED12** {% icon tool %})
   >    - *"Number of reads sampled "*: `200000`
   {: .hands_on}

   {% tool [Infer Experiment](toolshed.g2.bx.psu.edu/repos/nilesh/rseqc/rseqc_infer_experiment/5.0.3+galaxy0) %} tool genera un archivo con información sobre:
    - Biblioteca de extremo pareado o de extremo único
    - Fracción de lecturas que no se han podido determinar
    - 2 líneas
        - Para extremo único
            - `Fraction of reads explained by "++,--"`: fracción de lecturas asignadas a la cadena directa
            - `Fraction of reads explained by "+-,-+"`: fracción de lecturas asignadas a la cadena inversa
        - Para los extremos pareados
            - `Fraction of reads explained by "1++,1--,2+-,2-+"`: fracción de lecturas asignadas a la cadena directa
            - `Fraction of reads explained by "1+-,1-+,2++,2--"`: fracción de lecturas asignadas a la cadena inversa

   Si los dos números de "Fracción de lecturas explicadas por" están próximos entre sí, concluimos que la biblioteca no es un conjunto de datos de cadena específica (o no cadena).

   > <question-title></question-title>
   > 
   > 1. ¿Cuáles son los resultados de "Fracción de las lecturas explicada por" para `GSM461177_untreat_paired`?
   > 2. ¿Cree que el tipo de biblioteca de las 2 muestras es trenzada o no trenzada?
   > 
   > > <solution-title></solution-title>
   > > 
   > > 1. Resultados para `GSM461177_untreat_paired`:
   > > 
   > >    {% snippet faqs/galaxy-es/analysis_results_may_vary.md %}
   > > 
   > >    ```text
   > >    This is PairEnd Data
   > >    Fraction of reads failed to determine: 0.1013
   > >    Fraction of reads explained by "1++,1--,2+-,2-+": 0.4626
   > >    Fraction of reads explained by "1+-,1-+,2++,2--": 0.4360
   > >    ```
   > > 
   > >    de modo que el 46,26% de las lecturas se asignan a la cadena directa y el 43,60% a la cadena inversa.
   > > 
   > > 2. Se encuentran estadísticas similares para `GSM461180_treat_paired`, por lo que la biblioteca parece ser del tipo unstranded para ambas muestras.
   > > 
   > > > <comment-title>¿Cómo sería si la biblioteca estuviera trenzada?</comment-title>
   > > > 
   > > > Siguiendo con el ejemplo de 2 BAM, obtenemos para el no trenzado
   > > > 
   > > > ```text
   > > > This is PairEnd Data
   > > > Fraction of reads failed to determine: 0.0382
   > > > Fraction of reads explained by "1++,1--,2+-,2-+": 0.4847
   > > > Fraction of reads explained by "1+-,1-+,2++,2--": 0.4771
   > > > ```
   > > > 
   > > > Y para la cadena inversa:
   > > > 
   > > > ```text
   > > > This is PairEnd Data
   > > > Fraction of reads failed to determine: 0.0504
   > > > Fraction of reads explained by "1++,1--,2+-,2-+": 0.0061
   > > > Fraction of reads explained by "1+-,1-+,2++,2--": 0.9435
   > > > ```
   > > > 
   > > {: .comment}
   > {: .solution}
   {: .question}

> <details-title>Strandness y configuración del software</details-title>
> 
> Como a veces es bastante difícil averiguar qué ajustes corresponden a los de otros programas, la siguiente tabla puede ser útil para identificar el tipo de biblioteca:
> 
> | Library type         | **Infer Experiment** | **TopHat**       | **HISAT2**         | **HTSeq-count** | **featureCounts** |
> | -------------------- | -------------------- | ---------------- | ------------------ | --------------- | ----------------- |
> | Paired-End (PE) - SF | 1++,1--,2+-,2-+      | FR Second Strand | Second Strand F/FR | yes             | Forward (1)       |
> | PE - SR              | 1+-,1-+,2++,2--      | FR First Strand  | First Strand R/RF  | reverse         | Reverse (2)       |
> | Single-End (SE) - SF | ++,--                | FR Second Strand | Second Strand F/FR | yes             | Forward (1)       |
> | SE - SR              | +-,-+                | FR First Strand  | First Strand R/RF  | reverse         | Reverse (2)       |
> | PE, SE - U           | undecided            | FR Unstranded    | default            | no              | Unstranded (0)    |
> 
{: .details}

## Recuento de lecturas por genes


{% include _includes/cyoa-choices.html option1="featureCounts" option2="STAR" default="featureCounts" text="Para contar el número de lecturas por gen, ofrecemos un tutorial paralelo para los 2 métodos (STAR y featureCounts) que dan resultados muy similares. ¿Qué método prefiere utilizar?" disambiguation="tool"%}

<div class="featureCounts" markdown="1">

Como usted eligió utilizar la opción featureCounts del tutorial, ahora ejecutamos **featureCounts** para contar el número de lecturas por gen anotado.

> <hands-on-title>Contando el número de lecturas por gen anotado</hands-on-title>
> 
> 1. {% tool [featureCounts](toolshed.g2.bx.psu.edu/repos/iuc/featurecounts/featurecounts/2.0.3+galaxy2) %} con los siguientes parámetros para contar el número de lecturas por gen:
>    - {% icon param-collection %} *"Alignment file "*: `RNA STAR on collection N: mapped.bam` (salida de **RNA STAR** {% icon tool %})
>    - *"Specify strand information "*: `Unstranded`
>    - *"Gene annotation file "*: `A GFF/GTF file in your history`
>        - {% icon param-file %} *"Gene annotation file "*: `Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz`
>    - *"GFF feature type filter "*: `exon`
>    - *"GFF gene identifier "*: `gene_id`
>    - *"Output format "*: `Gene-ID "\t" read-count (MultiQC/DESeq2/edgeR/limma-voom compatible)`
>    - *"Create gene-length file "*: `Yes`
>    - *"¿"Does the input have read pairs*: `Yes, paired-end and count them as 1 single fragment`
>    - En *"Read filtering options "*:
>        - *"Minimum quality per read "*: `10`
> 
> 2. {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.11+galaxy1) %} para agregar los informes con los siguientes parámetros:
>    - En *"Results "*:
>        - *"Results "*
>            - *"Which tool was used generate logs?"*: `featureCounts`
>                - {% icon param-collection %} *"Output of FeatureCounts "*: `featureCounts on collection N: Summary` (salida de **featureCounts** {% icon tool %})
> 
>    > <question-title></question-title>
>    > 
>    > 1. ¿Cuántas lecturas se han asignado a un gen?
>    > 2. ¿Cuándo debemos preocuparnos por la tasa de asignación? ¿Qué debemos hacer?
>    > 
>    > > <solution-title></solution-title>
>    > > 
>    > > 1. Alrededor del 63% de las lecturas han sido asignadas a genes: esta cantidad es suficientemente buena.
>    > > 
>    > >    ![featureCounts assignment](../../images/ref-based/featureCounts_assignment_plot.png "Asignaciones con featureCounts")
>    > > 
>    > >    Algunas lecturas no han sido asignadas porque estaban mapeadas; otras no han sido asignadas a ninguna característica o han sido asignadas a características ambiguas.
>    > > 
>    > > 2. Si el porcentaje es inferior al 50%, debería investigar dónde se están mapeando sus lecturas (dentro de genes o no, con IGV) y comprobar que la anotación corresponde a la versión correcta del genoma de referencia.
>    > > 
>    > {: .solution}
>    {: .question}
{: .hands_on}

La salida principal de **featureCounts** es una tabla con los recuentos, es decir, el número de lecturas (o fragmentos en el caso de lecturas paired-end) mapeadas a cada gen (en filas, con su ID en la primera columna) en la anotación proporcionada. **FeatureCount** genera también los conjuntos de datos de salida **feature length**. Necesitaremos este archivo más adelante cuando ejecutemos la herramienta **goseq**.
</div>

<div class="STAR" markdown="1">

Como usted eligió usar la versión STAR del tutorial, usaremos **STAR** para contar las lecturas.

Como se ha escrito anteriormente, durante el mapeo, **STAR** contó las lecturas para cada gen proporcionado en el fichero de anotación de genes (esto se consiguió mediante la opción `Per gene read counts (GeneCounts)`). Sin embargo, esta salida proporciona algunas estadísticas al principio y los recuentos para cada gen dependiendo de la biblioteca (unstranded es la columna 2, stranded forward es la columna 3 y stranded reverse es la columna 4).

> <hands-on-title>Inspeccionar la salida de STAR</hands-on-title>
> 
> 1. Inspeccionar los recuentos de `GSM461177_untreat_paired` en la colección `RNA STAR on collection N: reads per gene`
> 
{: .hands_on}

> 
> <question-title></question-title>
> 
> 1. ¿Cuántas lecturas están sin mapear/multi-mapeadas?
> 2. ¿En qué línea comienza el recuento de genes?
> 3. ¿Cuáles son las diferentes columnas?
> 4. ¿Qué columnas son las más interesantes para nuestro conjunto de datos?
> 
> > <solution-title></solution-title>
> > 
> > 1. Hay 1.190.029 lecturas no mapeadas y 571.324 lecturas multimapeadas.
> > 2. Comienza en la línea 5 con el gen `FBgn0250732`.
> > 3. Hay 4 columnas:
> >    1. ID del gen
> >    2. Recuentos para RNA-seq no trenzado
> >    3. Recuentos de la primera cadena de lectura alineada con el ARN
> >    4. Recuentos de la 2ª cadena de lectura alineada con el ARN
> > 4. Necesitamos la columna Gene ID y la 2ª columna debido a la falta de coherencia de nuestros datos
> > 
> {: .solution}
{: .question}

Reformatearemos la salida de **STAR** para que sea similar a la salida de **featureCounts** (u otros programas de recuento) que sólo tiene 2 columnas, una con IDs y la otra con recuentos.

> <hands-on-title>Reformateo de la salida STAR</hands-on-title>
> 
> 1. {% tool [Select last](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_tail_tool/1.1.0) %} líneas de un conjunto de datos (tail) para eliminar las 4 primeras líneas con los siguientes parámetros:
>    - {% icon param-collection %} *"Text file"*: `RNA STAR on collection N: reads per gene` (salida de **RNA STAR** {% icon tool %})
>    - *"Operation "*: `Keep everything from this line on`
>    - *"Number of lines "*: `5`
> 
> 2. {% tool [Cut](Cut1) %} columnas de una tabla con los siguientes parámetros:
>    - *"Cut columns "*: `c1,c2`
>    - *"Delimited by "*: `Tab`
>    - {% icon param-collection %} *"De "*: `Select last on collection N` (salida de **Select last** {% icon tool %})
> 
> 3. Cambiar el nombre de la colección `FeatureCount-like files`
> 
{: .hands_on}

Más adelante en el tutorial necesitaremos obtener el tamaño de cada gen. Esta es una de las salidas de **FeatureCounts** pero también podemos obtenerla directamente del fichero de anotación del gen. Como esto es bastante largo, recomendamos lanzarlo ahora.

> <hands-on-title>Cómo obtener la longitud del gen</hands-on-title>
> 
> 1. {% tool [Gene length and GC content](toolshed.g2.bx.psu.edu/repos/iuc/length_and_gc_content/length_and_gc_content/0.1.2) %} con los siguientes parámetros:
>    - *"Select a built-in GTF file or one from your history "*: `Use a GTF from history`
>      - {% icon param-file %} *"Select a GTF file "*: `Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz`
>    - *"Analysis to perform "*: `gene lengths only`
> 
>    > <warning-title>Compruebe la versión de la herramienta a continuación</warning-title>
>    > 
>    > Esto sólo funcionará con la versión 0.1.2 o superior
>    > 
>    > {% snippet faqs/galaxy-es/tools_change_version.md %}
>    {: .warning}
> 
{: .hands_on}

</div>

> <question-title></question-title>
> 
> ¿Qué característica tiene más recuentos en ambas muestras? (Sugerencia: utilice la herramienta Ordenar)
> 
> > <solution-title></solution-title>
> > 
> > Para mostrar la característica detectada más abundante, debemos ordenar la tabla de recuentos. Esto se puede hacer de la siguiente manera:
> > 
> > 1. {% tool [Sort](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_sort_header_tool/1.1.1) %} con los siguientes parámetros:
> >    - {% icon param-collection %} *"Sort Query "*: <span class="featureCounts" markdown="1">`featureCounts on collection N: Counts` (salida de **featureCounts** {% icon tool %})</span><span class="STAR" markdown="1">Utilizar la colección `FeatureCount-like files`</span>
> >    - *"Number of header "*: <span class="featureCounts" markdown="1">`1`</span><span class="STAR" markdown="1">`0`</span>
> >    - En *"1: Selecciones de columna "*:
> >      - *"on column "*: `Column: 2`
> > 
> >        Esta columna contiene el número de lecturas = counts
> > 
> >      - *"in "*: `Descending order`
> > 
> > 2. Inspeccionar el resultado
> > 
> >    El resultado de ordenar la tabla en la columna 2 revela que FBgn0284245 es la característica con más recuentos (alrededor de 128.740 en `GSM461177_untreat_paired` y 127.400 en `GSM461180_treat_paired`).
> > 
> >    Comparar diferentes archivos de salida es más fácil si podemos ver más de un conjunto de datos simultáneamente. La función Scratchbook nos permite crear una colección de conjuntos de datos que se mostrarán juntos en la pantalla.
> > 
> >    > <hands-on-title>(Opcional) Ver los recuentos ordenados utilizando el Scratchbook</hands-on-title>
> >    > 
> >    > 1. El **Scratchbook** se activa haciendo clic en el icono de nueve bloques que aparece a la derecha de la barra de menú superior de Galaxy:
> >    > 
> >    >    ![icono scratchbook](../../images/ref-based/menubarWithScratchbook.png "Barra de menú con icono Scratchbook")
> >    > 
> >    > 2. Cuando el Scratchbook está **activado** los conjuntos de datos que se están visualizando (haciendo clic en el icono-ojo) se añaden a la vista Scratchbook:
> >    > 
> >    >    ![Icono de Scratchbook activado](../../images/ref-based/menubarWithScratchbookEnabled.png "Barra de menú con el icono de Scratchbook activado")
> >    > 
> >    > 3. Haga clic en el {% icon galaxy-eye %} (ojo) para ver uno de los archivos de **conteos clasificados**. En lugar de ocupar toda la barra central, la vista del conjunto de datos se muestra ahora superpuesta:
> >    > 
> >    >    ![Scratchbook one dataset shown](../../images/ref-based/scratchbookOneDataset.png "Scratchbook muestra un conjunto de datos superpuesto")
> >    > 
> >    > 4. A continuación, haga clic en el {% icon galaxy-eye %} (ojo) en el **segundo archivo de recuentos ordenados**. El segundo conjunto de datos pasa por encima del primero, pero puede mover la ventana para ver los dos conjuntos de datos uno al lado del otro:
> >    > 
> >    >    ![Scratchbook two datasets shown](../../images/ref-based/scratchbookTwoDatasetsShown.png "Scratchbook mostrando dos conjuntos de datos uno al lado del otro")
> >    > 
> >    > 5. Para **salir** del modo de selección Scratchbook, haga clic de nuevo en el icono **Scratchbook**. Puede decidir cerrar las ventanas o reducirlas para visualizarlas más tarde.
> >    > 
> >    {: .hands_on}
> {: .solution}
{: .question}

Aquí contamos las lecturas asignadas a genes de dos muestras. Es muy interesante volver a realizar el mismo procedimiento en los otros conjuntos de datos, especialmente para comprobar cómo difieren los parámetros dado el diferente tipo de datos (single-end frente a paired-end).

> <hands-on-title>(Opcional) Reejecutar en los otros conjuntos de datos</hands-on-title>
> 
> Puede realizar el mismo proceso en los otros archivos de secuencia disponibles en [Zenodo]({{ page.zenodo_link }}) y en la biblioteca de datos.
> 
> - Datos por pares
>   - `GSM461178_1` y `GSM461178_2` que se pueden etiquetar como `GSM461178_untreat_paired`
>   - `GSM461181_1` y `GSM461181_2` que se pueden etiquetar como `GSM461181_treat_paired`
> - Datos de extremo único
>   - `GSM461176` que puede etiquetar `GSM461176_untreat_single`
>   - `GSM461179` que puede etiquetar `GSM461179_treat_single`
>   - `GSM461182` que puede etiquetar `GSM461182_untreat_single`
> 
> Los enlaces a estos archivos se encuentran a continuación:
> 
> ```text
> {{ page.zenodo_link }}/files/GSM461178_1.fastqsanger
> {{ page.zenodo_link }}/files/GSM461178_2.fastqsanger
> {{ page.zenodo_link }}/files/GSM461181_1.fastqsanger
> {{ page.zenodo_link }}/files/GSM461181_2.fastqsanger
> {{ page.zenodo_link }}/files/GSM461176.fastqsanger
> {{ page.zenodo_link }}/files/GSM461179.fastqsanger
> {{ page.zenodo_link }}/files/GSM461182.fastqsanger
> ```
> 
> Para los datos de un solo extremo, no es necesario aplanar la colección antes del paso **Falco**. Los parámetros de todas las herramientas son los mismos excepto **STAR** para el que puede establecer `Length of the genomic sequence around annotated junctions` a 74 ya que un conjunto de datos tiene lecturas de 75bp (otros son 44bp y 45bp) y **FeatureCount** si sus datos ya no están emparejados.
> 
{: .hands_on}

# Análisis de la expresión génica diferencial

## Identificación de los rasgos expresados diferencialmente

Para poder identificar la expresión génica diferencial inducida por la depleción de PS, todos los conjuntos de datos (3 tratados y 4 no tratados) deben analizarse siguiendo el mismo procedimiento. Para ahorrar tiempo, hemos ejecutado los pasos anteriores por usted. Obtenemos entonces 7 ficheros con los recuentos de cada gen de *Drosophila* para cada muestra.

> <hands-on-title>Importar todos los ficheros de recuento</hands-on-title>
> 
> 1. Crear un **nuevo historial vacío**
> 
>    {% snippet faqs/galaxy-es/histories_create_new.md %}
> 
> 2. Importar los siete archivos de recuento de [Zenodo]({{ page.zenodo_link }}) o la biblioteca de Datos Compartidos:
> 
>    - `GSM461176_untreat_single_featureCounts.counts`
>    - `GSM461177_untreat_paired_featureCounts.counts`
>    - `GSM461178_untreat_paired_featureCounts.counts`
>    - `GSM461179_treat_single_featureCounts.counts`
>    - `GSM461180_treat_paired_featureCounts.counts`
>    - `GSM461181_treat_paired_featureCounts.counts`
>    - `GSM461182_untreat_single_featureCounts.counts`
> 
>    ```text
>    {{ page.zenodo_link }}/files/GSM461176_untreat_single_featureCounts.counts
>    {{ page.zenodo_link }}/files/GSM461177_untreat_paired_featureCounts.counts
>    {{ page.zenodo_link }}/files/GSM461178_untreat_paired_featureCounts.counts
>    {{ page.zenodo_link }}/files/GSM461179_treat_single_featureCounts.counts
>    {{ page.zenodo_link }}/files/GSM461180_treat_paired_featureCounts.counts
>    {{ page.zenodo_link }}/files/GSM461181_treat_paired_featureCounts.counts
>    {{ page.zenodo_link }}/files/GSM461182_untreat_single_featureCounts.counts
>    ```
> 
{: .hands_on}

Se podría pensar que podemos comparar directamente los valores de recuento de los archivos y calcular el grado de expresión génica diferencial. Sin embargo, no es tan sencillo.

Imaginemos que tenemos recuentos de RNA-Seq de 3 muestras para un genoma con 4 genes:

| Gene     | Sample 1 Counts | Sample 2 Counts | Sample 3 Counts |
| -------- | --------------- | --------------- | --------------- |
| A (2kb)  | 10              | 12              | 30              |
| B (4kb)  | 20              | 25              | 60              |
| C (1kb)  | 5               | 8               | 15              |
| D (10kb) | 0               | 0               | 1               |

La muestra 3 tiene más lecturas que las otras réplicas, independientemente del gen. Tiene una profundidad de secuenciación mayor que las otras réplicas. El gen B es el doble de largo que el gen A: podría explicar por qué tiene el doble de lecturas, independientemente de las réplicas.

Por lo tanto, el número de lecturas secuenciadas asignadas a un gen depende de:

- La **profundidad de secuenciación** de las muestras

  Las muestras secuenciadas con mayor profundidad tendrán más lecturas asignadas a cada gen

- La **longitud del gen**

  Los genes más largos tendrán más lecturas asignadas a ellos

Para comparar muestras o expresiones génicas, es necesario normalizar los recuentos de genes. Podríamos utilizar TPM (Transcritos por Kilobase Millón).

> <details-title>RPKM, FPKM y TPM?</details-title>
> 
> Estas tres métricas se utilizan para normalizar las tablas de recuento de:
> 
> - profundidad de secuenciación (la parte del "Millón")
> - longitud del gen (la parte "Kilobase")
> 
> Utilicemos el ejemplo anterior para explicar RPKM, FPKM y TPM.
> 
> Para **RPKM** (lecturas por kilobase millón),
> 
> 1. Calcule el factor de escala "por millón": sume el total de lecturas de una muestra y divida ese número por 1.000.000.
> 
>    | Gene               | Sample 1 Counts | Sample 2 Counts | Sample 3 Counts |
>    | ------------------ | --------------- | --------------- | --------------- |
>    | A (2kb)            | 10              | 12              | 30              |
>    | B (4kb)            | 20              | 25              | 60              |
>    | C (1kb)            | 5               | 8               | 15              |
>    | D (10kb)           | 0               | 0               | 1               |
>    | **Total reads**    | 35              | 45              | 106             |
>    | **Scaling factor** | 3.5             | 4.5             | 10.6            |
> 
>    *Debido a los pequeños valores del ejemplo, utilizamos "por decenas" en lugar de "por millones" y, por tanto, dividimos la suma por 10 en lugar de por 1.000.000.*
> 
> 2. Dividir el número de lecturas por el factor de escala "por millón
> 
>    Esto normaliza la profundidad de secuenciación, dando lecturas por millón (RPM)
> 
>    | Gene     | Sample 1 RPM | Sample 2 RPM | Sample 3 RPM |
>    | -------- | ------------ | ------------ | ------------ |
>    | A (2kb)  | 2.86         | 2.67         | 2.83         |
>    | B (4kb)  | 5.71         | 5.56         | 5.66         |
>    | C (1kb)  | 1.43         | 1.78         | 1.43         |
>    | D (10kb) | 0            | 0            | 0.09         |
> 
>    *En el ejemplo utilizamos el factor de escala "por decenas" y obtenemos lecturas por decenas*
> 
> 3. Divida los valores de RPM por la longitud del gen, en kilobases.
> 
>    | Gene     | Sample 1 RPKM | Sample 2 RPKM | Sample 3 RPKM |
>    | -------- | ------------- | ------------- | ------------- |
>    | A (2kb)  | 1.43          | 1.33          | 1.42          |
>    | B (4kb)  | 1.43          | 1.39          | 1.42          |
>    | C (1kb)  | 1.43          | 1.78          | 1.42          |
>    | D (10kb) | 0             | 0             | 0.009         |
> 
> **FPKM** (Fragmentos por Kilobase Millón) es muy similar a RPKM. RPKM se utiliza para RNA-seq de extremo único, mientras que FPKM se utiliza para RNA-seq de extremo pareado. En el caso del extremo único, cada lectura corresponde a un único fragmento secuenciado. Con RNA-seq paired-end, dos lecturas de un par se mapean a partir de un único fragmento, o si una lectura del par no se mapeó, una lectura puede corresponder a un único fragmento (en caso de que decidiéramos conservarlas). FPKM mantiene un registro de los fragmentos de forma que un fragmento con 2 lecturas se cuenta sólo una vez.
> 
> 
> **TPM** (Transcripciones por Kilobase Millón) es muy similar a RPKM y FPKM, excepto en el orden de la operación
> 
> 1. Divida el número de lecturas por la longitud de cada gen en kilobases
> 
>    Esto da las lecturas por kilobase (RPK).
> 
>    | Gene     | Sample 1 RPK | Sample 2 RPK | Sample 3 RPK |
>    | -------- | ------------ | ------------ | ------------ |
>    | A (2kb)  | 5            | 6            | 15           |
>    | B (4kb)  | 5            | 6.25         | 15           |
>    | C (1kb)  | 5            | 8            | 15           |
>    | D (10kb) | 0            | 0            | 0.1          |
> 
> 2. Calcular el factor de escala "por millón": sumar todos los valores RPK de una muestra y dividir este número por 1.000.000
> 
>    | Gene               | Sample 1 RPK | Sample 2 RPK | Sample 3 RPK |
>    | ------------------ | ------------ | ------------ | ------------ |
>    | A (2kb)            | 5            | 6            | 15           |
>    | B (4kb)            | 5            | 6.25         | 15           |
>    | C (1kb)            | 5            | 8            | 15           |
>    | D (10kb)           | 0            | 0            | 0.1          |
>    | **Total RPK**      | 15           | 20.25        | 45.1         |
>    | **Scaling factor** | 1.5          | 2.03         | 4.51         |
> 
>    *Como en el caso anterior, debido a los pequeños valores del ejemplo, utilizamos "por decenas" en lugar de "por millones" y, por lo tanto, dividimos la suma por 10 en lugar de por 1.000.000.*
> 
> 3. Dividir los valores RPK por el factor de escala "por millón
> 
>    | Gene     | Sample 1 TPM | Sample 2 TPM | Sample 3 TPM |
>    | -------- | ------------ | ------------ | ------------ |
>    | A (2kb)  | 3.33         | 2.96         | 3.33         |
>    | B (4kb)  | 3.33         | 3.09         | 3.33         |
>    | C (1kb)  | 3.33         | 3.95         | 3.33         |
>    | D (10kb) | 0            | 0            | 0.1          |
> 
> A diferencia de RPKM y FPKM, al calcular TPM, normalizamos primero la longitud del gen y luego la profundidad de secuenciación. Sin embargo, los efectos de esta diferencia son bastante profundos, como ya vimos con el ejemplo.
> 
> Las sumas de cada columna son muy diferentes:
> 
> 1. RPKM
> 
>    | Gene      | Sample 1 RPKM | Sample 2 RPKM | Sample 3 RPKM |
>    | --------- | ------------- | ------------- | ------------- |
>    | A (2kb)   | 1.43          | 1.33          | 1.42          |
>    | B (4kb)   | 1.43          | 1.39          | 1.42          |
>    | C (1kb)   | 1.43          | 1.78          | 1.42          |
>    | D (10kb)  | 0             | 0             | 0.009         |
>    | **Total** | 4.29          | 4.5           | 4.25          |
> 
> 2. TPM
> 
>    | Gene      | Sample 1 TPM | Sample 2 TPM | Sample 3 TPM |
>    | --------- | ------------ | ------------ | ------------ |
>    | A (2kb)   | 3.33         | 2.96         | 3.33         |
>    | B (4kb)   | 3.33         | 3.09         | 3.33         |
>    | C (1kb)   | 3.33         | 3.95         | 3.33         |
>    | D (10kb)  | 0            | 0            | 0.1          |
>    | **Total** | 10           | 10           | 10           |
> 
> La suma de todos los TPM en cada muestra es la misma. Esto facilita la comparación de la proporción de lecturas que corresponden a un gen en cada muestra. Por el contrario, con RPKM y FPKM, la suma de las lecturas normalizadas en cada muestra puede ser diferente, y esto hace más difícil comparar muestras directamente.
> 
> En el ejemplo, el TPM para el gen A en la muestra 1 es de 3,33 y en la muestra 2 es de 3,33. La misma proporción del total de lecturas se asigna entonces al gen A en ambas muestras (0,33 aquí). De hecho, la suma de los TPM en ambas muestras suma el mismo número (10 aquí), el denominador requerido para calcular las proporciones es entonces el mismo independientemente de la muestra, y por lo tanto la proporción de lecturas para el gen A (3.33/10 = 0.33) para ambas muestras.
> 
> Con RPKM o FPKM, es más difícil comparar la proporción de lecturas totales porque la suma de lecturas normalizadas en cada muestra puede ser diferente (4,29 para la Muestra 1 y 4,25 para la Muestra 2). Por lo tanto, si el RPKM para el gen A en la Muestra 1 es 1,43 y en la Muestra B es 1,43, no sabemos si la misma proporción de lecturas en la Muestra 1 corresponden al gen A que en la Muestra 2.
> 
> Dado que en RNA-Seq se trata de comparar la proporción relativa de lecturas, TPM parece más apropiado que RPKM/FPKM.
> 
{: .details}

RNA-Seq se utiliza a menudo para comparar un tipo de tejido con otro, por ejemplo, músculo frente a tejido epitelial. Y puede ser que haya muchos genes específicos del músculo transcritos en el músculo pero no en el tejido epitelial. A esto lo llamamos **diferencia en la composición de la biblioteca**.

También es posible observar una diferencia en la composición de la biblioteca en el mismo tipo de tejido tras la eliminación de un factor de transcripción.

Imaginemos que tenemos recuentos de RNA-Seq de 2 muestras (mismo tamaño de biblioteca: 635 reads), para un genoma con 6 genes. Los genes tienen la misma expresión en ambas muestras, excepto uno: sólo la Muestra 1 transcribe el gen D, a un nivel alto (563 reads). Como el tamaño de la biblioteca es el mismo para ambas muestras, la muestra 2 tiene 563 lecturas adicionales que se distribuirán entre los genes A, B, C, E y F.

| Gene      | Sample 1 | Sample 2 |     |
| --------- | -------- | -------- | --- |
| A         | 30       | 235      |     |
| B         | 24       | 188      |     |
| C         | 0        | 0        |     |
| D         | 563      | 0        |     |
| E         | 5        | 39       |     |
| F         | 13       | 102      |     |
| **Total** | 635      | 635      |     |

Como resultado, el recuento de lecturas para todos los genes excepto para los genes C y D es realmente alto en la Muestra 2. No obstante, el único gen expresado diferencialmente es el gen D.

TPM, RPKM o FPKM no tienen en cuenta estas diferencias en la composición de la biblioteca durante la normalización, pero herramientas más complejas, como DESeq2, sí lo hacen.

[**DESeq2**](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) ({% cite love2014moderated %}) es una gran herramienta para tratar con datos de ARN-seq y ejecutar análisis de Expresión Génica Diferencial (EGD). Toma archivos de conteo de lecturas de diferentes muestras, los combina en una gran tabla (con los genes en las filas y las muestras en las columnas) y aplica normalización para **profundidad de secuenciación** y **composición de la biblioteca**. No necesitamos tener en cuenta la normalización de la longitud de los genes porque estamos comparando los recuentos entre grupos de muestras para el mismo gen.

> <details-title>Normalización en DESeq2</details-title>
> 
> Tomemos un ejemplo para ilustrar cómo DESeq2 escala las distintas muestras:
> 
> Gen | Muestra 1 | Muestra 2 | Muestra 3
> A | 0 | 10 | 4
> B | 2 | 6 | 12
> C | 33 | 55 | 200
> 
> El objetivo es calcular un factor de escala para cada muestra, que tenga en cuenta la profundidad de lectura y la composición de la biblioteca.
> 
> 1. Tomar el log$$_e$$ de todos los valores:
> 
>    Gen | log(Muestra 1) | log(Muestra 2) | log(Muestra 3)
>    A | -Inf | 2.3 | 1.4
>    B | 0.7 | 1.8 | 2.5
>    C | 3.5 | 4.0 | 5.3
> 
> 2. Promedio de cada fila:
> 
>    Gen | Media de los valores logarítmicos
>    A | -Inf
>    B | 1,7
>    C | 4,3
> 
>    La media de los valores logarítmicos (también conocida como media geométrica) se utiliza aquí porque no se ve afectada fácilmente por valores atípicos (por ejemplo, el gen C con su valor atípico para la Muestra 3).
> 
> 3. Filtra los genes cuyo valor es infinito.
> 
>    Gen | Media de los valores logarítmicos
>     | 
>    B | 1,7
>    C | 4,3
> 
>    Aquí filtramos los genes sin recuento de lecturas en al menos 1 muestra, por ejemplo, genes que sólo se transcriben en un tejido como el gen D en el ejemplo anterior. Esto ayuda a centrar los factores de escala en genes transcritos a niveles similares, independientemente de la condición.
> 
> 4. Reste el valor logarítmico medio de los recuentos logarítmicos:
> 
>    Gen | log(Muestra 1) | log(Muestra 2) | log(Muestra 3)
>     | | 
>    B | -1.0 | 0.1 | 0.8
>    C | -0.8 | -0.3 | 1.0
> 
>    $$log(\textrm{cuentas del gen X}) - media(\textrm{valores logarítmicos de las cuentas del gen X}) = log(\frac{\textrm{cuentas del gen X}}{\textrm{media del gen X}})$$
> 
>    Este paso compara la relación entre los recuentos de cada muestra y la media de todas las muestras.
> 
> 5. Calcule la mediana de las proporciones de cada muestra:
> 
>    Gen | log(Muestra 1) | log(Muestra 2) | log(Muestra 3)
>     | |
>    B | -1.0 | 0.1 | 0.8
>    C | -0.8 | -0.3 | 1.0
>    **Mediana** | -0.9 | -0.1 | 0.9
> 
>    La mediana se utiliza aquí para evitar que los genes extremos (probablemente los menos frecuentes) influyan demasiado en el valor en una dirección. Ayuda a poner más énfasis en los genes moderadamente expresados.
> 
> 6. Calcule el factor de escala tomando la exponencial de las medianas:
> 
>    Gen | Muestra 1 | Muestra 2 | Muestra 3
>    **Mediana** | -0,9 | -0,1 | 0,9
>    **Factores de escala** | 0,4 | 0,9 | 2,5
> 
> 7. Calcular los recuentos normalizados: dividir los recuentos originales por los factores de escala:
> 
>    Gen | Muestra 1 | Muestra 2 | Muestra 3
>    A | 0 | 11.11 | 1.6
>    B | 5 | 6.67 | 4.8
>    C | 83 | 61.11 | 80
> 
> *Esta explicación es una transcripción y adaptación del [StatQuest video explaining Library Normalization in DESEq2](https://www.youtube.com/watch?v=UFB993xufUU&t=35s)*.
> 
{: .details}

DESeq2 también ejecuta el análisis de Expresión Génica Diferencial (GED), que tiene dos tareas básicas:

- Estimación de la varianza biológica utilizando las réplicas para cada condición
- Estimación de la importancia de las diferencias de expresión entre dos condiciones cualesquiera

Este análisis de expresión se estima a partir de recuentos de lecturas y se intenta corregir la variabilidad en las mediciones utilizando réplicas, que son absolutamente esenciales para obtener resultados precisos. Para su propio análisis, le aconsejamos que utilice al menos 3, pero preferiblemente 5 réplicas biológicas por condición. Es posible tener diferentes números de réplicas por condición.

> <details-title>Réplicas técnicas frente a biológicas</details-title>
> 
> Una réplica técnica es un experimento que se realiza una vez pero se mide varias veces (por ejemplo, secuenciación múltiple de la misma biblioteca). Una réplica biológica es un experimento realizado (y también medido) varias veces.
> 
> En nuestros datos, tenemos 4 réplicas biológicas (aquí llamadas muestras) sin tratamiento y 3 réplicas biológicas con tratamiento (gen *Pasilla* agotado por RNAi).
> 
> Recomendamos combinar las tablas de recuento para diferentes réplicas técnicas (pero no para réplicas biológicas) antes de un análisis de expresión diferencial (véase [documentación DESeq2](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#collapsing-technical-replicates))
> 
{: .details}

A continuación, se pueden incorporar al análisis múltiples factores con varios niveles que describan fuentes de variación conocidas (por ejemplo, tratamiento, tipo de tejido, sexo, lotes), con dos o más niveles que representen las condiciones de cada factor. Tras la normalización, podemos comparar la respuesta de la expresión de cualquier gen a la presencia de distintos niveles de un factor de forma estadísticamente fiable.

En nuestro ejemplo, tenemos muestras con dos factores variables que pueden contribuir a las diferencias en la expresión génica:

- Tratamiento (tratado o no tratado)
- Tipo de secuenciación (paired-end o single-end)

Aquí, el tratamiento es el factor principal que nos interesa. El tipo de secuenciación es información adicional que conocemos sobre los datos y que podría afectar al análisis. El análisis multifactorial nos permite evaluar el efecto del tratamiento, teniendo también en cuenta el tipo de secuenciación.

> <comment-title></comment-title>
> 
> Le recomendamos que añada todos los factores que crea que pueden afectar a la expresión génica en su experimento. Puede ser el tipo de secuenciación como aquí, pero también puede ser la manipulación (si hay diferentes personas involucradas en la preparación de la librería), otros efectos de lote, etc...
> 
{: .comment}

Si sólo tiene uno o dos factores con un número reducido de réplicas biológicas, la configuración básica de **DESeq2** es suficiente. En el caso de una configuración experimental compleja con un gran número de réplicas biológicas, las colecciones basadas en etiquetas son apropiadas. Ambos enfoques dan los mismos resultados. El enfoque basado en etiquetas requiere algunos pasos adicionales antes de ejecutar la herramienta **DESeq2**, pero valdrá la pena cuando se trabaje con una configuración experimental compleja.

{% include _includes/cyoa-choices.html option1="Basic" option2="Tag-based" option3="Collection split" default="Basic" text="¿Qué enfoque prefiere utilizar?" disambiguation="deseq"%}

<div class="Basic" markdown="1">

Ya podemos ejecutar **DESeq2**:

> <hands-on-title>Determinar características expresadas diferencialmente</hands-on-title>
> 
> 1. {% tool [DESeq2](toolshed.g2.bx.psu.edu/repos/iuc/deseq2/deseq2/2.11.40.8+galaxy0) %} con los siguientes parámetros:
>    - *"cómo "*: `Select datasets per level`
>        - En *"Factor "*:
>           - *"Specify a factor name, e.g. effects_drug_x or cancer_markers" "*: `Treatment`
>           - En *"1: Factor level "*:
>               - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'*: `treated`
>               - En *"Count file(s) "*: `Select all the treated count files (GSM461179, GSM461180, GSM461181)`
>           - En *"2: Factor level "*:
>               - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'*: `untreated`
>               - En *"Count file(s) "*: `Select all the untreated count files (GSM461176, GSM461177, GSM461178, GSM461182)`
>       - {% icon param-repeat %} *"Insertion factor
>           - *"Specify a factor name, e.g. effects_drug_x or cancer_markers" "*: `Sequencing`
>               - En *"Factor level "*:
>                    - {% icon param-repeat %} *"Insert Factor level "*
>                        - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'*: `PE`
>                        - En *"Count file(s) "*: `Select all the paired-end count files (GSM461177, GSM461178, GSM461180, GSM461181)`
>                    - {% icon param-repeat %} *"Insert Factor level "*
>                        - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'*: `SE`
>                        - En *"Count file(s) "*: `Select all the single-end count files (GSM461176, GSM461179, GSM461182)`
>    - *"Files have header? "*: `Yes`
>    - *"Choice of input data"*: `Count data (e.g. from HTSeq-count, featureCounts or StringTie)`
>    - En *"Advanced options "*:
>        - *"Use beta priors "*: `Yes`
>    - En *"Output options "*:
>        - *"Output selector "*: `Generate plots for visualizing the analysis results`, `Output normalised counts`
> 
{: .hands_on}

</div>

<div class="Tag-based" markdown="1">

DESeq2 requiere proporcionar para cada factor, recuentos de muestras en cada categoría. Por lo tanto, utilizaremos etiquetas en nuestra colección de recuentos para seleccionar fácilmente todas las muestras que pertenecen a la misma categoría. Para más información sobre formas alternativas de establecer etiquetas de grupo, consulte [este tutorial]({% link topics/galaxy-interface/tutorials/group-tags/tutorial.md %}).

> <hands-on-title>Añada etiquetas a su colección para cada uno de estos factores</hands-on-title>
> 
> 1. Cree una lista de colección con todos estos recuentos que etiquete como `all counts`. Nombre cada elemento de forma que sólo contenga el identificador GSM, el tratamiento y la biblioteca, por ejemplo, `GSM461176_untreat_single`.
> 
>    {% snippet faqs/galaxy-es/collections_build_list.md %}
> 
> 2. {% tool [Extract element identifiers](toolshed.g2.bx.psu.edu/repos/iuc/collection_element_identifiers/collection_element_identifiers/0.0.2) %} con los siguientes parámetros:
>    - {% icon param-collection %} *"Dataset collection "*: `all counts`
> 
>    Ahora extraeremos de los nombres los factores:
> 
> 3. {% tool [Replace Text in entire line](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_line/9.3+galaxy1) %}
>      - {% icon param-file %} *"Fichero a procesar "*: salida de **Extraer identificadores de elementos** {% icon tool %}
>      - En *"Replacement "*:
>         - En *"1: Replacement "*
>            - *"Find pattern "*: `(.*)_(.*)_(.*)`
>            - *"Replace with "*: `\1_\2_\3\tgroup:\2\tgroup:\3`
> 
>    Este paso crea 2 columnas adicionales con el tipo de tratamiento y secuenciación que pueden utilizarse con la herramienta {% tool [Tag elements](__TAG_FROM_FILE__) %} tool
> 
> 4. Cambie el tipo de datos a `tabular`
> 
>    {% snippet faqs/galaxy-es/datasets_change_datatype.md datatype="tabular" %}
> 
> 5. {% tool [Etiquetar elementos](__TAG_FROM_FILE__) %}
>      - {% icon param-collection %} *"Input Collection "*: `all counts`
>      - {% icon param-file %} *"Tag collection elements according to this file "*: salida de **Replace Text** {% icon tool %}
> 
> 6. Inspeccionar la nueva colección
> 
>    > <tip-title>¿No puede ver los cambios?</tip-title>
>    > 
>    > Puede que no lo vea a primera vista, ya que los nombres son los mismos. Sin embargo, si hace clic en uno de ellos y en {% icon galaxy-tags %} **Edit dataset tags**, deberías ver 2 etiquetas que empiezan por 'group:'. Esta palabra clave permitirá utilizar estas etiquetas en **DESeq2**.
>    > 
> > 
>     {: .tip}
> 
{: .hands_on}

Ya podemos ejecutar **DESeq2**:

> <hands-on-title>Determinar características expresadas diferencialmente</hands-on-title>
> 
> 1. {% tool [DESeq2](toolshed.g2.bx.psu.edu/repos/iuc/deseq2/deseq2/2.11.40.8+galaxy0) %} con los siguientes parámetros:
>    - *"how "*: `Select group tags corresponding to levels`
>        - {% icon param-collection %} *"Count file(s) collection "*: salida de **Elementos de etiqueta** {% icon tool %}
>        - En *"Factor "*:
>            - {% icon param-repeat %} *"Insertion factor
>                - *"Specify a factor name, e.g. effects_drug_x or cancer_markers" "*: `Treatment`
>                - En *"Factor level "*:
>                    - {% icon param-repeat %} *"Insert Factor level "*
>                        - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'*: `treated`
>                        - *"Select groups that correspond to this factor level "*: `Tags: treat`
>                    - {% icon param-repeat %} *"Insert Factor level "*
>                        - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'*: `untreated`
>                        - *"Select groups that correspond to this factor level "*: `Tags: untreat`
>            - {% icon param-repeat %} *"Insertion factor
>                - *"Specify a factor name, e.g. effects_drug_x or cancer_markers" "*: `Sequencing`
>                - En *"Factor level "*:
>                    - {% icon param-repeat %} *"Insert Factor level "*
>                        - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'*: `PE`
>                        - *"Select groups that correspond to this factor level "*: `Tags: paired`
>                    - {% icon param-repeat %} *"Insert Factor level "*
>                        - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'*: `SE`
>                        - *"Select groups that correspond to this factor level "*: `Tags: single`
>    - *"Files have header? "*: `Yes`
>    - *"Choice of input data"*: `Count data (e.g. from HTSeq-count, featureCounts or StringTie)`
>    - En *"Advanced options "*:
>        - *"Use beta priors "*: `Yes`
>    - En *"Output options "*:
>        - *"Output selector "*: `Generate plots for visualizing the analysis results`, `Output normalised counts`
> 
{: .hands_on}

</div>
<div class="Collection-split" markdown="1">

DESeq2 requiere proporcionar para cada factor, recuentos de muestras en cada categoría. Por lo tanto, utilizaremos patrones en el nombre de nuestras muestras para seleccionar fácilmente todas las muestras que pertenecen a la misma categoría.

> <hands-on-title>Generar una colección de cada categoría</hands-on-title>
> 
> 1. Cree una lista de colección con todos estos recuentos que etiquete como `all counts`. Nombre cada elemento de forma que sólo contenga el identificador GSM, el tratamiento y la biblioteca, por ejemplo, `GSM461176_untreat_single`.
> 
>    {% snippet faqs/galaxy-es/collections_build_list.md %}
> 
> 2. {% tool [Extract element identifiers](toolshed.g2.bx.psu.edu/repos/iuc/collection_element_identifiers/collection_element_identifiers/0.0.2) %} con los siguientes parámetros:
>    - {% icon param-collection %} *"Dataset collection "*: `all counts`
> 
>    Ahora dividiremos la colección por tratamiento. Tenemos que encontrar un patrón que sólo esté presente en una de las 2 categorías. Utilizaremos la palabra `untreat`:
> 
> 3. {% tool [Search in textfiles](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_grep_tool/9.3+galaxy1) %} (grep) con los siguientes parámetros:
>    - *"Select lines from "*: `Extract element identifiers on data XXX` (salida de **Extract element identifiers** {% icon tool %})
>    - *"that "*: `Match`
>    - *"Regular expression "*: `untreat`
> 
> 4. {% tool [Filter collecion](__FILTER_FROM_FILE__) %} con los siguientes parámetros:
>    - *"Input collection "*: `all counts`
>    - *"How should the elements to remove be determined"*: `Remove if identifiers are ABSENT from file`
>        - *"Filter out identifiers absent from "*: `Search in textfiles on data XXX` (salida de **Búsqueda en archivos de texto** {% icon tool %})
> 
> 5. Renombrar ambas colecciones `untreated` (la colección filtrada) y `treated` (la colección descartada).
> 
> Repetiremos el mismo proceso utilizando `single`
> 
> 6. {% tool [Search in textfiles](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_grep_tool/9.3+galaxy1) %} (grep) con los siguientes parámetros:
>    - *"Select lines from "*: `Extract element identifiers on data XXX` (salida de **Extract element identifiers** {% icon tool %})
>    - *"that "*: `Match`
>    - *"Regular expression "*: `single`
> 
> 7. {% tool [Filter collecion](__FILTER_FROM_FILE__) %} con los siguientes parámetros:
>    - *"Input collection "*: `all counts`
>    - *"How should the elements to remove be determined"*: `Remove if identifiers are ABSENT from file`
>        - *"Filter out identifiers absent from "*: `Search in textfiles on data XXX` (salida de **Búsqueda en archivos de texto** {% icon tool %})
> 
> 8. Renombrar ambas colecciones `single` (la colección filtrada) y `paired` (la colección descartada).
> 
{: .hands_on}

Ya podemos ejecutar **DESeq2**:

> <hands-on-title>Determinar características expresadas diferencialmente</hands-on-title>
> 
> 1. {% tool [DESeq2](toolshed.g2.bx.psu.edu/repos/iuc/deseq2/deseq2/2.11.40.8+galaxy0) %} con los siguientes parámetros:
>    - *"cómo "*: `Select datasets per level`
>        - En *"Factor "*:
>           - *"Specify a factor name, e.g. effects_drug_x or cancer_markers" "*: `Treatment`
>           - En *"1: Factor level "*:
>               - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'*: `treated`
>               - {% icon param-collection %} *"Contar fichero(s) "*: Seleccione la colección `treated`
>           - En *"2: Factor level "*:
>               - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'*: `untreated`
>               - {% icon param-collection %} *"Contar fichero(s) "*: Seleccione la colección `untreated`
>       - {% icon param-repeat %} *"Insertion factor
>           - *"Specify a factor name, e.g. effects_drug_x or cancer_markers" "*: `Sequencing`
>               - En *"Factor level "*:
>                    - {% icon param-repeat %} *"Insert Factor level "*
>                        - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'*: `PE`
>                        - {% icon param-collection %} *"Contar fichero(s) "*: Seleccione la colección `paired`
>                    - {% icon param-repeat %} *"Insert Factor level "*
>                        - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'*: `SE`
>                        - {% icon param-collection %} *"Contar fichero(s) "*: Seleccione la colección `single`
>    - *"Files have header? "*: `Yes`
>    - *"Choice of input data"*: `Count data (e.g. from HTSeq-count, featureCounts or StringTie)`
>    - En *"Advanced options "*:
>        - *"Use beta priors "*: `Yes`
>    - En *"Output options "*:
>        - *"Output selector "*: `Generate plots for visualizing the analysis results`, `Output normalised counts`
> 
{: .hands_on}

</div>

**DESeq2** generó 3 resultados:

- Una tabla con los recuentos normalizados para cada gen (filas) en las muestras (columnas)
- Resumen gráfico de los resultados, útil para evaluar la calidad del experimento:

    1. Un gráfico de las 2 primeras dimensiones de un análisis de componentes principales ([PCA](https://en.wikipedia.org/wiki/Principal_component_analysis)), ejecutado sobre los recuentos normalizados de las muestras

       > <details-title>¿Qué es un PCA?</details-title>
       > 
       > Imaginemos que tenemos algunas botellas de cerveza sobre la mesa. Podemos describir cada cerveza por su color, su espuma, su graduación, etcétera. Podemos componer toda una lista de características diferentes de cada cerveza en una cervecería. Pero muchas de ellas medirán propiedades relacionadas y, por tanto, serán redundantes. Si es así, deberíamos poder resumir cada cerveza con menos características. Esto es lo que hace el ACP o análisis de componentes principales.
       > 
       > Con el ACP, no nos limitamos a seleccionar algunas características interesantes y descartar las demás. En su lugar, construimos algunas características nuevas que resumen bien nuestra lista de cervezas. Estas nuevas características se construyen a partir de las anteriores. Por ejemplo, se puede calcular una nueva característica, como el tamaño de la espuma menos el pH de la cerveza. Se trata de combinaciones lineales.
       > 
       > De hecho, PCA encuentra las mejores características posibles, las que resumen la lista de cervezas. Estas características pueden utilizarse para encontrar similitudes entre las cervezas y agruparlas.
       > 
       > Volviendo a los recuentos de lecturas, el PCA se ejecuta en los recuentos normalizados de todas las muestras. Aquí, nos gustaría describir las muestras basándonos en la expresión de los genes. Así que las características son el número de lecturas mapeadas en cada gen. Las utilizamos, así como combinaciones lineales de las mismas, para representar las muestras y sus similitudes.
       > 
       > *La analogía de la cerveza ha sido adaptada de [una respuesta en StackExchange](https://stats.stackexchange.com/questions/2691/making-sense-of-principal-component-analysis-eigenvectors-eigenvalues)*.
       > 
       {: .details}

       Muestra las muestras en el plano 2D abarcado por sus dos primeros componentes principales. Cada réplica se representa como un punto de datos individual. Este tipo de gráfico es útil para visualizar el efecto global de las covariables experimentales y los efectos de lote.

       > <question-title></question-title>
       > 
       > ![DESeq PCA](../../images/ref-based/deseq2_pca.png "Diagrama de componentes principales de las muestras")
       > 
       > 1. ¿Qué separa la primera dimensión (PC1)?
       > 2. ¿Y la segunda dimensión (PC2)?
       > 3. ¿Qué podemos concluir sobre el diseño DESeq (factores, niveles) que elegimos?
       > 
       > > <solution-title></solution-title>
       > > 
       > > 1. La primera dimensión es separar las muestras tratadas de las no tratadas.
       > > 2. La segunda dimensión consiste en separar los conjuntos de datos de extremo único de los conjuntos de datos de extremo pareado.
       > > 3. Los conjuntos de datos se agrupan siguiendo los niveles de los dos factores. No parece haber ningún efecto oculto en los datos. Si hay una variación no deseada presente en los datos (por ejemplo, efectos de lote), siempre se recomienda corregirla, lo que puede lograrse en DESeq2 incluyendo en el diseño cualquier variable de lote conocida.
       > {: .solution}
       {: .question}

    2. Heatmap de la matriz de distancias muestra-muestra (con agrupación) basada en los recuentos normalizados.

       El mapa de calor ofrece una visión general de las similitudes y disimilitudes entre las muestras: el color representa la distancia entre las muestras. El azul oscuro significa una distancia más corta, es decir, muestras más cercanas dados los recuentos normalizados.

       > <question-title></question-title>
       > 
       > ![Mapa térmico de las distancias entre muestras](../../images/ref-based/deseq2_sample_sample_distance_heatmap.png "Mapa térmico de las distancias entre muestras")
       > 
       > ¿Cómo se agrupan las muestras?
       > 
       > > <solution-title></solution-title>
       > > 
       > > En primer lugar, se agrupan según el tratamiento (primer factor) y, en segundo lugar, según el tipo de secuenciación (segundo factor), como en el gráfico PCA.
       > > 
       > {: .solution}
       {: .question}

    3. Estimaciones de dispersión: estimaciones por genes (negro), los valores ajustados (rojo) y las estimaciones finales máximas a posteriori utilizadas en las pruebas (azul)

       Este gráfico de dispersión es típico, con las estimaciones finales reducidas desde las estimaciones por genes hacia las estimaciones ajustadas. Algunas estimaciones por genes se marcan como valores atípicos y no se reducen hacia el valor ajustado. La cantidad de contracción puede ser mayor o menor de lo que se ve aquí, dependiendo del tamaño de la muestra, el número de coeficientes, la media de la fila y la variabilidad de las estimaciones por genes.

    4. Histograma de valores *p* para los genes en la comparación entre los 2 niveles del 1er factor

    5. Un [MA plot](https://en.wikipedia.org/wiki/MA_plot):

       Muestra la vista global de la relación entre el cambio de expresión de las condiciones (log ratios, M), la fuerza de expresión media de los genes (media media, A) y la capacidad del algoritmo para detectar la expresión génica diferencial. Los genes que superaron el umbral de significación (valor p ajustado < 0,1) están coloreados en rojo.

- Un archivo de resumen con los siguientes valores para cada gen:

    1. Identificadores de genes
    2. Recuentos medios normalizados, promediados en todas las muestras de ambas condiciones
    3. Doble cambio en log2 (logaritmo base 2)

       Los cambios de pliegues log2 se basan en el factor primario de nivel 1 frente al factor de nivel 2, por lo que el orden de entrada de los niveles del factor es importante. Aquí, DESeq2 calcula los cambios de pliegues de las muestras "tratadas" frente a las "no tratadas" a partir del primer factor "Tratamiento", *es decir* los valores corresponden a la regulación al alza o a la baja de los genes en las muestras tratadas.

    4. Estimación del error estándar para la estimación del cambio de pliegue log2
    5. Estadística [Wald](https://en.wikipedia.org/wiki/Wald_test)
    6. Valor *p* para la significación estadística de este cambio
    7. valor *p* ajustado para pruebas múltiples con el procedimiento Benjamini-Hochberg, que controla la tasa de falsos descubrimientos ([FDR](https://en.wikipedia.org/wiki/False_discovery_rate))

  > <tip-title>¿Qué son los valores p y para qué se utilizan?</tip-title>
  > 
  > El valor p es una medida utilizada a menudo para determinar si una observación concreta posee o no significación estadística. En sentido estricto, el valor p es la probabilidad de que los datos hayan podido surgir al azar, suponiendo que la hipótesis nula sea correcta. En el caso concreto de RNA-Seq, la hipótesis nula es que no hay expresión génica diferencial. Por lo tanto, un valor p de 0,13 para un gen concreto indica que, para ese gen, suponiendo que no se exprese de forma diferencial, hay un 13% de posibilidades de que cualquier expresión diferencial aparente pueda producirse simplemente por una variación aleatoria en los datos experimentales.
  > 
  > el 13% sigue siendo bastante elevado, por lo que no podemos estar seguros de que se esté produciendo una expresión génica diferencial. La forma más común en que los científicos utilizan los valores p es establecer un umbral (normalmente 0,05, a veces otros valores como 0,01) y rechazar la hipótesis nula sólo para valores p por debajo de este valor. Así, para los genes con valores p inferiores a 0,05, podemos afirmar con seguridad que la expresión génica diferencial desempeña un papel. Hay que tener en cuenta que cualquier umbral de este tipo es arbitrario y que no hay ninguna diferencia significativa entre un valor p de 0,049 y 0,051, aunque sólo rechacemos la hipótesis nula en el primer caso.
  > 
  > Desgraciadamente, los p-valores a menudo se utilizan mal en la investigación científica, hasta el punto de que Wikipedia ofrece un [artículo dedicado](https://en.wikipedia.org/wiki/Misuse_of_p-values) sobre el tema. Véase también [este artículo](https://fivethirtyeight.com/features/not-even-scientists-can-easily-explain-p-values/) (dirigido a un público general, no científico).
  {: .tip}

Para más información sobre **DESeq2** y sus resultados, puede consultar la [documentación de **DESeq2**](https://www.bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf).

> <question-title></question-title>
> 
> 1. ¿Se expresa de forma diferencial el gen FBgn0003360 debido al tratamiento? En caso afirmativo, ¿en qué medida?
> 2. ¿Está el gen *Pasilla* (ps, FBgn0261552) regulado a la baja por el tratamiento de ARNi?
> 3. También podríamos estar hipotéticamente interesados en el efecto de la secuenciación (u otros factores secundarios en otros casos). ¿Cómo conoceríamos los genes expresados diferencialmente debido al tipo de secuenciación?
> 4. Nos gustaría analizar la interacción entre el tratamiento y la secuenciación. ¿Cómo podríamos hacerlo?
> 
> > <solution-title></solution-title>
> > 
> > 1. FBgn0003360 se expresa de forma diferencial debido al tratamiento: tiene un valor p ajustado significativo ($$2,8 \cdot 10^{-171} << 0,05$$). Se expresa menos (`-` en la columna log2FC) en las muestras tratadas en comparación con las no tratadas, en un factor ~8 ($$2^{log2FC} = 2^{2.99542318410271}$$).
> > 
> > 2. Puede comprobar manualmente si hay `FBgn0261552` en la primera columna o ejecutar {% tool [Filter data on any column using simple expressions](Filter1) %}
> >   - {% icon param-file %} *"Filter "*: el `DESeq2 result file` (salida de **DESeq2** {% icon tool %})
> >   - *"With condition"*: `c1 == "FBgn0261552"`
> > 
> > El cambio de pliegues log2 es negativo, por lo que está regulado a la baja, y el valor p ajustado es inferior a 0,05, por lo que forma parte de los genes con cambios significativos.
> > 
> > 3. DESeq2 en Galaxy devuelve la comparación entre los distintos niveles para el 1er factor, tras la corrección por la variabilidad debida al 2º factor. En nuestro caso actual, tratado contra no tratado para cualquier tipo de secuenciación. Para comparar los tipos de secuenciación, deberíamos ejecutar DESeq2 de nuevo cambiando los factores: el factor 1 (tratamiento) se convierte en el factor 2 y el factor 2 (secuenciación) se convierte en el factor 1.
> > 4. Para añadir la interacción entre dos factores (por ejemplo, tratado para datos paired-end vs no tratado para single-end), debemos ejecutar DESeq2 otra vez pero con un solo factor con los siguientes 4 niveles:
> >    - tratado-PE
> >    - no tratado-PE
> >    - tratado-SE
> >    - sin tratar-SE
> > 
> >    Seleccionando *"Output all levels vs all levels of primary factor (use when you have >2 levels for primary factor) "* a `Yes`, podemos comparar treated-PE vs untreated-SE.
> > 
> {: .solution}
{: .question}

## Anotación de los resultados de DESeq2

El ID de cada gen es algo así como FBgn0003360, que es un ID de la base de datos correspondiente, aquí Flybase ({% cite thurmond2018flybase %}). Estos IDs son únicos, pero a veces preferimos tener los nombres de los genes, incluso si pueden no hacer referencia a un gen único (por ejemplo, duplicado después de la re-anotación). Sin embargo, los nombres de los genes pueden indicar ya una función o ayudar a buscar los candidatos deseados. También nos gustaría mostrar la localización de estos genes dentro del genoma. Podemos extraer dicha información del archivo de anotación que hemos utilizado para el mapeo y el recuento.

> <hands-on-title>Anotación de los resultados de DESeq2</hands-on-title>
> 
> 1. Importar la anotación de genes Ensembl para *Drosophila melanogaster* (`Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz`) del historial anterior, o de la biblioteca de Datos Compartidos o de Zenodo:
> 
>    ```text
>    {{ page.zenodo_link }}/files/Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz
>    ```
> 
> 2. {% tool [Annotate DESeq2/DEXSeq output tables](toolshed.g2.bx.psu.edu/repos/iuc/deg_annotate/deg_annotate/1.1.0) %} with:
>    - {% icon param-file %} *"Tabular output of DESeq2/edgeR/limma/DEXSeq "*: the `DESeq2 result file` (output of **DESeq2** {% icon tool %})
>    - *"Input file type "*: `DESeq2/edgeR/limma`
>    - {% icon param-file %} *"Reference annotation in GFF/GTF format "*: imported gtf `Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz`
> 
{: .hands_on}

La salida generada es una extensión del archivo anterior:

1. Identificadores de genes
2. Recuentos medios normalizados de todas las muestras
3. Cambio de pliegue Log2
4. Estimación del error estándar para la estimación del cambio de pliegue log2
5. Estadística de Wald
6. Valor *p* para el estadístico de Wald
7. *valor p* ajustado para pruebas múltiples con el procedimiento Benjamini-Hochberg para la estadística Wald
8. Cromosoma
9. Inicio
10. Fin
11. Hebra
12. Característica
13. Nombre del gen

> <question-title></question-title>
> 
> 1. ¿Dónde se encuentra el gen más sobreexpresado?
> 2. ¿Cuál es el nombre del gen?
> 3. ¿Dónde se encuentra el gen *Pasilla* (FBgn0261552)?
> 
> > <solution-title></solution-title>
> > 
> > 1. FBgn0025111 (el gen mejor clasificado con el valor log2FC positivo más alto) se encuentra en el reverso del cromosoma X, entre 10.778.953 pb y 10.786.907 pb.
> > 2. De la tabla, obtenemos el símbolo del gen: Ant2. Después de algunas búsquedas en las [bases de datos biológicas en línea](https://www.ncbi.nlm.nih.gov/gene/32008), encontramos que Ant2 corresponde a la adenina nucleótido translocasa 2.
> > 3. El gen *Pasilla* se encuentra en la cadena anterior del cromosoma 3R, entre 9.417.939 pb y 9.455.500 pb.
> > 
> {: .solution}
{: .question}

La tabla de anotaciones no contiene nombres de columnas, lo que dificulta su lectura. Nos gustaría añadirlos antes de continuar.

> <hands-on-title>Añadir nombres de columnas</hands-on-title>
> 
> 1. Cree un nuevo archivo (`header`) a partir de lo siguiente (línea de encabezado de la salida DESeq2)
> 
>    ```text
>    GeneID	Base mean	log2(FC)	StdErr	Wald-Stats	P-value	P-adj	Chromosome	Start	End	Strand	Feature	Gene name
>    ```
> 
>    {% snippet faqs/galaxy-es/datasets_create_new_file.md name="header" format="tabular" %}
> 
> 2. {% tool [Concatenate datasets](cat1) %} para añadir esta línea de cabecera a la salida **Annotate**:
>    - {% icon param-file %} *"Concatenate Dataset "*: el conjunto de datos `header`
>    - *"Dataset "*
>       - Haga clic en {% icon param-repeat %} *"Insertar conjunto de datos "*
>         - {% icon param-file %} *"select "*: salida de **Annotate** {% icon tool %}
> 
> 3. Cambie el nombre de la salida a `Annotated DESeq2 results`
> 
{: .hands_on}

## Extracción y anotación de genes expresados diferencialmente

Ahora queremos extraer los genes más expresados diferencialmente debido al tratamiento con un cambio de pliegue > 2 (o < 1/2).

> <hands-on-title>Extraer los genes más expresados diferencialmente</hands-on-title>
> 
> 1. {% tool [Filter data on any column using simple expressions](Filter1) %} para extraer los genes con un cambio significativo en la expresión génica (valor *p* ajustado inferior a 0,05) entre las muestras tratadas y no tratadas:
>    - {% icon param-file %} *"Filter "*: `Annotated DESeq2 results`
>    - *"With following condition"*: `c7<0.05`
>    - *"Number of header lines to skip "*: `1`
> 
> 2. Renombrar la salida `Genes with significant adj p-value`.
> 
>    > <question-title></question-title>
>    > 
>    > ¿Cuántos genes presentan un cambio significativo en la expresión génica entre estas condiciones?
>    > 
>    > > <solution-title></solution-title>
>    > > 
>    > > Obtenemos 966 (967 líneas incluyendo una cabecera) genes (4,04%) con un cambio significativo en la expresión génica entre las muestras tratadas y no tratadas.
>    > > 
>    > {: .solution}
>    {: .question}
> 
> > <comment-title></comment-title>
> > 
> > El archivo con los resultados filtrados de forma independiente puede utilizarse para análisis posteriores, ya que excluye genes con sólo unos pocos recuentos de lecturas, ya que estos genes no se considerarán expresados de forma diferencial significativa.
> > 
> {: .comment}
> 
> Ahora seleccionaremos sólo los genes con un cambio de pliegue (FC) > 2 o FC < 0,5. Tenga en cuenta que el archivo de salida de DESeq2 contiene $$log_{2} FC$$, en lugar del propio FC, por lo que filtraremos por $$abs(log_{2} FC) > 1$$ (lo que implica FC > 2 o FC < 0,5).
> 
> 3. {% tool [Filter data on any column using simple expressions](Filter1) %} para extraer genes con un $$abs(log_{2} FC) > 1$$:
>    - {% icon param-file %} *"Filter "*: `Genes with significant adj p-value`
>    - *"With following condition"*: `abs(c3)>1`
>    - *"Number of header lines to skip "*: `1`
> 
> 4. Renombrar la salida `Genes with significant adj p-value & abs(log2(FC)) > 1`.
> 
>    > <question-title></question-title>
>    > 
>    > 1. ¿Cuántos genes se han conservado?
>    > 2. ¿Se puede encontrar el gen *Pasilla* (ps, FBgn0261552) en esta tabla?
>    > 
>    > > <solution-title></solution-title>
>    > > 
>    > > 1. Obtenemos 113 genes (114 líneas incluida la cabecera), es decir, el 11,79% de los genes con expresión diferencial significativa.
>    > > 2. El gen *Pasilla* puede encontrarse con una búsqueda rápida (o incluso utilizando {% tool [Filter data on any column using simple expressions](Filter1) %} )
>    > {: .solution}
>    {: .question}
{: .hands_on}

Ahora tenemos una tabla con 113 líneas y una cabecera correspondiente a los genes más expresados diferencialmente. Para cada gen, tenemos su ID, sus recuentos medios normalizados (promediados sobre todas las muestras de ambas condiciones), su $$log_{2} FC$$ y otra información que incluye el nombre y la posición del gen.

## Visualización de la expresión de los genes expresados diferencialmente

Podríamos trazar el $$log_{2} FC$$ para los genes extraídos, pero aquí nos gustaría ver un mapa de calor de la expresión de estos genes en las diferentes muestras. Así que tenemos que extraer los recuentos normalizados para estos genes.

Procedemos en varios pasos:

- Extraiga y trace los recuentos normalizados de estos genes para cada muestra con un mapa de calor, utilizando el archivo de recuentos normalizados generado por DESeq2
- Calcular, extraer y representar gráficamente la puntuación Z de los recuentos normalizados

> <comment-title>Tutoriales avanzados sobre visualización</comment-title>
> 
> En este tutorial, explicaremos rápidamente algunas posibilidades de visualización. Para más detalles, eche un vistazo a los tutoriales adicionales sobre visualización de resultados de RNA-Seq:
> 
> - [Visualización de los resultados de RNA-Seq con heatmap2]({% link topics/transcriptomics/tutorials/rna-seq-viz-with-heatmap2/tutorial.md %})
> - [Visualización de resultados de RNA-Seq con Volcano Plot]({% link topics/transcriptomics/tutorials/rna-seq-viz-with-volcanoplot/tutorial.md %})
> 
{: .comment}

### Visualización de los recuentos normalizados

Para extraer los recuentos normalizados de los genes interesantes, unimos la tabla de recuentos normalizados generada por DESeq2 con la tabla que acabamos de generar. Conservaremos entonces sólo las columnas correspondientes a los recuentos normalizados.

> <hands-on-title>Extraer los recuentos normalizados de los genes más diferencialmente expresados</hands-on-title>
> 
> 1. {% tool [Join two Datasets side by side on a specified field](join1) %} con los siguientes parámetros:
>    - {% icon param-file %} *"Join "*: el archivo `Normalized counts` (salida de **DESeq2** {% icon tool %})
>    - *"using column" "*: `Column: 1`
>    - {% icon param-file %} *"with "*: `Genes with significant adj p-value & abs(log2(FC)) > 1`
>    - *"and column "*: `Column: 1`
>    - *"Keep lines of first input that do not join with second input "*: `No`
>    - *"Keep header files "*: `Yes`
> 
>    El archivo generado tiene más columnas de las que necesitamos para el mapa de calor: recuentos medios normalizados, $$log_{2} FC$$ y otra información de anotación. Tenemos que eliminar las columnas adicionales.
> 
> 2. {% tool [Cut](Cut1) %} columnas de una tabla con los siguientes parámetros para extraer las columnas con los IDs de los genes y los recuentos normalizados:
>    - *"Cut columns "*: `c1-c8`
>    - *"Delimited by "*: `Tab`
>    - {% icon param-file %} *"From "*: el conjunto de datos unido (resultado de **Join two Datasets** {% icon tool %})
> 
> 3. Cambie el nombre de la salida a `Normalized counts for the most differentially expressed genes`
> 
{: .hands_on}

Ahora tenemos una tabla con 114 líneas (los 113 genes más expresados diferencialmente y una cabecera) y los recuentos normalizados de estos genes en las 7 muestras.

> <hands-on-title>Planifica el mapa de calor de los recuentos normalizados de estos genes para las muestras</hands-on-title>
> 
> 1. {% tool [heatmap2](toolshed.g2.bx.psu.edu/repos/iuc/ggplot2_heatmap2/ggplot2_heatmap2/3.1.3.1+galaxy0) %} para trazar el mapa de calor:
>    - {% icon param-file %} *"Input should have column headers "*: `Normalized counts for the most differentially expressed genes`
>    - *"Data transformation "*: `Log2(value+1) transform my data`
>    - *"Enable data clustering "*: `Yes`
>    - *"Labeling columns and rows "*: `Label columns and not rows`
>    - *"Type of colormap to use "*: `Gradient with 2 colors`
> 
{: .hands_on}

Debería obtener algo similar a:

![Mapa de calor con los recuentos normalizados de los genes más expresados diferencialmente](../../images/ref-based/heatmap2_normalized_counts.png "Recuentos normalizados de los genes más expresados diferencialmente")

> <question-title></question-title>
> 
> 1. ¿Qué representa el eje X del mapa de calor? ¿Y el eje Y?
> 2. ¿Observa algo en la agrupación de las muestras y los genes?
> 3. ¿Qué cambia si regenera el mapa de calor, esta vez seleccionando `Plot the data as it is` en *"Data transformation "*?
> 4. ¿Por qué no podemos utilizar `Log2(value) transform my data` en *"Data transformation "*?
> 5. ¿Cómo podría generar un mapa térmico de recuentos normalizados para todos los genes regulados al alza con un cambio de pliegue > 2?
> 
> > <solution-title></solution-title>
> > 
> > 1. El eje X muestra las 7 muestras, junto con un dendrograma que representa la similitud entre sus niveles de expresión génica. El eje Y muestra los 113 genes expresados diferencialmente, igualmente con un dendrograma que representa la similitud entre los niveles de expresión génica.
> > 2. Las muestras se agrupan por tratamiento.
> > 3. La escala cambia y sólo vemos unos pocos genes.
> > 4. Porque la expresión normalizada del gen `FBgn0013688` en `GSM461180_treat_paired` está en `0`.
> > 5. Extraer los genes con $$log_{2} FC$$ > 1 (filtro para genes con `c3>1` en el resumen de los genes expresados diferencialmente) y ejecute **heatmap2** {% icon tool %} en la tabla generada.
> > 
> {: .solution}
{: .question}

### Visualización de la puntuación Z

Para comparar la expresión génica entre muestras, también podríamos utilizar la puntuación Z, que a menudo se representa en las publicaciones.

La puntuación Z da el número de desviaciones estándar que un valor se aleja de la media de todos los valores del mismo grupo, aquí el mismo gen. Una puntuación Z de -2 para el gen X en la muestra A significa que este valor es 2 desviaciones estándar inferior a la media de los valores del gen X en todas las muestras (A, B, C, etc.).

La puntuación Z $$z_{i,j}$$ para un gen $$i$$ en una muestra $$j$$ dado el recuento normalizado $$x_{i,j}$$ se calcula como $$z_{i,j} = \frac{x_{i,j}- \overline{x_i}}{s_i}$$ con $$\overline{x_i}$$ la media y $$s_i$$ la desviación estándar de los recuentos normalizados para el gen $$i$$ en todas las muestras.

> <details-title>Calcular la puntuación Z de todos los genes</details-title>
> 
> A menudo necesitamos la puntuación Z para algunas visualizaciones. Para calcular la puntuación Z, dividimos el proceso en 2 pasos:
> 
> 1. Reste cada valor por la media de los valores de la fila (es decir, $$x_{i,j}- \overline{x_i}$$) utilizando la tabla de recuento normalizado
> 2. Dividir los valores anteriores por la desviación estándar de los valores de la fila, utilizando 2 tablas (los recuentos normalizados y la tabla calculada en el paso anterior)
> 
> > <hands-on-title>Calcular la puntuación Z de todos los genes</hands-on-title>
> > 
> > 1. {% tool [Table Compute](toolshed.g2.bx.psu.edu/repos/iuc/table_compute/table_compute/1.2.4+galaxy0) %} con los siguientes parámetros para > restar primero los valores medios por fila
> >    - *"Input Single or Multiple Tables "*: `Single Table`
> >      - {% icon param-file %} *"Table "*: `Normalized counts file on data ... and others` (salida de **DESeq2** {% icon tool %})
> >      - *"Type of table operation "*: `Perform a full table operation`
> >        - *"Operation "*: `Custom`
> >          - *"Custom expression on 'table', along 'axis' (0 or 1) "*: `table.sub(table.mean(1), 0)`
> > 
> >            La expresión `table.mean(1)` calcula la media de cada fila (aquí los genes) y `table.sub(table.mean(1), 0)` resta a cada valor la media de la fila (calculada con `table.mean(1)`)
> > 
> > 2. {% tool [Table Compute](toolshed.g2.bx.psu.edu/repos/iuc/table_compute/table_compute/1.2.4+galaxy0) %} con los siguientes parámetros:
> >    - *"Input Single or Multiple Tables "*: `Multiple Table`
> >      - Haga clic en {% icon param-repeat %} *"Insert tables "*
> >      - En *"1: Tablas "*:
> >        - {% icon param-file %} *"Table "*: `Normalized counts file on data ... and others` (salida de **DESeq2** {% icon tool %})
> >      - Haga clic en {% icon param-repeat %} *"Insert tables "*
> >      - En *"2: Tablas "*:
> >        - {% icon param-file %} *"Table "*: salida del primer **Table Compute** {% icon tool %}
> >      - *"Custom expression on 'tableN*: `table2.div(table1.std(1),0)`
> > 
> >        La expresión `table1.std(1)` calcula las desviaciones estándar de cada fila de la 1ª tabla (recuentos normalizados) y `table2.div` divide los valores de la 2ª tabla (calculados previamente) por estas desviaciones estándar.
> > 
> > 3. Cambie el nombre de la salida a `Z-scores`
> > 4. Inspeccione el archivo de salida
> > 
> {: .hands_on}
> 
> Ahora tenemos una tabla con la puntuación Z de todos los genes de las 7 muestras.
> 
> > <question-title></question-title>
> > 
> > 1. ¿Cuál es el rango para la puntuación Z?
> > 2. ¿Por qué algunas filas están vacías?
> > 3. ¿Qué podemos decir de las puntuaciones Z de los genes expresados diferencialmente (por ejemplo, `FBgn0037223`)?
> > 4. ¿Podemos utilizar la puntuación Z para estimar la intensidad de la expresión diferencial de un gen?
> > 
> > > <solution-title></solution-title>
> > > 
> > > 1. La puntuación Z oscila entre -3 desviaciones estándar y +3 desviaciones estándar. Puede situarse en una curva de distribución normal: -3 es el extremo izquierdo de la curva de distribución normal y +3 el extremo derecho de la curva de distribución normal
> > > 2. Si todos los recuentos son idénticos (normalmente a 0), la desviación estándar es 0, la puntuación Z no puede calcularse para estos genes.
> > > 3. Cuando un gen se expresa diferencialmente entre dos grupos (aquí tratado y no tratado), las puntuaciones Z para este gen serán (mayoritariamente) positivas para las muestras de un grupo y (mayoritariamente) negativas para las muestras del otro grupo.
> > > 4. La puntuación Z es una relación señal-ruido. Las puntuaciones Z absolutas grandes, es decir, los valores positivos o negativos grandes, no son una estimación directa del efecto, es decir, de la fuerza de la expresión diferencial. Una misma puntuación Z grande puede tener diferentes significados, dependiendo del ruido:
> > >    - con mucho ruido: un efecto muy grande
> > >    - con algo de ruido: un efecto bastante grande
> > >    - con poco ruido: un efecto más bien pequeño
> > >    - casi sin ruido: un efecto minúsculo
> > > 
> > >    El problema es que aquí el "ruido" no es sólo ruido de la medida. También puede estar relacionado con el "rigor" del control de la regulación de los genes. Los genes no estrechamente controlados, es decir, cuya expresión puede variar en un amplio rango a lo largo de las muestras, pueden ser inducidos o reprimidos considerablemente. Su puntuación Z absoluta será pequeña, ya que las variaciones entre muestras son grandes. Por el contrario, los genes que están estrechamente controlados pueden tener sólo cambios muy pequeños en su expresión, sin ningún impacto biológico. La puntuación Z absoluta será grande para estos genes.
> > > 
> > {: .solution}
> {: .question}
{: .details}

Ahora nos gustaría trazar un mapa de calor para las puntuaciones Z:

![Mapa de calor con las puntuaciones Z de los genes con mayor expresión diferencial](../../images/ref-based/z-score-heatmap.png "Puntuaciones Z de los genes con mayor expresión diferencial")

> <hands-on-title>Planificar la puntuación Z de los genes más expresados diferencialmente</hands-on-title>
> 
> 1. {% tool [heatmap2](toolshed.g2.bx.psu.edu/repos/iuc/ggplot2_heatmap2/ggplot2_heatmap2/3.1.3.1+galaxy0) %} para trazar el mapa de calor:
>    - {% icon param-file %} *"Input should have column headers "*: `Normalized counts for the most differentially expressed genes`
>    - *"Data transformation "*: `Plot the data as it is`
>    - *"Compute z-scores prior to clustering "*: `Compute on rows`
>    - *"Enable data clustering "*: `Yes`
>    - *"Labeling columns and rows "*: `Label columns and not rows`
>    - *"Type of colormap to use "*: `Gradient with 3 colors`
> 
{: .hands_on}

# Análisis de enriquecimiento funcional de los genes DE

Hemos extraído los genes que se expresan diferencialmente en las muestras tratadas (con genes PS agotados) en comparación con las muestras no tratadas. Ahora, nos gustaría saber si los genes expresados diferencialmente son transcritos enriquecidos de genes que pertenecen a categorías más comunes o específicas para identificar las funciones biológicas que podrían verse afectadas.

## Análisis de ontología génica

El análisis [Gene Ontology (GO)](http://www.geneontology.org/) se utiliza ampliamente para reducir la complejidad y resaltar los procesos biológicos en los estudios de expresión de todo el genoma. Sin embargo, los métodos estándar dan resultados sesgados en los datos de RNA-Seq debido a la sobredetección de la expresión diferencial para transcritos largos y altamente expresados.

[**goseq**](https://bioconductor.org/packages/release/bioc/vignettes/goseq/inst/doc/goseq.pdf) ({% cite young2010gene %}) proporciona métodos para realizar análisis GO de datos RNA-Seq teniendo en cuenta el sesgo de longitud. **goseq** también podría aplicarse a otros análisis basados en categorías de datos RNA-Seq, como el análisis de vías KEGG, como se discute en una sección posterior.

**goseq** necesita 2 archivos de entrada:

- Un archivo tabular con los genes expresados diferencialmente de todos los genes analizados en el experimento RNA-Seq con 2 columnas:
  - los Gene IDs (únicos dentro del archivo), en mayúsculas
  - un booleano que indica si el gen se expresa de forma diferencial o no (`True` si se expresa de forma diferencial o `False` si no)
- Un archivo con información sobre la longitud de un gen para corregir el posible sesgo de longitud en genes expresados diferencialmente

> <hands-on-title>Preparar el primer conjunto de datos para goseq</hands-on-title>
> 
> 1. {% tool [Compute](toolshed.g2.bx.psu.edu/repos/devteam/column_maker/Add_a_column1/2.0) %} en filas con los siguientes parámetros:
>    - {% icon param-file %} *"Input file"*: el `DESeq2 result file` (salida de **DESeq2** {% icon tool %})
>    - En *"Expresiones"*:
>      - {% icon param-text %} *"Add expression"*: `bool(float(c7)<0.05)`
>      - {% icon param-select %} *"Mode of the operation?"*: `Append`
>    - En *"Tratamiento de errores"*:
>      - {% icon param-toggle %} *"Autodetect column types"*: `No`
>      - {% icon param-select %} *"If an expression cannot be computed for a row"*: `Fill in a replacement value`
>      - {% icon param-select %} *"Replacement value"*: `False`
> 
> 2. {% tool [Cut](Cut1) %} columnas de una tabla con los siguientes parámetros:
>    - *"Cut columns"*: `c1,c8`
>    - *"Delimited by"*: `Tab`
>    - {% icon param-file %} *"From"*: la salida de **Compute** {% icon tool %}
> 
> 3. {% tool [Change Case](ChangeCase) %} con
>    - {% icon param-file %} *"From"*: la salida del **Cut** anterior {% icon tool %}
>    - *"Change case of columns"*: `c1`
>    - *"Delimited by"*: `Tab`
>    - *"To"*: `Upper case`
> 
> 4. Cambie el nombre de la salida a `Gene IDs and differential expression`
{: .hands_on}

Acabamos de generar la primera entrada para **goseq**. Como segunda entrada para **goseq** necesitamos las longitudes de los genes. Podemos usar aquí las longitudes de los genes generadas por **featureCounts** o **Gene length and GC content** y formatear los IDs de los genes.

> <hands-on-title>Preparar el archivo de longitud del gen</hands-on-title>
> 
> <div class="featureCounts" markdown="1">
> 1. Copy the feature length collection previously generated by **featureCounts** {% icon tool %} into this history
> 
> {% snippet faqs/galaxy-es/histories_copy_dataset.md %}
> 
> 2. {% tool [Extraer conjunto de datos](__EXTRACT_DATASET__) %} con:
>    - {% icon param-collection %} *"Input List "*: `featureCounts on collection N: Feature lengths`
>    - *"¿Cómo debe seleccionarse un conjunto de datos?*: `The first dataset`
> 
> </div>
> 
> <div class="STAR" markdown="1">
> 1. Copy the output of **Gene length and GC content** {% icon tool %} (`Gene length`) into this history
> 
> {% snippet faqs/galaxy-es/histories_copy_dataset.md %}
> </div>
> 
> 2. {% tool [Change Case](ChangeCase) %} con los siguientes parámetros:
> 
>    - {% icon param-file %} *"De "*: <span class="featureCounts" markdown="1">`GSM461177_untreat_paired` (salida de **Extract Dataset** {% icon tool %})</span><span class="STAR" markdown="1">`Gene length`</span>
>    - *"Change case of columns "*: `c1`
>    - *"Delimited by"*: `Tab`
>    - *"To "*: `Upper case`
> 
> 3. Cambie el nombre de la salida a `Gene IDs and length`
{: .hands_on}

Ya tenemos los dos archivos de entrada necesarios para goseq.

> <hands-on-title>Realizar análisis GO</hands-on-title>
> 
> 1. {% tool [goseq](toolshed.g2.bx.psu.edu/repos/iuc/goseq/goseq/1.50.0+galaxy0) %} with
>    - *"Differentially expressed genes file "*: `Gene IDs and differential expression`
>    - *"Gene lengths file "*: `Gene IDs and length`
>    - *"Gene categories "*: `Get categories`
>       - *"Select a genome to use "*: `Fruit fly (dm6)`
>       - *"Select Gene ID format "*: `Ensembl Gene ID`
>       - *"Select one or more categories "*: `GO: Cellular Component`, `GO: Biological Process`, `GO: Molecular Function`
>    - En *"Output options "*
>      - *"Output Top GO terms plot? "*: `Yes`
>      - *"Extract the DE genes for the categories (GO/KEGG terms)? "*: `Yes`
{: .hands_on}

**goseq** genera con estos parámetros 3 salidas:

1. Una tabla (`Ranked category list - Wallenius method`) con las siguientes columnas para cada término GO:

    1. `category`: Categoría GO
    2. `over_rep_pval`: *valor p* para la sobrerrepresentación del término en los genes expresados diferencialmente
    3. `under_rep_pval`: *valor p* para la infrarrepresentación del término en los genes expresados diferencialmente
    4. `numDEInCat`: número de genes expresados diferencialmente en esta categoría
    5. `numInCat`: número de genes en esta categoría
    6. `term`: detalle del término
    7. `ontology`: MF (Función Molecular - actividades moleculares de los productos génicos), CC (Componente Celular - donde los productos génicos son activos), BP (Proceso Biológico - vías y procesos más amplios compuestos por las actividades de múltiples productos génicos)
    8. `p.adjust.over_represented`: *valor p* para la sobrerrepresentación del término en los genes con expresión diferencial, ajustado para pruebas múltiples con el procedimiento Benjamini-Hochberg
    9. `p.adjust.under_represented`: *valor p* para la infrarrepresentación del término en los genes con expresión diferencial, ajustado para pruebas múltiples con el procedimiento Benjamini-Hochberg

   Para identificar categorías significativamente enriquecidas/no enriquecidas por debajo de algún valor p de corte, es necesario utilizar el valor *p* ajustado.

   > <question-title></question-title>
   > 
   > 1. ¿Cuántos términos GO están sobrerrepresentados con un valor P ajustado < 0,05? ¿Cuántos están infrarrepresentados?
   > 2. ¿Cómo se dividen los términos GO sobrerrepresentados en MF, CC y BP? ¿Y para los términos GO infrarrepresentados?
   > 
   > > <solution-title></solution-title>
   > > 
   > > 1. 60 términos GO (0,50%) están sobrerrepresentados y 7 (0,07%) infrarrepresentados.
   > > 
   > >    {% tool [Filter data on any column using simple expressions](Filter1) %} on c8 (adjusted p-value for over-represented GO terms) and c9 (adjusted p-value for under-represented GO terms)
   > > 
   > > 2. Para sobre-representados, 50 BP, 5 CC y 5 MF y para sub-representados, 5 BP, 2 CC y 0 MF
   > > 
   > >    {% tool [Agrupar datos](Grouping1) %} en la columna 7 (categoría) y recuento en la columna 1 (IDs)
   > > 
   > {: .solution}
   {: .question}

2. Un gráfico con los 10 términos GO más sobrerrepresentados

   > <question-title></question-title>
   > 
   > ![Top over-represented GO terms](../../images/ref-based/top_over-represented_go_terms.png)
   > 
   > ¿Qué es el eje x? ¿Cómo se calcula?
   > 
   > > <solution-title></solution-title>
   > > 
   > > El eje x es el porcentaje de genes de la categoría que se han identificado como expresados diferencialmente: $$100 \times \frac{numDEInCat}{numInCat}$$
   > > 
   > {: .solution}
   {: .question}

3. Una tabla con los genes expresados diferencialmente (de la lista que proporcionamos) asociados a los términos GO (`DE genes for categories (GO/KEGG terms)`)

> <comment-title>Tutorial avanzado sobre análisis de enriquecimiento</comment-title>
> 
> En este tutorial, hemos cubierto el análisis de enriquecimiento GO con **goseq**. Para aprender otros métodos y herramientas para el análisis de enriquecimiento de conjuntos de genes, por favor eche un vistazo al tutorial ["RNA-Seq genes to pathways"]({% link topics/transcriptomics/tutorials/rna-seq-genes-to-pathways/tutorial.md %}).
{: .comment}

## Análisis de rutas KEGG

**goseq** también puede utilizarse para identificar vías KEGG interesantes. La base de datos de vías KEGG es una colección de mapas de vías que representan el conocimiento actual de las redes de interacciones, reacciones y relaciones moleculares. Un mapa puede integrar muchas entidades, incluidos genes, proteínas, ARN, compuestos químicos, glicanos y reacciones químicas, así como genes de enfermedades y dianas farmacológicas.

Por ejemplo, la ruta `dme00010` representa el proceso de glucólisis (conversión de glucosa en piruvato con generación de pequeñas cantidades de ATP y NADH) para Drosophila melanogaster:

![dme00010 KEGG pathway](../../images/ref-based/dme00010_empty.png)

> <hands-on-title>Realizar análisis de vías KEGG</hands-on-title>
> 
> 1. {% tool [goseq](toolshed.g2.bx.psu.edu/repos/iuc/goseq/goseq/1.50.0+galaxy0) %} with
>    - *"Differentially expressed genes file "*: `Gene IDs and differential expression`
>    - *"Gene lengths file "*: `Gene IDs and length`
>    - *"Gene categories "*: `Get categories`
>       - *"Select a genome to use "*: `Fruit fly (dm6)`
>       - *"Select Gene ID format "*: `Ensembl Gene ID`
>       - *"Select one or more categories "*: `KEGG`
>    - En *"Output options "*
>      - *"Output Top GO terms plot? "*: `No`
>      - *"Extract the DE genes for the categories (GO/KEGG terms)? "*: `Yes`
> 
{: .hands_on}

**goseq** genera con estos parámetros 2 salidas:

1. Una tabla grande con los términos KEGG y algunas estadísticas

   > <question-title></question-title>
   > 
   > 1. ¿Cuántos términos de vías KEGG se han identificado?
   > 2. ¿Cuántos términos de vías KEGG están sobrerrepresentados con un valor P ajustado < 0,05?
   > 3. ¿Cuáles son los términos de vías KEGG sobrerrepresentados?
   > 4. ¿Cuántos términos de vías KEGG están infrarrepresentados con un valor P ajustado < 0,05?
   > 
   > > <solution-title></solution-title>
   > > 
   > > 1. El archivo tiene 128 líneas, incluida la cabecera, por lo que se han identificado 127 vías KEGG.
   > > 2. 2 vías KEGG (2,34%) están sobrerrepresentadas, utilizando {% tool [Filtrar datos en cualquier columna utilizando expresiones simples](Filter1) %} en c6 (valor p ajustado para vías KEGG sobrerrepresentadas)
   > > 3. Las dos vías KEGG sobrerrepresentadas son `01100` y `00010`. Buscándolas en la [base de datos KEGG](https://www.genome.jp/kegg/kegg2.html), podemos encontrar más información sobre estas vías: `01100` corresponde a todas las vías metabólicas y `00010` a la vía de la Glucólisis / Gluconeogénesis.
   > > 4. Ninguna vía KEGG está infrarrepresentada, utilizando {% tool [Filtrar datos en cualquier columna utilizando expresiones simples](Filter1) %} en c7 (valor p ajustado para vías KEGG infrarrepresentadas)
   > {: .solution}
   {: .question}

2. Una tabla con los genes expresados diferencialmente (de la lista que proporcionamos) asociados a las vías KEGG (`DE genes for categories (GO/KEGG terms)`)

Podríamos investigar qué genes están implicados en qué rutas consultando el segundo archivo generado por **goseq**. Sin embargo, esto puede ser engorroso y nos gustaría ver las rutas como se representa en la imagen anterior. **Pathview** ({% cite luo2013pathview %}) puede ayudar a generar automáticamente imágenes similares a la anterior, a la vez que añade información extra sobre los genes (por ejemplo, expresión) en nuestro estudio.

Esta herramienta necesita 2 entradas principales:

- ID de la(s) vía(s) a trazar, ya sea como un único ID o como un archivo con una columna con los ID de la vía
- Un archivo tabular con los genes del experimento RNA-Seq con 2 (o más) columnas:
  - los ID de los genes (únicos dentro del archivo)
  - información sobre los genes

    Puede ser, por ejemplo, un valor p o un cambio de pliegue. Esta información se añadirá al trazado de la ruta: el nodo del gen correspondiente se coloreará con el valor. Si hay diferentes columnas, las diferentes informaciones se trazarán una al lado de la otra en el nodo.

Aquí nos gustaría visualizar las 2 vías KEGG: la sobrerrepresentada `00010` (Glucólisis / Gluconeogénesis) y la más infrarrepresentada (aunque no significativamente) `03040` (Spliceosoma). Nos gustaría que los nodos de genes fueran coloreados por Log2 Fold Change para los genes diferencialmente expresados debido al tratamiento.

> <hands-on-title>Overlay log2FC on KEGG pathway</hands-on-title> (Superponer log2FC en la ruta KEGG)
> 
> 1. {% tool [Cut](Cut1) %} columnas de una tabla con los siguientes parámetros:
>    - *"Cut columns "*: `c1,c3`
>    - *"Delimited by"*: `Tab`
>    - {% icon param-file %} *"De "*: `Genes with significant adj p-value`
> 
> 2. Cambie el nombre a `Genes with significant adj p-value and their Log2 FC`
> 
>    Extraemos el ID y el Log2 Fold Change de los genes que tienen un valor p ajustado significativo.
> 
> 3. Cree un nuevo archivo tabular a partir de lo siguiente (IDs de las rutas a trazar) llamado `KEGG pathways to plot`
> 
>    ```text
>    00010
>    03040
>    ```
> 
> 4. {% tool [Pathview](toolshed.g2.bx.psu.edu/repos/iuc/pathview/pathview/1.34.0+galaxy0) %} with
>    - *"Number of pathways to plot "*: `Multiple`
>      - {% icon param-file %} *"KEGG pathways "*: `KEGG pathways to plot`
>      - *"Does the file have header (a first line with column names)? "*: `No`
>    - *"Especies a utilizar "*: `Fly`
>    - *"Provide a gene data file?*: `Yes`
>      - {% icon param-file %} *"Gene data "*: `Genes with significant adj p-value and their Log2 FC`
>      - *"Does the file have header (a first line with column names)? "*: `Yes`
>      - *"Format for gene data "*: `Ensembl Gene ID`
>    - *"Provide a compound data file? "*: `No`
>    - En *"Output options "*
>      - *"Output for pathway "*: `KEGG native`
>        - *"Plot on same layer?"*: `Yes`
> 
{: .hands_on}

**Pathview** genera una colección con la visualización KEGG: un archivo por ruta.

> <question-title></question-title>
> 
> `dme00010` Ruta KEGG de **Pathview**
> 
> ![KEGG pathway](../../images/ref-based/dme00010.png)
> 
> 1. ¿Qué son las cajas de color?
> 2. ¿Cuál es el código de color?
> 
> > <solution-title></solution-title>
> > 
> > 1. Los recuadros de color son genes de la vía que se expresan de forma diferencial
> > 2. Preste atención a que el código de colores es contrario a la intuición: el verde es para valores inferiores a 0, es decir, para genes con un log2FC < 0 y el rojo para genes con un log2FC > 0.
> {: .solution}
{: .question}

{% comment %}

# Inferencia del uso diferencial de exones

A continuación, nos gustaría conocer el uso diferencial de exones entre las muestras tratadas (PS depleted) y las no tratadas utilizando los recuentos de exones RNA-Seq. Reelaboraremos los resultados del mapeo que generamos anteriormente.

Utilizaremos [DEXSeq](https://www.bioconductor.org/packages/release/bioc/html/DEXSeq.html). DEXSeq detecta genes de alta sensibilidad, y en muchos casos exones, que están sujetos a un uso diferencial de exones. Pero primero, como para la expresión génica diferencial, necesitamos contar el número de lecturas que mapean a los exones.

## Contar el número de lecturas por exón

Este paso es similar al paso de [contar el número de lecturas por gen anotado](#count-the-number-of-reads-per-annotated-gene) excepto que, en lugar de HTSeq-count, estamos usando DEXSeq-Count.

> <hands-on-title>Cuento del número de lecturas por exón</hands-on-title>
> 
> 1. {% tool [DEXSeq-Count](toolshed.g2.bx.psu.edu/repos/iuc/dexseq/dexseq_count/1.28.1.0) %}: Utilice el **DEXSeq-Count** para preparar las anotaciones *Drosophila* para extraer sólo exones con ids de genes correspondientes
>     - *"Mode of operation "*: `Prepare annotation`
>       - {% icon param-file %} *"GTF file "*: `Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz`
> 
>    El resultado es de nuevo un archivo GTF listo para ser utilizado para el recuento
> 
> 2. {% tool [DEXSeq-Count](toolshed.g2.bx.psu.edu/repos/iuc/dexseq/dexseq_count/1.28.1.0) %}: Contar lecturas usando **DEXSeq-Count** con
>     - *"Mode of operation "*: `Count reads`
>       - {% icon param-files %} *"Input bam file"*: los archivos `BAM` generados por **RNA STAR**
>       - {% icon param-file %} *"DEXSeq compatible GTF file "*: el archivo GTF generado por **DEXSeq-Count**
>       - *"Is library paired end? "*: `Yes`
>       - *"Is library strand specific?*: `No`
>       - *"Skip all reads with alignment quality lower than the given minimum value "*: `10`
> 
{: .hands_on}

DEXSeq genera una tabla de recuento similar a la generada por featureCounts, pero con recuentos para los exones.

> <question-title></question-title>
> 
> 1. ¿Qué exón tiene más lecturas asignadas en ambas muestras?
> 2. ¿A qué gen pertenece este exón?
> 3. ¿Hay alguna relación con el resultado anterior obtenido con featureCounts?
> 
> > <solution-title></solution-title>
> > 
> > FBgn0284245:005 es el exón con más lecturas asignadas en ambas muestras. Forma parte de FBgn0284245, el rasgo con más lecturas mapeadas en él (a partir de featureCounts).
> > 
> {: .solution}
{: .question}

## Uso diferencial de exones

El uso de DEXSeq es similar al de DESeq2. Utiliza estadísticas similares para encontrar exones utilizados de forma diferencial.

Al igual que para DESeq2, en el paso anterior, sólo contamos las lecturas que correspondían a exones del cromosoma 4 y para una sola muestra. Para poder identificar el uso diferencial de exones inducido por la depleción de PS, todos los conjuntos de datos (3 tratados y 4 no tratados) deben analizarse siguiendo el mismo procedimiento. Para ahorrar tiempo, lo hemos hecho por usted. Los resultados están disponibles en [Zenodo]({{ page.zenodo_link }}):

- Resultados de la ejecución de DEXSeq-count en modo "Preparar anotación
- Siete archivos de recuento generados en modo 'Count reads

> <hands-on-title></hands-on-title>
> 
> 1. Crear un **nuevo historial vacío**
> 
>    {% snippet faqs/galaxy-es/histories_create_new.md %}
> 
> 2. Importar los siete archivos de recuento de [Zenodo]({{ page.zenodo_link }}) o la biblioteca de Datos Compartidos (si está disponible):
> 
>    - `Drosophila_melanogaster.BDGP6.87.dexseq.gtf`
>    - `GSM461176_untreat_single.exon.counts`
>    - `GSM461177_untreat_paired.exon.counts`
>    - `GSM461178_untreat_paired.exon.counts`
>    - `GSM461179_treat_single.exon.counts`
>    - `GSM461180_treat_paired.exon.counts`
>    - `GSM461181_treat_paired.exon.counts`
>    - `GSM461182_untreat_single.exon.counts`
> 
>    ```text
>    {{ page.zenodo_link }}/files/Drosophila_melanogaster.BDGP6.87.dexseq.gtf
>    {{ page.zenodo_link }}/files/GSM461176_untreat_single.exon.counts
>    {{ page.zenodo_link }}/files/GSM461177_untreat_paired.exon.counts
>    {{ page.zenodo_link }}/files/GSM461178_untreat_paired.exon.counts
>    {{ page.zenodo_link }}/files/GSM461179_treat_single.exon.counts
>    {{ page.zenodo_link }}/files/GSM461180_treat_paired.exon.counts
>    {{ page.zenodo_link }}/files/GSM461181_treat_paired.exon.counts
>    {{ page.zenodo_link }}/files/GSM461182_untreat_single.exon.counts
>    ```
> 
> 3. {% tool [DEXSeq](toolshed.g2.bx.psu.edu/repos/iuc/dexseq/dexseq/1.28.1+galaxy1) %}: Ejecutar **DEXSeq** con
>    - {% icon param-file %} *"Archivo GTF creado a partir de la herramienta DEXSeq-Count "*: `Drosophila_melanogaster.BDGP6.87.dexseq.gtf`
>    - En *"Factor "*:
>       - En "1: Factor"
>           - *"Specify a factor level "*: `condition`
>           - En *"Factor level "*:
>               - En *"1: Factor level "*:
>                   - *"Specify a Factor level "*: `treated`
>                   - {% icon param-files %} *"Counts file(s) "*: los 3 ficheros de recuento de exones con `treated` en su nombre
>               - En *"2: Factor level "*:
>                   - *"Specify a Factor level "*: `untreated`
>                   - {% icon param-files %} *"Counts file(s) "*: los 4 ficheros de recuento de exones con `untreated` en su nombre
>       - Haga clic en *"Insert Factor "* (no en "Insert Factor level")
>       - En "2: Factor"
>           - "Specify a factor level" a `sequencing`
>           - En *"Factor level "*:
>               - En *"1: Factor level "*:
>                   - *"Specify a Factor level "*: `PE`
>                   - {% icon param-files %} *"Counts file(s) "*: los 4 ficheros de recuento de exones con `paired` en su nombre
>               - En *"2: Factor level "*:
>                   - *"Specify a Factor level "*: `SE`
>                   - {% icon param-files %} *"Counts file(s) "*: los 3 ficheros de recuento de exones con `single` en su nombre
> 
>    > <comment-title></comment-title>
>    > 
>    > A diferencia de DESeq2, DEXSeq no permite nombres de factores primarios flexibles. Utilice siempre el nombre del factor primario como "condición"
>    > 
>    {: .comment}
{: .hands_on}

Al igual que DESeq2, DEXSeq genera una tabla con:

1. Identificadores de exón
2. Identificadores de genes
3. Identificadores de exón en el gen
4. Recuentos medios normalizados, promediados en todas las muestras de ambas condiciones
5. Logaritmo (en base 2) del cambio de pliegue

   Los cambios de pliegue log2 se basan en el Factor level primario 1 frente al Factor level 2. El orden de los niveles del factor es entonces importante. Por ejemplo, para el factor "Condición", DESeq2 calcula los cambios de pliegues de las muestras "tratadas" frente a las "no tratadas", *es decir* los valores corresponden a la regulación al alza o a la baja de los genes en las muestras tratadas.

6. Estimación del error estándar para la estimación del cambio de pliegue log2
7. Valor *p* para la significación estadística de este cambio
8. Valor *p* ajustado para pruebas múltiples con el procedimiento Benjamini-Hochberg que controla la tasa de falsos descubrimientos ([FDR](https://en.wikipedia.org/wiki/False_discovery_rate))

> <hands-on-title></hands-on-title>
> 
> 1. {% tool [Filter data on any column using simple expressions](Filter1) %} para extraer exones con un uso diferencial significativo (valor *p* ajustado igual o inferior a 0,05) entre muestras tratadas y no tratadas
> 
> > <question-title></question-title>
> > 
> > ¿Cuántos exones muestran un cambio de uso significativo entre estas condiciones?
> > 
> > > <solution-title></solution-title>
> > > 
> > > Obtenemos 38 exones (12,38%) con un cambio de uso significativo entre las muestras tratadas y no tratadas.
> > > 
> > {: .solution}
> {: .question}
{: .hands_on}

{% endcomment %}

# Conclusión


En este tutorial, hemos analizado datos reales de secuenciación de ARN para extraer información útil, como qué genes están regulados al alza o a la baja por la depleción del gen *Pasilla*, pero también en qué términos GO o vías KEGG están implicados. Para responder a estas preguntas, analizamos conjuntos de datos de secuencias de ARN utilizando un enfoque de análisis de datos de ARN-Seq basado en referencias. Este enfoque puede resumirse con el siguiente esquema:

![Resumen del proceso de análisis utilizado](../../images/ref-based/tutorial-scheme.png "Resumen del proceso de análisis utilizado")

