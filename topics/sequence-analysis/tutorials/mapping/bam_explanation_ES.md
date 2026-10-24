Un fichero BAM ([Binary Alignment Map](https://en.wikipedia.org/wiki/SAM_(file_format))) es un archivo binario comprimido que almacena las secuencias leídas, indicando si han sido alineadas a una secuencia de referencia (por ejemplo, un cromosoma) y, en caso afirmativo, la posición en la secuencia de referencia donde han sido alineadas.

> <hands-on-title>Inspeccionar un BAM/SAM file</hands-on-title>
>
> 1. Inspecciona el {% icon param-file %} output de **{{ include.mapper }}** {% icon tool %}
>
{: .hands_on}

Un fichero BAM (o SAM, la versión sin comprimir) consiste de:

- Una sección de cabecera (las líneas que comienzan con `@`) que contiene metadatos, en particular los nombres y longitudes de los cromosomas (líneas que comienzan con el símbolo `@SQ`).
- Una sección de alineamiento que consiste en una tabla con 11 campos obligatorios, además de un número variable de campos opcionales:

    Col | Campo | Tipo | Breve descripción
    --- | --- | --- | ---
    1 | QNAME | String | Query template NAME
    2 | FLAG | Integer | Bitwise FLAG
    3 | RNAME | String | References sequence NAME
    4 | POS | Integer | 1- based leftmost mapping POSition
    5 | MAPQ | Integer | MAPping Quality
    6 | CIGAR | String | CIGAR String
    7 | RNEXT | String | Ref. name of the mate/next read
    8 | PNEXT | Integer | Position of the mate/next read
    9 | TLEN | Integer | Observed Template LENgth
    10 | SEQ | String | Segment SEQuence
    11 | QUAL | String | ASCII of Phred-scaled base QUALity+33 

> <question-title></question-title>
>
> 1. ¿Qué información se encuentra en un fichero BAM/SAM?
> 2. ¿Cuál es la información adicional en comparación a un archivo FASTQ?
>
> > <solution-title></solution-title>
> > 1. Información de calidad y de secuencias, como un FASTQ
> > 2. Información de mapeo, localización de la lectura en el cromosoma, calidad del mapeo, etc.
> {: .solution }
{: .question}
