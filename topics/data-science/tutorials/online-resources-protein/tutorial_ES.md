---
layout: tutorial_hands_on
title: Una proteína a lo largo de la página UniProt
level: Introductory
zenodo_link: ''
questions:
- ¿Cómo se pueden buscar proteínas utilizando texto, nombres de genes o nombres de proteínas?
- ¿Cómo se interpreta la información en la parte superior de la página de entrada de UniProt?
- ¿Qué tipo de información se puede esperar de los distintos formatos de descarga, como FASTA y JSON?
- ¿Cómo se describe la función de una proteína como las opsinas en la sección "Function"?
- ¿Qué información estructurada se encuentra en las secciones "Names and Taxonomy", "Subcellular location", "Disease & Variants", "PTM/Processing"?
- ¿Cómo informarse sobre la expresión, interacciones, estructura, familia, secuencia y proteínas similares?
- ¿Cómo ayudan las pestañas "Variant viewer" y "Feature viewer" a mapear la información de la proteína a lo largo de su secuencia?
- ¿Qué muestra la pestaña "Publications" y cómo se pueden filtrar las publicaciones?
- ¿Cuál es la importancia de hacer un seguimiento de los cambios en la anotación de la entrada a lo largo del tiempo?
objectives:
- Explorando las entradas de proteínas en UniProtKB, interpreta la función de la proteína, su taxonomía, estructura, interacciones, variantes y más.
- Utiliza identificadores únicos para conectar bases de datos, descargar datos de genes y proteínas, visualizar y comparar características de secuencias.
time_estimation: 1H
key_points:
- Cómo navegar por las entradas de UniProtKB, accediendo a detalles exhaustivos sobre proteínas, como su función, taxonomía e interacciones.
- El Variant y el Feature viewer son herramientas útiles para explorar visualmente variantes, dominios, modificaciones y otras características clave de la secuencia.
- Amplía tu comprensión utilizando enlaces externos para cruzar referencias de datos y descubrir relaciones complejas.
- Explora la pestaña Historia para acceder a versiones anteriores de las anotaciones de la entrada.
requirements:
- type: internal
  topic_name: data-science
  tutorials:
  - online-resources-gene
contributions:
  authorship:
  - lisanna
  - bebatut
  translation:
  - hvelab
  - unode
  - ocaisa
  funding:
  - biont
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


Al realizar un análisis de datos biológicos, es posible que nos encontremos con algunas proteínas interesantes, que necesitemos explorar estos genes. ¿Pero cómo podemos hacerlo? ¿Cuáles son los recursos disponibles para ello? ¿Y cómo navegar por ellos?

El objetivo de este tutorial es familiarizarnos con ello, utilizando como ejemplo las opsinas humanas.

> <comment-title></comment-title>
>
> Este tutorial es un poco atípico: no trabajaremos en Galaxy sino principalmente fuera de ella, en las páginas de la base de datos [UniProt](https://uniprot.org).
{: .comment}

> <comment-title></comment-title>
>
> Este tutorial está diseñado para ser la continuación del tutorial ["Un gen a través de formatos de archivo"]({% link topics/data-science/tutorials/online-resources-gene/tutorial.md %}), pero también puede ser consultado como un módulo independiente.
{: .comment}

Las opsinas se encuentran en las células de la retina. Captan la luz e inician la secuencia de señales que dan lugar a la visión, y esa es la razón por la que, cuando se ven comprometidas, se asocian al daltonismo y a otras deficiencias visuales.

> <comment-title>Fuentes de información de este tutorial</comment-title>
> 
> El tutorial que usted está consultando fue desarrollado principalmente consultando los recursos de UniProtKB, en particular el tutorial [Explore UniProtKB entry](https://www.uniprot.org/help/explore_uniprotkb_entry). Algunas frases se reportan desde allí sin modificaciones.
> 
> Además, el tema se eligió basándose en el [Tutorial de Bioinformática] de Gale Rhodes (https://spdbv.unil.ch/TheMolecularLevel/Matics/index.html). Aunque el tutorial ya no se puede seguir paso a paso debido a cómo los recursos mencionados cambiaron con el tiempo, podría proporcionar ideas adicionales sobre las opsinas y, en particular, sobre cómo se pueden construir modelos estructurales de proteínas basados en información evolutiva.
{: .comment}


> <agenda-title></agenda-title>
> 
> En este tutorial trataremos:
> 
> 1. TOC
> {:toc}
{: .agenda}

# La página de entrada de UniProtKB

El portal que hay que visitar para obtener toda la información sobre una proteína es [UniProtKB](https://www.uniprot.org/). Podemos buscar en él utilizando una búsqueda de texto, o el nombre del gen o de la proteína. Probemos primero con un conjunto de palabras clave genéricas, como `Human opsin`.

> <hands-on-title>Buscar opsina humana en UniProtKB</hands-on-title>
> 
> 1. Abrir el [UniProtKB](https://www.uniprot.org/)
> 2. Escriba `Human opsin` en la barra de búsqueda
> 3. Lanza la búsqueda
> 
{: .hands-on}

> <question-title></question-title>
> 
> ¿Cuántos resultados obtuvimos?
> 
> > <solution-title></solution-title>
> > 
> > 410 resultados (en el momento de la preparación de este tutorial)
> > 
> {: .solution}
{: .question}

Estos 410 resultados nos dan la sensación de que tenemos que ser más específicos (aunque - spoiler - nuestro objetivo real está entre los primeros resultados).

Para ser suficientemente específicos, sugerimos utilizar un identificador único. Del [tutorial anterior]({% link topics/data-science/tutorials/online-resources-gene/tutorial.md %}) conocemos el nombre del gen de la proteína que estamos buscando, `OPN1LW`.

> <hands-on-title>Buscar OPN1LW en UniProtKB</hands-on-title>
> 
> 1. Escriba `OPN1LW` en la barra de búsqueda superior
> 2. Lanza la búsqueda
> 
{: .hands-on}

> <question-title></question-title>
> 
> 1. ¿Cuántos resultados obtuvimos?
> 2. ¿Qué debemos hacer para reducir este número?
> 
> > <solution-title></solution-title>
> > 
> > 1. 200+ resultados (en el momento de la preparación de este tutorial)
> > 2. Necesitamos aclarar lo que estamos buscando: OPN1LW humano
> > 
> {: .solution}
{: .question}

Necesitamos añadir `Human` para aclarar lo que estamos buscando.

> <hands-on-title>Búsqueda de OPN1LW humano en UniProtKB</hands-on-title>
> 
> 1. Escriba `Human OPN1LW` en la barra de búsqueda superior
> 2. Lanza la búsqueda
> 
{: .hands-on}

> <question-title></question-title>
> 
> 1. ¿Cuántos resultados obtuvimos?
> 2. ¿Tenemos un resultado que incluya OPN1LW como nombre de gen?
> 
> > <solution-title></solution-title>
> > 
> > 1. 7 resultados (en el momento de la preparación de este tutorial)
> > 2. El primer resultado está etiquetado con `Gene: OPN1LW (RCP)`
> > 
> {: .solution}
{: .question}

El primer resultado, etiquetado con `Gene: OPN1LW (RCP)`, es nuestro objetivo, `P04000 · OPSR_HUMAN`. Antes de abrir la página, hay que tener en cuenta dos cosas:

1. El nombre de la proteína `OPSR_HUMAN` es diferente del nombre del gen, así como sus IDs son.
2. Esta entrada tiene una estrella dorada, lo que significa que fue anotada y curada manualmente.

# Inspeccionar una entrada UniProt

> <hands-on-title>Abrir un resultado en UniProt</hands-on-title>
> 
> 1. Haga clic en `P04000 · OPSR_HUMAN`
> 
{: .hands-on}

![Captura de pantalla de la cabecera de la página de entrada de UniProt](./images/UniProt.png "UniProt page")

Para navegar por esta larga página, el menú (barra de navegación) de la izquierda será de gran utilidad. Sólo por él, entendemos que esta base de datos contiene información sobre la entrada en:
- las funciones conocidas,
- la taxonomía,
- la ubicación,
- variantes y enfermedades asociadas,
- Modificación postraduccional (PTM),
- la expresión,
- las interacciones,
- la estructura,
- los dominios y su clasificación,
- las secuencias
- proteínas similares.

La barra de navegación permanece en el mismo lugar de la pantalla a medida que se avanza y retrocede en una entrada, de modo que se puede navegar rápidamente a las secciones de interés. Consultaremos todas las secciones mencionadas por separado, pero centrémonos primero en las cabeceras de la izquierda.

En la parte superior de la página, puede ver la entrada UniProt y su nombre, el nombre de la proteína y del gen, el organismo, si la entrada de la proteína ha sido revisada manualmente por un curador de UniProt, su puntuación de anotación y el nivel de evidencia de su existencia.

Debajo de la cabecera principal, encontrará una serie de pestañas (*Entrada*, *Visor de variantes*, *Visor de características*, *Publicaciones*, *Enlaces externos*, *Historia*). Las pestañas permiten cambiar entre la entrada, una vista gráfica de las características de la secuencia (Visor de características), las publicaciones y los enlaces externos, pero por el momento las ignora y no se mueve de la pestaña *Entrada*.

## Entrada

El siguiente menú ya forma parte de la pestaña *Entrada*. Nos permite ejecutar una búsqueda de similitud de secuencias BLAST en la entrada, alinearla con todas sus isoformas, descargar la entrada en varios formatos, o añadirla a la cesta para guardarla para más tarde.

> <question-title></question-title>
> 
> 1. ¿Cuáles son los formatos disponibles en el menú desplegable *Descargar*?
> 2. ¿Qué tipo de información descargaríamos a través de estos formatos de archivo?
> 
> > <solution-title></solution-title>
> > 
> > 1. Los formatos son: `Text`, `FASTA (canonical)`, `FASTA (canonical & isoform`, `JSON`, `XML`, `RDF/XML`, `GFF`
> > 2. Los formatos`FASTA` deberían sonar familiares (después del tutorial preliminar), e incluyen la secuencia de la proteína, eventualmente con sus isoformas (en cuyo caso será un multi-FASTA). Aparte de éstos, todos los demás formatos no son específicos de una proteína, ni siquiera de una biología. Se trata de formatos de archivo generales muy utilizados por los sitios web para incluir la información que figura en la página. Por lo tanto, descargando el fichero `text` (o incluso mejor el `json`), descargaríamos la misma anotación a la que accedemos en esta página, pero en un formato más fácil de parsear programáticamente.
> > 
> {: .solution}
{: .question}

Desplacémonos ahora por la página de entrada, sección por sección.

### Función

Esta sección resume las funciones de esta proteína como sigue:

Los pigmentos visuales son las moléculas que absorben la luz y se encargan de la visión. Están formados por una apoproteína, la opsina, unida covalentemente al cis-retinal.

Independientemente del nivel de detalle que entiendas (dependiendo de tu formación), esto es impresionantemente corto y específico teniendo en cuenta la enorme cantidad de literatura y estudios que existen más allá de la determinación de la función de una proteína. De todos modos, alguien hizo el trabajo por nosotros, y esta proteína ya está completamente clasificada en la Ontología Genética (GO), que describe la función molecular, el proceso biológico y el componente celular de cualquier proteína clasificada.

GO es un ejemplo perfecto de base de datos / recurso que parte de un universo de conocimientos muy complejo y lo traduce a un gráfico más simple, a riesgo de perder detalles. Esto tiene la gran ventaja de organizar la información, hacerla contable y analizable y accesible mediante programación, lo que en última instancia nos permite disponer de estas largas páginas de resumen y Bases de Conocimiento.

> <question-title></question-title>
> 
> 1. ¿A qué funciones moleculares está anotada esta proteína?
> 2. ¿A qué componentes celulares está anotada esta proteína?
> 3. ¿A qué procesos biológicos está anotada esta proteína?
> 
> > <solution-title></solution-title>
> > 
> > 1. Proteína fotorreceptora, receptor acoplado a proteína G
> > 2. Membrana del disco fotorreceptor
> > 2. Transducción sensorial, Visión
> {: .solution}
{: .question}

### Nombres y taxonomía

Otros ejemplos de información estructurada están disponibles en la siguiente sección, por ejemplo, en la taxonomía. Esta sección también informa de otros identificadores únicos que se refieren a la misma entidad biológica o a entidades vinculadas (por ejemplo, enfermedades asociadas en el menú `MIM`).

> <question-title></question-title>
> 
> 1. ¿Cuál es el identificador taxonómico asociado a esta proteína?
> 2. ¿Cuál es el identificador del proteoma asociado a esta proteína?
> 
> > <solution-title></solution-title>
> > 
> > 1. 9606, es decir, Homo sapiens
> > 2. UP000005640, componente del cromosoma Xs
> {: .solution}
{: .question}

## Localización subcelular

Ya sabemos dónde está nuestra proteína en el cuerpo humano (en la retina, como se especifica en el resumen de la función), pero ¿dónde está en la célula?

> <question-title></question-title>
> 
> 1. ¿Dónde se encuentra nuestra proteína en la célula?
> 2. ¿Es coherente con la anotación GO observada anteriormente?
> 
> > <solution-title></solution-title>
> > 
> > 1. La sección explica que se trata de una "proteína de membrana multipasos", lo que significa que es una proteína que se inserta en la membrana celular y la atraviesa varias veces.
> > 2. La anotación GO de la parte superior menciona que nos referimos a la membrana (celular) del fotorreceptor en particular.
> > 
> {: .solution}
> 
{: .question}

La sección Localización subcelular incluye un área de *Características* que detalla qué secciones, a lo largo de la secuencia de la proteína, están insertadas en la membrana (Transmembrana) y cuáles no (Dominio topológico).

> <question-title></question-title>
> 
> ¿Cuántos dominios transmembrana y dominios topológicos existen?
> 
> > <solution-title></solution-title>
> > 
> > 8 dominios transmembrana y 7 topológicos
> > 
> {: .solution}
> 
{: .question}

### Enfermedad & Variantes

Como sabemos por el tutorial anterior, este gen/proteína está asociado a múltiples enfermedades. Esta sección detalla esta asociación y enumera las variantes específicas que se han detectado como relacionadas con la enfermedad.

> <question-title></question-title>
> 
> ¿Qué tipos de estudios científicos permiten evaluar la asociación de una variante genética a enfermedades?
> 
> > <solution-title></solution-title>
> > 
> > Tres métodos comúnmente utilizados para evaluar la asociación de una variante genética con una enfermedad son:
> > 
> > - Estudios de asociación del genoma completo (GWAS)
> > 
> >   Los GWAS se utilizan ampliamente para identificar variantes genéticas comunes asociadas a enfermedades. Consisten en escanear el genoma completo de un gran número de individuos para identificar variaciones vinculadas a una enfermedad o rasgo concreto.
> > 
> > - Estudios Caso-Control
> > 
> >   Los estudios de casos y controles se emplean con frecuencia para comparar individuos con una enfermedad con aquellos que no la padecen, centrándose en la presencia o frecuencia de variantes genéticas específicas en ambos grupos.
> > 
> > - Estudios Familiares
> > 
> >   Los estudios basados en familias implican el análisis de variantes genéticas dentro de familias en las que varios miembros están afectados por una enfermedad. Mediante el estudio de los patrones de herencia de las variantes genéticas y su asociación con la enfermedad dentro de las familias, los investigadores pueden identificar posibles genes asociados a la enfermedad.
> > 
> > Este tipo de estudios implicaría un uso extensivo de los tipos de archivos para gestionar datos genómicos, tales como: SAM (Sequence Alignment Map), BAM (Binary Alignment Map), VCF (Variant Calling Format) etc.
> > 
> {: .solution}
{: .question}

Esta sección también incluye un área *Features*, donde se mapean las variantes naturales a lo largo de la secuencia. Abajo, también se destaca que una vista más detallada de las características a lo largo de la secuencia se proporciona en la pestaña *Enfermedad & Variantes*, pero no la abramos por ahora.

### PTM/Procesamiento

Una modificación post-traduccional (PTM) es un evento de procesamiento covalente resultante de una escisión proteolítica o de la adición de un grupo modificador a un aminoácido.

> <question-title></question-title>
> 
> ¿Cuáles son las modificaciones postraduccionales de nuestra proteína?
> 
> > <solution-title></solution-title>
> > 
> > Cadena, glicosilación, enlace disulfuro, residuo modificado
> > 
> {: .solution}
> 
{: .question}

### Expresión

Ya sabemos dónde se encuentra la proteína en la célula, pero en el caso de las proteínas humanas a menudo disponemos de información sobre dónde se encuentra en el cuerpo humano, es decir, en qué tejidos. Esta información puede proceder del Human [ExpressionAtlas](https://www.ebi.ac.uk/gxa/home) o de otros recursos similares.

> <question-title></question-title>
> 
> ¿En qué tejido se encuentra la proteína?
> 
> > <solution-title></solution-title>
> > 
> > Los tres pigmentos colorantes se encuentran en las células fotorreceptoras de los conos.
> > 
> {: .solution}
> 
{: .question}

### Interacción

Las proteínas desempeñan su función a través de su interacción con el entorno, en particular con otras proteínas. Esta sección informa de los interactores de nuestra proteína de interés, en una tabla que también podemos filtrar por localización subcelular, enfermedades y tipo de interacción.

La fuente de esta información son bases de datos como STRING, y la página de entrada de nuestra proteína está directamente enlazada desde esta sección.

> <hands-on-title>Búsqueda de OPN1LW humano en UniProtKB</hands-on-title>
> 
> 1. Haga clic en el enlace [STRING](https://string-db.org/network/9606.ENSP00000358967) en una pestaña diferente
> 
{: .hands-on}

> <question-title></question-title>
> 
> 1. ¿Cuántos formatos de archivo diferentes se pueden descargar desde allí?
> 2. ¿Qué tipo de información se transmitirá en cada fichero?
> 
> > <solution-title></solution-title>
> > 
> > STRING proporciona datos en formatos de archivo descargables para apoyar análisis posteriores. El principal formato de archivo utilizado por STRING es el formato "TSV" (Tab-Separated Values), que presenta los datos de interacción de proteínas en un formato tabular estructurado. Este formato es idóneo para integrarlo fácilmente en diversas herramientas y programas de análisis de datos. Además, STRING ofrece datos en formato XML PSI-MI (Proteomics Standards Initiative Molecular Interactions), un estándar para representar datos de interacciones proteínicas que permite la compatibilidad con otras bases de datos de interacciones y plataformas de análisis. Estos formatos de archivo permiten a los investigadores aprovechar la gran cantidad de información sobre interacciones proteicas que contiene STRING para sus propios estudios y análisis. Los investigadores también pueden descargar representaciones visuales de las redes de proteínas en formatos de imagen como PNG y SVG, adecuados para presentaciones y publicaciones. Para análisis avanzados, STRING ofrece "archivos planos" que contienen información detallada sobre las interacciones, y archivos "MFA" (Multiple Alignment Format), útiles para comparar múltiples secuencias de proteínas. Estos diversos formatos de archivo descargables permiten a los investigadores aprovechar la riqueza de la información sobre interacciones proteicas de STRING para sus propios estudios y análisis.
> > 
> > 
> {: .solution}
> 
{: .question}

### Estructura

¿Siente curiosidad por las intrincadas estructuras tridimensionales de las proteínas? La sección *Estructura* de la página de entrada de UniProtKB es su puerta de entrada para explorar el fascinante mundo de la arquitectura de las proteínas.

En esta sección encontrará información sobre estructuras de proteínas determinadas experimentalmente. Estas estructuras proporcionan información crucial sobre el funcionamiento de las proteínas y su interacción con otras moléculas. Descubrirá vistas interactivas de la estructura de la proteína que puede explorar directamente dentro de la entrada UniProtKB. Esta característica proporciona una forma atractiva de navegar por los dominios de la proteína, los sitios de unión y otras regiones funcionales. Profundizando en la sección *Estructura*, comprenderá mejor la base física de la función de las proteínas y descubrirá la riqueza de información que pueden revelar los datos estructurales.

> <question-title></question-title>
> 
> 1. ¿Cuál es la variante asociada al daltonismo?
> 2. ¿Puedes encontrar ese aminoácido específico en la estructura?
> 3. ¿Puede formular una conjetura de por qué esta mutación es disruptiva?
> 
> > <solution-title></solution-title>
> > 
> > 1. En la sección *Enfermedad y Variantes*, descubrimos que el cambio de Glicina (G) a Ácido Glutámico (E) en la posición 338 de la secuencia de la proteína está asociado al Daltonismo.
> > 2. En el visor de estructura, podemos mover la molécula y pasar el ratón sobre la estructura para encontrar el AA en la posición 338. Puede llevar algún tiempo seguir las múltiples disposiciones helicoidales de estas estructuras. La glicina en 338 no está en una hélice, sino en lo que parece un bucle justo antes de una zona de baja confianza en la estructura.
> > 3. Basándonos en la información que hemos recopilado hasta ahora, podríamos formular una hipótesis de por qué es disruptiva. No está en una hélice (normalmente, en las proteínas transmembrana, las hélices se insertan en la membrana), por lo tanto, está en uno de los dominios más grandes que sobresalen de la membrana, dentro o fuera de la célula. Esta mutación probablemente no interrumpe la estructura en sus segmentos intramembrana, sino en uno de los dominios funcionales. Si quieres profundizar, puedes comprobar si se trata del segmento extra o intracelular en el **Visor de características**.
> > 
> {: .solution}
> 
{: .question}

¿De dónde procede la información del visor de estructura?

> <hands-on-title>Búsqueda de OPN1LW humano en UniProtKB</hands-on-title>
> 
> 1. Haga clic en el icono de descarga situado debajo de la estructura
> 2. Comprueba el archivo que se ha descargado
> 
{: .hands-on}

Se trata de un archivo PDB (Protein Data Bank), que permite visualizar y analizar la disposición de los átomos y aminoácidos de la proteína.

Sin embargo, no hay ninguna referencia a la base de datos PDB en los enlaces entre las *bases de datos de estructuras 3D*. En su lugar, el primer enlace hace referencia a la AlphaFoldDB. La base de datos AlphaFold es un recurso completo que proporciona estructuras 3D predichas para una amplia gama de proteínas. Utilizando técnicas de aprendizaje profundo e información evolutiva, AlphaFold predice con precisión la disposición espacial de los átomos dentro de una proteína, contribuyendo a nuestra comprensión de la función y las interacciones de las proteínas.

Por lo tanto, se trata de una *predicción* de la estructura, no de una estructura validada experimentalmente. Esta es la razón por la que está coloreada por confianza: las secciones en azul son las que tienen un valor de confianza alto, por lo que son aquellas para las que la predicción es muy fiable, mientras que las que están en naranja son menos reilizables o tienen una estructura desordenada (más flexible y móvil). No obstante, esta información se representa a través de un archivo PDB, porque sigue siendo estructural.

### Familia y dominios

La sección *Familia y Dominios* de la página de entrada de UniProtKB proporciona una visión completa de las relaciones evolutivas y los dominios funcionales dentro de una proteína. Esta sección ofrece información sobre la pertenencia de la proteína a familias de proteínas, superfamilias y dominios, arrojando luz sobre sus características estructurales y funcionales.

El área *Features* confirma efectivamente que al menos uno de los dos dominios que sobresalen de la membrana (el N-terminal) está desordenado. Esta área suele incluir información sobre regiones conservadas, motivos y características importantes de la secuencia que contribuyen al papel de la proteína en diversos procesos biológicos. La sección confirma una vez más que estamos ante una proteína transmembrana, y enlaza con varios recursos de datos filogenéticos, de familias de proteínas o de dominios, que nos guían en la comprensión de cómo las proteínas comparten ancestros comunes, evolucionan y adquieren funciones especializadas.

### Secuencia

Toda esta información sobre la evolución de la proteína, su función, su estructura, está codificada en última instancia en su secuencia. Una vez más, en esta sección tenemos la oportunidad de descargar el archivo FASTA que la transcribe, así como de acceder a la fuente de estos datos: los experimentos de secuenciación genómica que la evaluaron. En esta sección también se informa de cuándo se han detectado isoformas.

> <question-title></question-title>
> 
> ¿Cuántas isoformas potenciales están asignadas a esta entrada?
> 
> > <solution-title></solution-title>
> > 
> > 1: H0Y622
> > 
> > 
> {: .solution}
> 
{: .question}

### Proteínas similares

La última sección de la página de entrada de UniProt informa de proteínas similares (esto es básicamente el resultado de una agrupación, con umbrales de identidad del 100%, 90% y 50%).

> <question-title></question-title>
> 
> 1. ¿Cuántas proteínas similares al 100% de identidad?
> 2. ¿Cuántas proteínas similares al 90% de identidad?
> 3. ¿Cuántas proteínas similares al 50% de identidad?
> 
> > <solution-title></solution-title>
> > 
> > 1. 0
> > 2. 83
> > 3. 397
> > 
> {: .solution}
> 
{: .question}

Como habrá adivinado al ver esta página, gran parte del procesamiento de datos biológicos sobre una proteína consiste en mapear distintos tipos de información a lo largo de la secuencia y comprender cómo se influyen mutuamente. Un mapeo visual (y una tabla con la misma información) es proporcionado por las dos pestañas alternativas para ver esta entrada, que son el *Visor de Variantes* y el *Visor de Características*.

## Visor de variantes

> <hands-on-title>Visor de variantes</hands-on-title>
> 
> 1. Haga clic en la pestaña *Visor de variantes*
> 
{: .hands-on}

El *Visor de Variantes* mapea todas las versiones alternativas conocidas de esta secuencia. Para algunas de ellas se conoce el efecto (patógeno o benigno), para otras no.

> <question-title></question-title>
> 
> ¿Cuántas variantes son probablemente patógenas?
> 
> > <solution-title></solution-title>
> > 
> > Al alejar la vista de variantes, vemos que tenemos 5 puntos rojos, es decir, 5 variantes probablemente patógenas.
> > 
> > 
> {: .solution}
> 
{: .question}

El elevado número de variantes que se encuentran en esta sección sugiere que las "secuencias de proteínas" (así como las secuencias de genes, las estructuras de proteínas, etc.) son en realidad entidades menos fijas de lo que podríamos pensar.

## Visor de características

> <hands-on-title>Visor de características</hands-on-title>
> 
> 1. Haga clic en la pestaña *Visor de características*
> 
{: .hands-on}

El *Visor de características* es básicamente una versión fusionada de todas las áreas de *Características* que encontramos en la página *Entrada*, incluyendo *Dominios y sitios*, *Procesamiento de moléculas*, *PTMs*, *Topología*, *Proteómica*, *Variantes*. Si en el visor hace clic en cualquier característica, se enfocará la región correspondiente en la estructura, como la variante de interés

> <hands-on-title>Visor de variantes</hands-on-title>
> 
> 1. Expandir la parte *Variantes
> 2. Zoom out
> 3. Haga clic en nuestra variante de interés (el punto rojo en la posición 338)
> 
{: .hands-on}

> <question-title></question-title>
> 
> ¿Cuál es la topología en este lugar?
> 
> > <solution-title></solution-title>
> > 
> > Un dominio citoplasmático topológico
> > 
> > 
> {: .solution}
> 
{: .question}

Por último, echemos un vistazo rápido a las otras pestañas.

## Publicaciones

> <hands-on-title>Publicación</hands-on-title>
> 
> 1. Haga clic en la pestaña *Publicación*
> 
{: .hands-on}

En *Publicaciones* se enumeran las publicaciones científicas relacionadas con la proteína. Éstas se recopilan fusionando una lista completamente curada en UniProtKB/Swiss-Prot y otras importadas automáticamente. En esta pestaña, puede filtrar la lista de publicaciones por fuente y categorías que se basan en el tipo de datos que contiene una publicación sobre la proteína (como función, interacción, secuencia, etc.), o por el número de proteínas en el estudio correspondiente que describe ("pequeña escala" frente a "gran escala").

> <question-title></question-title>
> 
> 1. ¿Cuántas publicaciones hay asociadas a esta proteína?
> 2. ¿Cuántas publicaciones contienen información sobre su función?
> 
> > <solution-title></solution-title>
> > 
> > 1. 57
> > 2. 23
> > 
> {: .solution}
> 
{: .question}

## Enlaces externos

> <hands-on-title>Enlaces externos</hands-on-title>
> 
> 1. Haga clic en la pestaña *Enlaces externos*
> 
{: .hands-on}

La pestaña *Enlaces externos* reúne todas las referencias a bases de datos y recursos de información externos que encontramos en cada sección de la página de entrada. El texto de los enlaces suele informar de los identificadores únicos que representan a la misma entidad biológica en otras bases de datos. Para hacerse una idea de esta compexidad, consulte la siguiente imagen (que ya está parcialmente desfasada).

![Un gráfico que representa cómo todas las diferentes bases de datos están conectadas por IDs únicos, donde los nodos mayores son DBs, y las flechas los conectan con los IDs (nodos menores) que reportan. El mapa está muy abarrotado, especialmente alrededor del Nombre de entrada de UniProt, la ID de gen y la ID de gen de Ensembl](./images/complexDB.jpeg "Mejor ilustración de la ID compleja crossref: bioDBnet Network Diagram - [fuente](https://biodbnet-abcc.ncifcrf.gov/dbInfo/netGraph.php)")

## Historia

Por último, la pestaña *Historia* también es interesante. Informa y pone a disposición para su descarga todas las versiones anteriores de las anotaciones de esta entrada, es decir: toda la "evolución" de su anotación, que en este caso se remonta a 1988.

> <question-title></question-title>
> 
> ¿Esta entrada no ha sido anotada manualmente?
> 
> > <solution-title></solution-title>
> > 
> > Para responder a esta pregunta puede desplazarse hacia atrás en el tiempo por la tabla y comprobar la columna `Database`. ¿Estuvo alguna vez en TrEMBL en lugar de en SwissProt? No, por lo que esta entrada fue anotada manualmente desde su inicio.
> > 
> > 
> {: .solution}
> 
{: .question}


