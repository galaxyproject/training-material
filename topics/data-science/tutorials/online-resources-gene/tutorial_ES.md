---
layout: tutorial_hands_on
title: Aprendizaje sobre un gen a través de recursos y formatos biológicos
level: Introductory
draft: true
zenodo_link: https://zenodo.org/record/8304465
questions:
- ¿Cómo utilizar recursos de bioinformática para investigar una familia específica de proteínas (opsinas)?
- ¿Cómo navegar por el Genome Data Viewer para encontrar opsinas en el genoma humano?
- ¿Cómo identificar genes asociados con las opsinas y analizar su localización en los cromosomas?
- ¿Cómo explorar la literatura científica y los contextos clínicos relacionados con el gen OPN1LW?
- ¿Cómo usar archivos de secuencias proteicas para realizar búsquedas de similitud utilizando BLAST?
objectives:
- Partiendo de una búsqueda por texto, navega por múltiples recursos web para examinar distintos tipos de información sobre un gen, presentada en diversos formatos de archivo. 
time_estimation: 1H
key_points:
- Puedes buscar genes y proteínas usando texto específico en el genoma de NCBI.
- Una vez encuentres un gen o proteína relevante, puedes obtener su secuencia y anotación en varios formatos desde NCBI.
- También puedes conocer la localización cromosómica y la composición de exones‑intrones del gen de interés.
- NCBI ofrece una herramienta BLAST para realizar búsquedas de similitud con secuencias.
- Puedes explorar más recursos incluidos en este tutorial para aprender sobre condiciones asociadas al gen y sus variantes.
- Puedes introducir un archivo FASTA con una secuencia de interés para búsquedas BLAST.
contributions:
  authorship:
  - lisanna
  - bebatut
  - teresa-m
  translation:
  - hvelab
  - unode
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


Cuando realizamos un análisis bioinformático, por ejemplo RNA-seq, podemos acabar con una lista de nombres de genes. A continuación, necesitamos explorar estos genes. ¿Pero cómo podemos hacerlo? ¿Cuáles son los recursos disponibles para ello? ¿Y cómo navegar por ellos?

El objetivo de este tutorial es familiarizarnos con ello, utilizando como ejemplo las opsinas humanas.

Las opsinas humanas se encuentran en las células de la retina. Las opsinas captan la luz e inician la secuencia de señales que dan lugar a la visión. Procederemos formulando preguntas sobre las opsinas y los genes de opsina, y luego utilizaremos diferentes bases de datos y recursos bioinformáticos para responderlas.

> <comment-title></comment-title>
> Este tutorial es un poco atípico: no trabajaremos en Galaxy sino principalmente fuera de él, navegando por bases de datos y herramientas a través de sus propios interfaces web. El objetivo de este tutorial es ilustrar varias fuentes de datos biológicos en diferentes formatos de archivo, y que representan información diferente.
{: .comment}

> <agenda-title></agenda-title>
> 
> En este tutorial trataremos:
> 
> 1. TOC
> {:toc}
{: .agenda}

# Búsqueda de Opsinas Humanas

Para buscar Opsinas humanas, empezaremos por consultar el [Visor de Datos del Genoma del NCBI](https://www.ncbi.nlm.nih.gov/genome/gdv). El Visor de Datos del Genoma del NCBI (GDV) ({% cite rangwala2021accessing %}) es un navegador del genoma que permite la exploración y el análisis de ensamblajes de genomas eucariotas anotados. El navegador GDV muestra información biológica asignada a un genoma, incluida la anotación de genes, datos de variación, alineaciones BLAST y datos de estudios experimentales de las bases de datos NCBI GEO y dbGaP. Las notas de la versión GDV describen las nuevas características relacionadas con este navegador.

> <hands-on-title>Abrir el Visor de Datos del Genoma del NCBI</hands-on-title>
> 
> 1. Abra el visor de datos del genoma del NCBI en [www.ncbi.nlm.nih.gov/genome/gdv](https://www.ncbi.nlm.nih.gov/genome/gdv/)
{: .hands-on}

La página de inicio incluye un simple "árbol de la vida" donde el nodo humano está resaltado porque es el organismo por defecto para buscar. Podemos cambiar eso en la casilla *Buscar organismos* pero lo dejaremos por ahora ya que estamos interesados en las Opsinas Humanas.

![Captura de pantalla de la página de inicio del Visor de Datos del Genoma, se escribe la palabra "opsina" en la caja de búsqueda y se previsualiza el resultado](./images/GenomeDataViewerpage.png "Página de inicio del Visor de Datos del Genoma")

El panel de la derecha informa de múltiples ensamblajes del genoma de interés, y un mapa de los cromosomas en ese genoma. Allí podemos buscar Opsinas.

> <hands-on-title>Buscar Opsinas - Abrir el Visor de Datos del Genoma del NCBI</hands-on-title>
> 
> 1. Escriba `opsin` en la casilla *Buscar en genoma*
> 2. Haga clic en el icono de la lupa o pulse <kbd>Intro<kbd>
{: .hands-on}

Debajo del cuadro se muestra ahora una tabla con genes relacionados con la opsina junto con sus nombres y localización, es decir, el número de cromosoma, así como la posición inicial y final

En la lista de genes relacionados con el término de búsqueda opsina, aparecen el gen de la rodopsina (RHO), y tres pigmentos de los conos, opsinas sensibles a longitudes de onda corta, media y larga (para la detección de luz azul, verde y roja). Existen otras entidades, por ejemplo una -LCR (Locus Control region), genes putativos y receptores.

Múltiples coincidencias en el cromosoma X, uno de los cromosomas que determinan el sexo.

> <question-title></question-title>
> 
> 1. ¿Cuántos genes se han encontrado en el cromosoma X?
> 2. ¿Cuántos son los genes codificadores de proteínas?
> 
> > <solution-title></solution-title>
> > 
> > 1. Los resultados en ChrX son:
> > - OPSIN-LCR
> > - OPN1LW
> > - OP1MW
> > - OPN1MW2
> > - OPN1MW3
> > 
> > 2. Al pasar el ratón por encima de cada gen, se abre un cuadro y podemos hacer clic en *Detalles* para saber más sobre cada gen. Entonces nos enteramos de que el primero (OPSIN-LCR) no codifica proteínas, sino que es una región reguladora de genes, y que los otros son genes codificadores de proteínas. Así que hay 4 genes codificadores de proteínas relacionados con las opsinas en el cromosoma X. En particular, el cromosoma X incluye un gen de pigmento rojo (OPN1LW) y tres genes de pigmento verde (OPN1MW, OPN1MW2 y OPN1MW3 en el ensamblaje del genoma de referencia).
> {: .solution}
{: .question}

Centrémonos ahora en una opsina específica, el gen OPN1LW.

> <hands-on-title>Buscador del genoma abierto para el gen OPN1LW</hands-on-title>
> 
> 1. Haga clic en la flecha azul que aparece en la tabla de resultados al pasar el ratón por encima de la fila OPN1LW
{: .hands-on}

Debería haber aterrizado en [esta página](https://www.ncbi.nlm.nih.gov/genome/gdv/browser/genome/?id=GCF_000001405.40), que es la vista del genoma del gen OPN1LW.

![Visor de datos genómicos del gen OPN1LW, captura de pantalla de los dos paneles principales del visor, con los cromosomas a la izquierda y el visor de características a la derecha](./images/GenomeDataViewerofgeneOPN1LW.png "Visor de datos genómicos del gen OPN1LW")

Hay mucha información en esta página, centrémonos en una sección cada vez.

1. El Visor de Datos del Genoma, en la parte superior, nos dice que estamos viendo los datos del organismo `Homo sapiens`, ensamblaje `GRCh38.p14` y en particular en `Chr X` (Cromosoma X). Cada una de estas informaciones tiene un ID único.
2. Todo el cromosoma está representado directamente debajo, y las posiciones a lo largo de los brazos corto (`p`) y largo (`q`) están resaltadas.
3. A continuación, un recuadro azul resalta que ahora nos centramos en la Región correspondiente al Gen `OPN1LW`.

   Hay múltiples formas de interactuar con el visor de abajo. Pruebe, por ejemplo, a pasar el ratón por encima de los puntos que representan exones en el recuadro azul.

4. En el gráfico siguiente, la requencia génica es una línea verde con los exones (fragmentos codificantes de proteínas) representados por rectángulos verdes.

   Pase el ratón sobre la línea verde correspondiente a `NM_020061.6` (nuestro gen de interés) para obtener información más detallada.

   > <question-title></question-title>
   > 
   > 1. ¿Cuál es la ubicación del segmento OPN1LW?
   > 2. ¿Cuál es la longitud del segmento OPN1LW?
   > 3. ¿Qué son los intrones y los exones?
   > 4. ¿Cuántos exones e intrones hay en el gen OPN1LW?
   > 5. ¿Cuál es la longitud total de la región codificante?
   > 6. ¿Cuál es la distribución entre regiones codificantes y no codificantes? ¿Qué significa eso en el término de biología?
   > 7. ¿Cuál es la longitud de la proteína en número de aminoácidos?
   > 
   > > <solution-title></solution-title>
   > > 
   > > 1. De 154.144.243 a 154.159.032
   > > 2. 1.4790 nucleótidos, encontrados en *Span on 14790 nt, nucleotides)*
   > > 3. Los genes eucarióticos suelen estar interrumpidos por regiones no codificantes denominadas secuencias intermedias o intrones. Las regiones codificantes se denominan exones.
   > > 4. En este diagrama se puede ver que el gen OPN1LW consta de 6 exones y 5 intrones, y que los intrones son mucho más grandes que los exones.
   > > 5. La longitud de la CDS es de 1.095 nucleótidos.
   > > 6. De los 14790 nt del gen, sólo 1095 nt codifican proteínas, lo que significa que menos del 8% de los pares de bases contienen el código. Cuando este gen se expresa en las células de la retina humana, se sintetiza una copia de ARN de todo el gen. A continuación, se recortan las regiones de intrones y se unen las regiones de exones para producir el ARNm maduro (un proceso denominado splicing), que será traducido por los ribosomas al fabricar la proteína opsina roja. En este caso, se desecha el 92% del transcrito inicial de ARN, dejando el código proteico puro.
   > > 7. La longitud de la proteína resultante es de 364 aa, aminoácidos.
   > {: .solution}
   {: .question}

Pero, ¿cuál es la secuencia de este gen? Hay múltiples formas de recuperar esta información, vamos a repasar la que creemos que es una de las más intituitivas.

> <hands-on-title>Buscador del genoma abierto para el gen OPN1LW</hands-on-title>
> 
> 1. Haga clic en el icono {% icon tool %} *Herramientas* en la parte superior derecha del cuadro que muestra el gen
> 2. Haga clic en *Vista de secuencias de texto*
{: .hands-on}

Este panel informa de la secuencia de ADN de los intrones (en verde), así como de la de los exones (en rosa, incluyendo la secuencia de proteína traducida más abajo).

![Captura de pantalla de la vista de secuencia del recurso NHI, el texto aparece resaltado en diferentes colores](./images/SequenceTextView.png "Vista de texto de secuencia")

Este cuadro de secuencia no muestra el gen completo en este momento, sino una subsecuencia del mismo. Puede moverse aguas arriba y aguas abajo del código genético con las flechas *Página anterior* y *Página siguiente*, o empezar desde una posición específica con el botón *Ir a la posición*. Sugerimos empezar con el inicio de la parte codificante del gen, que como aprendimos antes está en la posición 154,144,243.

> <hands-on-title>Ir a una posición específica en la Vista de Secuencia</hands-on-title>
> 
> 1. Haga clic en *Ir a la posición*
> 2. Escriba en `154144243`
> 
>    Tenemos que quitar las comas para validar el valor
{: .hands-on}

La secuencia resaltada en púrpura señala aquí una región reguladora.

> <question-title></question-title>
> 
> 1. ¿Cuál es el primer aminoácido del producto proteico resultante?
> 2. ¿Cuál es la última?
> 3. ¿Puede anotar los tres primeros y los tres últimos AA de esta proteína?
> 
> > <solution-title></solution-title>
> > 
> > 1. La proteína correspondiente empieza por Metionina, M (todas lo hacen).
> > 2. El último AA del último exón (que se encuentra en la 2ª página) es Alanina (A). Después viene el codón de parada TGA, que no se traduce en un AA.
> > 3. Los tres primeros AA son: M,A,Q; los tres últimos: S,P,A.
> {: .solution}
{: .question}

Ahora podemos cerrar la *Vista Secuencia*.

Desde este recurso, también podemos obtener archivos, en diferente formato, que describen el gen. Están disponibles en la sección *Download*.

1. *Descargar FASTA* nos permitirá descargar el formato de archivo más sencillo para representar la secuencia de nucleótidos de todo el rango visible del genoma (más largo que el gen solamente).
2. *Download GenBank flat file* nos permitirá acceder a la anotación disponible en esta página (y más allá) en un formato de texto plano.
3. *Download Track Data* nos permite inpectar dos de los formatos de archivo que presentamos en las diapositivas: los formatos GFF (GFF3) y BED. Si se cambian las pistas, cada una puede o no estar disponible.

# Encontrar más información sobre nuestro gen

Veamos ahora la información que tenemos (en la literatura) sobre nuestro gen, utilizando los recursos del NCBI

> <hands-on-title>Ir a una posición específica en la Vista de Secuencia</hands-on-title>
> 
> 1. Abra la búsqueda NCBI en [www.ncbi.nlm.nih.gov/search](https://www.ncbi.nlm.nih.gov/search/)
> 2. Escriba `OPN1LW` en la casilla de búsqueda *Search NCBI
{: .hands-on}

![Captura de pantalla de la página de resultados del NIH, con tarjetas denominadas Literature, Genes, Proteins, genomes, Clinical y PubChem](./images/NIHresults.png "Results when searching for `OPN1LW` on NCBI").

## Literatura

Empecemos con la literatura y en particular con los resultados de *PubMed* o *PubMed Central

> <details-title>¿Cuál es la diferencia entre PubMed y PubMed Central? </details-title>
> 
> PubMed es una base de datos de literatura biomédica que contiene los resúmenes de las publicaciones de la base de datos.
> 
> PubMed Central es un repositorio de texto completo, que contiene el texto completo de las publicaciones de la base de datos.
> 
> Aunque el número exacto de aciertos puede variar en el tiempo con respecto a la captura de pantalla anterior, cualquier nombre de gen debería tener más aciertos en PubMed Central (búsqueda en los textos completos de las publicaciones) que en PubMed (búsqueda sólo en los resúmenes).
{: .details}

> <hands-on-title>Abrir PubMed</hands-on-title>
> 
> 1. Haga clic en *PubMed* en la casilla *Literatura
{: .hands-on}

Ha introducido en PubMed, una base de datos gratuita de literatura científica, los resultados de una búsqueda completa de artículos directamente asociados a este locus génico.

Haciendo clic en el título de cada artículo, podrá ver resúmenes del mismo. Si se encuentra en un campus universitario donde hay acceso en línea a revistas específicas, también podrá ver enlaces a artículos completos. PubMed es su punto de entrada a una amplia variedad de literatura científica en ciencias de la vida. En la parte izquierda de cualquier página de PubMed encontrará enlaces a una descripción de la base de datos, ayuda y tutoriales sobre búsquedas.

> <question-title></question-title>
> 
> 1. ¿Puede adivinar qué tipo de enfermedades están asociadas a este gen?
> 
> > <solution-title></solution-title>
> > 
> > 1. Responderemos a esta pregunta más adelante
> {: .solution}
{: .question}

> <hands-on-title>Volver a la página de búsqueda del NCBI</hands-on-title>
> 
> 1. Volver a la [página de búsqueda del NCBI](https://www.ncbi.nlm.nih.gov/search/all/?term=OPN1LW)
{: .hands-on}

## Clínica

Centrémonos ahora en el cuadro *Clinical*, y especialmente en *OMIM*. OMIM, el Online Mendeliam Inheritance in Man (¡y woman!), es un catálogo de genes humanos y trastornos genéticos.

> <hands-on-title>Abrir OMIM</hands-on-title>
> 
> 1. Haga clic en *OMIM* en la casilla *Clinical
{: .hands-on}

Cada entrada OMIM es un trastorno genético (aquí principalmente tipos de daltonismo) asociado con mutaciones en este gen.

> <hands-on-title>Lea tanto como le dicte su interés</hands-on-title>
> 
> 1. Siga los enlaces para obtener más información sobre cada entrada
{: .hands-on}

> <comment-title>Lea tanto como le dicte su interés</comment-title>
> 
> Para obtener más información sobre OMIM, haga clic en el logotipo de OMIM en la parte superior de la página. A través de OMIM, se dispone de abundante información sobre innumerables genes del genoma humano, y toda la información está respaldada por referencias a los artículos de investigación más recientes.
{: .comment}

¿Cómo afectan las variaciones en el gen al producto proteico y a sus funciones? Volvamos a la página de los NIH e investiguemos el acceso a la lista de Polimorfismos de Nucleótido Único (SNPs) que fueron detectados por estudios genéticos en el gen.

> <hands-on-title>Abrir dbSNP</hands-on-title>
> 
> 1. Volver a la [página de búsqueda del NCBI](https://www.ncbi.nlm.nih.gov/search/all/?term=OPN1LW)
> 2. Haga clic en *dbSNP* en la casilla *Clinical
{: .hands-on}

![Captura de pantalla de la página de dbSNPs sobre el gen OPN1LW. Tres paneles principales, el de la izquierda para filtrar la búsqueda basada en etiquetas, el central mostrando resultados, el de la derecha para una búsqueda más detallada y programática](./images/dbSNPs.png "dbSNP en OPN1LW")

> <question-title></question-title>
> 
> 1. ¿Cuál es el significado clínico de las rs5986963 y rs5986964 (las 2 primeras variantes listadas en el momento de creación de este tutorial)?
> 2. ¿Cuál es la consecuencia funcional de rs104894912?
> 3. ¿Cuál es la consecuencia funcional de rs104894913?
> 
> > <solution-title></solution-title>
> > 
> > 1. La significación clínica es `benign` por lo que parece que no tienen efecto sobre el producto proteico final
> > 2. la mutación rs104894912 conduce a una variante `stop_gained`, que trunca la proteína resultante demasiado pronto y es por tanto `pathogenic`
> > 3. la mutación rs104894913 conduce a un `missense_variant`, también `pathogenic`.
> {: .solution}
{: .question}

Investiguemos más sobre la variante rs104894913

> <hands-on-title>Aprenda más sobre una variante dbSNP</hands-on-title>
> 
> 1. Haga clic en `rs104894913` para abrir su [página dedicada](https://www.ncbi.nlm.nih.gov/snp/rs104894913)
> 2. Haga clic en *Clinical Significance* (Importancia clínica)
> 
>    > <question-title></question-title>
>    > 
>    > ¿Qué tipo de afección se asocia con la variante rs104894913?
>    > 
>    > > <solution-title></solution-title>
>    > > 
>    > > El nombre de la enfermedad asociada es "defecto de Protan". Una rápida búsqueda en Internet con su buscador le aclarará que se trata de un tipo de daltonismo.
>    > {: .solution}
>    {: .question}
> 
> 3. Haga clic en *Detalles de la variante*
> 
>    > <question-title></question-title>
>    > 
>    > 1. ¿Qué sustitución está asociada a esta variante?
>    > 2. ¿Cuál es el impacto de esta subtitución en términos de codón y aminoácido?
>    > 3. ¿En qué posición de la proteína se encuentra esta sustitución?
>    > 
>    > > <solution-title></solution-title>
>    > > 
>    > > 1. La sustitución `NC_000023.10:g.153424319G>A` corresponde al cambio de una Guanina (G) a una Adenina (A)
>    > > 2. Esta sustitución cambia el codón `GGG`, una Glicina, por `GAG`, un Glutatión
>    > > 3. `p.Gly338Glu` significa que la sustitución se encuentra en la posición 338 de la proteína.
>    > {: .solution}
>    {: .question}
{: .hands-on}

¿Qué significa esta sustitución para la proteína? Echemos un vistazo más profundo a esta proteína.

## Proteína

> <hands-on-title>Proteína abierta</hands-on-title>
> 
> 1. Volver a la [página de búsqueda del NCBI](https://www.ncbi.nlm.nih.gov/search/all/?term=OPN1LW)
> 2. Haga clic en *Proteína* en la casilla *Proteínas
> 3. Haga clic en `OPN1LW – opsin 1, long wave sensitive` en el cuadro superior
{: .hands-on}

![Captura de pantalla de la página de la proteína Opsin 1 NIH, dos paneles principales. El de la izquierda informa sobre el gen, el de la derecha es una tabla de contenido y una serie de enlaces a otros recursos](./images/Opsin1NIH.png "Opsin 1 NIH protein page")

Esta página presenta de nuevo algunos datos con los que estamos familiarizados (por ejemplo, la distribución de los exones a lo largo de la secuencia del gen).

> <hands-on-title>Descargar las secuencias de proteínas</hands-on-title>
> 
> 1. Haga clic en *Descargar conjuntos de datos*
> 2. Seleccionar
>    - `Gene Sequences (FASTA)`
>    - `Transcript sequences (FASTA)`
>    - `Protein sequences (FASTA)`
> 3. Haga clic en el botón *Descargar
> 4. Abra el archivo ZIP descargado
{: .hands-on}

> <question-title></question-title>
> 
> 1. ¿Qué contiene la carpeta?
> 2. ¿Cree que han aplicado buenas prácticas de datos?
> 
> > <solution-title></solution-title>
> > 
> > 1. La carpeta incluye
> >    - una carpeta `ncbi_datasets` con diferentes subcarpetas en ella leadig algunos archivos de datos (múltiples formatos),
> >    - un `README.md` (un archivo Markdown), que está diseñado para "viajar" junto con los datos y explicar cómo se recuperaron los datos, cuál es la estructura de la subcarpeta que contiene los datos, y dónde encontrar documentación extensa.
> > 2. Sin duda, es una buena práctica de gestión de datos guiar a los usuarios (no sólo a sus colaboradores, sino también a usted mismo en un futuro no muy lejano, cuando olvide de dónde procede ese archivo de su carpeta Descargas) hasta la fuente de datos y la estructura de datos.
> {: .solution}
{: .question}

# Búsqueda por secuencia

¿Qué podríamos hacer con estas secuencias que acabamos de descargar? Supongamos que acabamos de secuenciar los transcritos que hemos aislado mediante un experimento, de modo que conocemos la secuencia de nuestra entidad de interés, pero no sabemos cuál es. Lo que tenemos que hacer en este caso es buscar en toda la base de datos de secuencias conocidas por la ciencia y emparejar nuestra entidad desconocida con una entrada que tenga alguna anotación. Hagámoslo.

> <hands-on-title>Busca la secuencia de proteína contra todas las secuencias de proteína</hands-on-title>
> 
> 1. Abra (con el editor de texto más sencillo que tenga instalado) el archivo `protein.faa` que acaba de descargar.
> 2. Copie su contenido
> 3. BLAST abierto [blast.ncbi.nlm.nih.gov](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
> 4. Haga clic en `Protein BLAST, protein > protein`
> 
>    En efecto, utilizaremos una secuencia de proteínas para buscar en una base de datos de proteínas
> 
> 5. Pegue la secuencia de la proteína en el cuadro de texto grande
> 6. Compruebe el resto de parámetros
> 7. Haga clic en el botón azul `BLAST`
{: .hands-on}

Esta fase llevará algún tiempo, después de todo hay algún servidor en alguna parte que está comparando la totalidad de las secuencias conocidas con su objetivo. Cuando la búsqueda se haya completado, el resultado debe ser similar a la de abajo:

![Captura de pantalla de los resultados de BLAST, un gran encabezado en la parte superior y los resultados listados como una tabla en la parte inferior](./images/BLASTresults.png "BLAST results")

> <hands-on-title>Resumen gráfico de las secuencias de proteínas</hands-on-title>
> 
> 1. Haga clic en la pestaña *Resumen gráfico*
{: .hands-on}

Accedemos a un cuadro que contiene muchas líneas de colores. Cada línea representa un acierto de su búsqueda blast. Si hace clic en una línea roja, el recuadro estrecho que aparece justo encima ofrece una breve descripción del resultado.

> <hands-on-title>Descripciones de las secuencias de proteínas</hands-on-title>
> 
> 1. Haga clic en la pestaña *Descripciones*
{: .hands-on}

> <question-title></question-title>
> 
> 1. ¿Cuál es el primer resultado? ¿Se espera?
> 2. ¿Cuáles son los otros resultados? ¿Para qué organismos?
> 
> > <solution-title></solution-title>
> > 
> > 1. El primer acierto es nuestra opsina roja. Eso es alentador, porque la mejor coincidencia debe ser a la secuencia de consulta en sí, y usted tiene esta secuencia de esa entrada de genes.
> > 2. Otros hits son otras opsinas. Incluyen entradas de otros primates (por ejemplo, `Pan troglogytes`).
> {: .solution}
{: .question}

Los resultados son para nuestra opsina roja en humanos pero también para otras opsinas en otros primates. Podríamos querer eso, por ejemplo si quisiéramos utilizar estos datos para construir un árbol filogenético. Si en cambio estamos bastante seguros de que nuestra secuencia de interés es humana, también podríamos haber filtrado la búsqueda sólo en secuencias humanas.

> <hands-on-title>Filtrar una búsqueda BLAST</hands-on-title>
> 
> 1. Haga clic en *Editar búsqueda*
> 2. Escriba `Homo sapiens` en el campo *Organism*
> 3. Haga clic en el botón azul `BLAST`
{: .hands-on}

Con esta nueva búsqueda, encontramos las otras opsinas (verde, azul, pigmento de células bastón) en la lista. Los otros resultados tienen un número menor de residuos coincidentes. Si hace clic en cualquiera de las líneas de color en el *Resumen Gráfico*, abrirá más información sobre ese resultado, y podrá ver cuánta similitud tiene cada uno con la opsina roja, nuestra secuencia de consulta original. A medida que se desciende en la lista, cada secuencia sucesiva tiene menos en común con la opsina roja. Cada secuencia se muestra en comparación con la opsina roja en lo que se denomina un alineamiento de secuencias por pares. Más adelante, realizará alineaciones de secuencias múltiples a partir de las cuales podrá discernir las relaciones entre los genes.

> <details-title>Más detalles sobre las puntuaciones BLAST</details-title>
> 
> Las visualizaciones contienen dos medidas destacadas de la importancia del acierto:
> 
> 1. la puntuación BLAST - lableled Score (bits)
> 
>    La puntuación BLAST indica la calidad del mejor alineamiento entre la secuencia buscada y la secuencia encontrada (hit). Cuanto mayor sea la puntuación, mejor será el alineamiento. Las puntuaciones se ven reducidas por los desajustes y las lagunas en el mejor alineamiento. El cálculo de la puntuación es complejo e implica una matriz de sustitución, que es una tabla que asigna una puntuación a cada par de residuos alineados. La matriz más utilizada para la alineación de proteínas se conoce como BLOSUM62.
> 
> 2. el Valor de Expectativa (etiquetado Expect o E)
> 
>     El valor de la expectativa E de un acierto indica si el acierto es probablemente el resultado de la semejanza casual entre el acierto y la consulta, o de la ascendencia común del acierto y la consulta. ()
> 
>     > <comment-title>Filtrar una Búsqueda BLAST</comment-title>
>     > 
>     > Si E es menor que $$10\mathrm{e}{-100}$$, a veces se da como 0.0.
>     {: .comment}
> 
>     El valor de expectativa es el número de aciertos que esperarías que se produjeran por pura casualidad si buscaras tu secuencia en un genoma aleatorio del tamaño del genoma humano.
>
>      $$E = 25$$ significa que podría esperar encontrar 25 coincidencias en un genoma de este tamaño, puramente por azar. Por lo tanto, una coincidencia con $$E = 25$$ es probablemente una coincidencia casual, y no implica que la secuencia encontrada comparta ascendencia común con la secuencia buscada.
>
>      Los valores de expectativa de alrededor de 0,1 pueden o no ser biológicamente significativos (se necesitarían otras pruebas para decidirlo).
>
>      Pero valores muy pequeños de E significan que la coincidencia es biológicamente significativa. La correspondencia entre su secuencia de búsqueda y esta coincidencia debe surgir de la ascendencia común de las secuencias, porque las probabilidades son simplemente demasiado bajas de que la coincidencia pueda surgir por casualidad. Por ejemplo, $$E = 10\mathrm{e}{-18}$$ para una coincidencia en el genoma humano significa que sólo se esperaría una coincidencia casual en un billón de billones de genomas diferentes del mismo tamaño que el genoma humano.
>
>      La razón por la que creemos que todos procedemos de antepasados comunes es que la similitud masiva de secuencias en todos los organismos es sencillamente demasiado improbable para ser un hecho fortuito. Cualquier familia de secuencias similares en muchos organismos debe haber evolucionado a partir de una secuencia común en un antepasado remoto.
> 
{: .details}

> <hands-on-title>Descargando</hands-on-title>
> 
> 1. Haga clic en la pestaña *Descripciones
> 2. Haga clic en cualquier coincidencia de secuencia
> 3. Haga clic en *Descargar*
> 4. Select `FASTA (aligned sequences)`
> 
{: .hands-on}

Descargará un nuevo tipo de archivo, ligeramente diferente: un FASTA alineado. Si lo desea, explórelo antes de la siguiente sección.

Mientras que en las secciones anteriores de este tutorial hemos utilizado ampliamente las interfaces web de las herramientas (visores genómicos, exploración rápida de la literatura, lectura de anotaciones, etc.), esta búsqueda BLAST es un ejemplo de un paso que podría automatizar completamente con Galaxy.

> <hands-on-title> Búsqueda de similitud con BLAST en Galaxia </hands-on-title>
> 
> 1. Crear un nuevo historial para este análisis
> 
>    {% snippet faqs/galaxy-es/histories_create_new.md %}
> 
> 2. Cambiar el nombre del historial
> 
>    {% snippet faqs/galaxy-es/histories_rename.md %}
> 
> 3. Importar la secuencia de proteínas a través de enlace desde [Zenodo]({{ page.zenodo_link }}) o bibliotecas de datos compartidos Galaxy:
> 
>    ```text
>    {{ page.zenodo_link }}/files/protein.faa
>    ```
> 
>    {% snippet faqs/galaxy-es/datasets_import_via_link.md %}
> 
>    {% snippet faqs/galaxy-es/datasets_import_from_data_library.md %}
> 
> 1. {% tool [NCBI BLAST+ blastp](toolshed.g2.bx.psu.edu/repos/devteam/ncbi_blast_plus/ncbi_blastp_wrapper/2.10.1+galaxy2) %} con los siguientes parámetros:
>    - *"Protein query sequence(s) "*: `protein.faa`
>    - *"Subject database/sequences "*: `Locally installed BLAST database`
>    - *"Protein BLAST database "*: `SwissProt`
> 
>      Para buscar sólo secuencias anotadas en UniProt, necesitamos seleccionar la última versión de `SwissProt`
> 
>    - *"Set expectation value cutoff "*: `0.001`
>    - *"Output format "*: `Tabular (extended 25 columns)`
> 
{: .hands_on}

> <question-title></question-title>
> 
> ¿Cree que estamos viendo exactamente los mismos resultados que nuestra búsqueda original de `opsin` en [www.ncbi.nlm.nih.gov/genome/gdv](https://www.ncbi.nlm.nih.gov/genome/gdv/)? ¿Por qué?
> 
> > <solution-title></solution-title>
> > 
> > Los resultados pueden ser similares, pero sin duda hay algunas diferencias. De hecho, no sólo una búsqueda de texto es diferente a una búsqueda de secuencias en términos de método, sino que además en esta segunda ronda partimos de la secuencia de una opsina específica, es decir, de una rama de todo el árbol genealógico de proteínas. Algunos de los miembros de la familia son más similares entre sí, por lo que este tipo de búsqueda examina toda la familia desde una perspectiva bastante sesgada.
> {: .solution}
{: .question}

# Más información sobre nuestra proteína

Hasta ahora, hemos explorado esta información sobre las opsinas:
- cómo saber qué proteínas de un determinado tipo existen en un genoma,
- cómo saber dónde se encuentran a lo largo del genoma,
- cómo obtener más información sobre un gen de interés,
- cómo descargar sus secuencias en diferentes formatos,
- cómo utilizar estos archivos para realizar una búsqueda por similitud.

Puede que ahora sienta curiosidad por saber más sobre las proteínas que codifican. Ya hemos recopilado alguna información (por ejemplo, enfermedades asociadas), pero en los próximos pasos la cruzaremos con datos sobre la estructura de la proteína, localización, interactores, funciones, etc.

El portal a visitar para obtener toda la información sobre una proteína es [UniProt](https://www.uniprot.org/). Podemos buscar en él utilizando una búsqueda de texto, o el nombre del gen o de la proteína. Vamos a utilizar nuestra palabra clave habitual `OPN1LW`.

> <hands-on-title>Búsqueda en UniProt</hands-on-title>
> 
> 1. Abrir [UniProt](https://www.uniprot.org/)
> 2. Escriba `OPN1LW` en la barra de búsqueda
> 3. Seleccione la vista de tarjeta
{: .hands-on}

El primer resultado debe ser `P04000 · OPSR_HUMAN`. Antes de abrir la página, dos cosas a notar:

1. El nombre de la proteína `OPSR_HUMAN` es diferente del nombre del gen, así como sus IDs son.
2. Esta entrada tiene una estrella dorada, lo que significa que fue anotada y curada manualmente.

> <hands-on-title>Abrir un resultado en UniProt</hands-on-title>
> 
> 1. Haga clic en `P04000 · OPSR_HUMAN`
{: .hands-on}

![Captura de pantalla de la cabecera de la página de entrada de UniProt](./images/UniProt.png "UniProt page")

Esta es una página larga con mucha información, hemos diseñado un [tutorial entero]({% link topics/data-science/tutorials/online-resources-protein/tutorial.md %}) para repasarla.

