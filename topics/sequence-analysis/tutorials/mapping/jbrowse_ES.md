> <hands-on-title>Visualización de lecturas JBrowse</hands-on-title>
>
> 1. {% tool [JBrowse](toolshed.g2.bx.psu.edu/repos/iuc/jbrowse/jbrowse/1.16.11+galaxy1) %} con los siguientes parámetros:
>    - *"Reference genome to display"*: Use a built-in genome
>       - *"Select a reference genome"*: `mm10`
>    - *"JBrowse-in-Galaxy Action"*: `New JBrowse Instance`
>    - *"Insert Track Group"*:
>       - *"Insert Annotation Track"*:
>          - *"Track Type"*: BAM Pileups
>          - *"BAM Track Data"*: `aligned reads` (output de **Bowtie2** {% icon tool %})
>          - *"Autogenerate SNP Track"*: Yes
>          - *"Track Visibility"*: On for new users
> 2. Visualiza the dataset {% icon galaxy-eye %}
> 3. Zoom en el `{{ include.region_to_zoom }}`
{: .hands_on}

> <comment-title>Slow</comment-title>
>
> Esto puede tardar uno o dos minutos en ejecutarse, dependiendo de los recursos de tu instancia de entrenamiento. Toma tiempo porque el servidor construye un pequeño sitio web para ti y preprocesa el genoma de referencia en un formato más eficiente. Si quisieras compartir esto con tus colegas, podrías descargar este conjunto de datos y colocarlo directamente en tu propio servidor web.
{: .comment}

Las lecturas tienen una dirección: se alinean a la hebra directa o a la hebra inversa, respectivamente. Al hacer clic sobre una lectura, se muestra información adicional.

> <question-title></question-title>
>
> 1. ¿Qué significan la gota y la línea en el SNP autogenerado? 
> 2. ¿Qué significan las lecturas de diferentes colores?
>
> > <solution-title></solution-title>
> > 1. Si suficientes lecturas tienen un valor diferente, se marca con un ícono en forma de lágrima. La gráfica de cobertura se representa en altura según el porcentaje de lecturas con una llamada diferente en esa posición.
> > 2. Códigos de colores:
> >
> >    Colour                                                                         | Meaning
> >    ---                                                                            | ---
> >    <i style="background:#ec8b8b">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</i> Original red  | Forward strand
> >    <i style="background:#8f8fd8">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</i> Original blue | Reverse strand
> >    <i style="background:#d11919">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</i> Hard red      | Forward strand, missing mate
> >    <i style="background:#1919d1">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</i> Hard Blue     | Reverse strand, missing mate
> >    <i style="background:#ecc8c8">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</i> Light red     | Forward strand not proper
> >    <i style="background:#bebed8">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</i> Light blue    | Reverse strand, not proper
> >    <i style="background:#000000">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</i> Black         | Forward, diff chr
> >    <i style="background:#969696">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</i> Grey          | Reverse, diff chr
> >    <i style="background:#999999">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</i> Grey          | No strand
> {: .solution }
{: .question}
