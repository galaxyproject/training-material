> <question-title></question-title>
>
> Por qué ejecutamos la herramienta de recorte solo una vez en un conjunto de datos paired-end y no dos veces, una para cada conjunto de datos?
>
> > <solution-title></solution-title>
> >
> > La herramienta puede eliminar secuencias si se vuelven demasiado cortas durante el proceso de recorte. En archivos paired-end, elimina pares de secuencias completos si una (o ambas) de las dos lecturas quedan por debajo del umbral de longitud establecido. Las lecturas de un par que superan dicho umbral, pero cuya pareja ha quedado demasiado corta, pueden opcionalmente guardarse en archivos single-end. Esto garantiza que no se pierda por completo la información de un par de lecturas si solo una de ellas es de buena calidad.
> >
> {: .solution}
{: .question}
