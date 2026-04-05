> <question-title></question-title>
>
> ¿Cuál es la relación entre `{{ include.forward }}` y `{{ include.reverse }}` ?
>
> > <solution-title></solution-title>
> > Los datos han sido secuenciados utilizando secuenciación paired-end (de extremos emparejados).
> > 
> > La secuenciación paired-end se basa en la idea de que los fragmentos iniciales de ADN (más largos que la longitud de lectura real) se secuencian desde ambos extremos. Este enfoque da lugar a dos lecturas por fragmento: la primera en orientación directa y la segunda en orientación inversa complementaria. La distancia entre ambas lecturas es conocida, por lo que puede utilizarse como información adicional para mejorar el mapeo de las lecturas.
> > 
> > Con la secuenciación paired-end, cada fragmento queda más cubierto que con la secuenciación single-end (solo se secuencia en orientación directa):
> >
> > ```
> >     ------>
> >       ------>
> >         ------>
> >           ------>
> >
> >     ----------------------------- [fragment]
> >
> >     ------>         <------
> >       ------>         <------
> >         ------>         <------
> >           ------>         <------
> > ```
> > 
> > La secuenciación paired-end genera entonces dos ficheros:
> > - Un fichero con la secuenciación correspondiente a la orientación hacia adelante de todos los fragmentos
> > - Un fichero con la secuenciación correspondiente a la orientación inversa de todos los fragmentos
> >
> > Aquí, `{{ include.forward }}` corresponde a las lecturas hacia adelante y `{{ include.reverse }}` a las inversas.
> {: .solution }
{: .question}
