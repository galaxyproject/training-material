> <question-title></question-title>
>
> Perché lo strumento di trimming viene eseguito solo una volta su un set di dati paired-end e non due volte, una per ogni set di dati?
>
> > <solution-title></solution-title>
> >
> > Lo strumento può rimuovere le sequenze che diventano troppo corte durante il processo di trimming. Per i file paired-end rimuove intere coppie di sequenze se una (o entrambe) delle due letture diventa più corta del limite di lunghezza impostato. Le letture di una coppia di letture che sono più lunghe di una determinata soglia, ma per le quali la lettura partner è diventata troppo corta, possono essere facoltativamente scritte in file single-end. In questo modo si garantisce che le informazioni di una coppia di letture non vadano perse del tutto se solo una lettura è di buona qualità.
> >
> {: .solution}
{: .question}
