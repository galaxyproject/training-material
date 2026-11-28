> <question-title></question-title>
> 
> Qual è la relazione tra `{{ include.forward }}` e `{{ include.reverse }}`?
> 
> > <solution-title></solution-title>
> > I dati sono stati sequenziati utilizzando il sequenziamento paired-end.
> > 
> > Il sequenziamento paired-end si basa sull'idea che i frammenti di DNA iniziali (più lunghi della lunghezza effettiva della lettura) vengano sequenziati da entrambi i lati. Questo approccio dà luogo a due letture per frammento, con la prima lettura in orientamento in avanti e la seconda in orientamento inverso-complementare. La distanza tra le due letture è nota. Pertanto, può essere utilizzata come informazione aggiuntiva per migliorare la mappatura delle letture.
> > 
> > Con il sequenziamento paired-end, ogni frammento è più coperto che con il sequenziamento single-end (viene sequenziato solo l'orientamento forward):
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
> > Il sequenziamento paired-end genera quindi 2 file:
> > - Un file con le sequenze corrispondenti all'orientamento in avanti di tutti i frammenti
> > - Un file con le sequenze corrispondenti all'orientamento inverso di tutti i frammenti
> > 
> > Qui `{{ include.forward }}` corrisponde alle letture in avanti e `{{ include.reverse }}` alle letture inverse.
> {: .solution }
{: .question}
