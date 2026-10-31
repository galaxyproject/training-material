> <question-title></question-title>
> 
> Welche Beziehung besteht zwischen `{{ include.forward }}` und `{{ include.reverse }}`?
> 
> > <solution-title></solution-title>
> > Die Daten wurden mittels Paired-End-Sequenzierung sequenziert.
> > 
> > Die Paired-End-Sequenzierung basiert auf der Idee, dass die anfänglichen DNA-Fragmente (länger als die tatsächliche Leselänge) von beiden Seiten sequenziert werden. Dieser Ansatz führt zu zwei Reads pro Fragment, wobei der erste Read in Vorwärtsorientierung und der zweite Read in Rückwärtskomplementorientierung erfolgt. Der Abstand zwischen beiden Reads ist bekannt. Er kann daher als zusätzliche Information verwendet werden, um das Read-Mapping zu verbessern.
> > 
> > Bei der Paired-End-Sequenzierung ist jedes Fragment stärker abgedeckt als bei der Single-End-Sequenzierung (nur die Vorwärtsrichtung wird sequenziert):
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
> > Die Paired-End-Sequenzierung erzeugt dann 2 Dateien:
> > - Eine Datei mit den Sequenzen, die der Vorwärtsorientierung aller Fragmente entsprechen
> > - Eine Datei mit den Sequenzen, die der umgekehrten Orientierung aller Fragmente entsprechen
> > 
> > Hier entspricht `{{ include.forward }}` den Forward Reads und `{{ include.reverse }}` den Reverse Reads.
> {: .solution }
{: .question}
