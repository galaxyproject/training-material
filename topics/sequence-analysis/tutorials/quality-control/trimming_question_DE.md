> <question-title></question-title>
>
> Warum führen wir das Trimming-Tool nur einmal auf einem Paired-End-Datensatz aus und nicht zweimal, einmal für jeden Datensatz?
>
> > <solution-title></solution-title>
> >
> > Das Tool kann Sequenzen entfernen, wenn sie während des Trimming-Prozesses zu kurz werden. Bei Paired-End-Dateien entfernt es ganze Sequenzpaare, wenn einer (oder beide) der beiden Reads kürzer als der eingestellte Längen-Cutoff geworden sind. Reads eines Read-Paares, die länger als ein bestimmter Schwellenwert sind, bei denen aber der Partner-Read zu kurz geworden ist, können optional in Single-End-Dateien ausgeschrieben werden. Dadurch wird sichergestellt, dass die Informationen eines Lesepaares nicht vollständig verloren gehen, wenn nur ein Read von guter Qualität ist.
> {: .solution}
{: .question}
