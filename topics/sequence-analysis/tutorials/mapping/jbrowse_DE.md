> <hands-on-title>Visualisierung der Reads in JBrowse</hands-on-title>
> 
> 1. {% tool [JBrowse](toolshed.g2.bx.psu.edu/repos/iuc/jbrowse/jbrowse/1.16.11+galaxy1) %} browser, mit den folgenden Parametern:
>    - *"Referenzgenom zur Anzeige "*: Verwenden Sie ein eingebautes Genom
>       - *"Wählen Sie ein Referenzgenom "*: `mm10`
>    - *"JBrowse-in-Galaxy Action "*: `New JBrowse Instance`
>    - *"Insert Track Group "*:
>       - *"Insert Annotation Track "*:
>          - *"Track Type "*: BAM Pileups
>          - *"BAM Track Data "*: `aligned reads` (Ausgabe von **Bowtie2** {% icon tool %})
>          - *"Autogenerate SNP Track "*: Ja
>          - *"Track Visibility "*: Ein für neue Benutzer
> 2. Visualisieren Sie den Datensatz {% icon galaxy-eye %}
> 3. Zoom auf die `{{ include.region_to_zoom }}`
{: .hands_on}

> <comment-title>Langsam</comment-title>
> Die Ausführung kann je nach den Ressourcen Ihrer Trainingsinstanz ein bis zwei Minuten dauern. Der Server baut eine kleine Website für Sie auf und bereitet das Referenzgenom in einem effizienteren Format vor. Wenn Sie dies mit Ihren Kollegen teilen möchten, können Sie diesen Datensatz herunterladen und direkt auf Ihren Webserver stellen.
{: .comment}

Die Reads haben eine Richtung: Sie werden auf den Vorwärts- bzw. Rückwärtsstrang abgebildet. Wenn Sie auf einen Read klicken, werden zusätzliche Informationen angezeigt

> <question-title></question-title>
> 
> 1. Was bedeuten die Tropfenform und die Linie in der automatisch generierten SNP-Spur?
> 2. Was bedeuten unterschiedlich gefärbte Reads?
> 
> > <solution-title></solution-title>
> > 1. Wenn genügend Reads einen anderen Wert haben, wird dieser mit einem Tränensymbol markiert. Der Coverage Plot wird in der Höhe mit dem Prozentsatz der Reads mit einem anderen Call an dieser Position markiert
> > 2. Farbcodes:
> > 
> >    | Colour                                                | Meaning                      |
> >    | ----------------------------------------------------- | ---------------------------- |
> >    | <i style="background:#ec8b8b">     </i> Original red  | Forward strand               |
> >    | <i style="background:#8f8fd8">     </i> Original blue | Reverse strand               |
> >    | <i style="background:#d11919">     </i> Hard red      | Forward strand, missing mate |
> >    | <i style="background:#1919d1">     </i> Hard Blue     | Reverse strand, missing mate |
> >    | <i style="background:#ecc8c8">     </i> Light red     | Forward strand not proper    |
> >    | <i style="background:#bebed8">     </i> Light blue    | Reverse strand, not proper   |
> >    | <i style="background:#000000">     </i> Black         | Forward, diff chr            |
> >    | <i style="background:#969696">     </i> Grey          | Reverse, diff chr            |
> >    | <i style="background:#999999">     </i> Grey          | No strand                    |
> {: .solution }
{: .question}
