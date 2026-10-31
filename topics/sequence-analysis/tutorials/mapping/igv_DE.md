> <hands-on-title>Visualisierung der Reads in IGV</hands-on-title>
> 
> Es gibt zwei Möglichkeiten, es auszuführen:
> 
> 1. Wenn Sie IGV installiert haben (oder es installieren wollen):
>     1. Installieren Sie [IGV](https://software.broadinstitute.org/software/igv/download)
>     2. Starten Sie IGV auf Ihrem Computer
>     3. Erweitern Sie die {% icon param-file %} Ausgabe von **{{ include.tool }}** {% icon tool %}
>     4. Klicken Sie auf das `local` in `display with IGV`, um die Reads in den IGV-Browser zu laden
> 2. Wenn Sie nicht über IGV verfügen
>     1. Klicken Sie auf `Mouse mm10` (oder den richtigen Organismus) in `display with IGV`, um die Reads in den IGV-Browser zu laden
> 4. Zoom auf die `{{ include.region_to_zoom }}`
{: .hands_on}

Die Reads haben eine Richtung: Sie werden auf den Vorwärts- bzw. Rückwärtsstrang abgebildet. Wenn Sie den Mauszeiger über einen Read bewegen, werden zusätzliche Informationen angezeigt

> <question-title></question-title>
> 
> 1. Was könnte es bedeuten, wenn ein Balken in der Abdeckungsansicht farbig ist?
> 2. Was könnte der Grund dafür sein, dass eine Anzeige weiß statt grau ist?
> 
> > <solution-title></solution-title>
> > 1. Wenn ein Nukleotid in mehr als 20% der qualitätsgewichteten Reads von der Referenzsequenz abweicht, färbt IGV den Balken im Verhältnis zur Anzahl der Reads jeder Base.
> > 2. Sie haben eine Abbildungsqualität gleich Null. Die Interpretation dieser Mapping-Qualität hängt vom Mapping-Aligner ab, da einige häufig verwendete Aligner diese Konvention verwenden, um einen Read mit mehreren Alignments zu markieren. In einem solchen Fall wird der Read auch auf eine andere Stelle mit gleich guter Platzierung abgebildet. Es ist auch möglich, dass der Read nicht eindeutig platziert werden kann, aber die anderen Platzierungen nicht unbedingt gleichwertige Treffer ergeben.
> {: .solution }
{: .question}

> <comment-title>Tipps für IGV</comment-title>
> 1. Da die Anzahl der Reads in einer Region recht groß sein kann, zeigt der IGV-Browser standardmäßig nur die Reads an, die in ein kleines Fenster fallen. Dieses Verhalten kann im IGV von `view > Preferences > Alignments` geändert werden.
> 2. Wenn das Genom, das Sie interessiert, dort nicht zu finden ist, überprüfen Sie, ob es über **More...** verfügbar ist. Wenn dies nicht der Fall ist, können Sie es manuell über das Menü **Genome** -> **Genom laden von...** hinzufügen
> 
>    ![Genom in IGV auswählen](../../images/igv_select_genome.png "Genom in IGV auswählen")
> 
> Eine allgemeine Beschreibung der Benutzeroberfläche des IGV-Browsers finden Sie hier: [IGV Browser Beschreibung]({% link topics/introduction/tutorials/igv-introduction/tutorial.md %})
{: .comment}


