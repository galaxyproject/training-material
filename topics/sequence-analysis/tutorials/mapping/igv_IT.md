> <hands-on-title>Visualizzazione delle letture in IGV</hands-on-title>
> 
> Ci sono due modi per eseguirlo:
> 
> 1. Se avete installato IGV (o volete installarlo):
>     1. Installare [IGV](https://software.broadinstitute.org/software/igv/download)
>     2. Avvio di IGV sul computer
>     3. Espande il {% icon param-file %} output di **{{ include.tool }}** {% icon tool %}
>     4. Fare clic sul `local` in `display with IGV` per caricare le letture nel browser IGV
> 2. Se non si dispone di IGV
>     1. Fare clic su `Mouse mm10` (o sull'organismo corretto) in `display with IGV` per caricare le letture nel browser IGV
> 4. Zoom sul `{{ include.region_to_zoom }}`
> 
{: .hands_on}

Le letture hanno una direzione: sono mappate, rispettivamente, sul filamento avanti o sul filamento indietro. Quando si passa il mouse su una lettura, vengono visualizzate informazioni supplementari

> <question-title></question-title>
> 
> 1. Cosa può significare se una barra nella vista di copertura è colorata?
> 2. Quale potrebbe essere il motivo per cui una lettura è bianca invece che grigia?
> 
> > <solution-title></solution-title>
> > 1. Se un nucleotide differisce dalla sequenza di riferimento in più del 20% delle letture ponderate per la qualità, IGV colora la barra in proporzione al numero di letture di ciascuna base.
> > 2. Hanno una qualità di mappatura pari a zero. L'interpretazione di questa qualità di mappatura dipende dall'allineatore di mappatura, poiché alcuni allineatori comunemente utilizzati utilizzano questa convenzione per contrassegnare una lettura con allineamenti multipli. In questo caso, la lettura si mappa anche in un'altra posizione con un posizionamento altrettanto buono. È anche possibile che la lettura non sia posizionata in modo univoco, ma che gli altri posizionamenti non diano necessariamente risultati altrettanto buoni.
> {: .solution }
{: .question}

> <comment-title>Consigli per IGV</comment-title>
> 1. Poiché il numero di letture in una regione può essere molto grande, il browser IGV visualizza per impostazione predefinita solo le letture che rientrano in una piccola finestra. Questo comportamento può essere modificato in IGV da `view > Preferences > Alignments`.
> 2. Se il genoma di vostro interesse non è presente, controllate se è disponibile tramite **More...**. In caso contrario, è possibile aggiungerlo manualmente tramite il menu **Genomi** -> **Carica genoma da...**
> 
>    ![Seleziona genoma in IGV](../../images/igv_select_genome.png "Seleziona genoma in IGV")
> 
> Una descrizione generale dell'interfaccia utente del browser IGV è disponibile qui: [Descrizione del browser IGV]({% link topics/introduction/tutorials/igv-introduction/tutorial.md %})
{: .comment}
