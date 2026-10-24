> <hands-on-title>Visualizzazione delle letture in JBrowse</hands-on-title>
> 
> 1. {% tool [JBrowse](toolshed.g2.bx.psu.edu/repos/iuc/jbrowse/jbrowse/1.16.11+galaxy1) %} browser, con i seguenti parametri:
>    - *"Genoma di riferimento da visualizzare "*: Utilizzare un genoma incorporato
>       - *"Seleziona un genoma di riferimento "*: `mm10`
>    - *"Azione JBrowse-in-Galaxy "*: `New JBrowse Instance`
>    - *"Inserisci gruppo di tracce "*:
>       - *"Inserisci traccia annotazione "*:
>          - *"Tipo di traccia "*: BAM Pileup
>          - *"BAM Track Data "*: `aligned reads` (output di **Bowtie2** {% icon tool %})
>          - *"Autogenerazione traccia SNP "*: Sì
>          - *"Visibilità traccia "*: Attiva per i nuovi utenti
> 2. Visualizzare il dataset {% icon galaxy-eye %}
> 3. Zoom sul `{{ include.region_to_zoom }}`
> 
{: .hands_on}

> <comment-title>Slow</comment-title>
> L'esecuzione può richiedere uno o due minuti, a seconda delle risorse dell'istanza di formazione. Ci vuole tempo perché il server crea un piccolo sito web per voi e preelabora il genoma di riferimento in un formato più efficiente. Se si desidera condividere questo lavoro con i colleghi, è possibile scaricare questo set di dati e inserirlo direttamente nel proprio server web.
{: .comment}

Le letture hanno una direzione: sono mappate rispettivamente sul filamento avanti o sul filamento indietro. Quando si fa clic su una lettura, vengono visualizzate informazioni aggiuntive

> <question-title></question-title>
> 
> 1. Cosa significano la forma a goccia e la linea nella traccia SNP autogenerata?
> 2. Cosa significano letture di colore diverso?
> 
> > <solution-title></solution-title>
> > 1. Se un numero sufficiente di letture ha un valore diverso, viene contrassegnato con un'icona a goccia. Il grafico della copertura è contrassegnato in altezza con la percentuale di letture con una chiamata diversa in quella posizione
> > 2. Codici colore:
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
