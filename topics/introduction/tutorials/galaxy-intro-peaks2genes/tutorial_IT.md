---
layout: tutorial_hands_on
title: Dai picchi ai geni
zenodo_link: https://doi.org/10.5281/zenodo.1025586
level: Introductory
questions:
- Come usare Galaxy?
- Come derivare una lista di identificativi di geni da regioni *peaks*
objectives:
- Familiarizza con le basi di Galaxy
- Impara come ottenere dati da fonti esterno
- Impara come lanciare dei *tool*
- Impara come funzionano le cronologie
- Impara come creare dei flussi di lavoro
- Impara come condividere il tuo lavoro
time_estimation: 3H
key_points:
- Galaxy offre un'interfaccia utente grafica facile da usare per strumenti da riga di comando spesso complessi
- Galaxy conserva una cronologia completa delle analisi
- I flussi di lavoro consentono di ripetere l'analisi su dati diversi
- Galaxy può connettersi a fonti esterne per l'importazione e la visualizzazione dei dati
- Galaxy offre modalità per condividere risultati e metodi con altri
subtopic: next-steps
contributions:
  authorship:
  - pajanne
  - blankclemens
  - bebatut
  - bgruening
  - nsoranzo
  - dyusuf
  - sarah-peter
  - hexylena
  editing:
  - teresa-m
  - dadrasarmin
  translation:
  - lisanna
  - silviadg87
  - unode
  funding:
  - elixir-europe
  - deNBI
  - uni-freiburg
  - biont
lang: it
tags:
  - deutsch
  - english
  - español
translations:
  - de
  - en
  - es
---


Ci siamo imbattuti in un articolo ({% cite Li2012 %}) intitolato *"L'istone
acetiltransferasi MOF è un regolatore chiave della rete trascrizionale centrale delle
cellule staminali embrionali"*. L'articolo contiene l'analisi dei possibili geni
bersaglio di un'interessante proteina chiamata Mof. I target sono stati ottenuti
mediante ChIP-seq nei topi e i dati grezzi sono disponibili su [GEO]
(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37268). Tuttavia, l'elenco dei
geni non si trova né nei material supplementari dell'articolo, né sono parte dell'invio a GEO. La cosa più
simile che abbiamo trovato è un file in GEO contenente un elenco delle regioni in cui il
segnale è significativamente arricchito (i cosiddetti *peaks*):

1 | 3660676 | 3661050 | 375 | 210 | 62.0876250438913 | -2.00329386666667
1 | 3661326 | 3661500 | 175 | 102 | 28.2950833625942 | -0.695557142857143
1 | 3661976 | 3662325 | 350 | 275 | 48.3062708406486 | -1.29391285714286
1 | 3984926 | 3985075 | 150 | 93 | 34.1879823073944 | -0.816992
1 | 4424801 | 4424900 | 100 | 70 | 26.8023246007435 | -0.66282

**Tabella 1** Sottocampione del file disponibile

L'obiettivo di questo esercizio è di **trasformare questo elenco di regioni genomiche in
un elenco di possibili geni bersaglio**.

{% snippet faqs/galaxy-it/analysis_results_may_vary.md %}

> <agenda-title></agenda-title>
> 
> In questo tutorial, ci occuperemo di:
> 
> 1. TOC
> {:toc}
{: .agenda}

# Pretrattamenti

> <hands-on-title>Apri Galaxy</hands-on-title>
> 
> 1. Naviga verso un'istanza Galaxy: quella raccomandata dal tuo istruttore o una
>    nell'elenco **Istanze Galaxy** all'inizio di questa pagina
> 2. Accesso o registrazione (pannello superiore)
> 
>    ![Accedi o registrati nel pannello superiore](../../images/login_register.png)
{: .hands_on}

L'interfaccia di Galaxy è composta da tre parti principali. Gli strumenti disponibili
sono elencati a sinistra, la cronologia delle analisi è registrata a destra e il
pannello centrale mostra gli strumenti e i set di dati.

![schermata dell'interfaccia Galaxy che mostra il pannello della cronologia a destra, il pannello degli strumenti a sinistra e il pannello principale al centro](../../images/galaxy_interface.png "L'interfaccia Galaxy")

Cominciamo con una nuova cronologia.

> <hands-on-title>Creare la cronologia</hands-on-title>
> 
> 1. Assicurarsi di avere una cronologia di analisi vuota.
> 
>    {% snippet faqs/galaxy-it/histories_create_new.md %}
> 
> 2. **Rinomina la cronologia** per facilitarne il riconoscimento
> 
>    > <tip-title>Rinominare una cronologia</tip-title>
>    > 
>    > * Cliccare sul titolo della cronologia (per impostazione predefinita il titolo è
>    >   `Unnamed history`)
>    > 
>    >   ![Rinominare la cronologia]({% link shared/images/rename_history.png %})
>    > 
>    > * Digitare `Galaxy Introduction` come nome
>    > * Premere <kbd>Invio</kbd>
>    {: .tip}
>
{: .hands_on}

## Caricamento dei dati

> <hands-on-title>Caricamento dati</hands-on-title>
> 
> 1. scaricare l'elenco delle regioni di picco (il file [`GSE37268_mof3.out.hpeak.txt.gz`](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE37268&format=file&file=GSE37268%5Fmof3%2Eout%2Ehpeak%2Etxt%2Egz)) da [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37268) sul computer
> 2. Fare clic sul pulsante di caricamento in alto a sinistra dell'interfaccia
> 
>    ![Icona caricamento dati](../../../galaxy-interface/images/upload_button.png)
> 
> 3. Premere **Scegli file locali** e cercare il file sul computer
> 4. Selezionare `interval` come **Tipo**
> 5. premete **Avvio**
> 6. Premere **Chiudo**
> 7. Attendere il completamento del caricamento. Galaxy de-compatterà automaticamente il
>    file.
> 
> 8. Dopo questa operazione si vedrà il primo elemento della cronologia nel pannello
>    destro di Galaxy. Passerà attraverso gli stati grigio (preparazione/in attesa) e
>    giallo (esecuzione) per poi diventare verde (successo):
> 
>    ![Sezione storia](../../images/intro_01.png)
> 
>    Il caricamento diretto dei file non è l'unico modo per inserire i dati in Galaxy
> 
>    {% snippet faqs/galaxy-it/datasets_import_via_link.md format="interval" %}
> 
>    > <tip-title>Importazione di dati in Galaxy</tip-title> Ci sono [ulteriori opzioni]({% link topics/galaxy-interface/tutorials/get-data/slides.html %}) per gli utenti avanzati.
>    {: .tip}
>
{: .hands_on}


> <comment-title>Formato file Intervallo</comment-title>
> Il formato **Intervallo** è un formato Galaxy per rappresentare intervalli genomici. È separato da tabulazioni, ma ha il requisito aggiuntivo che tre delle colonne devono essere:
> - ID cromosoma
> - posizione iniziale (in base 0)
> - posizione finale (fine esclusiva)
> 
> è possibile specificare anche una colonna "filamento' opzionale e utilizzare una riga di
> intestazione iniziale per etichettare le colonne, che non devono essere in un ordine
> particolare. A differenza del formato BED (vedi sotto), possono essere presenti anche
> colonne aggiuntive arbitrarie.
> 
> Per ulteriori informazioni sui formati utilizzabili in Galaxy, consultare la pagina [Galaxy Data Formats](https://usegalaxy.org/static/formatHelp.html).
{: .comment}


> <hands-on-title>Controllare e modificare gli attributi di un file</hands-on-title>
> 
> 1. fare clic sul file nel pannello della cronologia
> 
>    vengono visualizzate alcune meta-informazioni (ad esempio, formato, database di
>    riferimento) sul file e l'intestazione del file, oltre al numero di righe del file
>    (48.647):
> 
>    ![File espanso nella cronologia](../../images/input_expanded_file.png)
> 
> 2. Fare clic sull'icona {% icon galaxy-eye %} (occhio) (**Visualizza dati**) nel set
>    di dati nella cronologia
> 
>    il contenuto del file è visualizzato nel pannello centrale
> 
> 3. Fare clic sull'icona {% icon galaxy-pencil %} (matita) (**Modifica attributi**)
>    nel vostro set di dati nella cronologia
> 
>    Nel pannello centrale viene visualizzato un modulo per modificare gli attributi del
>    set di dati
> 
> 4. Cerca `mm9` nell'attributo **Database/Build** e seleziona `Mouse July 2007
>    (NCBI37/mm9)` (la carta ci dice che i picchi sono da `mm9`)
> 
>    ![La versione del database può essere selezionata da un menu a discesa. Gli utenti possono iniziare a digitare il nome del database per filtrare l'elenco](../../images/Search-for-mm9.PNG)
> 
> 5. Cliccare su **Salva** in alto
> 6. Aggiungere un tag chiamato `#peaks` al set di dati per renderlo più facilmente
>    rintracciabile nella cronologia
> 
>    {% snippet faqs/galaxy-it/datasets_add_tag.md %}
> 
>    Il set di dati dovrebbe ora apparire come segue nella cronologia
> 
>    ![File picchi](../../images/input_tagged_file.png){: width="250px" height="300px"}
> 
{: .hands_on}

Per trovare i geni correlati a queste regioni di picco, abbiamo bisogno anche di un
elenco di geni nei topi, che possiamo ottenere dall'UCSC.

> <hands-on-title>Caricamento dati da UCSC</hands-on-title>
> 
> 1. Cercare `UCSC Main` nella barra di ricerca dello strumento (in alto a sinistra)
> 
>    ![strumento principale UCSC nella sezione strumenti](../../images/101_01.png)
> 
> 2. Cliccare su `UCSC Main` {% icon tool %}
> 
>    Verrà visualizzato il browser delle tabelle **UCSC**, che ha un aspetto simile a
>    questo:
> 
>    ![interfaccia browser tabella UCSC](../../images/intro_02.png)
> 
> 3. Impostare le seguenti opzioni:
>     - *"clade "*: `Mammal`
>     - *"genoma "*: `Mouse`
>     - *"assemblaggio "*: `July 2007 (NCBI37/mm9)`
>     - *"gruppo "*: `Genes and Gene Predictions`
>     - *"traccia "*: `RefSeq Genes`
>     - *"tabella "*: `refGene`
>     - *"regione "*: `genome`
>     - *"formato di uscita "*: `BED - browser extensible data`
>     - *"Invia l'output a "*: `Galaxy` (solo)
> 
> 4. Fare clic sul pulsante **ottenere l'output**
> 
>    Verrà visualizzata la schermata successiva:
> 
>    ![Impostazioni di output](../../images/intro_03.png)
> 
> 5. Assicurarsi che *"Crea un record BED per "* sia impostato su `Whole Gene`
> 6. Fare clic sul pulsante **Invia query a Galaxy**
> 7. Attendere che il caricamento sia terminato
> 8. Rinominare il nostro set di dati in qualcosa di più riconoscibile come `Genes`
> 
>    {% snippet faqs/galaxy-it/datasets_rename.md name="Genes" %}
> 
> 9. Aggiungere un tag chiamato `#genes` al set di dati per renderlo più facilmente rintracciabile nella cronologia
> 
{: .hands_on}

> <comment-title>Formato file BED</comment-title>
> Il formato **BED - Browser Extensible Data** fornisce un modo flessibile per codificare le regioni geniche. Le linee BED hanno tre campi obbligatori:
> - ID cromosoma
> - posizione iniziale (in base 0)
> - posizione finale (fine esclusiva)
> 
> Possono esserci fino a nove campi opzionali aggiuntivi, ma il numero di campi per riga
> deve essere coerente in ogni singolo set di dati.
> 
> è possibile trovare maggiori informazioni al riguardo su
> [UCSC](https://genome.ucsc.edu/FAQ/FAQformat#format1), compresa una descrizione dei
> campi opzionali.
{: .comment}

Ora abbiamo raccolto tutti i dati necessari per iniziare la nostra analisi.

# Parte 1: approccio ingenuo

Per prima cosa utilizzeremo un approccio "ingenuo" per cercare di identificare i geni a
cui sono associate le regioni di picco. Identificheremo i geni che si sovrappongono per
almeno 1bp alle regioni di picco.

## Preparazione del file

Diamo un'occhiata ai nostri file per vedere cosa abbiamo qui.

> <hands-on-title>Visualizza il contenuto del file</hands-on-title>
> 
> 1. Fare clic sull'icona {% icon galaxy-eye %} (occhio) del file di picco per
>    visualizzarne il contenuto (occhio) (**Visualizza dati**) del file di picco per
>    visualizzarne il contenuto
> 
>    Dovrebbe essere così:
> 
>    ![Contenuto del file di picco](../../images/intro_04.png)
> 
> 2. Visualizza il contenuto delle regioni dei geni da UCSC
> 
>    ![Contenuto del file UCSC](../../images/intro_05.png)
> 
{: .hands_on}

> <question-title></question-title>
> 
> Mentre il file dell'UCSC ha etichette per le colonne, il file del picco non le ha.
> Riuscite a indovinare il significato delle colonne?
> 
> 
> > <solution-title></solution-title>
> > 
> > Questo file di picco non ha un formato standard e, solo guardandolo, non è possibile
> > scoprire il significato dei numeri nelle diverse colonne. Nel documento gli autori
> > affermano di aver utilizzato il tool
> > [HPeak](https://www.ncbi.nlm.nih.gov/pubmed/20598134).
> > 
> > consultando il manuale di HPeak possiamo scoprire che le colonne contengono le
> > seguenti informazioni:
> > 
> >  - nome del cromosoma in base al numero
> >  - coordinata iniziale
> >  - coordinata finale
> >  - lunghezza
> >  - posizione all'interno del picco con la più alta copertura di frammenti di DNA >    ipotetico (vertice)
> >  - non rilevante
> >  - non rilevante
> {: .solution}
{: .question}

Per confrontare i due file, dobbiamo assicurarci che i nomi dei cromosomi seguano lo
stesso formato. Come si può vedere, nel file di picco manca `chr` prima di qualsiasi
numero di cromosoma. Ma cosa succede con i cromosomi 20 e 21? Saranno invece X e Y?
Controlliamo:

> <hands-on-title>Vedi la fine del file</hands-on-title>
> 
> 1. Cercare lo strumento {% tool [Select last lines from a dataset (tail)](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_tail_tool/9.3+galaxy1) %} ed eseguirlo con le seguenti impostazioni:
>     - *"File di testo "*: il nostro file di picco `GSE37268_mof3.out.hpeak.txt.gz`
>     - *"Operazione "*: `Keep last lines`
>     - *"Numero di righe "*: Scegliere un valore, ad es. `100`
> 2. fare clic su **Strumento di esecuzione**
> 3. Attendere che il lavoro sia terminato
> 4. Ispezionare il file attraverso l'icona {% icon galaxy-eye %} (occhio) icona (**Visualizza dati**)
> 
>    > <question-title></question-title>
>    > 
>    > 1. Come si chiamano i cromosomi?
>    > 2. Come si chiamano i cromosomi X e Y?
>    > 
>    > > <solution-title></solution-title>
>    > > 1. I cromosomi sono dati solo dal loro numero. Nel file dei geni dell'UCSC,
>    > >    iniziavano con `chr`
>    > > 2. i cromosomi X e Y sono denominati 20 e 21
>    > {: .solution }
>    {: .question}
{: .hands_on}

Per convertire i nomi dei cromosomi abbiamo quindi due cose da fare:

1. aggiungere `chr`
2. cambia 20 e 21 in X e Y

> <hands-on-title>Adegua i nomi dei cromosomi</hands-on-title>
> 
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_column/1.1.3) %} in una colonna specifica con le seguenti impostazioni:
>     - *"File da elaborare "*: il nostro file di picco `GSE37268_mof3.out.hpeak.txt.gz`
>     - *"in colonna "*: `1`
>     - *"Trova schema "*: `[0-9]+`
> 
>       Questo cercherà le cifre numeriche
> 
>     - *"Sostituisci con "*: `chr&`
> 
>       `&` è un segnaposto per il risultato della ricerca del modello
> 
> 2. Rinominare il file di output `chr prefix added`.
> 
> 3. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_column/1.1.3) %} : Eseguiamo nuovamente il tool con altre due sostituzioni
>    - *"File da elaborare "*: l'output dell'ultima esecuzione, `chr prefix added`
>    - *"in colonna "*: `1`
>    - {% icon param-repeat %} Sostituzione
>      - *"Trova schema "*: `chr20`
>      - *"Sostituisci con "*: `chrX`
>    - {% icon param-repeat %} Inserisci Sostituzione
>      - *"Trova schema "*: `chr21`
>      - *"Sostituisci con "*: `chrY`
> 
>    > <tip-title>Eseguire uno strumento</tip-title>
>    > 
>    > * Espandere le informazioni sul set di dati
>    > * Premere l'icona {% icon galaxy-refresh %} (**Eseguire di nuovo questo lavoro**)
>    {: .tip}
> 
> 4. ispezionare il file più recente attraverso l'icona {% icon galaxy-eye %} (occhio).
>    Abbiamo avuto successo?
> 
>    Ora abbiamo molti file e dobbiamo fare attenzione a selezionare quelli corretti a
>    ogni passo.
> 
>    > <question-title></question-title>
>    > 
>    > Quante regioni ci sono nel nostro file di output? È possibile fare clic sul nome
>    > dell'output per espanderlo e vedere il numero.
>    > 
>    > > <solution-title></solution-title> Dovrebbe essere uguale al numero di regioni
>    > > nel vostro primo file, `GSE37268_mof3.out.hpeak.txt.gz`: 48.647 Se il vostro
>    > > dice 100 regioni, allora avete eseguito il programma sul file `Tail` e dovete
>    > > eseguire nuovamente i passaggi.
>    > {: .solution }
>    {: .question}
> 
> 5. Rinomina il file in qualcosa di più riconoscibile, ad esempio `Peak regions`
{: .hands_on}

## Analisi

Il nostro obiettivo è confrontare i due file di regione (il file dei geni e il file dei
picchi/*peaks*) per sapere quali picchi sono correlati a quali geni. Se si
vuole sapere solo quali picchi si trovano **all'interno** dei geni (all'interno del
corpo del gene) si può saltare il passaggio successivo. Altrimenti, potrebbe essere
ragionevole includere la regione **promoter** dei geni nel confronto, ad esempio perché
si vogliono includere i fattori di trascrizione negli esperimenti ChIP-seq. Non esiste
una definizione rigorosa di regione promotrice, ma comunemente si utilizzano 2kb a monte
del TSS (inizio della regione). Utilizzeremo lo strumento **Get Flanks** per ottenere le
regioni da 2kb basi a monte dell'inizio del gene a 10kb basi a valle dell'inizio (12kb
di lunghezza). Per fare ciò, diciamo allo strumento Get Flanks che vogliamo regioni a
monte dell'inizio, con un offset di 10kb, che siano lunghe 12kb, come mostrato nel
diagramma seguente.

![Ottieni fianchi](../../images/intro_get_flanks.png)

> <hands-on-title>Aggiungi la regione del promotore ai record del gene</hands-on-title>
> 
> 1. {% tool [Get Flanks](toolshed.g2.bx.psu.edu/repos/devteam/get_flanks/get_flanks1/1.0.0) %} restituisce le regioni di affiancamento per ogni gene, con le seguenti impostazioni:
>     - *"Seleziona dati "*: file `Genes` da UCSC
>     - *"Regione "*: `Around Start`
>     - *"Posizione della/e regione/i di affiancamento "*: `Upstream`
>     - *"Offset "*: `10000`
>     - *"Lunghezza della/e regione/i affiancata/e "*: `12000`
> 
>    Questo strumento restituisce le regioni di affiancamento per ogni gene
> 
> 2. confronta le righe del file BED risultante con l'input per scoprire come sono
>    cambiate le posizioni di inizio e fine
> 
>    > <tip-title>Ispezione di diversi file usando lo scratchbook</tip-title>
>    > 
>    > * Fare clic su **Abilitazione/disabilitazione di Scratchbook** nel pannello
>    >   superiore
>    > 
>    >   ![Abilita/Disabilita Scratchbook](../../images/intro_scratchbook_enable.png)
>    > 
>    > * Fare clic sull'icona {% icon galaxy-eye %} (occhio) dei file da ispezionare
>    > * Cliccare su **Mostra/Nascondi quaderno**
>    > 
>    >   ![Mostra/Nascondi libro dei graffi](../../images/intro_scratchbook_show_hide.png)
>    {: .tip}
> 
> 3. Rinominare il set di dati per riflettere i risultati ottenuti (`Promoter regions`)
{: .hands_on}

l'output è costituito da regioni che partono da 2kb a monte del TSS e includono 10kb a
valle. Per le regioni di input sul filamento positivo, ad esempio `chr1 134212701
134230065`, si ottiene `chr1 134210701 134222701`. Per le regioni sul filamento
negativo, ad esempio `chr1 8349819 9289958`, si ottiene `chr1 9279958 9291958`.

Si sarà notato che il file UCSC è in formato `BED` e ha un database associato. Questo è
ciò che vogliamo anche per il nostro file di picco. Lo strumento **Intersect** che
utilizzeremo è in grado di convertire automaticamente i file di intervallo in formato
BED, ma convertiremo il nostro file di intervallo esplicitamente qui per mostrare come
si può ottenere questo risultato con Galaxy.

> <hands-on-title>Cambia formato e database</hands-on-title>
> 
> 1. Fare clic sull'icona {% icon galaxy-pencil %} (matita) nella voce della cronologia del file della regione di picco
> 2. passa alla scheda **Datatype**
> 3. Nella sezione **Convert to Datatype** sotto *"Target datatype "* selezionare: `bed (using 'Convert Genomic Interval To Bed')`
> 4. Premere **Crea set di dati**
> 5. Verificare che "Database/Build" sia `mm9` (la build del database per i topi utilizzata nel documento)
> 6. Rinomina il file in qualcosa di più riconoscibile, ad esempio `Peak regions BED`
{: .hands_on}

È il momento di trovare gli intervalli di sovrapposizione (finalmente!). Per farlo,
vogliamo estrarre i geni che si sovrappongono/intersecano con i nostri picchi.

> <hands-on-title>Trova sovrapposizioni</hands-on-title>
> 
> 1. {% tool [Intersect](toolshed.g2.bx.psu.edu/repos/devteam/intersect/gops_intersect_1/1.0.0) %} gli intervalli di due set di dati, con le seguenti impostazioni:
>     - *"Ritorno "*: `Overlapping Intervals`
>     - *"di "*: il file UCSC con le regioni dei promotori (`Promoter regions`)
>     - *"che intersecano "*: il nostro file della regione di picco da **Replace** (`Peak regions BED`)
>     - *"per almeno "*: `1`
> 
>    > <comment-title></comment-title>
>    > L'ordine degli input è importante! Vogliamo ottenere un elenco di **geni**, quindi il set di dati corrispondente con le informazioni sui geni deve essere il primo input (`Promoter regions`).
>    {: .comment}
>    ![Picchi di sovrapposizione dei geni](../../images/intro_overlapping_genes.png)
{: .hands_on}

Ora abbiamo l'elenco dei geni (colonna 4) che si sovrappongono alle regioni di picco,
come mostrato sopra.

Per avere una migliore visione d'insieme dei geni ottenuti, vogliamo esaminare la loro
distribuzione nei diversi cromosomi. Raggrupperemo la tabella per cromosoma e conteremo
il numero di geni con picchi su ciascun cromosoma

> <hands-on-title>conta i geni su diversi cromosomi</hands-on-title>
> 
> 1. {% tool [Group](Grouping1) %} dati in base a una colonna ed eseguire operazioni di aggregazione su altre colonne, con le seguenti impostazioni:
>     - *"Seleziona dati "* al risultato dell'intersezione
>     - *"Raggruppa per colonna "*:`Column 1`
>     - Premere **Operazione di inserimento** e scegliere:
>         - *"Tipo "*: `Count`
>         - *"Su colonna "*: `Column 1`
>         - *"Arrotondare il risultato al numero intero più vicino? "*: `No`
> 
>    > <question-title></question-title>
>    > 
>    > Quale cromosoma contiene il maggior numero di geni target?
>    > 
>    > > <solution-title></solution-title>
>    > > 
>    > > Il risultato varia a seconda delle impostazioni, ad esempio l'annotazione può
>    > > cambiare a causa degli aggiornamenti dell'UCSC. Se si segue il passaggio, con
>    > > la stessa annotazione, il risultato dovrebbe essere il cromosoma 11 con 2164
>    > > geni. Per garantire la riproducibilità, è necessario conservare tutti i dati di
>    > > input utilizzati nell'analisi. La ripetizione dell'analisi con lo stesso
>    > > insieme di parametri, memorizzati in Galaxy, può portare a un risultato diverso
>    > > se gli input sono cambiati, ad esempio l'annotazione di UCSC.
>    > {: .solution }
>    {: .question}
>
{: .hands_on}

## Visualizzazione

Abbiamo dei bei dati aggregati, quindi perché non disegnare un grafico a barre?

Prima di fare questo, però, dovremmo perfezionare i nostri dati raggruppati.

Si può notare che i cromosomi di topo non sono elencati nell'ordine corretto in questo
set di dati (lo strumento **Group** ha cercato di ordinarli, ma lo ha fatto in ordine
alfabetico).

Possiamo risolvere il problema eseguendo uno strumento dedicato all'ordinamento dei
dati.

> <hands-on-title>Correggere l'ordine della tabella dei conteggi dei geni</hands-on-title>
> 
> 1. {% tool [Sort](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_sort_header_tool/1.1.1) %} dati in ordine crescente o decrescente, con le seguenti impostazioni:
>     - *"Sort Query "*: risultato dell'esecuzione dello strumento Gruppo
>     - in {% icon param-repeat %} *"Selezioni di colonne "*
>       - *"su colonna "*: `Column 1`
>       - *"in "*: `Ascending order`
>       - *"Flavor "*: `Natural/Version sort (-V)`
> 
>    {% snippet faqs/galaxy-it/tools_similar.md toolname="Ordina i dati in ordine crescente o decrescente" toolversion="1.1.1" %}
> 
{: .hands_on}

Bene, siamo pronti a visualizzare!

> <hands-on-title>Disegna grafico a barre</hands-on-title>
> 
> 1. Fare clic sull'icona {% icon galaxy-barchart %} (visualizza) sull'output dello strumento **Sort**
> 2. selezionare `Bar diagram (NVD3)`
> 3. Fare clic sul pulsante **<<** nell'angolo in alto a destra
> 4. Scegliere un titolo in **Provvedere un titolo**, ad esempio `Gene counts per chromosome`
> 5. passare alla scheda {% icon galaxy-chart-select-data %} **Selezionare i dati** e
>    testare le impostazioni
> 6. Quando si è soddisfatti, fare clic sull'icona {% icon galaxy-save %} **Salva** in
>    alto a destra del *quadro principale*
> 
>    questo file verrà memorizzato nelle visualizzazioni salvate. In seguito sarà
>    possibile visualizzarla, scaricarla o condividerla con altri da **Dati ->
>    Visualizzazioni** nel menu superiore di Galaxy.
> 
{: .hands_on}

## Estrazione del flusso di lavoro / *workflow*

Osservando attentamente la cronologia, si può notare che contiene tutti i passaggi della
nostra analisi, dall'inizio alla fine. Costruendo questa cronologia, abbiamo costruito
un record completo della nostra analisi, con Galaxy che conserva tutte le impostazioni
dei parametri applicate in ogni fase. Non sarebbe bello convertire questa cronologia in
un flusso di lavoro da eseguire più volte?

Galaxy lo rende molto semplice con l'opzione `Extract workflow`. Ciò significa che ogni
volta che si desidera creare un flusso di lavoro, è possibile eseguirlo manualmente una
volta e poi convertirlo in un flusso di lavoro, in modo che la prossima volta sarà molto
meno faticoso eseguire la stessa analisi. Inoltre, consente di condividere o pubblicare
facilmente le analisi.

> <hands-on-title>Estrarre il flusso di lavoro</hands-on-title>
> 
> 1. **Pulisci** la tua cronologia: rimuovi tutti i lavori falliti (rossi) dalla tua cronologia facendo clic sul pulsante {% icon galaxy-delete %}.
> 
>    Questo faciliterà la creazione del flusso di lavoro.
> 
> 2. Fare clic su {% icon galaxy-gear %} (**Opzioni cronologia**) nella parte superiore del pannello cronologia e selezionare **Estrai flusso di lavoro**.
> 
>    ![`Extract Workflow` nel menu delle opzioni della cronologia](../../images/history_menu_extract_workflow.png)
> 
>    Il pannello centrale mostrerà il contenuto della cronologia in ordine inverso (il
>    più vecchio in cima) e sarà possibile scegliere quali passaggi includere nel flusso
>    di lavoro.
> 
> 3. Sostituire il **nome del flusso di lavoro** con qualcosa di più descrittivo, ad
>    esempio: `From peaks to genes`
> 
> 4. Se ci sono dei passaggi che non dovrebbero essere inclusi nel flusso di lavoro, è
>    possibile **deselezionarli** nella prima colonna di caselle.
> 
>    poiché abbiamo eseguito alcuni passaggi specifici per il nostro file di picco
>    personalizzato, potremmo voler escludere:
>    - **Seleziona per ultimo** {% icon tool %}
>    - tutti i passaggi **Sostituisci testo** {% icon tool %}
>    - **Convertire gli intervalli genomici in BED**
>    - **Prendere i fianchi** {% icon tool %}
> 
> 5. Fare clic sul pulsante **Crea flusso di lavoro** in alto.
> 
>    Verrà visualizzato un messaggio che indica che il flusso di lavoro è stato creato.
>    Ma dove è andato a finire?
> 
> 6. Fare clic su **Flusso di lavoro** nel menu a sinistra di Galaxy
> 
>    Qui è presente un elenco di tutti i flussi di lavoro
> 
> 7. Selezionare il flusso di lavoro appena generato e fare clic su **Modifica**
> 
>    Si dovrebbe vedere qualcosa di simile a questo:
> 
>    ![Interfaccia di modifica del flusso di lavoro](../../images/intro_06.png)
> 
>    > <comment-title>L'editor del flusso di lavoro</comment-title> Possiamo esaminare
>    > il flusso di lavoro nell'editor del flusso di lavoro di Galaxy. Qui è possibile
>    > visualizzare/modificare le impostazioni dei parametri di ogni fase, aggiungere e
>    > rimuovere strumenti e collegare l'uscita di uno strumento all'ingresso di un
>    > altro, il tutto in modo semplice e grafico. È inoltre possibile utilizzare questo
>    > editor per creare flussi di lavoro da zero.
>    {: .comment}
>
>     Sebbene abbiamo i nostri due input nel flusso di lavoro, manca loro la connessione con il primo strumento (**Intersect** {% icon tool %}), perché non abbiamo riportato alcuni dei passaggi intermedi.
> 
> 8. Collegate ogni set di dati di input allo strumento **Intersect** {% icon tool %} trascinando la freccia rivolta verso l'esterno a destra del suo riquadro (che denota un output) a una freccia rivolta verso l'interno a sinistra del riquadro **Intersect** (che denota un input)
> 9. Rinominare i set di dati di input in `Reference regions` e `Peak regions`
> 10. Premere **Auto Re-layout** per ripulire la nostra vista
>     ![Auto re-layouting](../../images/intro_07.png)
> 11. cliccare sull'icona {% icon galaxy-save %} *icona *Salva** (in alto) per salvare le modifiche
>     ![Pulsante Salva flusso di lavoro]({% link topics/contributing/images/save_workflow.png %}){: width="50%"}
> 
>    > <tip-title>Nascondere le fasi intermedie</tip-title>
>    > Quando si esegue un flusso di lavoro, di solito l'utente è interessato principalmente al prodotto finale e non a
>    > tutte le fasi intermedie. Per impostazione predefinita, tutti gli output di un
>    > flusso di lavoro vengono mostrati, ma è possibile indicare esplicitamente a Galaxy
>    > quali output mostrare e quali nascondere per un determinato flusso di lavoro. Questo
>    > comportamento è controllato dal piccolo asterisco accanto a ogni set di dati di
>    > output:
>    > 
>    > ![Workflow editor mark > output](../../../../shared/images/workflow_editor_mark_output.png)
>    > 
>    > se si fa clic su questo asterisco per uno qualsiasi dei set di dati di output,
>    > verranno mostrati *solo* i file con l'asterisco e tutti gli output senza asterisco
>    > verranno nascosti (si noti che fare clic su *tutti* gli output ha lo stesso effetto
>    > di fare clic su *nessuno* degli output, in entrambi i casi verranno mostrati tutti i
>    > set di dati).
>    {: .tip}
> 
{: .hands_on}

Ora è il momento di riutilizzare il nostro flusso di lavoro per un approccio più
sofisticato.

# Parte 2: approccio più sofisticato

Nella prima parte abbiamo utilizzato una definizione di sovrapposizione di 1 bp
(impostazione predefinita) per identificare i geni associati alle regioni di picco.
Tuttavia, i picchi potrebbero essere ampi, quindi, per ottenere una definizione più
significativa, potremmo identificare i geni che si sovrappongono nel punto in cui si
concentra la maggior parte delle letture, il **capo del picco**. Utilizzeremo le
informazioni sulla posizione del vertice del picco contenute nel file del picco
originale e controlleremo la sovrapposizione dei vertici con i geni.

## Preparazione

Abbiamo di nuovo bisogno del nostro file di picco, ma vorremmo lavorare in una
cronologia pulita. Invece di caricarlo due volte, possiamo copiarlo in una nuova
cronologia.

> <hands-on-title>Copia elementi della cronologia</hands-on-title>
> 
> 1. Creare una nuova cronologia e darle un nuovo nome come `Galaxy Introduction Part 2`
> 
>    {% snippet faqs/galaxy-it/histories_create_new.md %}
> 
> 2. Fare clic su **Opzioni cronologia** in alto a destra della cronologia. Fare clic su
>    **Mostra cronologia affiancata**
> 
>    Ora si dovrebbero vedere entrambe le cronologie affiancate
> 
> 3. trascinare e rilasciare il file del picco modificato (`Peak regions`, dopo i
>    passaggi di sostituzione), che contiene le informazioni sulla cima, nella nuova
>    cronologia.
> 4. Fare clic sul nome di Galaxy nella barra dei menu in alto a sinistra per
>    tornare alla finestra di analisi
> 
{: .hands_on}

## Creare il file del picco di vetta

Dobbiamo generare un nuovo file BED dal file di picco originale che contenga le
posizioni dei vertici dei picchi. L'inizio del picco è l'inizio del picco (colonna 2)
più la posizione all'interno del picco che ha la massima copertura del frammento di DNA
ipotetico (colonna 5, arrotondata al numero intero più piccolo perché alcuni picchi
cadono tra due basi). Come fine della regione del picco, definiremo semplicemente `start + 1`.

> <hands-on-title>Creare il file della cima del picco</hands-on-title>
> 
> 1. {% tool [Compute on rows](toolshed.g2.bx.psu.edu/repos/devteam/column_maker/Add_a_column1/2.0) %} con i seguenti parametri:
>   - *"File di input "*: il nostro file di picco `Peak regions` (il file in formato intervallo)
>   - *"L'input ha una riga di intestazione con nomi di colonne?": `No`
>   - In *"Espressioni "*:
>       - {% icon param-repeat %} *"Espressioni "*
>           - *"Aggiungi espressione "*: `c2 + int(c5)`
>           - *"Modalità dell'operazione "*: Aggiungi
>       - {% icon param-repeat %} *"Espressioni "*
>           - *"Aggiungi espressione "*: `c8 + 1`
>           - *"Modalità dell'operazione "*: Aggiungi
> 
>   Questo creerà un'ottava e una nona colonna nella nostra tabella, che utilizzeremo nel prossimo passo:
> 
> 2. Rinomina l'uscita `Peak summit regions`
> 
{: .hands_on}

Ora tagliamo solo il cromosoma più l'inizio e la fine del vertice:

> <hands-on-title>Tagliare le colonne</hands-on-title>
> 1. {% tool [Cut](Cut1) %} colonne di una tabella con le seguenti impostazioni:
>   - *"Tagliare colonne "*: `c1,c8,c9`
>   - *"Delimitato da Tab "*: `Tab`
>   - *"Da "*: `Peak summit regions`
> 
> L'output di **Cut** sarà in formato `tabular`.
> 
> 2. Cambiare il formato in `interval` (usare l'icona {% icon galaxy-pencil %}) poiché è quello che si aspetta lo strumento **Intersect**.
> 
>    {% snippet faqs/galaxy-it/datasets_change_datatype.md datatype="interval" %}
> 
>    L'output dovrebbe essere simile al seguente:
> 
>    ![Cime di picco](../../images/intro_summits.png){: width="200px"}
> 
{: .hands_on}

## Ottenere i nomi dei geni

I geni RefSeq che abbiamo scaricato dall'UCSC contenevano solo gli identificatori
RefSeq, ma non i nomi dei geni. Per ottenere un elenco di nomi di geni alla fine,
utilizziamo un altro file BED dalle Librerie di dati.

> <comment-title></comment-title>
> Ci sono diversi modi per inserire i nomi dei geni, se si ha bisogno di farlo da soli. Un modo è recuperare una mappatura attraverso Biomart
> e poi unire i due file (**Unisci due set di dati affiancati su un campo specificato**
> {% icon tool %}). Un altro metodo è quello di ottenere la tabella RefSeq completa da
> UCSC e convertirla manualmente in formato BED.
> 
{: .comment}

> <hands-on-title>Caricamento dei dati</hands-on-title>
> 
> 1. {% tool [Import](upload1) %} `mm9.RefSeq_genes_from_UCSC.bed` da [Zenodo](https://zenodo.org/record/1025586) o dalla libreria dati:
> 
>    ```
>    https://zenodo.org/record/1025586/files/mm9.RefSeq_genes_from_UCSC.bed
>    ```
> 
>    {% snippet faqs/galaxy-it/datasets_import_via_link.md genome="mm9" %}
> 
>    {% snippet faqs/galaxy-it/datasets_import_from_data_library.md path='Cliccare su "GTN - Material", "Introduction to Galaxy Analyses", "From peaks to genes", e poi "DOI: 10.5281/zenodo.1025586"' %}
> 
>    Per impostazione predefinita, Galaxy prende il link come nome, quindi rinominarli.
> 
> 2. Ispezionare il contenuto del file per verificare se contiene nomi di geni. Dovrebbe
>    essere simile al seguente: ![Nomi dei geni](../../images/intro_gene_names.png)
> 
> 3. Rinominalo `mm9.RefSeq_genes`
> 4. Applica il tag `#genes`
> 
{: .hands_on}

## Ripetizione del flusso di lavoro

È il momento di riutilizzare il flusso di lavoro creato in precedenza.

> <hands-on-title>Eseguire un flusso di lavoro</hands-on-title>
> 1. Aprire il menu del flusso di lavoro (barra dei menu a sinistra)
> 2. Trovare il flusso di lavoro creato nella sezione precedente e selezionare l'opzione **Esegui**
> 3. Scegliere come input il nostro file BED `mm9.RefSeq_genes` (`#genes`) e il risultato dello strumento **Cut** (`#peaks`)
> 4. Fare clic su **Eseguire flusso di lavoro**
> 
>    I risultati dovrebbero apparire nella cronologia, ma potrebbe volerci del tempo
>    prima che vengano completati.
{: .hands_on}

Abbiamo utilizzato il nostro flusso di lavoro per eseguire nuovamente l'analisi con i
picchi. Lo strumento **Group** ha nuovamente prodotto un elenco contenente il numero di
geni trovati in ciascun cromosoma. Ma non sarebbe più interessante conoscere il numero
di picchi in ogni singolo gene? Eseguiamo nuovamente il flusso di lavoro con
impostazioni diverse!

> <hands-on-title>Eseguire un flusso di lavoro con le impostazioni cambiate</hands-on-title>
> 1. Aprire il menu del flusso di lavoro (barra dei menu a sinistra)
> 2. Trovare il flusso di lavoro creato nella sezione precedente e selezionare l'opzione **Esegui**
> 3. Scegliere come input il nostro file BED `mm9.RefSeq_genes` (`#genes`) e il
>    risultato dello strumento **Cut** (`#peaks`)
> 4. Fare clic sul titolo dello strumento {% icon tool %} **Gruppo** per espandere le opzioni.
> 5. Modificare le seguenti impostazioni facendo clic sull'icona {% icon galaxy-pencil %} (matita) a sinistra:
>     - *"Raggruppa per colonna "*: `7`
>     - In *"Operazione "*:
>       - *"Su colonna "*: `7`
> 6. Fare clic su **Eseguire flusso di lavoro**
> 
{: .hands_on}

Congratulazioni! Dovresti avere un file con tutti i nomi unici dei geni e un
conteggio di quanti picchi contengono.

> <question-title></question-title>
> 
> L'elenco dei geni unici non è ordinato. Prova a ordinarlo da solo!
> 
> > <solution-title></solution-title>
> > È possibile utilizzare lo strumento "Ordina i dati in ordine crescente o decrescente" sulla colonna 2 e "ordinamento numerico veloce".
> {: .solution }
{: .question}


# Condividi il tuo lavoro

Una delle caratteristiche più importanti di Galaxy si concretizza alla fine di un'analisi.
Quando si pubblicano risultati eclatanti, è importante che altri ricercatori siano in
grado di riprodurre l'esperimento in silico. Galaxy consente agli utenti di condividere
facilmente i loro flussi di lavoro e le loro cronologie con altri.

Per condividere una cronologia, fare clic sulle opzioni della cronologia {% icon
galaxy-history-options %} e selezionare `Share or Publish`. In questa pagina si possono
fare 3 cose:


1. **Rendere accessibile tramite link**

   Questo genera un link che si può dare ad altri. Chiunque abbia questo link potrà
   visualizzare la propria cronologia.

2. **Mettere la cronologia a disposizione del pubblico in Cronologie pubblicate**

   Questo non solo crea un collegamento, ma pubblica anche la cronologia. Ciò significa
   che la cronologia sarà elencata sotto `Data → Histories → Published Histories` nel
   menu in alto.

3. **Condivisione con i singoli utenti**

   Questo condividerà la cronologia solo con utenti specifici dell'istanza Galaxy.

![Il menu per la condivisione delle cronologie comprende i pulsanti per rendere accessibile la cronologia, pubblicarla su questo server Galaxy e visualizzare un link condivisibile alla cronologia. In fondo c'è un pulsante per condividere la cronologia con i singoli utenti](../../images/publish.png)

> <hands-on-title>Condivisione della storia e del flusso di lavoro</hands-on-title>
> 
> 1. Condividi una delle tue storie con un tuo collega
> 2. Vedete se riuscite a fare lo stesso con il vostro flusso di lavoro!
> 3. Trova la cronologia e/o il flusso di lavoro condiviso dal tuo collega
> 
>    Le cronologie condivise con utenti specifici possono essere consultate da tali
>    utenti con `Data → Histories → Histories shared with me`.
> 
{: .hands_on}

# Conclusione


{% icon trophy %} Avete appena eseguito la vostra prima analisi in Galaxy. Avete anche
creato un flusso di lavoro dalla vostra analisi, in modo da poter ripetere facilmente la
stessa analisi su altri set di dati. Inoltre, avete condiviso i vostri risultati e
metodi con altri.

