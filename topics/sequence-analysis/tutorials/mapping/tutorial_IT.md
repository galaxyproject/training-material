---
layout: tutorial_hands_on
title: Mappatura
zenodo_link: https://doi.org/10.5281/zenodo.1324070
questions:
  - Che cos’è la mappatura?
  - Quali due elementi sono cruciali per una mappatura corretta?
  - Che cos’è un file BAM?
objectives:
  - Eseguire uno strumento per mappare le letture su un genoma di riferimento
  - Spiegare cos’è un file BAM e cosa contiene
  - Usare un browser genomico per comprendere i propri dati
time_estimation: 1h
key_points:
  - Conosci i tuoi dati!
  - La mappatura non è banale
  - Esistono molti algoritmi di mappatura – la scelta dipende dai tuoi dati
subtopic: basics
requirements:
- type: internal
  topic_name: sequence-analysis
  tutorials:
  - quality-control
follow_up_training:
- type: internal
  topic_name: transcriptomics
  tutorials:
  - ref-based
- type: internal
  topic_name: epigenetics
  tutorials:
  - formation_of_super-structures_on_xi
level: Introductory
edam_ontology:
- topic_0102
contributions:
  authorship:
  - joachimwolff
  - bebatut
  - hexylena
  translation:
  - lisanna
  - silviadg87
  - unode
  editing:
  - tflowers15
  funding:
  - biont
  - unimelb
  - melbournebioinformatics
  - AustralianBioCommons
recordings:
- captioners:
  - shiltemann
  date: '2021-02-15'
  galaxy_version: '21.01'
  length: 20M
  youtube_id: 1wm-62E2NkY
  speakers:
  - pvanheus
- youtube_id: B-CFWSkZzRc
  length: 24M
  galaxy_version: 24.1.2.dev0
  date: '2024-09-07'
  speakers:
  - DinithiRajapaksha
  captioners:
  - DinithiRajapaksha
  bot-timestamp: 1725707919
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



Il sequenziamento produce una raccolta di sequenze prive di contesto genomico.
Non sappiamo a quale parte del genoma corrispondano queste sequenze.
La mappatura delle letture di un esperimento su un genoma di riferimento è una fase fondamentale nell’analisi moderna dei dati genomici. Attraverso la mappatura, le letture vengono assegnate a una posizione specifica nel genoma, consentendo di ottenere informazioni come il livello di espressione dei geni.

Le letture non contengono informazioni sulla posizione, quindi non sappiamo da quale parte del genoma provengano.
Dobbiamo utilizzare la sequenza della lettura stessa per trovare la regione corrispondente nella sequenza di riferimento. Tuttavia, la sequenza di riferimento può essere molto lunga (circa 3 miliardi di basi per l’uomo), il che rende difficile individuare la regione corrispondente. Poiché le nostre letture sono brevi, potrebbero esserci più posizioni ugualmente probabili nella sequenza di riferimento da cui esse potrebbero derivare, questo è particolarmente vero per le regioni ripetitive.

In linea di principio, potremmo eseguire un’analisi BLAST per determinare dove si adattano meglio i frammenti sequenziati nel genoma noto. Tuttavia, sarebbe necessario eseguire questa operazione per ciascuna delle milioni di letture presenti nei dati di sequenziamento.
Allineare milioni di brevi sequenze in questo modo richiederebbe settimane di elaborazione. Inoltre, non ci interessa la corrispondenza esatta base per base (alignment), ma piuttosto da quale parte del genoma provengono le letture.
Questo approccio è chiamato **mappatura**.

Di seguito, elaboreremo un set di dati con il mappatore **Bowtie2** e visualizzeremo i dati con il programma **IGV**. 

> <agenda-title></agenda-title>
> 
> In questo tutorial, ci occuperemo di:
> 
> 1. TOC
> {:toc}
{: .agenda}

# Preparazione dei dati

> <hands-on-title>Caricamento dei dati</hands-on-title>
> 
> 1. Creare una nuova storia per questa esercitazione e assegnarle un nome appropriato
> 
>    {% snippet faqs/galaxy-it/histories_create_new.md %}
> 
>    {% snippet faqs/galaxy-it/histories_rename.md %}
> 
> 2. Importare `wt_H3K4me3_read1.fastq.gz` e `wt_H3K4me3_read2.fastq.gz` da
>    [Zenodo](https://zenodo.org/record/1324070) o dalla libreria dei dati (chiedere al
>    proprio docente)
> 
>    ```
>    https://zenodo.org/record/1324070/files/wt_H3K4me3_read1.fastq.gz
>    https://zenodo.org/record/1324070/files/wt_H3K4me3_read2.fastq.gz
>    ```
> 
>    {% snippet faqs/galaxy-it/datasets_import_via_link.md %}
> 
>    {% snippet faqs/galaxy-it/datasets_import_from_data_library.md %}
> 
>    Per impostazione predefinita, Galaxy assegna come nome il link stesso, quindi è necessario rinominarli.
> 
> 3. Rinominare i file in `reads_1` e `reads_2`
> 
>    {% snippet faqs/galaxy-it/datasets_rename.md %}
>
> 4. Crea una raccolta di coppie denominata `Paired Reads` 
>
>    {% snippet faqs/galaxy-it/collections_build_list_paired.md %}
> 
{: .hands_on}

Abbiamo appena importato in Galaxy i file FASTQ corrispondenti a dati paired-end, come quelli che si ottengono direttamente da un centro di sequenziamento.
Durante il sequenziamento possono essere introdotti errori, come nucleotidi chiamati in modo errato.
Tali errori possono influenzare l’analisi e portare a un’interpretazione sbagliata dei dati.
Il primo passo in qualsiasi analisi di dati di sequenziamento è sempre verificare la qualità delle letture.

Esiste un tutorial dedicato al [controllo di qualità]({% link topics/sequence-analysis/tutorials/quality-control/tutorial.md %}) dei dati di sequenziamento. NNon ripeteremo qui i passaggi descritti. Ti consigliamo di seguire il
[tutorial]({% link topics/sequence-analysis/tutorials/quality-control/tutorial.md %}) e di applicarlo ai tuoi dati prima di proseguire.

# Mappare le letture su un genoma di riferimento

La mappatura delle letture è il processo di allineamento delle sequenze lette a un genoma di riferimento. Un mapper prende in input un genoma di riferimento e un insieme di letture, con l’obiettivo di allineare ciascuna lettura al genoma di riferimento, consentendo mismatch, inserzioni/delezioni (indel) e clipping di piccoli frammenti alle due estremità delle letture:

![Spiegazione della mappatura](../../images/mapping/mapping.png "L'input è costituito da un insieme di letture e da un genoma di riferimento. Al centro sono mostrati i risultati della mappatura: le posizioni delle letture sul genoma di riferimento.
La prima lettura è allineata alla posizione 100 e presenta due mismatch. La seconda lettura è allineata alla posizione 114: si tratta di un allineamento locale con clipping (ritagli) alle estremità sinistra e destra. La terza lettura è allineata alla posizione 123 e mostra un'inserzione di 2 basi e una delezione di 1 base.")

Abbiamo bisogno di un genoma di riferimento su cui mappare le letture.

{% include topics/sequence-analysis/tutorials/mapping/ref_genome_explanation_IT.md answer_3="Questi dati provengono dal ChIP-seq dei topi, quindi utilizzeremo mm10 (*Mus musculus*)"%}

Attualmente esistono oltre 60 diversi mappatori, e il loro numero continua a crescere. In questo tutorial utilizzeremo [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/), uno strumento open source veloce ed efficiente in termini di memoria, particolarmente adatto per l’allineamento di letture di lunghezza compresa tra circa 50 e diverse migliaia di basi su genomi relativamente grandi.

> <hands-on-title>Mappatura con Bowtie2</hands-on-title>
> 1. {% tool [Bowtie2](toolshed.g2.bx.psu.edu/repos/devteam/bowtie2/bowtie2/2.4.2+galaxy0) %} con i seguenti parametri
>    - *"Is this single or paired library"*: `Paired-end`
>       - {% icon param-file %} *"FASTA/Q file #1"*: `reads_1`
>       - {% icon param-file %} *"FASTA/Q file #2"*: `reads_2`
>       - *"Do you want to set paired-end options?"*: `No`
> 
>         È comunque utile esaminare i parametri disponibili, in particolare l’orientamento dei mate, se conosciuto. Queste opzioni possono migliorare la qualità della mappatura paired-end.
> 
>     - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a built-in genome index`
>       - *"Select reference genome"*: `Mouse (Mus musculus): mm10`
>     - *"Select analysis mode"*: `Default setting only`
> 
>       È consigliabile consultare anche i parametri non predefiniti e comprenderne la funzione, poiché possono influenzare la qualità e l’efficienza della mappatura.
> 
>     - *"Save the bowtie2 mapping statistics to the history"*: `Yes`
> 
> 2. Ispezionare il file `mapping stats` cliccando sull'icona {% icon galaxy-eye %} (occhio)
> 
{: .hands_on}

> <question-title></question-title>
> 
> 1. Quali informazioni sono fornite in questo file?
> 2. Quante letture sono state mappate esattamente una volta?
> 3. Quante letture sono state mappate più di una volta? Come è possibile? Cosa si dovrebbe fare con esse?
> 4. Quante coppie di letture non sono state mappate? Quali possono essere le cause?
> 
> > <solution-title></solution-title>
> > 1. Il file fornisce informazioni quantitative: mostra quante sequenze sono state allineate, ma non dà indicazioni dirette sulla qualità dell’allineamento.
> > 2. Circa il 90% delle letture è stato allineato esattamente una volta
> > 3. Circa il 7% delle letture è stato allineato concordemente più di una volta. Queste sono dette multi-mapped reads. Ciò può accadere a causa di regioni ripetitive nel genoma di riferimento (ad esempio copie multiple di un gene), soprattutto quando le letture sono brevi. È difficile stabilire la loro origine esatta, perciò la maggior parte delle pipeline le ignora. È comunque importante verificare queste statistiche per assicurarsi di non escludere troppe informazioni nelle analisi successive.
> > 4. Circa il 3% delle coppie di letture non è stato mappato, perché:
> >     - entrambe le letture della coppia sono allineate, ma le loro posizioni non concordano (`aligned discordantly 1 time`)
> >     - le letture della coppia sono multi-mappate (`aligned >1 times` in `pairs
> >       aligned 0 times concordantly or discordantly`)
> >     - una delle due letture è mappata ma non la sua compagna (`aligned
> >       exactly 1 time` in `pairs aligned 0 times concordantly or discordantly`)
> >     - il resto non viene mappato affatto
> {: .solution }
{: .question}

Verificare le statistiche di mappatura è un passaggio cruciale prima di proseguire con qualsiasi analisi.
Esistono numerose possibili fonti di errore nella mappatura, tra cui (ma non solo): 

- **Artefatti della reazione a catena della polimerasi (PCR)**: Molti metodi di
  sequenziamento ad alta velocità (HTS) prevedono una o più fasi di PCR. Gli errori di
  PCR si manifestano come mismatch nell'allineamento e, in particolare, gli errori nei
  primi cicli di PCR si manifestano con letture multiple, suggerendo falsamente una
  variazione genetica nel campione. Un errore correlato è rappresentato dai duplicati di
  PCR, in cui la stessa coppia di letture si presenta più volte, alterando i calcoli di
  copertura nell'allineamento.
- **Errori di sequenziamento**: La macchina di sequenziamento può effettuare una
  chiamata errata per motivi fisici (ad esempio, olio su un vetrino Illumina) o a causa
  delle proprietà del DNA sequenziato (ad esempio, omopolimeri). Poiché gli errori di
  sequenziamento sono spesso casuali, possono essere filtrati come letture singleton
  durante la chiamata di variante.
- **Errori di mappatura**: L'algoritmo di mappatura può mappare una lettura nella
  posizione sbagliata del riferimento. Ciò accade spesso in prossimità di ripetizioni o
  di altre regioni a bassa complessità.

Pertanto, se le statistiche di mappatura non sono buone, è necessario indagare sulla
causa di questi errori prima di procedere con le analisi.

Dopo di che, è necessario dare un'occhiata alle letture e ispezionare il file BAM in cui
sono memorizzate le mappature delle letture.

# Ispezione di un file BAM

{% include topics/sequence-analysis/tutorials/mapping/bam_explanation_IT.md mapper="Bowtie2" %}

Il file BAM include molte informazioni, in particolare sulla qualità della mappatura.

> <hands-on-title>Riepilogo della qualità di mappatura</hands-on-title>
> 1. {% tool [Samtools Stats](toolshed.g2.bx.psu.edu/repos/devteam/samtools_stats/samtools_stats/2.0.2+galaxy2) %} con i seguenti parametri
>    - {% icon param-file %} *"BAM file"*: `aligned reads` (output of **Bowtie2** {% icon tool %})
>    - *"Use reference sequence"*: `Locally cached/Use a built-in genome`
>      - *"Using genome"*: `Mouse (Mus musculus): mm10 Full`
> 
> 2. Ispezionare il file {% icon param-file %} file `Stats`
{: .hands_on}

> <question-title></question-title>
> 
> 1. Qual è la percentuale di mismatch nelle letture mappate rispetto al genoma di riferimento?
> 2. Cosa rappresenta il tasso di errore?
> 3. Qual è la qualità media della mappatura e come viene rappresentata?
> 4. Qual è la dimensione media degli inserti?
> 5. Quante letture hanno un punteggio di qualità di mappatura inferiore a 20?
> 
> > <solution-title></solution-title>
> > 1. Ci sono ~21.900 mismatches per ~4.753.900 basi mappate, il che produce in media
> >    ~0,005 mismatches per basi mappate.
> > 2. Il tasso di errore è la proporzione di mismatch per basi mappate, quindi il
> >    rapporto calcolato subito prima.
> > 3. La qualità media è il punteggio medio di qualità della mappatura. Si tratta di un
> >    punteggio Phred come quello utilizzato nel file FASTQ per ciascun nucleotide. Ma
> >    qui il punteggio non è per nucleotide, ma per lettura e rappresenta la
> >    probabilità della qualità della mappatura.
> > 4. La dimensione dell'inserto è la distanza tra le due letture nelle coppie.
> > 5. Per ottenere le informazioni:
> >      1. {% tool [Filter BAM](toolshed.g2.bx.psu.edu/repos/devteam/bamtools_filter/bamFilter/2.5.2+galaxy2) %} con un filtro per mantenere solo le letture con una qualità di mappatura >= 20
> >      2. {% tool [Samtools Stats](toolshed.g2.bx.psu.edu/repos/devteam/samtools_stats/samtools_stats/2.0.5) %} sull'output di **Filter**
> > 
> >    Prima del filtraggio: 95.412 letture e dopo il filtraggio: 89.664 letture.
> {: .solution }
{: .question}

# Visualizzazione con il browser del genoma

## IGV

L’Integrative Genomics Viewer (IGV) è uno strumento di visualizzazione ad alte prestazioni per l’esplorazione interattiva di grandi insiemi di dati genomici integrati.
Supporta un’ampia varietà di tipi di dati, tra cui i dati di sequenziamento basati su microarray e di nuova generazione (NGS), oltre alle annotazioni genomiche.
Di seguito lo utilizzeremo per visualizzare le letture mappate.

{% include topics/sequence-analysis/tutorials/mapping/igv_IT.md tool="Bowtie2" region_to_zoom="chr2:98,666,236-98,667,473" %}

## JBrowse

{% tool [JBrowse](toolshed.g2.bx.psu.edu/repos/iuc/jbrowse/jbrowse/1.16.11+galaxy0) %} è
è un browser genomico alternativo basato sul web.
Mentre IGV è un software da scaricare ed eseguire localmente, JBrowse è accessibile tramite un’interfaccia web ospitata online, che consente di esplorare i dati genomici direttamente dal browser.
Lo useremo per visualizzare le letture mappate.

{% include topics/sequence-analysis/tutorials/mapping/jbrowse_IT.md tool="Bowtie2" region_to_zoom="chr2:98,666,236-98,667,473" %}

# Conclusione

Dopo il controllo di qualità, la mappatura rappresenta una fase fondamentale nella maggior parte delle analisi di dati di sequenziamento (RNA-Seq, ChIP-Seq, ecc.).
Serve a determinare l’origine delle letture nel genoma di riferimento e a utilizzare queste informazioni nelle analisi successive (downstream analyses).

