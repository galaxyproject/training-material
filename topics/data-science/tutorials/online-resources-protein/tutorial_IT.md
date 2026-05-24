---
layout: tutorial_hands_on
title: Una proteina lungo la pagina UniProt
level: Introductory
zenodo_link: ''
questions:
- Come si possono cercare proteine utilizzando testo, geni o nomi di proteine?
- Come si interpretano le informazioni in cima alla pagina di UniProt?
- Che tipo di informazioni ci si può aspettare da diversi formati di download, come FASTA e JSON?
- Come viene descritta la funzione di una proteina come le opsine nella sezione "Funzione"?
- Quali informazioni strutturate si trovano nelle sezioni "Nomi e tassonomia", "Localizzazione subcellulare", "Malattia e varianti", "PTM/Elaborazione"?
- Come si possono conoscere l'espressione proteica, le interazioni, la struttura, la famiglia, la sequenza e le proteine simili?
- In che modo le schede "Visualizzatore varianti" e "Visualizzatore caratteristiche" aiutano a mappare le informazioni proteiche lungo la sequenza?
- Cosa elenca la scheda "Pubblicazioni" e come si possono filtrare le pubblicazioni?
- Qual è l'importanza di monitorare le modifiche delle annotazioni delle voci nel tempo?
objectives:
- Esplorando le voci proteiche in UniProtKB, è possibile interpretare la funzione, la tassonomia, la struttura, le interazioni, le varianti e altro ancora delle proteine.
- Utilizzare identificatori univoci per collegare database, scaricare dati su geni e proteine, visualizzare e confrontare le caratteristiche delle sequenze.
time_estimation: 1H
key_points:
- Come navigare tra le voci di UniProtKB, accedendo a dettagli completi sulle proteine, come funzioni, tassonomia e interazioni.
- Il visualizzatore di varianti e caratteristiche è uno strumento utile per esplorare visivamente varianti proteiche, domini, modifiche e altre caratteristiche chiave della sequenza.
- Amplia le tue conoscenze utilizzando link esterni per incrociare i dati e scoprire relazioni complesse.
- Esplora la scheda Cronologia per accedere alle versioni precedenti delle annotazioni delle voci.
requirements:
- type: internal
  topic_name: data-science
  tutorials:
  - online-resources-gene
contributions:
  authorship:
  - lisanna
  - bebatut
  translation:
  - lisanna
  - silviadg87
  - unode
  funding:
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


Quando si esegue un'analisi di dati biologici, è possibile che ci si ritrovi con alcune proteine interessanti e che sia necessario esplorare questi geni. Ma come possiamo farlo? Quali sono le risorse disponibili per questo? E come navigare al loro interno?

Lo scopo di questo tutorial è di familiarizzare con questo aspetto, utilizzando come esempio le opsine umane.

> <comment-title></comment-title>
> Questo tutorial è un po' atipico: non lavoreremo in Galaxy ma soprattutto al di fuori di esso, nelle pagine del database [UniProt](https://uniprot.org).
> 
{: .comment}

> <comment-title></comment-title>
> Questo tutorial è stato progettato per essere la continuazione del tutorial ["Un gene attraverso i formati di file"]({% link topics/data-science/tutorials/online-resources-gene/tutorial.md %}), ma può essere consultato anche come modulo autonomo.
> 
{: .comment}

Le opsine si trovano nelle cellule della retina. Catturano la luce e iniziano la sequenza di segnali che porta alla visione, ed è per questo che, quando sono compromesse, sono associate al daltonismo e ad altri disturbi visivi.

> <comment-title>Fonti di informazione di questo tutorial</comment-title>
> 
> Il tutorial che state consultando è stato sviluppato principalmente consultando le risorse UniProtKB, in particolare il tutorial [Explore UniProtKB entry](https://www.uniprot.org/help/explore_uniprotkb_entry). Alcune frasi sono riportate da lì senza modifiche.
> 
> Inoltre, l'argomento è stato scelto sulla base del [Tutorial di bioinformatica] di Gale Rhodes (https://spdbv.unil.ch/TheMolecularLevel/Matics/index.html). Sebbene il tutorial non possa più essere seguito passo dopo passo a causa del modo in cui le risorse citate sono cambiate nel tempo, potrebbe fornire ulteriori spunti di riflessione sulle opsine e in particolare su come si possano costruire modelli strutturali di proteine basati su informazioni evolutive.
> 
{: .comment}


> <agenda-title></agenda-title>
> 
> In questo tutorial ci occuperemo di:
> 
> 1. TOC
> {:toc}
{: .agenda}

# La pagina di inserimento di UniProtKB

Il portale da visitare per ottenere tutte le informazioni su una proteina è [UniProtKB](https://www.uniprot.org/). Possiamo effettuare la ricerca utilizzando una ricerca testuale o il nome del gene o della proteina. Proviamo prima con una serie di parole chiave generiche, come `Human opsin`.

> <hands-on-title>Ricerca opsina umana su UniProtKB</hands-on-title>
> 
> 1. Aprire il file [UniProtKB](https://www.uniprot.org/)
> 2. Digitare `Human opsin` nella barra di ricerca
> 3. Avviare la ricerca
> 
{: .hands-on}

> <question-title></question-title>
> 
> Quanti risultati abbiamo ottenuto?
> 
> > <solution-title></solution-title>
> > 
> > 410 risultati (al momento della preparazione di questo tutorial)
> > 
> {: .solution}
{: .question}

Questi 410 risultati ci danno la sensazione di dover essere più specifici (anche se - spoiler - il nostro obiettivo effettivo è tra i primi risultati).

Per essere sufficientemente specifici, suggeriamo di usare un identificatore unico. Dal [tutorial precedente]({% link topics/data-science/tutorials/online-resources-gene/tutorial.md %}) conosciamo il nome del gene della proteina che stiamo cercando, `OPN1LW`.

> <hands-on-title>Ricerca di OPN1LW su UniProtKB</hands-on-title>
> 
> 1. Digitare `OPN1LW` nella barra di ricerca in alto
> 2. lanciare la ricerca
> 
{: .hands-on}

> <question-title></question-title>
> 
> 1. Quanti risultati abbiamo ottenuto?
> 2. Cosa dovremmo fare per ridurre questo numero?
> 
> > <solution-title></solution-title>
> > 
> > 1. 200+ risultati (al momento della preparazione di questo tutorial)
> > 2. Dobbiamo chiarire cosa stiamo cercando: OPN1LW umano
> {: .solution}
{: .question}

Dobbiamo aggiungere `Human` per chiarire cosa stiamo cercando.

> <hands-on-title>Ricerca di OPN1LW umana su UniProtKB</hands-on-title>
> 
> 1. Digitare `Human OPN1LW` nella barra di ricerca in alto
> 2. Avviare la ricerca
> 
{: .hands-on}

> <question-title></question-title>
> 
> 1. Quanti risultati abbiamo ottenuto?
> 2. Abbiamo un risultato che includa OPN1LW come nome del gene?
> 
> > <solution-title></solution-title>
> > 
> > 1. 7 risultati (al momento della preparazione di questo tutorial)
> > 2. Il primo risultato è etichettato con `Gene: OPN1LW (RCP)`
> {: .solution}
{: .question}

Il primo risultato, etichettato con `Gene: OPN1LW (RCP)`, è il nostro obiettivo, `P04000 · OPSR_HUMAN`. Prima di aprire la pagina, occorre notare due cose:

1. Il nome della proteina `OPSR_HUMAN` è diverso dal nome del gene, così come i loro ID sono.
2. Questa voce ha una stella dorata, il che significa che è stata annotata e curata manualmente.

# Ispezione di una voce UniProt

> <hands-on-title>Aprire un risultato in UniProt</hands-on-title>
> 
> 1. Clicca su `P04000 · OPSR_HUMAN`
{: .hands-on}

![Schermata dell'intestazione della pagina di inserimento di UniProt](./images/UniProt.png "pagina UniProt")

Per navigare in questa lunga pagina, il menu (barra di navigazione) a sinistra sarà estremamente utile. Solo da esso si capisce che questo database contiene informazioni sulla voce su:
- le funzioni note,
- la tassonomia,
- la posizione,
- varianti e malattie associate,
- Modifica post-traduzionale (PTM),
- l'espressione,
- le interazioni,
- la struttura,
- i domini e la loro classificazione,
- le sequenze
- proteine simili.

La barra di navigazione rimane nella stessa posizione sullo schermo mentre si sale e si scende in una voce, in modo da poter navigare rapidamente verso le sezioni di interesse. Consulteremo separatamente tutte le sezioni citate, ma prima concentriamoci sulle intestazioni a sinistra.

Nella parte superiore della pagina, è possibile vedere l'identificativo e il nome della voce UniProt, il nome della proteina e del gene, l'organismo, se la voce della proteina è stata rivista manualmente da un curatore UniProt, il suo punteggio di annotazione e il livello di evidenza della sua esistenza.

Sotto l'intestazione principale si trova una serie di schede (*Ingresso*, *Variant viewer*, *Feature viewer*, *Pubblicazioni*, *Link esterni*, *Storia*). Le schede consentono di passare dalla voce, alla visualizzazione grafica delle caratteristiche della sequenza (Feature viewer), alle pubblicazioni e ai collegamenti esterni, ma per il momento le ignorano e non si spostano dalla scheda *Entry*.

## Entrata

Il menu successivo fa già parte della scheda *Entry*. Permette di eseguire una ricerca di similarità di sequenza BLAST sulla voce, di allinearla con tutte le sue isoforme, di scaricare la voce in vari formati o di aggiungerla al carrello per salvarla in seguito.

> <question-title></question-title>
> 
> 1. Quali sono i formati disponibili nel menu a discesa *Download*?
> 2. Che tipo di informazioni scaricheremmo attraverso questi formati di file?
> 
> > <solution-title></solution-title>
> > 
> > 1. I formati sono: `Text`, `FASTA (canonical)`, `FASTA (canonical & isoform`, `JSON`, `XML`, `RDF/XML`, `GFF`
> > 2. I formati `FASTA` dovrebbero suonare familiari (dopo il tutorial preliminare), e includono la sequenza della proteina, eventualmente con le sue isoforme (nel qual caso si tratterà di un multi-FASTA). Oltre a questi, tutti gli altri formati non sono specifici per le proteine o per la biologia. Si tratta di formati di file generali ampiamente utilizzati dai siti web per includere le informazioni contenute nella pagina. Quindi, scaricando il file `text` (o meglio ancora il file `json`), si scaricherebbe la stessa annotazione a cui si accede in questa pagina, ma in un formato più facile da analizzare programmaticamente.
> > 
> {: .solution}
{: .question}

Scorriamo ora la pagina di ingresso, sezione per sezione.

### Funzione

Questa sezione riassume le funzioni di questa proteina come segue:

I pigmenti visivi sono le molecole che assorbono la luce e che sono alla base della visione. Sono costituiti da un'apoproteina, l'opsina, legata covalentemente alla cis-retina.

Indipendentemente dal livello di dettaglio che si comprende (a seconda del proprio background), questo è impressionantemente breve e specifico, considerando l'enorme quantità di letteratura e di studi che esistono al di là della determinazione della funzione di una proteina. Comunque, qualcuno ha fatto il lavoro per noi e questa proteina è già completamente classificata nella Gene Ontology (GO), che descrive la funzione molecolare, il processo biologico e il componente cellulare di ogni proteina classificata.

GO è un perfetto esempio di database/risorsa che si basa su un universo di conoscenze molto complesso e lo traduce in un grafico più semplice, con il rischio di perdere i dettagli. Questo ha il grande vantaggio di organizzare le informazioni, renderle conteggiabili, analizzabili e programmaticamente accessibili, permettendoci di ottenere queste lunghe pagine riassuntive e le KnowledgeBase.

> <question-title></question-title>
> 
> 1. A quali funzioni molecolari è associata questa proteina?
> 2. A quali componenti cellulari è associata questa proteina?
> 3. A quali processi biologici si riferisce questa proteina?
> 
> > <solution-title></solution-title>
> > 
> > 1. Proteina del fotorecettore, recettore accoppiato a proteine G
> > 2. Membrana del disco dei fotorecettori
> > 2. Trasduzione sensoriale, Visione
> > 
> {: .solution}
{: .question}

### Nomi e tassonomia

Altri esempi di informazioni strutturate sono disponibili nella sezione successiva, ad esempio nella tassonomia. Questa sezione riporta anche altri identificatori univoci che si riferiscono alla stessa entità biologica o a entità collegate (ad esempio, le malattie associate nel menu `MIM`).

> <question-title></question-title>
> 
> 1. Qual è l'identificatore tassonomico associato a questa proteina?
> 2. Qual è l'identificatore del proteoma associato a questa proteina?
> 
> > <solution-title></solution-title>
> > 
> > 1. 9606, cioè Homo sapiens
> > 2. UP000005640, componente del cromosoma Xs
> {: .solution}
{: .question}

## Posizione subcellulare

Sappiamo già dove si trova la nostra proteina nel corpo umano (nella retina, come specificato dal sommario della funzione), ma dove si trova nella cellula?

> <question-title></question-title>
> 
> 1. Dove si trova la nostra proteina nella cellula?
> 2. È coerente con l'annotazione GO osservata in precedenza?
> 
> > <solution-title></solution-title>
> > 
> > 1. La sezione spiega che si tratta di una "proteina di membrana a più passaggi", il che significa che è una proteina inserita nella membrana cellulare e che la attraversa più volte.
> > 2. L'annotazione GO in alto indica che ci stiamo riferendo in particolare alla membrana dei fotorecettori (cellule).
> {: .solution}
{: .question}

La sezione della localizzazione subcellulare comprende un'area *Caratteristiche* che specifica quali sezioni, lungo la sequenza della proteina, sono inserite nella membrana (Transmembrana) e quali no (Dominio topologico).

> <question-title></question-title>
> 
> Quanti domini transmembrana e domini topologici esistono?
> 
> > <solution-title></solution-title>
> > 
> > 8 domini transmembrana e 7 domini topologici
> {: .solution}
{: .question}

### Malattia e Varianti

Come sappiamo dal tutorial precedente, questo gene/proteina è associato a diverse malattie. Questa sezione illustra in dettaglio questa associazione, elencando anche le varianti specifiche che sono state individuate come correlate alla malattia.

> <question-title></question-title>
> 
> Quali tipi di studi scientifici permettono di valutare l'associazione di una variante genetica alle malattie?
> 
> > <solution-title></solution-title>
> > 
> > Tre metodi comunemente usati per valutare l'associazione di una variante genetica con una malattia sono:
> > 
> > - Studi di associazione a livello di genoma (GWAS)
> > 
> >   I GWAS sono ampiamente utilizzati per identificare le varianti genetiche comuni associate alle malattie. Comportano la scansione dell'intero genoma di un gran numero di individui per identificare le variazioni legate a una particolare malattia o tratto.
> > 
> > - Studi caso-controllo
> > 
> >   Gli studi caso-controllo sono spesso utilizzati per confrontare gli individui affetti da una malattia con quelli che non ne sono affetti, concentrandosi sulla presenza o sulla frequenza di specifiche varianti genetiche in entrambi i gruppi.
> > 
> > - Studi sulla famiglia
> > 
> >   Gli studi basati sulla famiglia prevedono l'analisi delle varianti genetiche all'interno delle famiglie in cui più membri sono affetti da una malattia. Studiando i modelli di ereditarietà delle varianti genetiche e la loro associazione con la malattia all'interno delle famiglie, i ricercatori possono identificare potenziali geni associati alla malattia.
> > 
> > Questo tipo di studi implica un uso estensivo dei tipi di file per la gestione dei dati genomici, come: SAM (Sequence Alignment Map), BAM (Binary Alignment Map), VCF (Variant Calling Format) ecc.
> > 
> {: .solution}
{: .question}

Questa sezione comprende anche un'area *Caratteristiche*, dove sono mappate le varianti naturali lungo la sequenza. Di seguito, si evidenzia anche che una vista più dettagliata delle caratteristiche lungo la sequenza è fornita nella scheda *Malattie e varianti*, ma non apriamola per ora.

### PTM/elaborazione

Una modifica post-traduzionale (PTM) è un evento di elaborazione covalente risultante da una scissione proteolitica o dall'aggiunta di un gruppo modificante a un amminoacido.

> <question-title></question-title>
> 
> Quali sono le modifiche post-traduzionali della nostra proteina?
> 
> > <solution-title></solution-title>
> > 
> > Catena, glicosilazione, legame disolfuro, residuo modificato
> {: .solution}
{: .question}

### Espressione

Sappiamo già dove si trova la proteina nella cellula, ma per le proteine umane abbiamo spesso informazioni su dove si trova nel corpo umano, cioè in quali tessuti. Queste informazioni possono provenire dallo Human [ExpressionAtlas](https://www.ebi.ac.uk/gxa/home) o da altre risorse simili.

> <question-title></question-title>
> 
> In quale tessuto si trova la proteina?
> 
> > <solution-title></solution-title>
> > 
> > I tre pigmenti colorati si trovano nelle cellule dei fotorecettori del cono.
> {: .solution}
{: .question}

### Interazione

Le proteine svolgono la loro funzione attraverso l'interazione con l'ambiente circostante, in particolare con altre proteine. Questa sezione riporta gli interagenti della nostra proteina di interesse, in una tabella che possiamo anche filtrare per localizzazione subcellulare, malattie e tipo di interazione.

La fonte di queste informazioni sono i database come STRING, e la pagina di ingresso della nostra proteina è direttamente collegata a questa sezione.

> <hands-on-title>Ricerca di OPN1LW umana su UniProtKB</hands-on-title>
> 
> 1. Fare clic sul [link STRINGA](https://string-db.org/network/9606.ENSP00000358967) in un'altra scheda
> 
{: .hands-on}

> <question-title></question-title>
> 
> 1. Quanti formati di file diversi si possono scaricare da lì?
> 2. Che tipo di informazioni saranno trasmesse in ogni file?
> 
> > <solution-title></solution-title>
> > 
> > STRING fornisce dati in formati di file scaricabili per supportare ulteriori analisi. Il formato di file principale utilizzato da STRING è il formato TSV (Tab-Separated Values), che presenta i dati di interazione proteica in un layout strutturato e tabellare. Questo formato è adatto a una facile integrazione in vari strumenti e software di analisi dei dati. Inoltre, STRING offre dati in formato XML PSI-MI (Proteomics Standards Initiative Molecular Interactions), uno standard per la rappresentazione dei dati di interazione proteica che consente la compatibilità con altri database di interazione e piattaforme di analisi. Questi formati di file consentono ai ricercatori di sfruttare la ricchezza di informazioni sulle interazioni proteiche di STRING per i propri studi e analisi. I ricercatori possono anche scaricare rappresentazioni visive delle reti proteiche in formati immagine come PNG e SVG, adatti per presentazioni e pubblicazioni. Per le analisi avanzate, STRING offre "file piatti" contenenti informazioni dettagliate sulle interazioni e file "MFA" (Multiple Alignment Format), utili per confrontare più sequenze di proteine. Questi diversi formati di file scaricabili consentono ai ricercatori di sfruttare la ricchezza di informazioni sulle interazioni proteiche di STRING per i propri studi e analisi.
> > 
> {: .solution}
{: .question}

### Struttura

Siete curiosi di conoscere le intricate strutture tridimensionali delle proteine? La sezione *Struttura* della pagina di ingresso di UniProtKB è la porta d'accesso per esplorare l'affascinante mondo dell'architettura delle proteine.

In questa sezione troverete informazioni sulle strutture proteiche determinate sperimentalmente. Queste strutture forniscono informazioni cruciali sul funzionamento delle proteine e sulla loro interazione con altre molecole. Si possono scoprire visualizzazioni interattive della struttura della proteina che possono essere esplorate direttamente all'interno della voce UniProtKB. Questa funzione offre un modo accattivante per navigare tra i domini, i siti di legame e altre regioni funzionali della proteina. Approfondendo la sezione *Struttura*, si potrà comprendere meglio la base fisica della funzione delle proteine e scoprire la ricchezza di informazioni che i dati strutturali possono svelare.

> <question-title></question-title>
> 
> 1. Qual è la variante associata al daltonismo?
> 2. Riesci a trovare quello specifico amminoacido nella struttura?
> 3. Puoi formulare un'ipotesi sul motivo per cui questa mutazione è dirompente?
> 
> > <solution-title></solution-title>
> > 
> > 1. Nella sezione *Malattie e varianti*, scopriamo che il cambiamento da glicina (G) ad acido glutammico (E) nella posizione 338 della sequenza proteica è associato al daltonismo.
> > 2. nel visualizzatore di strutture, possiamo muovere la molecola e passare il mouse sulla struttura per trovare l'AA in posizione 338. Potrebbe essere necessario un po' di tempo per seguire le molteplici disposizioni elicoidali di queste strutture. La glicina in posizione 338 non si trova in un'elica, ma in quella che sembra un'ansa appena prima di un'area a bassa confidenza nella struttura.
> > 3. Sulla base delle informazioni raccolte finora, potremmo formulare un'ipotesi sul motivo per cui questo è distruptivo. Non è in un'elica (di solito, nelle proteine transmembrana, le eliche sono inserite nella membrana), quindi è in uno dei domini più grandi che sporgono dalla membrana, dentro o fuori la cellula. Questa mutazione probabilmente non interrompe la struttura nei suoi segmenti intra-membrana, ma piuttosto uno dei domini funzionali. Se si desidera approfondire, è possibile verificare se si tratta del segmento extra- o intra-cellulare nel **Feature viewer**.
> > 
> {: .solution}
{: .question}

Da dove provengono le informazioni nel visualizzatore di strutture?

> <hands-on-title>Ricerca di OPN1LW umano su UniProtKB</hands-on-title>
> 
> 1. Fare clic sull'icona di download sotto la struttura
> 2. Controllare il file che è stato scaricato
> 
{: .hands-on}

Questo è un file PDB (Protein Data Bank) che consente di visualizzare e analizzare la disposizione degli atomi e degli amminoacidi della proteina.

Tuttavia, non c'è alcun riferimento al database PDB nei collegamenti tra i database delle *strutture 3D*. Invece, il primo link fa riferimento all'AlphaFoldDB. Il database AlphaFold è una risorsa completa che fornisce strutture 3D previste per un'ampia gamma di proteine. Utilizzando tecniche di deep learning e informazioni evolutive, AlphaFold predice la disposizione spaziale degli atomi all'interno di una proteina, contribuendo alla comprensione della funzione e delle interazioni proteiche.

Si tratta quindi di una *predizione* della struttura, non di una struttura convalidata sperimentalmente. Questo è il motivo per cui è colorato in base alla confidenza: le sezioni in blu sono quelle con un alto valore di confidenza, quindi quelle per cui la predizione è molto affidabile, mentre quelle in arancione sono meno affidabili o hanno una struttura disordinata (più flessibile e mobile). Tuttavia, queste informazioni sono rappresentate attraverso un file PDB, perché sono comunque strutturali.

### Famiglia e domini

La sezione *Famiglia e domini* nella pagina della voce UniProtKB fornisce una visione completa delle relazioni evolutive e dei domini funzionali di una proteina. Questa sezione offre informazioni sull'appartenenza della proteina a famiglie, superfamiglie e domini, facendo luce sulle sue caratteristiche strutturali e funzionali.

L'area *Caratteristiche* conferma effettivamente che almeno uno dei due domini che sporgono dalla membrana (quello N-terminale) è disordinato. Quest'area di solito include informazioni su regioni conservate, motivi e caratteristiche di sequenza importanti che contribuiscono al ruolo della proteina in vari processi biologici. La sezione conferma ancora una volta che siamo di fronte a una proteina transmembrana e rimanda a diverse risorse di dati filogenetici, famiglie di proteine o domini, che ci guidano nella comprensione di come le proteine condividano un'ascendenza comune, si evolvano e acquisiscano funzioni specializzate.

### Sequenza

Tutte queste informazioni sull'evoluzione, la funzione e la struttura della proteina sono codificate nella sua sequenza. Ancora una volta, in questa sezione abbiamo la possibilità di scaricare il file FASTA che la trascrive, nonché di accedere alla fonte di questi dati: gli esperimenti di sequenziamento genomico che li hanno valutati. Questa sezione riporta anche i casi in cui sono state rilevate isoforme.

> <question-title></question-title>
> 
> Quante potenziali isoforme sono mappate su questa voce?
> 
> > <solution-title></solution-title>
> > 
> > 1: H0Y622
> > 
> {: .solution}
{: .question}

### Proteine simili

L'ultima sezione della pagina UniProt Entry riporta le proteine simili (si tratta sostanzialmente del risultato di un clustering, con soglie di identità del 100%, 90% e 50%).

> <question-title></question-title>
> 
> 1. Quante proteine simili hanno un'identità del 100%?
> 2. Quante proteine simili hanno un'identità del 90%?
> 3. Quante sono le proteine simili al 50% di identità?
> 
> > <solution-title></solution-title>
> > 
> > 1. 0
> > 2. 83
> > 3. 397
> > 
> {: .solution}
{: .question}

Come avrete intuito consultando questa pagina, gran parte dell'elaborazione dei dati biologici su una proteina consiste nel mappare diversi tipi di informazioni lungo la sequenza e capire come si influenzano a vicenda. Una mappatura visiva (e una tabella con le stesse informazioni) è fornita dalle due schede alternative per visualizzare questa voce, cioè il *Variant viewer* e il *Feature viewer*.

## Visualizzazione delle varianti

> <hands-on-title>Visualizzazione delle varianti</hands-on-title>
> 
> 1. Cliccare sulla scheda *Variant viewer*
> 
{: .hands-on}

Il *Variant viewer* mappa tutte le versioni alternative note di questa sequenza. Per alcune di esse l'effetto (patogeno o benigno) è noto, per altre no.

> <question-title></question-title>
> 
> Quante varianti sono probabilmente patogene?
> 
> > <solution-title></solution-title>
> > 
> > ingrandendo la vista delle varianti, vediamo che abbiamo 5 punti rossi, quindi 5 varianti probabilmente patogene.
> > 
> {: .solution}
{: .question}

L'elevato numero di varianti che si trovano in questa sezione suggerisce che le "sequenze proteiche" (così come le sequenze geniche, le strutture proteiche ecc.) sono in realtà entità meno fisse di quanto si possa pensare.

## Visualizzazione delle caratteristiche

> <hands-on-title>Visualizzazione delle caratteristiche</hands-on-title>
> 
> 1. Fare clic sulla scheda *Feature viewer*
> 
{: .hands-on}

Il *Feature viewer* è fondamentalmente una versione fusa di tutte le aree *Features* che abbiamo trovato nella pagina *Entry*, tra cui *Domains & sites*, *Molecule processing*, *PTMs*, *Topology*, *Proteomics*, *Variants*. Se nel visualizzatore si fa clic su una qualsiasi caratteristica, la regione corrispondente nella struttura verrà messa a fuoco, come variante di interesse

> <hands-on-title>Variante visualizzatore</hands-on-title>
> 
> 1. espandere la parte *Varianti*
> 2. Zoom out
> 3. Fare clic sulla variante di interesse (il punto rosso in posizione 338)
> 
{: .hands-on}

> <question-title></question-title>
> 
> Qual è la topologia in questa posizione?
> 
> > <solution-title></solution-title>
> > 
> > Un dominio topologico citoplasmatico
> > 
> {: .solution}
{: .question}

Infine, diamo una rapida occhiata alle altre schede.

## Pubblicazioni

> <hands-on-title>Pubblicazione</hands-on-title>
> 
> 1. Fare clic sulla scheda *Pubblicazione*
> 
{: .hands-on}

La sezione *Pubblicazioni* elenca le pubblicazioni scientifiche relative alla proteina. Queste sono raccolte unendo un elenco completamente curato in UniProtKB/Swiss-Prot e un elenco importato automaticamente. In questa scheda è possibile filtrare l'elenco delle pubblicazioni in base alla fonte e alle categorie che si basano sul tipo di dati contenuti in una pubblicazione sulla proteina (come funzione, interazione, sequenza, ecc.) o sul numero di proteine nello studio corrispondente che descrive ("small scale" vs "large scale").

> <question-title></question-title>
> 
> 1. Quante pubblicazioni sono associate a questa proteina?
> 2. quante pubblicazioni contengono informazioni sulla sua funzione?
> 
> > <solution-title></solution-title>
> > 
> > 1. 57
> > 2. 23
> > 
> {: .solution}
{: .question}

## Collegamenti esterni

> <hands-on-title>Collegamenti esterni</hands-on-title>
> 
> 1. Fare clic sulla scheda *Collegamenti esterni*
> 
{: .hands-on}

La scheda *Collegamenti esterni* riunisce tutti i riferimenti a banche dati e risorse di informazione esterne che abbiamo trovato in ogni sezione della pagina di entrata. Il testo dei link riporta spesso gli identificatori unici che rappresentano la stessa entità biologica in altri database. Per avere un'idea di questa complessità, si veda l'immagine seguente (che è già parzialmente obsoleta).

![Un grafico che rappresenta come tutti i diversi database sono collegati da ID univoci, dove i nodi maggiori sono i DB e le frecce li collegano agli ID (nodi minori) che riportano. La mappa è molto affollata, soprattutto intorno a UniProt Entry Name, Gene ID e Ensembl Gene ID](./images/complexDB.jpeg "Best illustration of the complex ID crossref: bioDBnet Network Diagram - [source](https://biodbnet-abcc.ncifcrf.gov/dbInfo/netGraph.php)")

## Storia

Infine, anche la scheda *Storia* è interessante. Riporta e rende disponibili per il download tutte le versioni precedenti delle annotazioni di questa voce, cioè tutta l'"evoluzione" della sua annotazione, in questo caso risalente al 1988.

> <question-title></question-title>
> 
> Questa voce non è mai stata annotata automaticament?
> 
> > <solution-title></solution-title>
> > 
> > Per rispondere a questa domanda si può scorrere indietro nel tempo nella tabella e controllare la colonna `Database`. Questa voce è mai stata inserita in TrEMBL invece che in SwissProt? No, quindi questa voce è stata annotata manualmente fin dall'inizio.
> > 
> {: .solution}
{: .question}
