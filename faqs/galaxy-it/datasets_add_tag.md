---
title: Aggiunta di un tag
description: I tag possono aiutarti a organizzare meglio la cronologia e a tenere traccia dei dataset.
area: datasets
layout: faq
box_type: tip
contributions:
  authorship:
    - bebatut
    - wm75
    - hexylena
    - shiltemann
    - nekrut
  translation:
    - lisanna
    - silviadg87
    - unode
  funding:
    - biont
---


I dataset possono essere etichettati. Questo semplifica il tracciamento dei dataset nell'interfaccia di Galaxy. I tag possono contenere qualsiasi combinazione di lettere o numeri, ma non possono contenere spazi.

**Per etichettare un set di dati**:

1. fare clic sul set di dati per espanderlo
2. Cliccare su **Aggiungi tag** {% icon galaxy-tags %}
3. Aggiungere {% if include.tag %} un tag chiamato `{{include.tag}}`{% else %} tag text{% endif %}. I tag che iniziano con `#` saranno automaticamente propagati agli output degli strumenti che utilizzano questo set di dati (vedi sotto).
4. Premere <kbd>Invio</kbd>
5. verificare che il tag appaia sotto il nome del set di dati

**I tag che iniziano con `#` sono speciali!

Sono chiamati **Name tags**. La caratteristica unica di questi tag è che si *propagano*: se un set di dati è etichettato con un tag name, tutti i derivati (figli) di questo set di dati erediteranno automaticamente questo tag (vedi sotto). La figura seguente spiega perché questo è così utile. Si consideri la seguente analisi (i numeri tra parentesi corrispondono ai numeri dei set di dati nella figura sottostante):

1. un insieme di letture forward e reverse (set di dati 1 e 2) viene mappato rispetto a un riferimento utilizzando {% tool Bowtie2 %} generando il set di dati 3;
1. il dataset 3 è usato per calcolare la copertura delle letture usando {% tool BedTools Genome Coverage %} *separatamente* per i filamenti `+` e `-`. Questo genera due set di dati (4 e 5 per il più e il meno, rispettivamente);
1. i set di dati 4 e 5 sono utilizzati come input per i set di dati {% tool Macs2 broadCall %} che generano i set di dati 6 e 8;
1. gli insiemi di dati 6 e 8 vengono intersecati con le coordinate dei geni (insiemi di dati 9) usando {% tool BedTools Intersect %} generando gli insiemi di dati 10 e 11.

![Una storia senza name tag contro una storia con name tag]({% link shared/images/histories_why_nametags.svg %})

Ora si consideri che questa analisi è stata fatta senza tag dei nomi. Questo è mostrato sul lato sinistro della figura. È difficile individuare quali set di dati contengono dati "più" e quali "meno". Ad esempio, il set di dati 10 contiene dati "positivi" o "negativi"? Probabilmente "meno", ma ne siete sicuri? Nel caso di una storia di piccole dimensioni, come quella mostrata qui, è possibile tracciarla manualmente, ma con l'aumentare delle dimensioni di una storia diventa molto impegnativo.

La parte destra della figura mostra esattamente la stessa analisi, ma utilizzando i tag dei nomi. Quando è stata condotta l'analisi, i dataset 4 e 5 erano etichettati rispettivamente con `#plus` e `#minus`. Quando sono stati utilizzati come input per Macs2, i dataset 6 e 8 li hanno ereditati automaticamente e così via... Di conseguenza, è facile tracciare entrambi i rami (più e meno) di questa analisi.

Maggiori informazioni sono contenute in un [tutorial dedicato ai #nametag]({% link topics/galaxy-interface/tutorials/name-tags/tutorial.md %}).


<!-- Image is here = https://docs.google.com/drawings/d/1iiNsau6ddiE2MV9qMyekUq2mrpDHHcc02bXtcFEAnhY/edit?usp=sharing -->


