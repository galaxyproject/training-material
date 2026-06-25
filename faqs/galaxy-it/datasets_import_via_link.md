---
title: Importazione tramite link
area: data upload
box_type: tip
layout: faq
contributions:
  authorship:
    - bebatut
    - hexylena
    - mtekman
    - lecorguille
    - shiltemann
    - nomadscientist
    - nekrut
    - mblue9
  translation:
    - lisanna
    - silviadg87
    - hvelab
    - unode
  funding:
    - biont
optional_parameters:
  reset_form: Suggerisci all’utente di resettare prima il modulo
  collection: Aggiunge il passo per cliccare il pulsante raccolta
  collection_type: Suggerisce il tipo di raccolta
  link: (Deprecato) Al posto delle raccolte usa un input di tipo link
  link2: (Deprecato) Secondo formato di link suggerito da impostare
  genome: Suggerisce un genome da impostare per l’utente
  pairswaptext: Fornisce un filtro per il rilevamento delle coppie
  collection_name_convention: Suggerisce una convenzione di denominazione per la raccolta
examples:
  Import a list of files:
    collection: 'true'
    collection_type: List
    collection_name: raw_files
    format: thermo.raw
  Create a paired collection with a given naming scheme:
    collection: 'true'
    collection_type: Paired
    collection_name_convention: '''&lt;name&gt;_&lt;plate&gt;_&lt;batch&gt;'' to preserve
      the sample names, sequencing plate number and batch number.'
    collection_name: Here we will write 'C57_P1_B1'
    genome: GRCm38/mm10
    pairswaptext: '''SRR5683689_1'' and ''SRR5683689_2'''
  Import an interval file with specified genome, though providing the links in a code block is preferred:
    link: https://zenodo.org/record/3254917/files/chr21.gtf
    format: interval
    genome: mm9
---


* Copia la posizione del collegamento
* Fare clic su {% icon galaxy-upload %} **Carica i dati** nella parte superiore del pannello degli strumenti {% if include.reset_form %}
* Fare clic sul pulsante **Reset** in fondo al modulo. Se il pulsante è grigio -> passare al passo successivo. {% endif %} {% if include.collection %}
* Fare clic su **Collection** in alto {% endif %} {% if include.collection_type %}
* Fare clic su **Tipo di raccolta** e selezionare `{{ include.collection_type }}` {% endif %}
* Selezionare {% icon galaxy-wf-edit %} **Incollare/recuperare i dati**
* Incollare il/i link nel campo di testo {% if include.link %} `{{ include.link }}` {% endif %} {% if include.link2 %} `{{ include.link2 }}` {% endif %} {% if include.format %}
* Cambiare **Type (set all):** da "Auto-detect" a `{{ include.format }}` {% endif %} {% if include.genome %}
* Cambiare **Genome** in `{{ include.genome }}` {% endif %}
* Premere **Avvio** {% if include.collection %}
* Fare clic su **Costruisci** quando disponibile {% if include.pairswaptext %}
* Assicurarsi che le letture avanti e indietro siano impostate rispettivamente su {{ include.pairswaptext }}.
    * Fare clic su **Swap** altrimenti {% endif %}
* Inserire un nome per la collezione {% if include.collection_name_convention %}
    * Una convenzione di denominazione utile è usare {{ include.collection_name_convention }} {% endif %} {% if include.collection_name %}
    * {{ include.collection_name }} {% endif %}
* Fare clic su **Crea elenco** (e attendere un po') {% else %}
* **Chiude** la finestra {% endif %}

