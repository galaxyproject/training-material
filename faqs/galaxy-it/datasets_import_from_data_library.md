---
title: Importare i dati da una libreria di dati
area: data upload
box_type: tip
layout: faq
contributions:
  authorship:
    - bebatut
    - shiltemann
    - nsoranzo
    - hexylena
    - wm75
  translation:
    - lisanna
    - silviadg87
    - hvelab
    - unode
  funding:
    - biont
---



In alternativa al caricamento dei dati da un URL o dal proprio computer, i file possono essere resi disponibili da una libreria di dati condivisi:

1. Entrare in **Librerie** (pannello sinistro)
2. Navigare verso {% if include.path %}: *{{ include.path }}* o {% endif %} alla cartella corretta indicata dal vostro istruttore. {% unless include.path %} Nella maggior parte dei Galaxies i dati delle esercitazioni vengono forniti in una cartella denominata **GTN - Materiale --> Nome argomento -> Nome esercitazione**. {% endunless %}
3. selezionare i file desiderati
4. Fare clic su **Aggiungi alla cronologia** {% icon galaxy-dropdown %} vicino alla parte superiore e selezionare **{{ include.astype | default: "as Datasets" }}** dal menu a tendina
5. Nella finestra pop-up, scegliere {% if include.collection_type %}
   * *"Tipo di raccolta "*: **{{ include.collection_type }}** {% endif %}
   * *"Seleziona cronologia "*: {% if include.tohistory %}{{ include.tohistory }}{% else %} la cronologia in cui si desidera importare i dati (o crearne una nuova){% endif %}
6. Cliccare su {% if include.collection_type %}**Continue**{% else %}**Import**{% endif %} {% if include.collection_type %}
7. Nella finestra di dialogo successiva, dare un **nome**{% if include.collection_name %} adeguato, come `{{ include.collection_name }}`,{% endif %} alla nuova raccolta
8. Cliccare su **Crea raccolta** {% endif %}

