---
title: Creare un nuovo file
description: Galaxy permette di creare nuovi file dal menu di caricamento. È possibile fornire direttamente il contenuto del file.
area: datasets
box_type: tip
layout: faq
contributions:
  authorship:
    - bebatut
    - shiltemann
    - wm75
    - hexylena
    - simonbray
  translation:
    - lisanna
    - silviadg87
    - hvelab
    - unode
  funding:
    - biont
optional_parameters:
  name: The name of the dataset
  format: The format of the dataset
  genome: The genome of the dataset
  convertspaces: Ask the user to convert spaces to tabs in the settings
examples:
  Creating a new file: {}
  Creating a specific file:
    name: SARS-CoV-2 feature mapping
    format: tabular
    convertspaces: 'true'
---


* Fare clic su {% icon galaxy-upload %} **Caricamento dati** nella parte superiore del pannello degli strumenti
* Selezionare {% icon galaxy-wf-edit %} **Incolla/Recupera dati** in basso
* Incollare il contenuto del file nel campo di testo {% if include.name %}
* Cambiare il nome del dataset da "New File" a `{{ include.name }}` {% endif %} {%- if include.format -%}
* Cambiare **Type** da "Auto-detect" a `{{ include.format }}` {%- endif -%} {%- if include.genome -%}
* Cambiare **Genome** in `{{ include.genome }}` {%- endif -%} {%- if include.convertspaces -%}
* Dal menu Impostazioni ({% icon galaxy-gear %}) selezionare **Converti spazi in tabulazioni** {%- endif -%}
* Premere **Avvio** e **Chiudere** la finestra

