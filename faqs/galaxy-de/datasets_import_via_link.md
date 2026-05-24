---
title: Importieren über Links
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
    - Tillsa
    - unode
  funding:
    - biont
optional_parameters:
  reset_form: Den Benutzer zuerst anweisen, das Formular zurückzusetzen
  collection: Fügt einen Schritt hinzu, um die Schaltfläche für Sammlungen zu klicken
  collection_type: Schlägt einen Sammlungstyp vor
  link: (Veraltet) Statt Sammlungen einen Link‑Eingabetyp verwenden
  link2: (Veraltet) Zweites Linkformat
  genome: Schlägt ein Genom vor, das der Benutzer setzen sollte
  pairswaptext: Bietet einen Filter für die Paar‑Erkennung
  collection_name_convention: Schlägt eine Namenskonvention für die Sammlung vor
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


* Kopieren der Linkposition
* Klicken Sie auf {% icon galaxy-upload %} **Daten hochladen** am oberen Rand der Werkzeugleiste
{% if include.reset_form %}
* Klicken Sie auf die Schaltfläche **Zurücksetzen** am unteren Rand des Formulars. Wenn die Schaltfläche ausgegraut ist -> zum nächsten Schritt übergehen.
{% endif %}
{% if include.collection %}
* Klicken Sie oben auf **Sammlung**
{% endif %}
{% if include.collection_type %}
* Klicken Sie auf **Sammlungstyp** und wählen Sie `{{ include.collection_type }}`
{% endif %}
* Wählen Sie {% icon galaxy-wf-edit %} **Daten einfügen/holen**
* Fügen Sie den/die Link(s) in das Textfeld ein
{% if include.link %}
  `{{ include.link }}`
{% endif %}
{% if include.link2 %}
  `{{ include.link2 }}`
{% endif %}
{% if include.format %}
* Ändern Sie **Type (set all):** von "Auto-detect" auf `{{ include.format }}`
{% endif %}
{% if include.genome %}
* **Genom** auf `{{ include.genome }}` ändern
{% endif %}
* Drücken Sie **Start**
{% if include.collection %}
* Klicken Sie auf **Build** wenn verfügbar
{% if include.pairswaptext %}
* Stellen Sie sicher, dass die Vorwärts- und Rückwärtslesungen jeweils auf {{ include.pairswaptext }} gesetzt sind.
    * Klick auf **Tauschen** sonst
{% endif %}
* Geben Sie einen Namen für die Sammlung ein
{% if include.collection_name_convention %}
    * Eine nützliche Namenskonvention ist die Verwendung von {{ include.collection_name_convention }}
{% endif %}
{% if include.collection_name %}
    * {{ include.collection_name }}
{% endif %}
* Klicken Sie auf **Liste erstellen** (und warten Sie ein wenig)
{% else %}
* **Schließen** Sie das Fenster
{% endif %}
