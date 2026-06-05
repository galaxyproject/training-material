---
title: Erstellen einer neuen Datei
description: Galaxy ermöglicht es Ihnen, über das Upload-Menü neue Dateien zu erstellen. Sie können den Inhalt der Datei angeben.
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
    - Tillsa
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


* Klicken Sie auf {% icon galaxy-upload %} **Daten hochladen** am oberen Rand des Werkzeugfensters
* Wählen Sie {% icon galaxy-wf-edit %} **Einfügen/Daten holen** am unteren Rand
* Fügen Sie den Inhalt der Datei in das Textfeld ein
{% if include.name %}
* Ändern Sie den Namen des Datensatzes von "Neue Datei" in `{{ include.name }}`
{% endif %}
{%- if include.format -%}
* Ändern Sie **Type** von "Auto-detect" auf `{{ include.format }}`
{%- endif -%}
{%- if include.genome -%}
* Ändern Sie **Genom** in `{{ include.genome }}`
{%- endif -%}
{%- if include.convertspaces -%}
* Wählen Sie im Menü Einstellungen ({% icon galaxy-gear %}) **Leerzeichen in Tabulatoren umwandeln**
{%- endif -%}
* Drücken Sie **Start** und **Schließen** Sie das Fenster

