---
title: Umbenennen eines Datensatzes
area: datasets
box_type: tip
layout: faq
contributions:
  authorship:
    - bebatut
    - nsoranzo
    - shiltemann
    - wm75
    - hexylena
  translation:
    - Tillsa
    - unode
  funding:
    - biont
optional_parameters:
  name: Der neue Name für den Datensatz
examples:
  Without prescribing a name: {}
  Suggest the user rename the dataset genomes.fa:
    name: genomes.fa
---

- Klicken Sie auf das {% icon galaxy-pencil %} **Bleistift-Symbol** für den Datensatz, um seine Attribute zu bearbeiten
- Ändern Sie im zentralen Bereich das Feld **Name** {% if include.name %} in `{{ include.name }}` {% endif %}
- Klicken Sie auf die Schaltfläche **Speichern**
