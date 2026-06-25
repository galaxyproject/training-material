---
title: Importieren von Daten aus einer Datenbibliothek
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
    - Tillsa
    - unode
  funding:
    - biont
---



Als Alternative zum Hochladen der Daten von einer URL oder Ihrem Computer können die Dateien auch von einer *Shared Data Library* zur Verfügung gestellt werden:

1. Gehen Sie in **Bibliotheken** (linker Bereich)
2. Navigieren Sie zu {% if include.path %}: *{{ include.path }}* oder {% endif %} dem richtigen Ordner, wie von Ihrem Ausbilder angegeben.
   {% unless include.path %}- Auf den meisten Galaxies werden die Tutoriumsdaten in einem Ordner mit dem Namen **GTN - Material --> Topic Name -> Tutorial Name** bereitgestellt. {% endunless %}
3. Wählen Sie die gewünschten Dateien aus
4. Klicken Sie auf **Zur Historie hinzufügen** {% icon galaxy-dropdown %} am oberen Rand und wählen Sie **{{ include.astype | default: "as Datasets" }}** aus dem Dropdown-Menü
5. Wählen Sie im Pop-up-Fenster
   {% if include.collection_type %}
   * *"Sammlungstyp "*: **{{ include.collection_type }}**
   {% endif %}
   * *"Historie auswählen "*: {% if include.tohistory %}{{ include.tohistory }}{% else %}die Historie, in die Sie die Daten importieren möchten (oder erstellen Sie eine neue){% endif %}
6. Klicken Sie auf {% if include.collection_type %}**Fortfahren**{% else %}**Importieren**{% endif %}
{% if include.collection_type %}
7. Geben Sie der neuen Sammlung im nächsten Dialog einen passenden **Name**{% if include.collection_name %}, wie `{{ include.collection_name }}`,{% endif %}
8. Klicken Sie auf **Sammlung erstellen**
{% endif %}
