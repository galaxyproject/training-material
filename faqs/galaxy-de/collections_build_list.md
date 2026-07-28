---
title: Erstellen einer Datensatzsammlung
area: collections
box_type: tip
layout: faq
contributions:
  authorship:
    - shiltemann
    - hexylena
  translation:
    - Tillsa
    - unode
  funding:
    - biont
---


* Klicken Sie auf {% icon galaxy-selector %} **Elemente auswählen** am oberen Rand des Verlaufsfensters ![Schaltfläche Elemente auswählen]({% link topics/galaxy-interface/images/historyItemControls.png %})
* Überprüfen Sie {% if include.datasets_description %}{{ include.datasets_description }}{% else %}alle Datensätze in Ihrer Historie, die Sie aufnehmen möchten{% endif %}
* Klicken Sie auf **{% if include.n %}{{ include.n }}{% else %}n{% endif %} of N selected** und wählen Sie **Build Dataset List**

  ![Menüeintrag Liste erstellen]({% link topics/galaxy-interface/images/buildList.png %}){:width="15%"}

* You are in collection building wizard. Choose **Flat List** and click 'Next' button at the right bottom corner.

  ![Assistent zum Erstellen von Sammlungen – Flache Liste]({% link topics/galaxy-interface/images/flat_list_selector.png %}){:width="15%"}

* Doppelklicken Sie auf die Dateinamen, um sie zu bearbeiten. Entfernen Sie z. B. Dateiendungen oder gemeinsame Präfixe/Suffixe, um die Namen zu bereinigen.

  ![Liste bearbeiten und erstellen]({% link topics/galaxy-interface/images/flat_list_edit.png %}){:width="15%"}

* Geben Sie einen Namen für Ihre Sammlung ein {% if include.name %} bis {{ include.name }} {% endif %}
* Klicken Sie auf **Sammlung erstellen**, um Ihre Sammlung zu erstellen
* Klicken Sie erneut auf das Häkchensymbol oben in Ihrem Verlauf
