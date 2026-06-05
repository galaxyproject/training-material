---
title: Erstellen einer gepaarten Sammlung
area: collections
box_type: tip
layout: faq
contributions:
  authorship:
    - shiltemann
    - hexylena
    - lldelisle
  translation:
    - Tillsa
    - unode
  funding:
    - biont
---



* Klicken Sie auf {% icon galaxy-selector %} **Elemente auswählen** am oberen Rand des Verlaufsfensters ![Schaltfläche Elemente auswählen]({% link topics/galaxy-interface/images/historyItemControls.png %})

* Überprüfen Sie {% if include.datasets_description %}{{ include.datasets_description }}{% else %}alle Datensätze in Ihrem Verlauf, die Sie einschließen möchten{% endif %}
* Klicken Sie auf **{% if include.n %}{{ include.n }}{% else %}n{% endif %} of N selected** und wählen Sie **Liste der Datensatzpaare erstellen**

  ![Menüpunkt 'Paare-Liste erstellen']({% link topics/galaxy-interface/images/buildList.png %}){:width="15%"}

* Sie befinden sich im Assistenten zum Erstellen von Sammlungen. Wählen Sie **Liste gepaarter Datensätze** und klicken Sie auf die Schaltfläche 'Weiter' unten rechts.

  ![Assistent zum Erstellen von Sammlungen – gepaarte Liste]({% link topics/galaxy-interface/images/paired_list_selector.png %}){:width="15%"}

* Überprüfen und konfigurieren Sie das automatische Paaren. Gewöhnlich haben Mate-Paare die Endungen `_1 ` und `_2` oder `_R1` und `_R2`. Klicken Sie unten auf 'Weiter'.

  ![Bearbeiten und Erstellen einer Liste gepaarter Sammlungen]({% link topics/galaxy-interface/images/paired_list_edit.png %}){:width="15%"}

* Bearbeiten Sie den **Listenbezeichner** nach Bedarf.
* Geben Sie einen Namen für Ihre Sammlung ein
* Klicken Sie auf **Erstellen**, um Ihre Sammlung zu erstellen
* Klicken Sie erneut auf das Häkchen-Symbol oben in Ihrem Verlauf
