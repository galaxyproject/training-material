---
title: Kopieren eines Datensatzes zwischen Historien
description: Manchmal möchten Sie ein Dataset in mehreren Histories verwenden. Sie müssen die Daten nicht erneut hochladen, sondern können Datasets von einer History in eine andere kopieren.
area: histories
box_type: tip
layout: faq
contributions:
  authorship:
    - lecorguille
    - shiltemann
    - hexylena
    - bebatut
    - lldelisle
  translation:
    - Tillsa
    - unode
  funding:
    - biont
---


Es gibt 3 Möglichkeiten, Datensätze zwischen Historien zu kopieren

1. Aus der ursprünglichen Historie

    1. Klicken Sie auf das Symbol {% icon galaxy-gear %}, das sich am Anfang der Liste der Datensätze im Verlaufsfenster befindet
    2. Klicken Sie auf **Datensätze kopieren**
    3. Wählen Sie die gewünschten Dateien
    {% if include.history_name %}
    4. "Neuer Verlauf mit dem Namen:" `{{ include.history_name }}`
    {% else %}
    4. Geben Sie der "Neuen Historie" einen passenden Namen
    {% endif %}
    5. Validieren durch 'Kopieren von Verlaufsdaten'
    5. Klicken Sie auf den neuen Namen der Historie in dem grünen Feld, das gerade erschienen ist, um zu dieser Historie zu wechseln

2. Verwendung des {% icon galaxy-columns %} **Historien nebeneinander anzeigen**

    1. Klicken Sie auf den {% icon galaxy-dropdown %} Dropdown-Pfeil oben rechts im Verlaufsfenster (**History options**)
    2. Klicken Sie auf {% icon galaxy-columns %} **Historien nebeneinander anzeigen**
    3. Wenn Ihr Zielverlauf nicht vorhanden ist
        1. Klicken Sie auf 'Historien auswählen'
        2. Klicken Sie auf Ihren Zielverlauf
        3. Validieren durch 'Change Selected'
    3. Ziehen Sie den zu kopierenden Datensatz aus seiner ursprünglichen Historie
    4. Legen Sie ihn in der Zielhistorie ab

3. Aus der Zielhistorie

    1. Klicken Sie auf **Benutzer** in der oberen Leiste
    2. Klicken Sie auf **Datensätze**
    3. Suche nach dem zu kopierenden Datensatz
    4. Klicken Sie auf den Namen
    5. Klicken Sie auf **Kopieren in die aktuelle Historie**

