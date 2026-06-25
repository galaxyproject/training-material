---
title: Hinzufügen eines Tags
description: Tags können dir helfen, deinen Verlauf besser zu organisieren und Datensätze nachzuverfolgen.
area: datasets
layout: faq
box_type: tip
contributions:
  authorship:
    - bebatut
    - wm75
    - hexylena
    - shiltemann
    - nekrut
  translation:
    - Tillsa
    - unode
  funding:
    - biont
---


Datensätze können getaggt werden. Dies vereinfacht die Verfolgung von Datensätzen über die Galaxy-Schnittstelle. Tags können eine beliebige Kombination von Buchstaben oder Zahlen enthalten, dürfen aber keine Leerzeichen enthalten.

**Um einen Datensatz zu taggen**:

1. Klicken Sie auf den Datensatz, um ihn zu erweitern
2. Klicken Sie auf **Add Tags** {% icon galaxy-tags %}
3. Fügen Sie {% if include.tag %} ein Tag namens `{{include.tag}}`{% else %} tag text{% endif %} hinzu. Tags, die mit `#` beginnen, werden automatisch an die Ausgaben von Tools weitergegeben, die diesen Datensatz verwenden (siehe unten).
4. Drücken Sie <kbd>Enter</kbd>
5. Stellen Sie sicher, dass das Tag unter dem Namen des Datensatzes erscheint

**Tags, die mit `#` beginnen, sind speziell!**

Sie werden **Name-Tags** genannt. Das Besondere an diesen Tags ist, dass sie sich *vererben*: Wenn ein Datensatz mit einem Name-Tag gekennzeichnet ist, erben alle Ableitungen (Kinder) dieses Datensatzes automatisch dieses Tag (siehe unten). Die folgende Abbildung erklärt, warum dies so nützlich ist. Betrachten Sie die folgende Analyse (die Zahlen in Klammern entsprechen den Datensatznummern in der Abbildung unten):

1. Ein Satz von Vorwärts- und Rückwärts-Reads (Datensätze 1 und 2) wird mit {% tool Bowtie2 %} gegen eine Referenz gemappt, wodurch Datensatz 3 erzeugt wird;
1. Datensatz 3 wird verwendet, um die Read Coverage mit {% tool BedTools Genome Coverage %} *getrennt* für `+` und `-` Stränge. Dies erzeugt zwei Datensätze (4 und 5 für plus bzw. minus);
1. Die Datensätze 4 und 5 werden als Input für die {% tool Macs2 broadCall %}-Datensätze verwendet, die die Datensätze 6 und 8 erzeugen;
1. Die Datensätze 6 und 8 werden mit {% tool BedTools Intersect %} mit den Koordinaten von Genen (Datensatz 9) geschnitten, wodurch die Datensätze 10 und 11 entstehen.

![Eine Geschichte ohne Namens-Tags versus Geschichte mit Namens-Tags]({% link shared/images/histories_why_nametags.svg %})

Betrachten wir nun, dass diese Analyse ohne Namens-Tags durchgeführt wird. Dies ist auf der linken Seite der Abbildung zu sehen. Es ist schwer nachzuvollziehen, welche Datensätze "Plus"-Daten und welche "Minus"-Daten enthalten. Enthält zum Beispiel Datensatz 10 "Plus"-Daten oder "Minus"-Daten? Wahrscheinlich "minus", aber sind Sie sicher? Bei einem kleinen Verlauf wie dem hier gezeigten ist es möglich, dies manuell nachzuvollziehen, aber wenn der Umfang eines Verlaufs wächst, wird es sehr schwierig.

Die rechte Seite der Abbildung zeigt genau die gleiche Analyse, jedoch unter Verwendung von Namens-Tags. Als die Analyse durchgeführt wurde, waren die Datensätze 4 und 5 mit `#plus` bzw. `#minus` gekennzeichnet. Als sie als Eingaben in Macs2 verwendet wurden, übernahmen die Datensätze 6 und 8 automatisch diese Markierungen und so weiter... Daher ist es einfach, beide Zweige (plus und minus) dieser Analyse zu verfolgen.

Weitere Informationen finden Sie in einem [dedizierten #nametag-Tutorial]({% link topics/galaxy-interface/tutorials/name-tags/tutorial.md %}).


<!-- Image is here = https://docs.google.com/drawings/d/1iiNsau6ddiE2MV9qMyekUq2mrpDHHcc02bXtcFEAnhY/edit?usp=sharing -->


