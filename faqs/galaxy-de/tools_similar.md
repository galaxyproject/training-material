---
title: Multipile ähnliche Werkzeuge verfügbar
area: tools
box_type: tip
layout: faq
contributions:
  authorship:
    - shiltemann
  translation:
    - Tillsa
    - unode
  funding:
    - biont
---


Manchmal gibt es mehrere Werkzeuge mit sehr ähnlichen Namen. Wenn die Parameter im Tutorial nicht mit dem übereinstimmen, was Sie in Galaxy sehen, versuchen Sie bitte das Folgende:

1. Verwenden Sie den **Tutorial Mode** {% icon curriculum %} in Galaxy, und klicken Sie auf den blauen Tool-Button im Tutorial, um automatisch das richtige Tool und die richtige Version zu öffnen (noch nicht für alle Tutorials verfügbar)

   {% snippet faqs/galaxy-de/tutorial_mode.md %}

2. Überprüfen Sie, ob der **gesamte Werkzeugname** mit dem im Tutorial übereinstimmt.
   {% if include.toolname %} Bitte überprüfen Sie das:
   - **Vollständiger Werkzeugname**: `{{ include.toolname }}`
   {% endif %}
   {% if include.toolversion %}
   - **Werkzeugversion**: `{{include.toolversion}}` (wird hinter den Werkzeugnamen geschrieben)
   {% endif %}

