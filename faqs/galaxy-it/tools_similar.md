---
title: Strumenti simili alla compilazione multipla disponibili
area: tools
box_type: tip
layout: faq
contributions:
  authorship:
    - shiltemann
  translation:
    - lisanna
    - silviadg87
    - unode
  funding:
    - biont
---


A volte ci sono più strumenti con nomi molto simili. Se i parametri indicati nel tutorial non corrispondono a quelli visualizzati in Galaxy, provare con i seguenti:

1. Usare la modalità **Tutorial** {% icon curriculum %} in Galaxy e fare clic sul pulsante blu dello strumento nel tutorial per aprire automaticamente lo strumento e la versione corretti (non ancora disponibile per tutti i tutorial)

   {% snippet faqs/galaxy-it/tutorial_mode.md %}

2. Verificare che il **nome completo dello strumento** corrisponda a quello che si vede nel tutorial. {% if include.toolname %} Verificare che:
   - **Nome completo dello strumento**: `{{ include.toolname }}` {% endif %} {% if include.toolversion %}
   - **Versione dello strumento**: `{{include.toolversion}}` (scritto dopo il nome dello strumento) {% endif %}

