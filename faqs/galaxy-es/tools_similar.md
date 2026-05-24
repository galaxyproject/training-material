---
title: Múltiples herramientas similares disponibles
area: tools
box_type: tip
layout: faq
contributions:
  authorship:
    - shiltemann
  translation:
    - hvelab
    - unode
  funding:
    - biont
---


A veces hay varias herramientas con nombres muy parecidos. Si los parámetros del tutorial no coinciden con lo que ves en Galaxy, prueba lo siguiente:

1. Utiliza el **Modo Tutorial** {% icon curriculum %} en Galaxy, y haz click en el botón azul de la herramienta en el tutorial para abrir automáticamente la herramienta y versión correctas (aún no disponible para todos los tutoriales)

   {% snippet faqs/galaxy-es/tutorial_mode.md %}

2. Comprueba que el **nombre completo de la herramienta** coincide con lo que ves en el tutorial. {% if include.toolname %} Compruebe que:
   - **Nombre completo de la herramienta**: `{{ include.toolname }}` {% endif %} {% if include.toolversion %}
   - **Versión de la herramienta**: `{{include.toolversion}}` (escrita después del nombre de la herramienta) {% endif %}

