---
title: Crear una colección de conjuntos de datos
area: collections
box_type: tip
layout: faq
contributions:
  authorship:
    - shiltemann
    - hexylena
  translation:
    - hvelab
    - unode
  funding:
    - biont
---


* Haz clic en {% icon galaxy-selector %} **Seleccionar elementos** en la parte superior del panel del historial ![Botón Seleccionar elementos]({% link topics/galaxy-interface/images/historyItemControls.png %})
* Marque {% if include.datasets_description %}{{ include.datasets_description }}{% else %}todos los conjuntos de datos de su historial que desee incluir{% endif %}
* Haga clic en **{% if include.n %}{{ include.n }}{% else %}n{% endif %} of N selected** y elija **Build dataset list**

  ![elemento de menú crear colección de listas]({% link topics/galaxy-interface/images/buildList.png %}){:width="15%"}

* Introduzca un nombre para su colección {% if include.name %} to {{ include.name }} {% endif %}
* Haz clic en **Create collection** para crear tu colección
* Vuelve a hacer clic en el icono de la marca de verificación situado en la parte superior del historial


