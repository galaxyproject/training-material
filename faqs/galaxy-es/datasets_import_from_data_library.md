---
title: Importar datos de una biblioteca de datos
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
    - hvelab
    - unode
  funding:
    - biont
---



Como alternativa a cargar los datos desde una URL o desde su ordenador, los archivos también pueden estar disponibles desde una *biblioteca de datos compartidos*:

1. Entra en **Libraries** (panel izquierdo)
2. Navega a {% if include.path %}: *{{ include.path }}* o {% endif %} la carpeta correcta indicada por su instructor. {% unless include.path %}- En la mayoría de los tutoriales de Galaxias los datos serán proporcionados en una carpeta llamada **GTN - Material --> Topic name -> Tutorial name**. {% endunless %}
3. Seleccione los archivos deseados
4. Haz clic en **Add to History** {% icon galaxy-dropdown %} cerca de la parte superior y selecciona **{{ include.astype | default: "as Datasets" }}** en el menú desplegable
5. En la ventana emergente, elige {% if include.collection_type %}
   * *"Tipo de colección "*: **{{ include.collection_type }}** {% endif %}
   * *"Seleccionar historial "*: {% if include.tohistory %}{{ include.tohistory }}{% else %}el historial al que desea importar los datos (o crear uno nuevo){% endif %}
6. Haga clic en {% if include.collection_type %}**Continue**{% else %}**Import**{% endif %} {% if include.collection_type %}
7. En el siguiente cuadro de diálogo, asigne un **Name**{% if include.collection_name %} adecuado, como `{{ include.collection_name }}`,{% endif %} a la nueva colección
8. Haz clic en **Create collection** {% endif %}

