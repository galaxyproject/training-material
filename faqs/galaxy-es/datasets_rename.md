---
title: Cambiar el nombre de un conjunto de datos
area: datasets
box_type: tip
layout: faq
contributions:
  authorship:
    - bebatut
    - nsoranzo
    - shiltemann
    - wm75
    - hexylena
  translation:
    - hvelab
    - unode
  funding:
    - biont
optional_parameters:
  name: The new name for the dataset
examples:
  Without prescribing a name: {}
  Suggest the user rename the dataset genomes.fa:
    name: genomes.fa
---


- Haga clic en el {% icon galaxy-pencil %} **icono lápiz** del conjunto de datos para editar sus atributos
- En el panel central, cambie el campo **Name** {% if include.name %} a `{{ include.name }}` {% endif %}
- Haga clic en el botón **Save**

