---
title: Creación de un nuevo fichero
description: Galaxy permite crear nuevos ficheros desde el menú de carga. Puedes proporcionar el contenido del fichero.
area: datasets
box_type: tip
layout: faq
contributions:
  authorship:
    - bebatut
    - shiltemann
    - wm75
    - hexylena
    - simonbray
  translation:
    - hvelab
    - unode
  funding:
    - biont
optional_parameters:
  name: The name of the dataset
  format: The format of the dataset
  genome: The genome of the dataset
  convertspaces: Ask the user to convert spaces to tabs in the settings
examples:
  Creating a new file: {}
  Creating a specific file:
    name: SARS-CoV-2 feature mapping
    format: tabular
    convertspaces: 'true'
---


* Pulse {% icon galaxy-upload %} **Cargar Datos** en la parte superior del panel de herramientas
* Seleccione {% icon galaxy-wf-edit %} **Pegar/Obtener Datos** en la parte inferior
* Pegue el contenido del archivo en el campo de texto {% if include.name %}
* Cambie el nombre del conjunto de datos de "New File" a `{{ include.name }}` {% endif %} {%- if include.format -%}
* Cambie **Type** de "Auto-detect" a `{{ include.format }}` {%- endif -%} {%- if include.genome -%}
* Cambiar **Genome** a `{{ include.genome }}` {%- endif -%} {%- if include.convertspaces -%}
* En el menú Configuración ({% icon galaxy-gear %}) seleccione **Convertir espacios en tabuladores** {%- endif -%}
* Pulse **Iniciar** y **Cerrar** la ventana

