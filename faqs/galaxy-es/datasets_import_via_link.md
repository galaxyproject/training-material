---
title:  Importación por medio de enlaces
area: data upload
box_type: tip
layout: faq
contributors: [nomadscientist]
---

* Copia los enlaces
* Abre el manejador de carga de datos de Galaxy ({% icon galaxy-upload %} en la parte superior derecha del panel de herramientas)
{% if include.reset_form %}
* Selecciona **Reset** butón en la parte inferior del forma. Si el butón es gris, muevete a la próxima instrucción.
{% endif %}
{% if include.collection %}
* Selecciona **Collection** del encima
{% endif %}
{% if include.collection_type %}
* Selecciona **Collection Type** y selecciona `{{ include.collection_type }}`
{% endif %}
* Selecciona ‘Pegar/Traer datos’ **Paste/Fetch Data**
* Copia los enlaces en el campo de textos
{% if include.link %}
  `{{ include.link }}`
{% endif %}
{% if include.link2 %}
  `{{ include.link2 }}`
{% endif %}
{% if include.format %}
* Cambia **Type (set all):** de "Auto-detect" a `{{ include.format }}`
{% endif %}
{% if include.genome %}
* Cambia **Genome** a `{{ include.genome }}`
{% endif %}
* Presiona ‘Iniciar’ **Start**
{% if include.collection %}
* Selecciona **Build** cuando está disponible
{% if include.pairswaptext %}
* Asegúrate de que ambas lecturas, sentido y antisentido, tengan la configuración {{ include.pairswaptext }}, respectivamente.
* Ensure that the forward and reverse reads are set to {{ include.pairswaptext }}, respectively.
    * Si no, selecciona **Swap**
{% endif %}
* Pon un nombre por la colección
{% if include.collection_name_convention %}
    * Es útil usar la próxima manera de nombrar archivos: {{ include.collection_name_convention }}
{% endif %}
{% if include.collection_name %}
    * {{ include.collection_name }}
{% endif %}
* Selecciona **Create list** (y espera un momento)
{% else %}
* **Close** Cierra la ventana.
{% endif %}
{% if include.renaming == undefined or include.renaming == true %}
* Galaxy utiliza los URLs como nombres de forma predeterminada , así que los tendrás que cambiar a algunos que sean más útiles o informativos.
the window
{% endif %}
