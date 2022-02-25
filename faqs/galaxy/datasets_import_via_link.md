---
title: Importing via links
area: data upload
box_type: tip
layout: faq
contributors: [bebatut,hexylena,mtekman,lecorguille,shiltemann,nomadscientist,nekrut,mblue9]
---

* Copy the link location
* Open the Galaxy Upload Manager ({% icon galaxy-upload %} on the top-right of the tool panel)
{% if include.reset_form %}
* Click **Reset** button at the bottom of the form. If the button is greyed out -> skip to the next step.
{% endif %}
{% if include.collection %}
* Click on **Collection** on the top
{% endif %}
{% if include.collection_type %}
* Click on **Collection Type** and select `{{ include.collection_type }}`
{% endif %}
* Select **Paste/Fetch Data**
* Paste the link into the text field
{% if include.link %}
  `{{ include.link }}`
{% endif %}
{% if include.link2 %}
  `{{ include.link2 }}`
{% endif %}
{% if include.format %}
* Change **Type (set all):** from "Auto-detect" to `{{ include.format }}`
{% endif %}
{% if include.genome %}
* Change **Genome** to `{{ include.genome }}`
{% endif %}
* Press **Start**
{% if include.collection %}
* Click on **Build** when available
{% if include.pairswaptext %}
* Ensure that the forward and reverse reads are set to {{ include.pairswaptext }}, respectively.
    * Click **Swap** otherwise
{% endif %}
* Enter a name for the collection
{% if include.collection_name_convention %}
    * A useful naming convention is to use {{ include.collection_name_convention }}
{% endif %}
{% if include.collection_name %}
    * {{ include.collection_name }}
{% endif %}
* Click on **Create list** (and wait a bit)
{% else %}
* **Close** the window
{% endif %}
