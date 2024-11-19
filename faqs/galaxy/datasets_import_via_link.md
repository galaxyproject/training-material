---
title: Importing via links
area: data upload
box_type: tip
layout: faq
contributors: [bebatut,hexylena,mtekman,lecorguille,shiltemann,nomadscientist,nekrut,mblue9]
optional_parameters:
  reset_form: Tell the user to reset the form first
  collection: Adds step to click the collection button
  collection_type: Suggests a collection type
  link: (Deprecated) Instead of collections use a link-type input
  link2: (Deprecated) Second link
  format: Suggests a format for the user to set
  genome: Suggests a genome for the user to set
  pairswaptext: Provides a filter for the pair detection
  collection_name_convention: Suggests a naming convention for the collection
examples:
  Import a list of files:
    collection: 'true'
    collection_type: List
    collection_name: raw_files
    format: thermo.raw
  Create a paired collection with a given naming scheme:
    collection: 'true'
    collection_type: Paired
    collection_name_convention: "'&lt;name&gt;_&lt;plate&gt;_&lt;batch&gt;' to preserve the sample names, sequencing plate number and batch number."
    collection_name: "Here we will write 'C57_P1_B1'"
    genome: "GRCm38/mm10"
    pairswaptext: "'SRR5683689_1' and 'SRR5683689_2'"
  "Import an interval file with specified genome, though providing the links in a code block is preferred":
    link: "https://zenodo.org/record/3254917/files/chr21.gtf"
    format: interval
    genome: mm9
---

* Copy the link location
* Click {% icon galaxy-upload %} **Upload Data** at the top of the tool panel
{% if include.reset_form %}
* Click **Reset** button at the bottom of the form. If the button is greyed out -> skip to the next step.
{% endif %}
{% if include.collection %}
* Click on **Collection** on the top
{% endif %}
{% if include.collection_type %}
* Click on **Collection Type** and select `{{ include.collection_type }}`
{% endif %}
* Select {% icon galaxy-wf-edit %} **Paste/Fetch Data**
* Paste the link(s) into the text field
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
