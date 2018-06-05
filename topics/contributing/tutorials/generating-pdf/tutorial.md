---
layout: tutorial_hands_on
topic_name: contributing
tutorial_name: generating-pdf
---

# Introduction
{:.no_toc}

The website with the training material can be run locally. Sometimes, it is also interesting to freeze the tutorials or to get PDFs of the tutorials.

> ### Agenda
>
> In this tutorial, you will learn how to run a local instance of the GTN website:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Generate PDFs artifact

To generate the PDFs, a command `make pdf` is given. This command:

- Launches a detached Jekyll server to serve the website
- Generates the PDFs of the tutorials by calling Chrome via command line
- Generates the PDFs of the slide decks by calling decktage

> ### {% icon hands_on %} Hands-on: Checking the website generation locally
>
> 1. (If not done) Activate the conda environment: `source activate galaxy_training_material`
> 2. Install Chrome
>    - For OSX, install the [Chrome browser]()
>    - For Ubuntu, follow [these instructions]()
> 1. Generate the PDFs: `make pdf`
> 2. Check the generated PDFs in `_pdf` folder
{: .hands_on}

# Conclusion
{:.no_toc}

> ### Developing GTN training material
>
> This tutorial is part of a series to develop GTN training material, feel free to also look at:
>
> {% assign topic = site.data[page.topic_name] %}
> {% for material in topic.material %}
>  {% if material.enable != "false" and material.name != page.tutorial_name %}
>   {% if material.type == "introduction" %}
> 1. [{{ material.title }}]({{ site.baseurl }}/topics/{{ topic.name }}/slides/{{ material.name }}.html)
>  {% elsif material.type == "tutorial" %}
>   {% if material.hands_on %}
> 1. [{{ material.title }}]({{ site.baseurl }}/topics/{{ topic.name  }}/tutorials/{{ material.name }}/tutorial.html)
>   {% elsif material.slides %}
> 1. [{{ material.title }}]({{ site.baseurl }}/topics/{{ topic.name }}/tutorials/{{ material.name }}/slides.html)
>    {% endif %}
>   {% endif %}
>  {% endif %}
> {% endfor %}
{: .agenda}
