---
title: Importing and launching a GTN workflow
area: workflows
box_type: hands_on
layout: faq
contributors: [hexylena]
---

{% if include.path %}

<div class="show-when-galaxy-proxy-active">
<span class="workflow" data-workflow="{{ site.url }}{{ site.baseurl }}{{ include.path | convert_workflow_path_to_trs }}">Launch <strong>{{ include.title }}</strong> <i class="fas fa-share-alt" aria-hidden="true"></i></span>
(<a href="https://github.com/galaxyproject/training-material/blob/main/{{ include.path }}">View on GitHub</a>, <a href="https://training.galaxyproject.org/training-material/{{ include.path }}">Download workflow</a>)
workflow.
</div>

<div class="hide-when-galaxy-proxy-active">

Click to 
<a href="https://my.galaxy.training/?path=/workflows/trs_import%3Frun_form=true%26trs_url={{ site.url }}{{ site.baseurl }}{{ include.path | convert_workflow_path_to_trs }}">
    Launch <strong>{{ include.title }}</strong> {% icon workflow aria=false %}
</a>
(<a href="https://github.com/galaxyproject/training-material/blob/main/{{ include.path }}">View on GitHub</a>, <a href="https://training.galaxyproject.org/training-material/{{ include.path }}">Download workflow</a>)

</div>

{% capture import_url %}https://training.galaxyproject.org/training-material/{{ include.path }}{% endcapture %}
{% snippet faqs/galaxy/workflows_import.md override_title="If this does not work" import_url=import_url %}

{% else %}

1. Find the material you are interested in
2. View it's workflows, which can be found in the metadata box at the top of the tutorial
3. Click the button on any workflow to run it.

{% endif %}
