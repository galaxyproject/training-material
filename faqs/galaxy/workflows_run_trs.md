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
workflow.
</div>

<div class="hide-when-galaxy-proxy-active">

Click to 
<a href="https://my.galaxy.training/?path=/workflows/trs_import%3Frun_form=true%26trs_url={{ site.url }}{{ site.baseurl }}{{ include.path | convert_workflow_path_to_trs }}">
    Launch <strong>{{ include.title }}</strong> {% icon workflow aria=false %}
</a>

</div>

{% else %}

1. Find the material you are interested in
2. View it's workflows, which can be found in the metadata box at the top of the tutorial
3. Click the button on any workflow to run it.

{% endif %}
