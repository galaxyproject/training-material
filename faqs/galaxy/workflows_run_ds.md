---
title: Importing and Launching a Dockstore workflow
area: workflows
box_type: hands_on
layout: faq
contributors: [hexylena]
---

{% if include.dockstore_id %}

<div class="show-when-galaxy-proxy-active">

Click here to run <span class="workflow" data-workflow="https://dockstore.org/api/ga4gh/trs/v2/tools/#{{ include.dockstore_id }}"><strong>{{ include.title }}</strong> <i class="fas fa-share-alt" aria-hidden="true"></i></span>

</div>

<div class="hide-when-galaxy-proxy-active">

Click here to run <a href="https://my.galaxy.training/?path=/workflows/trs_import%3ftrs_server=dockstore.org%26run_form=true%26trs_id=%2523{{ include.dockstore_id }}"><strong>{{ include.title }}</strong> <i class="fas fa-share-alt" aria-hidden="true"></i></a>


</div>

{% else %}

1. Go to [Workflow → Import](https://my.galaxy.training/?path=/workflows/import) in your Galaxy
2. Switch tabs to TRS ID
3. Ensure the TRS server is set to "Dockstore.org"
4. Provide your workflow hub ID

{% endif %}
