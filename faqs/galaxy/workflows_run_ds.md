---
title: Importing and Launching a Dockstore workflow
area: workflows
box_type: hands_on
layout: faq
contributors: [hexylena]
optional_parameters:
  dockstore_id: The Dockstore ID of the workflow you want to import
  title: The title of the workflow you want to import, can be any text
examples:
  Load a workflow by Dockstore ID: {dockstore_id: "github.com/jmchilton/galaxy-workflow-dockstore-example-1/mycoolworkflow", title: "My Cool Workflow"}
---

{% if include.dockstore_id %}

<div class="show-when-galaxy-proxy-active">

<span class="workflow" data-workflow="https://dockstore.org/api/ga4gh/trs/v2/tools/#workflow/{{ include.dockstore_id }}">Launch <strong>{{ include.title }}</strong> <i class="fas fa-share-alt" aria-hidden="true"></i></span>
(<a href="https://dockstore.org/workflows/{{ include.dockstore_id }}">View on Dockstore</a>)


</div>

<div class="hide-when-galaxy-proxy-active">

<a href="https://my.galaxy.training/?path=/workflows/trs_import%3ftrs_server=dockstore.org%26run_form=true%26trs_id=%2523workflow/{{ include.dockstore_id }}">Launch <strong>{{ include.title }}</strong> <i class="fas fa-share-alt" aria-hidden="true"></i></a>
(<a href="https://dockstore.org/workflows/{{ include.dockstore_id }}">View on Dockstore</a>)

</div>

{% snippet faqs/galaxy/workflows_import_from_dockstore.md override_title="If this does not work" dockstore_id=include.dockstore_id %}

{% else %}

1. Go to [Workflow → Import](https://my.galaxy.training/?path=/workflows/import) in your Galaxy
2. Switch tabs to TRS ID
3. Ensure the TRS server is set to "dockstore.org"
4. Provide your workflow hub ID

{% endif %}
