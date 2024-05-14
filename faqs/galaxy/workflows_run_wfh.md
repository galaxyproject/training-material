---
title: Importing and Launching a WorkflowHub.eu Workflow
area: workflows
box_type: hands_on
layout: faq
contributors: [hexylena]
---

{% if include.wfhub_id %}

<div class="show-when-galaxy-proxy-active">

<span class="workflow" data-workflow="https://workflowhub.eu/ga4gh/trs/v2/tools/{{ include.wfhub_id }}/versions/{{ include.wfhub_version | default: 1 }}">Launch <strong>{{ include.title }}</strong> <i class="fas fa-share-alt" aria-hidden="true"></i></span>
(<a href="https://workflowhub.eu/workflows/{{ include.wfhub_id }}?version={{ include.wfhub_version | default: 1 }}">View on WorkflowHub</a>)

</div>

<div class="hide-when-galaxy-proxy-active">

<a href="https://my.galaxy.training/?path=/workflows/trs_import%3ftrs_server=workflowhub.eu%26run_form=true%26trs_id={{ include.wfhub_id }}%26trs_version={{ include.wfhub_version | default: 1}}">Launch <strong>{{ include.title }}</strong> <i class="fas fa-share-alt" aria-hidden="true"></i></a>
(<a href="https://workflowhub.eu/workflows/{{ include.wfhub_id }}?version={{ include.wfhub_version | default: 1 }}">View on WorkflowHub</a>)

</div>

{% else %}

1. Go to [Workflow → Import](https://my.galaxy.training/?path=/workflows/import) in your Galaxy
2. Switch tabs to TRS ID
3. Ensure the TRS server is set to "workflowhub.eu"
4. Provide your workflow hub ID

{% endif %}
