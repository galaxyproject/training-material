---
title: Importing and Launching a Dockstore Workflow
area: workflows
box_type: hands_on
layout: faq
contributors: [hexylena]
optional_parameters:
  dockstore_id: The Dockstore ID of the workflow you want to import, minus the "#workflow/" portion
  title: The title of the workflow you want to import, can be any text
  version: Version of the workflow you want to pin
examples:
  Load a workflow by Dockstore ID: {dockstore_id: "github.com/jmchilton/galaxy-workflow-dockstore-example-1/mycoolworkflow", title: "My Cool Workflow"}
  Load a v0.1.3 of a specific IWC workflow: {dockstore_id: "github.com/iwc-workflows/kmer-profiling-hifi-VGP1/main", title: "Kmer Profiling HiFi VGP1", version: "v0.1.3"}
---

<!-- GTN:IGNORE:011 _blank is used here due to iframeing sites that (can) set x-frame-options. The contexts in which _blank is set are only visible when iframe'd -->

{% if include.dockstore_id %}
{% capture external_page %}https://dockstore.org/workflows/{{ include.dockstore_id }}:{{ include.version }}{% endcapture %}
<div class="show-when-galaxy-proxy-active">


<a class="workflow" target="_top" href="/workflows/trs_import?trs_server=dockstore.org&run_form=true&trs_id=%23workflow/{{ include.dockstore_id }}&trs_version={{ include.version }}">
  Launch <strong>{{ include.title }} ({{ include.version }})</strong> <i class="fas fa-share-alt" aria-hidden="true"></i>
</a>
(<a target="_blank" href="{{ external_page }}">View on Dockstore</a>)

</div>

<div class="hide-when-galaxy-proxy-active">

<a href="https://my.galaxy.training/?path=/workflows/trs_import%3ftrs_server=dockstore.org%26run_form=true%26trs_id=%2523workflow/{{ include.dockstore_id }}%26trs_version={{ include.version }}">
  Launch <strong>{{ include.title }} ({{ include.version }})</strong> <i class="fas fa-share-alt" aria-hidden="true"></i>
</a>
(<a href="{{ external_page }}">View on Dockstore</a>)

</div>

{% snippet faqs/galaxy/workflows_import_from_dockstore.md override_title="If this does not work" dockstore_id=include.dockstore_id %}

{% else %}

1. Go to [Workflow → Import](https://my.galaxy.training/?path=/workflows/import) in your Galaxy
2. Switch tabs to TRS ID
3. Ensure the TRS server is set to "dockstore.org"
4. Provide your workflow hub ID

{% endif %}
