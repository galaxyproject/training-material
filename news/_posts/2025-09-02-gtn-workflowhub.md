---
title: "GTN is now integrated with WorkflowHub"
contributions:
  authorship: [supernord, hexylena, fbacall, bedroesb, elichad, dadrasarmin, shiltemann, bgruening, frederikcoppens, carolegoble]
  infrastructure: [hexylena, fbacall, bedroesb, uwwint, bgruening, frederikcoppens, carolegoble]
  funding: [eurosciencegateway, AustralianBioCommons]
tags: [gtn, workflowhub, workflows, fair]
layout: news
cover: news/images/workflowhub.png
coveralt: "WorkflowHub logo."

---

Thanks to a collaborative effort between the teams at [Galaxy Training Network (GTN)](https://training.galaxyproject.org/), [WorkflowHub](https://workflowhub.eu/), and [Australian BioCommons](https://www.biocommons.org.au/), GTN workflows are now registered automatically with WorkflowHub. The existing set of workflows can be viewed here: [**GTN on WorkflowHub**](https://workflowhub.eu/projects/12/workflows). For every new tutorial that is added to the GTN, any workflows that it contains will now also be pushed to WorkflowHub.

Workflows maintained by the [Intergalactic Workflow Commission](https://iwc.galaxyproject.org) are also automatically registered with WorkflowHub thanks to an earlier integration. Those workflows can be viewed here: [**IWC on WorkflowHub**](https://workflowhub.eu/projects/33#workflows)

## Why is it important?

The GTN encompasses a wide array of tutorials that guide researchers in the application of Galaxy for specific analyses, while also onboarding them to both the Galaxy platform itself and its vibrant community of contributors. The tutorials are often underpinned by a computational workflow that has been created to assist with the described training outcomes, and there are currently more than 330 such computational workflows. Each workflow is an investment of time, effort, and most importantly expertise. When combined with the scale, global reach, and impact of Galaxy, there is a demonstrable need for these workflows to be as findable as possible. To address this, every new tutorial that features a workflow will now be registered automatically with the WorkflowHub registry.

## About WorkflowHub

[WorkflowHub](https://workflowhub.eu/) is a public registry dedicated to sharing computational workflows across all scientific domains. It is deliberately agnostic to workflow language, management system or discipline, meaning a Galaxy workflow, a Snakemake pipeline or a custom Python script are equally welcome.

Each workflow entry contains rich metadata, version history, license information and ontology annotations. WorkflowHub supports citation infrastructure (such as DOIs and [ORCID](https://orcid.org/) IDs) so that authors receive scholarly credit.

[Galaxy’s existing integration with WorkflowHub]({% link news/2023/11/20/workflow-search.md %}) allows researchers to discover, import, and run workflows from WorkflowHub directly within Galaxy. Alternatively, from the WorkflowHub website, users can utilise the [“Run on Galaxy” button](https://galaxyproject.org/news/2023-11-13-run-in-galaxy-button-workflowhub/) to launch the workflow in Galaxy.
## Where can I learn more?

The existing set of registered GTN workflows can be [**viewed here**](https://workflowhub.eu/projects/12/workflows).

Read the new WorkflowHub publication to learn more about the registry: {% cite gustafsson2025workflowhub %}

The integration with WorkflowHub will be improved over time, enriching the metadata made discoverable through registration, and further improving the interoperability of WorkflowHub with the Galaxy ecosystem.
