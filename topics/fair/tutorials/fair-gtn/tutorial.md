---
layout: tutorial_hands_on
title: FAIR Galaxy Training Material
abbreviations:
  FAIR: Findable, Accessible, Interoperable, Reusable
  GTN: Galaxy Training Network
zenodo_link: ''
questions:
- What are the FAIR training materials?
- How to test, reproduce and share your content?
- How to collaborate and don’t duplicate?
objectives:
- Learn about metadata and findability
- Learn how to support system and content curation
time_estimation: "30M"
key_points:
- FAIR principles in Galaxy training development and content creation.
tags:
- fair
- gtn
- training
priority: 3
contributions:
  authorship:
    - kkamieniecka
    - poterlowicz-lab
  editing:
    - hexylena
  funding:
      - ELIXIR-UK-DaSH
subtopic: fair-data

requirements:
  - type: "internal"
    topic_name: fair
    tutorials:
      - fair-intro
---


Encouraging computational reproducibility in research, we will present a variety of data stewardship recommendations that we have found useful in the process of training development. As part of that process, we are exploring the application of the FAIR (Findable, Accessible, Interoperable, Reusable) guidelines to the Galaxy Training Network (GTN) materials, in order to improve their secondary use and adaptation.

This tutorial outlines how to set and use existing resources to make Galaxy training development and content creation FAIR.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

Here, we refer to a set of good practices as described in *"Ten simple rules for making training materials FAIR"* {% cite Garcia2020 %}.

![Ten simple rules for making training materials FAIR.]({% link topics/fair/images/fair_gtn.png %} "Ten simple rules for making training materials FAIR. The primary guideline is to share; the Findability rules are description, identity, and registration; the latter two, along with access rules, correspond to accessibility; and the first rule is to share; With the exception of the format rule, which stands alone for interoperability, the remaining four criteria all relate to various facets of reusability")

Image credit: Luc Wiegers and Celia van Gelder {% cite wiegers_luc_2019_3593258 %}

## Plan to share your training materials

The Galaxy Training Network (GTN) provides researchers with online training materials, connects them with local trainers, and helps promoting open data analysis practices worldwide. It provides a record of training content development and ensure materials curation via  cross-domain repositories (i.g., [GitHub](https://github.com/galaxyproject/training-material), [Zenodo](https://zenodo.org/)). Keeping your materials in the right places from the beginning will make it possible for you to more effectively and widely distribute your work without any duplication in contribution. Instruction where to start creating a new tutorial can be found at [GTN contributing tutorial]({% link topics/contributing/tutorials/create-new-tutorial/tutorial.md %}).

## Improve findability of your training materials by properly describing them

Describing your training materials with structured metadata is fundamental to making them FAIR resource. Creating metadata boost findability and preserve information and can be updated after publication. GTN  tutorial require several mandatory metadata information such as learning objectives following Bloom’s taxonomy {% cite chevron2014metacognitive %}, prerequisites, time estimate, or questions addressed by the tutorial. Content developers needs to know where to find all of the available metadata to reference it later. [Schema.org](https://schema.org/) is a collaborative project with a mission to create the addition of structured metadata to web pages. It describes data using shared vocabulary. It is in this spirit that the [BioSchemas](https://bioschemas.org/) community initiative was created to extend the Schema.org standard to life-science resources with training content specification. GTN offers schemas command line [tutorial]({% link topics/contributing/tutorials/schemas/tutorial.md %}) where metadata can be added or updated. This reduces the complexity of keywords and set the scene for controlled vocabularies environment. This increases the effectiveness of metadata filtering, which decreases ambiguity and makes it easier to find and retrieve information.

## Give your training materials a unique identity

To support teachers and trainers, the GTN tutorials rely on specific tags and identifiers. Adding persistent identifiers (PIDs) to training materials makes them easy to cite and aids citation counting in research metric systems {% cite mcmurry2017identifiers %}. Providing PIDs ensure that tag/identifier will continue to work even if webpage location changes. Data used for tutorials are required to be stored at [Zenodo](https://zenodo.org/) and associated with [Digital Object Identifiers](https://www.doi.org/) (DOI).

## Register your training materials online
Sharing and publishing with the GTN training helps minimize the amount of time and effort required for instructors to prepare for and run their training courses and workshops, by providing templates and a complete training infrastructure. Online registry makes your content more discoverable and accessible to wider community. GTN tutorials are automatically registered on [TeSS](https://tess.elixir-europe.org/) (ELIXIR Training e-support system).

## Define access rules for your training materials
The GTN materials are designed to be flexible in their use, open access and community-driven. It is important to follow the hosting website disclaimers and keep materials metadata in place.

## Use an interoperable format for your training materials
GTN Tutorials aim to follow best practices in course design, so that they can be used in different environments. Content of the tutorials and slides are written in [Markdown]({% link topics/contributing/tutorials/create-new-tutorial-content/tutorial.md %}) and supported by templates. Metadata are stored in YAML and workflows in JSON. Self-contained structure allow others to tune or reuse existing content.

## Make your training materials (re)usable for trainers
To help others determine whether the training materials are relevant and adaptable to their particular situations, metadata published alongside training materials should include context and sufficient description including: contributor details, license, description, learning outcomes, audience, requirements, tags/keywords, duration and last revision date Applying the proper licence and tagging training materials with metadata can make it simpler for others to (re)use and adapt them.

GTN provide strong technical support and set of [contributing self learning material]({% link topics/contributing/index.md %}).

## Make your training materials usable for trainees
Prerequisites and learning outcomes are particularly helpful metadata. To be effective, learning objectives must be written using active verbs that describe the expected trainee behaviours as well as the knowledge, skills, and expertise they will have received. Rich metadata requirements and SMART—Specific (Measurable, Attainable, Relevant, and Time-bound) learning outcomes following Bloom's taxonomy {% cite chevron2014metacognitive %} helps to clarify which trainees will benefit most from the training. Self-learning structure of the GTN materials supported by slides and video walkthroughs/tours adds another layer of usability.

## Make your training materials contribution friendly
GTN have clear [guidelines]({% link topics/contributing/tutorials/create-new-tutorial/tutorial.md %}) for contribution and involvement. Community supported CONTRIBUTING files and chat offer the chance to exchange information such as contact details, expectations for contributions, and more. All contributors should be listed and thanked in the acknowledgements; how to cite the tutorial and give credit to contributors can be found at the end of each tutorial.

## Keep your training materials up-to-date
It is crucial to keep your training materials up-to-date and so you are aware of any new features, trends, or improvements in the topic (such as updated tools or databases). Transparent peer-review and curation process in collaborative and open set up ensure materials quality.

# Conclusion
There are many ways to improve the training content. It is crucial that we work together to make training materials FAIR so that everyone can benefit from them. These simple suggestions are intended to encourage dialogue and cooperation among Galaxy training community and help bring the latest developments to users.

The Galaxy Training Network is an example of a robust, effective Community of Practice.
For more information please look at this great article {% cite hiltemann2023galaxy %}, the corresponding FAIR guidelines {% cite fair-training-materials %} and follow [short introduction to FAIR data stewardship](http://fellowship.elixiruknode.org/).
