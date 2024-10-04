---
title: Input Histories & Answer Keys
area: gtn
box_type: tip
layout: faq
contributors: [hexylena, nomadscientist]
---

Tutorials sometimes require significant amounts of data or data prepared in a very specific manner which often is shown to cause errors for learners that significantly affect downstream results. Input histories are an answer to that:

{% assign kid_key = "Input Histories Schema" %}
{% assign kid_val = site.data['schema-tutorial']['mapping']['input_histories'] %}
{% include _includes/schema-render.html key=kid_key value=kid_val %}

Additionally once the learner has gotten started, tutorials sometimes feature tools which produce stochastic outputs, or have very long-running steps. In these cases, the tutorial authors may provide answer histories to help learners verify that they are on the right track, or to enable them to catch up if they fall behind or something goes wrong.

{% assign kid_key = "Answer Histories Schema" %}
{% assign kid_val = site.data['schema-tutorial']['mapping']['answer_histories'] %}
{% include _includes/schema-render.html key=kid_key value=kid_val %}
