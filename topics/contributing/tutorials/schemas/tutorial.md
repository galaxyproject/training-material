---
layout: tutorial_hands_on
questions:
- What metadata is required or possible to set in a Tutorial, Slide, Topic, or FAQ
objectives:
- Know where to find all of the available metadata, so you can reference it later.
title: "GTN Metadata"
time_estimation: "10m"
contributors:
  - hexylena
---

## Training Materials

{% assign kid_key = "Tutorial Schema" %}
{% assign kid_val = site.data['schema-tutorial'] %}
{% include _includes/schema-render.html key=kid_key value=kid_val %}

{% assign kid_key = "Slides Schema" %}
{% assign kid_val = site.data['schema-slides'] %}
{% include _includes/schema-render.html key=kid_key value=kid_val %}

## GTN Resources

{% assign kid_key = "FAQ Schema" %}
{% assign kid_val = site.data['schema-faq'] %}
{% include _includes/schema-render.html key=kid_key value=kid_val %}

{% assign kid_key = "Topic Schema" %}
{% assign kid_val = site.data['schema-topic'] %}
{% include _includes/schema-render.html key=kid_key value=kid_val %}

{% assign kid_key = "Quiz Schema" %}
{% assign kid_val = site.data['schema-quiz'] %}
{% include _includes/schema-render.html key=kid_key value=kid_val %}

{% assign kid_key = "Learning Pathway Schema" %}
{% assign kid_val = site.data['schema-learning-pathway'] %}
{% include _includes/schema-render.html key=kid_key value=kid_val %}

{% assign kid_key = "News" %}
{% assign kid_val = site.data['schema-news'] %}
{% include _includes/schema-render.html key=kid_key value=kid_val %}

## Contributors

{% assign kid_key = "Contributor Schema" %}
{% assign kid_val = site.data['schema-contributors'] %}
{% include _includes/schema-render.html key=kid_key value=kid_val %}

{% assign kid_key = "Organisation Schema" %}
{% assign kid_val = site.data['schema-organisations'] %}
{% include _includes/schema-render.html key=kid_key value=kid_val %}

{% assign kid_key = "Funder Schema" %}
{% assign kid_val = site.data['schema-funders'] %}
{% include _includes/schema-render.html key=kid_key value=kid_val %}
