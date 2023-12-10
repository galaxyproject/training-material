---
title: GTN Stats
area: other
layout: faq
box_type: none
contributors: [hexylena]
---

<!-- tutorial stats -->
{% assign tutorials = site.pages | where:"layout", "tutorial_hands_on" | where_exp:"item","item.enable != false" %}
{% assign faqs = site.pages | where:"layout", "faq" %}
{% assign topics = site | list_topics_by_category: "science-technical" | to_keys %}
{% assign contributors = site.data['contributors'] | where_exp: "item", "item.halloffame != 'no'" | sort: "joined" %}


{% if include.compact %}
<div class="row" style="color: var(--text-color-boxtitle)">

<div class="col-md-3 col-sm-6 col-6">
 <div class="gtn-card color-comment">
   <div class="card-title-small">{{ contributors | size }}</div>
   <div class="card-text-small">Contributors</div>
 </div>
</div>

<div class="col-md-3 col-sm-6 col-6">
 <div class="gtn-card color-agenda">
   <div class="card-title-small">{{ topics | size }}</div>
   <div class="card-text-small">Topics</div>
 </div>
</div>

<div class="col-md-3 col-sm-6 col-6">
 <div class="gtn-card color-tip">
   <div class="card-title-small">{{ tutorials | size }}</div>
   <div class="card-text-small">Tutorials</div>
 </div>
</div>

<div class="col-md-3 col-sm-6 col-6">
 <div class="gtn-card color-handson">
   <div class="card-title-small">{{ site.age | round: 1}}</div>
   <div class="card-text-small">Years</div>
 </div>
</div>
</div>

{% else %}

<div class="row" style="color: var(--text-color-boxtitle)">

<div class="col-md-4">
 <div class="gtn-card color-agenda">
   <div class="card-title">{{ topics | size }}</div>
   <div class="card-text">Topics</div>
 </div>
</div>

<div class="col-md-4">
 <div class="gtn-card color-tip">
   <div class="card-title">{{ tutorials | size }}</div>
   <div class="card-text">Tutorials</div>
 </div>
</div>

<div class="col-md-4">
 <div class="gtn-card color-details">
   <div class="card-title">{{ faqs | size }}</div>
   <div class="card-text"><abbr title="Frequently Asked Questions">FAQs</abbr></div>
 </div>
</div>

<div class="col-md-6">
 <div class="gtn-card color-comment">
   <div class="card-title">{{ contributors | size }}</div>
   <div class="card-text">Contributors</div>
 </div>
</div>

<div class="col-md-6">
 <div class="gtn-card color-handson">
   <div class="card-title">{{ site.age | round: 1}}</div>
   <div class="card-text">Years</div>
 </div>
</div>
</div>
{% endif %}
