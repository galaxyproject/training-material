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
{% assign topics = site.data | where_exp: "item", "item.type" | where_exp:"item","item.enable != false" %}
{% assign contributors = site.data['contributors'] | where_exp: "item", "item.halloffame != 'no'" | sort: "joined" %}

<div class="row">
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
