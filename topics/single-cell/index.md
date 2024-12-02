---
layout: topic
topic_name: single-cell
---

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% snippet faqs/galaxy/analysis_troubleshooting.md sc=true %}

## Want to explore analysis beyond our tutorials?

<iframe src="https://training.galaxyproject.org/training-material/workflows/embed.html?query=single-cell" height="600px" width="100%" class="gtn-embed" frameborder="1"></iframe>

## News & Events - new testing
<section>
    <h2 class="mb-3">News and Events</h2>
      <iframe width="100%" height="600px" src="https://training.galaxyproject.org/training-material/feeds/single-cell-month.w.html"></iframe>
  </div>
</section>

## News & Events - ROW2

<div class="row">
  <!-- First Column: News and Events -->
  <div class="col-md-6 mb-4">
    <h2 class="mb-3">News and Events</h2>
    <iframe width="100%" height="600px" src="https://training.galaxyproject.org/training-material/feeds/single-cell-month.w.html"></iframe>
  </div>
  <!-- Second Column: How to contribute -->
<div class="col-md-6 mb-4">
    <h2 class="mt-4 mb-3">Want to Contribute?</h2>
    <p>If you want to help us behind the scenes, from testing workflows and tutorials to building tools, join our Galaxy Single-cell & Spatial Omics Community of Practice!</p>
    <ul class="contribute-list">
      <li>{% icon point-right %} <a href="https://galaxyproject.org/projects/singlecell/" target="_blank">Community of Practice</a></li>
      <li>{% icon feedback %} <a href="https://matrix.to/#/#spoc3:matrix.org" target="_blank">Matrix Chat Forum</a></li>
      <li>{% icon email %} <a href="https://lists.galaxyproject.org/lists/single-cell-cop.lists.galaxyproject.org/" target="_blank">Mailing List</a></li>
    </ul>
  </div>
</div>

## Try event horizon
<a class="Single-cell & sPatial Omics Event Horizon" href="{% if event.external %}{{event.external}}{% else %}{{site.baseurl}}{{event.url}}{% endif %}{% if include.campaign %}?utm_source=gtn&utm_medium=event-table&utm_campaign={{ include.campaign }}{% endif %}">{{event.title}}{% if event.draft %} (draft, will be hidden) {% endif %}</a>

## Try event horizon 2
{% assign upcoming_events = site | get_upcoming_events_for_this: new_material %}
{% assign upcoming_event_count = upcoming_events | size %}
{% if upcoming_event_count > 0 %}
<blockquote class="details hide-when-printing" id="upcoming-events">
              <div id="upcoming-events-c" class="box-title">
                <button type="button" aria-controls="upcoming-events-c" aria-expanded="false">
                  <i class="fas fa-calendar" aria-hidden="true"></i>
                  Upcoming events which cover this topic<span role="button" class="fold-unfold fa fa-minus-square" aria-hidden="true"></span>
                </button>
               </div>

   <p>Want to learn more with a live instructor? Check out these upcoming events:</p>
   {% include _includes/event-table.html events=upcoming_events campaign="via-tutorial" %}
</blockquote>
{% endif %}

## Want to contribute?

If you want to help us behind the scenes, from testing workflows and tutorials to building tools, join our Galaxy Single-cell & sPatial Omics Community of Practice!

 - {% icon point-right %}  [Community of Practice](https://galaxyproject.org/projects/singlecell/)
 - {% icon feedback %}  [Matrix Chat Forum](https://matrix.to/#/#spoc3:matrix.org)
 - {% icon email %}  [Mailing List](https://lists.galaxyproject.org/lists/single-cell-cop.lists.galaxyproject.org/)
