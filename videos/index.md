---
layout: base
---

{% include _includes/default-header.html %}

{% assign sorted_topics = "" | split: "," %}
{% assign sorted_topics_pre = site.data | sort | order: "title" %}

{% for topic in sorted_topics_pre %}
    {% if topic[0] == "introduction" %}
        {% assign sorted_topics = sorted_topics | unshift: topic %}
    {% else %}
        {% assign sorted_topics = sorted_topics | push: topic %}
    {% endif %}
{% endfor %}


<div class="container main-content">
<section>

<h1> GTN Video Slides </h1>

The GTN now generates videos for selected slide decks. Click on a topic below to jump to the video page for that topic!

<br/><br/>

<table class="table table-striped">
 <thead>
  <tr><th>Topic</th></tr>
 </thead>
 <tbody>
  {% for topic in sorted_topics %} {% unless topic[0] == 'contributors' %}
  <tr><td><a href="{{ site.baseurl }}/topics/{{ topic[1].name }}/videos/">{{ topic[1].title }}</a></td></tr>
  {% endunless %}{% endfor %}
 </tbody>
</table>

</section>
</div>
{% include _includes/default-footer.html %}
