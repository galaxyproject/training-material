---
layout: base
title: GTN Editorial Board
---

{% assign contributors = site.data['contributors'] %}
{% assign sorted_topics = site | list_topics_ids %}

<h1> GTN Editorial Board </h1>
<section>
 <p class="lead">
  A big <em>Thank You!</em> to the GTN Editorial Board members who ensure the quality of the tutotials in their respective topics.
 </p>

{% for t in sorted_topics %}
 {% assign topic = site.data[t] %}
 <h3> {{ topic.title }}</h3>

 {% include _includes/contributor-list.html contributors=topic.editorial_board badge=true %}
 {% unless topic.tag_based %}
 <p>View the <a href="{{site.baseurl}}/topics/{{topic.name}}/maintainer.html" >Topic maintainer page</a> for this topic.</p>
 {% endunless %}
{% endfor %}

</section>

<h2> Codeowners</h2>

Automatic CODEOWNERS file (just to make it easier to keep updated)


<!-- autogenerate a codeowners file based on this info -->

```
assets/     @bebatut @shiltemann @hexylena
bin/        @bebatut @shiltemann @hexylena
metadata/   @bebatut @shiltemann @hexylena
badges/     @hexylena
_layouts/   @bebatut @shiltemann @hexylena
_includes/  @bebatut @shiltemann @hexylena
_plugins/   @bebatut @hexylena

{% for t in sorted_topics %}{% assign topic = site.data[t] %}
topics/{{t}}/               {% for member in topic.editorial_board %} @{{member}} {% endfor %}{% endfor%}
```


