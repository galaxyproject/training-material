---
layout: page
---


# Frequently Asked Questions

Thanks for using the GTN! If you have any questions regarding the GTN, Galaxy, or any of the topics covered by the tutorials, you can check one of the FAQ pages below to see if your question is among them.


## Quick start

<a href="{{site.baseurl}}/faqs/galaxy/"><button type="button" class="btn btn-warning btn-info">Galaxy FAQs</button></a>
<a href="{{site.baseurl}}/faqs/gtn/"><button type="button" class="btn btn-warning btn-info">GTN FAQs</button></a>

## GTN Questions

Have questions about this website or the GTN in general? Please check the [GTN FAQs]({{site.baseurl}}/faqs/gtn/)

## Galaxy Questions

Have a question about Galaxy? Please check the [Galaxy FAQs]({{site.baseurl}}/faqs/galaxy/)


## Topic FAQ pages

Each topic in the GTN also has a dedicated FAQ page:

{% assign sorted_topics = "" | split: "," %}
{% assign sorted_topics_pre = site.data | sort | order: "title" %}

{% for topic in sorted_topics_pre %}
    {% if topic[0] == "introduction" %}
        {% assign sorted_topics = sorted_topics | unshift: topic %}
    {% else %}
        {% assign sorted_topics = sorted_topics | push: topic %}
    {% endif %}
{% endfor %}


**Galaxy for Scientists**

<ul>
{% for topic in sorted_topics %}
  {% if topic[1].enable != false and topic[1].type == 'use'  %}
<li><a href="{{site.baseurl}}/topics/{{topic[1].name}}/faqs/"> {{topic[1].title}}</a></li>
  {% endif %}
{% endfor %}
</ul>

**Galaxy for Instructors, Developers & System Administrators**

<ul>
{% for topic in sorted_topics %}
  {% if topic[1].enable != false and topic[1].type and topic[1].type != 'use'  %}
<li><a href="{{site.baseurl}}/topics/{{topic[1].name}}/faqs/"> {{topic[1].title}}</a></li>
  {% endif %}
{% endfor %}
</ul>


## Tutorial FAQ pages

Tutorials may also have FAQ pages associated with them, you can find the link to these in the overview box in the tutorial page itself.


## Ask the Community

Couldn't find the answer to your question here? Please feel free to ask in the [GTN Gitter channel]({{site.gitter_url}}), or on the [Galaxy Help Forum](https://help.galaxyproject.org) as well.
