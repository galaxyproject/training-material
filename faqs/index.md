---
layout: page
---


# Frequently Asked Questions

Thanks for using the GTN! If you have any questions regarding the GTN, Galaxy, or any of the topics covered by the tutorials, you can check one of the FAQ pages below to see if your question is among them.


## General FAQs

Common questions about the Galaxy platform, or about the GTN itself, can be found on the following pages:

<div markdown="0">

{% include _includes/button.html link="/faqs/galaxy" colour="yellow" label="Galaxy FAQs" buttonsize="large" %}

{% include _includes/button.html link="/faqs/gtn" colour="purple" label="GTN FAQs" buttonsize="large" %}

</div>


## Topic FAQ pages

Each topic in the GTN also has a dedicated FAQ page:

{% assign sorted_topics = site | list_topics: "all" %}

### Galaxy for Scientists

<div markdown="0">
{% assign buttoncolours = "red,orange,yellow,green,blue,purple,pink" | split: "," %}

{% assign buttonnum = 0 %}
{% for topic in sorted_topics %}
  {% if topic[1].enable != false and topic[1].type == 'use'  %}
    {% assign l = topic[1].name | prepend: "/topics/" | append: "/faqs/" %}
    {% assign lab = topic[1].title %}
    {% assign i = buttonnum | modulo: buttoncolours.size %}
    {% assign col = buttoncolours[i] %}

    {% include _includes/button.html link=l label=lab buttonsize="large" colour=col %}

    {% assign buttonnum = buttonnum | plus: 1 %}
  {% endif %}
{% endfor %}

</div>

### Galaxy for Instructors, Developers & System Administrators

<div markdown="0">
{% for topic in sorted_topics %}
  {% if topic[1].enable != false and topic[1].type and topic[1].type != 'use'  %}

    {% assign l = topic[1].name | prepend: "/topics/" | append: "/faqs/" %}
    {% assign lab = topic[1].title %}
    {% assign i = buttonnum | modulo: buttoncolours.size %}
    {% assign col = buttoncolours[i] %}

    {% include _includes/button.html link=l label=lab buttonsize="large" colour=col %}

    {% assign buttonnum = buttonnum | plus: 1 %}

  {% endif %}
{% endfor %}
</div>


## Tutorial FAQ pages

Tutorials may also have FAQ pages associated with them, you can find the link to these in the overview box in the tutorial page itself.


## Ask the Community

Couldn't find the answer to your question here? Please feel free to ask in the [GTN Gitter channel]({{site.gitter_url}}), or on the [Galaxy Help Forum](https://help.galaxyproject.org) as well.
