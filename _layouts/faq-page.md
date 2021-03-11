---
layout: faqs
---



{% assign faqs =  site.pages | where: "layout", "faq"  %}
{% assign topic_faqs = faqs | where_exp: "item", "item.path contains page.snippets_dir "%}
{% assign topic_faqs_sorted = topic_faqs | sort: "area" %}


{% assign area_prev = '' %}

{% for q in topic_faqs_sorted %}
{% unless q.area == area_prev %}
{% unless forloop.first %} <br><br> {% endunless %}
<h2> {{ q.area | capitalize}}</h2>
<hr/>
{% assign area_prev = q.area %}
{% endunless %}

<h3 class="faq-area"> {{q.title}} <a href="{{site.baseurl}}{{q.url}}">{% icon galaxy_instance %}</a> </h3>
{{q.description}}
{% snippet {{q.path}} %}



{% endfor %}

