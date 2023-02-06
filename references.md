---
layout: base
title: GTN Top Cited Papers
---

{% assign topcitations = site.citation_count | top_citations %}
<pre>
{{ site.citation_count | jsonify }}
</pre>

<ol>
{% for cite in topcitations %}
    <li>{{ cite[1]['count'] }} - {{cite[1]['text']}}</li>
{% endfor %}
</ol>
