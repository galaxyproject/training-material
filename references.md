---
layout: base
title: GTN Top Cited Papers
---

This page lists the top cited papers, based on citations used throughout the GTN.

<ol>
{% for cite in topcitations %}
    <li>{{ cite[1]['count'] }} - {{cite[1]['text']}}</li>
{% endfor %}
</ol>
