---
layout: page
title: GTN Events
---

This page lists all past and upcoming GTN events.

Would you like to add your own event here? See [instructions]()

# Upcoming & Past Training Events

{% assign events = site.pages |  where: "layout", "event" %}


{% for event in events %}

event: <a href="{{site.baseurl}}{{event.url}}">{{event.title}}</a><br>

{% endfor %}


