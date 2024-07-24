---
layout: page
title: GTN Event Horizon
description: |
  This page lists all past and upcoming GTN-related events organised by members of the community.
---

![Logo with a woman standing in front of a globe with multiple location pins in it]({% link assets/images/events.svg %}){: style="float:right; width:25%;"}

This page lists all past and upcoming GTN-related events organised by members of the community.

For a complete list of all Galaxy-related events, please also see the [Galaxy Event Horizon](https://galaxyproject.org/events/) and [ELIXIR's TeSS](https://tess.elixir-europe.org/events?q=galaxy).

Anybody is welcome to add their events here.

[Add your own event!]({% link faqs/gtn/gtn_event_create.md %}){: .btn.btn-success}

{% assign upcoming_events = site.pages| where_exp: "item", "item.layout == 'event' or item.layout == 'event-external' " | where_exp: "item", "item.event_state == 'upcoming' " | sort: 'date_start'  %}
{% assign ongoing_events = site.pages | where_exp: "item", "item.layout == 'event' or item.layout == 'event-external' " | where_exp: "item", "item.event_state == 'ongoing'  " | sort: 'date_start' | reverse %}
{% assign past_events = site.pages    | where_exp: "item", "item.layout == 'event' or item.layout == 'event-external' " | where_exp: "item", "item.event_state == 'ended'    " | sort: 'date_start' | reverse %}
{% assign upcoming_length = upcoming_events | size %}
{% assign ongoing_length = ongoing_events | size %}

{% if ongoing_length >0 %}

# Current Events

{% include _includes/event-table.html events=ongoing_events %}

{% endif %}

# Upcoming Events

{% if upcoming_length == 0 %}
  <p>No known upcoming events. Check out <a href="https://tess.elixir-europe.org/events?q=galaxy">ELIXIR's TeSS</a> or the <a href="https://galaxyproject.org/events/">Hub</a> for other Galaxy events going on!</p>
{% else %}
  {% include _includes/event-table.html events=upcoming_events %}
{% endif %}

# Past Events

{% include _includes/event-table.html events=past_events %}
