---
layout: page
title: GTN Event Horizon
---

![Logo with a woman standing in front of a globe with multiple location pins in it]({% link assets/images/events.svg %}){: style="float:right; width:25%;"}

This page lists all past and upcoming GTN-related events organised by members of the community.

For a complete list of all Galaxy-related events, please also see the [Galaxy Event Horizon](https://galaxyproject.org/events/)

Anybody is welcome to add their events here.

<a href="TODO"><button type="button" class="btn btn-success">Add your own event!</button></a>

{% assign upcoming_events = site.pages| where_exp: "item", "item.layout == 'event' or item.layout == 'event-external' " | where_exp: "item", "item.event_state == 'upcoming' " | sort: 'date_start' | reverse %}
{% assign ongoing_events = site.pages | where_exp: "item", "item.layout == 'event' or item.layout == 'event-external' " | where_exp: "item", "item.event_state == 'ongoing'  " | sort: 'date_start' | reverse %}
{% assign past_events = site.pages    | where_exp: "item", "item.layout == 'event' or item.layout == 'event-external' " | where_exp: "item", "item.event_state == 'ended'    " | sort: 'date_start' | reverse %}
{% assign ongoing_length = ongoing_events | size %}

{% if ongoing_length >0 %}

# Current Events

{% include _includes/event-table.html events=ongoing_events %}

{% endif %}

# Upcoming Events

{% include _includes/event-table.html events=upcoming_events %}

# Past Events

{% include _includes/event-table.html events=past_events %}
