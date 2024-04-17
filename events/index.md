---
layout: page
title: GTN Event Horizon
---

![Logo with a woman standing in front of a globe with multiple location pins in it]({% link assets/images/events.svg %}){: style="float:right; width:25%;"}

This page lists all past and upcoming GTN-related events organised by members of the community.

For a complete list of all Galaxy-related events, please also see the [Galaxy Event Horizon](https://galaxyproject.org/events/)

Anybody is welcome to add their events here.

<a href="TODO"><button type="button" class="btn btn-success">Add your own event!</button></a>

# Training Events

{% assign events = site.pages |  where: "layout", "event" | sort: 'date_start' | reverse %}


<table class="eventtable table table-striped">
 <thead>
  <tr>
   <th>Date</th>
   <th>Topic/Event</th>
   <th>Venue/Location</th>
   <th>Contact</th>
  </tr>
 </thead>
 <tbody>
 {% for event in events  %}
 <tr>
  <td class="eventtable-date"> {{event | collapse_date_pretty }} </td>
  <td>
   <a class="eventtable-title" href="{{site.baseurl}}{{event.url}}">{{event.title}}</a>
   <div class="eventtable-description"> {{event.description}} </div>
  </td>
  <td> {{event.location | format_location}} </td>
  <td> {% for org in event.contributions.organisers %}
			{% include _includes/contributor-badge-inline.html id=org %}
		{% endfor %}
  </td>
 </tr>
 {% endfor %}
 </tbody>
</table>


