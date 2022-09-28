---
layout: tutorial_hands_on
redirect_from:
- /topics/instructors/philosophies/
title: Teaching experiences
subtopic: practises
time_estimation: 1h
questions:
  - How to teach using Galaxy?
  - What are the different practises?
objectives:
  - Learn how others teach using Galaxy
key_points:
  - Every instructor has their prefered way to teach using Galaxy
contributors:
  - bebatut
  - fpsom
---

We all teach in different ways, share your training philosophy and help other instructors find the way that is best for them!

Some instructors have contributed their ideas. You can read their stories below.

# Teaching stories

{% assign stories = site.pages | where:"layout", "training_philosophy" %}
{% assign contributors = site.data['contributors'] %}

{% for p in stories %}
    {% assign username = p.username %}
<div id="{{ username }}">
    <h2>{{ contributors[username].name | default: username }}</h2>
    <div class="hall-of-fame-hero">
        <a href="{{ site.baseurl }}/hall-of-fame/{{ username }}/" class="thumbnail">
            <img src="https://avatars.githubusercontent.com/{{ username }}" alt="{{ username }}'s avatar">
        </a>
        <div class="contact-items">
            {% if contributors[username].github != false %}
                <a title="GitHub" href="https://github.com/{{ username }}">
                    {% icon github %}
                </a>
            {% endif %}
            {% if contributors[username].email %}
                <a title="E-mail" href="mailto:{{ contributors[username].email }}">
                    {% icon email %}
                </a>
            {% endif %}
            {% if contributors[username].gitter %}
                <a title="Gitter" href="https://gitter.im/{{ contributors[username].gitter }}">
                    {% icon gitter %}
                </a>
            {% endif %}
            {% if contributors[username].twitter %}
                <a title="Twitter" href="https://twitter.com/{{ contributors[username].twitter }}">
                    {% icon twitter %}
                </a>
            {% endif %}
            {% if contributors[username].linkedin %}
                <a title="LinkedIn" href="https://www.linkedin.com/in/{{ contributors[username].linkedin }}">
                    {% icon linkedin %}
                </a>
            {% endif %}
            {% if contributors[username].orcid %}
                <a title="ORCID" href="https://orcid.org/{{ contributors[username].orcid }}">
                    {% icon orcid %}
                </a>
            {% endif %}
        </div>
    </div>
    {{ p.content | markdownify }}
</div>
{% endfor %}