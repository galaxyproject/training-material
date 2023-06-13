---
title: Galaxy Admin Training Path
area: ansible
box_type: comment
layout: faq
contributors: [hexylena]
---

The yearly Galaxy Admin Training follows a specific ordering of tutorials. Use this timeline to help keep track of where you are in Galaxy Admin Training.

{% assign tutorials = "admin/ansible-galaxy admin/backup-cleanup admin/customization admin/tus admin/cvmfs admin/apptainer admin/tool-management admin/reference-genomes admin/data-library dev/bioblend-api admin/connect-to-compute-cluster admin/job-destinations admin/pulsar admin/celery admin/gxadmin admin/reports admin/monitoring admin/tiaas admin/sentry admin/ftp admin/beacon" | split: " " %}


{% assign seen_tuto = 0 %}
<ol id="git-gat-timeline">
{% for tutorial_id in tutorials %}
    {% assign parts = tutorial_id | split: "/" %}
    {% assign tp_topic = parts[0] %}
    {% assign tp_tutorial = parts[1] %}

    <li class="{% if include.tutorial == tp_tutorial %}active{% elsif seen_tuto == 0 %}disabled{% endif %}">
        <a href="{{ site.baseurl }}/topics/{{ tp_topic }}/tutorials/{{ tp_tutorial }}/tutorial.html">
            <div>Step {{ forloop.index }}</div>
            <div>{% if tp_topic != "admin" %}{{ tp_topic }}/{% endif %}{{ tp_tutorial }}</div>
        </a>
    </li>
    {% if include.tutorial == tp_tutorial %}{% assign seen_tuto = 1 %}{% endif %}
    {% unless forloop.last %}
    <span aria-hidden="true">
        <i class="fas fa-arrow-right" aria-hidden="true"></i>
    </span>
    {% endunless %}
{% endfor %}
</ol>
