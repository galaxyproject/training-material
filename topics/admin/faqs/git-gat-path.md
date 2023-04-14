---
title: Galaxy Admin Training Path
area: ansible
box_type: comment
layout: faq
contributors: [hexylena]
---

The yearly Galaxy Admin Training follows a specific ordering of tutorials. Use this timeline to help keep track of where you are in Galaxy Admin Training.

{% assign tutorials = "ansible-galaxy backup-cleanup tus cvmfs apptainer tool-management data-library connect-to-compute-cluster job-destinations pulsar celery gxadmin monitoring tiaas reports ftp beacon" | split: " " %}

{% assign seen_tuto = 0 %}
<ol id="git-gat-timeline">
{% for tutorial in tutorials %}
    <li class="{% if include.tutorial == tutorial %}active{% elsif seen_tuto == 0 %}disabled{% endif %}">
        <a href="{{ site.baseurl }}/topics/admin/tutorials/{{ tutorial }}/tutorial.html">
            <div>Step {{ forloop.index }}</div>
            <div>{{ tutorial }}</div>
        </a>
    </li>
    {% if include.tutorial == tutorial %}{% assign seen_tuto = 1 %}{% endif %}
    {% unless forloop.last %}
    <span aria-hidden="true">
        <i class="fas fa-arrow-right" aria-hidden="true"></i>
    </span>
    {% endunless %}
{% endfor %}
</ol>

<style type="text/css">
#git-gat-timeline {
    display: flex;
    flex-direction: row;
    flex-wrap: wrap;
}
#git-gat-timeline li  {
    display: flex;
    flex-direction: column;
    border: 1px solid black;
    border-radius: 5px;
    padding: 0.5em;
    margin: 0.5em;
}
#git-gat-timeline li.active {
    background: #a8ffa8;
    color: black;
}
#git-gat-timeline li.disabled {
    background: #eee;
}
#git-gat-timeline span {
    align-self: center;
}
</style>
