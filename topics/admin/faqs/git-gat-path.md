---
title: Galaxy Admin Training Path
area: ansible
box_type: comment
layout: faq
contributors: [hexylena]
---

The yearly Galaxy Admin Training follows a specific ordering of tutorials. Use this timeline to help keep track of where you are in Galaxy Admin Training.

{% assign tutorials = "ansible-galaxy tus cvmfs singularity tool-management data-library connect-to-compute-cluster job-destinations pulsar gxadmin monitoring tiaas reports ftp" | split: " " %}

{% assign seen_tuto = 0 %}
<ol id="git-gat-timeline">
{% for tutorial in tutorials %}
    <a href="{{ site.baseurl }}/topics/admin/tutorials/{{ tutorial }}/tutorial.html">
    <li class="{% if include.tutorial == tutorial %}active{% elsif seen_tuto == 0 %}disabled{% endif %}">
        <span>Step {{ forloop.index }}</span>
        <span>{{ tutorial }}</span>
    </li>
    {% if include.tutorial == tutorial %}{% assign seen_tuto = 1 %}{% endif %}
    </a>
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
    background: #86d486;
    color: black;
}
#git-gat-timeline li.disabled {
    background: #eee;
}
#git-gat-timeline span {
    align-self: center;
}
</style>
