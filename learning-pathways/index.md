---
layout: page
---

{% assign pathways_science = site.pages | where: "layout", "learning-pathway" | where: "type", "use" | sort: "priority", "last" %}

{% assign pathways_other = site.pages | where: "layout", "learning-pathway" | where_exp: "item", "item.type != 'use'" | sort: "priority", "last" %}


# Learning Pathways

Learning pathways are sets of tutorials curated for you by community experts to form a coherent set of lessons around a topic, building up knowledge as you go. We always recommend to follow the tutorials in the order they are listed in the pathway.

## For Scientists


<!-- list all available pathways as cards  -->
<div class="pathwaylist row">


{% for path in pathways_science %}

{% include _includes/pathway-card.html path=path %}

{% endfor %}

</div>


## For Teachers, Developers and System Administrators

<!-- list all available pathways as cards  -->
<div class="pathwaylist row">

{% for path in pathways_other %}

{% include _includes/pathway-card.html path=path %}

{% endfor %}

</div>
