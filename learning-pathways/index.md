---
layout: page
---

{% assign pathways = site.pages | where: "layout", "learning-pathway" | sort: "priority", "last" %}


# Learning Pathways

Learning pathways are sets of tutorials curated for you by community experts to form a coherent set of lessons around a topic, building up knowledge as you go. We always recommend to follow the tutorials in the order they are listed in the pathway.

<!-- list all available pathways as cards  -->
<div class="pathway">
{% for path in pathways %}

{% assign coverimage = path.coverimage | default: "/assets/images/GTN.png" %}
{% assign coverimagealt = path.coverimagealt | default: "GTN logo with a multi-coloured star and the words Galaxy Training Network"%}
<div class="card pathwayitem">
 <div class="card-header">
   <a href="{{site.baseurl}}{{path.url}}"><h3 class="card-title">{{path.title}}</h3></a>
 </div>
 <div class="row no-gutters">
  <div class="col-sm-5">
   {% if path.cover-image %}
   <img class="card-img pathwaycover" src="{{site.baseurl}}/{{path.cover-image}}" alt="{{ path.cover-image-alt }}" loading="lazy">
   {% else %}
   <img class="card-img pathwaycover" src="{{site.baseurl}}{{coverimage}}" alt="{{ coverimagealt }}" loading="lazy">
   {% endif %}
  </div>
  <div class="col-sm-7">
   <div class="card-body">
   <p></p>
   {% for tag in path.tags %}
   <span class="label label-default tutorial_tag" style="{{ tag | colour_tag }}">{{ tag  }}</span>
    {% endfor %}
   <hr/>
   <p class="card-text">
     {{ path.description }}
   </p>
   <a href="{{site.baseurl}}{{path.url}}" class="btn btn-primary">View Learning Pathway</a>
   </div>
  </div>
 </div>
</div>

{% endfor %}
</div>
