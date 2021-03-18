---
layout: page
---

{% assign contributors = site.data['contributors'] %}

# The Latest GTN News

Keep an eye on this page for the latest news around the GTN. New tutorials, GTN features, upcoming training events, and much much more!


<div class="newslist">
{% for n in site.categories['news'] %}{% unless n.external %}
{% if n.cover %}
 {% assign coverimage = n.cover %}
{% elsif n.tags contains "cofest" %}
 {% assign coverimage = "/assets/images/cofest.png" %}
{% elsif n.tags contains "galaxy" %}
 {% assign coverimage = "/assets/images/GalaxyNews.png" %}
{% else %}
 {% assign coverimage = "/assets/images/GTN.png" %}
{% endif %}


<div class="card newsitem">
 <div class="card-header">
     <a href="{{site.baseurl}}{{n.url}}"><h3 class="card-title">{{n.title}}</h3></a>
  </div>
 <div class="row no-gutters">
  <div class="col-sm-5">
   <img class="card-img newscover" src="{% unless coverimage contains 'http' %}{{site.baseurl}}/{% endunless %}{{coverimage}}" alt="cover image for news item">
  </div>
  <div class="col-sm-7">
   <div class="card-body">
        <!--<a href="{{site.baseurl}}{{n.url}}"><h4 class="card-title">{{n.title}}</h4></a>-->
        {% if n.contributors %}
        <div class="contributors-line"> {% include _includes/contributor-badge-list.html contributors=n.contributors %}</div>
        {% endif %}
        {% for tag in n.tags %}
<button class="label label-default tutorial_tag" id="{{ tag }}" style="{{ tag | colour_tag }}" title="Click to show all tutorials with this tag">{{ tag  }}</button>
{% endfor %}
    <hr/>
    <p class="card-text">{{n.excerpt}}</p>
    <a href="{{site.baseurl}}{{n.url}}" class="btn btn-primary">Full Story</a>
   </div>
  </div>
 </div>
</div>
{% endunless %}{% endfor %}
</div>

