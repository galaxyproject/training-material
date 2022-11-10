---
layout: page
---

{% assign contributors = site.data['contributors'] %}

# The Latest GTN News

Keep an eye on this page for the latest news around the GTN. New tutorials, GTN features, upcoming training events, and much much more!

Want to add your own news here (e.g. new tutorial, event, publication, anything else training related)? Check out how to do that in [this FAQ]({% link faqs/gtn/gtn_news_create_post.md %})

<div class="newslist">
{% for n in site.categories['news'] %}


{% include _includes/news-card.html news=n %}

{% endfor %}
</div>

