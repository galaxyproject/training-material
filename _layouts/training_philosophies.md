---
layout: base
---

{% assign stories = site.posts %}

{% include _includes/default-header.html %}

<div class="container main-content">
    <section>
        <h1>Training Stories</h1>

        <p class="lead">
            We all teach in different ways, share your training philosophy and help other instructors find the way that is best for them!
        </p>

    {{ posts }}
    {% for p in posts %}
    {{ p }}
    {% endfor %}
    </section>
</div>

{% include _includes/default-footer.html %}
