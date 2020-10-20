---
layout: base
---

{% assign sorted_topics = "" | split: "," %}
{% assign sorted_topics_pre = site.data | sort | order: "title" %}

{% for topic in sorted_topics_pre %}
    {% if topic[0] == "introduction" %}
        {% assign sorted_topics = sorted_topics | unshift: topic %}
    {% else %}
        {% assign sorted_topics = sorted_topics | push: topic %}
    {% endif %}
{% endfor %}

{% include _includes/default-header.html %}


<style type="text/css">
::cue {
   color:#333;
}
</style>
<div class="container main-content">
	<section>
		<h1>Videos</h1>

		{% for topic in sorted_topics %}

			{% assign t = site.data[topic[0]] %}
			{% assign topic_material = site.pages | topic_filter:topic[0] %}

			{% for material in topic_material %}
				{% if material.video %}
					<h2 id="video-{{ material.topic_name }}-{{ material.tutorial_name }}">{{ topic[1].title }} / {{ material.title }}</h2>
					<video controls preload="metadata" width="800">
						<source
							src="https://galaxy-training.s3.amazonaws.com/videos/topics/{{ material.topic_name }}/tutorials/{{ material.tutorial_name }}/slides.mp4"
							type="video/mp4">
						<track label="English" kind="captions" srclang="en"
							src="https://training.galaxyproject.org/videos/topics/{{ material.topic_name }}/tutorials/{{ material.tutorial_name }}/slides.en.vtt" default>
					</video>

				{% endif %}
			{% endfor %}

		{% endfor %}


	</section>
</div>

{% include _includes/default-footer.html %}
