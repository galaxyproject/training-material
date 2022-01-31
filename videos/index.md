---
title: GTN Videos
layout: page
redirect_from:
  - /topics/admin/videos/index
  - /topics/assembly/videos/index
  - /topics/climate/videos/index
  - /topics/computational-chemistry/videos/index
  - /topics/contributing/videos/index
  - /topics/dev/videos/index
  - /topics/ecology/videos/index
  - /topics/epigenetics/videos/index
  - /topics/galaxy-interface/videos/index
  - /topics/genome-annotation/videos/index
  - /topics/imaging/videos/index
  - /topics/instructors/videos/index
  - /topics/introduction/videos/index
  - /topics/metabolomics/videos/index
  - /topics/metagenomics/videos/index
  - /topics/proteomics/videos/index
  - /topics/sequence-analysis/videos/index
  - /topics/statistics/videos/index
  - /topics/transcriptomics/videos/index
  - /topics/variant-analysis/videos/index
  - /topics/visualisation/videos/index
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

The GTN now generates videos for selected slide decks. Click on a topic below to jump to the video page for that topic!

{% for topic in sorted_topics %}
{% assign topic_id = topic[0] %}
{% assign t = site.data[topic_id] %}

	{% assign has_video = false %}
	{% assign topic_material = site | topic_filter:topic[0] %}
	{% for material in topic_material %}
		{% if material.video %}
			{% assign has_video = true %}
		{% endif %}
	{% endfor %}

{% if has_video == true %}
<h2>{{ t.title }}</h2>
<div id="playlist">
	{% for material in topic_material %}
		{% if material.video %}
			{% capture vid %}{{ topic_id }}/{% if material.type == "introduction" %}slides/introduction{% else %}tutorials/{{ material.tutorial_name }}/slides{% endif %}{% endcapture %}
			<div class="pl-item">
				<a href="watch.html?v={{ vid }}">
					<div class="cover">
						<img src="https://training.galaxyproject.org/videos/topics/{{ vid }}.mp4.png" width="200px"/>
					</div>
					<div>
						<div class="title">{{ material.title }}</div>
					</div>
				</a>
			</div>
			{% for lang in material.translations.slides %}
			{% capture vid %}{{ topic_id }}/{% if material.type == "introduction" %}slides/introduction_{{ lang | upcase }}{% else %}tutorials/{{ material.tutorial_name }}/slides_{{ lang | upcase }}{% endif %}{% endcapture %}
			<div class="pl-item">
				<a href="watch.html?v={{ vid }}">
					<div class="cover">
						<img src="https://training.galaxyproject.org/videos/topics/{{ vid }}.mp4.png" width="200px"/>
					</div>
					<div>
						<div class="title">[{{ lang | upcase }}] {{ material.title }}</div>
					</div>
				</a>
			</div>
			{% endfor %}
		{% endif %}
	{% endfor %}
</div>
{% endif %}
{% endfor %}
