{% capture newLine %}
{% endcapture %}
{% for training in include.extra_trainings %}
    {% if training.type == "internal" %}
        {% assign extra_topic_metadata = site.data[training.topic_name] %}
        {% assign extra_topic = site | topic_filter:training.topic_name %}
            {% capture topic_desc %}[{{ extra_topic_metadata.title }}]({{ site.baseurl }}/topics/{{ training.topic_name }}){% endcapture %}
        {% if training.tutorials %}
            {% for extra_tuto in training.tutorials %}
                {% for topic_tuto in extra_topic %}
                    {% if extra_tuto == topic_tuto.tutorial_name %}
                        {% assign tuto_desc = topic_tuto.title | append: ": " | prepend: "    - " | prepend: newLine %}
                        {% if topic_tuto.slides %}
                            {% capture tuto_slide_desc %}[{% icon slides %} slides]({{ site.baseurl }}/topics/{{ training.topic_name }}/tutorials/{{ topic_tuto.tutorial_name }}/slides.html){% endcapture %}
                            {% assign tuto_desc = tuto_desc | append: tuto_slide_desc %}
                            {% assign sep = " - " %}
                        {% endif %}
                        {% if topic_tuto.hands_on %}
                            {% if topic_tuto.hands_on_url %}
                                {% capture tuto_hands_on_desc %}{{ sep }}[{% icon tutorial %} hands-on]({{ topic_tuto.hands_on_url }}){% endcapture %}
                            {% else %}
                                {% capture tuto_hands_on_desc %}{{ sep }}[{% icon tutorial %} hands-on]({{ site.baseurl }}/topics/{{ training.topic_name }}/tutorials/{{ topic_tuto.tutorial_name }}/tutorial.html){% endcapture %}
                            {% endif %}
                            {% assign tuto_desc = tuto_desc | append: tuto_hands_on_desc %}
                        {% endif %}
                    {% assign topic_desc = topic_desc | append: tuto_desc %}
                    {% endif %}
                {% endfor %}
            {% endfor %}
        {% endif %}
- {{ topic_desc }}
    {% elsif training.type == "none" %}
- {{ training.title }}
    {% else %}
- [{{ training.title }}]({{ training.link }})
    {% endif %}
{% endfor %}