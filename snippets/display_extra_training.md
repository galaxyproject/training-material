{% for training in include.extra_trainings %}
    <li>
    {% if training.type == "internal" %}
        {% assign extra_topic_metadata = site.data[training.topic_name] %}
        {% assign extra_topic = site.pages | topic_filter:training.topic_name %}
        <a href="{{ site.baseurl }}/topics/{{ training.topic_name }}">{{ extra_topic_metadata.title }}</a>
        {% if training.tutorials %}
            <ul>
                {% for extra_tuto in training.tutorials %}
                    {% for topic_tuto in extra_topic %}
                        {% if extra_tuto == topic_tuto.tutorial_name %}
                            {% assign sep = "" %}
                            <li> {{ topic_tuto.title }}:
                            {% if topic_tuto.slides %}
                                <a href="{{ site.baseurl }}/topics/{{ training.topic_name }}/tutorials/{{ topic_tuto.tutorial_name }}/slides.html">{% icon slides %} slides</a>
                                {% assign sep = "-" %}
                            {% endif %}
                            {% if topic_tuto.hands_on %}
                                {{ sep }} <a href="{{ site.baseurl }}/topics/{{ training.topic_name }}/tutorials/{{ topic_tuto.tutorial_name }}/tutorial.html">{% icon tutorial %} hands-on</a>
                            {% endif %}
                            </li>
                        {% endif %}
                    {% endfor %}
                {% endfor %}
            </ul>
        {% endif %}
    {% else %}
        <a href="{{ training.link }}">{{ training.title }}</a>
    {% endif %}
    </li>
{% endfor %}