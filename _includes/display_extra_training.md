{% for training in include.extra_trainings %}
    {% if training.type == "internal" %}
        {% assign extra_topic_metadata = site.data[training.topic_name] %}
        {% assign extra_topic = site | topic_filter:training.topic_name %}
        {% unless training.tutorials %}
        <li>
          <a href="{{ site.baseurl }}/topics/{{ training.topic_name }}">{{ extra_topic_metadata.title }}</a>
        </li>
        {% else training.tutorials %}
            {% for extra_tuto in training.tutorials %}
                {% for topic_tuto in extra_topic %}
                    {% if extra_tuto == topic_tuto.tutorial_name %}
                        {% if topic_tuto.slides %}
                            <li>
                              <a href="{{ site.baseurl }}/topics/{{ training.topic_name }}/tutorials/{{ topic_tuto.tutorial_name }}/slides.html">{% icon slides %} Slides: {{ topic_tuto.title }}</a>
                            </li>
                        {% endif %}
                        {% if topic_tuto.hands_on %}
                            {% if topic_tuto.hands_on_url %}
                                <li>
                                  <a href="{{ topic_tuto.hands_on_url }}">{% icon tutorial %} Hands-on: {{ topic_tuto.title }}</a>
                                </li>
                            {% else %}
                                <li>
                                  <a href="{{ site.baseurl }}/topics/{{ training.topic_name }}/tutorials/{{ topic_tuto.tutorial_name }}/tutorial.html">{% icon tutorial %} Hands-on: {{ topic_tuto.title }}</a>
                                </li>
                            {% endif %}
                        {% endif %}
                    {% endif %}
                {% endfor %}
            {% endfor %}
        {% endunless %}
    {% elsif training.type == "none" %}
      <li>
        {{ training.title }}
      </li>
    {% else %}
      <li>
        <a href="{{ training.link }}">{{ training.title }}</a>
      </li>
    {% endif %}
{% endfor %}
