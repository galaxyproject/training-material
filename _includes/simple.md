{% for day in site.data.training_sessions %}
 {% assign daynum = day[0] %}

 {% if day[1].subday %}
 <h4 class="daystart" style="margin-top:1em;"> {{day[1].title}} </h4>
 {% else %}
 <h3 id="{{ day[0] }}" class="daystart" style="margin-top:1em;"> {{day[1].title}} </h3>
 {% endif %}
 <p> {{day[1].description}} </p>

  {% assign daysessions = day[1].sessions %}
  {% assign selfstudy = day[1].selfstudy %}
  {% if daysessions %}
   {% include _includes/program_table.html sessions=daysessions fullprogram=true selfstudy=selfstudy %}
  {% endif %}
{% endfor %} <!-- end schedule -->

