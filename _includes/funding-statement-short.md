{% if include.contributions.funding %}
{% for id in include.contributions.funding %}
	{% assign name = site.data.contributors[id].name | default: id -%}
	<div>
		<div>
		{% if site.data.contributors[id].avatar %}
		<img class="funder-avatar" src="{{ site.data.contributors[id].avatar }}" alt="Logo">
		{% else %}
		<img class="funder-avatar" src="https://avatars.githubusercontent.com/{{ id }}" alt="Logo">
		{% endif %}
		</div>
		<div>
		{% if site.data.contributors[id].funding_id %}
		{{ site.data.contributors[id].funding_id }}
		{% endif %}
		</div>
	</div>
{% endfor %}
{% endif %}
