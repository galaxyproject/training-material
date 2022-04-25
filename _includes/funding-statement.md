{% if include.contributions.funding %}
<div markdown="1">
<h2 id="funding">{{locale['references']| default: "Funding" }}</h2>
<p>These individuals or organisations provided funding support for the development of this resource</p>

<div class="row">
{% for id in include.contributions.funding %}
	{% assign name = site.data.contributors[id].name | default: id -%}
	<div class="col-md-3 col-xs-12">
		{% if site.data.contributors[id].avatar %}
		<img class="funder-avatar" src="{{ site.data.contributors[id].avatar }}" alt="Logo">
		{% else %}
		<img class="funder-avatar" src="https://avatars.githubusercontent.com/{{ id }}" alt="Logo">
		{% endif %}
		<a href="{{ site.baseurl }}/hall-of-fame/{{ id }}/" class="btn btn-secondary">See Funder Profile</a>
	</div>
	<div class="col-md-9 col-xs-12">
		{{ site.data.contributors[id].funding_statement | markdownify }}
	</div>
{% endfor %}
</div>
</div>
{% endif %}
