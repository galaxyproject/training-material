{% if include.contributions.funding %}
<div markdown="1">
<h2 id="funding">{{locale['references']| default: "Funding" }}</h2>
<p>These individuals or organisations provided funding support for the development of this resource</p>

<div class="row">
{% for id in include.contributions.funding %}
	{% assign name = site.data.contributors[id].name | default: id -%}
	<div class="col-md-3 col-xs-12">
		<a href="{{ site.baseurl }}/hall-of-fame/{{ id }}/" class="funder-badge">
			<div>
				{% if site.data.contributors[id].avatar %}
				<img class="funder-avatar" src="{{ site.data.contributors[id].avatar }}" alt="Logo">
				{% else %}
				<img class="funder-avatar" src="https://avatars.githubusercontent.com/{{ id }}" alt="Logo">
				{% endif %}
			</div>

			<div class="info">
				<div class="name">{{ site.data.funders[id].short_name | default: site.data.funders[id].name | default: id }}</div>
				<div class="description">
				{{ site.data.funders[id].funding_statement | markdownify | strip_html }}
				</div>
			</div>
		</a>
	</div>
	<div class="col-md-9 col-xs-12">
		{{ site.data.contributors[id].funding_statement | markdownify }}
	</div>
{% endfor %}
</div>
</div>
{% endif %}
