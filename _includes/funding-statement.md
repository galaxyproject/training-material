{% if include.contributions.funding %}
<div markdown="1">
<h2 id="funding">{{locale['references']| default: "Funding" }}</h2>
<p>These individuals or organisations provided funding support for the development of this resource</p>

<div class="d-flex flex-wrap">
{% for id in include.contributions.funding %}
	{% assign name = site.data.contributors[id].name | default: id -%}
	<a href="{{ site.baseurl }}/hall-of-fame/{{ id }}/" class="funder-badge">
		{% assign pfo = site.data.funders[id] | default: site.data.organisations[id] | default: site.data.contributors[id] | default: nil %}
		<div class="avatar">
			{% if pfo.avatar %}
			<img class="funder-avatar" src="{{ pfo.avatar }}" alt="Logo">
			{% else %}
				{% unless pfo.github %}
				<img class="funder-avatar" src="https://avatars.githubusercontent.com/{{ id }}" alt="Logo">
				{% endunless %}
			{% endif %}
		</div>

		<div class="info">
			<div class="name">{{ pfo.short_name | default: pfo.name | default: id }}</div>
			<div class="description">
			{{ site.data.funders[id].funding_statement | markdownify | strip_html }}
			</div>
		</div>
	</a>
{% endfor %}
</div>
</div>
{% endif %}
