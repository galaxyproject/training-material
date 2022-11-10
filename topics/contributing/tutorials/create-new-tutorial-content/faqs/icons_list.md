---
title: Which icons are available to use in my tutorial?
box_type: tip
layout: faq
contributors: [shiltemann]
---


To use icons in your tutorial, take the name of the icon, 'details' in this example, and write something like this in your tutorial:

```markdown
{% raw %}{% icon details %}{% endraw %}
```

The following icons are currently available:

<div class="row">
{% for icon in site["icon-tag"] %}
	<div class="col-md-2 col-sm-3" style="text-align: center">
		<div style="font-size: 400%">{% icon_var icon[0] %}</div>
		<div>{{ icon[0] }}</div>
	</div>
{% endfor %}
</div>
