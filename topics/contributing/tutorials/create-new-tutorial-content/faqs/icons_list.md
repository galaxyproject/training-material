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

Some icons have multiple aliases, any may be used, but we'd suggest trying to choose the most semantically appropriate one in case Galaxy later decides to change the icon.

The following icons are currently available:

<div class="row">
{% assign icon_groups = site['icon-tag'] | group_icons %}
{% for icon in icon_groups %}
	<div class="col-md-2 col-sm-3" style="text-align: center">
		<div style="font-size: 400%">{% icon_var icon[0][0] %}</div>
		<div>{% for z in icon[0] %}{{ z }}{%unless forloop.last%},{%endunless%} {% endfor %}</div>
	</div>
{% endfor %}
</div>

New icons can be added in `_config.yaml`, and you can search for the corresponding icons at https://fontawesome.com/v4/icons/
