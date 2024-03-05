---
layout: page
title: Theme Customiser
---

This definitely doesn't need to exist but it's fun to play with.

<div class="row">
<div class="col-md-3">
<table id="theme">
</table>
</div>

<div class="col-md-9">
<p>Here's a preview of what your changes will look like.</p>

<iframe id="preview" src="{{ site.baseurl }}/topics/admin/tutorials/ansible-galaxy/tutorial.html" width="100%" height="800px"></iframe>

</div>

</div>

<script>
const css_var = [
  '--color-background',
  '--text-color',
  '--text-color-muted',
  '--text-color-inverted',
  '--text-color-boxtitle',
  '--slide-heading-color',
  '--code-foreground',
  '--code-background',
  '--border-light',
  '--hyperlink',
  '--hyperlink-hover',
// From Galaxy
  '--brand-color',
  '--brand-color-contrast',
  '--navbar-border-size',
  '--navbar-border-color',
  '--navbar-item-background',
  '--navbar-item-hover',
  '--secondary-color',
  '--secondary-color-text',
  '--outline-color',
  '--border',
  '--color-disabled',
];

const table = document.getElementById('theme');
css_var.forEach((v) => {
  const row = table.insertRow(-1);
  const cell1 = row.insertCell(0);
  const cell2 = row.insertCell(1);
  cell1.innerHTML = v.replace('--', '').replace(/-/g, ' ').replace(/\b\w/g, l => l.toUpperCase());
  const input = document.createElement('input');
  input.type = 'color';
  input.value = getComputedStyle(document.documentElement).getPropertyValue(v);
  input.oninput = () => {
    document.documentElement.style.setProperty(v, input.value);
    document.getElementById("preview").contentDocument.documentElement.style.setProperty(v, input.value);
  };
  cell2.appendChild(input);
});

</script>
