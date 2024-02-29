---
layout: page
title: Theme
redirect_from:
- /preferences
---

## Colour Scheme

We have attempted to decompose our themes into several axes depending on your accessibility needs:

This control responds to [prefers-color-scheme](https://developer.mozilla.org/en-US/docs/Web/CSS/@media/prefers-color-scheme)

<select class="form-control theme-control" id="brightness" onchange="savePrefs()">
	<option value="auto">Auto</option>
	<option value="light">Light</option>
	<option value="dark">Dark</option>
</select>


### Preferred Light Theme

<select class="form-control theme-control" id="light_theme" onchange="savePrefs()">
	<option value="white">White</option>
	<option value="yellow">Paper</option>
</select>

### Preferred Dark Theme

<select class="form-control theme-control" id="dark_theme" onchange="savePrefs()">
	<option value="night">Night</option>
	<option value="midnight">Midnight</option>
</select>

### UI Contrast

This control responds to [prefers-contrast](https://developer.mozilla.org/en-US/docs/Web/CSS/@media/prefers-contrast)

<select class="form-control theme-control" id="contrast" onchange="savePrefs()">
	<option value="auto">Auto</option>
	<option value="low">Low</option>
	<option value="high">High</option>
</select>


### UI Theme

<select class="form-control theme-control" id="theme" onchange="savePrefs()">
	<option value="default">Default</option>
	<option value="rainbow">🌈</option>
	<option value="blm">✊🏿</option>
	<option value="halloween">🎃</option>
	<option value="progress">🏳️‍🌈</option>
	<option value="trans">🏳️‍⚧️ </option>
	<option value="straya">🇦🇺</option>
</select>


<script>
function savePrefs() {
	// Convert this into a hash
	var prefs = {};
	[...document.querySelectorAll(".theme-control")]
		.map(x => { return [x.id, x.value]})
		.forEach(x => { prefs[x[0]] = x[1] })
	gtnLocalSet('theme2', JSON.stringify(prefs))
	processTheme2();


	if(prefs.theme === "straya"){
		document.body.classList.add('downunder');
		setTimeout(function(){
			document.body.classList.remove('downunder');
		}, 8000);
	}
}

function restorePrefs(){
	var prefs = JSON.parse(gtnLocalGet("theme2")) || {};
	Object.keys(prefs).forEach(k => {
		document.getElementById(k).value = prefs[k]
	})
	processTheme2();
}
restorePrefs();
</script>

## Fonts

While some of our fonts are available on Google Fonts, we **do not** use Google Fonts to serve them.

<select class="form-control font-control" id="font" onchange="saveFont()">
	<option value="default">Default (Atkinson Hyperlegible)</option>
	<option value="open-dyslexic">Open Dyslexic</option>
	<option value="comic-sans">Comic Sans</option>
</select>

Font for code blocks:

<select class="form-control font-control-code" id="font-code" onchange="saveFont()">
	<option value="default">Default</option>
	<option value="comic-sans">Comic Sans Mono</option>
</select>

<script>
function saveFont(){
	gtnLocalSet("fontMain", document.getElementById("font").value);
	gtnLocalSet("fontCode", document.getElementById("font-code").value);

	document.body.dataset["font_main"] = document.getElementById("font").value
	document.body.dataset["font_code"] = document.getElementById("font-code").value
}
document.getElementById("font").value = gtnLocalGet("fontMain");
document.getElementById("font-code").value = gtnLocalGet("fontCode");
</script>

Are you looking for Data Privacy settings? [They have moved to their own page]({% link user/privacy.md %})
