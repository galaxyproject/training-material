---
layout: page
title: Data Privacy
---

The GTN is invested in preserving your privacy. We attempt to balance your
rights and needs to private access of the GTN resources, with our needs to
continue receiving funding and provide a good experience for you, the user.

## Alternative Access

If you do not like the privacy options available to you, the GTN can be
[downloaded and hosted yourself, offline]({% link topics/contributing/tutorials/running-jekyll/tutorial.md %}),
for users needing completely anonymous access to the site and it's tutorials
and slides.

## Data Collection Systems

We use two data collection systems, for different purposes. Both, to the best
of our knowledge, render you effectively anonymous as they are currently
configured. Each system collects different data, and can be opted out of individually.

### Plausible

We collect the following information on visitors:

- Device information
- Location (country granularity)
- Page visited

We do not have information correlating visits across pages, we cannot identify
if you visit one page, or many pages in a single visit.

You can see the aggregate collected information in
[Plausible](https://plausible.galaxyproject.eu/training.galaxyproject.org/). 

**We do not have access to more granular data.**

However, if you wish you may opt out of these aggregate statistics.

> > <code-in-title>What you provide</code-in-title>
> > - Device information
> > - Location (country granularity)
> > - Page visited
> {: .code-in}
> 
> > <code-out-title>How we use it</code-out-title>
> > - Helps us request Grant Funding
> > - Helps us know who is using our tutorials
> > - Makes authors happy to know people love their tutorial!
> {: .code-out}
{: .code-2col}

Opt out from plausible:
<select class="form-control privacy-control" id="plausible-opt-out" onchange="savePrivacy()">
	<option value="opt-in">Opt-in</option>
	<option value="opt-out">Opt-out</option>
</select>

### Sentry

When something breaks on the GTN in the UI, we collect information about how that happened to help our developers fix problems within the GTN.

> > <code-in-title>What you provide</code-in-title>
> > - Device information
> > - Page
> > - Stack trace
> > - GTN Javascript variables
> {: .code-in}
>
> > <code-out-title>How we use it</code-out-title>
> > - Fixing JavaScript Bugs
> > - Especially from unusual platforms which we cannot test on.
> {: .code-out}
{: .code-2col}

Opt out from sentry:
<select class="form-control privacy-control" id="sentry-opt-out" onchange="savePrivacy()">
	<option value="opt-in">Opt-in</option>
	<option value="opt-out">Opt-out</option>
</select>


<script>
function savePrivacy() {
	gtnLocalSet('sentry-opt-out', document.getElementById("sentry-opt-out").value)
	gtnLocalSet('plausible-opt-out', document.getElementById("plausible-opt-out").value)

	if(document.getElementById("plausible-opt-out").value === "opt-in") {
		localStorage.removeItem("plausible_ignore")
	} else {
		localStorage.setItem("plausible_ignore", "true")
	}
}
// restore from prefs
document.getElementById("sentry-opt-out").value = gtnLocalGet("sentry-opt-out") || "opt-in";
document.getElementById("plausible-opt-out").value = gtnLocalGet("plausible-opt-out") || "opt-in";

if(navigator.doNotTrack === "1"){
	document.getElementById("sentry-opt-out").disabled = true
	document.getElementById("plausible-opt-out").disabled = true

	document.getElementById("sentry-opt-out").innerHTML = `<option value="opt-out">Opted-out (Do not track is set in your browser)</option>`
	document.getElementById("plausible-opt-out").innerHTML = `<option value="opt-out">Opted-out (Do not track is set in your browser)</option>`
}
</script>

## Settings

Here we expose the several types of data storage that browsers offer, so you can easily see what data this website creates and uses:

### LocalStorage

This data is **never** transferred to the server

<dl id="settings-data">
</dl>

<script>
let gtnSettingsKeys = Object.keys(window.localStorage);
gtnSettingsKeys.sort()
gtnSettingsKeys.forEach(k => {
	// Add a row to the table with this key/value
	var dt = document.createElement("dt");
	var dd = document.createElement("dd");
	dt.innerHTML = `${k}`;
	dd.innerHTML = `<code>${window.localStorage[k]}</code>`;
	document.getElementById("settings-data").appendChild(dt);
	document.getElementById("settings-data").appendChild(dd);
})
if(gtnSettingsKeys.length === 0){
	document.getElementById("settings-data").innerHTML = `There is no data.`;
}
</script>

### SessionStorage

This data is **never** transferred to the server, and is deleted when your browser window closes.

<dl id="session-data">
</dl>

<script>
let gtnSessionKeys = Object.keys(window.sessionStorage);
gtnSessionKeys.sort()
gtnSessionKeys.forEach(k => {
	// Add a row to the table with this key/value
	var dt = document.createElement("dt");
	var dd = document.createElement("dd");
	dt.innerHTML = `${k}`;
	dd.innerHTML = `<code>${window.sessionStorage[k]}</code>`;
	document.getElementById("session-data").appendChild(dt);
	document.getElementById("session-data").appendChild(dd);
})
if(gtnSessionKeys.length === 0){
	document.getElementById("session-data").innerHTML = `There is no data.`;
}
</script>

### Cookies

This data is automatically transferred to the server. You can clear it using your browser's tools:

<pre id="cookies-data">
</pre>

<script>
if(document.cookies !== undefined){
	document.getElementById('cookies-data').innerHTML = document.cookies;
} else {
	document.getElementById('cookies-data').innerHTML = "No cookies have been set.";
}
</script>
