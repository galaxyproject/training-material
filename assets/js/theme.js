function gtnLocalSet(key, value){
	// Work around for localStorage not being available in Safari private browsing
	// https://stackoverflow.com/a/14555361/358804
	try {
		localStorage.setItem(key, value);
	} catch (e) {
		console.log(`localStorage not available, cannot set ${key} to ${value}`);
	}
}

function gtnLocalGet(key){
	// Work around for localStorage not being available in Safari private browsing
	try {
		return localStorage.getItem(key);
	} catch (e) {
		console.log(`localStorage not available, cannot get ${key}`);
		return undefined;
	}
}

function getTheme2(){
	return JSON.parse(gtnLocalGet('theme2')) || {brightness: "auto", light_theme: "white", dark_theme: "night", contrast: "auto", theme: "default"};
}

function getThemePreference(){
	// If there's a theme specified in the URL
	params = (new URL(document.location)).searchParams;
	paramTheme = params.get('theme')
	// Immediately override our preferences with it
	if(paramTheme){
		setThemePreference(paramTheme, false);
	}
	return getTheme2().theme;
}


function setThemePreference(newValue, automatic) {
	var tmpprefs = getTheme2();
	tmpprefs['theme'] = newValue;
	gtnLocalSet('theme2', JSON.stringify(tmpprefs));

	// If this setTheme call was triggered "automatically", not by user choice, mark it as such
	if(automatic === true){
		gtnLocalSet('training-theme-automatic', 'true');
	} else {
		gtnLocalSet('training-theme-automatic', 'false');
	}
}

function setTheme(theme, automatic){
	let body_elem = document.getElementsByTagName("body")[0];
	var old_classes = body_elem.getAttribute("class")
	if (old_classes === undefined || old_classes === null){
		old_classes = "";
	}
	var new_classes = old_classes.split(' ').filter(function(x){ return x.substring(0, 6) !== "theme-" })
	new_classes.push('theme-' + theme)

	body_elem.setAttribute("class", new_classes.join(' '))

	setThemePreference(theme, automatic);
}

// If it's in the URL, or in a saved peference
training_theme_cookie = getThemePreference()
if(training_theme_cookie){
	// Then restore the theme to what's in the URL/preference.
	setTheme(training_theme_cookie, false);
	var prefs = getTheme2();
	prefs['theme'] = training_theme_cookie;
	gtnLocalSet('theme2', JSON.stringify(prefs));
}

// Theme2
function processTheme2(){
	// MediaQueryList
	const darkModePreference = window.matchMedia("(prefers-color-scheme: dark)");
	if (darkModePreference.matches) {
		document.body.dataset.brightness = 'dark';
	} else {
		document.body.dataset.brightness = 'light';
	}

	// recommended method for newer browsers: specify event-type as first argument
	darkModePreference.addEventListener("change", e => {
		console.log(e.matches);
		if (e.matches) {
			document.body.dataset.brightness = 'dark';
		} else {
			document.body.dataset.brightness = 'light';
		}
	});


	const contrastModePreferenceMore = window.matchMedia("(prefers-contrast: more)");
	const contrastModePreferenceLess= window.matchMedia("(prefers-contrast: less)");
	if (contrastModePreferenceMore.matches) {
		document.body.dataset.contrast = 'high';
	}
	if (contrastModePreferenceLess.matches) {
		document.body.dataset.contrast = 'low';
	}

	// recommended method for newer browsers: specify event-type as first argument
	contrastModePreferenceMore.addEventListener("change", e => {
		if (e.matches) {
			document.body.dataset.contrast = 'high';
		}
	});
	contrastModePreferenceLess.addEventListener("change", e => {
		if (e.matches) {
			document.body.dataset.contrast = 'low';
		}
	});

	// Do this last so it overrides the above
	if (gtnLocalGet('theme2') !== null){
		var prefs = getTheme2();
		console.log("prefs", prefs);
		Object.keys(prefs).forEach(function(key) {
			if (key === "brightness" && prefs[key] === "auto"){
				return;
			}
			document.body.dataset[key] = prefs[key];
		})
	}

	if(gtnLocalGet('fontMain') !== null){
		document.body.dataset["font_main"] = gtnLocalGet("fontMain");
	}
	if(gtnLocalGet('fontCode') !== null){
		document.body.dataset["font_code"] = gtnLocalGet("fontCode");
	}

}
processTheme2();



// However we have some additional things:
var gtnThemeCurrentDate = new Date();
var gtnThemeCurrentMonth = gtnThemeCurrentDate.getMonth();
// If we had automatically set the theme in the past, or never automatically set one before
var gtnThemeIsAutomatic = (gtnLocalGet('training-theme-automatic') === "true" || gtnLocalGet('training-theme-automatic') === undefined)
// If the user had manually chosen a theme, then gtnThemeIsAutomatic === false

if(gtnThemeCurrentMonth === 1 && window.location.pathname === "/training-material/"){ // Just for February
	// If the user hasn't manually chosen a theme
	if(gtnThemeIsAutomatic || getThemePreference() === null || getThemePreference() === "undefined" || getThemePreference() === undefined){
		setTheme("blm", true);
	}
}
else if(gtnThemeCurrentMonth === 5 && window.location.pathname === "/training-material/"){ // Just for June
	// If the user hasn't manually chosen a theme
	if(gtnThemeIsAutomatic || getThemePreference() === null || getThemePreference() === "undefined" || getThemePreference() === undefined){
		setTheme("progress", true);
	}
}
else { // Not one of the "special" months
	if (gtnThemeIsAutomatic){
		// Then we mark it as automatic
		gtnLocalSet('training-theme-automatic', 'true');
		// And set the theme back to default
		setTheme("default", true);
	}
}

// URL overrides
// If there's a theme specified in the URL
var theme2Override = getTheme2();
var params = (new URL(document.location)).searchParams;
['brightness', 'light_theme', 'dark_theme', 'contrast', 'theme'].forEach(function(key) {
	if(params.get(key)){
		console.log(`Processing URL override for ${key}=${params.get(key)}`);
		document.body.dataset[key] = params.get(key);
		theme2Override[key] = params.get(key);
	}
})
if(Object.keys(theme2Override).length > 0){
	gtnLocalSet('theme2', JSON.stringify(theme2Override));
}
