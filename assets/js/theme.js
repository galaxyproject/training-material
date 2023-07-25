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

function getThemePreference(){
	// If there's a theme specified in the URL
	params = (new URL(document.location)).searchParams;
	paramTheme = params.get('theme')
	// Immediately override our preferences with it
	if(paramTheme){
		setThemePreference(paramTheme, false);
	}
	return gtnLocalGet('training-theme');
}


function setThemePreference(newValue, automatic) {
	gtnLocalSet('training-theme', newValue);

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
}

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
