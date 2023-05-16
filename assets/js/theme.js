function getThemePreference(){
	// If there's a theme specified in the URL
	params = (new URL(document.location)).searchParams;
	paramTheme = params.get('theme')
	// Immediately override our preferences with it
	if(paramTheme){
		setThemePreference(paramTheme, false);
	}
	return localStorage.getItem('training-theme');
}

function setThemePreference(newValue, automatic) {
	localStorage.setItem('training-theme', newValue);

	// If this setTheme call was triggered "automatically", not by user choice, mark it as such
	if(automatic === true){
		localStorage.setItem('training-theme-automatic', 'true');
	} else {
		localStorage.setItem('training-theme-automatic', 'false');
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
const currentDate = new Date();
const currentMonth = currentDate.getMonth();

if(currentMonth === 1 && window.location.pathname === "/"){ // Just for February
	// If the user hasn't chosen a theme
	if(getThemePreference() === null || getThemePreference() === "undefined" || getThemePreference() === undefined){
		setTheme("blm", true);
	}
}
else if(currentMonth === 6 && window.location.pathname === "/"){ // Just for June
	// If the user hasn't chosen a theme
	if(getThemePreference() === null || getThemePreference() === "undefined" || getThemePreference() === undefined){
		setTheme("progress", true);
	}
}
else { // Not one of the "special" months
	// If we had automatically set the theme in the past, or never automatically set one before
	if localStorage.getItem('training-theme-automatic') === "true" || localStorage.getItem('training-theme-automatic') === undefined {
		// Then we mark it as non-automatic
		localStorage.setItem('training-theme-automatic', 'false');
		// And set the theme back to default
		setTheme("default", false);
	}
}
