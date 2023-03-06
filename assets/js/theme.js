function getThemePreference(){
	params = (new URL(document.location)).searchParams;
	paramTheme = params.get('theme')

	if(paramTheme){
		setThemePreference(paramTheme);
	}
	return localStorage.getItem('training-theme');
}

function setThemePreference(newValue) {
	localStorage.setItem('training-theme', newValue);
}


function setTheme(theme){
	var old_classes = document.getElementsByTagName("body")[0].getAttribute("class")
	if (old_classes === undefined || old_classes === null){
		old_classes = "";
	}
	var new_classes = old_classes.split(' ').filter(function(x){ return x.substring(0, 6) !== "theme-" })
	new_classes.push('theme-' + theme)

	document.getElementsByTagName("body")[0].setAttribute("class", new_classes.join(' '))
	setThemePreference(theme);
}

training_theme_cookie = getThemePreference()
if(training_theme_cookie){
	setTheme(training_theme_cookie);
}

const d = new Date();
let month = d.getMonth();

// Just for February
if(month === 1 && window.location.pathname === "/"){
	if(getThemePreference() === null || getThemePreference() === "undefined" || getThemePreference() === undefined){
		setTheme("blm");
	}
}

// Just for July
if(month === 5 && window.location.pathname === "/"){
	if(getThemePreference() === null || getThemePreference() === "undefined" || getThemePreference() === undefined){
		setTheme("progress");
	}
}
