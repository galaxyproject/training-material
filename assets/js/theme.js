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
	var old_classes = $("body").attr("class");
	if (old_classes === undefined){
		old_classes = "";
	}
	var new_classes = old_classes.split(' ').filter(function(x){ return x.substring(0, 6) !== "theme-" })
	new_classes.push('theme-' + theme)

	$("body").attr('class', new_classes.join(' '))
	setThemePreference(theme);
}

(function (window, document) {
	function onDocumentReady(fn) {
		if (document.attachEvent ? document.readyState === "complete" : document.readyState !== "loading") {
			fn();
		} else {
			document.addEventListener('DOMContentLoaded', fn);
		}
	}

	onDocumentReady(function () {
		$("#theme-selector button").click(function(evt){
			var theme = $(evt.target).data('theme');
			setTheme(theme);
			if(theme === "straya"){
				$("body").addClass('downunder');
				setTimeout(function(){
					$("body").removeClass('downunder');
				}, 8000);
			}
		})

		training_theme_cookie = getThemePreference()
		if(training_theme_cookie){
			setTheme(training_theme_cookie);
		}

                // Just for february
                if(getThemePreference() === null || getThemePreference() === "undefined" || getThemePreference() === undefined){
                        setTheme("blm");
                }

	});
})(window, document);
