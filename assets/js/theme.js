(function (window, document) {
	function onDocumentReady(fn) {
		if (document.attachEvent ? document.readyState === "complete" : document.readyState !== "loading") {
			fn();
		} else {
			document.addEventListener('DOMContentLoaded', fn);
		}
	}

	onDocumentReady(function () {
		function getCookie(){
			return document.cookie.replace(/(?:(?:^|.*;\s*)training-theme\s*\=\s*([^;]*).*$)|^.*$/, "$1");
		}

		function setCookie(newValue) {
			document.cookie = "training-theme=; expires=Thu, 01 Jan 1970 00:00:00 GMT";
			document.cookie = "training-theme=" + newValue + ";path=/";
		}


		function setTheme(theme){
			var old_classes = $("body").attr("class");
			if (old_classes === undefined){
				old_classes = "";
			}
			var new_classes = old_classes.split(' ').filter(function(x){ return x.substring(0, 6) !== "theme-" })

			if (theme !== "default"){
				new_classes.push('theme-' + theme)
			}

			$("body").attr('class', new_classes.join(' '))
			setCookie(theme);
		}

		$("#theme-selector").click(function(evt){
			var theme = $(evt.target).data('value');
			setTheme(theme);
		})

		training_theme_cookie = getCookie()
		if(training_theme_cookie){
			setTheme(training_theme_cookie);
		}
	});
})(window, document);

