// make boxes collapsible
$("blockquote.solution>.box-title>button,blockquote.details>.box-title>button,blockquote.tip>.box-title>button,blockquote.question>.box-title>button,blockquote.hands_on>.box-title>button,blockquote.comment>.box-title>button").click(function(event) {
	// If they click the icons, we need to make sure we have a reference to the button properly
	var button;
	if(event.target.nodeName === "BUTTON"){
		button = $(event.target)
	} else {
		button = $(event.target).parents("button");
	}
	// Then we can collapse specifically things based on the button's ancestry
	var parentBlockquote = button.parents("blockquote");

	// Collapse every child of the blockquote, that is NOT a box-title
    $(">*:not(.box-title)", parentBlockquote).toggle();
	// And toggle our icon
    $(">span.fold-unfold", button).toggleClass("fa-plus-square fa-minus-square");

	// Naturally we also need to toggle the aria-expanded attribute to make sure we're accessible
    $(this).attr("aria-expanded",
		// if it's collapsed
		$(">span.fold-unfold", this).hasClass("fa-plus-square") ? 
		// mark as collapsed
		"true" : "false"
	);
});

// NEW
$("blockquote.solution,blockquote.details,blockquote.tip").each(function() {
    $(">.box-title", this).click();
});

$("section.tutorial .hands_on,section.tutorial .hands-on").append('<p class="text-muted" style="text-align:right;font-size:0.9rem;"><i class="far fa-question-circle" aria-hidden="true"></i> <a href="./faqs/">FAQs</a> | <a href="https://gitter.im/Galaxy-Training-Network/Lobby">Gitter Chat</a> | <a href="https://help.galaxyproject.org">Help Forum</a></p>');

// CYOA Support
function cyoaChoice(text){
	if(text !== undefined && text !== null){
		localStorage.setItem('gtn-cyoa', text);

		var inputs = document.querySelectorAll(".gtn-cyoa input"),
			options = [...inputs].map(x => x.value),
			nonMatchingOptions = options.filter(x => x !== text);

		nonMatchingOptions.forEach(value => {
			document.querySelectorAll(`.${value}`).forEach(el => el.classList.add("gtn-cyoa-hidden"));
		})

		document.querySelectorAll(`.${text}`).forEach(el => el.classList.remove("gtn-cyoa-hidden"));

		// Just in case we mark it as checked (e.g. if default/from URL)
		document.querySelector(`input[value="${text}"]`).checked = true
	}
}

function cyoaDefault(defaultOption){
	// Start with the URL parameter
	var urlOption = (new URL(document.location)).searchParams.get("gtn-cyoa");
	if(urlOption){
		cyoaChoice(urlOption);
		return;
	}

	// Otherwise fall back to local storage (survives refreshes)
	var lsOption = localStorage.getItem('gtn-cyoa');
	if(lsOption !== null){
		cyoaChoice(lsOption);
		return;
	}

	// Otherwise if the browser is remembering for us, use that.
	var currentlySelected = [...document.querySelectorAll("input[name='cyoa']")].filter(x => x.checked)[0];
	if(currentlySelected){
		cyoaChoice(currentlySelected);
		return;
	}

	// And failing that, use the default.
	cyoaChoice(defaultOption);
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
        // If you pass `?with-answers` to an URL, it will automatically open
        // the `<details>` blocks (i.e. the Q&A sections most of the time).
        if (window.location.search.match(/\?with-answers/gi)) {
            // Same as above selector
            $(".solution>h3,.details>h3,.tip>h3").click();
            $(".solution>.box-title,.details>.box-title,.tip>.box-title").click();

        }
        if (window.location.search.match(/\?with-answers/gi)) {
            // Same as above selector
            $(".solution>h3,.details>h3,.tip>h3,.hands_on>h3,.question>h3").click();
            $(".solution>.box-title,.details>.box-title,.tip>.box-title,.hands_on>.box-title,.question>.box-title").click();

        }
        // collapse all boxes on the faq overview pages
        if (window.location.href.indexOf("faqs") > -1) {
            $(".hands_on>h3,.question>h3,.comment>h3").click();
            $(".hands_on>.box-title,.question>.box-title,.comment>.box-title").click();
        }
    });

})(window, document);


function fixDiffPresentation(codeBlock){
	codeBlock.childNodes.forEach(x => {
		if(x.nodeName == '#text'){
			x.textContent = x.textContent.split('\n').map(q => { return q.startsWith(" ") ? q.slice(1) : q }).join('\n')
		} else {
			if(!(x.nodeName.toLowerCase() === 'span' && x.classList[0] === 'notranslate')){
				var fixed = $(x).text().split('\n').map(q => { return q.slice(1) }).join('\n');
				$(x).text(fixed);
			}
		}
	})
}

<!--  For admin training -->
document.querySelectorAll("section.tutorial.topic-admin div.language-diff pre code").forEach(codeBlock => fixDiffPresentation(codeBlock))
document.querySelectorAll("section.tutorial.topic-data-science div.language-diff pre code").forEach(codeBlock => fixDiffPresentation(codeBlock))
