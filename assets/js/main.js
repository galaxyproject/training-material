// make boxes collapsible
//LEGACY
$(".solution>h3,.details>h3,.tip>h3,.question>h3,.hands_on>h3,.comment>h3").click(function(event) {
    $(">*:not(h3)", $(this).parent()).toggle(400);
    $(">span.fold-unfold", this).toggleClass("fa-plus-square fa-minus-square");
});

//NEW
$("blockquote.solution>.box-title>button,blockquote.details>.box-title>button,blockquote.tip>.box-title>button,blockquote.question>.box-title>button,blockquote.hands_on>.box-title>button,blockquote.comment>.box-title>button").click(function(event) {
	// If they click the icons, we need to make sure we have a reference to the button properly
	var button;
	if(event.target.nodeName === "BUTTON"){
		button = $(event.target)
	} else {
		button = $(event.target).parents("button");
	}
	// Then we can collapse specifically things based on the button's ancestry
	var parentBlockquote = button.parents("blockquote")[0];

	// Collapse every child of the blockquote, that is NOT a box-title
    $(">*:not(.box-title)", parentBlockquote).toggleClass("box-collapsed");
	// And toggle our icon
    $(">span.fold-unfold", button).toggleClass("fa-plus-square fa-minus-square");

	// Naturally we also need to toggle the aria-expanded attribute to make sure we're accessible
    $(this).attr("aria-expanded",
		// if it's collapsed (i.e. showing plus icon indicating expanding)
		$(">span.fold-unfold", this).hasClass("fa-plus-square") ?
		// mark as collapsed
		"false" : "true"
	);
});

// collapse some box types by default
// LEGACY
 $(".solution>h3,.details>h3,.tip>h3").each(function() {
    $(">*:not(h3)", $(this.parent)).toggle("box-collapsed");
    $(this).append("<span role='button' class='fold-unfold fa fa-plus-square'></span>");
});


// NEW
$("blockquote.solution,blockquote.details,blockquote.tip").each(function() {
    $(">.box-title>button", this).click();
});

$("section.tutorial .hands_on,section.tutorial .hands-on").each((idx, el) => {
	var box_id = $(".box-title", el).attr("id");
	$(el).append(`
		<p class="text-muted" style="text-align:right;font-size:0.9rem;">
			<a href="#${box_id}">Link to here</a> |
			<i class="far fa-question-circle" aria-hidden="true"></i> <a href="./faqs/">FAQs</a> |
			<a href="https://gitter.im/Galaxy-Training-Network/Lobby">Gitter Chat</a> |
			<a href="https://help.galaxyproject.org">Help Forum</a>
		</p>`
	);
})

// CYOA Support
function cyoaChoice(text){
	if(text !== undefined && text !== null){
		var loc = new URL(document.location)
		try {
			localStorage.setItem(`gtn-cyoa-${loc.pathname}`, text);
		} catch(e) {
			// Helaas pindakaas
		}

		var inputs = document.querySelectorAll(".gtn-cyoa input"),
			options = [...inputs].map(x => x.value),
			nonMatchingOptions = options.filter(x => x !== text);

		nonMatchingOptions.forEach(value => {
			document.querySelectorAll(`.${value}`).forEach(el => el.classList.add("gtn-cyoa-hidden"));
		})

		document.querySelectorAll(`.${text}`).forEach(el => el.classList.remove("gtn-cyoa-hidden"));

		// Just in case we mark it as checked (e.g. if default/from URL)
		var input_el = document.querySelector(`input[value="${text}"]`)
		// Can be undefined
		if(input_el) {
			input_el.checked = true;
		}
	}
}

function cyoaDefault(defaultOption){
	// Start with the URL parameter
	var loc = new URL(document.location)
	var urlOption = loc.searchParams.get("gtn-cyoa");
	if(urlOption){
		cyoaChoice(urlOption);
		return;
	}

	// Otherwise fall back to local storage (survives refreshes)
	var lsOption;
	try {
		lsOption = localStorage.getItem(`gtn-cyoa-${loc.pathname}`);
	} catch(e) {
		// Helaas pindakaas
	}
	if(lsOption !== null && lsOption !== undefined){
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
        var withAnswers = (new URL(document.location)).searchParams.get("with-answers");
        if (withAnswers !== null) {
            // Same as above selector
            $(".solution>.box-title button,.details>.box-title button").click();
        }

        var expandAll = (new URL(document.location)).searchParams.get("expand-all");
        if (expandAll !== null) {
            $(".solution>.box-title button,.details>.box-title button,.tip>.box-title button").click();

        }
        // collapse all boxes on the faq overview pages
        if (window.location.href.indexOf("faqs") > -1) {
            $(".hands_on>.box-title,.question>.box-title,.comment>.box-title").click();
        }

        var handsOnOnly = (new URL(document.location)).searchParams.get("only-hands-on");
		if(handsOnOnly !== null) {
			$(".tutorial .container .col-sm-10>:not(.hands_on)").hide()
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

// For admin training
document.querySelectorAll("article.topic-admin section#tutorial-content div.language-diff pre code").forEach(codeBlock => fixDiffPresentation(codeBlock))
document.querySelectorAll("article.topic-data-science section#tutorial-content div.language-diff pre code").forEach(codeBlock => fixDiffPresentation(codeBlock))

// Redirects
if(window.location.hostname === "galaxyproject.github.io") {
	// Redirect
	var redirect = "https://training.galaxyproject.org" + window.location.pathname + window.location.search;
	$('div.container.main-content').prepend("<div class='alert alert-warning'><strong>Note: </strong>This content has a new home at <a href=\"" + redirect + "\">" + redirect + "</a>, which you will be redirected to in 5 seconds.</div>");

	window.setTimeout(function(){
	window.location.href = redirect;
	}, 5000)
}

// Copy paste buttons
document.querySelectorAll('div.highlight').forEach((snippet) => {
	// Google translate has additional #text nodes mixed in with
	// the pre for some reason.
	var gtn_snippet_pres = [...snippet.childNodes].filter(x => x.tagName == "PRE")
	if(gtn_snippet_pres && gtn_snippet_pres.length > 0){
		gtn_snippet_pres[0].insertAdjacentHTML('beforebegin','<button class="btn btn-light" data-clipboard-snippet tabindex="0"><i class="fa fa-copy"></i>&nbsp;Copy</button>');
	}
});

var clipboardSnippets=new ClipboardJS('[data-clipboard-snippet]',{
    target:function(trigger){return trigger.nextElementSibling;
}});
