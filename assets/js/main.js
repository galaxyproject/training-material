// make boxes collapsible
$(".solution>h3,.details>h3,.tip>h3,.question>h3,.hands_on>h3,.comment>h3").click(function(event) {
    $(">*:not(h3)", $(this).parent()).toggle(400);
    $(">span.fold-unfold", this).toggleClass("fa-plus-square fa-minus-square");
});

// collapse some box types by default
$(".solution,.details,.tip").each(function() {
    $(">*:not(h3)", this).toggle();
    var h3 = $("h3:first", this);
    h3.append("<span role='button' class='fold-unfold fa fa-plus-square'></span>");
});

// don't collapse the others by default
$(".question,.hands_on,.comment").each(function() {
    var h3 = $("h3:first", this);
    h3.append("<span role='button' class='fold-unfold fa fa-minus-square'></span>");
});

$("section.tutorial .hands_on").append('<p class="text-muted" style="text-align:right;font-size:0.9rem;"><i class="far fa-question-circle" aria-hidden="true"></i> <a href="./faqs/">FAQs</a> | <a href="https://gitter.im/Galaxy-Training-Network/Lobby">Gitter Chat</a> | <a href="https://help.galaxyproject.org">Help Forum</a></p>');

// CYOA Support
function cyoaChoice(text){
	if(text !== undefined && text !== null){
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
	var currentlySelected = [...document.querySelectorAll("input[name='cyoa']")].filter(x => x.checked)[0];
	var urlOption = (new URL(document.location)).searchParams.get("gtn-cyoa");

	// If there's a URL parameter, use that over other options
	if(urlOption !== null){
		cyoaChoice(urlOption)
		return
	}

	// Chrome forgets
	if(currentlySelected === undefined) {
		cyoaChoice(defaultOption)
	} else {
		// Firefox doesn't.
		cyoaChoice(currentlySelected.value)
	}
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

        }
        if (window.location.search.match(/\?with-answers/gi)) {
            // Same as above selector
            $(".solution>h3,.details>h3,.tip>h3,.hands_on>h3,.question>h3").click();

        }
        // collapse all boxes on the faq overview pages
        if (window.location.href.indexOf("faqs") > -1) {
            $(".hands_on>h3,.question>h3,.comment>h3").click();
        }
    });

})(window, document);


<!--  For admin training -->
document.querySelectorAll("section.tutorial.topic-admin div.language-diff pre code").forEach(codeBlock => {
	codeBlock.childNodes.forEach(x => {
		if(x.nodeName == '#text'){
			x.textContent = x.textContent.split('\n').map(q => { return q.slice(1) }).join('\n')
		} else {
			var fixed = $(x).text().split('\n').map(q => { return q.slice(1) }).join('\n');
			$(x).text(fixed);
		}
	})
})

document.querySelectorAll("section.tutorial.topic-data-science div.language-diff pre code").forEach(codeBlock => {
	codeBlock.childNodes.forEach(x => {
		if(x.nodeName == '#text'){
			x.textContent = x.textContent.split('\n').map(q => { return q.slice(1) }).join('\n')
		} else {
			var fixed = $(x).text().split('\n').map(q => { return q.slice(1) }).join('\n');
			$(x).text(fixed);
		}
	})
})
