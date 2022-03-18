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
$("section.tutorial.topic-admin div.language-diff pre code .gi,section.tutorial.topic-admin div.language-diff pre code .gd").each((x, e) => {
  var fixed = $(e).text().split('\n').map(q => { return q.slice(1) }).join('\n');
  $(e).text(fixed);
})
$("section.tutorial.topic-data-science div.language-diff pre code .gi,section.tutorial.topic-data-science div.language-diff pre code .gd").each((x, e) => {
  var fixed = $(e).text().split('\n').map(q => { return q.slice(1) }).join('\n');
  $(e).text(fixed);
})

