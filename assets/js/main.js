// Handle foldable challenges and solutions (on click and at start).
$(".solution>h3,.details>h3").click(function(event) {
    $(">*:not(h3)", $(this).parent()).toggle(400);
    $(">span.fold-unfold", this).toggleClass("fa-plus-square fa-minus-square");
});

$(".solution,.details").each(function() {
    $(">*:not(h3)", this).toggle();
    var h3 = $("h3:first", this);
    h3.append("<span class='fold-unfold fa fa-plus-square'></span>");
});

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
            // Loop over all the `<details>` tags and open them.
            var details = document.querySelectorAll('details');
            Array.prototype.forEach.call(details, function (detail) {
                detail.querySelector('summary').textContent = 'Answers';
                detail.setAttribute('open', true);
            });
        }
    });

})(window, document);