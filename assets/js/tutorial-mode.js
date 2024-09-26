/* GTN Tutorial Mode
 * Setup messaging to share events with any Galaxy that is embedding this content.
 *
 * Originally we suggested admins to setup a proxy of the GTN at
 * /training-material/ ensuring that Galaxy and the GTN were on the same
 * domain, the same origin.
 *
 * This let us intermingle the site's content a bit, and call javascript from both sides safely.
 *
 * We replace it now with postMessage, a widely available and supported API for sharing messages and data between domains.
 * https://caniuse.com/mdn-api_window_postmessage
 *
 * **Messages**
 *
 * The following communication occurs
 *
 * GTN                | Direction | Galaxy
 * ------------       | --------- | ------
 * Navigation         | ->        | Galaxy stores this in localStorage
 * Scroll             | ->        | Galaxy stores this in localStorage
 * Click to Load Tool | ->        | Galaxy loads the tool at the specified version
 * Click to Load WF   | ->        | Galaxy activates the TRS API.
 * Scrolls            | <-        | User opens the GTN in Tutorial Mode and Galaxy restores the scroll position
 *
 * These are the *only* supported messages.
 *
 * All GTN messages are sent as `action@@data@@other` where `@@` is the
 * delimiter, chosen since unlikely to occur naturally
 *
 * The only Galaxy message sent to us is scroll, and it is just a number. Any
 * non-numeric characters are stripped before a scroll is attempted.
 *
 * **Security**
 * precautions are taken where possible and useful. On the Galaxy
 * side we'll need to check the origin of the request and ensure it is really
 * training.galaxyproject.org (if we ever change domains again, like we've done
 * in the past, this should be acceptable as long as the old domain continues
 * to load and does not redirect automatically. If it does redirect, tutorial
 * mode will stop working (messages will not be processed) until the
 * corresponding Galaxy update is made.)
 */

// Wrapper for postMessage to ensure we always send the same format
function gtnPostMessage(category, message,other){
	parent.postMessage(`${category}@@${message}@@${other}`, "*")
}

// On load we inform any parent of which url we're at.
// So when the user e.g. navigates pages this can be stored on the galaxy side.
gtnPostMessage('navigate', window.location)

// We'll also track the scroll position to allow restoring that..
//
// https://developer.mozilla.org/en-US/docs/Web/API/Document/scroll_event
//
// Code samples added on or after August 20, 2010 are in the public domain CC0.
// No licensing notice is necessary but if you need one, you can use:
//
// Any copyright is dedicated to the Public Domain:
// https://creativecommons.org/publicdomain/zero/1.0/
let lastKnownScrollPosition = 0;
let ticking = false;
document.addEventListener("scroll", (event) => {
	lastKnownScrollPosition = window.scrollY;
	if (!ticking) {
		window.requestAnimationFrame(() => {
			gtnPostMessage('scroll', lastKnownScrollPosition)
			ticking = false;
		});
		ticking = true;
	}
});

// If the user clicks on a tool button/link, we'll inform the parent (Galaxy)
// of that so it can load said tool.
document.querySelectorAll("button.tool,span.tool").forEach((button) => {
	button.addEventListener("click", (event) => {
		gtnPostMessage(
			'loadTool',
			event.target.getAttribute("data-tool"),
			event.target.getAttribute("data-version")
		)
	});
});

// Same for workflows, no trs-version/trs-server attribute? TODO.
document.querySelectorAll("button.workflow,span.workflow").forEach((button) => {
	button.addEventListener("click", (event) => {
		gtnPostMessage(
			'loadWorkflow',
			event.target.getAttribute("data-workflow")
		)
	});
});


// Lastly, we listen to messages from the parent which enable us to restore our Scroll offset
//
// Every article you read about postMessage is concerned about security, and the origin of the message.
// We here clearly ignore it, so we are playing with potentially dangerous code.
//
// If anyone is able to inject a message into the page, and we are careless
// with this data, they could compromise the GTN.
//
// So we need to be extremely careful with the data we receive.
//
window.addEventListener('message', e => {
	// Cross domain messages are allowed from anywhere, since we don't know
	// who will embed the GTN and we would like to be able to support this
	// scroll behaviour for improved user experience.

	// It *MUST* be an integer.
	// So we serialise whatever we get to a string (we don't care about any
	// data inside, just, get it as quickly as possible to a string. And
	// make sure it contains only numbers.
	var safeScrollText = `${e.data}`.replace(/[^0-9]/g, '');
	var scrollPosition = parseInt(safeScrollText, 10);
	if (scrollPosition) {
		window.scrollTo(0, scrollPosition);
	}
	return;
},false);
