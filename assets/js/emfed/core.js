import { getToots } from "./client.js";
import DOMPurify from "../dompurify/purify.es.js";
/**
 * Mark a string as safe for inclusion in HTML.
 */
function safe(s) {
    return Object.assign(new String(s), { __safe: null });
}
/**
 * Format a value as a string for templating.
 */
function flat(v) {
    if (typeof v === "undefined" || v === null) {
        return "";
    }
    else if (typeof v === "string" || v instanceof String) {
        if (v.hasOwnProperty("__safe")) {
            return v;
        }
        else {
            // Escape strings for inclusion in HTML.
            return v
                .replaceAll("&", "&amp;")
                .replaceAll("<", "&lt;")
                .replaceAll(">", "&gt;")
                .replaceAll('"', "&quot;")
                .replaceAll("'", "&#039;");
        }
    }
    else {
        return v.map(flat).join("");
    }
}
/**
 * The world's dumbest templating system.
 */
function html(strings, ...subs) {
    let out = strings[0];
    for (let i = 1; i < strings.length; ++i) {
        out += flat(subs[i - 1]);
        out += strings[i];
    }
    return safe(out);
}
/**
 * Render a single toot object as an HTML string.
 */
function renderToot(toot) {
    // Is this a boost (reblog)?
    let boost = null;
    if (toot.reblog) {
        boost = {
            avatar: toot.account.avatar,
            username: toot.account.username,
            display_name: toot.account.display_name,
            user_url: toot.account.url,
        };
        toot = toot.reblog; // Show the "inner" toot instead.
    }
    const date = new Date(toot.created_at).toLocaleString();
    const images = toot.media_attachments.filter((att) => att.type === "image");
    return html `<li class="toot">
    <a class="permalink" href="${toot.url}">
      <time datetime="${toot.created_at}">${date}</time>
    </a>
    ${boost &&
        html ` <a class="user boost" href="${boost.user_url}">
      <img class="avatar" width="23" height="23" src="${boost.avatar}" />
      <span class="display-name">${boost.display_name}</span>
      <span class="username">@${boost.username}</span>
    </a>`}
    <a class="user" href="${toot.account.url}">
      <img class="avatar" width="46" height="46" src="${toot.account.avatar}" />
      <span class="display-name">${toot.account.display_name}</span>
      <span class="username">@${toot.account.username}</span>
    </a>
    <div class="body">${safe(DOMPurify.sanitize(toot.content))}</div>
    ${images.map((att) => html ` <a
          class="attachment"
          href="${att.url}"
          target="_blank"
          rel="noopener noreferrer"
        >
          <img
            class="attachment"
            src="${att.preview_url}"
            alt="${att.description}"
          />
        </a>`)}
  </li>`.toString();
}
/**
 * Get the toots for an HTML element and replace that element with the
 * rendered toot list.
 */
export async function loadToots(element) {
    // Fetch toots based on the element's `data-toot-*` attributes.
    const el = element;
    const toots = await getToots(el.href, el.dataset.tootAccountId, Number(el.dataset.tootLimit ?? 5), el.dataset.excludeReplies === "true");
    // Construct the HTML content.
    const list = document.createElement("ol");
    list.classList.add("toots");
    el.replaceWith(list);
    for (const toot of toots) {
        const html = renderToot(toot);
        list.insertAdjacentHTML("beforeend", html);
    }
}
/**
 * Transform all links on the page marked with the `mastodon-feed` class.
 */
export function loadAll() {
    document.querySelectorAll("a.mastodon-feed").forEach(loadToots);
}
