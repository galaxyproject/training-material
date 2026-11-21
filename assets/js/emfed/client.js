/**
 * Fetch recent toots for a user, given their Mastodon URL.
 */
export async function getToots(userURL, accountId, limit, excludeReplies) {
    const url = new URL(userURL);
    // Either use the account id specified or look it up based on the username
    // in the link.
    const userId = accountId ??
        (await (async () => {
            // Extract username from URL.
            const parts = /@(\w+)$/.exec(url.pathname);
            if (!parts) {
                throw "not a Mastodon user URL";
            }
            const username = parts[1];
            // Look up user ID from username.
            const lookupURL = Object.assign(new URL(url), {
                pathname: "/api/v1/accounts/lookup",
                search: `?acct=${username}`,
            });
            return (await (await fetch(lookupURL)).json())["id"];
        })());
    // Fetch toots.
    const tootURL = Object.assign(new URL(url), {
        pathname: `/api/v1/accounts/${userId}/statuses`,
        search: `?limit=${limit ?? 5}&exclude_replies=${!!excludeReplies}`,
    });
    return await (await fetch(tootURL)).json();
}
