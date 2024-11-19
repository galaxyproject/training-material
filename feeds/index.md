---
layout: page
title: GTN Feeds
---

We offer a wide range of [RSS feeds]({{ site.baseurl }}/news/2024/06/04/gtn-standards-rss.html) to help you keep up to date with the latest training materials and events. You can subscribe to these feeds in your favourite RSS reader, or embed them in your own website.

## Feed Directory

{% if jekyll.environment == "production" %}

- [GTN News]({{ site.baseurl }}/feed.xml)
- [GTN Events]({{ site.baseurl }}/events/feed.xml)
- Topic Feeds, which include all *new* tutorials, slides, FAQs, workflows, and events.
    - [Single Cell]({{ site.baseurl }}/topics/single-cell/feed.xml)
    - [Admin Training]({{ site.baseurl }}/topics/admin/feed.xml)
    - ...and every other topic
- Monthly/Weekly/Daily digests, which include all *new* tutorials, slides, workflows, FAQs, events, and contributors.
    - [GTN Firehose]({{ site.baseurl }}/feeds/matrix-all.xml) (every change as it happens!)
    - [GTN Monthly]({{ site.baseurl }}/feeds/matrix-month.xml)
    - [GTN Weekly]({{ site.baseurl }}/feeds/matrix-week.xml)
    - [GTN Daily]({{ site.baseurl }}/feeds/matrix-day.xml)
- Per-tag Monthly digests
    - [Single Cell Month]({{ site.baseurl }}/feeds/single-cell-month.xml)
    - [One Health Month]({{ site.baseurl }}/feeds/one-health-month.xml)
    - ...and every other topic / tag based topic (i.e. topics linked from the home page)

These are available as an [OPML file as well]({{ site.baseurl }}/feeds/gtn.opml).

{% else %}
GTN Feed listing is not available in development mode. (This is done so we don't need to generate the feed pages or add an exception to our URL checking, while keeping CI times fast.)
{% endif %}

## Embedding Feeds

Any[^1] of the above feeds can be embedded anywhere you like. Simply replace
`.xml` with `.w.html` in the URL and it'll produce a feed preview that can be embedded in an iframe easily. (`.w.xml` is also available but does not work for Safari users, ~10% of our traffic.)

[^1]: minus the main news feed currently as that is produced by a third party plugin
amenable to embedding.

<iframe width="340px" height="600px" src="/training-material/events/feed.w.html"></iframe>
<iframe width="340px" height="600px" src="/training-material/topics/admin/feed.w.html"></iframe>
<iframe width="340px" height="600px" src="/training-material/topics/one-health/feed.w.html"></iframe>

```html
<iframe width="340px" height="600px" src="{{ site.url }}{{ site.baseurl }}/events/feed.w.html"></iframe>
<iframe width="340px" height="600px" src="{{ site.url }}{{ site.baseurl }}/topics/admin/feed.w.html"></iframe>
<iframe width="340px" height="600px" src="{{ site.url }}{{ site.baseurl }}/topics/one-health/feed.w.html"></iframe>
```

### Digests

We created these 'digests' based on the bot that posted updates to our Matrix channel, they are simply digests of recent changes to help keep community members up to date.

<iframe width="340px" height="600px" src="/training-material/feeds/matrix-month.w.html"></iframe>
<iframe width="340px" height="600px" src="/training-material/feeds/matrix-week.w.html"></iframe>
<iframe width="340px" height="600px" src="/training-material/feeds/matrix-day.w.html"></iframe>

These can be embedded like so:

```html
<iframe width="340px" height="600px" src="{{ site.url }}{{ site.baseurl }}/feeds/matrix-month.w.html"></iframe>
<iframe width="340px" height="600px" src="{{ site.url }}{{ site.baseurl }}/feeds/matrix-week.w.html"></iframe>
<iframe width="340px" height="600px" src="{{ site.url }}{{ site.baseurl }}/feeds/matrix-day.w.html"></iframe>
```

### Community Specific Digests

Some digests were created for individual communities:

<iframe width="340px" height="600px" src="/training-material/feeds/single-cell-month.w.html"></iframe>
<iframe width="340px" height="600px" src="/training-material/feeds/one-health-month.w.html"></iframe>

These are easily embedded, note the `.w.html` ending, indicating a widget. (There is a `.w.xml` version with an atom feed and separate XSLT style sheet, but we've generated the HTML from that automatically as well, which is more compatible with Safari users..)

```html
<iframe width="340px" height="600px" src="{{ site.url }}{{ site.baseurl }}/feeds/single-cell-month.w.html"></iframe>
<iframe width="340px" height="600px" src="{{ site.url }}{{ site.baseurl }}/feeds/one-health-month.w.html"></iframe>
```

### Differences between feeds

Having two feeds with the same data might seem a bit odd but we have two separate user stories we want to address:

- **Feed Reader**: This is the feed you'd subscribe to in your feed reader, and more importantly, the sort of feed that you'd send to someone else if they were curious how to follow updates to the GTN.
- **Feed Widget (XML)**: This is the Atom feed version of the widget that you could subscribe to. You can embed this directly in an iframe, if you do not care about Safari users.
- **Feed Widget (HTML)**: This is a pre-XSLT'd version of the Atom feed, ready for user in an iframe for all users!
