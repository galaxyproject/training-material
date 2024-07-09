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
- Monthly/Weekly/Daily Rollups, which include all *new* tutorials, slides, workflows, FAQs, events, and contributors.
    - [GTN Monthly]({{ site.baseurl }}/feeds/matrix-month.xml)
    - [GTN Weekly]({{ site.baseurl }}/feeds/matrix-week.xml)
    - [GTN Daily]({{ site.baseurl }}/feeds/matrix-day.xml)
- Per-tag Monthly Rollups
    - [Single Cell Month]({{ site.baseurl }}/feeds/single-cell-month.xml)
    - [One Health Month]({{ site.baseurl }}/feeds/one-health-month.xml)
    - ...and every other topic / tag based topic (i.e. topics linked from the home page)

These are available as an [OPML file as well]({{ site.baseurl }}/feeds/gtn.opml).

{% else %}
GTN Feed listing is not available in development mode. (This is done so we don't need to generate the feed pages or add an exception to our URL checking, while keeping CI times fast.)
{% endif %}

## Embedding Feeds

Any[^1] of the above feeds can be embedded anywhere you like. Simply replace
`.xml` with `.w.xml` in the URL and it'll produce a feed preview that is more

[^1]: minus the main news feed currently as that is produced by a third party plugin
amenable to embedding.

<iframe width="340px" height="600px" src="/training-material/events/feed.w.xml"></iframe>
<iframe width="340px" height="600px" src="/training-material/topics/admin/feed.w.xml"></iframe>
<iframe width="340px" height="600px" src="/training-material/topics/one-health/feed.w.xml"></iframe>

```html
<iframe width="340px" height="600px" src="/training-material/events/feed.w.xml"></iframe>
<iframe width="340px" height="600px" src="/training-material/topics/admin/feed.w.xml"></iframe>
<iframe width="340px" height="600px" src="/training-material/topics/one-health/feed.w.xml"></iframe>
```

### Rollups

We created these 'rollups' based on the bot that posted updates to our Matrix channel, they are simply digests of recent changes to help keep community members up to date.

<iframe width="340px" height="600px" src="/training-material/feeds/matrix-month.w.xml"></iframe>
<iframe width="340px" height="600px" src="/training-material/feeds/matrix-week.w.xml"></iframe>
<iframe width="340px" height="600px" src="/training-material/feeds/matrix-day.w.xml"></iframe>

These can be embedded like so:

```html
<iframe width="340px" height="600px" src="/training-material/feeds/matrix-month.w.xml"></iframe>
<iframe width="340px" height="600px" src="/training-material/feeds/matrix-week.w.xml"></iframe>
<iframe width="340px" height="600px" src="/training-material/feeds/matrix-day.w.xml"></iframe>
```

### Community Specific Rollups

Some rollups were created for individual communities:

<iframe width="340px" height="600px" src="/training-material/feeds/single-cell-month.w.xml"></iframe>
<iframe width="340px" height="600px" src="/training-material/feeds/one-health-month.w.xml"></iframe>

These are easily embedded, note the `.w.xml` ending, indicating a widget. (This is simply used to provide an alternate XSLT that renders better in an `iframe`.

```html
<iframe width="340px" height="600px" src="/training-material/feeds/single-cell-month.w.xml"></iframe>
<iframe width="340px" height="600px" src="/training-material/feeds/one-health-month.w.xml"></iframe>
```

### Differences between feeds

Having two feeds with the same data might seem a bit odd but we have two separate user stories we want to address:

- **Feed Reader**: This is the feed you'd subscribe to in your feed reader, and more importantly, the sort of feed that you'd send to someone else if they were curious how to follow updates to the GTN.
- **Feed Widget**: This is the feed you'd embed in your website, or in a widget on a page. It's designed to be more visually appealing and easier to read in a small space.
