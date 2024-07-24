---
title: Standards compliant™ training news via RSS/Atom
layout: news

tags:
- gtn infrastructure
contributions:
  authorship:
    - hexylena
  infrastructure:
    - hexylena
---

<!-- GTN:IGNORE:003 -->

The GTN loves leveraging existing standards for both old and new purposes. One of our favourite standards that we adhere to is providing an RSS/Atom news feed for sharing news and events and updates!

**Update 2024 Jun 18**: We've added a new home page for all of our feeds: [GTN feeds home]({% link feeds/index.md %})

## Why RSS/Atom?

RSS feeds provide access to the latest news and updates from websites. They are a great way to keep up to date with your favourite websites without having to visit them every day, and without having to sign up for a mailing list or social media account. You can subscribe to feeds using a [feed reader](https://aboutfeeds.com/), which will automatically check for updates and display them for you.

There are some other benefits of RSS feeds however:

### Privacy

Our feed subscriptions are completely private! We do not need a mailing list or any other personal information to send you updates. You can subscribe to our feed and receive updates without us knowing anything about you. You can subscribe or unsubscribe at any time.

### Standards

[Atom feeds](https://en.wikipedia.org/wiki/Atom_(web_standard)) (similar to the more well known RSS feeds) are a well-established standard ([IETF RFC 4287](https://datatracker.ietf.org/doc/html/rfc4287)) for 'news' type feeds. This means that we can easily integrate with other systems that understand Atom feeds! Any feed reader you choose will be able to understand our feed.

<a href="http://validator.w3.org/feed/check.cgi?url=http%3A//training.galaxyproject.org/training-material/feed.xml">
  <img src="{% link assets/images/valid-atom.png %}" alt="[Valid Atom 1.0]" title="Validate my Atom 1.0 feed" />
</a>

### Environmental Friendliness

Feeds are small, lightweight, and easily cached by servers. Our feeds are generated at least once a day and cached aggressively with support for features like ETags and Last-Modified headers. This means that you can subscribe to our feed and receive updates without causing any additional load on our servers. Additionally checking if our feed is updated does not require any computation on our part! We just serve a static file to you.

### Syndication

This is one of the *coolest* benefits of standard RSS feeds, we believe: how many other systems interoperate with them! For instance if you wanted to subscribe to GTN news in your Slack channel, you could use the existing [RSS integration](https://slack.com/help/articles/218688467-Add-RSS-feeds-to-Slack) to do so without any additional coding or customisation required. Just run

```
/feed subscribe http://training.galaxyproject.org/training-material/feed.xml
```

If you wanted to see the same news in your Matrix channel, you can use one of the [RSS bridges](https://gitlab.com/imbev/matrix-rss-bridge) to do so. You'll get notified of every new post!

We even use this feed to post our updates to the Galaxy Hub

## What's in a feed?

If you've never looked at our feed before, we recommend you do, [it's a bit fancier than the average RSS feed](https://training.galaxyproject.org/training-material/feed.xml)! Each news item provides a short preview to let you decide if you're interested in reading more. If you are, you can click through to the full post on our website.

## GTN Feeds

We actually produce a large number of RSS feeds which can let you subscribe to topics you're interested in:

- [News](https://training.galaxyproject.org/training-material/feed.xml)
- [Events](https://training.galaxyproject.org/training-material/events/feed.xml)

We also produce per-topic feeds if you are just interested in a single topic and want to be informed about new tutorials or slides added to that topic:

- [Galaxy Interface Updates](https://training.galaxyproject.org/training-material/topics/galaxy-interface/feed.xml)
- [Single Cell Updates](https://training.galaxyproject.org/training-material/topics/single-cell/feed.xml)
- [Fair Training Updates](https://training.galaxyproject.org/training-material/topics/fair/feed.xml)

These feeds are updated every time we add new content to the GTN, so you can be sure you're always up to date! You can find their links on every tutorial page.

## More Feeds in the Future

We're exploring the possibility for 'tag' based feeds, where you can subscribe to anything added to our site with a given tag. Say you were interested in genome annotation, with a tag based feed you'd get to hear about any new tutorials, slides, events, FAQs, or news posts added to the GTN with that tag.

On top of that we're considering producing a 'daily digest' feed and potentially a 'weekly digest' like we currently create manually for posting to our chat channels in matrix. This would give us a lot more flexibility in how to share our news and updates with you. Perhaps we would also create those digests for specific topics, so you could subscribe to a weekly digest of all new single cell tutorials, for example.

Would you use [these features](https://github.com/galaxyproject/training-material/discussions/4999)? Let us know in [the discussions](https://github.com/galaxyproject/training-material/discussions/5000)!
