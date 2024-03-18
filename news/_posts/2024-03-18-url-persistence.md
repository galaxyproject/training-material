---
title: "Cool URLs Don't Change, GTN URLs don't either."
layout: news
tags:
- gtn infrastructure
- new feature
contributions:
  authorship:
  - hexylena
  infrastructure:
  - hexylena
---

At the Galaxy Training Network we are *really* committed to ensuring our training materials are [Findable, Accessible, Interoperable, and Reusable]({% link faqs/gtn/fair_training.md %}).
This means that we want to make sure that the URLs to our training materials are persistent and don't change.
The GTN wants you to be able to rely on our URLs once you've added them to a poster or training material, without having to worry about them breaking in the future.

For a long time we relied on contributors ensuring that when files are merged, we add appropriate redirects to each moved file, however this isn't a very reliable system.
We'd recently also introduced [Persistent URLs (PURLs) to our lessons as well]({{ site.baseurl }}/news/2023/04/19/shortlinks.html) but that only helps our users going forward, it doesn't ensure we are meeting our earlier commitments.

So now we've added a new test to each GTN merge that checks URLs from the last 3 years to ensure that they are *all* still working.
If they aren't, the merge will fail and the contributor will need to fix the URLs before the changes can be accepted.

We'll be expanding how far back we check URLs in the future, but for now, this will help us ensure that our URLs are completely persistent!

By implementing this we discovered only 50 pages (out of ~5.3k GTN pages) that had been moved without proper redirections, and we've fixed them all!
