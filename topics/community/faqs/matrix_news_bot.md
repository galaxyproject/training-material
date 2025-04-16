---
title: How do I add a news feed to a Matrix channel?
box_type: tip
layout: faq
contributors: [nomadscientist, hexylena]
---

1. You must be an *Admin* in the channel. Find this out by going to the channel and selecting **Room info** --> **People**, or clicking on the little circle images of people in a channel. Admins can make other admins.

2. Go to **Room info** --> **Extensions** --> **Add extension** --> **Feeds**

3. Under *Subscribe to a feed*, add a URL from this [GTN feeds listing]({% link feeds/index.md %}). Make sure that it ends in `.xml`.  For example, `https://training.galaxyproject.org/training-material/topics/community/feed.xml` would provide updates on any community-tagged GTN materials into the Matrix channel.  

4. Under *Template*, change the existing text to the following: `$LINK: $SUMMARY`

5. Provide a reasonable name, and then hit **Subscribe**!

Details from Matrix are here: https://ems-docs.element.io/books/element-cloud-documentation/page/migrate-to-the-new-github-and-feeds-bots
