---
title: Creating a GTN News post
area: contributors
layout: faq
box_type: tip
contributors: [shiltemann]
---

If you have created a new tutorial, running an event, published a paper around training, or have anything else interesting to share with the GTN community, we encourage you to write a **News item** about it!

News items will show up on the [GTN homepage]({% link index.md %}) and in the [GTN news feed]({% link news.md %}).

**Creating the news post**

Have a look at the existing news items in the [`news/_posts/` folder](https://github.com/galaxyproject/training-material/tree/main/news/_posts) of the GTN repository for some examples.

A news post is a markdown file that looks as follows:


```markdown
---
layout: news

title: "New Tutorial: My tutorial title"
tags:
  - new tutorial
  - transcriptomics
contributors:
  - shiltemann
  - hexylena

tutorial: "topics/introduction/tutorials/data-manipulation-olympics/tutorial.html"
cover: "path/to/cover-image.jpg"  # usually an image from your tutorial
coveralt: "description of the cover image"

---

A bit of text containing your news, this is all markdown formatted,
so you can do **bold** and *italic* text like this, and links look
like [this](https://example.com) etc.

Describe everything you want to convey here, can be as long as you
need.
```

Make sure the filename is structured as follows: `year-month-day-title.md`, so for example: `2022-10-28-my-new-tutorial.md`
