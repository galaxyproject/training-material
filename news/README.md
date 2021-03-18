# News items

To create a new item for the GTN newsfeed, add a markdown file in the `_posts` directory of this folder


```
---
title: The headline for your news item
contributors: [] #authors of the post or tutorial etc
tags: [keyword1, keyword2]  # optional
cover: <link to cover image> # optional
tutorial: topics/<topic>/tutorials/<tutorial>/tutorial.html # if this newsitem is about a GTN tutorial
link: "https://<yourlink>" #optional, if you would like to link to an external post or page
layout: news
---


The content of your news post goes here




```

- `tags`: for example: `new tutorial`, `transcriptomics`, `galaxy`, etc
- `cover`: if you dont supply a cover image, we have defined some defaults. E.g. the "galaxy news" logo if one of your tags is `galaxy`, or the GTN logo otherwise

