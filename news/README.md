# News items

To create a new item for the GTN newsfeed, add a markdown file in the `_posts` directory of this folder


```
---
title: The headline for your news item
contributors: [contributor1,contributor2]    # authors of the post or tutorial etc
tags: [keyword1, keyword2]   # optional
cover: <link to cover image> # optional
coveralt: # Alt text describing your cover image.
          # Mandatory if you set a cover image.
external: true|false         # optional, set to "true" if you just want to link to an external website
link: "https://<yourlink>"   # optional, if you would like to link to an external post or page
                             # mandataory if you set external: true above
tutorial: topics/<topic>/tutorials/<tutorial>/tutorial.html # if this newsitem is about a GTN tutorial
layout: news
---


The content of your news post goes here (unless you set external: true in the metadata)




```

- `tags`: for example: `new tutorial`, `transcriptomics`, `galaxy`, etc
- `cover`: if you dont supply a cover image, we have defined some defaults. E.g. the "galaxy news" logo if one of your tags is `galaxy`, or the GTN logo otherwise

