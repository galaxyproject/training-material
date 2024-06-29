---
title: Time to git commit
area: gat
box_type: hands_on
layout: faq
contributors: [hexylena]
required_parameters:
  page: The parent page object since we reference the page's title in the commit message
examples:
  Using this snippet:
    page: page
---

It's time to commit your work! Check the status with

```
git status
```

Add your changed files with

```
git add ... # any files you see that are changed
```

And then commit it!

```
git commit -m 'Finished {{ include.page.title }}'
```
