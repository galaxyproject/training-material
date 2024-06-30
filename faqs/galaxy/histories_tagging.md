---
title: History tagging
description: Explains how to add tags to a history
area: histories
box_type: tip
layout: faq
contributors: [nekrut]
---

Tags are short pieces of text used to describe the thing they're attached to and many things in Galaxy can be tagged.
Each item can have many tags and you can add new tags or remove them at any time. Tags can be another useful way to
organize and search your data. For instance, you might tag a history with the type of analysis you did in it: `assembly`
or `variants`. Or you may tag them according to data sources or some other metadata: `long-term-care-facility` or
`yellowstone-park:2014`.

**To tag a history:**

1. Click on {% icon galaxy-pencil %} (**Edit**) next to the history name (which by default is "Unnamed history").
2. Click on **Add tags** {% icon galaxy-tags %} and start typing. Any tags that you've used previously will show below your partial entry -
  allowing you to use this 'autocomplete' data to re-use your previous tags without typing them in full.
3. Click on **Save** {% icon galaxy-save %}.
4. To cancel, click the {% icon galaxy-undo %} "Cancel" button.

> <warning-title>Do not use spaces</warning-title>
> It is strongly recommended to replace spaces in tags with `_` or `-`, as spaces will automatically be removed when the tag is saved.
{: .warning}

![UI for tagging histories]({% link shared/images/history_tags.png %})