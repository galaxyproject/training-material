---
title: I can't find the tutorial I just added in GitPod/Locally.
area: jekyll
box_type: tip
layout: faq
contributors: [hexylena]
---

**Problem**: You've just added a new `tutorial.md` or `slides.html` to the GTN, but it doesn't show up when you browse the website which you're running locally on your laptop or in the cloud with GitPod.

**Solution**: `touch topics/<your-topic>/index.md`, and wait for the site to rebuild. Your page will show up

**Explanation**: when you add a new file to the site, Jekyll attempts to decide which pages need to be re-compiled. When running in `--incremental` mode (the default in GitPod), this is extremely important. However due to our use of plugins, and [the current inability of plugins to modify regeneration](https://github.com/jekyll/jekyll/issues/6418), the topic pages (which include mostly plugin-generated content) aren't regenerated correctly, and you need to manually touch the index file to trigger a rebuild.
