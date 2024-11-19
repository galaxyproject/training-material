---
title: What is my.galaxy.training
area: introduction
layout: faq
box_type: tip
contributors: [hexylena]
---

The [`my.galaxy.training`](https://my.galaxy.training) is part of the GTN. We found that often need to direct our learners to specific pages within Galaxy, but which Galaxy? Should we add three links, one for each of the current bigger UseGalaxy.* servers? That would be really annoying for users who aren't using one of those servers.

E.g. how do we link to [/user](https://my.galaxy.training/?path=/user), the user preferences page which is available on every Galaxy Instance? This service handles that in a private and user-friendly manner.

### (Learners) How to Use It

When you access a my.galaxy.training page you'll be prompted to select a server, simply select one and you're good to go!

If you want to enter a private Galaxy instance, perhaps, behind a firewall, that's also an option! Just select the 'other' option and provide your domain. Since the redirection happens in your browser with no servers involved, as long as *you* can access the server, you'll get redirected to the right location.

### (Tutorial Authors) How to use it

If you want to link to a specific page within Galaxy, simple construct the URL: `https://my.galaxy.training/?path=/user` where everything after `?path` is the location they should be redirected to on Galaxy. That example link will eventually redirect the learner to something like `https://usegalaxy.eu/user`.

### Technical Background

So we took inspiration from [Home Assistant](https://my.home-assistant.io/) which had the same problem, how to redirect users to pages on their own servers. The `my.galaxy.training` service is a very simple static page which looks in the user's [`localStorage`](https://developer.mozilla.org/en-US/docs/Web/API/Window/localStorage) for their preferred server.
If it's not set, the user can click one of the common domains, and be redirected. When they access another link, they'll be prompted to use a button that remembers which server they chose.

### Data Privacy

Any domain selected is not tracked nor communicated to any third party. Your preferred server is stored in your browser, and never transmitted to the GTN. That's why we use localStorage instead of cookies.
