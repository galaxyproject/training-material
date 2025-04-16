---
layout: tutorial_hands_on
priority: 1
subtopic: sig
title: "Creating community content"
questions:
  - "What content does my Topic automatically get from the GTN?"
  - "How can I make sure this data is collected well?"
objectives:
  - "Signpost community leads and users to useful resources"
  - "Explain why metadata is key for such community resources"
  - "Provide a reference, rather than a tutorial"
time_estimation: "30M"
key_points:
  - "The GTN has worked hard to provide automated metrics and resources to highlight and acknowledge the efforts of communities and community leads"
  - "This only works if we contribute with effective tagging."
contributions:
  authorship:
  - nomadscientist
  - shiltemann
requirements:
  -
      type: "internal"
      topic_name: community
      tutorials:
          - sig_define
---

Galaxy *[Special Interest Group](https://galaxyproject.org/community/sig)* (**SIG**)s work hard to build and maintain training resources. The GTN has worked hard to acknowledge this and offer nice impact pages to communities.

Here is a list of resources that you can use!

> <comment-title></comment-title>
> - We want this material to grow (similar to the 'Creating content in markdown' tutorial) so please do add further resources that can help communities!
{: .comment}

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Maintainer Home

Maintaining resources is **fundamental** to the quality and usefulness of any software resource. We do *not* throw code over the wall!

To help topic maintainers to quickly recognise what materials need updating and fixing, topic **Maintainer Homes** were built.

> <hands-on-title>Go to your topic Maintainer Home</hands-on-title>
>
> 1. Go to any GTN training topic of interest.
> 2. Scroll down past the tutorials, and click on the "Maintainer Home" button
>
>    {% snippet topics/community/faqs/maintainer_home.md %}
>
> 3. Explore the Maintainer Home!
>    - For example, the [Single Cell Maintainer Home]({% link topics/single-cell/maintainer.md %})
>
{: .hands_on}

You may instantly see some key information missing from tutorials, or how long its been since someone checked it! Time to update some materials!

You can see an example from the Single-cell topic below:

<iframe
    src="/training-material/topics/single-cell/maintainer.html"
    width="100%"
    height="400"
    style="border: 2px solid #00008B;"
    title="Single-Cell Topic Maintainer Home">
</iframe>

# Community Home

Where the **Maintainer Home** helps you sustain your community, the **Community Home** helps you show off your community. An end-of-year gift in 2024, this page will sift through news, events, and GTN contributions for your community tag of interest (example: single-cell) and provide a beautiful visualization of your efforts.

> <hands-on-title>Go to your topic Community Home</hands-on-title>
>
> 1. Go to any GTN training topic of interest.
> 2. Scroll down past the tutorials, and click on the "Community Home" button
>
>    {% snippet topics/community/faqs/community_home.md %}
>
> 3. Explore the **Community Home**!
>    - For example, the [Single Cell Community Home]({% link topics/single-cell/community.md %})
>
{: .hands_on}

You can see an example from the Single-cell topic below:

<iframe
    src="/training-material/topics/single-cell/community.html"
    width="100%"
    height="400"
    style="border: 2px solid #00008B;"
    title="Single-Cell Topic Community Home">
</iframe>

# Topic usage statistics

Next up, you might want to know how many people are actually using your materials? Welcome to your **Topic usage statistics**! You may have already found this, actually, as it's (currently) at the bottom of the Maintainer Home.

> <hands-on-title>Go to your topic usage statistics</hands-on-title>
>
> 1. Go to any GTN training topic of interest.
> 2. Scroll down past the tutorials, and click on the "Maintainer Home" button
>
>    {% snippet topics/community/faqs/maintainer_home.md %}
>
> 3. Scroll down to the section "Statistics For Your Materials"
> 4. Explore the usage statistics!
>    - For example, the [Single Cell Usage Statistics]({% link topics/single-cell/maintainer.md %}#statistics-for-your-materials)
>
{: .hands_on}

You can see an example from the Single-cell topic below.

<iframe
    src="/training-material/topics/single-cell/maintainer.html#statistics-for-your-materials"
    width="100%"
    height="400"
    style="border: 2px solid #00008B;"
    title="Single-Cell Topic Usage Statistics">
</iframe>

# News widgets

You can also embed news into your pages, subdomains/ Galaxy Labs, or even your Matrix channels.

Follow this documentation to learn how:

1. [GTN Feeds]( {% link feeds/index.md %} )
2. Bot integration into matrix

{% snippet topics/community/faqs/matrix_news_bot.md %}

You can see an example from the Single-cell topic below.

<h3 class="mb-3">News and Events</h3>
<iframe width="100%" height="300px" src="/training-material/feeds/single-cell-month.w.html"></iframe>

# Workflow search

Want to see all the workflows tagged with your community tag across public servers? Look no further!

Follow this documentation to learn how:

1. [Galaxy Pan-Galactic Workflow Search]( {% link workflows/list.md %} )
2. [Workflow Search Querying]( {% link news/2023/11/20/workflow-search.md %} )

You can see an example from the Single-cell topic below.

<h3 class="mb-3">Public workflows</h3>
<iframe src="/training-material/workflows/embed.html?query=single-cell" height="400px" width="100%" class="gtn-embed" frameborder="1"></iframe>

# Galaxy Community Activities calendar

We host a GoogleCalendar on the [Galaxy Special Interest Group Hub page](https://galaxyproject.org/community/sig/).

{% snippet topics/community/faqs/community_activites_calendar.md %}

# Conclusion

{% icon congratulations %} Congratulations! You've made it to the end! Hopefully you think these resources are brilliant, and are making sure to tag everything (news, events, training materials, workflows, FAQs, you name it, you should tag it!) with your community tag!
