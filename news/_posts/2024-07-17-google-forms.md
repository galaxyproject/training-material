---
title: Simplifying GTN contribution with Google Forms
layout: news
tags:
- gtn infrastructure
- new feature
- automation
contributions:
  authorship:
  - hexylena
  - shiltemann
  infrastructure:
  - shiltemann
  - hexylena
cover: news/images/2024-google-forms.png
coveralt: A screenshot of the GTN Google Form for news contributions. The form is titled 'GTN News' and has fields for 'Title', with a screenshot of a galaxy single cell news post as the header image.
---

Since the last time we [announced GTN news posting via Google Form]({{ site.baseurl }}/news/2024/01/29/simplified-gtn-news-submission-via-google-form.html), we've found that this has been an excellent fit for the community and decided to greatly expand the use of Google Forms for GTN contributions.

Our users were clear: *GitHub, pull requests, commits, and YAML files are not always a pleasant or user-friendly way to contribute.*

After the success of news contributions via form we have expanded that process significantly. Now you can contribute the following all via Google Form:

- [News](https://forms.gle/TqGTr6y46wrJDri7A)
- [FAQs](https://forms.gle/2JVMfd1AgtenZPvv9)
- [Recordings](https://forms.gle/qNG8FkTN1yRZPNZY6)
- [Events](https://forms.gle/M6ECp1e3pZoFGYnV8)

## The Details

Each of these works similarly: you fill out the form, and then in the background a GTN GitHub Action will create a pull request for you.
The GTN maintainers will review the pull request, and when everything is in order and tests pass, merge it. Your contribution will be in the GTN without you needing to learn development tooling.
