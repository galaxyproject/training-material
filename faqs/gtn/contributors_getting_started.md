---
title: How can I get started with contributing?
area: contributors
layout: faq
box_type: tip
contributors: [bebatut,shiltemann]
---

If you would like to get involved in the project but are unsure where to start, there are some easy ways to contribute which will also help you familiarize yourself with the project.

#### 1. Checking existing tutorials

A great way to help out the project is to test/edit existing tutorials. Pick a tutorial and check the contents. Does everything work as expected? Are there things that could be improved?

Below is a checklist of things to look out for to help you get started. If you feel confident in making changes yourself, please open a pull request, otherwise please file an issue with any problems you run into or suggestions for improvements.

*Basic*
- **Test** the tutorial on a running Galaxy instance
   - For example [Galaxy Main](https://usegalaxy.org), [Galaxy Europe](https://usegalaxy.eu), or [Galaxy Australia](https://usegalaxy.org.au)
   - Report any issues you run into
- **Language** editing
  - Fix spelling and grammar mistakes
  - Simplify the English (to make it more accessible)

*Intermediate*
- **Metadata**
  - Are the objectives, keypoints and time estimate filled in?
  - Do they fit with the contents of the tutorial?
- **Content**
  - Is there enough background information provided in the introduction section and throughout the tutorial?
  - **Question boxes**
    - Add questions or question boxes where you think they might be useful (make people think about results they got, test their understanding, etc)
    - Check that answers are still up-to-date
  - **Screenshots and Videos**
    - Make sure there is also a textual description of the image/video contents
    - Does the screenshot add value to the tutorial or can it be removed?

*Advanced*
- **Workflows**
  - Add a workflow definition file `.ga` if none is present
  - Check that the existing workflow is up-to-date with the tutorial contents
  - Enable [workflow testing](https://github.com/usegalaxy-eu/workflow-testing)
- **Tours**
  - Add a tour if none exists
  - Run the existing tour and check that it is up-to-date with the tutorial contents
- **Datasets**
  - Check that all datasets used in the tutorial are present in Zenodo
  - Add a data-library.yaml file if none exists


#### 2. Reviewing pull requests

Another great way to help out the project is by reviewing [open pull requests]({{ site.github_repository }}/pulls?q=is%3Apr+is%3Aopen+sort%3Aupdated-desc). You can use the above checklist as a guide for your review. Some documentation about how to add your review in the GitHub interface can be found [here](https://help.github.com/articles/about-pull-request-reviews/)


