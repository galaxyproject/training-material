# FAQs

Formerly known as snippets, FAQs provide answers to short, commonly asked questions. These can then be included inside tutorials where appropriate, and easily updated when needed. We also use them to create FAQ pages listing all FAQs on a given subject.

## Directory structure

FAQs may exist in different places in the repository, depending on the question:

- Project-level FAQs: `faqs/`
  - `faqs/galaxy/` for general Galaxy questions
  - `faqs/gtn/` for questions regarding the GTN website itself

- Topic-level FAQs: `topics/<topic>/faqs/`
  - for questions pertaining to that specific topic

- Tutorial-level FAQs: `topics/<topic>/tutorials/<tutorial>/faqs/`
  - for questions pertaining to that specific tutorial
  - if this is present, it is linked to from the tutorial overview box at the top, and from the end of the tutorial


## FAQ file

Each snippet (question) is a separate file, with some metadata

```yaml
---
title: How do I run a workflow?
area: workflows
box_type: tip        # tip/comment/hands_on; optional, if you want the content to be in a box
layout: faq          # if you set this the snippet will get its own page and be included in the FAQs page
---

Here you can write the answer in Markdown

- Go to `Workflows` on the top menu bar
- Click on ..
- ..

```


These FAQs can be included as snippets inside tutorials as well: `{% snippet faqs/galaxy/workflows_run.md %}`

If you would like to override the default box type, supply an include variable: `{% snippet faqs/galaxy/workflows_run.md box_type="hands_on" %}`
or to render without a box: `{% snippet faqs/galaxy/workflows_run.md box_type="none"  %}`


## FAQ page

We automatically generate a page containing all FAQs in that directory (and any possible subdirectories)

To do this, create a file named `index.md` inside the faq folder:

```yaml
---
layout: faq-page
---
```

If you would like to enforce an order of the faq areas, you can do so:

```yaml
---
layout: faq-page
area_order: [introduction, learners, instructors, contributors, other]
---
```

(just make sure you list all existing areas in the folder)



