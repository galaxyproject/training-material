---
layout: tutorial_hands_on

title: "Generating PDF artefacts of the website"
questions:
  - "How to generate PDF of the different tutorials and slides?"
objectives:
  - "Generating PDFs"
time_estimation: "10m"
key_points:
  - "PDFs can be easily generated for the different tutorials to share with learnees or to keep a fixed version of a tutorial"
contributors:
  - bebatut
---

# Introduction


The website with the training material can be run locally. Sometimes, it is also interesting to freeze the tutorials or to get PDFs of the tutorials.

> <agenda-title></agenda-title>
>
> In this tutorial, you will learn how to run a local instance of the GTN website:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Generate PDFs artifact

To generate the PDFs, a command `make pdf` is given. This command:

- Launches a detached Jekyll server to serve the website
- Generates the PDFs of the tutorials by calling Chrome via command line
- Generates the PDFs of the slide decks by calling decktape

> <hands-on-title>Checking the website generation locally</hands-on-title>
>
> 1. Install a browser
> 2. Generate the PDFs: `make pdf`
> 3. Check the generated PDFs in `_pdf` folder
{: .hands_on}
