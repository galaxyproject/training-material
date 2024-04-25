---
layout: tutorial_hands_on

title: "Creating Slides in Markdown"
questions:
  - "How to format slides?"
  - "How do we add presenter notes?"
  - "How to use the features of the slide show tool?"
  - "What sort of content should be included in slides?"
objectives:
  - "Create a new set of slides"
  - "Add presenter comments"
time_estimation: "15m"
key_points:
  - "Slides are often presented before the hands-on portion, provide relevant background information"
  - "Provide examples that are relevant to your audience"
  - "Consider writing speaker notes"
subtopic: slides
contributions:
  authorship:
  - bebatut
  - hexylena
  editing:
  - bgruening
  - shiltemann
---

Once we have set up the infrastructure, we are ready to write a slide deck.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

The tutorial's content should be placed in the file `slides.html`. Its syntax and structure are simple, and will have the following structure:

```markdown
---
layout: tutorial_slides
logo: "GTN" # You can also provide an image path, e.g. topics/ai4life/images/AI4Life-logo_giraffe-nodes.png

title: Title of the Slide Deck
questions:
- Which biological questions are addressed by the tutorial?
- Which bioinformatics techniques are important to know for this type of data?
objectives:
- The learning objectives are the goals of the tutorial
- They will be informed by your audience and will communicate to them and to yourself
  what you should focus on during the course
- They are single sentences describing what a learner should be able to do once they
  have done completed tutorial
- You can use Bloom's Taxonomy to write effective learning objectives
time_estimation: ''
key_points:
- The take-home messages
- They will appear at the end of the tutorial

contributions:
  authorship:
  - contributor1
  editing:
  - contributor2
---

# Slides

- Introduce a subject
- Quick overview
- Not too much text
- Images are Nice

???

- Slides a great way to introduce an audience to your subject.

---

## Another Slide

And some contents

---

## A third slide

More content!

???

- Here go speaker notes
```

# Metadata

The `slides.html` needs to start with some metadata at the top:

- `layout: tutorial_slides`: please keep the default unless 
- `title`: title of the tutorial (it will appear on the tutorial page and the topic page)
- `contributions`: everybody who has contributed to this tutorial (usernames must match those in `CONTRIBUTORS.yaml` file)

> <hands-on-title>Fill the basic metadata</hands-on-title>
>
> 1. Update the slide information in the header section of your slides:
>
>     ```
>     title: "Similarity search with BLAST"
>     ```
>
{: .hands_on}

This information is used to display the data from the topic and tutorial page. They are also used to check which information is missing for the tutorials.

We also define metadata related to the pedagogical content of the tutorial, which will appear at the top ("Overview" box) and bottom of the online tutorial:

{% assign kid_key = "Slides Schema" %}
{% assign kid_val = site.data['schema-slides'] %}
{% include _includes/schema-render.html key=kid_key value=kid_val %}


> <hands-on-title>Fill out the pedagogical metadata</hands-on-title>
>
> 1. Define 2 questions that will be addressed during the tutorial and add them to the metadata
> 2. Define 2 learning objectives for the tutorial and add them to the metadata
{: .hands_on}

## Listing contributors

All tutorials and slides must give credit to all contributors. This can be any type of contribution, adding them in GitHub, creating images for it, etc.

1. Make sure all contributors are listed in the [`CONTRIBUTORS.yaml`](https://github.com/galaxyproject/training-material/blob/main/CONTRIBUTORS.yaml) file.
   Each contributor is defined in this file like:

   ```yaml
   contributor-username:                 # GitHub username (if the contributor has one)
     name: Full Name                     # mandatory
     joined: 2020-06                     # mandatory
     email: saskia.hiltemann@gmail.com   # optional
     twitter: shiltemann                 # optional
     linkedin: shiltemann                # optional
     gitter: shiltemann                  # optional
     orcid: 0000-0003-3803-468X          # optional
     bio: Researcher at EMC              # optional
     elixir_node: AU                     # optional
     contact_for_training: false         # optional
     affiliations: [gallantries]         # optional
   ```

2. Add all contributors to the metadata of the slide deck.

   ```yaml
   contributions:
     authorship:
       - contributor-username
     editing:        # Optional
       - hexylena
     funding:        # Optional
       - carpentries
     testing:        # Optional
       - userX
     ux:             # Optional
       - userY
     infrastructure: # Optional
       - userZ
   ```

   To define a funding body in the `CONTRIBUTORS.yaml` there are a few extra fields available:

   ```yaml
   gallantries:
     name: Gallantries Project
     joined: 2020-09
     avatar: "https://www.erasmusplus.nl/assets/images/logo.png"
     github: false
     funder: true
     funding_id: 2020-1-NL01-KA203-064717
     funding_statement: |
        This project ([`2020-1-NL01-KA203-064717`](https://erasmus-plus.ec.europa.eu/projects/search/details/2020-1-NL01-KA203-064717)) is funded with the support of the Erasmus+ programme of the European Union. Their funding has supported a large number of tutorials within the GTN across a wide array of topics.
        ![eu flag with the text: with the support of the erasmus programme of the european union](https://gallantries.github.io/assets/images/logosbeneficaireserasmusright_en.jpg)
   ```

   Funding bodies will be credited at the bottom of the tutorial with the appropriate funding statement, and will get a page in the hall of fame listing all tutorials that list them as a funder.

   For an example of how this all looks, see the [R basics tutorial]({% link topics/data-science/tutorials/r-basics/tutorial.md %}) (top and bottom of the tutorial).

{% assign kid_key = "Contributions Schema" %}
{% assign kid_val = site.data['schema-tutorial']['mapping']['contributions'] %}
{% include _includes/schema-render.html key=kid_key value=kid_val %}

{% assign kid_key = "CONTRIBUTORS Schema" %}
{% assign kid_val = site.data['schema-contributors'] %}
{% include _includes/schema-render.html key=kid_key value=kid_val %}

# Content

The tutorial's content is written directly after the section of metadata. This is written in Markdown, a simple markup language.

> <comment-title>Markdown</comment-title>
>
> Check [this cheatsheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet) to learn more how to use Markdown.
{: .comment}

The Markdown content is then transformed into a user friendly webpage through a templating system. With this approach, there is no need to add the name of every tutorial each time, since they are automatically added based on the tutorial's metadata.

## Formatting

Please see the [associated slide deck]({% link topics/contributing/tutorials/create-new-tutorial-slides/slides.html %}) for details on the formatting commands available. In general the markdown you already know will continue to work, but there are additional slide formatting commands available (e.g. `.pull-left[]` for classic two column layouts.)

## Automated Video Recordings

Please see the [dedicated tutorial on adding video support]({% link topics/contributing/tutorials/slides-with-video/tutorial.md %}) on adding automatic video recordings to your slides, if you're interested in that feature.

## Adding images with captions

To add an image in Markdown file, we need to use the markdown syntax for this: {% raw %}`!​[proper alt text describing the image for visually impaired learners](../../images/image.png)`{% endraw %}.

We have also added a small plugin to handle captions for each image:

![A textual description of the image](../../images/image_caption_screenshot.png "Example of an image with a caption ")<!-- Adding a space to the caption to not trigger figurigy skip_titles -->

The prefix "Figure 1." is automatically added before its caption. This is done with the following Markdown syntax:

{% raw %}
```markdown
!​[A textual description of the image](../images/image.png "Example of an image with a caption")
```
{% endraw %}

### Guidelines on Alt vs Figcaption Text

> While both the alt attribute and the figcaption element provide a way to
> describe images, the way we write for them is different. **`alt` descriptions
> should be functional; `figcaption` descriptions should be editorial or
> illustrative.**
>
{: .quote cite="https://thoughtbot.com/blog/alt-vs-figcaption" author="thoughtbot.com"}

As an example for this image:

![alt text]({% link topics/microbiome/images/plasmid-metagenomics-nanopore/sequence_method.jpg %} "Example of an image with a caption ")

{% raw %}
```markdown
!​[Alt text (shown when image cannot be displayed)]({% link path/to/image.png %} "Example of an image with a caption")
```
{% endraw %}

Field          | Appropriate Contents
----           | -----
alt text       | Image of cell membrance with an embedded protein with central pore. DNA is shown splitting and entering the pore, an electrical signal comes out reading A C T or G.
figure caption | Using nanopore sequencing, a single molecule of DNA or RNA can be sequenced without the need for PCR amplification or chemical labeling of the sample. (Image from: <a href="https://nanoporetech.com/sites/default/files/s3/white-papers/WGS_Assembly_white_paper.pdf?submissionGuid=40a7546b-9e51-42e7-bde9-b5ddef3c3512">Nanopore sequencing: The advantages of long reads for genome assembly</a>)

## Writing mathematical expressions

Mathematical expressions can be written in LaTeX, and are automatically rendered with [MathJax](https://www.mathjax.org/).

Surround your math expression with two `$` signs on each side (like in LaTeX math blocks):

- inline expressions, *e.g.* `$$ 5 + 5 $$` will be rendered as $$ 5 + 5 $$
- block expressions, *e.g.* `$$ 5 + 5 $$` will be rendered in its own line block as

   $$ 5 + 5 $$

Dollar signs are therefore *reserved characters* for instructing the templating system to open/close LaTeX math blocks. If you want to use a `$` within your expression, you will need to *escape* it: `$$ a + 3\$ = 5\$ $$` will be rendered as: $$ a + 3\$ = 5\$ $$


LaTeX code that uses the pipe symbol `|` in inline math statements may lead to a line being recognized as a table line by the templating system.
This can be avoided by using the `\vert` command instead of `|`

## Tables and Matrices

Tables can be generated using markdown by using the `|` symbol to indicate column dividers, and `--` for table headers:

{% raw %}
```markdown
|       | Obs1 | Obs2 | Obs3 |
|------ |--------------------|
| Feat1 | 0    | 1    | 2    |
| Feat2 | 1    | 2    | 3    |
| Feat3 | 2    | 3    | 4    |
```
{% endraw %}

When rendered, they will take the full width of the page:

|       | Obs1 | Obs2 | Obs3 |
|------ |--------------------|
| Feat1 | 0    | 1    | 2    |
| Feat2 | 1    | 2    | 3    |
| Feat3 | 2    | 3    | 4    |

This does not appear to be visually appealing when representing matrices, which is why a matrix box can be used instead:

{% raw %}
```markdown
> |       | Obs1 | Obs2 | Obs3 |
> | ----- |--------------------|
> | Feat1 | 0    | 1    | 2    |
> | Feat2 | 1    | 2    | 3    |
> | Feat3 | 2    | 3    | 4    |
{: .matrix}
```
{% endraw %}

The rendered table is then given as a minimum-width and centred matrix:

> |       | Obs1 | Obs2 | Obs3 |
> | ----- |--------------------|
> | Feat1 | 0    | 1    | 2    |
> | Feat2 | 1    | 2    | 3    |
> | Feat3 | 2    | 3    | 4    |
{: .matrix}

## Internally linking to other training material

If you want to link to other training material within your text, please use the {%raw%}`{​% link path/to/file.ext %​}`{%endraw%} tag:

{%raw%}
```markdown
[link text]( {% link topics/single-cell/tutorials/scrna-case_monocle3-trajectories/tutorial.md %} )
```
{%endraw%}

(Note the `.md` extension, and not `.html`, always provide the file name here, it will automatically be converted to the correct link)

If you want to link to a specific section in a tutorial using an anchor (e.g. `#getting-started`), place it outside of the {%raw%}`{​% link %​}`{%endraw%} tag:

{%raw%}
```markdown
[link text]({% link topics/single-cell/tutorials/scrna-case_monocle3-trajectories/tutorial.md %}#section-name)
```
{%endraw%}


# Improving the learning experience with Boxes

To improve the learning experience in our tutorial, we define some boxes to highlight content.

These boxes are defined always with the same structure:

{% raw %}
```markdown
> <type​-title>Name of the box</type​-title>
> text goes here!
{: .type}
```
{% endraw %}

where type is something like `tip`, so `<tip-title>` and `{: .tip}`. You must follow this structure **exactly** for it to be rendered correctly.

# Conclusion

If you have created a new slide deck, please also consider writing a [GTN news post]({% link faqs/gtn/gtn_news_create_post.md %}) about it to let people know about it (and make it easy for us to tweet about)!


