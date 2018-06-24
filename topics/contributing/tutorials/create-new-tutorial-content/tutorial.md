---
layout: tutorial_hands_on
topic_name: contributing
tutorial_name: create-new-tutorial-content
---

# Introduction
{:.no_toc}

Once we have set up the infrastructure, we are ready to write the tutorial.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

The tutorial's content should be placed in the file `tutorial.md`. Its syntax and structure are simple, and will have the following structure:

```markdown
---
layout: tutorial_hands_on
topic_name: training
tutorial_name: create-new-tutorial
---
# Introduction
{:.no_toc}

blabla

# Section 1

blabla

## Subsection 1

blabla

# Section 2

blabla

## Subsection 2

blabla

# Conclusion
{:.no_toc}

blabla
```

# Metadata

The `tutorial.md` needs to start with some metadata at the top:

- `layout: tutorial_hands_on`: keep the default
- `topic_name: training`: replace 'training' the name of the topic
- `tutorial_name: create-new-tutorial`: replace 'create-new-tutorial' with the `name` of tutorial that you used in the topic level metadata

These metadata are there to help the templating system linking between the tutorial's file and the global [metadata]({{site.baseurl}}/topics/contributing/tutorials/create-new-tutorial-metadata/tutorial.html).
If this is not correctly defined, the tutorial will not be shown on the website.

> ### {% icon hands_on %} Hands-on: Fix the top metadata
>
> 1. Change the `tutorial_name` and the `topic_name` to fit to the ones defined in the metadata
> 2. (Optional) Build the website locally by following the [Jekyll tutorial]({{ site.baseurl }}/topics/contributing/tutorials/running-jekyll/tutorial.html) and check that the tutorial is reachable
{: .hands_on}

# Content

The tutorial's content is written directly after the short section of metadata. This is written in Markdown, a simple markup language.

> ### {% icon tip %} Tip: Markdown
>
> Check [this cheatsheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet) to learn more how to use Markdown.
{: .tip}

The Markdown content is then transformed into a user friendly webpage throughout a templating system. With this approach, there is no need to add the name of every tutorial each time, since they are automatically added based on the tutorial's metadata.

We recommend to structure the tutorials as follows:

- An introduction, to bring an overview of the tutorial with its use cases, data, and methods
- Multiple sections, representing the steps of the analysis, complete with their hands-on parts (practicing is an important part of the learning process)
- A conclusion to summarize what has been done in the tutorial (with a graphic)

> ### {% icon hands_on %} Hands-on: Structuring the tutorial
>
> 1. Add a small description of the dataset
> 2. Add one or two sections with ideas for the tutorial
> 3. Add a small conclusion
{: .hands_on}

## Adding images with caption

To add an image in Markdown file, we need to use the markdown syntax for this: `![](../../images/image.png)`.

We have also added a small plugin to add a caption for each image:

![This figure shows an example of an image with a caption](../../images/image_caption_screenshot.png "Example of an image with a caption")

The prefix "Figure 1." is automatically added before its caption. This is done with the following Markdown syntax:

```markdown
![A textual description of the image](../images/image.png "This is my super caption")
```

We can also cross-reference images inside our Markdown with an anchor. For example, we can link to [the previous figure](#figure-1) using `[the display text](#figure-nb)` (changing `nb` with the image's number).

## Writing mathematical expressions

Mathematical expressions can be written in LaTeX, and are automtaically rendered with [MathJax](https://www.mathjax.org/).

Surround your math expression with two `$` signs (like in LaTeX math blocks):

- inline expressions, *e.g.* `$$ 5 + 5 $$` will be rendered as $$ 5 + 5 $$
- block expressions, *e.g.* `$$ 5 + 5 $$` will be rendered in its own line block as

   $$ 5 + 5 $$

Dollar signs are therefore *reserved characters* for instructing the templating system to open/close LaTeX math blocks. If you want to use a `$` within your expression, you will need to *escape* it: `$$ a + 3\$ = 5\$ $$` will be rendered as: $$ a + 3\$ = 5\$ $$

> ### {% icon comment %} Comments
> LaTeX code that uses the pipe symbol `|` in inline math statements may lead to a line being recognized as a table line by the templating system.
> This can be avoided by using the `\vert` command instead of `|`
{: .comment}

# Improving the learning experience

To improve the learning experience in our tutorial, we define some boxes to highlight content.

These boxes are defined always with the same structure:

{% raw %}
```markdown
> ### {% icon an_icon %} Type of box: Name of the box
> list
{: .type_of_box}
```
{% endraw %}

You must follow this structure exactly for it to be rendered correctly.

## **Overview** box

This box at the top of each tutorial is automatically generated using the metadata we defined in the topic's metadata file

> ### {% icon hands_on %} Hands-on: Checking the metadata
>
> 1. Check that the metadata added previously are correctly filling the overview box
>
>    > ### {% icon question %} Questions
>    >
>    > What metadata hasn't been added to this box?
>    >
>    >    > ### {% icon solution %} Solution
>    >    >
>    >    > The take-home messages are not added to this box but into the last box of the tutorial
>    >    {: .solution}
>    {: .question}
{: .hands_on}

## **Agenda** box

In most tutorials, the second box is the agenda box, placed at the end of the introduction. It shows the table of contents for the tutorial

{% raw %}
```markdown
> ### Agenda
>
> In this tutorial we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}
```
{% endraw %}

There is no need to fill out the list; this will be done automatically based off of your tutorial's section title.

To avoid adding the "Introduction" and "Conclusion" sections in the agenda, you can add `{:.no_toc}` below the section name. This will be rendered as follows:

> ### Agenda
>
> In this tutorial we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## **Hands-on** box

We find that having users walk through the tutorial, doing all of the steps is important for learning the concepts taught. We therefore emphasize this by regularly adding hands-on sections, where trainees are encouraged to do the analysis by themselves. We have designed some special boxes to make these sections easier to find.

{% raw %}
```markdown
> ### {% icon hands_on %} Hands-on: Spliced mapping
>
> 1. **RNA STAR** {% icon tool %}: Map your reads on the reference genome with
>    - *"Single-end or paired-end reads"*:  `Paired-end (as individual datasets)`
>    - *"RNA-Seq FASTQ/FASTA file, forward reads"*: the generated `trimmed reads pair 1` files (multiple datasets)
>    - *"RNA-Seq FASTQ/FASTA file, reverse reads"*: the generated `trimmed reads pair 2` files (multiple datasets)
>    - *"Custom or built-in reference genome"*: `Use a built-in index`
>    - *"Reference genome with or without an annotation"*: `use genome reference without builtin gene-model`
>    - *"Select reference genome"*: `Drosophila Melanogaster (dm6)`
>    - *"Gene model (gff3,gtf) file for splice junctions"*: the imported `Drosophila_melanogaster.BDGP6.87.gtf`
>    - *"Length of the genomic sequence around annotated junctions"*: `36`
>
>        This parameter should be length of reads - 1
>
> 2. **MultiQC** {% icon tool %}: Aggregate the STAR logs with
>      - *"Which tool was used generate logs?"*: `STAR`
>      - *"Type of FastQC output?"*: `Log`
>      - *"STAR log output"*: the generated `log` files (multiple datasets)
{: .hands_on}
```
{% endraw %}

For consistency please use:

- {% raw %}`{% icon hands_on %}`{% endraw %} emoji to define that is an hands-on
- Short imperative sentences to make it easy to identify the tasks
- Name of the tool in bold with the {% raw %}`{% icon tool %}`{% endraw %} emoji to make it easy to identify a Galaxy tool
- Parameters for the tool as a sublist

This will be rendered like:


> ### {% icon hands_on %} Hands-on: Spliced mapping
>
> 1. **RNA STAR** {% icon tool %}: Map your reads on the reference genome with
>    - *"Single-end or paired-end reads"*:  `Paired-end (as individual datasets)`
>    - *"RNA-Seq FASTQ/FASTA file, forward reads"*: the generated `trimmed reads pair 1` files (multiple datasets)
>    - *"RNA-Seq FASTQ/FASTA file, reverse reads"*: the generated `trimmed reads pair 2` files (multiple datasets)
>    - *"Custom or built-in reference genome"*: `Use a built-in index`
>    - *"Reference genome with or without an annotation"*: `use genome reference without builtin gene-model`
>    - *"Select reference genome"*: `Drosophila Melanogaster (dm6)`
>    - *"Gene model (gff3,gtf) file for splice junctions"*: the imported `Drosophila_melanogaster.BDGP6.87.gtf`
>    - *"Length of the genomic sequence around annotated junctions"*: `36`
>
>        This parameter should be length of reads - 1
>
> 2. **MultiQC** {% icon tool %}: Aggregate the STAR logs with
>      - *"Which tool was used generate logs?"*: `STAR`
>      - *"Type of FastQC output?"*: `Log`
>      - *"STAR log output"*: the generated `log` files (multiple datasets)
{: .hands_on}

There are also some predefined **parameter icons** available which can be used to indicate the *type* of parameter. These are optional,
but may be helpful in some cases (for example to distinguish between single file inputs and collection inputs).

The available icons are:

{% raw %}
```markdown
> ### {% icon hands_on %} Hands-on: My Step
>
> 1. **My Tool** {% icon tool %} with the following parameters
>  - {% icon param-text %} *"My text parameter"*: `my value`
>  - {% icon param-file %} *"My input file"*: `my file`
>  - {% icon param-files %} *"My multiple file input or collection"*: `my collection`
>  - {% icon param-select %} *"My select menu"*: `my choice`
>  - {% icon param-check %} *"My check box"*: `yes`
{: .hands_on}
```
{% endraw %}

which, when rendered, look like:

> ### {% icon hands_on %} Hands-on: My Step
>
> 1. **My Tool** {% icon tool %} with the following parameters
>  - {% icon param-text %} *"My text parameter"*: `my value`
>  - {% icon param-file %} *"My input file"*: `my file`
>  - {% icon param-files %} *"My multiple file input or collection"*: `my collection`
>  - {% icon param-select %} *"My select menu"*: `my choice`
>  - {% icon param-check %} *"My check box"*: `yes`
{: .hands_on}


## **Questions** and **solution** boxes

Questions can be added to force trainees to think about what they are currently doing, and to put things in perspective.
They can also help the instructors by exposing and clarifying common scenarios, errors, or applications.

{% raw %}
```markdown
> ### {% icon question %} Questions
>
> 1. Why are some tests filtered?
> 2. Does it improve the *p*-value distribution?
>
>    > ### {% icon solution %} Solution
>    >
>    > 1. Sol for the first question
>    > 2. Sol for the second question
>    >
>    {: .solution}
{: .question}
```
{% endraw %}

Which will be rendered as:

> ### {% icon question %} Questions
>
> 1. Why are some tests filtered?
> 2. Does it improve the *p*-value distribution?
>
>    > ### {% icon solution %} Solution
>    >
>    > 1. Sol for the first question
>    > 2. Sol for the second question
>    >
>    {: .solution}
{: .question}

Questions should be quick to answer. You can directly ask a question and expect an answer, or you can provide some answers and create multiple choices questions (MCQs).
With well chosen wrong answers, MCQs can do much more than just measure how much someone knows, such as exposing common misconceptions and mistakes.

In the box below, initially hidden, we add the correct answer and possibly any additional explanation. Self-trainees can then check the solution and its explanation.


## **Tips** box

{% raw %}
```markdown
> ### {% icon tip %} Tip: Importing data via links
>
> * Copy the link location
> * Open the Galaxy Upload Manager
> * Select **Paste/Fetch Data**
> * Paste the link into the text field
> * Press **Start**
{: .tip}
```
{% endraw %}

Rendered:

> ### {% icon tip %} Tip: Importing data via links
>
> * Copy the link location
> * Open the Galaxy Upload Manager
> * Select **Paste/Fetch Data**
> * Paste the link into the text field
> * Press **Start**
{: .tip}

## **Comments** boxes

{% raw %}
```markdown
> ### {% icon comment %} Comments
> - Edit the "Database/Build" to select "dm3"
> - Rename the datasets according to the samples
{: .comment}
```
{% endraw %}

Rendered:

> ### {% icon comment %} Comments
> - Edit the "Database/Build" to select "dm3"
> - Rename the datasets according to the samples
{: .comment}

## **Details** box

The detail box is used to give more background explanation on the subject. By default the box is collapsed, trainees can expand it if they wish to know extra information about a topic.

{% raw %}
```markdown
> ### {% icon details %} More details on the ....
>
> Add more details in Markdown...
>
{: .details}
```
{% endraw %}

Rendered:

> ### {% icon details %} More details on the ....
>
> Add more details in Markdown...
>
{: .details}

## **Key points** box

This last box of the tutorial is automatically created with the take-home messages defined in the topic's metadata

To render the boxes correctly, the syntax needs to be correct. If it doesn't work, have a look at similar tutorials and get inspiration.

## Nested boxes

Boxes can be nested, *e.g.* for having tips inside a hands-on:

{% raw %}
```markdown
> ### {% icon hands_on %} Hands-on: Defining the topic for the tutorial
>
> 1. Search for NCBI Blast+ on the [ToolShed](https://toolshed.g2.bx.psu.edu/)
> 2. Check in which category it is
>
>    > ### {% icon question %} Questions
>    >
>    > In which topic will you put the tutorial?
>    >
>    >    > ### {% icon solution %} Solution
>    >    >
>    >    > If we search for [NCBI Blast+ in the ToolShed](https://toolshed.g2.bx.psu.edu/view/devteam/ncbi_blast_plus/7538e2bfcd41), it is attributed to 2 categories (bottom): "Next Gen Mappers" and "Sequence Analysis".
>    >    > We decided to put it in "Sequence analysis" because this is the most general one for this tutorial.
>    >    {: .solution}
>    {: .question}
{: .hands_on}
```
{% endraw %}

# Conclusion
{:.no_toc}
