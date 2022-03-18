---
title: Knitting RMarkdown documents in RStudio
area: interactive tools
box_type: hands_on
layout: faq
contributors: [hexylena]
---

One of the other nice features of RMarkdown documents is making lovely presentation-quality worthy documents. You can take, for example, a tutorial and produce a nice report like output as HTML, PDF, or `.doc` document that can easily be shared with colleagues or students.

![Screenshot of the metadata with html_notebook and word_document being visible and a number of options controlling their output. TOC, standing for table of contents, has been set to true for both.]({% link topics/data-science/images/rstudio/r-meta-output.png %})

Now you're **ready to preview the document**:

![screenshot of preview dropdown with options like preview, knit to html, knit to pdf, knit to word]({% link topics/data-science/images/rstudio/r-preview.png %})

**Click Preview**. A window will popup with a preview of the rendered verison of this document.

![screenshot of rendered document with the table of contents on left, title is in a large font, and there are coloured boxes similar to GTN tutorials offering tips and more information]({% link topics/data-science/images/rstudio/r-preview.png %})

The preview is really similar to the GTN rendering, no cells have been executed, and no output is embedded yet in the preview document. But if you have run cells (e.g. the first few loading a library and previewing the `msleep` dataset:

![screenshot of the rendered document with a fancy table browser embedded as well as the output of each step]({% link topics/data-science/images/rstudio/r-preview-output.png %})

**When you're ready to distribute** the document, you can instead use the `Knit` button. This runs every cell in the entire document fresh, and then compiles the outputs together with the rendered markdown to produce a nice result file as HTML, PDF, or Word document.

![screenshot of the console with 'chunks' being knitted together]({% link topics/data-science/images/rstudio/r-knit.png %})

> ### {% icon tip %} Tip: PDF + Word require a LaTeX installation
> You might need to install additional packages to compile the PDF and Word document versions
{: .tip}

And at the end you can see a pretty document rendered with all of the output of every step along the way. This is a fantastic way to e.g. distribute read-only lesson materials to students, if you feel they might struggle with using an RMarkdown document, or just want to read the output without *doing* it themselves.

![screenshot of a PDF document showing the end of the tutorial where a pretty plot has been rendered and there is some text for conclusions and citations]({% link topics/data-science/images/rstudio/r-pdf.png %})
