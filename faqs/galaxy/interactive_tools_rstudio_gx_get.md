---
title: Use gx_get
area: interactive tools
box_type: hands_on
layout: faq
contributors: [Camila-goclowski]
---

><tip-title>How to import datasets from Galaxy to RStudio</tip-title>
> RStudio in Galaxy comes with a `gx_get()` function. This is critical for moving datasets from your Galaxy history into Galaxy RStudio. The function outputs the file path from which you access your data via RStudio.
>
> To use it,  use the numbered position of the dataset you are looking to import. For example:
>
> If we want to find the first dataset in your history, run the following command:
> ```r
> gx_get(1)
> ```
>The result of this command will be a file path to the first dataset in your Galaxy history. Use that file path for future importing commands.
{: .tip}
