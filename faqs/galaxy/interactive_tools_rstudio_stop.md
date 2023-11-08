---
title: Stop RStudio
area: interactive tools
box_type: hands_on
layout: faq
contributors: [yvanlebras,shiltemann]
---


When you have finished your R analysis, it's time to stop RStudio.

1. First, save your work into Galaxy, to ensure reproducibility:
    1. You can use `gx_put(filename)` to save individual files by supplying the filename
    2. You can use `gx_save()` to save the entire analysis transcript and any data objects loaded into your environment.

2. Once you have saved your data, you can proceed in 2 different ways:
     - Deleting the corresponding history dataset named `RStudio` and showing a "in progress state", so yellow, OR
     - Clicking on the "User" menu at the top and go to "Active InteractiveTools" and locate the RStudio instance you started, selecting the corresponding box, and finally clicking on the "Stop" button at the bottom.
