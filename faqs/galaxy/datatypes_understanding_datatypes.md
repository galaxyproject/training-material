---
title: Understanding Datatypes
area: datatypes
layout: faq
box_type: tip
contributors: [jennaj, garimavs]
---

- Allow Galaxy to detect the datatype during Upload, and adjust from there if needed.
- Tool forms will filter for the appropriate datatypes it can use for each input. 
- [Directly changing](https://training.galaxyproject.org/training-material/faqs/galaxy/datasets_change_datatype.html) a datatype can lead to errors. Be intentional and consider [converting](https://training.galaxyproject.org/training-material/faqs/galaxy/datasets_convert_datatype.html) instead when possible. 
- Dataset content can also be adjusted (tools: Data manipulation) and the expected [datatype detected](https://training.galaxyproject.org/training-material/faqs/galaxy/#detecting-the-datatype-file-format). Detected datatypes are the most reliable in most cases.
- If a tool does not accept a dataset as valid input, it is not in the correct format with the correct datatype.
- Once a dataset’s content matches the datatype, and that dataset is repeatedly used (example: Reference annotation) use that same dataset for all steps in an analysis or expect problems. This may mean rerunning prior tools if you need to make a correction. 
- Tip: Not sure what datatypes a tool is expecting for an input?
    1. Create a new empty history
    2. Click on a tool from the tool panel
    3. The tool form will list the accepted datatypes per input
- _Warning_: In some cases, tools will transform a dataset to a new datatype at runtime for you.
    - This is generally helpful, and best reserved for smaller datasets.
    - Why? This can also unexpectedly create hidden datasets that are near duplicates of your original data, only in a different format. 
    - For large data, that can quickly consume working space (quota). 
    - Deleting/purging any hidden datasets can lead to errors if you are still using the original datasets as an input. 
    - Consider [converting](https://training.galaxyproject.org/training-material/faqs/galaxy/datasets_convert_datatype.html) to the expected datatype yourself when data is large.
    - Then test the tool directly on converted data. If it works, purge the original to recover space.