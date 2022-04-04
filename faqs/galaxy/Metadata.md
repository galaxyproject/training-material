---
title: Metadata
area: Getting inputs right    # FAQs will be grouped by these areas on the FAQ page
description:  How do I find, adjust, and/or correct metadata?
layout: faq         

---

### When and how are metadata assigned?
- During **Upload** on new datasets, either autodetected (default) or directly assigned by the user
- With the **Edit Attributes functions** on existing datasets, either autodetected or directly assigned by the user
- By **upstream tools** based on the those tool's input metadata assignments (inherited) or as new defined outputs types

### Changing or updating dataset metadata
Click on the the dataset's pencil icon {% icon galaxy-pencil %} to reach the _Edit Attributes_ forms in the center panel to change or update metadata. Metadata attributes can include: database (dbkey, genome, transcriptome, custom genome), datatype (including compression format, e.g. gzipped versus not), column assignments (chrom, start, end, and more), row assignments (specify comments lines), permissions, and others.

### Tool errors related to metadata
Many tools can error when the wrong dataset attributes (metadata) are assigned. Some assignments are automatic without user action (example: tool outputs), others are or can be "autodetected" if incorrect or need modification (example: most datatypes), and others required action by the user if the attribute cannot be autodetected (example: assigning database to newly uploaded data).

### Example: Dataset metadata missing or incomplete

#### How to notice the problem
- The dataset will not download when using the disk icon
- Tools error when using a specific dataset that has been used before successfully
- Tools error with a message that ends with: ``OSError: [Errno 2] No such file or directory``.

### Common solutions
**Reset the dataset(s) metadata.**

This could be a newly uploaded dataset or one created by prior tools. If a job failed, it is one or more of the input datasets to the tool that errored. Click on the the dataset's pencil icon to reach the Edit Attributes forms in the center panel. Then do one of the following as applies:
- **_Directly reset metadata_** Find the tab for the metadata you want to change, make the change, and **save**.
-  **_Autodetect metadata_** Click on the **_Auto-detect_** button found on the first **_Edit Attributes_** tab. The dataset will turn yellow in the history while the job is processing.
The dataset will turn yellow in the history while the job is processing. Once compeleted, other attributes can be changed and/or the dataset can be used.

#### If problems persist
Your error may not really be related to metadata, or changing metadat may not be enough to resolve the input/tool usage problems. Review the other [FAQ solutions](https://galaxyproject.org/support/#troubleshooting) for more data input help. Review the [Galaxy Tutorials](https://galaxyproject.org/learn/) for tool usage help.


