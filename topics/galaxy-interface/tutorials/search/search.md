# Basic Searching

You can filter what datasets are shown and search for datasets using the search bar at the top of the panel. Enter
any text that a dataset you'd be looking for would contain, including:

* the name or part of the name
* any text (or partial text) from the info field
* the file format or reference database
* any text or partial text from the annotation or tags of a dataset

For example:

* To find all vcf files you might enter: `vcf` alone.
* To find all files whose names contain data 1, you can enter: `data 1`
* To search for a VCF file named 'VCF filter on data 1' and tagged with 'experiment_1', you could enter: `vcf filter on data 1 experiment_1`

## Clearing a Search

You can clear a search and show all visible datasets by clicking the round 'X' button in the right of the search bar
or - while entering text in the search bar - hitting the escape key ('Esc').

# Advanced Searching

You can also specify dataset properties that you want to filter on. If you search with multiple properties, these are connected with ANDs, so datasets must match all provided attributes.

Query                                              | Results
-----                                              | ------
`name:'FASTQC on'`                                 | Any datasets with "FASTQC on" in the title, but avoids items which have "FASTQC on" in other fields like the description or annotation.
`extension:vcf`                                       | Datasets with a specific format. Some formats are hierarchical, e.g. searching for `fastq` will find fastq files but also fastqsanger and fastqillumina files. You can see more formats in the upload dialogue.
`tag:experiment1 tag:to_publish`                   | for searching on (a partial) dataset tag. You can repeat to search for more tags.
`related:10`                                       | A specific history item ID (based on the ordering in the history)
`state:error`                                      | To show only datasets in a given state. Other options include `ok`, `running`, `paused`, and `new`.

If you find normal searching is showing too many datasets, and not what you're looking for, try the advanced search! Just use the {% icon galaxy-advanced-search %} button next to the search field to show the advanced selector.
