---
layout: tutorial_hands_on

title: "Combining single cell datasets after pre-processing"
subtopic: single-cell-CS
priority: 2

zenodo_link: 'https://zenodo.org/records/15090813'

questions:
  - I have some AnnData files from different samples that I want to combine into a single file. How can I combine these and label them within the object?

redirect_from:
  - /topics/transcriptomics/tutorials/scrna-case_alevin-combine-datasets/tutorial

answer_histories:
  - label: "UseGalaxy.eu"
    history: https://singlecell.usegalaxy.eu/u/wendi.bacon.training/h/locked-combining-single-cell-datasets-after-pre-processing-answer-key
    date: 2025-03-20
  - label: "UseGalaxy.eu-ARCHIVE2"
    history: https://usegalaxy.eu/u/j.jakiela/h/combining-datasets-key-history
    date: 2024-03-26
  - label: "UseGalaxy.eu-ARCHIVE1"
    history: https://usegalaxy.eu/u/wendi.bacon.training/h/combining-datasets-400k-key-history
    date: 2024-12-10
  - label: "All total samples - processed after Alevin into single object (UseGalaxy.eu)"
    history: https://usegalaxy.eu/u/j.jakiela/h/all-total-samples-processed-after-alevin-into-single-object
    date: 2024-03-26

input_histories:
  - label: "UseGalaxy.eu"
    history: https://singlecell.usegalaxy.eu/u/wendi.bacon.training/h/locked-combining-single-cell-datasets-after-pre-processing-input
    date: 2025-03-20
  - label: "UseGalaxy.eu-ARCHIVE2"
    history: https://usegalaxy.eu/u/j.jakiela/h/combining-datasets-input
  - label: "UseGalaxy.eu-ARCHIVE1"
    history: https://usegalaxy.eu/u/wendi.bacon.training/h/combining-datasets-input
    date: 2024-12-10
  - label: "UseGalaxy.org"
    history: https://usegalaxy.org/u/juliajot/h/combining-datasets-input

extra:
  example_histories:
    - label: "UseGalaxy.eu_Whole_Separate-AnnData"
      history: https://singlecell.usegalaxy.eu/u/wendi.bacon.training/h/locked-combining-single-cell-datasets-after-pre-processing-input-history-whole-samples
      date: 2025-03-26
    - label: "UseGalaxy.eu_Whole_Combined-AnnData"
      history: https://singlecell.usegalaxy.eu/u/wendi.bacon.training/h/locked-whole-samples-combining-single-cell-datasets-after-pre-processing-answer-key
      date: 2025-03-26

objectives:
  - Combine data matrices from different samples in the same experiment
  - Label the metadata for downstream processing

time_estimation: 1H

key_points:
  - Create a single scanpy-accessible AnnData object from multiple AnnData files, including relevant cell metadata according to the study design
  - Retreive partially analysed data from a public repository

tags:
  - paper-replication
  - MIGHTS

contributions:
  authorship:
    - nomadscientist


  editing:
    - hexylena
    - pinin4fjords

  testing:
    - wee-snufkin

requirements:
  - type: "internal"
    topic_name: single-cell
    tutorials:
        - scrna-intro
        - scrna-umis
        - scrna-case_alevin


gitter: Galaxy-Training-Network/galaxy-single-cell

follow_up_training:
  -
    type: "internal"
    topic_name: single-cell
    tutorials:
        - scrna-case_basic-pipeline

recordings:
- youtube_id: 22t-4qvHnow
  length: 11M
  date: '2023-05-09'
  speakers:
  - hrukkudyr
  captioners:
  - hrukkudyr
- captioners:
  - nomadscientist
  date: '2021-02-15'
  galaxy_version: '21.09'
  length: 11M
  youtube_id: U8pVa6csmUE
  speakers:
  - nomadscientist

---


This tutorial will take you from the multiple AnnData outputs of the [previous tutorial](https://humancellatlas.usegalaxy.eu/training-material/topics/transcriptomics/tutorials/scrna-case_alevin/tutorial.html) to a single, combined  AnnData object, ready for all the fun downstream processing. We will also look at how to add in metadata (for instance, SEX or GENOTYPE) for analysis later on.

> <details-title>Where am I?</details-title>
>
> You are in one of the four tutorials associated with a Case Study, which replicates and expands on the analysis performed in a manuscript {% cite Bacon2018 %}.
>
> ![case study overview](../../images/scrna-casestudy/CS2-case_study_overview_2.png "Overview of the four parts of the case study, represented by boxes connected with noodles. There is a signpost specifying that you are currently in the second part.")
>
{: .details}

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


## Get Data
The sample data is a subset of the reads from each of the seven samples in a mouse dataset of fetal growth restriction {% cite Bacon2018 %} (see the [study in Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-6945/results/tsne) and the [project submission](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6945/)).

> <details-title>Fun facts about the dataset</details-title>
>
> We are particularly keen for learners to be able to go from raw FASTQ files all the way through analysis. We aren't handing you a curated dataset that we specially modified in order for this tutorial to work.
> - This tutorial's input dataset is the full dataset generated from the {% icon level %} [previous tutorial in this case study]({% link topics/single-cell/tutorials/scrna-case_alevin/tutorial.md %}).
> - The only difference is that in that previous tutorial, we only analysed one dataset. Our input data, however, is the result of applying that workflow to all seven datasets.
> - You can find this dataset in the [study in Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-6945/results/tsne) and the [project submission](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6945/).
>
{: .details}

{% include _includes/cyoa-choices.html option1="Import History" option2="Zenodo" default="Import-History"
       text="If you're on the EU or ORG/USA server, then the quickest way to Get the Data for this tutorial is via importing a history. Otherwise, you can also import from Zenodo - it just might take a moment longer if you're in a live course and everyone is importing the same dataset at the same time!" %}
<div class="Import-History" markdown="1">

> <hands-on-title>Import History</hands-on-title>
>
> 1. Import the {% icon galaxy-history-input %} *Input history* by following the link below
>
>    {% for h in page.input_histories %}
>      [ {{h.label}} Input History]( {{h.history}} )
>    {% endfor %}
>
>    {% snippet faqs/galaxy/histories_import.md %}
>
> If you want to import the history to another Galaxy server, check how to do it below!
>
> {% snippet faqs/galaxy/histories_transfer_entire_histories_from_one_galaxy_server_to_another.md %}
>
{: .hands_on}

</div>

<div class="Zenodo" markdown="1">

> <hands-on-title>Data upload for 7 files</hands-on-title>
>
> 1. Create a new history for this tutorial (if you're not importing the history above)
> 2. Import the different AnnData files and the experimental design table from [Zenodo](https://zenodo.org/records/10852529).
>
>    ```
>    {{ page.zenodo_link }}/files/N701-400k-AnnData.h5ad
>    {{ page.zenodo_link }}/files/N702-400k-AnnData.h5ad
>    {{ page.zenodo_link }}/files/N703-400k-AnnData.h5ad
>    {{ page.zenodo_link }}/files/N704-400k-AnnData.h5ad
>    {{ page.zenodo_link }}/files/N705-400k-AnnData.h5ad
>    {{ page.zenodo_link }}/files/N706-400k-AnnData.h5ad
>    {{ page.zenodo_link }}/files/N707-400k-AnnData.h5ad
>    {{ page.zenodo_link }}/files/Experimental_Design.tabular
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Rename the datasets
> 4. Check that the AnnData files datatype is `h5ad`, otherwise you will need to change each file to `h5ad`!
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}

</div>

{% icon galaxy-eye %} Inspect the {% icon param-file %} `Experimental Design` text file. This shows you how each `N70X` corresponds to a sample, and whether that sample was from a male or female. This will be important metadata to add to our sample, which we will add very similarly to how you added the `gene_name` and `mito` metadata previously!

# Important tips for easier analysis

{% snippet faqs/galaxy/tutorial_mode.md %}

{% snippet topics/single-cell/faqs/single_cell_omics.md %}

{% snippet faqs/galaxy/analysis_troubleshooting.md sc=true %}

{% snippet faqs/gtn/gtn_example_histories.md %}

# Concatenating objects

> <hands-on-title>Concatenating AnnData objects</hands-on-title>
>
> 1. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.10.9+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `N701-400k`
>    - *"Function to manipulate the object"*: `Concatenate along the observations axis`
>    - {% icon param-file %} *"Annotated data matrix to add"*: `Select all the other matrix files from bottom to top, N702 to N707`
>
>    > <comment-title></comment-title>
>    > If you imported files from Zenodo instead of using the input history, yours might not be in the same order as ours. Since the files will be concatenated in the order that you click, it will be helpful if you click them in the same order, from N702 to N707. This will ensure your samples are given the same batch numbers as we got in this tutorial, which will help when we're adding in metadata later!
>    {: .comment}
>
>    > <warning-title>Don't add N701!</warning-title>
>    > You are adding files to N701, so do not add N701 to itself!
>    {: .warning}
>
>    - *"Join method"*: `Intersection of variables`
>    - *"Key to add the batch annotation to obs"*: `batch`
>    - *"Separator to join the existing index names with the batch category"*: `-`
> 2. Rename {% icon galaxy-pencil %} output `Combined_Object`
{: .hands_on}

Now let's look at what we've done! You can *peek* at the AnnData object in your {% icon galaxy-history %} history. You can also inspect the dataset for further details using a tool.

> <hands-on-title>Inspecting AnnData Objects</hands-on-title>
>
> 1. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.10.9+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Combined_object`
>    - *"What to inspect?"*: `General information about the object`
> 2. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.10.9+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Combined_object`
>    - *"What to inspect?"*: `Key-indexed observations annotation (obs)`
> 3. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.10.9+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Combined_object`
>    - *"What to inspect?"*: `Key-indexed annotation of variables/features (var)`
{: .hands_on}

Now have a look at the three {% icon tool %} **Inspect AnnData** outputs.

> <question-title></question-title>
>
> 1. How many cells do you have now?
> 2. Where is `batch` information stored?
>
> > <solution-title></solution-title>
> >
> > 1. If you peek at your dataset, or look at the **General information** {% icon tool %} output, you will find there are `316 cells`, as the matrix is now 316 cells (n_obs) x 35734 genes (n_var). You will also find these numbers in the **obs** {% icon tool %} (cells) and **var** {% icon tool %} (genes) file sizes.
> > 2. Batch information is stored under **Key-indexed observations annotation (obs)**. Different versions of the Manipulate tool might put the `batch` columns in different locations. The tool version in this tutorial puts `batch` in the `8th` column. Batch refers to the order in which the matrices were added. The files are added from the bottom of the history upwards, so be careful how you set up your histories when running this (i.e. if your first dataset is N703 and the second is N701, the `batch` will call N703 `0` and N701 `1`!)
> {: .solution}
>
{: .question}

# Adding cell metadata

I set up the example history with the earliest indices at the bottom.

![The files are numbered such that dataset #1 is N701, dataset #2 is N702, etc., up through #7 as N707. This puts N707 at the top and N701 at the bottom in the Galaxy history.](../../images/scrna-casestudy/CS2-history-files-ascending.png "Correct history ordering for combining datasets in order"){ width=400px }

Therefore, when it is all concatenated together, the `batch` appears as follows:

| Index | Batch | Genotype | Sex |
|------ |--------------------|
| N701 | 0    | wildtype    | male    |
| N702 | 1    | knockout   | male    |
| N703 | 2    | knockout   | female    |
| N704 | 3    | wildtype    | male    |
| N705 | 4    | wildtype    | male    |
| N706 | 5    | wildtype    | male    |
| N707 | 6    | knockout    | male    |

If you used Zenodo to import files, they may not have imported in order (i.e. N701 to N707, ascending). In that case, you will need to tweak the parameters of the next tools appropriately to label your batches correctly!

The two critical pieces of metadata in this experiment are **sex** and **genotype**. I will later want to color my cell plots by these parameters, so I want to add them in now!

## Sex & Genotype metadata

> <hands-on-title>Labelling sex</hands-on-title>
>
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_column/9.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: output of **Inspect AnnData: Key-indexed observations annotation (obs)** {% icon tool %})
>    - *"1. Replacement"*
>
>         - *"in column"*: `Column: 8` - or whichever column `batch` is in
>         - *"Find pattern"*: `0|1|3|4|5|6`
>         - *"Replace with"*: `male`
>    - **+  Insert Replacement**
>    - *"2. Replacement"*
>
>         - *"in column"*: `Column: 8`
>         - *"Find pattern"*: `2`
>         - *"Replace with"*: `female`
>    - **+  Insert Replacement**
>    - *"3. Replacement"*
>
>         - *"in column"*: `Column: 8`
>         - *"Find pattern"*: `batch`
>         - *"Replace with"*: `sex`
>
>    The output of the {% icon param-file %} **Replace Text** {% icon tool %} tool has many columns. However, we want only the column containing the **sex** information, in order by cell barcode, to add to our AnnData later.
>
> 2. {% tool [Cut columns from a table](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c8`
>    - *"Delimited by"*: `Tab`
>    - {% icon param-file %} *"From"*: output of **Replace text** {% icon tool %}
>
> 3. Rename {% icon galaxy-pencil %} output `Sex_metadata`
{: .hands_on}

That was so fun, let's do it all again but for genotype!

> <hands-on-title>Labelling genotype</hands-on-title>
>
> 1. {% tool [Replace Text in a specific column](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_column/9.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: output of **Inspect AnnData: Key-indexed observations annotation (obs)** {% icon tool %}
>    - *"1. Replacement"*
>
>         - *"in column"*: `Column: 8`
>         - *"Find pattern"*: `0|3|4|5`
>         - *"Replace with"*: `wildtype`
>    - **+ Insert Replacement**
>    - *"2. Replacement"*
>
>         - *"in column"*: `Column: 8`
>         - *"Find pattern"*: `1|2|6`
>         - *"Replace with"*: `knockout`
>    - **+ Insert Replacement**
>    - *"3. Replacement"*
>
>         - *"in column"*: `Column: 8`
>         - *"Find pattern"*: `batch`
>         - *"Replace with"*: `genotype`
>
>    Now we want only the column containing the genotype information - we will ultimately add this into the cell annotation in the AnnData object.
>
> 2. {% tool [Cut columns from a table](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c8`
>    - *"Delimited by"*: `Tab`
>    - {% icon param-file %} *"From"*: output of **Replace text** {% icon tool %}
>
> 3. Rename {% icon galaxy-pencil %} output `Genotype_metadata`
{: .hands_on}

You might want to do this with all sorts of different metadata - which labs handled the samples, which days they were run, etc. Once you've created and cut all your metadata columns, we can paste them together before adding them into the AnnData object itself.

> <hands-on-title>Combining metadata columns</hands-on-title>
>
> 1. {% tool [Paste two files side by side](Paste1) %} with the following parameters:
>    - {% icon param-file %} *"Paste"*: `Genotype_metadata`
>    - {% icon param-file %} *"and"*: `Sex_metadata`
>    - *"Delimit by"*: `Tab`
> 2. Rename {% icon galaxy-pencil %} output `Cell_Metadata`
{: .hands_on}

Let's add this metadata to the AnnData object!

> <hands-on-title>Adding metadata to AnnData object</hands-on-title>
>
> 1. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.10.9+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Combined_object`
>    - *"Function to manipulate the object"*: `Add new annotation(s) for observations or variables`
>    - *"What to annotate?"*: `Observations (obs)``
>    - {% icon param-file %} *"Table with new annotations"*: `Cell_Metadata`
{: .hands_on}

Woohoo! We're there! You can *peek* at your new AnnData object in your history {% icon galaxy-history %} to see the additional **Obs** categories to make sure this worked. You should now find a `sex` and `genotype` in the **Obs** listing.

## Add batch metadata

I want to clean up this AnnData object just a bit more first. It would be a lot nicer if 'batch' meant something, rather than 'the order in which the Manipulate AnnData tool added my datasets'.

> <hands-on-title>Labelling batches</hands-on-title>
>
> 1. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.10.9+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: output of **Manipulate AnnData - Add new annotations** {% icon tool %}
>    - *"Function to manipulate the object"*: `Rename categories of annotation`
>    - *"Key for observations or variables annotation"*: `batch`
>    - *"Comma-separated list of new categories"*: `N701,N702,N703,N704,N705,N706,N707`
>
> 2. Rename {% icon galaxy-pencil %} output `Batched_Object`
>
{: .hands_on}

{% icon congratulations %} Well done! *Peek* at your final {% icon param-file %} `Batched_Object` to see the wealth of information that has been added. You are now ready to move along to further filtering!

# Conclusion

You might find the {% icon galaxy-history-answer %} *Answer Key Histories* helpful to check or compare with:
  - {% for h in page.answer_histories %}
      [ {{h.label}} ]( {{h.history}} )
    {% endfor %}
<!-- Only currently want to iterate through the first history, but might want others in the future (different servers!) -->

You can also run this entire tutorial via a {% icon galaxy-workflows-activity %} *Workflow*, after performing the **Get data** step initially.
 - [Tutorial Workflow]({% link topics/single-cell/tutorials/scrna-case_alevin-combine-datasets/workflows/index.md %})

<iframe title="Galaxy Workflow Embed" style="width: 100%; height: 700px; border: none;" src="https://singlecell.usegalaxy.eu/published/workflow?id=d1a1e6070c0a36ca&embed=true&buttons=true&about=false&heading=false&minimap=true&zoom_controls=true&initialX=0&initialY=-20&zoom=0.5"></iframe>

Finally, you may remember that the datasets you analysed in this tutorial were downsampled. The idential analysis can be performed on whole samples, which you can find in the example histories below.
- {% icon warning %} Remember that *if* you are in a course, time for exploring these histories will not be factored into the schedule. Explore these outside of course time!
- {% for h in page.example_histories %}
    [ {{h.label}} ]( {{h.history}} )
  {% endfor %}

{% snippet topics/single-cell/faqs/user_community_join.md %}
