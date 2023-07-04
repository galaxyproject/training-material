---
layout: tutorial_hands_on

title: 'Mass spectrometry: GC-MS analysis with apLCMS'
zenodo_link: ''
questions:
- Which biological questions are addressed by the tutorial?
- Which bioinformatics techniques are important to know for this type of data?
objectives:
- The learning objectives are the goals of the tutorial
- They will be informed by your audience and will communicate to them and to yourself
  what you should focus on during the course
- They are single sentences describing what a learner should be able to do once they
  have completed the tutorial
- You can use Bloom's Taxonomy to write effective learning objectives
time_estimation: 3H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- contributor1
- contributor2

---


# Introduction
{:.no_toc}

`recetox-aplcms` is a software package designed for the processing of LC/MS based metabolomics data, in particular for peak detection in high resolution mass spectrometry (HRMS) data. 
It supports reading `.mzml` files in raw profile mode and uses a bi-Gaussian chromatographic peak shape ({% cite yu2010quantification %}) for feature detection and quantification. `recetox-aplcms` is based on the apLCMS R package ({% cite 10.1093/bioinformatics/btp291 %}) and includes various software updates and is actively developed and maintained on GitHub.

There are two major routes of data analysis. The first, which we call **unsupervised** analysis, does not use existing knowledge. It detects peaks de novo from the data based on the data itself. The second, which we call **hybrid** analysis, combines de novo peak detection with existing knowledge. The existing knowledge can come from two sources - known metabolites and historically detected features from the same machinery.

![recetox-aplcms overview](../../images/aplcms_scheme.png "Overview of individual steps of the recetox-aplcms package that can be combined in two separate workflows processing HRMS data in an unsupervised manner or by including a-priori knowledge.")


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Title for your first section

Give some background about what the trainees will be doing in the section.
Remember that many people reading your materials will likely be novices,
so make sure to explain all the relevant concepts.

## Title for a subsection
Section and subsection titles will be displayed in the tutorial index on the left side of
the page, so try to make them informative and concise!

# Hands-on Sections
Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

# Common part

It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:


## Sub-step with **recetox-aplcms - remove noise**

intensity is droping over time in the experiment, we want to normalise this and remove noise

Parameters
- `min_pres` - in how many consequent scans do we want to have the signal present
  - e.g. data point show up only every 3 scans for 10 seconds, this setting can filter it out
  - we can't set the maximum at this point
- `min_run` - length of elution time (seconds) - what has to be the minimum width of a peak to it be considered a peak
- `mz_tol` - ppm tolerance, helps us what is still considered as one signal, width of m/z peak
- `baseline` - intensity cutoff you want to use

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [recetox-aplcms - remove noise](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_remove_noise/recetox_aplcms_remove_noise/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input spectra data"*: `output` (Input dataset collection)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Can we provide some numerical intervals for input parameters?
> 2. Is the data filtered also on retention time?
>
> > ### {% icon solution %} Solution
> >
> > 1. Not really, because they are very instrument and data specific.
> > 2. Yes, and it can be adjusted with the `min_run` parameter.
> >
> {: .solution}
>
{: .question}

Output is in the `.parquet` format, just a binary (more accurate yet smaller in size) representation of tabular format. We would loose accuracy by storing it in text format. In the table we have our filtered data, already preliminary grouped on m/z basis - within the ppm tolerance we chose. 

## Sub-step with **recetox-aplcms - generate feature table**

Peaks detection
- we want to fit peak shapes to our data, this allows us to compute precise intensities
- this also resolves peak overlaps, by fitting them both as separate peaks
- this then ofc does work with centroid data (there we dont have the peak shapes anymore, just some "averages" of them)

This tool takes the grouped features from the previous step and computes the peak shape in retention time domain and then integrates the peak area


Parameters:
- `sd_cut_min/max` - here we can now specify the maximum and minimum peak width by selecting allowed range for standard deviation (either side)
- `sigma_ratio_lim_min/max` - ratio (sd1 divided by sd2) between standard deviations
  **TODO** can we add a picture explainig its effects?
- `bandwidth` - to improve the peak shape by smoothing

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [recetox-aplcms - generate feature table](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_generate_feature_table/recetox_aplcms_generate_feature_table/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input profile data"*: `output_file` (output of **recetox-aplcms - remove noise** {% icon tool %})
>    - *"shape_model"*: `bi-Gaussian`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Why do we need to fit our peaks to a (in this case Gaussian) shape?
> 2. Why do we get two standard deviations?
>
> > ### {% icon solution %} Solution
> >
> > 1. This allows us to precisely compute the area under the curve and get precise intensity.
> > 2. two standard deviations (sd1, sd2) because we use bi-Gausian shape (the sides of the shape can be different). if they would be same, then its Gaussian
> >
> {: .solution}
>
{: .question}

## Sub-step with **recetox-aplcms - compute clusters**

pre-alignment step - we put all peaks from all samples to one table and group them based on m/z and rt... for that we can specify parameters which parametrise the "size" of buckets (clusters)
- all peaks get assignd a cluster_id shared across samples

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [recetox-aplcms - compute clusters](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_compute_clusters/recetox_aplcms_compute_clusters/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input data"*: `output_file` (output of **recetox-aplcms - generate feature table** {% icon tool %})
>    - *"Tolerances input method"*: `direct`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Can we influence the clustering sensitivity w.r.t. retention time?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Yes, use `mz_tol_absolute` parameter.
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

Output is again separated to individual tables, but with assigned cluster ID.

tolerances - these are values entered in the clustering step - if not provided, they are estimated instead

## Sub-step with **recetox-aplcms - compute template**

We need a template to which we align the data - this takes peak time which has the highest number of features... gives us the highest number of reference points we fit our curve to. This step can be potentially skipped if you want to select a template manually or provide a custom file used in further steps.

**TODO** add hint how to select single file from collection as input for a tool

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [recetox-aplcms - compute template](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_compute_template/recetox_aplcms_compute_template/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input data"*: `clustered_feature_tables` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **recetox-aplcms - correct time**

We need our previously clustered features, selected template and tolerances.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [recetox-aplcms - correct time](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_correct_time/recetox_aplcms_correct_time/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input clustered features table"*: `clustered_feature_tables` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>    - {% icon param-file %} *"Input template features table"*: `output_file` (output of **recetox-aplcms - compute template** {% icon tool %})
>    - {% icon param-file %} *"Input tolerances values"*: `tolerances` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **recetox-aplcms - compute clusters**

second round of grouping... comment on how you can iterativelly combine this steps multiple times, since you still get the same format of outputs etc.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [recetox-aplcms - compute clusters](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_compute_clusters/recetox_aplcms_compute_clusters/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input data"*: `output_file` (output of **recetox-aplcms - correct time** {% icon tool %})
>    - *"Tolerances input method"*: `file`
>        - {% icon param-file %} *"Input tolerances values"*: `tolerances` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **recetox-aplcms - align features**

kernel density alignment

- `min_occurrence` - how many samples does the feature need to be present in... This way we can filter peaks that appear really consistenly across the samples.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [recetox-aplcms - align features](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_align_features/recetox_aplcms_align_features/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Clustered features"*: `clustered_feature_tables` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>    - {% icon param-file %} *"Input tolerances values"*: `tolerances` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

Outputs:
- rt_table - rt for all samples for our features
- intensity_table - intensities for all samples
- peak_metadata - all data related to specific peaks that have been detected (enumerate columns), number of peaks that have been grouped together... 

# Unsupervised

**TODO** somehow inspect tables, comment on many gaps

## Sub-step with **recetox-aplcms - recover weaker signals**

Our tables have many gaps, some features weren't detected in some samples... but it doesn't mean they actually aren't there... so we revisit the data trying to recover them, we do another round of peak picking without any noise filtering, but only on specific place (specified by m/z and rt ranges from already analysed data)

- `recover_min_count` - how many data points need to be there to consider them a peak
- `bandwidth` - again the same for smoothing

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [recetox-aplcms - recover weaker signals](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_recover_weaker_signals/recetox_aplcms_recover_weaker_signals/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input spectra data"*: `output` (Input dataset collection)
>    - {% icon param-file %} *"Input extracted feature samples collection"*: `output_file` (output of **recetox-aplcms - generate feature table** {% icon tool %})
>    - {% icon param-file %} *"Input corrected feature samples collection"*: `output_file` (output of **recetox-aplcms - correct time** {% icon tool %})
>    - {% icon param-file %} *"Metadata table"*: `metadata_file` (output of **recetox-aplcms - align features** {% icon tool %})
>    - {% icon param-file %} *"RT table"*: `rt_file` (output of **recetox-aplcms - align features** {% icon tool %})
>    - {% icon param-file %} *"Intensity table"*: `intensity_file` (output of **recetox-aplcms - align features** {% icon tool %})
>    - {% icon param-file %} *"Input tolerances values"*: `tolerances` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **recetox-aplcms - compute clusters**

We might have added new features, so we do the clustering again.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [recetox-aplcms - compute clusters](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_compute_clusters/recetox_aplcms_compute_clusters/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input data"*: `output_file` (output of **recetox-aplcms - recover weaker signals** {% icon tool %})
>    - *"Tolerances input method"*: `direct`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **recetox-aplcms - align features**

Features can now appear in more samples, so we also repeat the alignment step.

These steps can be potentially again repeated and combined (as well as for example combined with retention time correction) in arbitrary number of iterations.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [recetox-aplcms - align features](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_align_features/recetox_aplcms_align_features/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Clustered features"*: `clustered_feature_tables` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>    - {% icon param-file %} *"Input tolerances values"*: `tolerances` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

# Hybrid

## Sub-step with **recetox-aplcms - merge known table**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [recetox-aplcms - merge known table](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_merge_known_table/recetox_aplcms_merge_known_table/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Metadata table"*: `metadata_file` (output of **recetox-aplcms - align features** {% icon tool %})
>    - {% icon param-file %} *"RT table"*: `rt_file` (output of **recetox-aplcms - align features** {% icon tool %})
>    - {% icon param-file %} *"Intensity table"*: `intensity_file` (output of **recetox-aplcms - align features** {% icon tool %})
>    - {% icon param-file %} *"Table of known features"*: `output` (Input dataset)
>    - {% icon param-file %} *"Input tolerances values"*: `tolerances` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>    - *"Tables merge direction"*: `Merge known table to features`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **recetox-aplcms - recover weaker signals**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [recetox-aplcms - recover weaker signals](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_recover_weaker_signals/recetox_aplcms_recover_weaker_signals/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input spectra data"*: `output` (Input dataset collection)
>    - {% icon param-file %} *"Input extracted feature samples collection"*: `output_file` (output of **recetox-aplcms - generate feature table** {% icon tool %})
>    - {% icon param-file %} *"Input corrected feature samples collection"*: `output_file` (output of **recetox-aplcms - correct time** {% icon tool %})
>    - {% icon param-file %} *"Metadata table"*: `output_metadata_file` (output of **recetox-aplcms - merge known table** {% icon tool %})
>    - {% icon param-file %} *"RT table"*: `output_rt_file` (output of **recetox-aplcms - merge known table** {% icon tool %})
>    - {% icon param-file %} *"Intensity table"*: `output_intensity_file` (output of **recetox-aplcms - merge known table** {% icon tool %})
>    - {% icon param-file %} *"Input tolerances values"*: `tolerances` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **recetox-aplcms - compute clusters**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [recetox-aplcms - compute clusters](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_compute_clusters/recetox_aplcms_compute_clusters/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input data"*: `output_file` (output of **recetox-aplcms - recover weaker signals** {% icon tool %})
>    - *"Tolerances input method"*: `file`
>        - {% icon param-file %} *"Input tolerances values"*: `tolerances` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **recetox-aplcms - compute template**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [recetox-aplcms - compute template](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_compute_template/recetox_aplcms_compute_template/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input data"*: `clustered_feature_tables` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **recetox-aplcms - correct time**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [recetox-aplcms - correct time](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_correct_time/recetox_aplcms_correct_time/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input clustered features table"*: `clustered_feature_tables` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>    - {% icon param-file %} *"Input template features table"*: `output_file` (output of **recetox-aplcms - compute template** {% icon tool %})
>    - {% icon param-file %} *"Input tolerances values"*: `tolerances` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **recetox-aplcms - compute clusters**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [recetox-aplcms - compute clusters](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_compute_clusters/recetox_aplcms_compute_clusters/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input data"*: `output_file` (output of **recetox-aplcms - correct time** {% icon tool %})
>    - *"Tolerances input method"*: `file`
>        - {% icon param-file %} *"Input tolerances values"*: `tolerances` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **recetox-aplcms - align features**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [recetox-aplcms - align features](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_align_features/recetox_aplcms_align_features/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Clustered features"*: `clustered_feature_tables` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>    - {% icon param-file %} *"Input tolerances values"*: `tolerances` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **recetox-aplcms - merge known table**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [recetox-aplcms - merge known table](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_merge_known_table/recetox_aplcms_merge_known_table/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Metadata table"*: `metadata_file` (output of **recetox-aplcms - align features** {% icon tool %})
>    - {% icon param-file %} *"RT table"*: `rt_file` (output of **recetox-aplcms - align features** {% icon tool %})
>    - {% icon param-file %} *"Intensity table"*: `intensity_file` (output of **recetox-aplcms - align features** {% icon tool %})
>    - {% icon param-file %} *"Table of known features"*: `output` (Input dataset)
>    - {% icon param-file %} *"Input tolerances values"*: `tolerances` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>    - *"Tables merge direction"*: `Merge features to known table`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.