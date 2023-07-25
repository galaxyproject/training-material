---
layout: tutorial_hands_on

title: 'Mass spectrometry: GC-MS analysis with apLCMS'
zenodo_link: 'https://zenodo.org/record/7890956'
level: Intermediate
questions:
- TBA
objectives:
- TBA
time_estimation: 2H
key_points:
- The processing of untargeted GC-MS metabolomic data can be done using open-source tools.
- This processing workflow is complementary to XCMS.
contributions:
  authorship:
    - xtrojak
    - hechth

---


# Introduction
{:.no_toc}

`recetox-aplcms` is a software package designed for the processing of LC/MS based metabolomics data, in particular for peak detection in high resolution mass spectrometry (HRMS) data. 
It supports reading `.mzml` files in raw profile mode and uses a bi-Gaussian chromatographic peak shape ({% cite yu2010quantification %}) for feature detection and quantification. `recetox-aplcms` is based on the apLCMS R package ({% cite 10.1093/bioinformatics/btp291 %}) and includes various software updates and is actively developed and maintained on GitHub.

There are two major routes of data analysis. The first, which we call **unsupervised** analysis, does not use existing knowledge. It detects peaks de novo from the data based on the data itself. The second, which we call **hybrid** analysis, combines de novo peak detection with existing knowledge. The existing knowledge can come from two sources - known metabolites and historically detected features from the same machinery. While the unsupervised approach allows for unbiased exploration and discovery of patterns, the hybrid approach integrates prior knowledge or supervised techniques to enhance targeted analysis and interpretation. The choice between these approaches depends on the research objectives, available prior knowledge, and the specific questions being addressed in the metabolomics study.

![recetox-aplcms overview](../../images/aplcms_scheme.png "Overview of individual steps of the recetox-aplcms package that can be combined in two separate workflows processing HRMS data in an unsupervised manner or by including a-priori knowledge.")

The workflows consist of the following building blocks:

- **remove noise** - denoise the raw data and extract the EIC
- **generate feature table** - group features in EIC into peaks using peak-shape model
- **compute clusters** - compute mz and rt clusters across samples
- **compute template** - find the template for rt correction
- **correct time** - correct the rt across samples using splines
- **align features** - align identical features across samples
- **recover weaker signals** - recover missed features in samples based on the aligned features
- **merge known table** - add known features to the detected features table and vice versa

We use three GC-[EI+] high-resolution mass spectrometry data files generated from quality control seminal plasma samples to demonstrate data analysis in this tutorial.

> <details-title> Seminal plasma samples </details-title>
> 
> The seminal plasma samples were analyzed according to the standard operating procedure [(SOP) for metabolite profiling of seminal plasma via GC Orbitrap](https://zenodo.org/record/5734331). The three samples used in this training are pooled quality control (QC) samples coming from about 200 samples. The pooled samples were analyzed in dilution series to test the system suitability and the quality of the assay.
>
{: .details}

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data preparation and prepocessing

Prior to commencing the analysis pipeline, the dataset must be downloaded and prepared. Several preprocessing steps can be executed in parallel on individual samples. To facilitate this process, utilizing Galaxy's Dataset collections is recommended. When uploading data into Galaxy, you can select to create the dataset collection directly.

## Import the data into Galaxy

> <hands-on-title> Upload data </hands-on-title>
>
> 1. Create a new history for this tutorial
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) into a collection:
>
>    ```
>    https://zenodo.org/record/7890956/files/8_qc_no_dil_milliq.raw
>    https://zenodo.org/record/7890956/files/21_qc_no_dil_milliq.raw
>    https://zenodo.org/record/7890956/files/29_qc_no_dil_milliq.raw
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md collection=true format="mzml" collection_name="input" renaming=false %}
>
> 3. Make sure your data is in a **collection**. You can always manually create the collection from separate files:
>
>    {% snippet faqs/galaxy/collections_build_list.md %}
>
>    In the further steps, this dataset collection will be referred to as `input` (and we recommend naming this collection like that to avoid confusion).
>
> 4. Import the following extra file from [Zenodo]({{ page.zenodo_link }}): **TODO**
>
>    ```
>    https://gitlab.ics.muni.cz/umsa/umsa-files/-/raw/master/testdata/recetox-aplcms/hybrid/known_table.parquet
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
> 
>    Please pay attention to the format of the uploaded file, and make sure it is correctly imported.
>
>    {% snippet faqs/galaxy/datatypes_understanding_datatypes.md %}
>
>    > <comment-title> The known table </comment-title>
>    >
>    > **TODO**
>    >
>    {: .comment}
>
{: .hands_on}

As a result of this step, we should have in our history a green Dataset collection with three `.raw` files as well as the table of known metabolites.

## Convert the raw data to mzML

Our input data are in `.raw` format, which is not suitable for the downstream tools in this tutorial. We use the tool {% tool [msconvert](toolshed.g2.bx.psu.edu/repos/galaxyp/msconvert/msconvert/3.0.20287.2) %} to convert our samples to the appropriate format (`.mzML` in this case).

> <hands-on-title> Convert the raw data to mzML </hands-on-title>
>
> 1. {% tool [msconvert](toolshed.g2.bx.psu.edu/repos/galaxyp/msconvert/msconvert/3.0.20287.2) %} with the following parameters:
>    - {% icon param-collection %} *"Input unrefined MS data"*: `input` (Input dataset collection)
>    - *"Do you agree to the vendor licenses?"*: `Yes`
>    - *"Output Type"*: `mzML`
>
{: .hands_on}

# Common part

In this first part, `recetox-aplcms` integrates noise filtering, peak detection and alignment, and statistical analysis to process and extract meaningful information from LC-MS data. To enhance the quality of the data, the tool employs noise filtering techniques to remove false-positive peaks caused by background noise. It applies statistical methods or threshold-based approaches to distinguish true peaks from the noise. 

Then, `recetox-aplcms` detects and extracts individual peaks from the noise-free data. It uses an adaptive algorithm that iteratively identifies peaks by considering the local intensity distributions. Once the peaks are detected, they are grouped based on their chromatographic behaviour across multiple samples. It aligns the peaks by accounting for variations in retention time, which can occur due to instrument drift or other factors. This is achieved by retention time correction, when the tool estimates retention time shifts based on peak intensities and their alignment patterns, iteratively adjusting the retention time values to minimize misalignment and maximize the alignment of peaks with similar chromatographic behaviour. Finally, by considering retention time, m/z values, and peak intensities, `recetox-aplcms` matches corresponding features, ensuring their accurate alignment and enabling meaningful comparisons.

## Remove noise

This step extracts an ion chromatogram (EIC) showing the intensity of only particular ions of interest over time.
Since the intensity is dropping over time in the experiment, we want to normalise this and remove noise from the raw data.
It also performs a first clustering step of points with close m/z values into the extracted ion chromatograms (EICs).

> <details-title> Key parameters </details-title>
> 
> A precise tuning of input parameters is crucial in this step since it can potentially lead to the elimination of some of the data of interest or, on the other extreme, the preservance of some noisy background data:
>
> - **Minimal elution time** - minimal length of elution time of a peak to be actually recognised as a peak. It is closely related to the chromatography method used. Only peaks with greater elution length are kept.
> - **Minimal signal presence** - determines in how many consequent scans do we want to have the signal present. Sometimes the signal from a real feature is not present 100% of the time along the feature's retention time. This parameter sets the threshold for an ion trace to be considered a feature. This parameter is best determined by examining the raw data. For example, if we know a data point shows up only every three scans for 10 seconds, this setting can be used to filter it out.
> - **m/z tolerance** - tolerance (in ppm) to determine how far along m/z do two points need to be separated for them to be considered different a different peak. It can be seen as the width of the m/z peak.
> - **Baseline correction** - intensity cutoff to be used. After grouping the observations, the highest intensity in each group is found. If the highest is lower than this value, the entire group will be deleted.
>
> ![recetox-aplcms noise parameters](../../images/aplcms_explain_min_run_and_pres.png "Graphical explanation of effect of minimal elution time (min_run) and minimal signal presence (min_pres) parameters on input data.")
>
{: .details}

> ### {% icon hands_on %} Hands-on: Remove noise
>
> 1. {% tool [recetox-aplcms - remove noise](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_remove_noise/recetox_aplcms_remove_noise/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input spectra data"*: `mzML input (msconvert)` (output of **msconvert** {% icon tool %})
>    - *"Minimal signal presence [fraction of scans]"*: `0.8`
>    - *"Minimal elution time [unit corresponds to the retention time]"*: `0.5`
>    - *"m/z tolerance [ppm]"*: `30.0`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. Are there any numerical intervals available for the input parameters?
> 2. Can the data be filtered also on the retention time axis?
>
> > ### {% icon solution %} Solution
> >
> > 1. All input parameters are instrument-specific. Therefore we cannot provide recommended intervals for the values.
> > 2. Indeed, for this purpose, the `Minimal elution time` parameter can be used.
> >
> {: .solution}
>
{: .question}

> <details-title> Parquet format </details-title>
> 
> Output is in the `.parquet` format, which is a binary representation of tabular format. 
This format is used to increase the accuracy of stored values, which would be significantly lower when stored in text format.
>
{: .details}

## Generate feature table

This step takes the features grouped by m/z from the previous step and detects peaks. The goal is to fit peak shapes in the retention time domain to our data, which allows computing precise intensities by integrating the peak area. This step also resolves peak overlaps by fitting them both as separate peaks. As a consequence of this approach, `recetox-aplcms` does not work with centroid data since there are no peak shapes anymore, just some "averages" of them.

> <details-title> Key parameters </details-title>
> 
> - **Minimal/maximal standard deviation** - specify the maximum and minimum peak width by selecting the allowed range for the standard deviation (both $$\sigma_1$$ and $$\sigma_2$$).
> - **Minimal/maximal sigma ratio** - the lower and upper limit of the ratio between left-standard deviation and the right-standard deviation $$\frac{\sigma_1}{\sigma_2}$$. It represents the relative skewness of the peak.
> - **Bandwidth factor** - parameter used to scale down the overall range of retention times (the bandwidth) assumed in the kernel smoother used for peak identification. The value is between zero and one. The minimal and maximal bandwidth can be limited by explicit values. It is used to improve the peak shape by smoothing.
>
> ![recetox-aplcms sigma parameters](../../images/aplcms_explain_bi_gaussian.png "The picture shows the bi-Gaussian model, characterised by two standard deviations. The ratio of standard deviations influences the degree of skewness, thus producing a different shape of the peak.")
>
{: .details}

> ### {% icon hands_on %} Hands-on: Generate feature table
>
> 1. {% tool [recetox-aplcms - generate feature table](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_generate_feature_table/recetox_aplcms_generate_feature_table/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input profile data"*: `noise-free data` (output of **recetox-aplcms - remove noise** {% icon tool %})
>    - *"Minimal sigma ratio"*: `0.05`
>    - *"Maximal sigma ratio"*: `20.0`
>    - *"Standard deviations boundaries"*: `Yes`
>       - *"Minimal standard deviation"*: `0.1`
>       - *"Maximal standard deviation"*: `10.0`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. What is the purpose of fitting the peaks to a shape (in this case bi-Gaussian)?
> 2. Why are there two standard deviations?
>
> > ### {% icon solution %} Solution
> >
> > 1. It allows precise computation of the area under the curve and estimating the intensity.
> > 2. Two standard deviations ($$\sigma_1$$ and $$\sigma_2$$) come from a bi-Gaussian shape, where both sides of the shape can be different. When these values are equal, we obtain a Gaussian shape.
> >
> {: .solution}
>
{: .question}

## Compute clusters

The pre-alignment step, where we put all peaks from all samples into a single table and group them based on both m/z and rt. The tool takes a collection of all detected features and computes the clusters over a global feature table, adding the `sample_id` and `cluster` (shared across samples) columns to the table. This process is parametrised by influencing the "size" of buckets (clusters) using relative m/z tolerance and retention time tolerance.

> <details-title> Clustering algorithm details </details-title>
> 
> Features are first grouped in m/z dimension based on the relative m/z tolerance. Then, the absolute tolerance is computed for each feature and a new group is separated once the difference between consecutive features is above this threshold. The same process is then repeated for the retention time dimension. The individual indices are then combined into a single index in the `cluster` columns.
> 
> ![recetox-aplcms clustering](../../images/aplcms_clustering.png "Illustrative example of clustering algorithm with highlighted m/z groups (red) and rt groups (blue). Cluster is the same only when these two groups overlap (green).")
> 
{: .details}

> ### {% icon hands_on %} Hands-on: Compute clusters
>
> 1. {% tool [recetox-aplcms - compute clusters](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_compute_clusters/recetox_aplcms_compute_clusters/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input data"*: `feature tables` (output of **recetox-aplcms - generate feature table** {% icon tool %})
>    - *"Relative m/z tolerance [ppm]"*: `5.0`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. Can we influence the clustering sensitivity w.r.t. retention time?
>
> > ### {% icon solution %} Solution
> >
> > 1. Yes, use retention time tolerance parameter can be used for this purpose.
> >
> {: .solution}
>
{: .question}

Output is again separated into a collection of individual tables but with assigned `cluster`.

## Compute template

To continue with further steps, we need a template into which we can align the data - this means a peak table which has the highest number of features, consequently giving us the highest number of reference points we fit our curve to. This step can be potentially skipped if you want to select a template manually or provide a custom file used in further steps.

**TODO** add hint how to select single file from collection as input for a tool

> ### {% icon hands_on %} Hands-on: Compute template
>
> 1. {% tool [recetox-aplcms - compute template](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_compute_template/recetox_aplcms_compute_template/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input data"*: `clustered features` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>
{: .hands_on}

## Correct time

We need our previously clustered features and selected template.
Apply spline-based retention time correction to a feature table given the template table and the mz and rt tolerances.

> <details-title> Retention time correction effects </details-title>
> 
> The retention time alignment starts by grouping features based on their m/z values using a lenient m/z tolerance. Kernel density estimation is applied to the m/z values within each group to identify multiple modes and split the group accordingly. Next, a retention-time cutoff is determined using a model-based search. Retention time differences are calculated between features within each m/z group. The differences are combined, and a density estimation is performed.
>
> A threshold is determined by finding the largest distance where the observed density exceeds 1.5 times the fitted value. If any spacing in retention time is larger than the threshold, feature groups with similar m/z values are further divided.
> 
> The alignment is performed by selecting one profile as the template (one with the most number of features) and aligning all other profiles against it using a kernel smoother. Adjustments are made to the retention times based on the alignment results, and for features outside the range, their retention times are adjusted with the same amount as the nearest endpoint of the alignment.
>
> <figure>
>    <iframe src="data/time-correction-plot.html" width="800px" height="600px" frameBorder="0"></iframe>
>    <figcaption>The effect of time correction on a sample. The retention time of targeted sample (green) is shifted to corrected values (red) based on the reference template (blue). The data in the plot is reduced by 75% for efficiency reasons.</figcaption>
>  </figure>
>
{: .details}

> ### {% icon hands_on %} Hands-on: Correct time
>
> 1. {% tool [recetox-aplcms - correct time](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_correct_time/recetox_aplcms_correct_time/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input clustered features table"*: `clustered features` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>    - {% icon param-file %} *"Input template features table"*: `reference template` (output of **recetox-aplcms - compute template** {% icon tool %})
>    - *"Relative m/z tolerance [ppm]"*: `5.0`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. Can we determine in advance the value of the retention time shift applied to the sample used as a reference template?
>
> > ### {% icon solution %} Solution
> >
> > 1. The reference template is used to adjust other samples accordingly. However, the template sample itself is not influenced by the correction step. Therefore, the shift value is zero.
> >
> {: .solution}
>
{: .question}

## Compute clusters (2nd round) 

After we have aligned the retention time of our samples, we need to run the second round of clustering to reflect the introduced changes. 

> ### {% icon comment %} Steps iteration
>
> Note that the outputs and inputs of the previous steps are compatible, making it possible to iteratively combine these steps multiple times until desired quality of results is achieved.
{: .comment}

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [recetox-aplcms - compute clusters](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_compute_clusters/recetox_aplcms_compute_clusters/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input data"*: `corrected features` (output of **recetox-aplcms - correct time** {% icon tool %})
>    - *"Relative m/z tolerance [ppm]"*: `3.0`
>    - *"Retention time tolerance [unit corresponds to the retention time]"*: `1.0`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. What already used step could we use next? What would be its effects?
>
> > ### {% icon solution %} Solution
> >
> > 1. Running **recetox-aplcms - compute template** {% icon tool %}, followed by **recetox-aplcms - correct time** {% icon tool %}, has the potential to enhance features' retention time alignment and consequently improve the results of subsequent steps.
> >
> {: .solution}
>
{: .question}

## Features alignment

This step performs feature alignment after clustering and retention time correction. The peaks clustered across samples are grouped based on the given tolerances to create an aligned feature table, connecting identical features across samples. Among tolerances, the `Minimal occurrence in samples` parameter can be used to control at least in how many samples a feature has to be detected in order to be included in the aligned feature table. This allows us to preserve only peaks that appear really consistently across the samples.

> <details-title> Algorithm details </details-title>
> 
> To properly align features, we repeat the model-based search to determine the retention time cutoff and group features accordingly. Within each feature group, we employ kernel density estimators to assess the m/z and retention time dimensions, potentially leading to further divisions within the group. The resulting groups represent aligned features across all profiles. From each group, we extract the median m/z and median retention time as representative characteristic values for the features. This information is recorded in an aligned feature table (for practical purposes, separated into three tables), which includes the median m/z, median retention time, m/z range, and intensities in each profile for every feature.
>
> **TODO** add picture? - taking medians from values in a group
> 
{: .details}

> ### {% icon hands_on %} Hands-on: Align features
>
> 1. {% tool [recetox-aplcms - align features](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_align_features/recetox_aplcms_align_features/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Clustered features"*: `clustered features` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>    - *"Relative m/z tolerance [ppm]"*: `3.0`
>    - *"Retention time tolerance [unit corresponds to the retention time]"*: `1.0`
>
{: .hands_on}

> <details-title> Output files </details-title>
> 
> The output consists of three tables. All tables share an `id` column.
>
> #### Metadata Table
>
> Contains all quantitative data related to specific peaks that have been detected (mean m/z and rt with their maximal and minimal values). The `npeaks` column denotes the number of peaks which have been grouped into this feature. The columns with the sample names indicate whether this feature is present in the sample.
> 
> > |  id   | mz           |  mzmin       |  mzmax        |  rt            |  rtmin        |  rtmax        |   npeaks  |  21_qc_no_dil_milliq   |  29_qc_no_dil_milliq   |  8_qc_no_dil_milliq    |
> > |-------|--------------|--------------|---------------|----------------|---------------|---------------|-----------|------------------------|------------------------|------------------------|
> > |  1    | 70.03707021  |  70.037066   |  70.0370750   |  294.1038014   |  294.0634942  |  294.149985   |   3       |  1                     |  1                     |  1                     |
> > |  2    | 70.06505677  |  70.065045   |  70.0650676   |  141.9560055   |  140.5762528  |  143.335758   |   2       |  1                     |  0                     |  1                     |
> > |  57   | 78.04643252  |  78.046429   |  78.0464325   |  294.0063397   |  293.9406777  |  294.072001   |   2       |  1                     |  1                     |  0                     |
> > |  ...  | ...          |   ...        |  ...          |  ...           |  ...          |  ...          |   ...     |  ...                   |  ...                   |  ...                   |
> {: .matrix}
> 
> #### Intensity Table
> 
> This table contains the peak area for aligned features in all samples.
> 
> > |  id   |  21_qc_no_dil_milliq   |  29_qc_no_dil_milliq   |  8_qc_no_dil_milliq    |
> > |-------|------------------------|------------------------|------------------------|
> > |  1    |  13187487.20482895     |  7957395.699119729     |  11700594.397257797    |
> > |  2    |  2075168.6398983458    |  0                     |  2574362.159289044     |
> > |  57   |  2934524.4406785755    |  1333044.5065971944    |  0                     |
> > |  ...  |  ...                   |  ...                   |  ...                   |
> {: .matrix}
> 
> #### Retention Time Table
> 
> This table contains the retention times for all aligned features in all samples.
>
> > |  id   |  21_qc_no_dil_milliq   |  29_qc_no_dil_milliq   |  8_qc_no_dil_milliq    |
> > |-------|------------------------|------------------------|------------------------|
> > |  1    |  294.09792478513236    |  294.1499853056912     |  294.0634942428341     |
> > |  2    |  140.57625284242982    |  0                     |  143.33575827589172    |
> > |  57   |  294.07200187644435    |  293.9406777222317     |  0                     |
> > |  ...  |  ...                   |  ...                   |  ...                   |
> {: .matrix}
>
{: .details}

> ### {% icon question %} Questions
>
> 1. Why is the outputted aligned feature table separated into three tables?
> 2. How can I find retention times and intensities corresponding to a feature from the metadata table?
>
> > ### {% icon solution %} Solution
> >
> > 1. The purpose of separating the information into multiple tables is to simplify data inspection, eliminate the overhead of too many details in one location, and make the outputs compatible with other tools (such as XCMS).
> > 2. All three tables share the `id` column. The particular value of the `id` column can be used to find corresponding retention times and intensities.
> >
> {: .solution}
>
{: .question}

> ### {% icon comment %} Next steps
>
> At this point, there are two alternative routes how to continue - you can use [unsupervised]({{ site.baseurl }}/topics/metabolomics/tutorials/gc_ms_with_aplcms/tutorial.html#unsupervised) approach (use no existing knowledge, detects peaks de novo from the data based on the data itself) or [hybrid]({{ site.baseurl }}/topics/metabolomics/tutorials/gc_ms_with_aplcms/tutorial.html#hybrid) approach (combine de novo peak detection with existing knowledge).
> 
{: .comment}

# Unsupervised

The unsupervised approach focuses on exploring and discovering patterns in the data without prior knowledge or assumptions. This method employs techniques such as dimensionality reduction and clustering to reveal inherent structures and relationships within the metabolite profiles. The unsupervised approach is useful when the specific patterns or classes within the data are not known beforehand, allowing for unbiased exploration of the data and potential identification of novel patterns or subgroups. However, the unsupervised approach may require additional validation and annotation steps to assign biological meaning to the discovered patterns.

## Recover weaker signals

This step recovers features which are present in a sample but might have been filtered out initially as noise due to low signal intensity. It runs the second stage of peak detection based on the aligned feature table from the feature alignment step. If a feature is contained in the aligned feature table, this step revisits the raw data and searches for this feature at the retention time obtained by mapping the corrected retention time back to the original sample.

> <details-title> Gaps in the feature table </details-title>
>
> After inspecting the obtained metadata table, we can observe that many columns indicating the presence of the feature in the sample are empty (resp., there is a zero). The same holds for respective positions in intensity and retention time tables. That basically means we lost some information for metabolites that did not pass the noise filter and were not included in the feature table during LC/MS profiling. However, that does not necessarily mean they are not present in the samples. To double-check this, we revisit the data and try to recover them.
> 
> The recovery process involves sequentially examining each LC/MS profile. Firstly, an interpolating spline is used to adjust the retention times in the aligned feature table to match the specific profile. Secondly, for features in the aligned table with zero values in the profile, a search is conducted within a smaller region around their m/z and retention time locations. This region is narrower than the range defined by the m/z and retention time across other profiles. Noise filtering is not applied during this step. If one or more features are detected within this region, the feature most consistent with the parameters in the feature table is considered the missed feature, and its intensity replaces the zero value in the aligned feature table.
>
> **TODO** add table and picture: metadata table with many gaps, comment on their meaning... next to it a 3D plot where is a feature in sample but with very low intensity
>
{: .details}

> <details-title> Key parameters </details-title>
> 
> Since most of the previous steps are run again, also their parameters are repeated, giving the opportunity to set them more or less strictly.
>
> Besides them, there is a parameter `Minimal count to recover`, allowing to set the minimum number of raw data points to be considered as a true feature.
>
{: .details}

> ### {% icon hands_on %} Hands-on: Recover weaker signals
>
> 1. {% tool [recetox-aplcms - recover weaker signals](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_recover_weaker_signals/recetox_aplcms_recover_weaker_signals/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input spectra data"*: `mzML input (msconvert)` (output of **msconvert** {% icon tool %})
>    - {% icon param-collection %} *"Input extracted feature samples collection"*: `feature tables` (output of **recetox-aplcms - generate feature table** {% icon tool %})
>    - {% icon param-collection %} *"Input corrected feature samples collection"*: `corrected features` (output of **recetox-aplcms - correct time** {% icon tool %})
>    - {% icon param-file %} *"Metadata table"*: `aligned metadata table` (output of **recetox-aplcms - align features** {% icon tool %})
>    - {% icon param-file %} *"RT table"*: `aligned rt table` (output of **recetox-aplcms - align features** {% icon tool %})
>    - {% icon param-file %} *"Intensity table"*: `aligned intensity table` (output of **recetox-aplcms - align features** {% icon tool %})
>    - *"Relative m/z tolerance [ppm]"*: `10`
>    - *"Retention time tolerance [unit corresponds to the retention time]"*: `5.0`
>    - *"m/z tolerance [ppm]"*: `10`
>    - *"Minimal count to recover"*: `3`
>    - *"Bandwidth factor"*: `0.5`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. For the parameters that are repeated from previous steps, is it neccesary to use the same values for them?
>
> > ### {% icon solution %} Solution
> >
> > 1. **TODO** what is the answer? what is recommended?
> >
> {: .solution}
>
{: .question}

## Compute clusters (3rd round)

We might have added new features, so we do the clustering again.

> ### {% icon hands_on %} Hands-on: Compute clusters
>
> 1. {% tool [recetox-aplcms - compute clusters](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_compute_clusters/recetox_aplcms_compute_clusters/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input data"*: `recovered features` (output of **recetox-aplcms - recover weaker signals** {% icon tool %})
>    - *"Relative m/z tolerance [ppm]"*: `10`
>    - *"Retention time tolerance [unit corresponds to the retention time]"*: `5.0`
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

## Features alignment (2nd round)

Features can now appear in more samples than before, so we also need to repeat the alignment step.

> ### {% icon comment %} Steps iteration
>
> These steps can be again repeated and combined (as well as, for example, combined with retention time correction if necessary) in an arbitrary number of iterations.
>
{: .comment}

> ### {% icon hands_on %} Hands-on: Align features
>
> 1. {% tool [recetox-aplcms - align features](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_align_features/recetox_aplcms_align_features/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Clustered features"*: `clustered features` (output of **recetox-aplcms - compute clusters** {% icon tool %}))
>    - *"Relative m/z tolerance [ppm]"*: `10`
>    - *"Retention time tolerance [unit corresponds to the retention time]"*: `5.0`
>    - *"Minimal occurrence in samples"*: `2`
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

***TODO*** more explanation of the outputs? comparision with pre-recovered data?

# Hybrid

The hybrid approach combines unsupervised techniques with supervised or targeted methods. This approach incorporates external information, such as known metabolic pathways or class labels, to guide the analysis and interpretation of the data. The hybrid approach may leverage prior knowledge to enhance the detection and interpretation of specific metabolite classes, pathways, or biomarkers of interest. The hybrid approach is particularly valuable when there is prior knowledge available or when targeted analysis is desired.

## Merge known table

This step allows us to incorporate the knowledge of known metabolites or historically detected features on the same machinery to help detect and quantify lower-intensity peaks. The features we have found so far are merged with the data frame of known features. Then we will recover the weaker signals and repeat all the previous steps. Finally, the data frame of known metabolites and historical features is updated with the new information from the measured data.

> <details-title> Known table with example </details-title>
> 
> The known table can contain these columns: *chemical_formula*, *HMDB_ID*, *KEGG_compound_ID*, *mass*, *ion.type* (the ion form), *m.z* (either theoretical or mean observed m/z value of previously found features), *Number_profiles_processed* (the total number of processed samples to build this database), *Percent_found* (the percentage of historically processed samples in which the feature appeared), *mz_min* (minimum observed m/z value), *mz_max* (maximum observed m/z value), *RT_mean* (mean observed retention time), *RT_sd* (standard deviation of observed retention time), *RT_min* (minimum observed retention time), *RT_max* (maximum observed retention time), *int_mean.log* (mean observed log intensity), *int_sd.log* (standard deviation of observed log intensity), *int_min.log* (minimum observed log intensity), *int_max.log* (maximum observed log intensity).
> 
> Let us have the known table with the following contents:
> 
> > chemical_formula | HMDB_ID | KEGG_compound_ID | mass | ion.type | m.z |
> > -----------------|---------|------------------|------|----------|-----|
> > ... | ... | ... | ... | ... | ... |
> > C5H10O5 | HMDB00098 | C00181 | 150.05282343 | M+H   | 151.06009942999998 |
> > C8H10N4O3 | HMDB02123 | NA | 210.075290206 | M+H   | 211.082566206 |
> > C18H20N2S | HMDB15038 | C07175 | 296.13471934 | M+H   | 297.14199534 |
> > C16H13ClN2O2 | HMDB14376 | C07125 | 300.066555377 | M+H   | 301.073831377 |
> > C20H28O8 | HMDB32844 | NA | 396.178417872 | M+H   | 397.185693872 |
> > ... | ... | ... | ... | ... | ... |
> {: .matrix}
>
> Then the information about m/z will be added to our feature table for further processing: 
>
> > |  id   | mz           |  mzmin       |  mzmax        |  rt            |  rtmin        |  rtmax        |   npeaks  |  21_qc_no_dil_milliq   |  29_qc_no_dil_milliq   |  8_qc_no_dil_milliq    |
> > |-------|--------------|--------------|---------------|----------------|---------------|---------------|-----------|------------------------|------------------------|------------------------|
> > |  1    | 70.03707021  |  70.037066   |  70.0370750   |  294.1038014   |  294.0634942  |  294.149985   |   3       |  1                     |  1                     |  1                     |
> > |  2    | 70.06505677  |  70.065045   |  70.0650676   |  141.9560055   |  140.5762528  |  143.335758   |   2       |  1                     |  0                     |  1                     |
> > |  57   | 78.04643252  |  78.046429   |  78.0464325   |  294.0063397   |  293.9406777  |  294.072001   |   2       |  1                     |  1                     |  0                     |
> > |  ...  | ...          |   ...        |  ...          |  ...           |  ...          |  ...          |   ...     |  ...                   |  ...                   |  ...                   |
> > |  1871  | 151.06009942999998 |   NA       |  NA          |  NA          |  NA     |  NA     |   NA |  NA              |  NA              |  NA              |
> > |  1872  | 211.082566206 |   NA   |  NA     |  NA      |  NA     |  NA     |   NA |  NA              |  NA              |  NA              |
> > |  1873  | 297.14199534 |   NA   |  NA     |  NA      |  NA     |  NA     |   NA |  NA              |  NA              |  NA              |
> > |  1874  | 301.073831377 |   NA   |  NA     |  NA      |  NA     |  NA     |   NA |  NA              |  NA              |  NA              |
> > |  1875  | 397.185693872 |   NA   |  NA     |  NA      |  NA     |  NA     |   NA |  NA              |  NA              |  NA              |
> > |  ...  | ...          |   ...        |  ...          |  ...           |  ...          |  ...          |   ...     |  ...                   |  ...                   |  ...                   |
> {: .matrix}
>
> Note that the data values are only illustrative.
>
{: .details}

> <details-title> Key parameters </details-title>
> 
> - **Match tolerance [ppm]** - The ppm tolerance to match identified features to known metabolites/features depends on the mass accuracy of your machine.
>
{: .details}

> ### {% icon hands_on %} Hands-on: Merge known table
>
> 1. {% tool [recetox-aplcms - merge known table](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_merge_known_table/recetox_aplcms_merge_known_table/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Metadata table"*: `aligned metadata table` (output of **recetox-aplcms - align features** {% icon tool %})
>    - {% icon param-file %} *"RT table"*: `aligned rt table` (output of **recetox-aplcms - align features** {% icon tool %})
>    - {% icon param-file %} *"Intensity table"*: `aligned intensity table` (output of **recetox-aplcms - align features** {% icon tool %})
>    - {% icon param-file %} *"Table of known features"*: `known table` (Input dataset)
>    - *"Relative m/z tolerance [ppm]"*: `10`
>    - *"Retention time tolerance [unit corresponds to the retention time]"*: `5.0`
>    - *"Tables merge direction"*: `Merge known table to features`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. Can I use my custom database of known metabolites?
>
> > ### {% icon solution %} Solution
> >
> > 1. Yes, you can use a custom database with a similar structure as described above.
> >
> {: .solution}
>
{: .question}

## Recover weaker signals

This step recovers features which are present in a sample but might have been filtered out initially as noise due to low signal intensity, as well as those possibly obtained by merging the known table into our feature table. 

> ### {% icon hands_on %} Hands-on: Recover weaker signals
>
> 1. {% tool [recetox-aplcms - recover weaker signals](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_recover_weaker_signals/recetox_aplcms_recover_weaker_signals/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input spectra data"*: `mzML input (msconvert)` (output of **msconvert** {% icon tool %})
>    - {% icon param-collection %} *"Input extracted feature samples collection"*: `feature tables` (output of **recetox-aplcms - generate feature table** {% icon tool %})
>    - {% icon param-collection %} *"Input corrected feature samples collection"*: `corrected features` (output of **recetox-aplcms - correct time** {% icon tool %})
>    - {% icon param-file %} *"Metadata table"*: `aligned metadata table` (output of **recetox-aplcms - align features** {% icon tool %})
>    - {% icon param-file %} *"RT table"*: `aligned rt table` (output of **recetox-aplcms - align features** {% icon tool %})
>    - {% icon param-file %} *"Intensity table"*: `aligned intensity table` (output of **recetox-aplcms - align features** {% icon tool %})
>    - *"Relative m/z tolerance [ppm]"*: `10`
>    - *"Retention time tolerance [unit corresponds to the retention time]"*: `5.0`
>    - *"m/z tolerance [ppm]"*: `10`
>    - *"Minimal count to recover"*: `3`
>    - *"Bandwidth factor"*: `0.5`
>
{: .hands_on}

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

## Compute clusters (3rd round)

Now we repeat the whole workflow since we might have added new features, and we start with the clustering.

> ### {% icon hands_on %} Hands-on: Compute clusters
>
> 1. {% tool [recetox-aplcms - compute clusters](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_compute_clusters/recetox_aplcms_compute_clusters/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input data"*: `recovered features` (output of **recetox-aplcms - recover weaker signals** {% icon tool %})
>    - *"Relative m/z tolerance [ppm]"*: `10`
>    - *"Retention time tolerance [unit corresponds to the retention time]"*: `5.0`
>
{: .hands_on}

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

## Compute template (2nd round)

With new clusters, we need to prepare a new template for the time correction step.

> ### {% icon hands_on %} Hands-on: Compute template
>
> 1. {% tool [recetox-aplcms - compute template](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_compute_template/recetox_aplcms_compute_template/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input data"*: `clustered features` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>
{: .hands_on}

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

## Correct time (2nd round)

Adding features from the known table might require applying the retention time correction again.

> ### {% icon hands_on %} Hands-on: Correct time
>
> 1. {% tool [recetox-aplcms - correct time](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_correct_time/recetox_aplcms_correct_time/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input clustered features table"*: `clustered features` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>    - {% icon param-file %} *"Input template features table"*: `reference template` (output of **recetox-aplcms - compute template** {% icon tool %})
>    - *"Relative m/z tolerance [ppm]"*: `10`
>    - *"Retention time tolerance [unit corresponds to the retention time]"*: `5.0`
>
{: .hands_on}

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

## Compute clusters (4th round)

Correcting retention time requires updating the clustering.

> ### {% icon hands_on %} Hands-on: Compute clusters
>
> 1. {% tool [recetox-aplcms - compute clusters](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_compute_clusters/recetox_aplcms_compute_clusters/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input data"*: `corrected features` (output of **recetox-aplcms - correct time** {% icon tool %})
>    - *"Relative m/z tolerance [ppm]"*: `10`
>    - *"Retention time tolerance [unit corresponds to the retention time]"*: `5.0`
>
{: .hands_on}

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

## Features alignment (2nd round)

And this all leads to the final computation of features alignment.

> ### {% icon hands_on %} Hands-on: Align features
>
> 1. {% tool [recetox-aplcms - align features](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_align_features/recetox_aplcms_align_features/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Clustered features"*: `clustered features` (output of **recetox-aplcms - compute clusters** {% icon tool %})
>    - *"Relative m/z tolerance [ppm]"*: `10`
>    - *"Retention time tolerance [unit corresponds to the retention time]"*: `5.0`
>    - *"Minimal occurrence in samples"*: `2`
>
{: .hands_on}

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

## Merge known table (2nd round)

After finishing all steps, we can include all new knowledge into the existing database.

> <details-title> Key parameters </details-title>
> 
> - **Minimal occurence of feature** - The number of spectra a new feature must be present for it to be added to the database. We recommend setting this parameter in a stringent manner.
>
{: .details}

> ### {% icon hands_on %} Hands-on: Merge known table
>
> 1. {% tool [recetox-aplcms - merge known table](toolshed.g2.bx.psu.edu/repos/recetox/recetox_aplcms_merge_known_table/recetox_aplcms_merge_known_table/0.10.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Metadata table"*: `aligned metadata table` (output of **recetox-aplcms - align features** {% icon tool %})
>    - {% icon param-file %} *"RT table"*: `aligned rt table` (output of **recetox-aplcms - align features** {% icon tool %})
>    - {% icon param-file %} *"Intensity table"*: `aligned intensity table` (output of **recetox-aplcms - align features** {% icon tool %})
>    - {% icon param-file %} *"Table of known features"*: `known table` (Input dataset)
>    - *"Relative m/z tolerance [ppm]"*: `10`
>    - *"Retention time tolerance [unit corresponds to the retention time]"*: `5.0`
>    - *"Tables merge direction"*: `Merge features to known table`
>       - *"new_feature_min_count"*: `2`
>
{: .hands_on}

> <details-title> Example of updated known table </details-title>
> 
> Example of how can be new knowledge included in the known table:
> 
> > chemical_formula | HMDB_ID | KEGG_compound_ID | mass | ion.type | m.z | Number_profiles_processed | Percent_found | mz_min | mz_max | RT_mean | RT_sd | RT_min | RT_max | int_mean(log) | int_sd(log) | int_min(log) | int_max(log)
> > -----------------|---------|------------------|------|----------|-----|---------------------------|---------------|--------|--------|---------|-------|--------|--------|---------------|-------------|---------------|--------------|
> > ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
> > C5H10O5 | HMDB00098 | C00181 | 150.05282343 | M+H   | 151.06009942999998 | 3.0 | 1.0 | 151.06009942999998 | 151.06172515032375 | 191.2551471245085 | 0.4430276002417259 | 190.8004248337492 | 191.68547177224977 | 9.37396208102855 | 0.0247741820684487 | 9.347567941998712 | 9.396712797272118
> > C8H10N4O3 | HMDB02123 | NA | 210.075290206 | M+H   | 211.082566206 | 3.0 | 1.0 | 211.07901081624345 | 211.082566206 | 167.9441831090941 | 0.4031343080704094 | 167.54668517707515 | 168.3527267744501 | 7.199170317946007 | 0.0254833388951392 | 7.1698777762386285 | 7.216237500092092
> > C18H20N2S | HMDB15038 | C07175 | 296.13471934 | M+H   | 297.14199534 | 3.0 | 1.0 | 297.1371550623409 | 297.14199534 | 109.97665012187495 | 95.2526270699826 | 0.0 | 166.34891329963853 | 6.332251446676835 | 0.2113801719274527 | 6.182783093698555 | 6.481719799655116
> > C16H13ClN2O2 | HMDB14376 | C07125 | 300.066555377 | M+H   | 301.073831377 | 3.0 | 1.0 | 301.0736065580818 | 301.07921477804985 | 115.87415505877158 | 100.34999045612813 | 0.0 | 173.88690693535776 | 6.828260757688546 | 0.0399063569267865 | 6.800042702093164 | 6.856478813283927
> > C20H28O8 | HMDB32844 | NA | 396.178417872 | M+H   | 397.185693872 | 3.0 | 1.0 | 397.1817069960474 | 397.185693872 | 109.22118480492205 | 94.5892574861945 | 0.0 | 164.2527571907555 | 6.610149642385979 | 0.0322116360186999 | 6.587372576124044 | 6.632926708647915
> > ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
> {: .matrix}
>
> Note that the data values are only illustrative.
>
{: .details}

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

# Conclusion
{:.no_toc}

**TODO** comment on resulting tables
