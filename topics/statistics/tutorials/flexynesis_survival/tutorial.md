---
layout: tutorial_hands_on

title: Identifing Survival Markers of Brain tumor with Flexynesis
zenodo_link: https://zenodo.org/records/16287482
questions:
- How can multi-modal genomic data be integrated to identify survival markers in brain tumors?
- How can deep learning approaches improve survival prediction in cancer patients?
- Which genomic features are most predictive of patient survival in glioma subtypes?
objectives:
- Apply multi-modal data integration techniques to combine mutation, gene expression, and clinical data
- Perform survival analysis to identify prognostic biomarkers in brain tumors
- They are single sentences describing what a learner should be able to do once they
  have completed the tutorial
- Implement deep learning models for survival prediction using genomic data
time_estimation: 2H
key_points:
- Multi-modal data integration improves survival prediction accuracy compared to single data types
- Deep learning approaches can effectively handle high-dimensional genomic data for survival analysis
contributors:
- Nilchia
- bgruening
subtopic: flexynesis
answer_histories:
    - label: "usegalaxy.eu"
      history: https://usegalaxy.eu/u/nilchia/h/final-survival-markers-of-lower-grade-gliomas-test-wf
      date: 2025-08-01
---

Here, we use Flexynesis tool suite on a multi-omics dataset of 506 Brain Lower Grade Glioma (LGG) and 288 Glioblastoma Multiforme (GBM) samples with matching mutation and copy number alteration. This data were downloaded from the [cBioPortal](https://www.cbioportal.org) ({% cite cbioportal-website %}). The data was split into train (70% of the samples) and test (30% of the samples).

We will work on CNA (copy number alteration) and MUT (mutation data which is a binary matrix of genexsample) data.

This training is inspired from the original flexynesis analysis notebook: [survival_subtypes_LGG_GBM.ipynb](https://github.com/BIMSBbioinfo/flexynesis/blob/main/examples/tutorials/survival_subtypes_LGG_GBM.ipynb).

> <comment-title> Getting data from cBioPortal </comment-title>
>
> If you want to download data directly from cBioPortal, you can go through "**Prepare data from CbioPortal for Flexynesis integration**" training
{: .comment}

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

> <warning-title>LICENSE</warning-title>
>
>Flexynesis is only available for NON-COMMERCIAL use. Permission is only granted for academic, research, and educational purposes. Before using, be sure to review, agree, and comply with the license.
>For commercial use, please review the flexynesis license on GitHub and contact the [copyright holders](https://github.com/BIMSBbioinfo/flexynesis)
{: .warning}

# Get data

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}):
>
>    ```
>    https://zenodo.org/records/16287482/files/train_mut_lgggbm.tabular
>    https://zenodo.org/records/16287482/files/train_clin_lgggbm.tabular
>    https://zenodo.org/records/16287482/files/test_cna_lgggbm.tabular
>    https://zenodo.org/records/16287482/files/test_mut_lgggbm.tabular
>    https://zenodo.org/records/16287482/files/test_clin_lgggbm.tabular
>    https://zenodo.org/records/16287482/files/train_cna_lgggbm.tabular
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 4. Check that the datatype is `tabular`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tabular" %}
>
> 5. Add to each database a tag corresponding to their modality (#mut, #cna, and #clinical)
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

# Train Flexynesis model with a survival task

Now it is time import train and test datasets. We rank genes by Laplacian Scores and pick top 10% of the genes, while removing highly redundant genes with a correlation score threshold of 0.8 and a variance threshold of 50%. We will also use the intermediate fusion of omic layers.
For hyperparameter optimization we use the following parameters.


*`Model class`: We pick DirectPred (a fully connected network)
* `Column name in the train clinical data to use as survival event` and `Column name in the train clinical data to use as survival time`. Survival status (consists of 0s and 1s) and time since last followup. It is **important** that the clinical data contains both variables as ***numerical*** values.
* `Column name in the train clinical data to use for predictions, multiple targets are allowed`: We can concurrently train the same network to be able to predict other variables such as histological diagnosis, however, here we just focus on the survival endpoints, so we pass an empty list.
* `Number of iterations for hyperparameter optimization.`: We do 5 iterations of hyperparameter optimization. This is a reasonable number of demonstration purposes, but it could be beneficial to increase this value in order to discover even better models.
* `How many epochs to wait when no improvements in validation loss are observed. `: If a training does not show any signs of improving the performance on the validation part of the train_dataset for at least 10 epochs, we stop the training. This not only significantly decreases the amount spent on training by avoiding unnecessary continuation of unpromising training runs, but also helps avoid over-fitting the network on the training data.


> <hands-on-title> Survival task </hands-on-title>
>
> 1. {% tool [Flexynesis](toolshed.g2.bx.psu.edu/repos/bgruening/flexynesis/flexynesis/0.2.20+galaxy3) %} with the following parameters:
>    - *"I certify that I am not using this tool for commercial purposes."*: `Yes`
>    - *"Type of Analysis"*: `Supervised training`
>        - {% icon param-file %} *"Training clinical data"*: `train_clin_lgggbm.tabular`
>        - {% icon param-file %} *"Test clinical data"*: `test_clin_lgggbm.tabular`
>        - {% icon param-file %} *"Training omics data"*: `train_mut_lgggbm.tabular`
>        - {% icon param-file %} *"Test omics data"*: `test_mut_lgggbm.tabular`
>        - *"What type of assay is your input?"*: `mut`
>        - In *"Multiple omics layers?"*:
>            - {% icon param-repeat %} *"Insert Multiple omics layers?"*
>                - {% icon param-file %} *"Training omics data"*: `train_cna_lgggbm.tabular`
>                - {% icon param-file %} *"Test omics data"*: `test_cna_lgggbm.tabular`
>                - *"What type of assay is your input?"*: `cna`
>        - *"Model class"*: `DirectPred`
>        - *"Column name in the train clinical data to use as survival event"*: `c8`
>        - *"Column name in the train clinical data to use as survival time"*: `c7`
>        - In *"Advanced Options"*:
>            - *"Fusion method"*: `intermediate`
>            - *"Variance threshold (as percentile) to drop low variance features."*: `0.5`
>            - *"Correlation threshold to drop highly redundant features."*: `0.8`
>            - *"Minimum number of features to retain after feature selection."*: `1000`
>            - *"Top percentile features (among the features remaining after variance filtering and data cleanup) to retain after feature selection."*: `10.0`
>            - *"Number of iterations for hyperparameter optimization."*: `5`
>        - In *"Visualization"*:
>            - *"Generate embeddings plot?"*: `No`
>            - *"Generate kaplan meier curves plot?"*: `Yes`
>            - *"Generate hazard ratio plot?"*: `Yes`
>                - *"Omics layer to use for cox input"*: `mut`
>                - *"Number of top important features to include in Cox model"*: `5`
>                - *"Performs K-fold cross-validation?"*: `No`
>            - *"Generate scatter plot?"*: `No`
>            - *"Generate concordance heatmap plot?"*: `No`
>            - *"Generate precision-recall curves plot?"*: `No`
>            - *"Generate ROC curves plot?"*: `No`
>            - *"Generate boxplot?"*: `No`
>
{: .hands_on}

> <question-title></question-title>
>
> 1. What are main outputs of Flexynesis?
> 2. What does the generated plots show?
> 3. What information is in `job.stats` output?
> 4. What are top 10 features in `job.feature_importance.GradientShap`?
>
> > <solution-title></solution-title>
> >
> > 1. The results collection contains the following data:
> > * job.embeddings_test (latent space of test data)
> > * job.embeddings_train (latent space of train data)
> > * job.feature_importance.GradientShap (feature importance calculated by GradientShap method)
> > * job.feature_importance.IntegratedGradients (feature importance calculated by IntegratedGradients method)
> > * job.feature_logs.cna
> > * job.feature_logs.mut
> > * job.predicted_labels (Prediction of the created model)
> > * job.stats
> >
> > 2. There are two figures in plot collection:
> > * The first plot is the Kaplan Meier Curves of the risk subtypes.
> > * The second plot is a forest plot of calculated Cox-PH model of top 5 markers.
> >
> > 3. It shows that we achieved a reasonable performance of Harrell's Concordance Index (C-index) on test data.
> >
> > 4. *IDH2*, *IDH1*, *HIVEP3*, *PIK3CA*, *EGFR*, *RELN*, *ATRX*, *INSRR*, *TP53*, *COL6A3*. *IDH1*, *TP53*, and *ATRX* are extensively studied and been shown to be relevant in gliomas ({% cite Koschmann2016 %}).
> >
> {: .solution}
>
{: .question}

We can also visualize the high and low risk groups with the embeddings.

# Survival-risk subtypes

## Extract predicted labels and calculate median

Let's group the samples by predicted survival risk scores into 2 groups and visualize the sample embeddings colored by risk subtypes.

> <hands-on-title> Calculate median of predicted survival score</hands-on-title>
>
> First we should extract the `job.predicted_labels` and then calculate the median.
>
> 1. {% tool [Extract dataset](__EXTRACT_DATASET__) %} with the following parameters:
>    - {% icon param-file %} *"Input List"*: `results` (output of **Flexynesis** {% icon tool %})
>    - *"How should a dataset be selected?"*: `Select by element identifier`
>        - *"Element identifier:"*: `job.predicted_labels`
>
> 2. {% tool [Datamash](toolshed.g2.bx.psu.edu/repos/iuc/datamash_ops/datamash_ops/1.9+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input tabular dataset"*: `job.predicted_labels` (output of **Extract dataset** {% icon tool %})
>    - *"Input file has a header line"*: `Yes`
>    - In *"Operation to perform on each group"*:
>        - {% icon param-repeat %} *"Insert Operation to perform on each group"*
>            - *"Type"*: `Median`
>            - *"On column"*: `Column: 6` (*"predicted_label)
>
{: .hands_on}

> <question-title></question-title>
>
> 1. What is the median?
>
> > <solution-title></solution-title>
> >
> > 1. -0.66496020555496
> >
> {: .solution}
>
{: .question}

Now we add a column to our `job.predicted_labels` data indicating high and low risk groups (samples with predicted_label > median are high risk groups)

## Assign high and low risk groups

> <hands-on-title> Assign high/low groups </hands-on-title>
>
> 1. {% tool [Compute](toolshed.g2.bx.psu.edu/repos/devteam/column_maker/Add_a_column1/2.1) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `job.predicted_labels` (output of **Extract dataset** {% icon tool %})
>    - *"Input has a header line with column names?"*: `Yes`
>        - In *"Expressions"*:
>            - {% icon param-repeat %} *"Insert Expressions"*
>                - *"Add expression"*: `float(c6)>=-0.66496020555496`
>                - *"The new column name"*: `risk_groups`
>    - In *"Error handling"*:
>        - *"Autodetect column types"*: `No`
>
>
{: .hands_on}

For better visualization, let's change True to High_risk and False to Low_risk

> <hands-on-title> Rename labels </hands-on-title>
>
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_column/9.5+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `table` (output of **Compute** {% icon tool %})
>    - In *"Replacement"*:
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"in column"*: `Column: 9`
>            - *"Find pattern"*: `True`
>            - *"Replace with"*: `High_risk`
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"in column"*: `Column: 9`
>            - *"Find pattern"*: `False`
>            - *"Replace with"*: `Low_risk`
>
{: .hands_on}

Now we visualize the embedding.

> <hands-on-title> Dimension reduction plot </hands-on-title>
>
> 1. {% tool [Extract dataset](__EXTRACT_DATASET__) %} with the following parameters:
>    - {% icon param-file %} *"Input List"*: `results` (output of **Flexynesis** {% icon tool %})
>    - *"How should a dataset be selected?"*: `Select by element identifier`
>        - *"Element identifier:"*: `job.embeddings_test`
>
> 2. {% tool [Flexynesis plot](toolshed.g2.bx.psu.edu/repos/bgruening/flexynesis_plot/flexynesis_plot/0.2.20+galaxy3) %} with the following parameters:
>    - *"I certify that I am not using this tool for commercial purposes."*: `Yes`
>    - *"Flexynesis plot"*: `Dimensionality reduction`
>        - {% icon param-file %} *"Predicted labels"*: `table` (output of **Replace Text** {% icon tool %})
>        - {% icon param-file %} *"Embeddings"*: `job.embeddings_test` (output of **Extract dataset** {% icon tool %})
>        - *"Column in the labels file to use for coloring the points in the plot"*: `Column: 9`
>
{: .hands_on}

# Conclusion

In this tutorial, we have successfully applied the Flexynesis framework to analyze survival markers in brain tumors using multi-modal genomic data from 506 Lower Grade Glioma (LGG) and 288 Glioblastoma Multiforme (GBM) samples. Through hands-on analysis of mutation and copy number alteration data, we have separated our samples to high and low risk groups and ranked the genes by their importance from which some are already known to affect survival in glioma.