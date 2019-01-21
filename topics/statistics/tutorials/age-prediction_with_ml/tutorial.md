---
layout: tutorial_hands_on

title: Age prediction using machine learning
zenodo_link: https://zenodo.org/record/2545213#.XEWTJ9-YVa0
questions:
- How to use machine learning to create predictive models from biological datasets (RNA-seq and DNA methylation)?
objectives:
- Learn aging biomarkers from RNA-seq and DNA methylation datasets
- Apply regression based machine learning algorithms
- Learn feature selection and hyperparameter optimisation
time_estimation: 1H
contributors:
- polkhe
- anuprulez

---


# Introduction
{:.no_toc}

[Machine Learning](https://en.wikipedia.org/wiki/Machine_learning) is used to create predictive models by learning features from datasets. In this tutorial, we will apply a couple of (scikit-learn) machine learning tools to [RNA-seq](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1599-6#Sec9) and [DNA methylation](https://www.sciencedirect.com/science/article/pii/S1872497317301643?via%3Dihub) datasets to predict chronological age of humans. Using these tools in Galaxy, we can achieve comparable prediction scores as achieved by these analyses. The RNA-seq gene expression ([FPKM](https://www.ebi.ac.uk/training/online/glossary/fpkm)) dataset is generated using fibroblast cell lines from 133 healthy humans and their ages range from 1 to 94 years. The biomarkers of age are identified by a machine learning algorithm to create a age prediction model. Within each individual, DNA methylation changes with age. This knowledge is used to create an age prediction model using DNA methylation dataset. The CpGs sites with highest correlation to age are selected as biomarkers/features. These are used by a machine learning algorithm to create a predictive model. This tutorial is divided into two parts - one with RNA-seq and another with DNA methylation datasets.


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Analyze RNA-seq data 

The data on which we perform our first task of hyperparameter estimation is a RNAseq data of firoblast cell lines belonging to 133 healthy patients
of age from 1 to 94 years. On this data we perform an exhaustive search (known as grid search) for finding the best features in the dataset and then apply ElasticNet regressor with 5-fold crossvalidation. The R2 regression score is compared to the predictions found in the [original paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1599-6#Sec9).

## Get the raw data

We proceed to the analysis with uploading the data.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]()
>
>    ```
>    https://zenodo.org/record/2545213/files/training_data_normal.tsv?download=1
>    ```
>
>    {% include snippets/import_via_link.md %}
>
> 3. Rename the dataset to `training_data_normal`.
>
>    {% include snippets/rename_dataset.md %}
>
> 4. Check that the datatype is `tabular`
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
{: .hands_on}

## Pre-processing

We can see that the RNA-seq dataset is high-dimensional. There are over `27,000` columns/features. Generally, not all the columns in the dataset are useful for prediction. We need only those columns/features which increases the predictive ability of the model. To filter these columns, we perform feature selection and retain only those columns which are useful. To do that, we use `SelectKBest` module in the suite of data preprocessors. Again, we are not sure of how many of these columns we will need. To find the right number of columns, we do a hyperparameter search by setting different number of features and find out the best number. To create this preprocessor, we will use **pipeline builder** tool. This tool defines a sequential processing of datasets. After preprocessing step, we should add a regressor algorithm (ElasticNet) which analyzes the preprocessed dataset. This tool gives a `zip` file as output containing the specifications of different steps of the entire analysis. 

> ### {% icon hands_on %} Hands-on: Data pre-processing
>
> 1. **Pipeline Builder** {% icon tool %} with the following parameters:
>    - In *"1: Pre-processing step:"*:
>        - *"Choose the type of transformation:"*: `Feature Selection`
>            - *"Select a feature selection algorithm:"*: `SelectKBest - Select features according to the k highest scores`
>                - *"Select a score function:"*: `f_regression - Univariate linear regression tests`
>    - In *"Final Estimator:"*:
>        - *"Choose the module that contains target estimator:"*: `sklearn.linear_model`
>            - *"Choose estimator class:"*: `ElasticNet`
>
>    > ### {% icon comment %} Comment
>    >
>    > *ElasticNet* is a regularization method that combines lasso and ridge regression approaches.
>    {: .comment}
>
{: .hands_on}

## Hyperparameter optimisation

In any machine learning algorithm, there are many parameters (hyperparameters) whose values we are not sure of. There are default values given for these parameters but they may not be optimal for different datasets. To find the best combination of the values of different parameters, [hyperparameter optimisation](https://en.wikipedia.org/wiki/Hyperparameter_optimization) is performed. There are different techniques to do hyperparameter optimisation:
- Grid search
- Random search

For our analyses, we will use grid search approach. It is an exhaustive search which tries out all the combinations of different hyperparameters and ranks these combinations based on a scoring metric.

> ### {% icon details %} Cross-validation
>
> *Cross-validation* is a model validation technique which estimates the performance of a predictive model on an unseen data. A dataset is divided into `k` folds and these folds are divided over a training and validation sets. The performance is averaged over `k` folds.
>
{: .details}

In the pipeline builder, we added two steps - preprocessing (feature selection) and an estimator (regressor). There are different hyperparameters for these two steps and their best combination should be found out. We will perform grid search to estimate the best values for these parameters: **k**, **normalize**, **alpha**. For each parameter, we need to specify a set of values to choose from:
- **k**: [5880, 5890, 5895, 5900]
- **normalize**: [True, False]
- **alpha**: [0.00001, 0.0001, 0.001]

For these three parameters, we have 24 different combinations (4 x 2 x 3) of values and we will verify the performance of each combination. You might have noticed that the parameter **k** is used for feature selection and parameters **normalize** and **alpha** are used for regressor.

> ### {% icon hands_on %} Hands-on: Hyperparameter search
>
> 1. **Hyperparameter Search** {% icon tool %} with the following parameters:
>    - *"Select a model selection search scheme:"*: `GridSearchCV - Exhaustive search over specified parameter values for an estimator ` 
>        - *"Choose the dataset containing pipeline object"*: `Pipeline builder` zipped file
>        - In *"Search parameters Builder"*:
>            - In *"Parameter setting for search:"*:
>                - {% icon param-repeat %} *"Insert Parameter setting for search:"*
>                    - *"Choose the transformation the parameter belongs to"*: `Pre-processing step #1`
>                        - *"Pre_processing component #1  parameter:"*: `k: [5880, 5890, 5895, 5900]`
>                - {% icon param-repeat %} *"Insert Parameter setting for search:"*
>                    - *"Choose the transformation the parameter belongs to"*: `Final estimator`
>                        - *"Estimator parameter:"*: `normalize: [True, False]`
>                - {% icon param-repeat %} *"Insert Parameter setting for search:"*
>                    - *"Choose the transformation the parameter belongs to"*: `Final estimator`
>                        - *"Estimator parameter:"*: `alpha: [0.00001, 0.0001, 0.001]`
>        - In *"Advanced Options for SearchCV"*:
>            - *"Select the primary metric (scoring):"*: `Regression -- 'r2'`
>            - *"Select the cv splitter:"*: `KFold`
>                - *"n_splits"*: `5`
>                - *"Whether to shuffle data before splitting"*: `Yes`
>                - *"Random seed number"*: `3111696`
>            - *"Raise fit error:"*: `Yes`
>    - *"Select input type:"*: `tabular data`
>        - *"Training samples dataset:"*: `training_data_normal` tabular file
>        - *"Does the dataset contain header:"*: `Yes`
>        - *"Choose how to select data by column:"*: `All columns BUT by column header name(s)`
>            - *"Type header name(s):"*: `age`
>        - *"Dataset containing class labels or target values"*: `training_data_normal` tabular file
>        - *"Does the dataset contain header:"*: `Yes`
>        - *"Choose how to select data by column:"*: `Select columns by column header name(s)`
>            - *"Type header name(s):"*: `age`
>
>    > ### {% icon comment %} Comment
>    >
>    > The tool returns two outputs, one of which is a table with numerical results. Inspect it carefully: the last column shows the ranking of settings
based on performance. The ranking is based on the the 14th column where you can find mean test score and the values for parameters which output this result
are visible in the 8th column.
>    {: .comment}
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. What is the 'best' possible mean_test_score estimated by the tool for these parameters?
> 2. Which combination of parameter settings gives it?
> 3. How many combinations of possible parameters the tool estimated?
>
> > ### {% icon solution %} Solution
> >
> > 1. 0.7269799945732331
> > 2. alpha: 0.001, normalize: True, k: 5880
> > 3. 24
> >
> {: .solution}
>
{: .question}

## Parallel coordinates plot

We will visualize the tabular output of hyperparameter search tool from the previous step using **Parallel Coordinates Plot of tabular data**.

> ### {% icon hands_on %} Hands-on: Parallel coordinates plot
>
> 1. **Parallel Coordinates Plot** {% icon tool %} with the following parameters:
>    - *"Select data file:"*: Tabular output of hyperparameter search tool
>    - *"Select the columns for dimensions:"*: `c[5, 6, 7, 14]`
>    - *"Select a column containing the values for coloring:"*: `c14`
>
>    > ### {% icon comment %} Comment
>    >
>    > The output plot has the following legend: the colour-coding is based on the `mean_test_score` column. You can follow the line leading
to the score along every column with parameters' settings.
>
>
>    ![data](images/plotting_output.png "The vizualization of the hyperparameter optimization tool output.")
>
>    {: .comment}
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. What can you notice about the least performing (let's say least four) hyperparameters' settings (judging by the plot)?
>
> > ### {% icon solution %} Solution
> >
> > 1. The four 'worst' settings are:
> > - alpha: 0.00001, normalize: False, k: 5880
> > - alpha: 0.00001, normalize: False, k: 5890
> > - alpha: 0.00001, normalize: False, k: 5895
> > - alpha: 0.00001, normalize: False, k: 5900
> >
> {: .solution}
>
{: .question}

## Conclusion
In the plot shown above, we achieved an R2 score of `0.73` (last column) with 5-fold cross-validation. In the [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1599-6#Sec9) as well, a similar R2 score is mentioned for linear regressors. It showcases one of the usecases of using the scikit-learn machine learning tools in Galaxy to reproduce the results published in a paper.

# Preparing for the prediction

This second part of the analysis is covering the age prediction task. We start with repeating the same already familiar steps from the first part of the tutorial to estimate hyperparameters. The next stage is the actual training, using the tool for ensemble methods for classification and regression. Then we
compare the results with available test labels in order to calculate residuals.

## Get the train and test data

We proceed to the analysis with uploading new data. You might want to create a new history first.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]()
>
>    ```
>    https://zenodo.org/api/files/0d468136-5025-4c0f-bf8b-a8277a513a93/test_rows.csv
>    https://zenodo.org/api/files/0d468136-5025-4c0f-bf8b-a8277a513a93/test_rows_labels.csv
>    https://zenodo.org/api/files/0d468136-5025-4c0f-bf8b-a8277a513a93/train_rows.csv
>    ```
>
>    {% include snippets/import_via_link.md %}
>
> 3. Rename the datasets accordingly.
>
>    {% include snippets/rename_dataset.md %}
>
> 4. Check that the datatypes.
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
>    > ### {% icon comment %} Comment
>    >
>    > The `train_rows.csv` contains an additional column with ages, which is used for the training. We will estimate our model on
>    > `test_rows.csv` and compare the predicted age with labels in `test_rows_labels.csv`
>
>    {: .comment}
{: .hands_on}

## Pre-processing with *Pipeline Builder*

We move on to pre-processing with re-running **Pipeline Builder** on the new data to setup the [Gradient Boosting Regressor](https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.GradientBoostingRegressor.html).

> ### {% icon hands_on %} Hands-on: Data pre-processing
>
> 1. **Pipeline Builder** {% icon tool %} with the following parameters:
>    - In *"Final Estimator:"*:
>        - *"Choose the module that contains target estimator:"*: `sklearn.ensemble`
>            - *"Choose estimator class:"*: `GradientBoostingRegressor`
>
>    > ### {% icon comment %} Comment
>    >
>    > [*Ensemble*](https://en.wikipedia.org/wiki/Ensemble_learning) methods allow to use several learning models for better predictions. 
>    {: .comment}
>
{: .hands_on}

## *Hyperparameter Search*: the training

Before we can start testing, we need to train the model on the input data. At the same time the familiar **Hyperparameter Search** tool allows to estimate the optimal number of learners for the ensemble, which is a hyperparameter.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Hyperparameter Search** {% icon tool %} with the following parameters:
>    - *"Select a model selection search scheme:"*: `GridSearchCV - Exhaustive search over specified parameter values for an estimator `
>        - In *"Search parameters Builder"*:
>            - In *"Parameter setting for search:"*:
>                - {% icon param-repeat %} *"Insert Parameter setting for search:"*
>                    - *"Choose the transformation the parameter belongs to"*: `Final estimator`
>                        - *"Estimator parameter:"*: `n_estimators: [25, 50, 75, 100, 200]`
>        - In *"Advanced Options for SearchCV"*:
>            - *"Select the primary metric (scoring):"*: `Regression -- 'r2'`
>            - *"Select the cv splitter:"*: `KFold`
>                - *"n_splits"*: `5`
>                - *"Whether to shuffle data before splitting"*: `Yes`
>                - *"Random seed number"*: `3111696`
>            - *"Raise fit error:"*: `Yes`
>    - *"Select input type:"*: `tabular data`
>        - *"Does the dataset contain header:"*: `Yes`
>        - *"Choose how to select data by column:"*: `All columns BUT by column header name(s)`
>            - *"Type header name(s):"*: `Age`
>        - *"Does the dataset contain header:"*: `Yes`
>        - *"Choose how to select data by column:"*: `Select columns by column header name(s)`
>            - *"Type header name(s):"*: `Age`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. What is the predicted optimal number of estimators?
>
> > ### {% icon solution %} Solution
> >
> > 1. 75
> >
> {: .solution}
>
{: .question}

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Parallel Coordinates Plot** {% icon tool %} with the following parameters:
>    - *"Select the columns for dimentions:"*: `c[5, 12]`
>    - *"Select a column containg the values for coloring:"*: `c12`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. What hyperparameter value returned the lowest mean_test_score (according to the plot)?
>
> > ### {% icon solution %} Solution
> >
> > 1. n = 25
> >
> {: .solution}
>
{: .question}

# Ensemble prediction

Now when we prepared the trained model, the actual prediction can be launched.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Ensemble methods for classification and regression** {% icon tool %} with the following parameters:
>    - *"Select a Classification Task:"*: `Load a model and predict`
>        - *"Models"*: to the `output from Hyperparameter Search`
>        - *"Data (tabular)"*: `test_rows.csv`
>        - *"Does the dataset contain header:"*: `Yes`
>        - *"Select the type of prediction:"*: `Predict class labels`
>
{: .hands_on}

Let's plot the predictions and compare with the test labels.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Plot actual vs predicted curves and residual plots of tabular data** {% icon tool %} with the following parameters:
>    - *"Select input data file :"*: `test_rows_labels.csv`
>    - *"Select predicted data file :"*: to the `output from the previous tool`
>
>    > ### {% icon comment %} Comment
>    >
>    > The tool outputs three html files with the interactive plots.
>    {: .comment}
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. Inspect the plots. What can you say about our predictions?
>
> > ### {% icon solution %} Solution
> >
> > 1. The predictions (marked orange) show good results.
> >
> {: .solution}
>
{: .question}


# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
