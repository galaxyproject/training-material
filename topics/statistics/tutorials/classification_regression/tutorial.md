---
layout: tutorial_hands_on

title: "Machine learning: classification and regression"
questions:
  - "what are classification and regression techniques?"
  - "How they can be used for prediction?"
  - "How visualizations can be used to analyze predictions?"
objectives:
  - "Explain the types of supervised machine learning - classification and regression."
  - "Learn how to make predictions using the training and test data."
  - "Visualize the predictions."
requirements:
time_estimation: "1H"
key_points:
  - "Learn machine learning's supervised approaches - classification and regression."
  - "In supervised approaches, the target for each sample is known."
  - "For classification and regression tasks, data is divided into training and test sets."
  - "Using classification, the categories of rows are learned using the training set and predicted using the test set."
  - "Using regression, real-valued targets are learned using the training set and predicted using the test set."
contributors:
  - anuprulez

---

# Introduction
{:.no_toc}

Supervised learning methods in machine learning have targets/classes/categories defined in the datasets. These targets can either be discreet values or real-values (continuous). When the targets are discreet, the learning task is called as classification. When the targets are real-values, the task becomes regression. Classification is assigning a distinct category to each sample in the dataset. Regression assigns a real-valued output to a sample in the dataset. In the image below, the "green" line is a boundary which separates the blue balls from the red ones. The task of a classification method is to learn this boundary which can be used to differentiate between unseen blue and red balls.

![Dataset](images/classification_1.png)


The following image shows how a (regression) curve is fit which explains most of the data points. Here, the curve is a straight line. The regression task is to learn this curve which explains the underlying data distribution.

![Dataset](images/regression_1.png)


> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Classification

Classification task assigns a category/class to a sample by learning a decision boundary using a dataset. This dataset is called a training dataset and contains a class/category for each sample. The algorithm which performs this task is called a classifier. The training dataset contains "features" as columns and a mapping between these features and the target is learned for each sample. The performance of mapping is evaluated using test dataset. The test dataset contains only the feature columns and not the target column. The target column is predicted using the mapping learned on the training dataset. In this tutorial, we will use a classifier to train a model using a training dataset, predict the targets for test dataset and visualize the results using plots.

## Data upload

The datasets required for this tutorial contain 9 features of breast cancer which include the thickness of clump, cell-size, cell-shape and so on ([more information](https://github.com/EpistasisLab/penn-ml-benchmarks/tree/master/datasets/classification/breast-w)). In addition to these features, the training dataset contains one more column as `target`. It has a binary value (0 or 1) for each row. `0` indicates no breast cancer and `1` indicates breast cancer. The test dataset does not contain the `target` column. The third dataset contains the all the samples from test dataset but also the `target` column which would be needed to create plot showing comparison between actual and predicted targets.


> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the following datasets and choose the type of data as `tabular`.
> 
>    ```
>    https://zenodo.org/record/1401230/files/breast-w_train.tsv
>    https://zenodo.org/record/1401230/files/breast-w_test.tsv
>    https://zenodo.org/record/1401230/files/breast-w_targets.tsv
>    ```
> 
>    {% include snippets/import_via_link.md %}
>
> 3. Rename datasets to `breast-w_train`, `breast-w_test` and `breast-w_targets`.
>
>    {% include snippets/rename_dataset.md %}
>
{: .hands_on}

## Train a classifier
In this step, `SVM (Support vector machine)` classifier is trained using `breast-w_train` dataset. The last column of this dataset assigns a category for each row. The classifier learns a mapping between each row and its category. This mapping is called a trained model. It is used to predict the categories of unseen data (`breast-w_test`).

> ### {% icon hands_on %} Hands-on: Train a classifier
> 
> **SVM Classifier (Support vector machine)** {% icon tool %} with the following parameters
> 1. {% icon param-select %} *"Select a Classification Task"*: `Train a model`
> 2. {% icon param-select %} *"Classifier type"*: `Linear Support Vector Classification`
> 3. {% icon param-select %} *"Select input type"*: `tabular data`
> 4. {% icon param-file %} *"Training samples dataset"*: `breast-w_train`
> 5. {% icon param-check %} *"Does the dataset contain header"*: `Yes`
> 6. {% icon param-select %} *"Choose how to select data by column"*: `All columns but by column header name(s)`
> 7. {% icon param-text %} *"Type header name(s)"*: `target`
> 8. {% icon param-file %} *"Dataset containing class labels"*: `breast-w_train`
> 9. {% icon param-check %} *"Does the dataset contain header"*: `Yes`
> 10. {% icon param-select %} *"Choose how to select data by column"*: `Select columns by column header name(s)`
> 11. {% icon param-text %} *"Select target column(s)"*: `target`
> 12. `Execute` the classifier to train
{: .hands_on}


## Predict using a trained model
The previous step produces a model file of type `zip`. Rename this file to `model.zip` by using `edit` dataset property. The trained model is used to predict the categories of each row in `breast-w_test` dataset.

> ### {% icon hands_on %} Hands-on: Predict using a trained model
> 
> **SVM Classifier (Support vector machine)** {% icon tool %} with the following parameters
> 
> 1. {% icon param-select %} *"Select a Classification Task"*: `Load a model and predict`
> 2. {% icon param-file %} *"Models"*: `model.zip`
> 3. {% icon param-file %} *"Data (tabular)"*: `breast-w_test`
> 4. {% icon param-check %} *"Does the dataset contain header"*: `Yes`
> 5. {% icon param-select %} *"Select the type of prediction"*: `Predict class labels`
> 6. `Execute` to predict categories
{: .hands_on}


## Visualize the predictions
The previous step produces a predicted file of type `tabular`. Rename this file to `predicted_labels` by using `edit` dataset property. This tool will create three output files, one each for confusion matrix, precision, recall and f1-score and roc curves. These files are zipped.

> ### {% icon hands_on %} Hands-on: Visualize the predictions
> 
> **Plot confusion matrix, precision, recall and ROC and AUC curves** {% icon tool %} with the following parameters
> 
> 1. {% icon param-file %} *"Select input data file"*: `breast-w_targets`
> 2. {% icon param-file %} *"Select predicted data file"*: `predicted_labels`
> 3. {% icon param-file %} *"Select trained model"*: `model.zip`
> 4. `Execute` to create visualizations
>
>    ![Dataset](images/confusion_matrix.png)
>    ![Dataset](images/precision_recall_f1.png)
>    ![Dataset](images/roc.png)
{: .hands_on}


# Regression

Regression is also a supervised learning task where the target is a real number (continuous) instead of discreet like classification. The algorithms which are used for regression tasks are called regressors. A regressor learns mapping between features of a data row and its target value. Inherently, it tries to fit a curve for the targets. This curve can be linear (straight line curve) or non-linear. 


## Data upload

The datasets required for this tutorial contain 9 features of breast cancer which include the thickness of clump, cell-size, cell-shape and so on ([more information](https://github.com/EpistasisLab/penn-ml-benchmarks/tree/master/datasets/classification/breast-w)). In addition to these features, the training dataset contains one more column as `target`. It has a binary value (0 or 1) for each row. `0` indicates no breast cancer and `1` indicates breast cancer. The test dataset does not contain the `target` column. The third dataset contains the all the samples from test dataset but also the `target` column which would be needed to create plot showing comparison between actual and predicted targets.


> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the following datasets and choose the type of data as `tabular`.
> 
>    ```
>    https://zenodo.org/record/1401230/files/breast-w_train.tsv
>    https://zenodo.org/record/1401230/files/breast-w_test.tsv
>    https://zenodo.org/record/1401230/files/breast-w_targets.tsv
>    ```
> 
>    {% include snippets/import_via_link.md %}
>
> 3. Rename datasets to `breast-w_train`, `breast-w_test` and `breast-w_targets`.
>
>    {% include snippets/rename_dataset.md %}
>
{: .hands_on}


## Train a regressor
In this step, `Gradient boosting` regressor is used for the regression task. The last column of this dataset assigns a category for each row. The classifier learns a mapping between each row and its category. This mapping is called a trained model. It is used to predict the categories of unseen data (`breast-w_test`).

> ### {% icon hands_on %} Hands-on: Train a classifier
> 
> **Gradient boosting regressor** {% icon tool %} with the following parameters
> 1. {% icon param-select %} *"Select a Classification Task"*: `Train a model`
> 2. {% icon param-select %} *"Select an ensemble method"*: `Gradient Boosting Regressor`
> 3. {% icon param-select %} *"Select input type"*: `tabular data`
> 4. {% icon param-file %} *"Training samples dataset"*: `breast-w_train`
> 5. {% icon param-check %} *"Does the dataset contain header"*: `Yes`
> 6. {% icon param-select %} *"Choose how to select data by column"*: `All columns but by column header name(s)`
> 7. {% icon param-text %} *"Type header name(s)"*: `target`
> 8. {% icon param-file %} *"Dataset containing class labels"*: `breast-w_train`
> 9. {% icon param-check %} *"Does the dataset contain header"*: `Yes`
> 10. {% icon param-select %} *"Choose how to select data by column"*: `Select columns by column header name(s)`
> 11. {% icon param-text %} *"Select target column(s)"*: `target`
> 12. `Execute` the classifier to train
{: .hands_on}


## Predict using a trained model
The previous step produces a model file of type `zip`. Rename this file to `model.zip` by using `edit` dataset property. The trained model is used to predict the categories of each row in `breast-w_test` dataset.

> ### {% icon hands_on %} Hands-on: Predict using a trained model
> 
> **SVM Classifier (Support vector machine)** {% icon tool %} with the following parameters
> 
> 1. {% icon param-select %} *"Select a Classification Task"*: `Load a model and predict`
> 2. {% icon param-file %} *"Models"*: `model.zip`
> 3. {% icon param-file %} *"Data (tabular)"*: `breast-w_test`
> 4. {% icon param-check %} *"Does the dataset contain header"*: `Yes`
> 5. {% icon param-select %} *"Select the type of prediction"*: `Predict class labels`
> 6. `Execute` to predict categories
{: .hands_on}


## Visualize the predictions
The previous step produces a predicted file of type `tabular`. Rename this file to `predicted_labels` by using `edit` dataset property. This tool will create three output files, one each for confusion matrix, precision, recall and f1-score and roc curves. These files are zipped.

> ### {% icon hands_on %} Hands-on: Visualize the predictions
> 
> **Plot confusion matrix, precision, recall and ROC and AUC curves** {% icon tool %} with the following parameters
> 
> 1. {% icon param-file %} *"Select input data file"*: `breast-w_targets`
> 2. {% icon param-file %} *"Select predicted data file"*: `predicted_labels`
> 3. {% icon param-file %} *"Select trained model"*: `model.zip`
> 4. `Execute` to create visualizations
>
>    ![Dataset](images/confusion_matrix.png)
>    ![Dataset](images/precision_recall_f1.png)
>    ![Dataset](images/roc.png)
{: .hands_on}



