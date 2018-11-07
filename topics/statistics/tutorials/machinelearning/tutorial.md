---
layout: tutorial_hands_on

title: "Basics of machine learning"
zenodo_link: https://zenodo.org/record/1468039#.W8zyxBRoSAo
questions:
  - "What is machine learning?"
  - "Why is it useful?"
  - "What are its different approaches?"
objectives:
  - "Provide the basics of machine learning and its variants."
  - "Learn how to do classification using the training and test data."
  - "Learn how to use Galaxy's machine learning tools."
requirements:
time_estimation: "30M"
key_points:
  - "Machine learning algorithms learn features from data."
  - "It is used for multiple tasks like classification, regression, clustering and so on."
  - "Multiple learning tasks can be performed using Galaxy's machine learning tools."
  - "For the classification and regression tasks, data is divided into training and test sets."
  - "Each sample/record in the training data has a category/class/label."
  - "A machine learning algorithm learns features from the training data and do predictions on the test data."

contributors:
  - anuprulez

---

# Introduction
{:.no_toc}

[Machine learning](https://en.wikipedia.org/wiki/Machine_learning) uses the techniques from statistics, mathematics and computer science to make computer programs learn from data. It is one of the most popular fields of computer science and finds applications in multiple streams of data analysis like [classification](https://en.wikipedia.org/wiki/Statistical_classification), [regression](https://en.wikipedia.org/wiki/Regression_analysis), [clustering](https://en.wikipedia.org/wiki/Cluster_analysis), [dimensionality reduction](https://en.wikipedia.org/wiki/Dimensionality_reduction), [density estimation](https://en.wikipedia.org/wiki/Density_estimation) and many more. Some real-life applications are spam filtering, medical diagnosis, autonomous driving, recommendation systems, facial recognition, stock prices prediction and many more. The following image shows a basic flow of any machine learning task. A user has data and it is given to a machine learning algorithm for analysis.

>    ![data](images/ml_basics.png "Flow of a machine learning task.")

There are multiple ways in which machine learning can be used to perform data analysis. They depend on the nature of data and the kind of data analysis. The following image shows the most popular ones. In [supervised learning](https://en.wikipedia.org/wiki/Supervised_learning) techniques, the categories of data records are known beforehand. But in [unsupervised learning](https://en.wikipedia.org/wiki/Unsupervised_learning), the categories of data records are not known.

>    ![data](images/variants_ml.png "Different types of machine learning.")

In general, machine learning can be used in multiple real-life tasks by using applying its variants as depicted in the following image.

>    ![data](images/usage_ml.png "Real-life usage of machine learning.")

The following image shows how a classification task is performed. The complete data is divided into training and test sets. The training set is used by a classifier to learn features. It results in a trained model and its robustness (of learning) is evaluated using the test set (unseen by the classifier during the training).

>    ![data](images/prediction.png "Supervised learning.")

This tutorial shows how to use a machine learning module implemented as a Galaxy tool. The data used in this tutorial is available at [Zenodo](https://zenodo.org/record/1468039#.W8zyxBRoSAo).

> ### Agenda
>
> Performing a machine learning task (classification) using a tool involves the following steps:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Data upload

The datasets required for this tutorial contain 9 features of breast cancer which include the thickness of clump, cell-size, cell-shape and so on ([more information](https://github.com/EpistasisLab/penn-ml-benchmarks/tree/master/datasets/classification/breast-w)). In addition to these features, the training dataset contains one more column as `target`. It has a binary value (0 or 1) for each row. `0` indicates no breast cancer and `1` indicates breast cancer. The test dataset does not contain the `target` column.


> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the following datasets and choose the type of data as `tabular`.
> 
>    ```
>    https://zenodo.org/record/1401230/files/breast-w_train.tsv
>    https://zenodo.org/record/1401230/files/breast-w_test.tsv
>    ```
> 
>    {% include snippets/import_via_link.md %}
>
> 3. Rename datasets to `breast-w_train` and `breast-w_test`.
>
>    {% include snippets/rename_dataset.md %}
>
> 4. The datasets should look like these:
>
>
>    ![data](images/train_data.png "Training data (breast-w_train) with targets (9 features and one target).")
>
>
>    ![data](images/test_data.png "Test data (breast-w_test) (9 features and no target).")
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


## See predictions
The last column of the predicted dataset shows the category of each row. A row either gets `0` (no breast cancer) or `1` (breast cancer) as its predicted category. 

> ### {% icon hands_on %} Hands-on: See the predicted column
> 1. Click on `view data` link of the dataset created after executing the previous step.
> 2. The last column of the `tabular` data shows the predicted category (`target`) for each row.
> 
{: .hands_on}


> ### {% icon tip %} Additional resources:
>
> Read more about **machine learning using scikit-learn** [here](http://scikit-learn.org/stable/).
{:.tip}
