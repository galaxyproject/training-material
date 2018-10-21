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
  - "Learn how to make predictions using the training and test data."
  - "Classify data using a Galaxy workflow."
requirements:
time_estimation: "30M"
key_points:
  - "Machine learning algorithms learn features from data."
  - "It is used for multiple tasks like classification, regression, clustering and so on."
  - "For the classification task, data is divided into training and test sets."
  - "Each data sample in training and test sets has a category/class."
  - "Many learning tasks can be performed on datasets using Galaxy tools for machine learning."
contributors:
  - anuprulez

---

# Introduction
{:.no_toc}

Machine learning uses the techniques from statistics, mathematics and computer science to make computer programs learn from data. It is one of the most popular fields of computer science and finds applications in multiple streams of data analysis like classification, regression, clustering, dimensionality reduction, density estimation and many more. Some real-life applications are spam filtering, medical diagnosis, autonomous driving, recommendation systems, facial recognition, stock prices prediction and many more. The following image shows a basic flow of any machine learning task. A user has data and it is given to a machine learning algorithm for analysis.

![Dataset](images/ml_basics.png)

There are multiple ways in which machine learning can be used to perform data analysis. They depend on the nature of data and the kind of data analysis. The following image shows the most popular ones.

![Dataset](images/variants_ml.png)

The following image shows how a classification task is performed. The complete data is divided into training and test sets. The training set is used by a classifier to learn features. It results in a trained model and it is evaluated using the test set (unseen by the classifier during the training).

![Dataset](images/prediction.png)

This tutorial shows how to use machine learning modules implemented as Galaxy tools. Few machine learning tools are present in the tools collection under the header "statistics". These tools can be used to create workflows to perform a machine learning task.

The data used in this tutorial is available at [Zenodo](https://doi.org/10.5281/zenodo.1404173).

> ### Agenda
>
> Performing a machine learning task (classification) using a tool which involves the following steps:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Upload training data

This dataset contains 9 features of breast cancer. These features include the thickness of clump, cell-size, cell-shape and so on ([more information](https://github.com/EpistasisLab/penn-ml-benchmarks/tree/master/datasets/classification/breast-w)). The last column in the tabular data contains binary value (0 or 1) where a row (a set of features) is classified into breast cancer or not.

> ### {% icon hands_on %} Hands-on: Data upload
> 1. Download and import the following dataset in a new Galaxy history. Choose the type of data as `tabular`.
     - `breast-w_test.tsv`
>
>    ```
>    https://zenodo.org/record/1468039/files/breast-w_train.tsv
>    ```
>
>    {% include snippets/import_via_link.md %}
>
>    ![Dataset](images/train_data.png)
>
{: .hands_on}


# Choose a classifier and update its parameters

> ### {% icon hands_on %} Hands-on: **SVM Classifier (Support vector machine)** {% icon tool %}
> 1. Choose **SVM Classifier (Support vector machine)** {% icon tool %} tool.
> 2. Choose `train a model` as a classification task and choose the type of classifier as `linear support vector classification`.
> 4. Choose the type of data as `tabular`.
> 
{: .hands_on}


# Execute the classifier

> ### {% icon hands_on %} Hands-on: Execute the classifier
> 1. Click on `execute` button to execute the classifier.
> 2. See the resulting model (new dataset) in the history.
> 
{: .hands_on}


# Upload test data

> ### {% icon hands_on %} Hands-on: Upload test data
> 1. Upload test dataset `breast-w_test.tsv` from [Zenodo](https://zenodo.org/record/1468039/files/breast-w_test.tsv).
> 2. This dataset does not contain the target column. This column will be predicted in the following step using the trained model.
> 
> ![Dataset](images/test_data.png)
{: .hands_on}


# Choose trained model to predict

> ### {% icon hands_on %} Hands-on: Choose trained model to predict
> 1. Choose **SVM Classifier (Support vector machine)** {% icon tool %} tool.
> 2. Set the classification task as `load a model and predict`. 
> 3. Choose the trained model. 
> 4. Choose the test data from previous step.
> 5. Execute the tool.
> 
{: .hands_on}


# See predictions

> ### {% icon hands_on %} Hands-on: See the predicted columns
> 1. Click on `view data` link of the dataset created after executing the previous step.
> 2. The last column of the `tabular` data shows the predicted category (target) for each row.
> 
{: .hands_on}


> ### {% icon tip %} Additional resources:
>
> Read more about **machine learning using scikit-learn** [here](http://scikit-learn.org/stable/).
{:.tip}
