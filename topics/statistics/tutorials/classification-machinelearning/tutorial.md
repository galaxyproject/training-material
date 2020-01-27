---
layout: tutorial_hands_on

title: 'Classification in Machine Learning'
zenodo_link: https://zenodo.org/record/...
questions:
- What is classification and how we can use classification techniques?
objectives:
- Learn classification background
- Apply classification based machine learning algorithms
- Learn what a quantitative structure-analysis relationship (QSAR) model is and how it can be constructed in Galaxy
- Learn how visualizations can be used to analyze the classification results
key_points:
- Classification is a supervised approach in machine learning.
- For classification tasks, data is divided into training and test sets.
- Using classification, the samples are learned using the training set and predicted using the test set.
- For each classification algorithm, it parameters should be optimised based on the dataset.
- Machine learning algorithms can be applied to chemical datasets to predict important properties.
time_estimation: 2H
contributors:
- khanteymoori
- anuprulez
- simonbray
---

# Introduction
{:.no_toc}

In this tutorial you will learn how to apply Galaxy tools to solving [classification](ttps://en.wikipedia.org/wiki/Statistical_classification) problems. First, we will introduce classification briefly, and then examine the logistic regression which is the linear classifier. Next, we will discuss the nearest neighbor classifier, which is a simple but nonlinear classifier. Then advanced classifiers, such as support vector machines, random forest and ensemble classifiers will be introduced and applied. Furthermore, we will show how to visualize the results in each step. 
Finally, we will discuss how to train the classifiers by finding the values of their parameters that minimize a cost function. We will work through a real problem to learn how the classifiers and learning algorithms work.
Classification is a [supervised learning](https://en.wikipedia.org/wiki/Supervised_learning) method in machine learning and the algorithm which is used for this learning task is called a classifier. In this tutorial we will build a classifier which can predict whether a chemical substance is biodegradable or not. Substances which degrade quickly are preferable to those which degrade slowly, as they do not accumulate and pose a risk to the environment. Therefore, it is useful to be able to predict easily in advance whether a substance is biodegradable prior to production and usage in consumer products.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Classification

Classification is the process of assigning every object from a collection to exactly one class from a known set of classes by learning a decision boundary in a dataset. This dataset is called a training dataset and contains samples and desired class for each sample. The training dataset contains "features" as columns and a mapping between these features and the class label is learned for each sample. The performance of mapping is evaluated using a test dataset which is separate from training dataset. The test dataset contains only the feature columns and not the class column. The class column is predicted using the mapping learned on the training dataset. Examples of classification task is assigning a patient (the object) to a group of healthy or ill (the classes) people on the basis of his or her medical record. In this tutorial, we will use a classifier to train a model using a training dataset, predict the targets for test dataset and visualize the results using plots.

![classification](images/classification.png "Classification of samples belonging to different classes.")

In figure [1](#figure-1), the line is a boundary which separates a class from another class (for example from tumor to no tumor). The task of a classifier is to learn this boundary, which can be used to classify or categorize an unseen/new sample. The line is the decision boundary. There are different ways to learn this decision boundary. If the dataset is linearly separable, linear classifiers can produce good classification results. But, when the dataset is complex and requires non-linear decision boundaries, more powerful classifiers like `support vector machine` or `ensemble` based classifiers may prove to be beneficial. 

The Data Classification process includes two steps:
1. Building the Classifier or Model: This step is the learning step or the learning phase and in this step the classification algorithms build the classifier. The classifier is built from the training set made up of database samples and their associated class labels. Each sample that constitutes the training set is referred to as a class. 

2. Using Classifier for Classification: In this step, the classifier is used for classification. Here the test data is used to estimate the accuracy of classification rules. The classification rules can be applied to the new data samples if the accuracy is considered acceptable. 


# Quantitative Structure - Activity Relationship biodegradation

The classification problem we will study in this tutorial is related to biodegradation. Chemical substances which decay slowly will accumulate over time, which poses a threat to the environment. Therefore, it is useful to be able to predict in advance whether a substance will break down quickly or not.

Quantitative structure-activity relationip (QSAR) and quantitative structure-property relationship (QSPR) models attempt to predict the activity or property of chemicals based on their chemical structure. To achieve this, a database of compounds is collected for which the property of interest is known. For each compound, molecular descriptors are collected which describe the structure (for example: molecular weight, number of nitrogen atoms, number of carbon-carbon double bonds). Using these descriptors, a model is constructed which is capable of predicting the property of interest for a new, unknown molecule. In this tutorial we will use a database assembled from experimental data of the Japanese Ministry of International Trade and Industry to create a classification model for biodegradation. We then will be able to use this model to classify new molecules into one of two classes: biodegradable or non-biodegradable.

As a benchmark, we will use the [dataset](https://pubs.acs.org/doi/10.1021/ci4000213) assembled by Mansouri et al. using data from the National Institute of Technology and Evaluation of Japan. This database contains 1055 moelecules, together with precalculated molecular descriptors.

In this tutorial, we will apply a couple of ([scikit-learn](https://scikit-learn.org/stable/)) machine learning tools to dataset created by Mansouri et al. to predict whether a molecule is biodegradable or not.
In the following part, we will perform classification on the biodegradability dataset using a linear classifier and then will analyze the results with plots.

## Get train and test datasets

We have two datasets available; the training dataset contains 837 molecules, while the test dataset contains 218 molecules.

Let's begin by uploading the necessary datasets.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo](https://zenodo.org/record/....)
>
>    ```
>    https://zenodo.org/record/.../files/train_rows.csv
>    https://zenodo.org/record/.../files/test_rows_labels.csv
>    https://zenodo.org/record/.../files/test_rows.csv
>    ```
>
>    {% include snippets/import_via_link.md %}
>
> 3. Rename the datasets as `train_rows`, `test_rows_labels` and `test_rows` respectively.
>
>    {% include snippets/rename_dataset.md %}
>
> 4. Check that the datatype of all the three datasets is `tabular`.
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
{: .hands_on}

The `train_rows` contains a column `Class` which is the class label or target. We will evaluate our model on `test_rows` and compare the predicted class with the true class value in `test_rows_labels`
{: .comment}

> ### {% icon details %} Preparing the data for Classification
>
> Preparing the data involves these following major activities: 
> 1. Data Cleaning: involves removing the noise and treatment of missing values. The noise is removed by applying noise filtering techniques and the problem of missing values is solved by replacing a missing value with different techniques. 
> 2. Relevance Analysis: Database may also have the irrelevant attributes. Correlation analysis is used to know whether any two given attributes are related.
> 3. Normalization: The data is transformed using normalization. Normalization involves scaling all values for given attribute in order to make them fall within a small specified range. Normalization is used when in the learning step, the neural networks or the methods involving measurements are used.
>
{: .details}

# Learn the logistic regression classifier

At the first step, to learn the mapping between several features and the classes, we will apply linear classifier. It learns features from training dataset and maps all the rows to their respective class. The process of mapping gives a trained model. Logistic regression is an instance of supervised classification in which we know the correct label of the class for each sample and the algorithm estimate of the true class. We want to learn parameters (weight and bias for the line) that make estimated calss for each training observation as close as possible to the true class label. This requires two components, the first is a metric for how close the current class label is to the true label. Rather than measure similarity, we usually talk about the opposite of this, the distance between the classifier output and the desired output, and we call this distance, the loss function or the cost function. 
The second thing we need is an optimization algorithm for iteratively updating the weights so as to minimize this loss function. The standard algorithm for this is gradient descent. So, the dataset is divided into two parts - training and test sets. The training set is used to train a classifier and the test set is used to evaluate the performance of the trained model.

> ### {% icon hands_on %} Hands-on: Train logistic regression classifier
>
> 1. **Generalized linear models for classification and regression** {% icon tool %} with the following parameters to train the regressor:
>    - *"Select a Classification Task"*: `Train a model`
>       - *"Select a linear method"*: `Logistic Regression model`
>          - *"Select input type"*: `tabular data`
>             - {% icon param-file %} *"Training samples dataset"*: `train_rows.csv`
>             - {% icon param-check %} *"Does the dataset contain header"*: `Yes`
>             - {% icon param-select %} *"Choose how to select data by column"*: `All columns EXCLUDING some by column header name(s)`
>                - {% icon param-text %} *"Type header name(s)"*: `Class`
>             - {% icon param-file %} *"Dataset containing class labels"*: `train_rows.csv`
>             - {% icon param-check %} *"Does the dataset contain header"*: `Yes`
>             - {% icon param-select %} *"Choose how to select data by column"*: `Select columns by column header name(s)`
>                - {% icon param-text %} *"Select target column(s)"*: `Class`
> 2. Rename the generated file to `LogisticRegression_model`
{: .hands_on}

> ### {% icon question %} Question
>
> What is learned by the logistic regression model?
>
> > ### {% icon solution %} Solution
> >
> > In the logistic regressoion model, the coefficients of the best classifier line, is learned.
> > 
> {: .solution}
>
{: .question}

## Predict class using test dataset

After learning on the training dataset, we should evaluate the performance on the test dataset to know whether the learning algorithm learned a good classifier from the training dataset or not. This classifier is used to predict a new sample and a similar accuracy is expected. 
Now, we will predict class in the test dataset using this classifier in order to see if the classifier has learned important features which can generalize on a new dataset. The test dataset (`test_rows`) contains the same number of features but does not contain the `Class` column. This is predicted using the trained classifier.


> ### {% icon hands_on %} Hands-on: Predict class using the logistic regression classifier
>
> 1. **Generalized linear models for classification and regression** {% icon tool %} with the following parameters to predict targets of test dataset using the trained model:
>    - *"Select a Classification Task"*: `Load a model and predict`
>       - {% icon param-file %} *"Models"*: `LogisticRegression_model`
>       - {% icon param-file %} *"Data (tabular)"*: `test_rows.csv`
>       - {% icon param-check %} *"Does the dataset contain header"*: `Yes`
>       - {% icon param-select %} *"Select the type of prediction"*: `Predict class labels`
> 2. Rename the generated file to `classification_linear`
{: .hands_on}

## Visualize the logistic regression classification results

We will evaluate the classification by comparing them to the expected classes. In the previous step, we classify the test dataset (`classification_linear`). We have one more dataset (`test_rows_labels`) which contains the true class label of the test set. Using the true and predicted class labels in the test set, we will verify the performance by analyzing the plots. As you can see, `classification_linear` has no header, so first we should remove the header from `test_rows_labels` to compare. 

> ### {% icon hands_on %} Hands-on: Remove the header
>
> 1. **Remove beginning of a file** {% icon tool %} with the following parameters:
>       - {% icon param-file %} *"Remove first"*: `1`
>       - {% icon param-file %} *"from"*: `test_rows_labels.csv`
> 2. Rename the generated file to `test_rows_labels_without_header.csv`
{: .hands_on}


Now we visualize and analyze the classification using "Plot confusion matrix, precision, recall and ROC and AUC curves" tool in galaxy.

> ### {% icon hands_on %} Hands-on: Check and visualize the classification
> 1. **Plot confusion matrix, precision, recall and ROC and AUC curves** {% icon tool %} with the following parameters to visualize the classification:
>    - {% icon param-file %} *"Select input data file"*: `test_rows_labels_without_header.csv`
>    - {% icon param-file %} *"Select predicted data file"*: `classification_linear`
>    - {% icon param-file %} *"Select trained model"*: `LogisticRegression_model`
{: .hands_on}

The visualization tool creates the following plots:

1. [Confusion matrix](https://en.wikipedia.org/wiki/Confusion_matrix): Confusion matrix summarizes the classification performance of a classifier with respect to test data. It is a two-dimensional matrix, the horizontal axis (x-axis) shows the predicted labels and the vertical axis (y-axis) shows the true labels. Each rectangular box shows a count of samples falling into the four output combinations (true class, predicted class) - (1, 0), (1, 1), (0, 1) and (0, 0). In Figure 2, confusion matrix of the predictions is a heatmap. For a good prediction, the diagonal running from top-left to bottom-right should contain less number of samples, because it shows the counts of incorrectly predicted samples. Hovering over each box in Galaxy shows the true and predicted class labels and the count of samples.

    ![confusion_matrix](images/confusion_matrix_linear.png "Confusion matrix for the logistic regression classifier. ")

2. [Precision, recall and F1 score](https://en.wikipedia.org/wiki/Precision_and_recall): Precision, recall and F1 score. These scores determine the robustness of classification. It is important to analyze the plot for any classification task to verify the accuracy across different classes which provides more information about the balanced or imbalanced accuracy across multiple classes present in the dataset.

    ![prf1_scores](images/precision_recall_linear.png "Precision, recall and F1 score for the logistic regression classifier.")

3. [Receiver operator characteristics (ROC) and area under ROC (AUC)](https://towardsdatascience.com/understanding-auc-roc-curve-68b2303cc9c5): Receiver operator characteristics (ROC) and area under ROC (AUC). The ROC curve is shown in blue. For a good prediction, it should be more towards the top-left of this plot. For a bad prediction, it is close to the orange line (y = x).

    ![roc_scores](images/roc_linear.png "Receiver operator characteristics (ROC) and area under ROC (AUC) for the logistic regression classifier.")

These plots are important to visualize the quality of classifier and the true and predicted classes.


> ### {% icon question %} Question
>
> Inspect the plots. What can you say about the classification?
>
> > ### {% icon solution %} Solution
> >
> > Figures 2,3 and 4 show that the classification is acceptable, but as you will see in the next steps, the reults can be improved. 
> >
> {: .solution}
{: .question}

# Using k-nearest neighbor classifier (KNN)

At the second step, we will use k-nearest neighbor classifier. In the k-nearst neighbor classifier, a sample is classified by a majority vote of its neighbors.  The sample is assigned to the class which is most common among it's k nearest neighbors.  k is a positive integer and typically it is small. For example, if k = 1, then the sample is simply assigned to the class of that single nearest neighbor. Surprisingly, when the number of data points is large, this classifieris not that bad. Choosing the best value of k is very important. If k is too small, the classifier will be sensitive to noise points and If k is too large, neighborhood may include points from other classes and cause errors. To select the k that’s right for your data, we recoomend that run the KNN algorithm several times with different values of k and choose the k that reduces the number of errors. We encounter while maintaining the algorithm’s ability to accurately make predictions when it’s given data it hasn’t seen before.

> ### {% icon question %} Question
>
> What are advantages and disadvantages about this model?
>
> > ### {% icon solution %} Solution
> > Advantages:
> > - It is very simple algorithm to understand and interpret.
> >
> > - It is very useful for nonlinear data because there is no assumption about data in this algorithm.
> >
> > - It is a versatile algorithm as we can use it for classification as well as regression.
> >
> > - It has relatively high accuracy but there are much better supervised learning models than KNN.
> >
> > - It works very well in low dimensions for complex decision surfaces.
> >
> > Disadvantages:
> >
> > - Classification is slow, because it stores all the training data.
> >
> > - High memory storage required as compared to other supervised learning algorithms.
> >
> > - Prediction is slow in case of big training samples.
> >
> > - It is very sensitive to the scale of data as well as irrelevant features.
> >
> > - It suffers a lot from the curse-of-dimensionality.
> >
> {: .solution}
{: .question}

> ### {% icon hands_on %} Hands-on: Train k-nearest neighbor classifier
>
> 1. **Nearest Neighbors Classification** {% icon tool %} with the following parameters to train the regressor:
>    - *"Select a Classification Task"*: `Train a model`
>       - *"Classifier type"*: `Nearest Neighbors`
>          - *"Select input type"*: `tabular data`
>             - {% icon param-file %} *"Training samples dataset"*: `train_rows.csv`
>             - {% icon param-check %} *"Does the dataset contain header"*: `Yes`
>             - {% icon param-select %} *"Choose how to select data by column"*: `All columns EXCLUDING some by column header name(s)`
>                - {% icon param-text %} *"Type header name(s)"*: `Class`
>             - {% icon param-file %} *"Dataset containing class labels"*: `train_rows.csv`
>             - {% icon param-check %} *"Does the dataset contain header"*: `Yes`
>             - {% icon param-select %} *"Choose how to select data by column"*: `Select columns by column header name(s)`
>                - {% icon param-text %} *"Select target column(s)"*: `Class`
>             - {% icon param-select %} *"Neighbor selection method"*: `k-nearest neighbors`
> 2. Rename the generated file to `NearestNeighbors_model`
{: .hands_on}


Now, we should evaluate the performance on the test dataset to know whether the KNN classifier is a good model from the training dataset or not. 

> ### {% icon hands_on %} Hands-on: Predict class using the k-nearest neighbor classifier
>
> 1. **Nearest Neighbors Classification** {% icon tool %} with the following parameters to predict targets of test dataset using the trained model:
>    - *"Select a Classification Task"*: `Load a model and predict`
>       - {% icon param-file %} *"Models"*: `NearestNeighbors_model`
>       - {% icon param-file %} *"Data (tabular)"*: `test_rows.csv`
>       - {% icon param-check %} *"Does the dataset contain header"*: `Yes`
>       - {% icon param-select %} *"Select the type of prediction"*: `Predict class labels`
> 2. Rename the generated file to `classification_NearestNeighbors`
{: .hands_on}


Now we visualize and analyze the classification:

> ### {% icon hands_on %} Hands-on: Check and visualize the classification
> 1. **Plot confusion matrix, precision, recall and ROC and AUC curves** {% icon tool %} with the following parameters to visualize the classification:
>    - {% icon param-file %} *"Select input data file"*: `test_rows_labels_without_header.csv`
>    - {% icon param-file %} *"Select predicted data file"*: `classification_NearestNeighbors`
>    - {% icon param-file %} *"Select trained model"*: `NearestNeighbors_model`
{: .hands_on}

The visualization tool creates the Confusion matrix, Precision, recall and F1 score, Receiver operator characteristics (ROC) and area under ROC (AUC) as follows:
1.
    ![confusion_matrix](images/confusion_matrix_NN.png "Confusion matrix for the k-nearest neighbor classifier.")
2.
    ![prf1_scores](images/precision_recall_NN.png "Precision, recall and F1 score for the k-nearest neighbor classifier.")
3.
    ![roc_scores](images/roc_linear.png "Receiver operator characteristics (ROC) and area under ROC (AUC) for the k-nearest neighbor classifier.")


# Support Vector Machines (SVMs)

Support Vector Machines(SVMs) have been extensively researched in the data mining and machine learning communities for the last decade and actively applied to applications in various domains such as bioinformatics. SVM is a generalization of a classifier called maximal margin classifier and is introduced as a binary classifier intended to separate two classes when obtaining the optimal hyperplane and decision boundary. SVMs are based on the assumption that the input data can be linearly separable in a geometric space. The maximal margin classifier is simple, but it cannot be applied to the majority of datasets, since the classes must be separated by a linear boundary and this is often not the case when working with real word data.. That is why the support vector classifier was introduced as an extension of the maximal margin classifier, which can be applied in a broader range of cases.
To solve this problem SVM using kernel functions to map the input to a high dimension feature space, i.e hyperplane, where a linear decision boundary is constructed in such a manner that the boundary maximises the margin between two classes. The kernel approach is simply an efficient computational approach for accommodating a non-linear boundary between classes.
Without going into technical details, a kernel is a function that quantifies the similarity of two observations.
Two special properties of SVMs are that SVMs achieve (1) high generalization by maximizing the margin and (2) support an efficient learning of nonlinear functions by
kernel trick. In the next step, we will build a SVM classifier with our data. 

> ### {% icon hands_on %} Hands-on: Train a SVM classifier
>
> 1. **Support vector machines (SVMs) for classification** {% icon tool %} with the following parameters to train the regressor:
>    - *"Select a Classification Task"*: `Train a model`
>       - *"Select a linear method"*: `Linear Support Vector Classification`
>          - *"Select input type"*: `tabular data`
>             - {% icon param-file %} *"Training samples dataset"*: `train_rows.csv`
>             - {% icon param-check %} *"Does the dataset contain header"*: `Yes`
>             - {% icon param-select %} *"Choose how to select data by column"*: `All columns EXCLUDING some by column header name(s)`
>                - {% icon param-text %} *"Type header name(s)"*: `Class`
>             - {% icon param-file %} *"Dataset containing class labels"*: `train_rows.csv`
>             - {% icon param-check %} *"Does the dataset contain header"*: `Yes`
>             - {% icon param-select %} *"Choose how to select data by column"*: `Select columns by column header name(s)`
>                - {% icon param-text %} *"Select target column(s)"*: `Class`
> 2. Rename the generated file to `SVM_model`
{: .hands_on}

> ### {% icon question %} Question
>
> What is learned by the support vector machines?
>
> > ### {% icon solution %} Solution
> >
> > The coefficients of the line with the maximal margin in the kernel space is learned in training phase.
> > 
> {: .solution}
>
{: .question}


No we will evaluate the performance of the SVM classifier:

> ### {% icon hands_on %} Hands-on: Predict class SVM classifier
>
> 1. **Support vector machines (SVMs) for classification** {% icon tool %} with the following parameters to predict targets of test dataset using the trained model:
>    - *"Select a Classification Task"*: `Load a model and predict`
>       - {% icon param-file %} *"Models"*: `SVM_model`
>       - {% icon param-file %} *"Data (tabular)"*: `test_rows.csv`
>       - {% icon param-check %} *"Does the dataset contain header"*: `Yes`
>       - {% icon param-select %} *"Select the type of prediction"*: `Predict class labels`
> 2. Rename the generated file to `classification_SVM`
{: .hands_on}


Now lets visualize the resluts:

> ### {% icon hands_on %} Hands-on: Check and visualize the classification
> 1. **Plot confusion matrix, precision, recall and ROC and AUC curves** {% icon tool %} with the following parameters to visualize the classification:
>    - {% icon param-file %} *"Select input data file"*: `test_rows_labels_without_header.csv`
>    - {% icon param-file %} *"Select predicted data file"*: `classification_SVM`
>    - {% icon param-file %} *"Select trained model"*: `SVM_model`
{: .hands_on}

The visualization tool creates the following plots:
1.
    ![confusion_matrix](images/confusion_matrix_svm.png "Precision, recall and F1 score for the SVM classifier.")
2.
    ![prf1_scores](images/precision_recall_svm.png "Precision, recall and F1 score for the SVM classifier.")
3. 
    ![roc_scores](images/roc_svm.png "Receiver operator characteristics (ROC) and area under ROC (AUC) for the SVM classifier.")


# Using random forest for classification

Random forest is an ensemble of decision trees, and usually trained with the “bagging” method. [*Ensemble*](https://scikit-learn.org/stable/modules/ensemble.html#ensemble) method uses multiple learning models internally for better predictions and the general idea of the bagging method is that a combination of learning models increases the overall result. it uses multiple decision tree regressors internally and predicts by taking the collective performances of the predictions by multiple decision trees. It has a good predictive power and is robust to the outliers. It creates an ensemble of weak learners (decision trees) and iteratively minimizes error. One big advantage of random forest is that it can be used for both classification and regression problems. The main idea behind the random forest is adding additional randomness to the model, while growing the trees and instead of searching for the most important feature while splitting a node, it searches for the best feature among a random subset of features. This results better model because of wide diversity. Generally, the more trees in the forest the more robust the forest looks like. In the same way in the random forest classifier, the higher the number of trees in the forest gives the high accuracy results.the more trees in the forest the more robust the forest looks like. In the same way in the random forest classifier, the higher the number of trees in the forest gives the high accuracy results.
There are two stages in Random Forest algorithm, one is random forest creation, the other is to make a prediction from the random forest classifier created in the first stage.

> ### {% icon hands_on %} Hands-on: Train random forest
>
> 1. **Ensemble methods for classification and regression** {% icon tool %} with the following parameters to train the regressor:
>    - *"Select a Classification Task"*: `Train a model`
>       - *"Select an ensemble method"*: `Random forest classifier`
>          - *"Select input type"*: `tabular data`
>             - {% icon param-file %} *"Training samples dataset"*: `train_rows.csv`
>             - {% icon param-check %} *"Does the dataset contain header"*: `Yes`
>             - {% icon param-select %} *"Choose how to select data by column"*: `All columns EXCLUDING some by column header name(s)`
>                - {% icon param-text %} *"Type header name(s)"*: `Class`
>             - {% icon param-file %} *"Dataset containing class labels"*: `train_rows.csv`
>             - {% icon param-check %} *"Does the dataset contain header"*: `Yes`
>             - {% icon param-select %} *"Choose how to select data by column"*: `Select columns by column header name(s)`
>                - {% icon param-text %} *"Select target column(s)"*: `Class`
> 2. Rename the generated file to `RandomForestClassifier`
{: .hands_on}

> ### {% icon question %} Question
>
> What are the advantages of random forest classifier compare with classifiers?
>
> > ### {% icon solution %} Solution
> > 1. The overfitting problem will never come when we use the random forest algorithm in any classification problem.
> > 2. The same random forest algorithm can be used for both classification and regression task.
> > 3. The random forest algorithm can be used for feature engineering, Which means identifying the most important features out of the available features from the training dataset.
> {: .solution}
>
{: .question}


After learning on the training dataset, we should evaluate the performance on the test dataset.

> ### {% icon hands_on %} Hands-on: Predict targets using the random forest
>
> 1. **Ensemble methods for classification and regression** {% icon tool %} with the following parameters to predict targets of test dataset using the trained model:
>    - *"Select a Classification Task"*: `Load a model and predict`
>       - {% icon param-file %} *"Models"*: `RandomForestClassifier`
>       - {% icon param-file %} *"Data (tabular)"*: `train_rows_test.csv`
>       - {% icon param-check %} *"Does the dataset contain header"*: `Yes`
>       - {% icon param-select %} *"Select the type of prediction"*: `Predict class labels`
> 2. Rename the generated file to `predicted_class_random_forest`
{: .hands_on}


1. 
    ![roc_scores](images/roc_rf.png "Receiver operator characteristics (ROC) and area under ROC (AUC) for the random forest classifier.")

These plots are important to visualize the quality of regression and the true and predicted targets - how close or far they are from each other. The closer they are, the better is the prediction.



> ### {% icon question %} Question
>
> Inspect the plots. What can you say about the classification?
>
> > ### {% icon solution %} Solution
> >
> > Figures show that we achieved an AUC score of `1.0`  for the test set using random forest. It means the prediction is very good and without any error.
> {: .solution}
{: .question}


# Create data processing pipeline

At the final step, we will create a pipeline with **Pipeline builder** tool but this time, we just specify the classifier. The **Pipeline builder** tool will wrap this classifier and return a zipped file. We will use this zipped file with **Estimator attributes** tool set the search space of hyperparameters.

> ### {% icon hands_on %} Hands-on: Create pipeline
>
> 1. **Pipeline builder** {% icon tool %} with the following parameters:
>    - In *"Final Estimator"*:
>        - *"Choose the module that contains target estimator"*: `sklearn.ensemble`
>            - *"Choose estimator class"*: `BaggingClassifier`
>    - In *"Output the final estimator instead?"*: `Final Estimator`
> 
>      We choose `Final Estimator` as we have only the estimator and no preprocessor and need the parameters of only the estimator.
>
{: .hands_on}


## Extract hyperparameters

We use the **Estimator attributes** tool to get a list of different hyperparameters of the estimator. This tool creates a tabular file with a list of all the different hyperparameters of the preprocessors and estimators. This tabular file will be used in the **Hyperparameter search** tool to populate the list of hyperparameters with their respective (default) values.

> ### {% icon hands_on %} Hands-on: Estimator attributes
>
> 1. **Estimator attributes** {% icon tool %} with the following parameters:
>    - {% icon param-files %} *"Choose the dataset containing estimator/pipeline object"*:  `final estimator builder` file (output of **Pipeline builder** {% icon tool %})
>    - *"Select an attribute retrieval type"*: `Estimator - get_params()`
>
{: .hands_on}

## Search for the best values of hyperparameters

After extracting the parameter names from the **Pipeline builder** file, we will use the **Hyperparameter search** tool to find the best values for each hyperparameter. These values will lead us to create the best model based on the search space chosen for each hyperparameter. We use only one parameter `n_estimators` of `BaggingClassifier` for this task. This parameter specifies the number of bagging stages the learning process has to go through. The default value of `n_estimators` for this regressor is `10`. But, we are not sure if this gives the best accuracy. Therefore, it is important to set this parameter to different values to find the optimal one. We choose some values which are less than `10` and a few more than `10`. The hyperparameter search will look for the optimal number of estimators and gives the best-trained model as one of the outputs. This model is used in the next step to classify the test dataset.


> ### {% icon details %} 5-fold cross-validation
>
> It is a model validation technique which estimates the performance of a predictive model on an unseen data. A dataset is divided into `5` folds and these folds are categorised into training and validation sets. The idea of cross-validation is shown in figure [3](#figure-3). The complete dataset is divided into `5` equal parts. 80% of the dataset is used for training and the remaining 20% is used for validating the performance of training. This is done for `5` folds/iterations, each time the validation set (20% of the dataset) is different. In all five folds, the complete dataset is used for training and validation. The final validation performance is averaged over `5` folds.
>
> ![5fold_cv](../../images/age-prediction-with-ml/5fold_cv.png "5-fold cross validation. ")
>The image shows how the 5-fold cross-validation works. The complete dataset is divided into 5 equal parts/folds. 4 parts (80%) of the data (training set shown in yellow) are used for training the model and the remaining one part is used for evaluating (validation set shown in blue) the trained model. This is repeated for 5 times till every part/fold is used as the validation set. The accuracies computed for different validation folds are averaged to give 5-fold cross-validation accuracy.
>
{: .details}

> ### {% icon hands_on %} Hands-on: Hyperparameter search
>
> 1. **Hyperparameter search** {% icon tool %} with the following parameters:
>    - *"Select a model selection search scheme"*: `GridSearchCV - Exhaustive search over specified parameter values for an estimator `
>        - {% icon param-files %} *"Choose the dataset containing pipeline/estimator object"*: `zipped` file (output of **Pipeline builder** {% icon tool %})
>        - {% icon param-files %} *"Is the estimator a deep learning model?"*: `NO` {% icon tool %})
>        - In *"Search parameters Builder"*:
>             - {% icon param-files %} *"Choose the dataset containing parameter names"*: `tabular` file (output of **Estimator attributes** {% icon tool %})
>             - In *"Parameter settings for search"*:
>                 - {% icon param-repeat %} *"1: Parameter settings for search"*
>                    - *"Choose a parameter name (with current value)"*: `n_estimators: 10`
>                    - *"Search list"*: `[5,10,20,50]`
>    - *"Select input type"*: `tabular data`
>        - {% icon param-files %} *"Training samples dataset"*: `train_rows` tabular file
>        - *"Does the dataset contain header"*: `Yes`
>        - *"Choose how to select data by column"*: `All columns BUT by column header name(s)`
>            - *"Type header name(s)"*: `Class`
>        - {% icon param-files %} *"Dataset containing class labels or target values"*: `train_rows` tabular file
>        - *"Does the dataset contain header"*: `Yes`
>        - *"Choose how to select data by column"*: `Select columns by column header name(s)`
>            - *"Type header name(s)"*: `Class`
>    - *"Whether to hold a portion of samples for test exclusively?"*: `Nope`
>    - *"Save best estimator?"*: `Fitted best estimator or Detailed cv_results_from nested CV`
>
{: .hands_on}

> ### {% icon question %} Question
>
> What is the optimal number of estimators for the given dataset?
>
> Hint: Please look at the `mean_test_score` column in the tabular result from the **Hyperparameter search** tool.
>
> > ### {% icon solution %} Solution
> >
> > 20 (Even though the default value of the number of estimators for Bagging Classifier is `10`, `20` gives the best accuracy. That's why it is important to perform hyperparameter search to tune these parameters for any dataset). 
> >
> {: .solution}
>
{: .question}

Using the **Hyperparameter search** tool, we found the best model, based on the training data. Now, we will predict age in the test dataset using this model.

> ### {% icon hands_on %} Hands-on: Predict age
>
> 1. **Ensemble methods for classification and regression** {% icon tool %} with the following parameters:
>    - *"Select a Classification Task"*: `Load a model and predict`
>        - {% icon param-files %} *"Models"*: `zipped` file (output of **Hyperparameter search** {% icon tool %})
>        - {% icon param-files %} *"Data (tabular)"*: `test_rows` tabular file
>        - *"Does the dataset contain header"*: `Yes`
>
{: .hands_on}


Now we will verify the performance by creating and analysing the plots:


1. 
    ![confusion_matrix](images/confusion_matrix_bagging.png "Confusion matrix for the bagging classifier.")

2. 
    ![prf1_scores](images/precision_recall_bagging.png "Precision, recall and F1 score for the bagging classifier.")

3. 
    ![roc_scores](images/roc_bagging.png "Residual plot between residual (predicted - true) and predicted targets. The plot shows a random pattern of points.")


Figure ... shows that we again achieved  AUC `1.00` which shows that our model is highly effective at predicting whether or not a molecule is biodegradable.


# Conclusion
By following these steps, we learned how to build classifiers and visualize the classification results using Galaxy's machine learning and plotting tools. The features of the training dataset are mapped to the classes. This mapping is used to make predictions on an unseen (test) dataset. The quality of classifiers is visualized using a plotting tool. There are multiple other classification algorithms, few are simpler to use (with fewer parameters) and some are powerful, which can be tried out on this dataset and on other datasets as well. Different datasets can also be analyzed using these classifiers. The classifiers have lots of parameters which can be altered while performing the analyzes to see if they affect the classification accuracy. It may be beneficial to perform hyperparameter search to tune these parameters of classifiers for different datasets. In addition, we learned the relevance of machine algorithms for QSAR analyses and constructed a model which successfully predicted an important chemical property - the biodegradability of a substance.
