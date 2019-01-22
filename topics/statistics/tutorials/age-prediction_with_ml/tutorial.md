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
key_points:
    - Various machine learning algorithms should be used to find the best ones
    - For each machine learning algorithm, it hyperparameters should be optimised based on the dataset
    - Feature selection should be done for high-dimensional datasets
time_estimation: 2H
contributors:
- polkhe
- anuprulez

---


# Introduction
{:.no_toc}

[Machine Learning](https://en.wikipedia.org/wiki/Machine_learning) is used to create predictive models by learning features from datasets. In the studies [Jason G. Fleischer et al. 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1599-6#Sec9) and [Jana Naue et al. 2017](https://www.sciencedirect.com/science/article/pii/S1872497317301643?via%3Dihub), biomarkers were examined to predict choronological age of humans by analysing the RNA-seq gene expression levels and DNA methylation pattern respectively. Different machine learning algorithms were used in these studies to select specific biomarkers to make age prediction. The RNA-seq gene expression ([FPKM](https://www.ebi.ac.uk/training/online/glossary/fpkm)) dataset was generated using fibroblast cell lines from 133 healthy humans and their ages range from 1 to 94 years. Within each individual, [DNA methylation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3174260/) changes with age. This knowledge was used to select useful biomarkers from DNA methylation dataset. The [CpGs sites](https://en.wikipedia.org/wiki/CpG_site) with the highest correlation to age were selected as biomarkers/features. In both these studies, specific biomarkers were analysed by machine learning algorithms to create an age prediction model. [Regression](https://en.wikipedia.org/wiki/Regression_analysis) is used as the machine learning task in both the studies because the targets/labels/ages are real numbers.

Using machine learning tools in Galaxy, we can achieve comparable prediction scores as achieved by these analyses. In this tutorial, we will apply a couple of ([scikit-learn](https://scikit-learn.org/stable/)) machine learning tools to [RNA-seq](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1599-6#Sec9) and [DNA methylation](https://www.sciencedirect.com/science/article/pii/S1872497317301643?via%3Dihub) datasets to predict the chronological age of humans. These datasets have been taken directly from these studies. The training is done over 5-fold cross-validation and [r2](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.r2_score.html) (coefficient of determination) scoring metric is used to evaluate the performance of the trained model. The closer it is to 1.0, the better it is. If it is negative, then the trained model is not good. To infer how its values exhibit model performance, we can compare the figures 1 and 2. In both the plots, the true and predicted targets are plotted in a scatter plot. For a good model, most of the points should lie along the `x = y` line as the true and predicted targets are close to each other. In figure 1, we can see that the points are scattered and do not show any pattern. Therefore, the r2 score is `-0.06`. But, figure 2 shows a better pattern as most of the points lie along the line and the r2 score is almost `1.0`.

>    ![model_bad](images/model_bad.png "This shows an example of a bad model as most of the points are scattered.")

>    ![model_good](images/model_good.png "This shows an example of a good model as most of the points lie along the x = y line.")

This tutorial is divided into two parts - one with RNA-seq dataset and another with DNA methylation dataset. We will use the r2 metric to evaluate the performance of the trained model. For RNA-seq dataset, we will compute cross-validated r2 score using the training set and for DNA methylation dataset, we will compute the r2 score for the test set.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Analyze RNA-seq dataset 

The RNA-seq dataset is collected from fibroblast cell lines belonging to 133 healthy patients of age from 1 to 94 years. Using this dataset, we perform an [exhaustive search](https://en.wikipedia.org/wiki/Hyperparameter_optimization) (also known as grid search) for finding the best features/genes and then apply [ElasticNet regressor](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.ElasticNet.html#sklearn.linear_model.ElasticNet) with 5-fold cross-validation. The r2 score achieved is comparable to the predictions found in the study [Jason G. Fleischer et al. 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1599-6#Sec9).

> ### {% icon details %} Cross-validation
>
> It is a model validation technique which estimates the performance of a predictive model on an unseen data. A dataset is divided into `k` folds (parts) and these folds are categorised into training and validation sets. For example in 5-fold cross-validation, 4 parts of the dataset are used for training and 1 part is used for validating the performance. The final performance is averaged over `k` folds.
>
{: .details}

## Get the raw data

We proceed to the analysis by uploading the dataset. The dataset has `133` rows corresponding to humans and over `27,000` columns specifying genes. The last column `age` contains the chronological age of humans. For the validation set, this `age` column is predicted and r2 metric is computed.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo](https://zenodo.org/record/2545213#.XEWTJ9-YVa0)
>
>    ```
>    https://zenodo.org/record/2545213/files/training_data_normal.tsv
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

## Create a data processing pipeline

We can see that this RNA-seq dataset is high-dimensional. There are over `27,000` columns/features. Generally, not all the features in the dataset are useful for prediction. We need only those features which increase the predictive ability of the model. To filter these features, we perform feature selection and retain only those which are useful. To do that, we use `SelectKBest` module in the suite of data preprocessors. Again, we are not sure of how many of these features we will need. To find the right number of features, we do a hyperparameter search by setting a different number of features and find out the best number. To create this preprocessor, we will use **pipeline builder** tool. This tool defines sequential processing of datasets. After the preprocessing step, we should add a regressor algorithm (ElasticNet) which analyzes the preprocessed dataset. This tool gives a `zip` file as an output containing the specifications of different steps of the entire analysis.

> ### {% icon hands_on %} Hands-on: Create pipeline
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

In any machine learning algorithm, there are many parameters (hyperparameters). We are not sure which values of these parameters will give an optimal prediction. The default values given for these parameters may not be optimal for different datasets. To find the best combination of the values of different parameters for a dataset, [hyperparameter optimisation](https://en.wikipedia.org/wiki/Hyperparameter_optimization) is performed. There are different techniques to optimise the hyperparameters of any algorithm given a dataset:
- [Grid search](https://scikit-learn.org/0.16/modules/generated/sklearn.grid_search.GridSearchCV.html)
- [Random search](https://scikit-learn.org/0.16/modules/generated/sklearn.grid_search.RandomizedSearchCV.html#sklearn.grid_search.RandomizedSearchCV)

For our analyses, we will use the grid search approach. It is an exhaustive search which tries out all the combinations of different hyperparameters and ranks these combinations based on a scoring metric. In the random search, the values of a parameter are selected randomly from a given range and the best one is found. Grid search works well for parameters taking categorical as well as numerical values but for the random search, it becomes difficult for parameters which take categorical values.

In the pipeline builder, we added two steps - preprocessing (feature selection) and an estimator (regressor). There are different hyperparameters for these two steps and their best combination should be found out. We will perform grid search to estimate the best values for these parameters: **k** (number of features), **normalize** (subtract the mean and divide by the l2-norm of the dataset) and **alpha** (a constant which is multiplied to the regularisation term). For each parameter, we need to specify a set of values to choose from:
- **k**: [5880, 5890, 5895, 5900]
- **normalize**: [True, False]
- **alpha**: [0.00001, 0.0001, 0.001]

For these three parameters, we have 24 different combinations (4 x 2 x 3) of values and we will verify the performance of each combination. You might have noticed that the parameter **k** is used for feature selection and parameters **normalize** and **alpha** are used for regressor.

> ### {% icon hands_on %} Hands-on: Hyperparameter search
>
> 1. **Hyperparameter Search** {% icon tool %} with the following parameters:
>    - *"Select a model selection search scheme:"*: `GridSearchCV - Exhaustive search over specified parameter values for an estimator ` 
>        - {% icon param-files %} *"Choose the dataset containing pipeline object"*: `Zipped` file (output of **Pipeline Builder** {% icon tool %})
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
>            - *"Raise fit error:"*: `No`
>    - *"Select input type:"*: `tabular data`
>        - {% icon param-files %} *"Training samples dataset:"*: `training_data_normal` tabular file
>        - *"Does the dataset contain header:"*: `Yes`
>        - *"Choose how to select data by column:"*: `All columns BUT by column header name(s)`
>            - *"Type header name(s):"*: `age`
>        - {% icon param-files %} *"Dataset containing class labels or target values"*: `training_data_normal` tabular file
>        - *"Does the dataset contain header:"*: `Yes`
>        - *"Choose how to select data by column:"*: `Select columns by column header name(s)`
>            - *"Type header name(s):"*: `age`
>
>    > ### {% icon comment %} Comment
>    >
>    > The tool returns two outputs, one of which is a table with numerical results. Please inspect it carefully: the `rank_test_score` column shows the ranking of different combinations based on the values in `mean_test_score` column.
>    {: .comment}
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. What is the best `mean_test_score` value estimated by the hyperparameter search tool?
> 2. Which combination of parameters gives the best result?
> 3. How many possible combinations of parameters the hyperparameter search tool estimated?
>
> > ### {% icon solution %} Solution
> >
> > 1. 0.73
> > 2. alpha: 0.001, normalize: True, k: 5880
> > 3. 24
> >
> {: .solution}
>
{: .question}

## Parallel coordinates plot

We will visualize the tabular output of hyperparameter search tool from the previous step using **parallel coordinates plot of tabular data**. There are multiple columns in the tabular output, but we will focus on only a few of them.

> ### {% icon hands_on %} Hands-on: Parallel coordinates plot
>
> 1. **Parallel Coordinates Plot** {% icon tool %} with the following parameters:
>    - {% icon param-files %} *"Select data file:"*: 'Tabular' file (output of **Hyperparameter Search** {% icon tool %})
>    - *"Select the columns for dimensions:"*: `c5, c6, c7, c14`
>    - *"Select a column containing the values for coloring:"*: `c14`
>
>    > ### {% icon comment %} Comment
>    >
>    > The output plot has the following legend: the colour-coding is based on the `mean_test_score` column. You can follow the line leading
>    > to the score along every column with parameters' settings.
>    {: .comment}
>    ![data](images/plotting_output.png "The visualization of the hyperparameter optimisation tool output.")
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
In the plot shown above (figure 3), we achieved an r2 score of `0.73` (last column) with 5-fold cross-validation on the training set. In the study [Jason G. Fleischer et al. 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1599-6#Sec9) as well, a similar r2 score is mentioned for linear regressors. Moreover, the study also included a customised ensemble regressor which achieved better performance (`r2 = 0.81`). However, our analysis showcases the use of machine learning tools in Galaxy to reproduce the results published in the paper.

# Analyze DNA methylation dataset
In the second part of the analysis, we will use DNA methylation dataset to predict chronological age. The dataset is divided into two parts - training and test sets. The training set is used to train a regressor and the test set is used to evaluate the performance of the trained model using the r2 scoring metric. 5-fold cross-validation is used for training.

## Get the train and test datasets

We proceed with the analysis by uploading new datasets. You might want to create a new history first. The training set contains `208` rows corresponding to humans and `14` features. The last column is `age`. The test set contains `104` rows and the same number of features as the training set. The `age` column in the test set is predicted after training on the training set. Another dataset `test_rows_labels` contains the true age values of the test set which is used to compute r2 scores between true and predicted age.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo](https://zenodo.org/record/2545213#.XEWTJ9-YVa0)
>
>    ```
>    https://zenodo.org/record/2545213/files/train_rows.csv
>    https://zenodo.org/record/2545213/files/test_rows_labels.csv
>    https://zenodo.org/record/2545213/files/test_rows.csv
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
>    > ### {% icon comment %} Comment
>    > The `train_rows` contains a column `Age` which is the label or target. We will evaluate our model on
>    > `test_rows` and compare the predicted age with the true age in `test_rows_labels`
>    {: .comment}
{: .hands_on}

## Create a data processing pipeline

We will create a pipeline with **pipeline builder** tool but this time, we just specify the regressor, [gradient boosting regressor](https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.GradientBoostingRegressor.html). This is an ensemble-based regressor which uses multiple tree-based regressors internally and predicts by taking an ensemble of the predictions. The tool will output a zipped file.

> ### {% icon hands_on %} Hands-on: Create pipeline
>
> 1. **Pipeline Builder** {% icon tool %} with the following parameters:
>    - In *"Final Estimator:"*:
>        - *"Choose the module that contains target estimator:"*: `sklearn.ensemble`
>            - *"Choose estimator class:"*: `GradientBoostingRegressor`
>
>    > ### {% icon comment %} Comment
>    >
>    > [*Ensemble*](https://en.wikipedia.org/wiki/Ensemble_learning) method uses multiple learning models internally for better predictions. 
>    {: .comment}
>
{: .hands_on}

## Hyperparameter optimisation

For this analysis as well, we will use **hyperparameter search** tool to estimate the best values of parameters for the given dataset. We use only one parameter `n_estimators` of `Gradient Boosting` regressor for this task. This parameter specifies the number of boosting stages the learning process has to go through. The hyperparameter search will look for the optimal number of estimators and gives the best-trained model as one of the outputs. This model is used in the next step to predict age in the test dataset.

> ### {% icon hands_on %} Hands-on: Hyperparameter optimisation
>
> 1. **Hyperparameter Search** {% icon tool %} with the following parameters:
>    - *"Select a model selection search scheme:"*: `GridSearchCV - Exhaustive search over specified parameter values for an estimator `
>        - {% icon param-files %} *"Choose the dataset containing pipeline object"*: `Zipped` file (output of **Pipeline Builder** {% icon tool %})
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
>            - *"Raise fit error:"*: `No`
>    - *"Select input type:"*: `tabular data`
>        - {% icon param-files %} *"Training samples dataset:"*: `train_rows` tabular file
>        - *"Does the dataset contain header:"*: `Yes`
>        - *"Choose how to select data by column:"*: `All columns BUT by column header name(s)`
>            - *"Type header name(s):"*: `Age`
>        - {% icon param-files %} *"Dataset containing class labels or target values"*: `train_rows` tabular file
>        - *"Does the dataset contain header:"*: `Yes`
>        - *"Choose how to select data by column:"*: `Select columns by column header name(s)`
>            - *"Type header name(s):"*: `Age`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. What is the optimal number of estimators for the given dataset?
>
> > ### {% icon solution %} Solution
> >
> > 1. 75
> >
> {: .solution}

## Predict age

Using the hyperparameter search, we found the best model based on the training data. Now, we will predict age in the test dataset using this model.

> ### {% icon hands_on %} Hands-on: Prediction
>
> 1. **Ensemble methods for classification and regression** {% icon tool %} with the following parameters:
>    - *"Select a Classification Task:"*: `Load a model and predict`
>        - {% icon param-files %} *"Models"*: `Zipped` file (output of **Hyperparameter Search** {% icon tool %})
>        - {% icon param-files %} *"Data (tabular)"*: `test_rows` tabular file
>        - *"Does the dataset contain header:"*: `Yes`
>        - *"Select the type of prediction:"*: `Predict class labels`
>
{: .hands_on}

We will plot the true and predicted age and verify the performance.

> ### {% icon hands_on %} Hands-on: Plots
>
> 1. **Plot actual vs predicted curves and residual plots of tabular data** {% icon tool %} with the following parameters:
>    - {% icon param-files %} *"Select input data file :"*: `test_rows_labels` tabular file
>    - {% icon param-files %} *"Select predicted data file :"*: `Tabular` file (output of **Ensemble methods for classification and regression** {% icon tool %})
>
>    > ### {% icon comment %} Comment
>    >
>    > The tool outputs three html files with the interactive plots.
>    {: .comment}
>    > ![Scatter plot](images/scatter_plot.png "Scatter plot for true and predicted age")
>    > We can see in the scatter plot (figure 4) that most of the points lie along the x=y curve. It means that the true and predicted ages are 
>    > close to each other. The root mean square error (`RMSE`) is `3.76` and r2 score is `0.94`.
>    > ![Residuals](images/residual_plot.png "Residuals")
>    > The residual plot (figure 5) is generated to see if there is any visible pattern between residual (predicted age - true age)
>    > and predicted  age. For a good model, there should not be any visible pattern with the plotted points.
>    >
>    > ![True vs predicted age](images/true_vs_predicted_plot.png "True vs predicted age")
>    > The plot in figure 6 shows the true and predicted ages. It can be seen that the points are close.
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. Inspect the plots. What can you say about the predictions?
>
> > ### {% icon solution %} Solution
> >
> > 1. Figures 4 and 6 show that the prediction is good because the predicted age lie close to the true age.
> >
> {: .solution}
>
{: .question}

## Conclusion
We can see in figure 4 that we have achieved an r2 score of `0.94` and root mean square score of `3.76` for the test set using gradient boosting regressor. In the study [Jana Naue et al. 2017](https://www.sciencedirect.com/science/article/pii/S1872497317301643?via%3Dihub) as well, a similar root mean square score (`3.93`) is mentioned using random forest regressor. The root mean square score shows the difference in the true and predicted age of humans. The r2 score (`0.94`) is close to the best achievable score of `1.0` which shows that the trained model is good. Overall, the second part of the analysis also shows that using the machine learning tools in Galaxy, we can achieve state-of-the-art predictions mentioned in the recent scientific studies.
> 
{:.no_toc}

