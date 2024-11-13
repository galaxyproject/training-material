---
layout: tutorial_hands_on
level: Intermediate
title: Building the LORIS LLR6 Model Using Galaxy-PyCaret Tool
zenodo_link: https://zenodo.org/records/13885908
questions:
- How can I reproduce the LORIS LLR6 model published by Chang et al., 2024?
- Which tools can I use in Galaxy to obtain a logistic regression model?
- How can I evaluate the model to confirm its performance?
objectives:
- Use a large dataset of immune checkpoint blockade (ICB)-treated and non-ICB-treated patients across 18 solid tumor types, encompassing a wide range of clinical, pathologic and genomic features to build a Machine Learning Model.
- Build a Machine Learing model using PyCaret tool available in Galaxy.
- Evaluate the models for reproducibility and robustness by comparing them with the original LORIS LLR6 model published by Chang et al., 2024.

time_estimation: 1H
key_points:
- Use Galaxy tools to build similar LORIS LLR6 model published by Chang et al., 2024.
- Use the PyCaret tool to train a new model and confirm the robustness of the study published by Chang et al., 2024.
contributors:
- paulocilasjr 
- qchiujunhao 
- jgoecks
tags:
- LORIS Score Model
- Machine Learning
- PyCaret
- Galaxy-ML tools
---

> <comment-title>PyCaret Model Comparison Tool</comment-title>
>
> The PyCaret Model Comparison tool described in this tutorial is only available (for now) at: 
> [Cancer-Galaxy](https://cancer.usegalaxy.org)
>
> Galaxy-ML tools > PyCaret Model Comparison 
>
> As soon as it gets incorporated to the main galaxy project, this tutorial will be updated.
{:  .comment}

Using a comprehensive dataset of patients treated with immune checkpoint blockade (ICB) and non-ICB-treated patients across 18 solid tumor types, we will develop LORIS (logistic regression-based immunotherapy-response score). The goal is to accurately predict patient responses to the treatment.

To achieve this, we will follow three essential steps: (i) upload patient data training file to Galaxy, (ii) set upt and run the PyCaret Model Comparison Tool to training the best model (iii) evaluate the predictive performance of the model comparing with LORIS model published ({% cite Chang2024 %}).

![schema of the whole process of training model and test.](../../images/loris_tutorial/tutorial_schema.png "Overview of the process steps to obtain the model from the LORIS dataset.")

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

> <comment-title>Background</comment-title>
>
> [LORIS dataset](https://github.com/rootchang/LORIS/blob/main/02.Input/AllData.xlsx) comprises clinical, pathologic, and genomic features
> from a diverse cohort of patients, enabling robust analysis and model development.
> There are 10 different cohort dataset available in the raw data (xlsx): 
> 1)Chowell_train, 2) Chowell_test, 3) MSK1, 4) MSK2, 5) Kato_panCancer, 6) Shim_NSCLC, 7) Vanguri_NSCLC, 8) Ravi_NSCLC, 9) Pradat_panCancer, 10) MSK_nonICB
>
> For the proporse of being conscise in this tutorial, we are going to focus only on 1 cohort: Chowell_train.
>
{:  .comment}

# Dataset composition to train the model 

Before we begin the hands-on session, here’s a brief explanation of the features we’ll be using. These features were selected based on the findings of Chang et al. (2024), identifying them as the most important for training the model.

## TMB (Tumor mutation burden)

Tumor mutational burden (TMB) is defined as the total number of somatic mutations within a specified region of the tumor genome. It has been recognized as a biomarker for predicting the efficacy of immune checkpoint blockade (ICB) in solid tumors. The FDA has approved a threshold of 10 mutations per megabase (Mb) as a biomarker for response to ICB treatment.

In this dataset, TMB values range from 0 to over 368 mutations per megabase (mut/Mb), with several extreme values like 368.6 and 93.5. To mitigate the influence of these outliers, TMB values will be truncated at 50 mut/Mb, meaning any value exceeding 50 will be capped at 50. This is crucial because extreme TMB values can disproportionately skew the model's learning process, leading to unreliable predictions.

## Systemic Therapy History

This feature is a binary variable that indicates whether a patient received chemotherapy or targeted therapy prior to immunotherapy. It is coded as 1 if the patient had undergone such treatments before starting immunotherapy and 0 if they had not. The Systemic_therapy_history feature is a binary variable that indicates whether a patient received chemotherapy or targeted therapy prior to immunotherapy. It is coded as 1 if the patient had undergone such treatments before starting immunotherapy and 0 if they had not. 

## Albumin

The albumin feature represents the albumin levels measured in patients, which is an important biomarker often associated with nutritional status and liver function. The values are measured in grams per deciliter (g/dL) and typically range between 2.1 and 4.9 g/dL in the dataset. Higher levels of albumin are generally associated with better overall health and can serve as an indicator of a patient's ability to recover or respond to treatments like immunotherapy.

## Cancer Type

The CancerType feature represents the type of cancer diagnosed in each patient, which can vary significantly across the dataset. Common types include Non-Small Cell Lung Cancer (NSCLC), Small Cell Lung Cancer (SCLC), Melanoma, Endometrial cancer, and other cancer types such as Gastric, Colorectal, Renal, and Breast cancer. This feature is critical for understanding the heterogeneity of the patient cohort and may influence treatment decisions, response rates, and outcomes. 

Incorporating this feature into a machine learning model requires translating the categorical CancerType into one-hot encoded variables. Each cancer type will be represented as a binary feature (0 or 1), with each type becoming a separate column in the dataset. This enables the model to interpret the presence or absence of a specific cancer type for each patient.

## NLR (blood neutrophil–lymphocyte ratio)

The neutrophil–lymphocyte ratio (NLR), a biomarker derived from the ratio of neutrophils to lymphocytes in the blood, is increasingly used in cancer research due to its association with inflammation and immune response. It can serve as a prognostic factor in various cancer types. Higher NLR values often indicate a poorer prognosis, potentially reflecting a more aggressive disease or impaired immune response.

In this dataset, NLR values range, for example from 0.8 to 88, with several extreme outliers. To address this, NLR values will be truncated at 25, meaning any value above 25 will be set to 25. This truncation is important for preventing extreme outliers from disproportionately influencing the machine learning model.

## Age

In predictive models for patient outcomes, age is a crucial feature because it is often correlated with various health factors and disease risks. As people age, their immune systems, metabolism, and ability to recover from illnesses may change, influencing how they respond to treatments, medications, or disease progression. Including age as a feature helps models account for the biological changes that occur over time and can improve the accuracy of predictions across different age groups.

However, there are limits to how predictive age might be, particularly for extreme values. For example, patients over a certain age may share similar health characteristics, and further increases in age may not significantly add predictive value. Truncating age to a maximum value (like 85) helps to avoid overemphasizing small differences between very old patients, where the added predictive power might be negligible.

## Response

The feature Response is a categorical target variable indicating whether patients benefited from immune checkpoint blockade (ICB) therapy, classified as 0 (no benefit) or 1 (benefit). The model is trained using all previous features to predict patient outcomes when ICB is chosen as the treatment.

# Prepare environment and get the data 

> <comment-title>Preprocessing the raw data</comment-title>
>
> The raw data published by ({% cite Chang2024 %}) can be found here: 
> [LORIS raw dataset](https://github.com/rootchang/LORIS/blob/main/02.Input/AllData.xlsx)
> 
> We preprocessed the raw data using a Python script to:
> 1) Extract the `Chowell_Train` tab from the excel file.
> 2) Select the 7 features (`TMB`, `Systemic Therapy History`, `Albumin`, `Cancer Type`, `NLR`, `Age`, and `Response`) important for building the model.
> 3) Truncate the values for `Age`, `NLR`, and `TMB`.
> 4) Encode `Cancer Type` using one-hot encoding. 
> 5) Save the dataset as a .tsv file.
>
> A jupyter notebook can be found at Dockstore: [LORIS_preprocessing](https://dockstore.org/notebooks/github.com/paulocilasjr/pycaret-use-case/preprocessing:main?tab=info)
>
{:  .comment}

> <hands-on-title> Environment and Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial. If you are not inspired, you can name it *LORIS model classifier*.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the files from Zenodo or from the shared data library
>
>    ```
>    https://zenodo.org/records/13885908/files/Chowell_train_Response.tsv
>    https://zenodo.org/records/13885908/files/Chowell_test_Response.tsv
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Check that the data format assigned for the file is **tsv**.
>    If it is not, follow the Changing the datatype tip.
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> 5. Add a tag (`LORIS model dataset`) to the dataset corresponding to `Chowell_train_Response.tsv` and `Chowell_test_Response.tsv`
>    This is important to trace back on what dataset the model was built on.
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

# Using PyCaret Model Comparison Tool

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [PyCaret Model Comparison](https://toolshed.g2.bx.psu.edu/repos/paulo_lyra_jr/pycaret_model_comparison/PyCaret_Model_Comparison/2024.3.3.2+0) %} with the following parameters:
>    - {% icon param-file %} *"Input Dataset (CSV or TSV)"*: `Chowell_train_Response.tsv`
>    - {% icon param-file %} *"Test Dataset (CSV or TSV)"*: `Chowell_test_Response.tsv`
>    - {% icon param-file %} *"Select the target column"*: `C22: Response`
>    - {% icon param-file %} *"Task"*: `Classification`
> Run the tool 
{: .hands_on}

# Tool output files

After your model is trained and tested, you should see two new files in your history list:

- PyCaret Model Comparasion Best Model: The PyCaret model pickle file. This file allows the model to be reused without requiring retraining, ensuring consistent predictions.

- PyCaret Model Report: The file containing all the plots for the models trained and the best model selected.

As the proporse of this tutorial, we will focus on the PyCaret Model Report.

# PyCaret Model Report 

The PyCaret HTML report output provides a comprehensive and interactive overview of the trained model’s performance in an accessible, browser-ready format. This HTML report documents key aspects of the model’s training and evaluation process, supporting deeper insight into how well the model performed on both training and test datasets. There are four tabs: Setup & Best Model; Best Model Plots; Feature Importance; Explainer. 

We are going to provide a brief explanation about the content of each tab present in the report.

> <tip-title>Setup & Best Model Tab</tip-title>
>- Setup Parameters: Documents the initial configurations used in the PyCaret Model Comparison Tool.
>- Best Model Class and Hyperparameters: Specifies the model selected as the best performer from the model comparison, listing the model’s hyperparameters as chosen through tuning.
>- Performance Metrics: Summarizes key evaluation metrics, including Accuracy, ROC-AUC, Recall, Precision, F1-Score, Cohen’s Kappa, Matthews Correlation Coefficient (MCC), and Training Time (TT in seconds).
>![alt](... "label")
>
{: .tip}

> <tip-title>Best Model Plots Tab</tip-title>
>- Training and Test Evaluation Plots: Displays visualizations of the model’s performance, including an ROC-AUC curve for binary classification, Confusion Matrix, Precision-Recall (PR) curve, Error Analysis, and a Classification Report for detailed class-level performance.
>- Additional Model Insights: Includes diagnostic plots like the Learning Curve, Calibration Plot, Validation Curve (VC), and various dimensionality reduction plots (Manifold, RFE).
>- Feature Importance: Shows the contribution of each feature to the model, both individually and collectively (All Features), providing insight into the factors influencing model decisions.
>![alt](... "label")
>
{: .tip}

> <tip-title>Feature Importance Tab</tip-title>
>- Feature importance analysis from atrained Random Forest: While these metrics are not directly related to LASSO logistic regression, they can offer additional insights into how features behave in alternative models, potentially validating or constrasting with the model's results.
>- SHAP summary from a trained lightGBM model: SHAP (SHapley Additive exPlanations) is particularly useful for understanding feature impact. SHAP values provide a unified measure of feature contribution of each feature's influence on predictions.
>![alt](... "label")
>
{: .tip}

> <tip-title>Explainer Tab</tip-title>
>- Mean Absolute SHAP Value (Average Impact on Predicted Response): This shows the overall contribution of each feature to predictions. Higher SHAP values indicate features with a larger influence on the model, helping identify which inputs are most impactful.
>- Permutation Importance (Decrease in ROC AUC with Randomized Feature): By randomly shuffling each feature and observing the drop in ROC AUC, you can see how crucial each feature is for model performance. A large drop indicates high importance.
>- Partial Dependence Plot (PDP): PDPs show how a feature affects predictions on average. For each feature, observe the plot's curve to understand if it has a linear, nonlinear, or interaction effect on the prediction.
>- Percentage 1 vs. Predicted Probability: This plot compares the true proportion of positive cases (label=1) to predicted probabilities, helping assess the model’s calibration—whether predicted probabilities match observed outcomes.
>- Cumulative Percentage per Category (Top X%): This measures the cumulative impact of each category when sampling the top percentage. It's helpful for understanding feature value distribution or concentration.
>- Percentage Above and Below Cutoff: This evaluates the model's performance above and below a chosen threshold, offering insight into sensitivity and specificity at different probability cutoffs.
>- Confusion Matrix: A matrix of True Positives, False Positives, True Negatives, and False Negatives, allowing you to assess the model's accuracy, error types, and balance of predictions.
>- Lift Curve: The lift curve shows the improvement over random selection for each decile of predictions, with the lift value indicating the effectiveness of the model in capturing positive instances.
>- ROC AUC Curve: This curve plots True Positive Rate vs. False Positive Rate across thresholds, with the AUC score representing the model’s discriminative ability. Higher AUC means better distinction between classes.
>- PR AUC Curve: Plots Precision vs. Recall, focusing on the model’s performance on the positive class. Useful for imbalanced data, with a higher AUC indicating better handling of positive instances without sacrificing precision.
>![alt](... "label")
>
{: .tip}

# LORIS PanCancer LLR6 Model Robustness:

It's crucial to understand the concept we want to drive into this analyse. Since we are looking to generate a model similar to what was published by {% cite Chang2024 %}, we are going to use the model's metrics from the paper to compare the model we generate through Galaxy-PyCaret.

> <comment-title>Robustness definition </comment-title>
>
> Some evidence is robust across reasonable variations in analysis, and some evidence is fragile, meaning that support for the finding is contingent on specific decisions such as which observations are excluded and which covariates are included.
> 
> Thus, Robustness refers to testing the reliability of a prior finding using the same data and a different analysis strategy.
>
>![alt](... "label")
>
{:  .comment}

## Classification Algorithms

One of PyCaret's key features is the ability to train and compare multiple models with minimal code. By default, PyCaret evaluates a broad range of algorithms—including linear models, tree-based models, and ensemble methods—and then ranks these models based on their performance. It uses the Area Under the Curve (AUC) metric as the primary criterion to determine which model is the best.

![alt](... "label")

In this case, the best algorithm to perform in the datase matches what was reported in the article: Logistic Regression model.

## Hyperparameters

{% cite Chang2024 %}'s model for Pan-Cancer LLR6 model has the following hyperparameters set:
C = 0.1, Class Weight = Balanced, l1 ratio = 1, max iter = 100, Penalty = Elasticnet, Solver = Saga.

> <comment-title>Hyperparameters meaning</comment-title>
>
> *Penalty* - Defines the type of regularization applied to the model to prevent overffiting. Options for linear and logistic regression are: `L1`, `L2`, and `Elasticnet`. 
> Briefily, L1 (Lasso Regression) removes not important features, overcomming overffiting as well as in dimension reduction, but when most of the features (variables) in the model are useful, L2 (Ridge Regression) is used. Elasticnet regularization combines both L1 and L2, addressing multicolinearity while also enabling feature selection. When Elasticnet is selected, opens the L1 Ratio parameter.
>
> *L1-Ratio* - Controls the balance between L1 and L2 penalties. A value of `1` uses purely L1 regularization, which encourages sparsity. A value of `0` uses purely L2.
>
> *C* - is the inverse of the regularization strength. A smaller C (like `0.1`) implies stronger regularization, prevents overfitting by penalizing large coefficients, also makes the model simple (smaller coefficients). The contrary, higher values make the model more complex.
>
> *Solver* - Specifies the optmization algorithm used for fitting the model. `SAGA` is well-suited for large datasets and supports the elastinet penalty. Also effective for sparse data and is fast for L1. `LBFGS` is a quasi-newton optimization algorithm, efficient for smaller dataset and suppport L2 regularization but does not support elasticnet.
>
> *Class Weight* - is used to handle imbalanced classes by adjusting the weigth associated with each class. `Balanced` adjusts these weights inversely proportional to class frequencies in the data, more weight to the minority class.
>
> *Max Iter* - Specifies the maximum `number of iterations` the solver will run before stopping. If converge is achieved earlier, it will stop, if not, might need to increase the value.
>
{:  .comment}

Our PyCaret Best Model has the hyperparameters 

# Conclusion

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
