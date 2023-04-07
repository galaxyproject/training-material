---
layout: tutorial_hands_on

title: "Supervised Learning with Hyperdimensional Computing"
zenodo_link: "https://doi.org/10.5281/zenodo.6467875"
questions:
  - "How to encode data into vectors in a high-dimensional space?"
  - "What kind of operations can be performed on these vectors?"
  - "What is a vector-symbolic architecture?"
  - "How to build a classification model out of this architecture?"
objectives:
  - "Learn how to encode data into high-dimensional vectors"
  - "Build a vector-symbolic architecture"
  - "Use the architecture to build a classification model"
requirements:
time_estimation: "30m"
level: Intermediate
key_points:
  - "Every kind of data can be encoded into high-dimensional vectors"
  - "A vector-symbolic architecture is composed of vectors and a limited set of arithmentic operators"
  - "Classification models build according to the hyperdimensional computing paradigm can scale on datasets with massive amounts of features and data points"
contributors:
  - cumbof

---

# Introduction

`chopin2` ({% cite Cumbo2020 %}) implements a domain-agnostic supervised classification method based on the hyperdimensional (HD) computing paradigm. It is ano open-source tool and its code is available on GitHub at [https://github.com/cumbof/chopin2](https://github.com/cumbof/chopin2).

In this tutorial, we are going to work on a dataset with microbial relative abundances (RA) and presence/absence information (BIN) computed on metagenomic samples collected from individual affected by the colorectal cancer (CRC) in a case/control scenario.
The main goal is to build a supervised classification model over the RA profiles with `chopin2` in order to discriminate case and control samples with a high accuracy.

The dataset with RA profiles has been produced with `MetaPhlAn3` ({% cite beghini2021integrating %}) and it is available throught the `curatedMetagenomicData` package for R ({% cite pasolli2017accessible %}). The studies describing the analysis of the metagenomic samples are available in Nature Medicine ({% cite thomas2019metagenomic %} {% cite wirbel2019meta %}).
In particular, the dataset contains the profiles of 241 microbial species detected in 101 stool samples of patients affected by the CRC, and 92 samples collected from the stool of healthy individuals.

The data we use in this tutorial is also available on [Zenodo](https://doi.org/10.5281/zenodo.6467875).

Please note that both the RA and BIN datasets in Zenodo have also been stratified according to the age category and sex of both case and control individuals. However, in this tutorial we are going to analyze the unstratified datasets, called `RA__ThomasAM__species.csv` and `BIN__ThomasAM__species.csv`.

> <agenda-title></agenda-title>
>
> `chopin2` is composed of three steps:
>
> 1. [Get the data](#get-the-data)
> 2. [Build a classification model](#build-a-classification-model)
> 3. [Feature selection](#feature-selection)
{: .agenda}

# Get the data

The first step consists in importing the RA and BIN datasets into the Galaxy history. As previously mentioned, RA refers to the dataset with microbial relative abundance profiles computed with `MetaPhlAn3`, while BIN refers to the dataset with presence/absence information about microbes in samples.

> <hands-on-title>Get the data</hands-on-title>
> 1. Create a new history for this tutorial
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import datasets from [Zenodo](https://zenodo.org/record/7806264):
>    - RA dataset (`RA__ThomasAM__species.csv`)
>    - BIN dataset (`BIN__ThomasAM__species.csv`)
>
>    ```
>    https://zenodo.org/record/7806264/files/RA__ThomasAM__species.csv
>    https://zenodo.org/record/7806264/files/BIN__ThomasAM__species.csv
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    Both the datasets contain samples on the columns and microbial species on the rows. The RA dataset contains relative abundance profiles, thus values ranging from 0.0 to 100.0, while the BIN dataset contains presence/absence (i.e., 0 or 1) information computed on the relative abundance. It means that the value of the cell `ij` in the BIN matrix is 1 if the relative abundance of the species `i` in the sample `j` is greater than 0.0 in the RA matrix. Otherwise, it is 0.
>
> 3. Rename the files {% icon galaxy-pencil %}
>    - The names of files imported in the previous step are the full URL to those files;
>    - Let's change their name to something a little clearer like `RA__ThomasAM__species.csv` and `BIN__ThomasAM__species.csv`.
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
{: .hands_on}

# Build a classification model

Hyperdimensional computing is an emerging brain-inspired computing paradigm that deals with vectors in a high-dimensional space. Every kind of information is thus encoded into random 10,000-length vectors (usually binary or bipolar) that are combined together to represent complex concepts. The length of the vectors is usually in the order of 10,000 to guarantee their quasi-orthogonality.

There are usually three types of arithmetic operators that can be applied to combine vectors: bundling, binding, and permutation. They all have unique properties that must be taken into account when used ({% cite kanerva2009hyperdimensional %}).

For what concerns the bundling operator: (i) the resulting vector is similar to the input vectors, (ii) the more vectors are involved in bundling, the harder it is to determine the component vectors, (iii) if several copies of any vector are included in bundling, the resulting vector is closer to the dominant vector than to the other components.

The binding operation is instead: (i) invertible (unbinding), (ii) it distributes over bundling, (iii) it preserves the distance, (iv) and the resulting vector is dissimilar to the input vectors.

Finally, the permutation operation is: (i) invertible, (ii) it distributes over bundling and any elementwise operation, (iii) it rpeserves the distance, and (iv) the resulting vector is dissimilar to the input vectors.

The set of vectors and arithmetic operators used for combining vectors is called vector-symbolic architecture.

`chopin2` implements a supervised classification model that works by encoding every observation in a dataset into high-dimensional vectors by combining values under their features. The model is built by collapsing the vector representation of the observations in the training set into a vector for each class of observations (e.g., in the case of the datasets retrieved in the previous step, classes are `CRC` and `control`). The classification model os then tested by computing the inner product between the vector representations of the observations in the test set against the two classes vectors. Vectors are thus classified based on their closeness to the classes.

> <hands-on-title>Build a classification model with `chopin2`</hands-on-title>
> **chopin2** {% icon tool %}: Run **chopin2** with the following configuration:
>    * Select the `RA__ThomasAM__species.csv` dataset imported into the history during the **Get the data** step;
>    * Set a dimensionality of 10,000;
>    * Set 100 levels;
>    * Set 10 retraining iterations;
>    * Set the number of folds to 5 to cross-validate the model;
>    * Finally, disable the feature selection.
>
> The number of levels is the number of random vectors that `chopin2` generates for encoding encode data. This must be carefully selected since the accuracy of the resulting model is strictly correlated to this parameter ({% cite Cumbo2020 %}).
>
> The tool will end up creating a `summary` file in your history containing a few basic information about the generated classification model and, more importantly, it's accuracy.
>
> We suggest you to repeat the same steps for building a classification model on the binary presence/absence profiles in `BIN__ThomasAM__species.csv`. In this case there is no need to use so many levels since the dataset contains only two possible values (i.e., 0 and 1). Thus, the number of levels must be changed to 2.
>
>    > <question-title></question-title>
>    >
>    > Compare the `summary` files generate by `chopin2` as the result of building a classification model over the `RA__ThomasAM__species.csv` and `BIN__ThomasAM__species.csv` datasets.
>    > 1. Is there a difference in the accuracy of the two models?
>    > 2. Does the difference in the number of levels have an impact on the running time and the final accuracy of the models?
>    {: .question}
{: .hands_on}

# Feature selection

`chopin2` also implements a feature selection method based on the backward variable elimination strategy. It means that the tool will produce a classification model starting with the whole set of features in the dataset and iteratively remove those features that do not contribute to the accuracy of the model.

Be aware that this specific type of feature selection method could lead to the generation of thousands of classification models in order to determine the best features and discard those ones that do not significantly contribute to a good accuracy. However, despite the huge amount of computational resources required to run the algorithm, it can easily handle datasets with massive amounts of features.

> <hands-on-title>Identify the best features</hands-on-title>
> **chopin2** {% icon tool %}: Run **chopin2** with the following configuration:
>    * Select the `BIN__ThomasAM__species.csv` dataset imported into the history during the **Get the data** step;
>    * Set a dimensionality of 10,000;
>    * Set 2 levels;
>    * Set 10 retraining iterations;
>    * Set the number of folds to 5 to cross-validate the model;
>    * Finally, enable the feature selection.
>
> We are going to focus on the `BIN__ThomasAM__species.csv` dataset only this time since `chopin2` produced a classification model with a better accuracy compared to that of the model built over the `RA__ThomasAM__species.csv` dataset.
>
> The tool will end up creating a `summary` file in your history containing a line for each step of the backward variable elimination method. Every line corresponds to a classification model. Please note that the last column reports the features excluded in each step.
>
> A `selection` file will also appear in your history with the list of best features selected by `chopin2`.
>
>    > <question-title></question-title>
>    >
>    > Look at the list of microbial species (features) selected by `chopin2` reported in the `selection` file. Can you recognise any microbial species known to be linked in some way to the genesis and development of colorectal cancer?
>    {: .question}
{: .hands_on}

> <details-title>Feature selection with Hyperdimensional Computing</details-title>
> `chopin2` is the first tool that implements a feature selection method based on the hyperdimensional computing paradigm.
{: .details}
