---
layout: tutorial_hands_on

title: Introduction to deep learning
zenodo_link: https://zenodo.org/record/3706539#.XmjDYHVKg5k
questions:
- What is deep learning and neural networks?
- Why is it useful?
- How to create a neural network architecture for classification?
objectives:
- Learn about how to create an end-to-end neural network architecture
- Learn about Galaxy deep learning tools
- Learn how to interpret predictions
key_points:
- Multiple tools to constitute a neural network architecture
- Interpretation of predictions using visualisation tools
time_estimation: 1H
contributors:
- anuprulez
---

## Introduction

### Deep learning and neural networks
Deep learning is a branch of artificial intelligence which recognises patterns in large volumes of data. Using these learned patterns on the existing data, new data can be categorised. These patterns in data are learned by a computational model based on multiple architectures of neural networks. A neural network is a web of artificial neurons which is also called processing units. The idea of neural networks is inspired from mamalian cerebral cortex where neuronal circuits are used to learn. A neural network is structured into multiple layers where each layer contains several neurons. The neurons from adjacent layers are interconnected allowing the exchange of information. An artificial neuron is shown in Figure 1. The neuron, shown in orange, takes inputs (x1 and x2) and computes output (y). The entities w1 and w2 are the weights of the connections (between inputs and neuron). The weights and inputs are combined following the basic principles of mathematics.

x1 = x1 * w1

x2 = x2 * w2

![data](../../images/neuron.svg "An artificial neuron.")

The weights denote the significance of a particular input to produce the observed output. When it is large, the input is significant and whnen small, the input less significant to produce the output. These weights can be initialised randomly and they are modified over the course of learning by a neural network. Using the updated inputs (as shown in above equations), the output is computed using:

y = f(x1 + x2)

where *f* is an activation function. An activation function is a mathematical function which translates the combination of inputs to an output. The choices of these functions are many - sigmoid, linear, tanh, ReLU and so on. For example, sigmoid is:

f(x) = 1 / (1 + exp(-x))

The above equation will return a real number between 0 and 1.

ReLU is:

f(x) = max(0, x)

Neurons make the building blocks of a neural network and are arranged in several layers. A usual neural network will look like as shown in Figure 2.

![data](../../images/neural_network.svg "A neural network.")



### Relevance of deep learning in Bioinformatics

## Get training and test datasets

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
>
>    {% include snippets/create_new_history.md %}
>
> 2. Import the files from [Zenodo](https://zenodo.org/record/3706539#.XmjDYHVKg5k)
>
>    ```
>    https://zenodo.org/record/3706539/files/X_test.tsv
>    https://zenodo.org/record/3706539/files/X_train.tsv
>    https://zenodo.org/record/3706539/files/y_test.tsv
>    https://zenodo.org/record/3706539/files/y_train.tsv
>    ```
>
>    {% include snippets/import_via_link.md %}
>
> 3. Check that the datatype is `tabular`.
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
{: .hands_on}


## Neural network architecture

### Create architecture: Choose layers

> ### {% icon hands_on %} Hands-on: Create a deep learning model architecture using Keras
>
> 1. **Create a deep learning model architecture using Keras** {% icon tool %} with the following parameters:
>    - *"Select keras model type"*: `Sequential`
>    - *"input_shape"*: `(7129, )`
>       
>    - In *"LAYER"*:
>        - In *"1: LAYER"*:
>            - *"Choose the type of layer"*: `Core -- Dense`
>                - *"units"*: `16`
>                - *"Activation functions"*: `elu`
>        - In *"2: LAYER"*:
>            - *"Choose the type of layer"*: `Core -- Dense`
>                - *"units"*: `16`
>                - *"Activation functions"*: `elu`
>        - In *"3: LAYER"*:
>            - *"Choose the type of layer"*: `Core -- Dense`
>                - *"units"*: `1`
>                - *"Activation functions"*: `sigmoid`
> 
>
{: .hands_on}

The tool returns a JSON output file containing information about the nueral network layers and their attributes like their types, number of units they have and their activation functions.


### Create architecture: Add training parameters

> ### {% icon hands_on %} Hands-on: Create deep learning model with an optimizer, loss function and fit parameters
>
> 1. **Create deep learning model with an optimizer, loss function and fit parameters** {% icon tool %} with the following parameters:
>    - *"Choose a building mode"*: `Build a training model`
>    - *"Select the dataset containing model configurations (JSON)"*: `Keras model config` (output of **Create a deep learning model architecture using Keras** {% icon tool %})
>    - *"Do classification or regression?"*: `KerasGClassifier`
>    - In *"Compile Parameters"*:
>        - *"Select a loss function"*: `binary_crossentropy`
>        - *"Select an optimizer"*: `RMSprop - RMSProp optimizer`
>    - In *"Fit Parameters"*:
>        - *"epochs"*: 10`
>        - *"batch_size"*: 4
>
{: .hands_on}

The tool returns a zipped file containing the object of classifier. The classifier object will be used for training.

### Deep learning training

> ### {% icon hands_on %} Hands-on: Deep learning training and evaluation conduct deep training and evaluation either implicitly or explicitly
>
> 1. **Deep learning training and evaluation conduct deep training and evaluation** {% icon tool %} with the following parameters:
>    - *"Select a scheme"*: `Train and validate`
>    - *"Choose the dataset containing pipeline/estimator object"*: `Keras model builder` (output of **Create deep learning model** {% icon tool %})
>    - *"Select input type"*: `tabular data`
>    - *"Training samples dataset"*: `X_train.tsv`
>        - *"Does the dataset contain header"*: `Yes`
>        - *"Choose how to select data by column"*: `All columns`
>    - *"Dataset containing class labels or target values"*: `y_train.tsv`
>        - *"Does the dataset contain header"*: `Yes`
>        - *"Choose how to select data by column"*: `All columns`
>
{: .hands_on}

The tool returns 3 files - a `tabular` file containing output (accuracy of cross-validation) of training, a `zipped` file with the trained model and an `H5` file containing the weights of neural network layers.

### Prediction on test data

> ### {% icon hands_on %} Hands-on: Model Prediction predicts on new data using a preffited model
> 
> 1. **Model Prediction predicts on new data using a preffited model** {% icon tool %} with the following parameters:
>    - *"Choose the dataset containing pipeline/estimator object"*: `Fitted estimator or estimator skeleton` (output of **Deep learning training and evaluation** {% icon tool %})
>    - *"Choose the dataset containing weights for the estimator above"*: `Weights trained` (output of **Create deep learning model** {% icon tool %})
>    - *"Select invocation method"*: `predict`
>    - *"Select input data type for prediction"*: `tabular data`
>        - *"Training samples dataset"*: `X_test`
>        - *"Does the dataset contain header"*: `Yes`
>        - *"Choose how to select data by column"*: `All columns`
> 
>
{: .hands_on}

The tool returns the predictions for the test data in a tabular dataset.

## Visualisation

> ### {% icon hands_on %} Hands-on: Machine Learning Visualization Extension includes several types of plotting for machine learning
> 
> 1. **Machine Learning Visualization Extension includes several types of plotting for machine learning** {% icon tool %} with the following parameters:
>    - *"Select a plotting type"*: `Confusion matrix for classes` 
>    - *"Select dataset containing true labels"*: `y_test.tsv`
>    - *"Does the dataset contain header"*: `Yes`
>    - *"Choose how to select data by column"*: `All columns`
>    - *"Select dataset containing predicted labels"*: `Model prediction` (output of **Model Prediction predicts on new data using a preffited model** {% icon tool %})
> 
>
{: .hands_on}

![data](../../images/confusion_matrix_dl.png "Confusion matrix for true and predicted classes")

## Summary


## Conclusion


{:.no_toc}

