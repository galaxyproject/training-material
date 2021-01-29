---
layout: tutorial_hands_on

title: Introduction to recurrent neural networks (RNN)
zenodo_link: ''
questions:
- What is a recurrent neural network (RNN)?
- What are some applications of RNN?
objectives:
- Understand the difference between feedforward neural networks (FNN) and RNN
- Learn various RNN types and acrchitectures
- Solve a sentiment analysis problem on IMDB movie review dataset using RNN in Galaxy
time_estimation: 2H
contributors:
- Kaivan Kamali

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

Artificial neural networks are a machine learning discipline roughly inspired by how neurons in a 
human brain work. In the past decade, there has been a huge resurgence of neural networks thanks 
to the vast availability of training data and enormous increases in computing capacity (Successfully 
training complex neural networks in some domains requires lots of data and compute capacity). There 
are various types of neural networks (Feedforward, recurrent, etc). In this tutorial, we discuss 
the recurrent neural networks and explain how they differ from feedforward variants. We also describe 
various RNN architectures and solve a sentiment analysis problem using RNN in Galaxy.       

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Review of feedforward neural networks (FNN)

In feedforward neural networks (FNN) all the fields of a training example are presented
to the network at once, after which the the network generates an output. For example, a
lung X-ray image is passed to a FNN, and the network predicts tumor or no tumor. By contrast,
in RNNs the training example fields are presented to the network one at a time. For example,
a sequence of English words is passed to an RNN, one at a time, and the network generates a
sequence of Persian words, one at a time. RNNs handle sequential data, whether its temporal or ordinal.

## Single layer FNN

![Alternative text](../../images/FFNN_no_hidden.png "Single layer feedforward neural network")

Figure 1 shows a single layer FNN, where the input is 3 dimensional. Each input field is multiplied by a
weight. Afterwards, the results are summed up, along with a bias, and passed to an activation function.

![Alternative text](../../images/activation.gif "Activation of the output neuron o1. Activation function f could be Sigmoid, Tanh, ReLU, etc.")

The activation function can have many forms (sigmoid, tanh, ReLU, linear, step function, sign function, etc.).
Output layer neurons usually have sigmoid or tanh functions.

![Alternative text](../../images/sigmoid.gif "Sigmoid activation function")

## FFN with hidden layer

Minsky and Papert published a book in 1969 that showed that a single layer FNN cannot solve problems in which 
the data is not linearly separable (such as the XOR problem). Adding one (or more) hidden layers to FNN enables 
it solve problems in which data is non-linearly separable. In theory a FNN with one hidden layer can represent 
any function (add reference), although in practice training such a model is very difficult (if not impossible), 
hence, we usually add multiple hidden layers to solve complex problems.   

![Alternative text](../../images/FFNN.png "Feedforward neural network with a hidden layer. Biases to hidden/output layer neurons are ommitted for clarity")

## Learning algorithm 

The learning algorithm incrementally tweaks the network weights, so that the network error on the training 
set is minimized. Training set is composed of many training examples. Each training set is fed to the network 
and the network output is compared to the expected output. We need to define a **loss function** to objectively 
measure how much the predicted output is off of the expected output. For classification problems we use the 
**cross entropy** loss function.     

![Alternative text](../../images/CrossEntropy.gif "Cross entropy loss function")

Loss function is calculated for each training example. The average of the calculated loss functions over all training 
examples in the training set is the **Cost function**. The goal of the learning algorithm is to minimize the cost 
function. The cost function is a function of network weightsand biases of all neurons in all layers. The **backpropagation** 
learning algorithm iteratively computes the gradient of cost function relative to each weight and bias, then updates the weights 
and biases in the opposite direction of the gradient, to find the local minimum.

![Alternative text](../../images/CostFunction.gif "Cross entropy cost function")
  
# Recurrent neural networks

## Possible inputs/outputs
## RNN architectures

Section and subsection titles will be displayed in the tutorial index on the left side of
the page, so try to make them informative and concise!

# Text representation schemes

## Bag of words
## Text frequency inverse document frequency
## Word2Vec


You may want to cite some publications; this can be done by adding citations to the
bibliography file (`tutorial.bib` file next to your `tutorial.md` file). These citations
must be in bibtex format. If you have the DOI for the paper you wish to cite, you can
get the corresponding bibtex entry using [doi2bib.org](https://doi2bib.org).

With the example you will find in the `tutorial.bib` file, you can add a citation to
this article here in your tutorial like this:
{% raw %} `{% cite Batut2018 %}`{% endraw %}.
This will be rendered like this: {% cite Batut2018 %}, and links to a
[bibliography section](#bibliography) which will automatically be created at the end of the
tutorial.


**Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

# Hands-on Sections
Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]() or from the shared data library
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}

# Title of the section usually corresponding to a big step in the analysis

It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:


## Sub-step with **Create deep learning model**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Create deep learning model** {% icon tool %} with the following parameters:
>    - *"Choose a building mode"*: `Build a training model`
>        - In *"Compile Parameters"*:
>            - *"Select an optimizer"*: `Adam - Adam optimizer `
>            - *"Select metrics"*: ``
>        - In *"Fit Parameters"*:
>            - *"epochs"*: `2`
>            - *"batch_size"*: `128`
>            - In *"callback"*:
>                - {% icon param-repeat %} *"Insert callback"*
>                    - *"Choose a callback"*: `None`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Deep learning training and evaluation**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Deep learning training and evaluation** {% icon tool %} with the following parameters:
>    - *"Select a scheme"*: `Train and Validate`
>        - In *"Validation holdout"*:
>            - *"Select the splitting method"*: `ShuffleSplit`
>        - In *"Metrics for evaluation"*:
>            - *"Select the primary metric (scoring):"*: `default with estimator`
>    - *"Select input type:"*: `tabular data`
>        - *"Choose how to select data by column:"*: `All columns`
>        - *"Choose how to select data by column:"*: `All columns`
>    - *"Save the fitted model"*: ``
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Model Prediction**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Model Prediction** {% icon tool %} with the following parameters:
>    - *"Select input data type for prediction"*: `tabular data`
>        - *"Choose how to select data by column:"*: `All columns`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Machine Learning Visualization Extension**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Machine Learning Visualization Extension** {% icon tool %} with the following parameters:
>    - *"Select a plotting type"*: `Confusion matrix for classes`
>        - *"Choose how to select data by column:"*: `All columns`
>        - *"Does the dataset contain header:"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
