---
layout: tutorial_hands_on

title: Introduction to recurrent neural networks (RNN)
zenodo_link: https://zenodo.org/record/4477881
questions:
- What is a recurrent neural network (RNN)?
- What are some applications of RNN?
objectives:
- Understand the difference between feedforward neural networks (FNN) and RNN
- Learn various RNN types and acrchitectures
- Solve a sentiment analysis problem on IMDB movie review dataset using RNN in Galaxy
requirements:
-
  type: internal
  topic_name: statistics
  tutorials:
    - intro_deep_learning
time_estimation: 2H
contributors:
- kxk302

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

## Multi-layer FFN

Minsky and Papert showed that a single layer FNN cannot solve problems in which the data is not linearly separable, 
such as the XOR problem ({% cite Newell780 %}). Adding one (or more) hidden layers to FNN enables it solve problems 
in which data is non-linearly separable. Per Universal Approximation Theorem, a FNN with one hidden layer can represent 
any function ({% cite Cybenko1989 %}), although in practice training such a model is very difficult (if not impossible), 
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

Unlike FNN, in RNN the output of the network at time t is used as network input at time t+1. RNN handle sequential data (e.g. temporal or ordinal). 

## Possible RNN inputs/outputs

There are 4 possible input/output combinations for RNN and each have a specific application. One-to-one is basically a FNN. One-to-many, 
where we have one input and a variable number of output. One example application is image captioning, where a single image is provided 
as input and a variable number of words (which caption the image) is returned as output (See Figure 7).   

![Alternative text](../../images/RNN_1_to_n.png "One-to-many RNN")

Many-to-one RNN, on the other hand, have a variable number of inputs and a single output. One example application is document sentiment 
classification, where a variable number of words in a document are presented as input, and a single output predicts whether the document
has a positive or negative sentiment regarding a topic (See Figure 8).

![Alternative text](../../images/RNN_n_to_1.png "Many-to-one RNN")

There are two types of many-to-many RNN. One in which the number of inputs and outputs match, e.g., in labeling the video frames the number 
of frames matches the number of labels, and the other in which the number of inputs and outputs do not match, e.g., in language translation 
we pass in n words in English and get m words in Italian (See Figure 9).

![Alternative text](../../images/RNN_n_to_m.png "Many-to-many RNN")

## RNN architectures

Mainly, there are three types of RNN: 1) Vanilla RNN, 2) LSTM, and 3) GRU. A Vanilla RNN, simply combines the state information from the 
previous timestamp with the input from the current timestamp to generate the state information for current timestamp. The problem with Vanilla 
RNN is that training deep RNN networks is impossible due to the **vanishing gradient** problem. Basically, starting from the output layer, 
in order to determine weights/biases updates, we need to calculate the derivative of the loss function relative to the layers input, which is 
usually a small number. This is not a problem for the output layer, but for the previous layers, this process must be repeated recursively, 
resulting in very small updates in weights/biases of the initial layers of the RNN, lting the learning process.

LSTM and GRU are two RNN architectures that address vanishing gradient problem. Full description of LSTM/GRU is beyond the scope of this 
tutorial (Please refer to ref1 and ref2), but in a nutshell both LSTM and GRU used **gates** such that the weights/biases updates in previous 
layers are calculated via a series of additions (not multiplications). Hence, these architectures can learn even when the RNN has hundreds or 
thousands of layers. 

# Text representation schemes

In this tutorial we perform sentiment analysis on IMDB movie reviews dataset. We train our RNN on the training dataset, which is made up of 25000 
movie reviews, some positive and some negative. We then test our RNN on the test set, which is also made up of 25000 movie reviews, again some 
positive and some negative. The training and test sets have no overlap. Since we are dealing with text data, its a good idea to review various 
mechanisms for representing text data.
 
## Bag of words and TF-IDF

If you don't care about the order of the words in a document, you can use bag of words (BoW) or text frequency inverse document frequency (TF-IDF).
In these models we have a 2 dimensional array. The rows represent the documents (in our example, the movie reviews) and the columns
represent the words in our voabulary (all the unique words in all the documents). If a word is not present in a document, we have a zero 
at the corresponding row and column as the entry. If a word is present in the document, we have a one as the entry -- Alternatively, we could use 
the word count or frequency.

![Alternative text](../../images/BoW.png "Bag of words (BoW) representation")

Suppose we have the following 2 documents: 1) Magic passed the basketball to Kareem, and 2) Lebron stole the basketball from Curry. The BoW 
representation of these documents is given in Figure 10. 

BoW's advantage is its simplicity, yet it does not take into account the rarity of a word across documents, which unlike common words are
important for document classification. 

In TF-IDF, similar to BoW we have an entry for each document-word pair. In TD-IDF, the entry is the product of 1) Text frequency, the 
frequency of a word in a document, and 2) Inverse document frequency, the inverse of the number of documents that have the word divided 
by the total number of documents (we usually use logarithm of the IDF).

TF-IDF takes into account the rarity of a word across documents, but like BoW does not capture word order or word meaning in documents. BoW 
and TF-IDF are suitable representations for when word order is not important. They are used in document classification problems like spam detection.

## One hot encoding (OHE)

OHE is a technique to convert categorical variables such as words into a vector. Suppose our vocabulary has 3 words: orange, apple, banana. 
Each word for this vocabulary is represented by a vector of size 3. Orange is represented by a vector whose first element is 1 and other 
elements are 0; Apple is represented by a vector whose second element is 1 and other elements are 0; And banana is represented by a 
vector whose third element is 1 and other elements are 0. As you can see only one element in the vector is 1 and the rest are 0's. The same 
concept applies if the size of the vocabulary is N.    

![Alternative text](../../images/OHE.gif "One hot encoding (OHE) representation")

The problem with OHE is that for very large vocabulary sizes (say, 100,000 words) it requires tremendous amount of storage. Also, it has no 
concept of word similarity.   

## Word2Vec

In Word2Vec, each word is represented as an n dimensional vector (n being much smaller than vocabulary size), such that the words that have 
similar meanings are closer to each other in the vector space, and words that don't have a similar meaning are farther apart. Words are 
considered to have a similar meaning if they co-occur often in documents. There are 2 Word2Vec architectures, one that predicts the probability 
of a word given the surrounding words (Continous BOW), and one that given a word predicts the probability of the surrounding words (Continous skip-gram).

In this tutorial, we find an n dimensional representation of the IMDB movie review words, not based on word meanings, but based on how they
improve the sentiment classification task.    

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
