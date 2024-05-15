---
layout: tutorial_hands_on

title: Fine tune large protein AI model using HuggingFace
zenodo_link: https://doi.org/10.5281/zenodo.10986248
questions:
- How to load large protein AI models?
- How to fine-tune such models on downstream tasks such as post-translational site prediction?
objectives:
- Learn to load and use large protein models from HuggingFace
- Learn to fine-tune them on specific tasks such as predicting dephosphorylation sites
requirements:
  -
    type: internal
    topic_name: galaxy-interface
    tutorials:
      - jupyterlab
  -
    type: internal
    topic_name: statistics
    tutorials:
      - fine_tune_large_protein_models
time_estimation: 1H
tags:
- interactive-tools
- machine-learning
- deep-learning
- jupyter-lab
- fine-tuning
- dephosphorylation-site-prediction
contributors:
- sheetkulkarni01
- anuprulez

---

The advent of large language models has transformed the field of natural language processing, enabling machines to comprehend and generate human-like language with unprecedented accuracy. Pre-trained language models, such as BERT, RoBERTa, and their variants, have achieved state-of-the-art results on a wide range of tasks, from sentiment analysis and question answering to language translation and text classification. Moreover, the emergence of transformer-based models, such as Generative Pre-trained Transformer (GPT) and its variants, has enabled the creation of highly advanced language models that can generate coherent and context-specific text. The latest iteration of these models, ChatGPT, has taken the concept of conversational AI to new heights, allowing users to engage in natural-sounding conversations with machines. However, despite their impressive capabilities, these models are not yet perfect, and their performance can be significantly improved through fine-tuning. Fine-tuning involves adapting the pre-trained model to a specific task or domain by adjusting its parameters to optimize its performance on a target dataset. This process allows the model to learn task-specific features and relationships that may not be captured by the pre-trained model alone, resulting in highly accurate and specialized language models that can be applied to a wide range of applications. In this article, we will delve into the world of fine-tuning large language models, exploring the benefits and challenges of this approach, as well as the various techniques and strategies that can be employed to achieve optimal results. Protein large language models (LLMs) represent a groundbreaking advancement in bioinformatics, leveraging the power of deep learning to understand and predict the behavior of proteins at an unprecedented scale. These models, exemplified by the ProtTrans suite, are inspired by natural language processing (NLP) techniques, applying similar methodologies to biological sequences. ProtTrans models, including BERT and T5 adaptations, are trained on vast datasets of protein sequences, enabling them to capture the complex patterns and functions encoded within amino acid sequences. By interpreting these sequences much like language, protein LLMs offer transformative potential in drug discovery, disease understanding, and synthetic biology, bridging the gap between computational predictions and experimental biology.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Open JupyterLab in Galaxy Europe


### Open JupyterLab

> <hands-on-title>GPU-enabled Interactive Jupyter Notebook for Machine Learning</hands-on-title>
>
> - {% tool [GPU-enabled Interactive Jupyter Notebook for Machine Learning](interactive_tool_ml_jupyter_notebook) %}
>    - *"Do you already have a notebook?"*: `Start with default notebooks`
>    - Click *"Run Tool"*
>
>    > <comment-title></comment-title>
>    >  If you do not have access to this resource in Galaxy Europe, please apply for it at: [Access GPU-JupyterLab](http://usegalaxy.eu/gpu-request). It may take a day or two to receive access.
>    >
>    {: .comment}
{: .hands_on}


### Fetch notebook and protein sequences

> <hands-on-title>Fetch data from Zenodo</hands-on-title>
>
> 1. Create a new folder named `fine-tuning` alongside other folders such as "data", "outputs", "elyra" or you can use your favourite folder name.
> 2. Inside the created folder, clone a code repository by clicking on "Git" icon.
> 3. In the shown popup, provide the repository path as shown below and then, click on "clone":
>    ```
>    https://github.com/anuprulez/fine-tune-protTrans-repository
>    ```
>
{: .hands_on}


## Fine-tuning notebook

From the cloned repository, open the `fine-tune-protTrans-dephophorylation.ipynb` notebook. The notebook contains all the necessary script for processing protein sequences, creating and configuring protein large language model, training it on the protein sequences and evaluating them on the test protein sequences and visualising results. Let's look at these key steps of fine-tuning.

#### Install necessary Python packages


#### Fetch and split data

#### Define configurations for LoRA and transformer (T5) model

#### Create T5 model

#### Create model training method

#### Train model

#### Analyse results 

## Conclusion


