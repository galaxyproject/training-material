---
layout: tutorial_hands_on

title: COVID CT scan images segmentation prediction using uNet neural network in GPU-powered Jupyterlab
zenodo_link: https://zenodo.org/record/6091361#.Ygu4gIzMI5k
questions:
- How to use Jupyterlab and it several features?
- How to use it for creating input datasets and writing artificial intelligence (AI) algorithms?
objectives:
- Learn to use Jupyterlab - an online Python editor designed for developing AI algorithms
- Explore several of its features such as Git integration, workflow of jupyter notebook, integration to Galaxy
- Develop AI algorithms using scikit-learn and tensorflow
- Send long-running jobs to Galaxy's cluster and save results in its history
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
      - intro_deep_learning
time_estimation: 1H
contributors:
- anuprulez

---

# Introduction
Jupyterlab [https://jupyterlab.readthedocs.io/en/stable/] is a popular integrated development environment (IDE) for a variety of tasks in data science such as prototyping analyses, creating meaningful plots, data manipulation and preprocessing. Python is one of mostly used languages in such environment. Given the usefulness of Jupyterlab, more importantly in online platforms, a robust Jupyterlab notebook application has been developed that is powered by GPU acceleration and contains numerous packages such as Pandas, Numpy, Scipy, Scikit-learn, Tensorflow, ONNX for modern data science. It has been developed as an interactive Galaxy tool that runs on an isolated docker container [https://github.com/anuprulez/ml-jupyter-notebook]. The docker container has been built using "jupyter/tensorflow-notebook:tensorflow-2.6.0" as the base container. Moreover, With the use of Bioblend [https://bioblend.readthedocs.io/], a Galaxy tool [https://github.com/bgruening/galaxytools/pull/1157] can be executed to make use of Galaxy remote job processing for long-running deep learning training and the finished datasets (such a trained models, tabular files, ...) are saved in a Galaxy history. 


## Image segmentation and COVID CT scan


https://www.sciencedirect.com/science/article/pii/S2666990021000069

### Unet neural network

{:.no_toc}

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


## Custom Jupyterlab features features programs Git version control Galaxy 

### AI programs on GPU(s)

### Git version control

### Shareable AI models using ONNX

### Accessible Galaxy tools via Bioblend

### Miscellaneous - Elyra AI - workflow of notebooks, GPU utilization dashboards,


## Get data

## Dataset description

## Neural network architecture description

## Segmentation prediction using uNet in Jupyterlab

### Open Jupyterlab editor

### Pull code

### Attach dataset

### Train model in the notebook

### Train model remotely and make inference

### Predict masks using trained model

# Conclusion
{:.no_toc}
