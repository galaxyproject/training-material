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


## Image segmentation
In an image classification task, a label is assigned to each image. For example, an image containing a cat gets a label as "cat" ("cat" is given an integer to be used in programs). In an image segmentation task, each pixel of an image gets a label. For example, in an image containing a cat, all the pixels occupied by cat in the image are denoted by violet color (see image below [https://www.tensorflow.org/tutorials/images/segmentation]). Image segmentation is widely used in many fields such as medical imaging to find out abnormal regions in tissue images, self-driving cars to detect different objects on highways and so on.

![Image segmentation for cat](../../images/cat_segmentation.png "Image segmentation for cat. The outline of the violet section denoted by yellow color is cat's mask.")


## COVID CT scan
Segmenting lung CT scans to locate infected regions has been actively proposed to augment the RT-PCR testing for initial screening of COVID infection in humans. Deep learning has been used to predict these regions with high accuracy in works such as "Inf-net: Automatic COVID-19 lung infection segmentation from CT images" [https://ieeexplore.ieee.org/document/9098956], "JCS: An explainable COVID-19 diagnosis system by joint classification and segmentation" [https://arxiv.org/abs/2004.07054] and so on. In the image shown below, the differences between the CT scans of a normal person and suffering from COVID can be seen. The regions marked by white patches are infected by COVID [https://www.sciencedirect.com/science/article/pii/S2666990021000069]. 

![Normal and Covid CT scans](../../images/normal_covid_ct_scans.jpg "(Left) Lungs CT scan for normal person and (right) lungs CT scan for person having COVID-19.")

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
