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

{:.no_toc}

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Introduction
Jupyterlab [https://jupyterlab.readthedocs.io/en/stable/] is a popular integrated development environment (IDE) for a variety of tasks in data science such as prototyping analyses, creating meaningful plots, data manipulation and preprocessing. Python is one of mostly used languages in such environment. Given the usefulness of Jupyterlab, more importantly in online platforms, a robust Jupyterlab notebook application has been developed that is powered by GPU acceleration and contains numerous packages such as Pandas, Numpy, Scipy, Scikit-learn, Tensorflow, ONNX for modern data science. It has been developed as an interactive Galaxy tool that runs on an isolated docker container [https://github.com/anuprulez/ml-jupyter-notebook]. The docker container has been built using "jupyter/tensorflow-notebook:tensorflow-2.6.0" as the base container. Moreover, With the use of Bioblend [https://bioblend.readthedocs.io/], a Galaxy tool [https://github.com/bgruening/galaxytools/pull/1157] can be executed to make use of Galaxy remote job processing for long-running deep learning training and the finished datasets (such a trained models, tabular files, ...) are saved in a Galaxy history. 


## Image segmentation
In an image classification task, a label is assigned to each image. For example, an image containing a cat gets a label as "cat" ("cat" is given an integer to be used in programs). In an image segmentation task, each pixel of an image gets a label. For example, in an image containing a cat, all the pixels occupied by cat in the image are denoted by violet color (see image below [https://www.tensorflow.org/tutorials/images/segmentation]). Image segmentation is widely used in many fields such as medical imaging to find out abnormal regions in tissue images, self-driving cars to detect different objects on highways and so on.

![Image segmentation for cat](../../images/cat_segmentation.png "Image segmentation for cat. The outline of the violet section denoted by yellow color is cat's mask.")


## COVID CT scan
Segmenting lung CT scans to locate infected regions has been actively proposed to augment the RT-PCR testing for initial screening of COVID infection in humans. Deep learning has been used to predict these regions with high accuracy in works such as "Inf-net: Automatic COVID-19 lung infection segmentation from CT images" [https://ieeexplore.ieee.org/document/9098956], "JCS: An explainable COVID-19 diagnosis system by joint classification and segmentation" [https://arxiv.org/abs/2004.07054] and so on. In the image shown below, the differences between the CT scans of a normal person and suffering from COVID can be seen. The regions marked by white patches are infected by COVID [https://www.sciencedirect.com/science/article/pii/S2666990021000069].

![Normal and Covid CT scans](../../images/normal_covid_ct_scans.jpg "(Left) Lungs CT scan for normal person and (right) lungs CT scan for person having COVID-19.")


## CT scans and masks
The picture below shows CT scans of infected lungs (top row) and borders around the infected regions have been drawn with red color (middle row). In the last row respective masks, from the CT scans, have been taken out denoting the infected regions. These masks are the "segmented" regions from the corresponding CT scans. While creating the dataset for training the deep learning model (Unet), the CT scans are the data points and its "known" respective masks of infected regions become their labels [https://www.sciencedirect.com/science/article/pii/S2666990021000069]. 

![Masks](../../images/covid_ct_scan_masks.png "Picture shows CT scans and its masks of infected region (drawn with red color in images of the middle row).")

### Unet neural network
Unet neural network (Unet) is widely used for segmentation tasks in images. The name "Unet" resembles the shape of its architecture as "U" (see image below). It has two parts - encoder and decoder. The left half of the "U" shape is the encoder that learns the feature representation (downsampling) at lower level. The right half is the decoder that maps the low resolution feature representation on the higher resolution (upsampling) pixel space. For further reading, please refer to: [https://arxiv.org/pdf/1505.04597.pdf]

![Unet architecture](../../images/Unet.jpg "Architecture of Unet neural network for image segmentation")

In this tutorial, we will use the dataset of CTs scans and their respective masks to train a Unet model. The model learns to map the infected regions in the CT scans to their masks. For prediction, the trained model is given CT scans and it predicts infected regions. For this experiment, we will use Jupyterlab for pulling the notebooks that contain training and prediction scripts. The data (images and also the trained model) required for this notebook can be downloaded from [https://zenodo.org/record/6091361#.Ygu4gIzMI5k]. The model can either be trained in the Jupyterlab or can be sent to Galaxy's cluster for remote processing. After remote processing, the generated datasets such as trained model become available in a new Galaxy history.


## Custom Jupyterlab features
Jupyterlab notebook has been augmented with several useful features that makes it ready-to-use for quick prototying of AI projects. Feature such as **available online** makes it really convenient to share it with other researchers and users. GPUs have accelerated AI research, especially deep leanring. Therefore, the backend of the Jupyterlab is powered by GPU to make long running AI training programs finish faster by parallelizing matrix multiplications. Galaxy tool for remote job processing also runs on GPU. Jupyterlab is also integrated with **Git version control** that makes it easy to pull, push and maintain codebase directly into the notebook. Repositories from Github can be easily clone, updated and maintained. In addition, a standard model format **Open Neural Network Exchange(ONNX)**, has been added to transform scikit-learn and tensorflow models to "onnx" files. These files can then be conveniently shared and used for inference. Galaxy also supports "onnx" file format to make it easier to create, save and share such models. Using Bioblend, the notebook can be connected to Galaxy and different histories and tools can be accessed. Using this features, Galaxy tools can be executed with the correct input datasets such as the Galaxy tool for processing long running training. Many notebooks can be created serving different purposes. These notebooks can be knit together to form one pipeline where each notebook transforms data taking a form of data from its previous notebook and pass on the transformed data to its next nextbook. This feature 


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
