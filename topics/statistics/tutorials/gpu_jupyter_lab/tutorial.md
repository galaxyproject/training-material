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

## Jupyterlab features
Jupyterlab notebook has been augmented with several useful features that makes it ready-to-use for quick prototying of AI projects. Feature such as **available online** makes it really convenient to share it with other researchers and users. GPUs have accelerated AI research, especially deep leanring. Therefore, the backend of the Jupyterlab is powered by GPU to make long running AI training programs finish faster by parallelizing matrix multiplications. Galaxy tool for remote job processing also runs on GPU. Jupyterlab is also integrated with **Git version control** that makes it easy to pull, push and maintain codebase directly into the notebook. Repositories from Github can be easily clone, updated and maintained. In addition, a standard model format **Open Neural Network Exchange(ONNX)**, has been added to transform scikit-learn and tensorflow models to "onnx" files. These files can then be conveniently shared and used for inference. Galaxy also supports "onnx" file format to make it easier to create, save and share such models. Using Bioblend, the notebook can be connected to Galaxy and different histories and tools can be accessed. Using this features, Galaxy tools can be executed with the correct input datasets such as the Galaxy tool for processing long running training. Many notebooks can be created serving different purposes. These notebooks can be knit together to form one pipeline where each notebook transforms data taking a different form of data from its previous notebook and pass on the transformed data to its next nextbook.  
An example workflow created using Elyra AI [https://elyra.readthedocs.io/en/stable/] can be found at [<<open notebook>>/elyra/METABRIC_ML.pipeline]. These pipelines can also be executed remotely on different cluster using its "runtimes" features that pull different docker containers. An example pipeline can be seen in the picture below.

![Elyra AI pipeline](../../images/elyra_ai_pipeline.png "An example ML pipeline created using Elyra AI").

There are many other features such as GPU utilization dashboards for monitoring the GPU usage and system memory utilization, voila for rendering output cells of a notebook in a separate tab hiding all code cells, interactive bqplots, and many more. There are several packages suited for performing ML tasks such as Open-CV and Scikit-Image for image processing, NiBabel package for processing images files.

## Security
Security benefits of executing code in isolated environments such as docker containers 


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

# Image segmentation using Jupyterlab
In this tutorial, we will a few features of Jupyterlab in Galaxy to create and train a deep learning model and predict segmented regions using the trained model from COVID CT scans.

## Open Jupyterlab editor

> ### {% icon hands_on %} Hands-on: GPU enabled Interactive Jupyter Notebook for Machine Learning
>
> - {% tool [Create a deep learning model architecture](interactive_tool_ml_jupyter_notebook) %}
>    - *"Do you already have a notebook?"*: `Start with a fresh notebook`
>    - Click *"Execute"*
{: .hands_on}

Now, we should wait for a few minutes until Galaxy creates the required compute environment for opening a new Jupyterlab. The progress can be checked by clicking on the "User>Active Interactive tools". On "Active Interactive Tools" page, there is a list of all open interactive tools. We should find our interactive tool by "GPU enabled Interactive Jupyter Notebook for Machine Learning" name. When the job info starts showing "running", then the name of the interactive tool gets associated with a link to the running Jupyterlab. On clicking this link, we can open the running Jupyterlab.


## Run image segmentation analysis
> ### {% icon hands_on %} Hands-on: Pull code
> Several features of Jupyterlab editor running in Galaxy can be found out by opening the "home_page.ipynb" notebook. Folder "notebooks" contain several notebooks showing the running usecases of different packages and features. To use Git version control for pulling any codebase from GitHub, following steps should be performed
> 1. Create a new folder named "covid_ct_segmentation"
> 2. Inside this folder, clone a code repository by clicking on "Git" icon as shown in the picture below
> 3. In the popup that asks the repository path, enter "https://github.com/anuprulez/gpu_jupyterlab_ct_image_segmentation" and click on "clone"
> 4. The repository "anuprulez/gpu_jupyterlab_ct_image_segmentation" gets immediately cloned
> 5. A few notebooks can be found inside "gpu_jupyterlab_ct_image_segmentation"
>    ![Clone repository](../../images/git_clone.png "Clone a code repository using Git").
>
{: .hands_on}

Now, we have all the notebooks available for performing image segmentation. The entire analysis can be completed in two ways - one by running all the notebooks in the Jupyterlab itself and another by sending the notebook for training to a remote Galaxy cluster invoking another Galaxy tool from the notebook. First, let's look at how we can run this analysis in the notebook itself.

> ### {% icon hands_on %} Hands-on: Run notebooks in Jupyterlab
>
> 1. Download and save datasets in the notebook using **"fetch_datasets.ipynb"** notebook. It will also create all the necessary folders. It downloads two datasets, one "h5" file containing many matrices as sub-datasets belonging to training data, training labels, validation data, validation labels, test data and test labels. These sub-datasets are stored in different variables after reading the original "h5" file once. For training we only need these datasets/matrices - training data, training labels, validation data and validation labels. The matrices, test data and test labels, are used for prediction. We use "h5" format for storing and retrieving datasets as all AI algorithms need input datasets in the form of matrices. Since, in the field of AI, there are many different types of datasets such as images, sequences, real numbers and so on. To coverge all these different forms of datasets, we use "h5" format. In an analysis, any type of dataset that can be used with AI algorithms can also be saved as "h5" file. For image segmentation tasks, we also saved all the input datasets/matrices to the deep learning model as sub-datasets in one "h5" file that can be easily created, stored and downloaded.
>  
> 2. Create and train a deep learning Unet model using **"create_model_and_train.ipynb"** notebook. This will read the input "h5" datasets containing all images and train the model after creating deep learning model architecture. This notebook first creates a deep learning architecture based on Unet including custom loss functions such as total variation and binary cross-entropy losses. After creating the deep learning architecture, all training datasets such as training data, training labels, validation data, validation labels are loaded from the combined "h5" file. In the next step, all the datasets and deep learning architecture are compiled together and training starts for 10 epochs (10 iteration over entire training dataset). The training is really fast as it runs of GPU and finishes in a few minutes created a trained model. In the last step, the trained model containing several files are converted to one "Onnx" file.
> 
> 3. Predict unseen masks using the trained model in **"predict_masks.ipynb"** notebook. This will read test datasets and trained model and make prediction and print prediction as an image. To predict masks of unseen CT scans, we use "predict_masks.ipynb" notebook. First, it reads test datasets from the combined "h5" file and then, loads the "Onnx" model. Using this model, it predicts masks of unseen CT scans and then plots the ground truth masks and predicted masks in one image (see below).
> ![True and predicted masks of CTs scans](../../images/true_pred_masks.png "True and predicted masks of CT scans")
>
> The above picture shows origin lungs CT scans in the first column. The second column shows the corresponding true infected regions in white. The third and fourth columns show the infected regions predicted with different loss methods, one is the binary cross-entropy loss and another is the combination of binary cross-entropy loss and total variation loss. Binary cross-entropy calculates loss between two corresponding pixels of true and predicted images. It measure how close two corresponding pixels are. But, total variation loss measures the amount of noise in predicted images as noisy regions usually show large variations in their neighbourhood. The noise in the predicted images gets minimized and the connectivity of predicted masks improves as well when minimizing this loss.   
>
{: .hands_on}


> ### {% icon hands_on %} Hands-on: 2. Run notebooks in Jupyterlab as well in remote Galaxy cluster
>
> 1. Download and save datasets in the notebook using **"fetch_datasets.ipynb"** notebook in the same way as before
>  
> 2. Execute **"create_model_and_train_remote.ipynb"** notebook in a remote cluster using another Galaxy tool by running **"run_remote_training.ipynb"** notebook
> 
> 3. Predict unseen masks using the trained model in **"predict_masks.ipynb"** notebook. This will read test datasets and trained model and make prediction and print prediction as an image. To predict masks of unseen CT scans, we use "predict_masks.ipynb" notebook. First, it reads test datasets from the combined "h5" file and then, loads the "Onnx" model. Using this model, it predicts masks of unseen CT scans and then plots the ground truth masks and predicted masks in one image (see below).
> ![True and predicted masks of CTs scans](../../images/true_pred_masks.png "True and predicted masks of CT scans")
>
> The above picture shows origin lungs CT scans in the first column. The second column shows the corresponding true infected regions in white. The third and fourth columns show the infected regions predicted with different loss methods, one is the binary cross-entropy loss and another is the combination of binary cross-entropy loss and total variation loss. Binary cross-entropy calculates loss between two corresponding pixels of true and predicted images. It measure how close two corresponding pixels are. But, total variation loss measures the amount of noise in predicted images as noisy regions usually show large variations in their neighbourhood. The noise in the predicted images gets minimized and the connectivity of predicted masks improves as well when minimizing this loss.   
>
{: .hands_on}



# Conclusion
In this tutorial, we have discussed several features of Jupyterlab editor in Galaxy and used a few of them to perform an image segmentation task by training a deep learning model. The ready-to-use infrastructure that we used here for developing an AI algorithm and training model proves very convenient for making quick prototypes, accelerated by GPU, in data science without thinking about the storage and compute resources. The usage of "h5" files for input datasets to a neural network solves the problem of varying data formats across multiple fields that use AI algorithms extensively for predictive tasks. Using onnx for creating, storing and sharing AI models makes it easier to handle such model as one compact file. Git integration makes it easy to maintain an entire repository in the notebook. The entire environment can be shared with other researchers only via sharing the URL. Bioblend connects this infrastructure to Galaxy which makes its histories, workflows and tools accessible directly via the notebook. Overall, this infrastructure would prove useful for researchers working to create quick prototypes of their application.



{:.no_toc}
