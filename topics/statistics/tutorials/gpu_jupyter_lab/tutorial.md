---
layout: tutorial_hands_on

title: A Docker-based interactive Jupyterlab powered by GPU for artificial intelligence in Galaxy
zenodo_link: https://zenodo.org/record/6091361
questions:
- How to use Jupyterlab and it several features?
- How to use it for creating input datasets and writing artificial intelligence (AI) algorithms?
objectives:
- Learn to use Jupyterlab - an online Python editor designed for developing AI algorithms
- Explore several of its features such as Git integration, workflow of jupyter notebook, integration to Galaxy
- Develop AI algorithms using scikit-learn and tensorflow
- Send long-running jobs to Galaxy's cluster and save results in its history
- Reproduce results from recent scientific papers - COVID CT scan segmentation and 3D protein structure prediction
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
tags:
- interactive-tools
- machine-learning
- deep-learning
- jupyter-lab
- image-segmentation
contributors:
- anuprulez

---

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Introduction
[Jupyterlab](https://jupyterlab.readthedocs.io/en/stable/) is a popular integrated development environment (IDE) for a variety of tasks in data science such as prototyping analyses, creating meaningful plots, data manipulation and preprocessing. Python is one of the most used languages in such an environment. Given the usefulness of Jupyterlab, more importantly in online platforms, a robust Jupyterlab notebook application has been developed that is powered by GPU acceleration and contains numerous packages such as Pandas, Numpy, Scipy, Scikit-learn, Tensorflow, ONNX for modern data science. It has been developed as an interactive Galaxy tool that runs on an isolated [docker container](https://github.com/anuprulez/ml-jupyter-notebook). The docker container has been built using [jupyter/tensorflow-notebook:tensorflow-2.6.0](https://hub.docker.com/layers/jupyter/tensorflow-notebook/tensorflow-2.6.0/images/sha256-23ac04200f8f6ee4a224f7f8cb48c810e24bc18cdd20a3614d1914e01f950969?context=explore) as the base container. Moreover, a Galaxy [tool](https://github.com/bgruening/galaxytools/pull/1157) ( `run_jupyter_job`) can be executed using [Bioblend](https://bioblend.readthedocs.io/) which makes use of the Galaxy's remote job handling for long-running deep learning training. The long-running training happens on a Galaxy cluster and the outcome datasets such as trained models, tabular files and so on are saved in Galaxy history.


## Jupyterlab features
[Jupyterlab](https://jupyterlab.readthedocs.io/en/stable/) notebook has been augmented with several useful features and together they make it ready-to-use for quick prototyping of AI projects. Being able to **serve online** makes it convenient to share it with other researchers and users. GPUs have accelerated AI research, especially deep learning. Therefore, the backend of the Jupyterlab is **powered by GPU** to make long-running AI training programs finish faster by parallelizing matrix multiplications. Galaxy tool for remote job processing can also run on GPU. It is also integrated with **Git version control** making it easy to clone, commit and maintain codebase directly in the notebook. For example, repositories from Github can be easily cloned, updated and maintained. In addition, a standard model format, **open neural network exchange (ONNX)**, has been supported in the notebook to transform scikit-learn and tensorflow models to "onnx" files. These model files can be conveniently shared and used for inference. Galaxy also supports "onnx" file format which makes it easier to interoperate such models from a notebook to Galaxy and vice-versa. Using Bioblend, the notebook can be **connected to Galaxy** and its histories, tools and workflows can be accessed directly in the notebook. They can easily be operated via commands from the notebook. Many notebooks are created serving different purposes. These notebooks can be knit together to form one pipeline where each notebook transforms data taking a different form of data from its previous notebook and passing on the transformed data to its next notebook. An example workflow created using **[Elyra AI](https://elyra.readthedocs.io/en/stable/)** can be found by name "METABRIC\_ML.pipeline" in the "elyra" folder of the opened Jupyterlab instance. These pipelines can also be executed remotely on different clusters too using its "runtimes" features that pull different docker containers. An example pipeline created in the notebook can be seen in the picture below.

![Elyra AI pipeline](../../images/elyra_ai_pipeline.png "An example ML pipeline created using Elyra AI in a notebook. Each connected block of the pipeline is a different notebook")

There are several additional packages suited for performing machine learning tasks such as [Open-CV](https://github.com/opencv/opencv) and [Scikit-Image](https://scikit-image.org/) for image processing, [NiBabel](https://nipy.org/nibabel/) package for processing images files. These packages are available to researchers and users in this notebook. Many other features such as GPU utilization dashboards for monitoring the GPU usage and system memory utilization, Voila for rendering output cells of a notebook in a separate tab hiding all code cells, interactive bqplots are integrated into the Jupyterlab notebook for researchers and users to manage their projects more efficiently.

In the following sections, we will look at an example machine learning task, image segmentation, and go through each step of how to perform such analysis in a Jupyterlab notebook. In addition, we will also look at how to predict 3D structure of proteins using the same infrastructure.

## Image segmentation
In an image classification task, a label is assigned to each image. For example, an image containing a cat gets a label as "cat" ("cat" is given an integer to be used in a machine learning program). But, in an image segmentation task, each pixel of an image gets a label. For example, in an image containing a cat, all the pixels occupied by the cat in the image are denoted by violet color [see Figure 2](https://www.tensorflow.org/tutorials/images/segmentation). Image segmentation is widely used in many fields such as medical imaging to find out abnormal regions in images of tissues, self-driving cars to detect different objects on highways and many more.

![Image segmentation for a cat](../../images/cat_segmentation.png "Image segmentation for cat. The violet region of the middle plot denotes the cat and the outline of the violet section denoted by the yellow color is the cat's mask")


## COVID CT scan
Segmenting lung CT scans task to locate infected regions has been actively proposed to augment the RT-PCR testing for initial screening of COVID infection in humans. Deep learning has been used to predict these regions with high accuracy in studies such as {% cite InfNet %}, {% cite JCSWu2021 %} and many more. In Figure 3, the differences between the CT scans of a normal person and a person suffering from COVID-19 can be seen. The regions of lungs marked by white patches are infected by {% cite SAEEDIZADEH2021100007 %}.

![Normal and Covid CT scans](../../images/normal_covid_ct_scans.jpg "(Left) Lungs CT scan for a normal person and (right) lungs CT scan for a person having COVID-19")


## CT scans and masks
Figure 4 shows CT scans of infected lungs (top row) and borders around the infected regions have been drawn with red color (middle row). In the last row respective masks, from the CT scans, have been taken out denoting the infected regions. These masks are the "segmented" regions from the corresponding CT scans. While creating the dataset for training the deep learning model (Unet), the CT scans are the data points and its "known" respective masks of infected regions become their labels ({% cite SAEEDIZADEH2021100007 %}). 

![Masks](../../images/covid_ct_scan_masks.png "Picture shows CT scans and their masks of infected regions that are drawn with red color in middle row plots")

### Unet neural network
Unet neural network is widely used for segmentation tasks in images ({% cite Ronnebergerunet %}). The name "Unet" resembles the shape of its architecture as "U" (see Figure 5). It has two parts - encoder and decoder. The left half of the "U" shape is the encoder that learns the feature representation (downsampling) at a lower level. The right half is the decoder that maps the low-resolution feature representation on the higher resolution (upsampling) pixel space.

![Unet architecture](../../images/Unet.jpg "Architecture of Unet neural network for image segmentation. The name "Unet" resembles the shape of its architecture as 'U'")

In this tutorial, we will use the dataset of CTs scans and their respective masks to train an Unet neural network model. The model learns to map the infected regions in the CT scans to their masks. For prediction, the trained model is given unseen CT scans and for each CT scan, it predicts infected regions or masks. For this experiment, we will use Jupyterlab for cloning the notebooks from Github that contain scripts for downloading data, creating and training Unet model and prediction tasks. Data (CT scans and also the trained model) required for the notebooks can be downloaded from [Zenodo](https://zenodo.org/record/6091361). The model can either be trained in the Jupyterlab or can be sent to Galaxy's cluster for remote processing. After remote processing, the created datasets such as the trained model become available in the new Galaxy history.

## Protein 3D structure
Understanding the structure of proteins provides insights into their functions. But, only a few proteins have known structures out of billions of known proteins. To advance into the direction of predicting their 3D structures only from their sequences, AlphaFold2 ({% cite AlphaFold2 %}) has made a breakthrough to predict their structures with high accuracy. However, the database that it uses are large, approximately 2 TBs, and are hard to manage. Therefore, to make it more accessible, a few approaches such as ColabFold ({% cite Mirdita2021 %}), has been developed that consumes significantly less memory, much faster and produces 3D structures with similar accuracy. We will look at, in later sections, how we can use ColabFold to predict 3D structure of a protein sequence using the Jupyterlab infrastructure in Galaxy.

# Reproduce results from recent publications using Jupyterlab infrastructure in Galaxy
In this tutorial, we will discuss a few features of Jupyterlab to create and train a Unet deep learning model and then, predict segmented regions in COVID CT scans using the Unet trained model. In addition, we will predict 3D structure of protein sequence using ColabFold, a faster implementation of AlphaFold2.

## Open Jupyterlab editor

> ### {% icon hands_on %} Hands-on: GPU enabled Interactive Jupyter Notebook for Machine Learning
>
> - {% tool [GPU enabled Interactive Jupyter Notebook for Machine Learning](interactive_tool_ml_jupyter_notebook) %}
>    - *"Do you already have a notebook?"*: `Start with a fresh notebook`
>    - Click *"Execute"*
>
{: .hands_on}

Now, we should wait for a few minutes until Galaxy creates the required compute environment for opening a new Jupyterlab. Usually, this task takes around 10-15 minutes. The progress can be checked by clicking on the "User>Active Interactive tools". On "Active Interactive Tools" page, there is a list of all open interactive tools. The interactive tool can be found by "GPU enabled Interactive Jupyter Notebook for Machine Learning" name. When the job info starts showing "running" on the "Active Interactive Tools" page, then the name of the interactive tool gets associated with a URL to the running Jupyterlab. On clicking this URL, the running Jupyterlab can be opened. Several features of Jupyterlab editor running in Galaxy can be learned by going through the "home_page.ipynb" notebook. The folder "notebooks" contain several notebooks that show multiple use-cases of features and packages. In the following sections, we will look at two use-cases of reproducing results from recent papers - one is COVID CT scan segmentation ({% cite SAEEDIZADEH2021100007 %}) using deep learning and another one is prediction 3D structure of protein sequences using ColabFold ({% cite Mirdita2021 %}).


## Use-case 1: Run image segmentation analysis

### Clone Github repository
To use Git version control for cloning any codebase from GitHub, the following steps should be performed.

> ### {% icon hands_on %} Hands-on: Pull code
>
> 1. Create a new folder named "covid_ct_segmentation" alongside other folders such as "data", "outputs", "elyra"
> 2. Inside the created folder, clone a code repository by clicking on "Git" icon as shown in Figure 6
> 3. In the shown popup, provide the repository path as "https://github.com/anuprulez/gpu_jupyterlab_ct_image_segmentation" and then, click on "clone"
> 4. The repository "anuprulez/gpu_jupyterlab_ct_image_segmentation" gets immediately cloned
> 5. Move inside the created folder `gpu_jupyterlab_ct_image_segmentation`. A few notebooks can be found inside that are numbered
>    ![Clone repository](../../images/git_clone.png "Clone a code repository using Git")
>
{: .hands_on}

Now, we have all the notebooks available for performing image segmentation. The entire analysis can be completed in two ways - one by running all the notebooks in the Jupyterlab itself and another by sending the notebook for training to a remote Galaxy cluster invoking another Galaxy tool from the notebook. First, let's look at how we can run this analysis in the notebook itself.


### Run notebooks in Jupyterlab

> ### {% icon hands_on %} Hands-on: Run notebooks in Jupyterlab
>
> 1. Download and save datasets in the notebook using **"1_fetch_datasets.ipynb"** notebook. 
>
> It will also create all the necessary folders. It downloads two datasets, one "h5" file containing many matrices as sub-datasets belonging to training data, training labels, validation data, validation labels, test data and test labels. These sub-datasets are stored in different variables after reading the original "h5" file once. For training, we only need these datasets/matrices - training data, training labels, validation data and validation labels. The matrices, test data and test labels, are used for prediction. We use "h5" format for storing and retrieving datasets as all AI algorithms need input datasets in the form of matrices. Since, in the field of AI, there are many different types of datasets such as images, sequences, real numbers and so on, therefore, to converge all these different forms of datasets to one format, we use "h5" to store matrices. In any AI analysis, different forms of datasets can be stored as `h5` files. 

> For image segmentation task, we also saved all the input datasets/matrices to the deep learning model as sub-datasets in one "h5" file so that it can be easily created, stored, downloaded and used. This step may take a few minutes as it downloads around 450 MB of data from Zenodo. Once datasets are downloaded to the notebook, we can move to the next step of creating a deep learning model and start training it.
> 
> 2. Create and train a deep learning Unet model using **"2_create_model_and_train.ipynb"** notebook. 
> 
>This will read the input "h5" datasets containing all images and train the model after creating deep learning model architecture. This notebook first creates a deep learning architecture based on Unet including custom loss functions such as total variation and binary cross-entropy losses. After creating the deep learning architecture, all training datasets such as training data, training labels, validation data, validation labels are loaded from the combined "h5" file. In the next step, all the datasets and deep learning architecture are compiled together and training starts for 10 epochs (10 iterations over the entire training dataset).
> 
> The training is fast as it runs on GPU and finishes in a few minutes and creates a trained model. In the last step, the trained model containing several files are converted to one "onnx" file. Once a trained model is ready, we can move to the next step to make predictions on unseen CT scan masks.
> 
> 3. Predict unseen masks using the trained model in **"3_predict_masks.ipynb"** notebook. First, it reads test datasets from the combined "h5" file and then, loads the "onnx" model. Using this model, it predicts masks of unseen CT scans and then plots the ground truth and predicted masks in one plot (see Figure 7).
> 
> ![True and predicted masks of CTs scans](../../images/true_pred_masks.png "Ground truth and predicted masks of COVID CT scans. Predicted masks in the third and fourth columns (computed using slightly different loss functions) are very close to the corresponding original masks in the second column")
> The above picture shows origin lungs CT scans in the first column. The second column shows the corresponding true infected regions in white. The third and fourth columns show the infected regions predicted with different loss functions, one is the binary cross-entropy loss and another is the combination of binary cross-entropy loss and total variation loss. Binary cross-entropy calculates loss between two corresponding pixels of true and predicted images. It measures how close two corresponding pixels are. But, total variation loss measures the amount of noise in predicted images as noisy regions usually show large variations in their neighbourhood. The noise in the predicted images gets minimized and the connectivity of predicted masks improves as well when minimizing this loss.   
>
{: .hands_on}


### Run Jupyterlab notebooks remotely on a Galaxy cluster

The training task completed in the notebook above can also be sent to a Galaxy cluster by executing a Galaxy tool in the notebook itself. Using [Bioblend](https://bioblend.readthedocs.io/) APIs, datasets in the form of `h5` files and notebook are uploaded to a Galaxy history and then a Galaxy tool `run_jupyter_job` executes the notebook using the uploaded dataset on a Galaxy cluster and creates a trained model in the Galaxy history. Let's look at how to execute a notebook remotely on a Galaxy cluster.

> ### {% icon hands_on %} Hands-on: Run Jupyterlab notebooks remotely on a Galaxy cluster
>
> 1. Download and save datasets in the notebook using **"1_fetch_datasets.ipynb"** notebook in the same way as before. Ignore this step if this notebook has already been executed.
> 2. Execute **"5_run_remote_training.ipynb"** notebook to dynamically execute code inside **"4_create_model_and_train_remote.ipynb"** notebook on a remote cluster using a different Galaxy tool (`run_jupyter_job`). The notebook **"5_run_remote_training.ipynb"** provides the path of a notebook to be executed remotely along with the datasets to the Galaxy tool by calling a custom function `run_script_job` which is part of the Jupyterlab notebook. Executing **"5_run_remote_training.ipynb"** uploads datasets and dynamic Python script, extracted from the **"4_create_model_and_train_remote.ipynb"** notebook, to a newly created Galaxy history (Figure 8).
> 
> When the task of uploading dataset and dynamic code is finished, the Galaxy tool (`run_jupyter_job`) executes the dynamically uploaded script with the uploaded dataset on a remote Galaxy cluster which is similar to running any other Galaxy tool. When the Galaxy tool (`run_jupyter_job`) finishes its execution, the resulting models and other datasets appear in the created Galaxy history. While the job is running on the Galaxy cluster, the Jupyter notebook can be closed as the model training task gets decoupled from the notebook and is entirely transferred to the Galaxy cluster. Specific history running this job can be accessed in Galaxy. The name of the Galaxy history created by (`run_jupyter_job`) method can be changed by passing a new history name (parameter `new_history_name`).
> 
> ![Galaxy history](../../images/finished_history_remote_ai.png "Galaxy history showing finished datasets after remote training on a Galaxy cluster")
> 
> **Note**: The training may take longer depending on how busy Galaxy's queueing is as it sends the training task to be done on a Galaxy cluster. Therefore, this feature should be used when the training task is expected to run for at least a couple of hours. The training time is higher because a large Docker container is downloaded on the assigned cluster and only then, the training task can proceed.
>
{: .hands_on}


When the Galaxy job finishes, it creates a history having few datasets such as input `h5` dataset, dynamic Python script extracted from the notebook, a trained model saved as `onnx` file, all saved lists and NumPy arrays and a zipped file containing all the intermediate files created and saved in the dynamic script. One such example of a Galaxy history is available [here](https://usegalaxy.eu/u/kumara/h/ctsegmentationmarch18). Datasets from this history can be downloaded into the Jupyterlab notebook and further analyses can be carried on. Using the notebook **"6_predict_masks_remote_model.ipynb"**, the trained model can be downloaded from Galaxy history to Jupyterlab notebook and then, it can be used to make predictions of unseen CT scan masks in the same way as shown in **"3_predict_masks.ipynb"** notebook.


## Use-case 2: Run ColabFold to predict 3D structure of protein

Google Deepmind’s AlphaFold2 ({% cite AlphaFold2 %}) has made a breakthrough in predicting the 3D structure of proteins with outstanding accuracy. However, due to their large database size (a few TBs), it is not easily accessible to researchers. Therefore, a few approaches have been developed that replace the time-consuming steps of AlphaFold2 with slightly different steps to create input features and predict the 3D structure with similar accuracy as AlphaFold2. One such approach is ColabFold ({% cite Mirdita2021 %}) which replaces the large database search in AlphaFold2 for finding homologous sequences by a significantly faster [MMseqs2 API](https://github.com/soedinglab/MMseqs2) call to generate input features based on the query protein sequence. ColabFold is approximately 16 times faster than AlphaFold2 in predicting 3D structures of protein sequences. It has been integrated into the Docker container by adding two packages - ColabFold and the GPU enabled [JAX](https://github.com/google/jax) which is just-in-time compiler from Google used for mathematical transformations. In  the following hands-on section, we will run a notebook and predict a 3D structure of a protein sequence. 

> ### {% icon hands_on %} Hands-on: Predict 3D structure of protein using ColabFold
>
> 1. Clone a Github repository by following steps described in the previous hands-on section `Clone Github repository`. Ignore if already done.
> 2. Open **7_ColabFold_MMseq2.ipynb** notebook and execute it.
> 3. It needs a few minutes (around 5-10 minutes) to download AlphaFold2 weigths (of size approximately 4 GB) and uses these weights for prediction. An example prediction is shown in figure 9. 
> !["3D structure of a protein sequence of length 300 with side-chains predicted by ColabFold"](../../images/3D_300_L_protein.png "3D structure of a protein sequence of length 300 with side-chains predicted by ColabFold")
>
{: .hands_on}


## Docker's security
For remote training, a dynamic Python script is sent to a Galaxy tool for execution on a Galaxy cluster that may pose security risks of containing malicious code. To minimize the security risks that comes along with executing dynamic Python scripts, the Galaxy tool `run_jupyter_job` is executed inside a [Docker container](https://github.com/anuprulez/ml-jupyter-notebook). It is the same container inside which the Jupyterlab notebook runs. All the packages available inside the notebook are also available for the Galaxy tool. Programs run in isolated environments inside Docker containers and have their connections limited to necessary libraries, and network connections. Due to these restrictions, they have a reduced number of interactions with the host operating system.


# Conclusion
In this tutorial, we have discussed several features of the Jupyterlab notebook in Galaxy and used a few of them to perform an image segmentation task on COVID CT scans by training Unet, a deep learning model. In addition, we have used it for predicting 3D structure of a protein sequences which is a memory-intensive and time-consuming task. The ready-to-use infrastructure that we used here for developing an AI algorithm and training a model proves convenient for making quick prototypes, accelerated by GPU, in the data science field without caring about the storage and compute resources. 

The possibility of remote execution of long-running training tasks obviates keeping the notebook or computer open till the training finishes and at the same time, the datasets and models stay safe in Galaxy history. The usage of `h5` files for input datasets to a neural network solves the problem of varying data formats across multiple fields that use AI algorithms extensively for predictive tasks. Using `onnx` for creating, storing and sharing AI models makes it easier to handle such a model as one compact file. Git integration makes it easy to maintain an entire repository directly in the notebook. The entire environment can be shared with other researchers only via sharing the URL. Bioblend connects this infrastructure to Galaxy which makes its histories, workflows and tools accessible directly via the notebook. Overall, this infrastructure would prove useful for researchers working to create quick prototypes of their applications.
