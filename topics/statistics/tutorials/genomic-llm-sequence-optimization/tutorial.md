---
layout: tutorial_hands_on
title: Optimizing DNA Sequences for Biological Functions using a DNA LLM
level: Intermediate
draft: true
requirements:
-
  type: "internal"
  topic_name: data-science
  tutorials:
  - python-basics
  - python-warmup-stat-ml
-
  type: "internal"
  topic_name: statistics
  tutorials:
  - intro-to-ml-with-python
  - neural-networks-with-python
  - deep-learning-without-gai-with-python
  - genomic-llm-pretraining
  - genomic-llm-finetuning
  - genomic-llm-zeroshot-prediction
  - genomic-llm-sequence-generation
questions:
- to do
objectives:
- pretraining LLM for DNA
- finetuning LLM
- zeroshot prediction for DNA variants and synthetic DNA sequence generation.
time_estimation: 3H
key_points:
- To be added
contributions:
  authorship:
  - raphaelmourad
  - bebatut
tags:
- elixir
- ai-ml
- Large Language Model
subtopic: gai-llm
priority: 5
notebook:
  language: python
  pyolite: true
---


# Prepare resources 

## Install dependencies


```python
!pip install datasets==3.0.1
!pip install torch==2.5.0
!pip install transformers -U
!pip install accelerate==1.1.0
!pip install peft==0.13.2
!pip install bitsandbytes==0.44.1
!pip install flash-attn==2.6.3
!pip install Bio==1.7.1
!pip install orfipy
```

## Import Python libraries

- `torch`
- `flash_attn`
- `numpy`
- `transformers`
  - `AutoTokenizer`
  - `EarlyStoppingCallback`
  - `Trainer`
  - `TrainingArguments`
  - `AutoModelForCausalLM`
  - `AutoConfig`
  - `DataCollatorForLanguageModeling`
- `datasets`
- `accelerate`


```python
# 
# During the class, 
import os

import sys
import time
from os import path
import gc


import flash_attn
import torch
import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt

import transformers

from transformers import AutoTokenizer
from transformers import EarlyStoppingCallback, Trainer, TrainingArguments
from transformers import AutoModelForCausalLM, AutoConfig
from transformers import DataCollatorForLanguageModeling
from datasets import load_dataset

import accelerate
```

## Check versions

Numpy version > 1.26.4

```python
np.__version__
```

transformers version > 4.47.1

```python
transformers.__version__
```

flash_attn > 2.6.0.post1 and 2.7.0.post2

```python
flash_attn.__version__
```

accelerate > 0.32.1

```python
# Tested with accelerate==0.32.1

accelerate.__version__
```

## Prepare GPU


```python
# CHECK GPU
# We can see how many VRAM is used and how much the GPU is used.
!nvidia-smi
```

    Thu Feb  6 07:59:23 2025       
    +-----------------------------------------------------------------------------------------+
    | NVIDIA-SMI 550.54.15              Driver Version: 550.54.15      CUDA Version: 12.4     |
    |-----------------------------------------+------------------------+----------------------+
    | GPU  Name                 Persistence-M | Bus-Id          Disp.A | Volatile Uncorr. ECC |
    | Fan  Temp   Perf          Pwr:Usage/Cap |           Memory-Usage | GPU-Util  Compute M. |
    |                                         |                        |               MIG M. |
    |=========================================+========================+======================|
    |   0  Tesla T4                       Off |   00000000:00:04.0 Off |                    0 |
    | N/A   36C    P8              9W /   70W |       0MiB /  15360MiB |      0%      Default |
    |                                         |                        |                  N/A |
    +-----------------------------------------+------------------------+----------------------+
                                                                                             
    +-----------------------------------------------------------------------------------------+
    | Processes:                                                                              |
    |  GPU   GI   CI        PID   Type   Process name                              GPU Memory |
    |        ID   ID                                                               Usage      |
    |=========================================================================================|
    |  No running processes found                                                             |
    +-----------------------------------------------------------------------------------------+


```python
# LOOK AT GPU USAGE AND RAM
!nvidia-smi
```

    Thu Feb  6 16:49:41 2025       
    +-----------------------------------------------------------------------------------------+
    | NVIDIA-SMI 550.54.15              Driver Version: 550.54.15      CUDA Version: 12.4     |
    |-----------------------------------------+------------------------+----------------------+
    | GPU  Name                 Persistence-M | Bus-Id          Disp.A | Volatile Uncorr. ECC |
    | Fan  Temp   Perf          Pwr:Usage/Cap |           Memory-Usage | GPU-Util  Compute M. |
    |                                         |                        |               MIG M. |
    |=========================================+========================+======================|
    |   0  Tesla T4                       Off |   00000000:00:04.0 Off |                    0 |
    | N/A   54C    P8             10W /   70W |       2MiB /  15360MiB |      0%      Default |
    |                                         |                        |                  N/A |
    +-----------------------------------------+------------------------+----------------------+
                                                                                             
    +-----------------------------------------------------------------------------------------+
    | Processes:                                                                              |
    |  GPU   GI   CI        PID   Type   Process name                              GPU Memory |
    |        ID   ID                                                               Usage      |
    |=========================================================================================|
    |  No running processes found                                                             |
    +-----------------------------------------------------------------------------------------+



Control the use of ram by CUDA

```python
torch.backends.cudnn.benchmark=True
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:32 "
```

```python
# Check GPU
import torch
torch.device('cuda' if torch.cuda.is_available() else 'cpu')
```

    device(type='cuda')

# Get the model 

Mistral-DNA-v0.1 was derived from Mixtral-8x7B for the human genome. Mixtral-8x7B was simplified for DNA: the number of layers and the hidden size were reduced. The model was pretrained using the human genome hg38 with 200B DNA sequences.

The model can be downloaded on HuggingFace: https://huggingface.co/RaphaelMourad/Mistral-DNA-v0.1

```python
!git clone https://github.com/raphaelmourad/Mistral-DNA.git
```


```python
!ls
```

    Mistral-DNA  sample_data


```python
# SET DIRECTORY
os.chdir("Mistral-DNA/")
print(os.getcwd())
```

/content/Mistral-DNA

