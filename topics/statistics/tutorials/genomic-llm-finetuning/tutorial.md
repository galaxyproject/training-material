---
layout: tutorial_hands_on
title: Fine-tuning a LLM for DNA Sequence Classification
level: Intermediate
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
questions:
- How to classify a DNA sequence depending on if it binds a protein or not (transcription factor)?
objectives:
- Load a pre-trained model and modify its architecture to include a classification layer.
- Prepare and preprocess labeled DNA sequences for fine-tuning.
- Define and configure training parameters to optimize the model's performance on the classification task.
- Evaluate the fine-tuned model's accuracy and robustness in distinguishing between different classes of DNA sequences.
time_estimation: 3H
key_points:
- Fine-tuning pre-trained LLMs reduces training time and computational needs, making advanced research accessible.
- Techniques like LoRA enable fine-tuning on modest hardware, broadening access to powerful models.
- Rigorous testing on unseen data confirms a model's practical applicability and reliability.
contributions:
  authorship:
  - raphaelmourad
  - bebatut
tags:
- elixir
- ai-ml
- Large Language Model
subtopic: gai-llm
priority: 2
notebook:
  language: python
  pyolite: true
---

After preparing, training, and utilizing a language model for DNA sequences, we can now fine-tune a pre-trained Large Language Model (LLM) for specific DNA sequence classification tasks. Here, we will use a pre-trained model from Hugging Face, specifically the [Mistral-DNA-v1-17M-hg38](https://huggingface.co/RaphaelMourad/Mistral-DNA-v1-17M-hg38), and adapt it to classify DNA sequences based on their biological functions. Our objective is to classify sequences according to whether they bind to transcription factors.

> <comment-title>Transcription factors</comment-title>
>
> Transcription factors are proteins that play a crucial role in regulating gene expression by binding to specific DNA sequences, known as enhancers or promoters. These proteins act as molecular switches, turning genes on or off in response to various cellular signals and environmental cues. By binding to DNA, transcription factors either promote or inhibit the recruitment of RNA polymerase, the enzyme responsible for transcribing DNA into RNA, thereby influencing the rate of transcription.
> 
> ![Diagram illustrating DNA binding with CTCF. The left panel, outlined in red, shows a DNA sequence 'CCACCAGGGGGCGC' labeled as 'DNA binding CTCF,' with an oval labeled 'CTCF' above it. The right panel, outlined in blue, shows a different DNA sequence 'GTGGCTAGTAGGTAG' labeled as 'DNA not binding CTCF,' indicating that this sequence does not interact with CTCF.](images/two_dna_sequences.png "Two types of DNA sequences. On the left, a DNA sequence that binds the transcription factor CTCF. On the right, a DNA sequence that does not bind CTCF.")
>
> Transcription factors are essential for numerous biological processes, including cell differentiation, development, and response to external stimuli. Their ability to recognize and bind specific DNA sequences allows them to orchestrate complex gene expression programs, ensuring that the right genes are expressed at the right time and in the right place within an organism. Understanding the function and regulation of transcription factors is vital for deciphering the molecular mechanisms underlying health and disease, and it opens avenues for developing targeted therapeutic interventions.
>
{: .comment}

This classification task is crucial for understanding gene regulation, as transcription factors play a vital role in controlling which genes are expressed in a cell. By training a model to predict whether a DNA sequence binds to a transcription factor, we can gain insights into regulatory mechanisms and potentially identify novel binding sites or understand the impact of genetic variations on transcription factor binding.

By fine-tuning the model, we aim to leverage its pre-trained knowledge of DNA sequences to achieve high accuracy in this classification task. This tutorial will guide you through the necessary steps, from data preparation to model evaluation, ensuring you can apply these techniques to your own research or projects.

We will use [`Mistral-DNA-v1-17M-hg38`](https://huggingface.co/RaphaelMourad/Mistral-DNA-v1-1M-hg38), a mixed model that was pre-trained on the entire Human Genome. It contains approximately 17 million parameters and was trained using the Human Genome assembly GRCh38 on sequences of 10,000 bases (10K):

```python
model_name="RaphaelMourad/Mistral-DNA-v1-17M-hg38"
```

> <comment-title>Pretraining a LLM</comment-title>
>
> To learn how to pretrain a LLM on DNA, please follow the dedicated ["Pretraining a Large Language Model (LLM) from Scratch on DNA Sequences"]({% link topics/statistics/tutorials/genomic-llm-pretraining/tutorial.md %}) tutorial
>
{: .comment}


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Prepare resources

## Install dependencies

The first step is to install the required dependencies:

```python
!pip install accelerate==1.1.0
!pip install peft==0.13.2
!pip install torch==2.5.0
!pip install transformers -U
!pip install progressbar
!pip install bitsandbytes
```

> <question-title></question-title>
>
> 1. What is `accelerate`?
> 2. What is `peft`?
> 3. What is `torch`?
> 4. What is `transformers`?
>
> > <solution-title></solution-title>
> >
> > 1. `accelerate` is a library by [Hugging Face](https://huggingface.co/) -- a platform that provides tools and resources for building, training, and deploying machine learning models -- designed to simplify the process of training and deploying machine learning models across different hardware environments. It provides tools to optimize performance on GPUs, TPUs, and other accelerators, making it easier to scale models efficiently.
> >
> > 2. The PEFT (Parameter-Efficient Fine-Tuning) Python library, developed by Hugging Face, is a tool designed to efficiently adapt large pretrained models to various downstream tasks without the need to fine-tune all of the model's parameters. By focusing on a small subset of parameters, PEFT significantly reduces computational and storage costs, making it feasible to fine-tune large language models (LLMs) on consumer-grade hardware. The library integrates seamlessly with the Hugging Face ecosystem, including Transformers, Diffusers, and Accelerate, enabling streamlined model training and inference. PEFT supports techniques like LoRA (Low-Rank Adaptation) and prompt tuning, and it can be combined with quantization to further optimize resource usage. Its open-source nature fosters collaboration and accessibility, allowing developers to customize models for specific applications quickly and efficiently.
> > 
> > 3. `torch`, also known as PyTorch, it is an open-source machine learning library developed by Facebook's AI Research lab. It provides a flexible platform for building and training neural networks, with a focus on tensor computations and automatic differentiation.
> > 
> > 4. `transformers` is a library by Hugging Face that provides implementations of state-of-the-art transformer models for natural language processing (NLP). It includes pre-trained models and tools for fine-tuning, making it easier to apply transformers to various NLP tasks.
> >
> {: .solution}
>
{: .question}

## Import Python libraries

Let's now import them.

```python
import os

import accelerate
import flash_attn
import numpy as np
import pandas as pd
import torch
import transformers
from accelerate import FullyShardedDataParallelPlugin, Accelerator
from pathlib import Path
from peft import (
    LoraConfig,
    get_peft_model,
    get_peft_model_state_dict,
    prepare_model_for_kbit_training,
)
from progressbar import ProgressBar
from random import randrange
from torch.utils.data import TensorDataset, DataLoader
from torch.distributed.fsdp.fully_sharded_data_parallel import FullOptimStateDictConfig, FullStateDictConfig
from transformers import (
    AutoTokenizer,
    AutoModel,
    BitsAndBytesConfig,
    EarlyStoppingCallback,
    set_seed,
)
```

> <comment-title>Versions</comment-title>
>
> This tutorial has been tested with following versions:
> - `numpy` = 1.19 (and not 1.2)
> - `transformers` > 4.47.1
> 
> You can check the versions with:
>
> ```python
> np.__version__
> transformers.__version__
> ```
{: .comment}

# Configure fine-tuning

## Check and configure available resources

We select the appropriate device (CUDA-enabled GPU if available) for running PyTorch operations

```python
torch.device('cuda' if torch.cuda.is_available() else 'cpu')
```

Let's check the GPU usage and RAM:

```python
!nvidia-smi
```

We now set an environment variable that configures how PyTorch manages CUDA memory allocations

```python
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:32"
```

## Specify settings for quantization

Quantization is a technique used in machine learning and signal processing to reduce the precision of numerical values, typically to decrease memory usage and computational requirements. This process is particularly useful when working with large models as it allows them to be deployed on hardware with limited resources without significantly sacrificing performance.

Here, we use `BitsAndBytesConfig` to configure a 4-bit quantization. Using 4-bit precision reduces the memory footprint of the model, which is particularly useful for very large models that might not fit into GPU memory otherwise:

```python
bnb_config = BitsAndBytesConfig(
    load_in_4bit=True,
    bnb_4bit_use_double_quant=True,
    bnb_4bit_compute_dtype=torch.bfloat16
)
```

> <question-title></question-title>
>
> What do the parameters?
>
> 1. `load_in_4bit=True`
> 2. `bnb_4bit_use_double_quant=True`
> 3. `bnb_4bit_compute_dtype=torch.bfloat16`
>
> > <solution-title></solution-title>
> >
> > 1. `load_in_4bit=True`: Specifies that the model should be loaded with 4-bit quantization. Using 4-bit precision reduces the memory footprint of the model, which is particularly useful for very large models that might not fit into GPU memory otherwise.
> >
> > 2. `bnb_4bit_use_double_quant=True`: enables double quantization, which means that the quantization constants from the first quantization are quantized again. This further reduces the memory footprint, although it may introduce additional computational overhead.
> >
> > 3. `bnb_4bit_compute_dtype=torch.bfloat16`: sets the compute data type to bfloat16 (Brain Floating Point 16-bit format). Using bfloat16 can provide a good balance between computational efficiency and numerical stability, especially on hardware that supports this format, such as certain GPUs and TPUs.
> >
> {: .solution}
>
{: .question}

## Configure Accelerate

Now, we will configure the [Hugging Face Accelerate library](https://huggingface.co/docs/accelerate/en/index) to optimize the training process for large models using Fully Sharded Data Parallel (FSDP). This setup is crucial for efficiently utilizing GPU resources and enabling distributed training across multiple devices.

First, we need to configure the FSDP plugin, which will manage how model parameters and optimizer states are sharded across GPUs. This configuration helps in reducing memory usage and allows for the training of larger models.

```python
fsdp_plugin = FullyShardedDataParallelPlugin(
    state_dict_config=FullStateDictConfig(offload_to_cpu=True, rank0_only=False),
    optim_state_dict_config=FullOptimStateDictConfig(offload_to_cpu=True, rank0_only=False),
)
```
> <question-title></question-title>
>
> What do the parameters?
> 
> 1. `state_dict_config=FullStateDictConfig(offload_to_cpu=True, rank0_only=False)`?
> 2. `optim_state_dict_config=FullOptimStateDictConfig(offload_to_cpu=True, rank0_only=False)`?
>
> > <solution-title></solution-title>
> >
> > 1. `state_dict_config=FullStateDictConfig(offload_to_cpu=True, rank0_only=False)`
> >    - `FullStateDictConfig`: Configures how the model's state dictionary (parameters) is managed.
> >    - `offload_to_cpu=True`: Specifies that the model's parameters should be offloaded to CPU memory when not in use. This helps free up GPU memory, especially useful when working with large models.
> >    - `rank0_only=False`: Indicates that the state dictionary operations (like saving and loading) are not restricted to the rank 0 process. This allows all processes to participate in these operations, which can be beneficial for distributed training setups.
> >
> > 2. `optim_state_dict_config=FullOptimStateDictConfig(offload_to_cpu=True, rank0_only=False)`
> >    - `FullOptimStateDictConfig`: Configures how the optimizer's state dictionary is managed.
> >    - `offload_to_cpu=True`: Similar to the model's state dictionary, this setting offloads the optimizer states to CPU memory when not in use, further reducing GPU memory usage.
> >    - `rank0_only=False`: Allows all processes to handle the optimizer state dictionary operations, ensuring that the optimizer states are managed efficiently across the distributed setup.
> >
> {: .solution}
>
{: .question}

Next, we initialize the Accelerator from the Hugging Face Accelerate library, integrating the FSDP plugin for seamless distributed training:

```python
accelerator = Accelerator(fsdp_plugin=fsdp_plugin)
```

By passing the FSDP plugin to the `Accelerator`, we enable sharded data parallelism, which efficiently manages model and optimizer states across multiple GPUs.

With this configuration, the `Accelerator` will handle the complexities of distributed training, allowing us to focus on developing and experimenting with our models. This setup is particularly beneficial when working with large-scale models and limited GPU resources, as it optimizes memory usage and enables faster training times.

## Configure LoRA for Parameter-Efficient Fine-Tuning

We will configure the LoRA (Low-Rank Adaptation) settings for parameter-efficient fine-tuning of a large language model. LoRA is a technique that allows us to fine-tune only a small number of additional parameters while keeping the original model weights frozen, making it highly efficient for adapting large models to specific tasks.

We use the `LoraConfig` class to define the settings for LoRA. This configuration specifies how the low-rank adaptations are applied to the model.

```python
peft_config = LoraConfig(
    r=16,
    lora_alpha=16,
    lora_dropout=0.05,
    bias="none",
    task_type="SEQ_CLS",
    target_modules=["q_proj", "k_proj", "v_proj", "o_proj", "gate_proj"]
)
```

> <question-title></question-title>
>
> What do the parameters?
> 
> 1. `r=16`?
> 2. `lora_alpha=16`?
> 3. `lora_dropout=0.05`?
> 4. `bias="none"`?
> 5. `task_type="SEQ_CLS"`?
> 6. `target_modules=["q_proj", "k_proj", "v_proj", "o_proj", "gate_proj"]`?
>
> > <solution-title></solution-title>
> >
> > 1. `r=16`: This parameter specifies the rank of the low-rank matrices used in the adaptation. A higher rank allows the model to capture more complex patterns but also increases the number of trainable parameters.
> >
> > 2. `lora_alpha=16`: This scaling factor controls the magnitude of the updates applied by the low-rank matrices. It helps balance the influence of the adaptations relative to the original model weights.
> >
> > 3. `lora_dropout=0.05`: Dropout is applied to the low-rank matrices during training to prevent overfitting. A dropout rate of 0.05 means that 5% of the elements are randomly set to zero during each training step.
> >
> > 4. `bias="none"`: This setting specifies that no bias parameters are added to the low-rank adaptations. Other options include "all" to add biases to all layers or "lora_only" to add biases only to the LoRA layers.
> >
> > 5. `task_type="SEQ_CLS"`: This indicates that the model is being fine-tuned for a sequence classification task. Other task types might include "CAUSAL_LM" for causal language modeling or "SEQ_2_SEQ_LM" for sequence-to-sequence tasks.
> >
> > 6. `target_modules=["q_proj", "k_proj", "v_proj", "o_proj", "gate_proj"]`: This list specifies the modules within the model architecture to which the LoRA adaptations will be applied. These modules are typically the attention layers in transformer models:
> >    - `"q_proj"`: query projections
> >    - `"k_proj"`: key projections
> >    - `"v_proj"`: value projections
> >    - `"o_proj"`: output projections
> >    - `"gate_proj"`: gating projections in some architectures.
> >
> {: .solution}
>
{: .question}

By configuring LoRA in this way, we can efficiently adapt a large pretrained model to a specific task with minimal computational overhead, making it feasible to fine-tune on consumer-grade hardware. This approach is particularly useful for tasks like text classification, sentiment analysis, or any other application where we need to specialize a general-purpose language model.

## Configure Training Arguments

Let's now set up the training arguments using the `TrainingArguments` class from the Hugging Face Transformers library. These arguments define the training configuration, including hyperparameters and settings for saving and evaluating the model.

```python
training_args = transformers.TrainingArguments(
    output_dir="./results",
    evaluation_strategy="epoch",
    save_strategy="epoch",
    learning_rate=1e-5,
    per_device_train_batch_size=16,
    per_device_eval_batch_size=16,
    num_train_epochs=5,
    weight_decay=0.01,
    bf16=True,
    report_to="none",
    load_best_model_at_end = True,
)
```

> <question-title></question-title>
>
> What do the parameters?
>
> 1. `output_dir="./results"`
> 2. `evaluation_strategy="epoch"`
> 3. `save_strategy="epoch"`
> 4. `learning_rate=1e-5`
> 5. `per_device_train_batch_size=16`
> 6. `per_device_eval_batch_size=16`
> 7. `num_train_epochs=5`
> 8. `weight_decay=0.01`
> 9. `bf16=True`
> 10. `report_to="none"`
> 11. `load_best_model_at_end=True`
> 
> > <solution-title></solution-title>
> >
> > 1. `output_dir="./results"`: Specifies the directory where the model predictions and checkpoints will be saved.
> >
> > 2. `evaluation_strategy="epoch"`: The model will be evaluated at the end of each epoch. This allows for monitoring the model's progress and adjusting the training process as needed.
> >
> > 3. `save_strategy="epoch"`: The model checkpoints will be saved at the end of each epoch.  This ensures that checkpoints are available for each complete pass through the dataset.
> >
> > 4. `learning_rate=1e-5`: Sets the initial learning rate for the optimizer. This rate determines how much the model's weights are updated during training.
> >
> > 5. `per_device_train_batch_size=16`: The number of samples per device (e.g., GPU) to load for training. 
> >
> > 6. `per_device_eval_batch_size=16`: The number of samples per device to load for evaluation.
> >
> > 7. `num_train_epochs=5`: The total number of training epochs. An epoch is one complete pass through the training dataset.
> >
> > 8. `weight_decay=0.01`: Applies L2 regularization to the model weights to prevent overfitting.
> >
> > 9. `bf16=True`: Enables mixed precision training using bfloat16, which can speed up training and reduce memory usage on compatible hardware.
> >
> > 10. `report_to="none"`: Disables reporting to external services like WandB or TensorBoard. If you want to track metrics, you can set this to "wandb", "tensorboard", etc.
> >
> > 11. `load_best_model_at_end=True`: Ensures that the best model based on evaluation metrics is loaded at the end of training.
> >
> {: .solution}
>
{: .question}

These settings provide a balanced configuration for training a model efficiently while ensuring that the best version of the model is saved and can be used for further evaluation or deployment. Adjust these parameters based on your specific use case and available computational resources.



# Prepare the tokenizer

We will now set up the tokenizer to convert DNA sequences into numerical tokens that the model can process. The tokenizer is a crucial component in preparing the data for model training and inference: it transforms raw text into a format that can be processed by machine learning models. 

We use the `AutoTokenizer` class from the Hugging Face Transformers library to load a pre-trained tokenizer. We specify the pre-trained model from which to load the tokenizer. This should match the model you plan to use for training or inference. This tokenizer will be configured to handle DNA sequences efficiently.

```python
tokenizer = transformers.AutoTokenizer.from_pretrained(
    model_name,
    model_max_length=200,
    padding_side="right",
    use_fast=True,
    trust_remote_code=True,
)
```

> <question-title></question-title>
>
> What do the parameters?
>
> 1. `model_max_length=200`
> 2. `padding_side="right"`
> 3. `use_fast=True`
> 4. `trust_remote_code=True`
> 
> > <solution-title></solution-title>
> >
> > 1. `model_max_length=200`: Sets the maximum length of the tokenized sequences. Sequences longer than this will be truncated, and shorter ones will be padded.
> > 
> > 2. `padding_side="right"`: Specifies that padding should be added to the right side of the sequences. This ensures that all sequences in a batch have the same length.
> >
> > 3. `use_fast=True`: Enables the use of the fast tokenizer implementation, which is optimized for speed and is suitable for most use cases.
> >
> > 4. `trust_remote_code=True`: Allows the tokenizer to execute custom code from the model repository, which may be necessary for some models that require specific preprocessing steps.
> >
> {: .solution}
>
{: .question}

By configuring the tokenizer in this way, we ensure that our DNA sequences are properly tokenized and formatted for input into the model. This step is essential for preparing our data for efficient and effective model training and evaluation.

Let's now tailor the tokenizer to better suit our specific use case, ensuring that the model processes sequences accurately and efficiently. Special tokens play a crucial role in defining how sequences are processed and interpreted by the model. Here, we sets:
- the end-of-sequence (EOS) token, which indicates the end of a sequence. It is essential for tasks where the model needs to generate sequences or understand where a sequence ends.
- the padding (PAD) token, which is used to pad sequences to a uniform length within a batch. Padding ensures that all sequences in a batch have the same length, which is necessary for efficient processing during training and inference

```
tokenizer.eos_token = "[EOS]"
tokenizer.pad_token = "[PAD]"
```

# Prepare data

To finetune the model, we must provide a dataset to train the model. We will the data with the 1st transcription factor (`tf0`) in mouse from {% cite zhou2024dnabert2efficientfoundationmodel %}. The data is stored on [GitHub](https://github.com/raphaelmourad/Mistral-DNA).

## Get data

Let's get the data for from GitHub:

```python
!git clone https://github.com/raphaelmourad/Mistral-DNA.git
```

We now need to uncompress the labeled data:

```python
!tar -xf Mistral-DNA/data/GUE.tar.xz -C Mistral-DNA/data/
```

We change the current working directory to the `Mistral-DNA` folder.

```python
os.chdir("Mistral-DNA/")
print(os.getcwd())
```

Let's define experience and path to data variables

```python
expe = "tf/0"
data_path = f"data/GUE/{ expe }" 
```




## Prepare Datasets for Training and Validation

We now need to set up the datasets required for training and validating. Properly preparing these datasets is crucial for ensuring that the model finetunes effectively and generalizes well to new data.

We will use the files `data_path` folder we just defined:
- `train.csv` for training
- `dev.csv` for validation

> <question-title></question-title>
>
> How is the content of each file?
>
> 1. `train.csv`
> 2. `dev.csv`
> 
> > <solution-title></solution-title>
> >
> > The 2 files are CSV files with 2 columns (`sequence` and `label`) and different number of rows:
> > 1. `train.csv`: 32,379 rows. 
> > 2. `dev.csv`: 1,000 rows
> > 
> > Values in `label` are:
> > - `0`: The DNA sequence in `sequence` column does not bind to the 1st transcription factor.
> > - `1`: The DNA sequence in `sequence` column binds to the transcription factor.
> >
> {: .solution}
>
{: .question}


Before we proceed we import some classes and functions from `scriptPython/function.py`:

```python
### LOAD FUNCTIONS MODULE
import sys
sys.path.append("scriptPython/")
from functions import *
```

We use the `SupervisedDataset` class to load and prepare the datasets. This class handles the tokenization and formatting of the data, making it ready for model training and evaluation.

```python
train_dataset = SupervisedDataset(
    tokenizer=tokenizer,
    data_path=Path(data_path) / "train.csv",
    kmer=-1,
)
val_dataset = SupervisedDataset(
    tokenizer=tokenizer,
    data_path=Path(data_path) / "dev.csv",
    kmer=-1,
)
```

> <question-title></question-title>
>
> What does `kmer=-1`?
>
> > <solution-title></solution-title>
> >
> > This parameter is used to specify the length of k-mers (substrings of length k) to be considered in the dataset. A value of -1 typically means that no k-mer splitting is applied, and the sequences are processed as they are.
> >
> {: .solution}
>
{: .question}

## Configure Data Collation

A data collator ensures that sequences are properly padded and formatted, which is crucial for optimizing the training process.

We'll use the `DataCollatorForSupervisedDataset` class to handle the collation of tokenized data. This collator will manage padding and ensure that all sequences in a batch are of uniform length.

```python
data_collator = DataCollatorForSupervisedDataset(tokenizer=tokenizer)
```

# Load and Configure the Model for Sequence Classification

Let's now load the pre-trained model, a model originally trained for large language modeling tasks, not specifically for classification. To adapt it for our binary classification task, we will add a new classification head on top of the existing architecture. This head will consist of a single neuron that connects to the output of the language model, enabling it to classify whether a DNA sequence binds to a transcription factor (label `1`) or not (label `0`).

This additional layer, or **classification head**, is a simple neural network layer that takes the high-level features extracted by the language model and maps them to our binary classification output. It learns to weigh these features appropriately to make accurate predictions for our specific task.

We use the `AutoModelForSequenceClassification` class from the Hugging Face Transformers library to load the pre-trained model and set it up for our specific classification task:

```python
model=transformers.AutoModelForSequenceClassification.from_pretrained(
    model_name,
    num_labels=2,
    output_hidden_states=False,
    quantization_config=bnb_config,
    device_map="auto",
    trust_remote_code=True,
)
```

> <question-title></question-title>
>
> What do the parameters?
>
> 1. `num_labels=2`
> 2. `output_hidden_states=False`
> 3. `quantization_config=bnb_config`
> 4. `device_map="auto"`
> 5. `trust_remote_code=True`
> 
> > <solution-title></solution-title>
> >
> > 1. `num_labels=2`: Sets the number of output labels to 2, corresponding to the binary classification task (binding or not binding to transcription factors).
> >
> > 2. `output_hidden_states=False`: Indicates that the model should not output hidden states. This is typically set to False unless you need access to the intermediate representations for further analysis.
> > 
> > 3. `quantization_config=bnb_config`: Applies predefined quantization configuration to the model, which helps reduce memory usage and enables efficient training on consumer-grade hardware.
> > 
> > 4. `device_map="auto"`: Automatically determines the best device placement for the model's layers, optimizing for available hardware (e.g., GPUs). If it finds a GPU, it will use a GPU. If there's no GPU, it will not use the GPU
> > 
> > 5. `trust_remote_code=True`: Allows the model to execute custom code from the model repository, which may be necessary for certain architectures or preprocessing steps.
> >
> {: .solution}
>
{: .question}

To ensure that the model correctly handles padding tokens, we need to align the padding token configuration between the model and the tokenizer. This step is crucial for maintaining consistency during training and inference, especially when dealing with sequences of varying lengths:

```python
model.config.pad_token_id = tokenizer.pad_token_id
```

# Initialize the Trainer

We can now set up the `Trainer` to manage the training and evaluation process of our model. The `Trainer` class simplifies the training loop, handling many of the complexities involved in training deep learning models.

We first need to attach the LoRA adapter to the model:

```python
model.add_adapter(peft_config, adapter_name="lora_1")
```

Let's now set up the `Trainer`:

```python
trainer = transformers.Trainer(
    model=model,
    args=training_args,
    compute_metrics=compute_metrics,
    train_dataset=train_dataset,
    eval_dataset=val_dataset,
    data_collator=data_collator,
    callbacks = [EarlyStoppingCallback(early_stopping_patience=3)]
)
```

> <question-title></question-title>
>
> What do the `callbacks = [EarlyStoppingCallback(early_stopping_patience=3)]` parameter?
> 
> > <solution-title></solution-title>
> >
> > It adds an early stopping mechanism to the training process. This mechanism is designed to halt training when the model's performance on the validation set stops improving, helping to prevent overfitting and conserve computational resources.
> > 
> > How Early Stopping Works?
> > 
> > **Purpose**: The primary goal of early stopping is to capture the model parameters when the loss reaches its minimum value during training. This is crucial because, after a certain point, continued training may lead to overfitting, where the model starts to perform worse on unseen data.
> > 
> > **Patience Parameter**: The `early_stopping_patience=3` setting specifies that training should continue for three additional epochs after the model's performance on the validation set stops improving. This "patience" period helps mitigate the effects of noise in the training process. Noise can cause temporary fluctuations in the loss, making it seem like the model has reached a local minimum when further training might yield better results.
> >
> > **Process**: During training, the loss is monitored at each epoch. If the loss does not decrease for three consecutive epochs, training is stopped. However, if a better model with a lower loss is found within those three epochs, training continues. This approach ensures that the model has truly reached a robust local minimum, rather than being prematurely halted due to noise.
> >
> > By incorporating early stopping with a patience of three epochs, you balance the need to find an optimal model with the risk of overfitting, ultimately leading to more efficient and effective training outcomes.
> > 
> {: .solution}
>
{: .question}

Ffor distributed training, where multiple GPUs or nodes are used to accelerate the training process, it is essential to do:

```python
trainer.local_rank=training_args.local_rank
```

The `local_rank` parameter identifies the rank of the current process within its local node, enabling coordinated communication and synchronization between processes. This setup is crucial for managing tasks such as gradient synchronization and data partitioning, ensuring that each process operates on the correct portion of the model or dataset. By assigning the local rank from `training_args` to the `Trainer`, we facilitate efficient and scalable training, leveraging the full computational power of multi-GPU environments.

# Start the training

Let's start the training process for our model using the trainer.train() method:

```python
trainer.train()
```

After launching `trainer.train()`, we can notice that the training process is significantly faster compared to training a model from scratch seen in ["" tutorial]({% link topics/statistics/tutorials/genomic-llm-pretraining/tutorial.md %}). This efficiency is due to the use of a pre-trained model, which has already undergone extensive training on large datasets using powerful computational resources. For example, pre-training a model on even a small portion of the human genome can take dozens of hours, but fine-tuning this model on a specific task, such as classifying DNA sequences, is much quicker. Fine-tuning leverages the pre-trained model's foundational knowledge, allowing you to adapt it to new tasks with a smaller, labeled dataset. This approach not only saves time but also reduces the need for extensive computational power. By downloading a pre-trained model from platforms like Hugging Face and fine-tuning it on a local machine with a modest GPU, we can achieve high performance with minimal overhead, making advanced modeling techniques accessible for a wide range of applications.

# Evaluate Model Performance

After successfully training the model, the next essential step is to evaluate its performance on a test dataset. This evaluation process is crucial for understanding how well the model generalizes to new, unseen data and for assessing its readiness for real-world applications.

> <comment-title></comment-title>
>
> If finetuning is too long, you can stop the training.
>
{: .comment}

The test data is stored in `data_path/test.csv`, we prepare it as for training and validation data.

```python
test_dataset = SupervisedDataset(
    tokenizer=tokenizer,
    data_path=Path(data_path) / "test.csv",
    kmer=-1,
)
```

We then use the `trainer.evaluate()` method. This methods is designed to assess the model's performance on a specified dataset, typically the test dataset, which contains data that the model has not encountered during training.

```python
results = trainer.evaluate(eval_dataset=test_dataset)
```

The method computes various evaluation metrics, such as accuracy, precision, recall, and F1 score, depending on the task and the configuration specified in `compute_metrics`. These metrics provide a comprehensive view of the model's performance, highlighting its strengths and weaknesses.

The Trainer uses the `data_collator` to ensure that the test data is properly formatted and padded, maintaining consistency with the training process. This consistency is crucial for accurate evaluation.

The evaluation results are stored in the `results` variable, which contains the computed metrics. We can analyze these `results` to gain insights into the model's performance and make informed decisions about further improvements or deployment.

> <question-title></question-title>
>
> What is stored in `results`? How do you interpret this information?
> 
> > <solution-title></solution-title>
> >
> > `results` provides a comprehensive overview of the model's performance on the evaluation dataset with:
> > 1. **eval_loss (0.424961)**: This metric represents the loss value calculated on the evaluation dataset. Lower values indicate better model performance.
> >
> >    A loss of 0.425 suggests that the model is reasonably well-fitted to the data, though the specific interpretation depends on the context and the loss function used (e.g., cross-entropy for classification tasks).
> > 
> > 2. **eval_accuracy (0.804000)**: Accuracy measures the proportion of correctly predicted instances out of the total instances.
> > 
> >    An accuracy of 80.4% indicates that the model correctly predicted the class for 80.4% of the samples in the evaluation dataset.
> > 
> > 4. **eval_f1 (0.800838)**: The F1 score is the harmonic mean of precision and recall, providing a single metric that balances both concerns.
> >    
> >    An F1 score of 0.801 suggests a good balance between precision and recall, indicating that the model performs well in both identifying positive cases and minimizing false positives and negatives.
> > 
> > 4. **eval_matthews_correlation (0.628276)**: The Matthews Correlation Coefficient (MCC) is a measure of the quality of binary classifications, taking into account true and false positives and negatives.
> >    
> >    An MCC of 0.628 indicates a moderate to strong correlation between the predicted and actual classes, suggesting the model is performing better than random guessing.
> > 
> > 6. **eval_precision (0.824614)**: Precision is the ratio of correctly predicted positive observations to the total predicted positives.
> > 
> >    A precision of 82.5% means that out of all the instances predicted as positive, 82.5% were actually positive.
> > 
> > 6. **eval_recall (0.804000)**: Recall (or sensitivity) is the ratio of correctly predicted positive observations to all observations in the actual class.
> >    
> >    A recall of 80.4% indicates that the model correctly identified 80.4% of all actual positive cases.
> > 
> > 8. **eval_runtime (6.548800)**: The total time taken to evaluate the model on the dataset.
> >    
> >    A runtime of 6.55 seconds provides insight into the computational efficiency of the evaluation process.
> > 
> > 8. **eval_samples_per_second (152.699000)**: The number of samples processed per second during evaluation.
> > 
> >    Processing 152.7 samples per second indicates the efficiency of the evaluation pipeline.
> > 
> > 9. **eval_steps_per_second (9.620000)**: The number of evaluation steps completed per second.
> >     
> >    Completing 9.62 steps per second reflects the speed of the evaluation process.
> > 
> > 11. **epoch (3.000000)**: The number of training epochs completed before this evaluation.
> > 
> >    The evaluation was conducted after 3 epochs of training, providing context for the model's learning progress.
> > 
> {: .solution}
>
{: .question}


# Conclusion

In this tutorial, we explored the process of fine-tuning a large language model (LLM) for DNA sequence classification. By following the steps outlined, you have learned how to leverage pre-trained models to achieve efficient and effective classification of DNA sequences, specifically focusing on their binding affinity to transcription factors.

We began by configuring the fine-tuning process, ensuring that available computational resources were optimally utilized. This included specifying settings for quantization, configuring Accelerate for distributed training, and implementing LoRA for parameter-efficient fine-tuning. These steps were crucial for maximizing performance and minimizing computational overhead.

Next, we prepared the tokenizer and data, ensuring that DNA sequences were properly tokenized and formatted for model input. We created datasets for training and validation, and configured data collation to handle batch processing efficiently.

We then loaded and configured the model for sequence classification, adding a classification head to adapt the pre-trained model to our specific task. With the model and data prepared, we initialized the `Trainer`, which streamlined the training process by managing the training loop, evaluation, and checkpointing.
