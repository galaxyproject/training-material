---
layout: tutorial_hands_on
title: Predicting Mutation Impact with Zero-shot Learning using a pretrained DNA LLM
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
  - genomic-llm-finetuning
questions:
- How does zero-shot learning differ from traditional supervised learning, and what advantages does it offer in the context of predicting DNA mutation impacts?
- What steps are involved in computing embeddings for DNA sequences using a pre-trained LLM, and how do these embeddings capture the semantic meaning of the sequences?
- Why is the L2 distance used as a metric to quantify the impact of mutations, and how does a higher L2 distance indicate a more significant mutation effect?
objectives:
- Explain the concept of zero-shot learning and its application in predicting the impact of DNA mutations using pre-trained large language models (LLMs).
- Utilize a pre-trained DNA LLM from Hugging Face to compute embeddings for wild-type and mutated DNA sequences.
- Compare the embeddings of wild-type and mutated sequences to quantify the impact of mutations using L2 distance as a metric.
- Interpret the results of the L2 distance calculations to determine the significance of mutation effects and discuss potential implications in genomics research.
- Develop a script to automate the process of predicting mutation impacts using zero-shot learning, enabling researchers to apply this method to their own datasets efficiently.
time_estimation: 3H
key_points:
- Pre-trained DNA language models are powerful tools for analyzing genetic sequences, enabling efficient and accurate assessment of mutation impacts without extensive computational resources.
- SNPs in exons generally have a more significant impact on gene function compared to those in introns, highlighting the importance of focusing on coding regions for understanding genetic variations.
- Utilizing embeddings to quantify the effects of mutations provides a robust method for comparing sequence variations, offering insights into the functional consequences of SNPs and other genetic modifications.
- Employing statistical tests, such as t-tests and Wilcoxon rank-sum tests, is crucial for validating observed differences in genetic data, ensuring that findings are supported by empirical evidence.
contributions:
  authorship:
  - raphaelmourad
  - bebatut
tags:
- elixir
- ai-ml
- Large Language Model
subtopic: gai-llm
priority: 3
notebook:
  language: python
  pyolite: true
---

Predicting the impact of mutations is a critical task in genomics, as it provides insights into how genetic variations influence biological functions and contribute to diseases. Traditional methods for assessing mutation impact often rely on extensive experimental data or computationally intensive simulations. However, with the advent of large language models (LLMs) and zero-shot learning, we can now predict mutation impacts more efficiently and effectively.

Zero-shot learning is a technique that allows a pre-trained model to make predictions on tasks it wasn't explicitly trained for, leveraging its existing knowledge. This approach is particularly valuable when labeled data is scarce or when rapid predictions are needed. By using a pre-trained DNA LLM, we can compute embeddings for both wild-type and mutated DNA sequences and compare them to quantify the impact of mutations.

This tutorial focuses on this innovative method, utilizing a [pre-trained DNA LLM available on Hugging Face](https://huggingface.co/RaphaelMourad/Mistral-DNA-v1-17M-hg38) to assess the impact of mutations. This approach opens new avenues for bioinformatics research, particularly in genomics and personalized medicine, by enabling researchers to gain insights into the functional impact of DNA mutations efficiently.

Building upon this foundation, our new tutorial focuses on predicting the impact of mutations using zero-shot learning with a pre-trained DNA LLM. Zero-shot learning allows us to utilize the pre-trained model directly, without additional training, to make predictions on new, unseen tasks. Specifically, we will use a pre-trained model available on Hugging Face, designed for DNA sequences, to assess the impact of mutations.

We will use [`Mistral-DNA-v1-17M-hg38`](https://huggingface.co/RaphaelMourad/Mistral-DNA-v1-1M-hg38), a mixed model that was pre-trained on the entire Human Genome. It contains approximately 17 million parameters and was trained using the Human Genome assembly GRCh38 on sequences of 10,000 bases (10K):

```python
model_name="RaphaelMourad/Mistral-DNA-v1-17M-hg38"
```

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
!pip install Bio==1.7.1
!pip install transformers -U
```

## Import Python libraries

Let's now import them.

```python
import os

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import torch
import tensorflow as tf
import gzip
from Bio import SeqIO
from transformers import (
  AutoConfig,
  AutoModelForCausalLM,
  AutoTokenizer,
  EarlyStoppingCallback,
  Trainer,
  TrainingArguments,
)
```

> <comment-title>Versions</comment-title>
>
> This tutorial has been tested with following versions:
> - `transformers` = 4.48.3
> 
> You can check the versions with:
>
> ```python
> transformers.__version__
> ```
{: .comment}

## Check and configure available resources

We select the appropriate device (CUDA-enabled GPU if available) for running PyTorch operations

```python
torch.device('cuda' if torch.cuda.is_available() else 'cpu')
```

Let's check the GPU usage and RAM:

```python
!nvidia-smi
```

Let's configure PyTorch and the CUDA environment -- software and hardware ecosystem provided by NVIDIA to enable parallel computing on GPU -- to optimize GPU memory usage and performance:

1. Enables CuDNN benchmarking in PyTorch:

    ```python
    torch.backends.cudnn.benchmark=True
    ```

2. Set an environment variable that configures how PyTorch manages CUDA memory allocations

    ```python
    os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:32"
    ```


# Tokenizing DNA Sequences

We now set up the tokenizer to convert raw DNA sequences into a format that the model can process, enabling it to understand and analyze the sequences effectively:

```python
tokenizer = transformers.AutoTokenizer.from_pretrained(
    model_name,
    use_fast=True,
    trust_remote_code=True,
)
```

> <question-title></question-title>
>
> What do the parameters?
> 
> 1. `use_fast=True`?
> 2. `trust_remote_code=True`?
>
> > <solution-title></solution-title>
> >
> > 1. `use_fast=True`: Enables the use of a fast tokenizer implementation, which is optimized for speed and efficiency. This is particularly useful when working with large datasets or when performance is a priority.
> > 2. `trust_remote_code=True`: Allows the tokenizer to execute custom code from the model repository. This may be necessary for certain architectures or preprocessing steps that require additional functionality.
> >
> {: .solution}
>
{: .question}

# Load and Configure the Pre-trained Model

We will now load the pre-trained DNA large language model (LLM) and configure it for our specific task of predicting the impact of DNA mutations. 

```python
model=transformers.AutoModelForCausalLM.from_pretrained(
    model_name,
)
```

We would like to ensure that the model correctly handles padding tokens, which are used to standardize the length of sequences within a batch:

```python
model.config.pad_token_id = tokenizer.pad_token_id
```

Aligning the padding token ID between the model and tokenizer is crucial for maintaining consistency during training and inference.

Let's look at the model:

```python
model
```

```
MixtralForCausalLM(
  (model): MixtralModel(
    (embed_tokens): Embedding(4096, 256)
    (layers): ModuleList(
      (0-7): 8 x MixtralDecoderLayer(
        (self_attn): MixtralAttention(
          (q_proj): Linear(in_features=256, out_features=256, bias=False)
          (k_proj): Linear(in_features=256, out_features=256, bias=False)
          (v_proj): Linear(in_features=256, out_features=256, bias=False)
          (o_proj): Linear(in_features=256, out_features=256, bias=False)
        )
        (block_sparse_moe): MixtralSparseMoeBlock(
          (gate): Linear(in_features=256, out_features=8, bias=False)
          (experts): ModuleList(
            (0-7): 8 x MixtralBlockSparseTop2MLP(
              (w1): Linear(in_features=256, out_features=256, bias=False)
              (w2): Linear(in_features=256, out_features=256, bias=False)
              (w3): Linear(in_features=256, out_features=256, bias=False)
              (act_fn): SiLU()
            )
          )
        )
        (input_layernorm): MixtralRMSNorm((256,), eps=1e-05)
        (post_attention_layernorm): MixtralRMSNorm((256,), eps=1e-05)
      )
    )
    (norm): MixtralRMSNorm((256,), eps=1e-05)
    (rotary_emb): MixtralRotaryEmbedding()
  )
  (lm_head): Linear(in_features=256, out_features=4096, bias=False)
)
```

> <question-title></question-title>
>
> What do the parameters?
> 
> 1. How does the vocabulary size of 4,096 relate to k-mers in DNA sequences?
> 2. What is the role of the embedding layer in `MixtralForCausalLM`?
> 3. How many layers does the `MixtralForCausalLM` model have, and what is their purpose?
> 4. What components make up the self-attention mechanism in MixtralAttention?
> 5. How does the Mixture of Experts (MoE) mechanism reduce computational load?
> 6. Why is layer normalization important in `MixtralForCausalLM`?
> 7. What advantage do rotary embeddings offer in understanding DNA sequences?
>
> > <solution-title></solution-title>
> >
> > 1. The vocabulary size corresponds to the number of unique "words" or k-mers (subsequences of DNA) the model can recognize, similar to using k-mers of size six, enhanced by byte-pair encoding for more nuanced patterns.
> > 2. The embedding layer converts DNA sequences into numerical vectors, enabling the model to process and analyze them.
> > 3. The model has 8 layers, each containing a `MixtralDecoderLayer`. These layers process the embedded input sequences through a series of transformations, including self-attention and mixture of experts, to capture complex patterns in the data.
> > 4. The self-attention mechanism includes query, key, value, and output projections, which weigh the importance of different tokens in the sequence.
> > 5. MoE activates only a subset of parameters during each forward pass, using a routing mechanism to direct sequences to specific experts.
> > 6. Layer normalization stabilizes and accelerates training by ensuring consistent scaling of inputs to the attention mechanism and subsequent layers.
> > 7. Rotary embeddings enhance the model's understanding of positional information, providing a more nuanced representation of sequence order compared to traditional methods.
> >
> {: .solution}
>
{: .question}

> <question-title></question-title> 
>
> For a DNA sequence "ACGTAGCATCGGATCTATCTATCGACACTTGGTTATCGATCTACGAGCATCTCGTTAGC"
>
> 1. How can we get the hidden states?
> 2. How can we compute mean of the hidden states accross the sequence length dimension?
> 3. What is the shape of the output?
> 4. What does the mean of the hidden states accross the sequence length dimension represent?
>
> > <solution-title></solution-title>
> >
> > Let's start by defining the DNA sequence:
> > 
> > ```python
> > dna = "ACGTAGCATCGGATCTATCTATCGACACTTGGTTATCGATCTACGAGCATCTCGTTAGC"
> > ```
> > 
> > 1. To get the hidden state:
> >    1. Tokenize the DNA sequence using the tokenizer
> >      
> >        ```python
> >        tokenized_dna = tokenizer(dna, return_tensors = "pt")
> >        ```
> >
> >    2. Extract the tensor containing the token IDs from the tokenized output
> >      
> >        ```python
> >        inputs = tokenized_dna["input_ids"]
> >        ```
> >
> >    3. Pass the tokenized input through the model.
> >      
> >        ```python
> >        model_outputs = model(inputs)
> >        ```
> >
> >    4. Extract the hidden states from the model's output.
> >      
> >        ```python
> >        hidden_states = model_outputs[0].detach()
> >        ```
> >
> > 2. To compute mean of the hidden states accross the sequence length dimension:
> >    
> >    ```python
> >    embedding_mean = torch.mean(hidden_states[0], dim=0)
> >    ```
> >
> > 3. The shape is 4,096, the number of possible tokens
> > 
> > 4. It represents the average embedding of the DNA sequence. This fixed-size representation can be used for various downstream tasks, such as classification, clustering, or similarity comparisons.
> {: .solution}
>
{: .question}

# Compare the effect of mutations with or without amino acid modification

Let's look at the portion of the Cystic fibrosis transmembrane conductance regulator (CFTR) gene where a mutation responsible of the Cystic fibrosis appears.

In this portion, we can observe several cases:

Cases | Sequences
--- | ---
Wild-type without any mutation | ATTAAAGAAAATATCATCTTTGGTGTTTCCTAT
Mutation ATT to ATA | ATAAAAGAAAATATCATCTTTGGTGTTTCCTAT
Deletion of TCT codon | ATTAAAGAAAATATCA---TTGGTGTTTCCTAT

In the second case, the amino acid does not change (silent mutation). In the last case, an amino acid is removed. Let's look if the mutation and the deletion have the same distance to the wild-type when computing using the DNA embedding.

First, we define the sequences:

```python
dna_wt= "ATTAAAGAAAATATCATCTTTGGTGTTTCCTAT"
dna_mut="ATAAAAGAAAATATCATCTTTGGTGTTTCCTAT"
dna_del="ATTAAAGAAAATATCATTGGTGTTTCCTAT"
```

Let's compute the hidden states for all DNA sequences:

```python
tokenized_dna = tokenizer(
  [dna_wt, dna_mut, dna_del],
  return_tensors="pt",
  padding=True,
)
inputs_seqs = tokenized_dna["input_ids"]
model_outputs = model(inputs_seqs)
hidden_states = model_outputs[0].detach()
```

We now compute the maximum of the hidden states accross the sequence length dimension:

```python
embedding_max = torch.max(hidden_states, dim=1)[0]
```


To compare the effects of silent mutation and amino acid deletion, we will compute the distance between the wild-type embeddings and the mutation / deletion embedding using the L2 (Euclidean) distance

> <details-title>L2 (Euclidean) distance</details-title>
>
> The L2 distance, also known as the Euclidean distance, is a measure of the straight-line distance between two points in Euclidean space. It is commonly used to quantify the difference between two vectors, such as embeddings in machine learning.
>
> For two vectors $$a$$ and $$b$$ in an $$n$$-dimensional space, where $$a=[a_{1},a_{2},...,a_{n}]$$ and $$b=[b_{1},b_{2},...,b_{n}]$$, the L2 distance is calculated as:
>
> \\(L2 = \sqrt{\sum_{i} (a_{i}-b_{i})^{2}\\)
>
{: .details}
          
```python
wt_mut_L2 = torch.norm(embedding_max[0] - embedding_max[1])
print(wt_mut_L2)
wt_del_L2 = torch.norm(embedding_max[0] - embedding_max[2])
print(wt_del_L2)
```

```
tensor(27.4657)
tensor(145.7797)
```

> <question-title></question-title>
>
> 1. Is it the silent mutation or the amino acid deletion that have the lowest L2 distance to the wild-type?
> 2. How to interpret the result?
>
> > <solution-title></solution-title>
> >
> > 1. The silent mutation has a lower L2 distance to the wild-type compared to the mutation with amino acid deletion.
> > 2. A smaller L2 distance indicates that the sequences are more similar, while a larger distance suggests greater dissimilarity.
> {: .solution}
>
{: .question}


# Compare effects of SNPs in exons and introns

Single Nucleotide Polymorphisms (SNPs) are the most common type of genetic variation, involving a change in a single nucleotide within a DNA sequence. These variations can occur in different regions of the genome, including exons (the coding regions of genes) and introns (the non-coding regions within genes). Understanding the impact of SNPs in these regions is crucial for assessing their role in genetic disorders, phenotypic traits, and overall genome function.

We will now leverage the pre-trained DNA language model to compare the effects of SNPs in exons and introns. By utilizing embeddings generated from DNA sequences, we can quantify the impact of these variations and gain insights into their functional consequences.

## Load sequences

To begin our analysis, we need to load the sequence data, which includes sequences with Single Nucleotide Polymorphisms (SNPs) and their corresponding reference sequences (wild-type) without SNPs. These sequences are derived from both introns and exons and are stored in compressed FASTA format on [GitHub](https://github.com/raphaelmourad/Mistral-DNA/tree/main/data/SNP):

| Wild-type (without SNPs) | Mutated (with SNP)
--- | --- | ---
Intron | [`SNPintron_ref_201b.fasta.gz`](https://github.com/raphaelmourad/Mistral-DNA/raw/refs/head/master/data/SNP/SNPintron_ref_201b.fasta.gz)  |  [`SNPintron_alt_201b.fasta.gz`](https://github.com/raphaelmourad/Mistral-DNA/raw/refs/head/master/data/SNP/SNPintron_alt_201b.fasta.gz)
Exon | [`SNPexon_ref_201b.fasta.gz`](https://github.com/raphaelmourad/Mistral-DNA/raw/refs/head/master/data/SNP/SNPexon_ref_201b.fasta.gz)  |  [`SNPexon_alt_201b.fasta.gz`](https://github.com/raphaelmourad/Mistral-DNA/raw/refs/head/master/data/SNP/SNPexon_alt_201b.fasta.gz)

```python
baseurl = "https://github.com/raphaelmourad/Mistral-DNA/raw/refs/head/master/data/SNP/"
win=201
exon_wt_snp_fp = f"{ baseurl }/SNPexon_ref_{ win }b.fasta.gz"
exon_mut_snp_fp = f"{ baseurl }/SNPexon_alt_{ win }b.fasta.gz"
intron_wt_snp_fp = f"{ baseurl }/SNPintron_ref_{ win }b.fasta.gz"
intron_mut_snp_fp = f"{ baseurl }/SNPintron_alt_{ win }b.fasta.gz"
```

We need to get files and read the sequences from them:

```python
import requests
def downloadReadFastaFile(fasta_file):
  response = requests.get(fasta_file)
  # Check if the request was successful
  if response.status_code == 200:
      # Open a local file in binary write mode and save the content
      with open('file.gz', 'wb') as file:
          file.write(response.content)
      print("File downloaded successfully.")
  else:
      print(f"Failed to download file. HTTP Status code: {response.status_code}")
  # Read the file
  seql_list=[]
  with gzip.open("file.gz", "rt") as handle:
      for record in SeqIO.parse(handle, "fasta"):
          seqj=str(record.seq)
          seql_list.append(seqj)
  return seql_list

exon_wt_seqs = downloadReadFastaFile(exon_wt_snp_fp)
exon_mut_seqs = downloadReadFastaFile(exon_mut_snp_fp)
intron_wt_seqs = downloadReadFastaFile(intron_wt_snp_fp)
intron_mut_seqs = downloadReadFastaFile(intron_mut_snp_fp)
```

> <question-title></question-title>
>
> 1. How many sequences are in the data?
> 2. How long are the sequences?
>
> > <solution-title></solution-title>
> >
> > 1. ```python 
> > len( exon_wt_seqs ) , len(exon_mut_seqs) , len(intron_wt_seqs) , len( intron_mut_seqs ) 
> > ```
> >  `(10000, 10000, 10000, 10000)`
> > 2. the sequences are 201 nucleotides long.
> {: .solution}
>
{: .question}

We will keep only the 100 first sequences:

```python
kseq = 100
exon_wt_seqs = exon_wt_seqs[0:kseq]
exon_mut_seqs = exon_mut_seqs[0:kseq]
intron_wt_seqs = intron_wt_seqs[0:kseq]
intron_mut_seqs = intron_mut_seqs[0:kseq]
```

> <question-title></question-title>
>
> How many differences (SNPs) between reference and alternative sequences are there for first exon sequence?
>
> > <solution-title></solution-title>
> >
> > ```python
> > n_diff = 0
> > for i in range( len( exon_wt_seqs[0] ) ) :
> >     n_diff += exon_mut_seqs[0][i] == exon_wt_seqs[0][i]
> > print(n_diff)
> > ```
> > 
> {: .solution}
>
{: .question}

## Compute effect of SNPs

To compute the effect of SNPs, we need :

1. Compute the embdedding of DNA sequences

  ```python
  def computeEmbedding(seqs):
    tokenized_dna = tokenizer(
      seqs,
      return_tensors="pt",
      padding=True,
    )
    inputs_seqs = tokenized_dna["input_ids"]
    hidden_states = model(inputs_seqs)[0].detach().cpu()
    return torch.max(hidden_states, dim=1)[0]
  ```

2. Compute the L2 distance between reference (without SNP) and alternative (with SNP):

  ```python
  def computeMutationEffect(wt_seqs, mut_seqs):
    wt_embedding = computeEmbedding(wt_seqs)
    mut_embedding = computeEmbedding(mut_seqs)
    return torch.norm(mut_embedding - wt_embedding, dim=1)
  ```


```python
exon_SNP_distL2 = computeMutationEffect(exon_wt_seqs, exon_mut_seqs)
intron_SNP_distL2 = computeMutationEffect(intron_wt_seqs,intron_mut_seqs)
```

> <question-title></question-title>
>
> What are the dimensions of `exon_SNP_distL2` and `intron_SNP_distL2`?
>
> > <solution-title></solution-title>
> >
> > 10,000 each: one value per pair of wt/mutated sequence
> > 
> {: .solution}
>
{: .question}


## Compare SNPs effects between introns and exons

We can now quantify the impact of SNPs and determine if the differences are statistically significant.

To visualize the predicted effects of SNPs, we can use a boxplot to compare the L2 distances for SNPs in exons and introns. The L2 distance serves as a metric for the impact of mutations, with higher distances indicating a more significant effect.

Boxplot of predicted SNP effects

```python
SNPs_distL2 = {"exons": exon_SNP_distL2, "introns": intron_SNP_distL2}

fig, ax = plt.subplots()
ax.boxplot(SNPs_distL2.values())
ax.set_xticklabels(SNPs_distL2.keys())
```

> <question-title></question-title>
>
> ![Box plot comparing two categories: exonSNPs and intronSNPs. The plot displays the distribution of values for each category, with the median indicated by an orange line within each box. The boxes represent the interquartile range, while the whiskers extend to show the range of the data, excluding outliers, which are depicted as individual points above the whiskers. The exonSNPs category shows a higher median and more variability compared to the intronSNPs category.](./images/boxplot.png)
>
> What can we conclude from the boxplot?
>
> > <solution-title></solution-title>
> >
> > From the plot, we observe that the L2 distance is generally higher for SNPs in exons compared to those in introns. This suggests that SNPs in exons have a stronger predicted effect, which aligns with their role in coding regions of genes.
> {: .solution}
>
{: .question}

To determine if the observed differences in L2 distances are statistically significant, we can perform hypothesis tests to compare the distributions:

- T-test:

  ```python
  sp.stats.ttest_ind(exon_SNP_distL2, intron_SNP_distL2)
  ```

  > <question-title></question-title>
  >
  > What can we conclude from the t-test?
  >
  > > <solution-title></solution-title>
  > >
  > > The t-test yields a p-value of approximately $$7.56 \cdot 10^{−7}$$, indicating a statistically significant difference between the L2 distances of SNPs in exons and introns.
  > {: .solution}
  >
  {: .question}

- Wilcoxon rank-sum test:

  ```python
  sp.stats.wilcoxon(exon_SNP_distL2, intron_SNP_distL2)
  ```

  > <question-title></question-title>
  >
  > What can we conclude from the Wilcoxon test?
  >
  > > <solution-title></solution-title>
  > >
  > > The Wilcoxon rank-sum test also shows a significant p-value of approximately $$2.87 \cdot 10^{−6}$$, confirming that the distributions are significantly different.
  > {: .solution}
  >
  {: .question}

The analysis demonstrates that SNPs in exons have a more substantial impact on the sequence embeddings compared to SNPs in introns, as evidenced by higher L2 distances. This finding is statistically significant, highlighting the importance of exonic SNPs in potentially altering gene function and expression. By understanding these differences, researchers can gain insights into the functional consequences of genetic variations in different genomic regions.

# Conclusion

Throughout this tutorial, we explored a comprehensive workflow for analyzing DNA sequences and assessing the impact of mutations using a pre-trained DNA language model. We began by tokenizing DNA sequences, converting them into numerical representations that the model could process and analyze effectively. Following this, we loaded and configured a pre-trained model to handle these sequences, ensuring seamless integration with the tokenizer. We then delved into comparing the effects of mutations, both with and without amino acid modifications, showcasing the model's ability to discern subtle differences in sequence impacts. Furthermore, we focused on analyzing the effects of Single Nucleotide Polymorphisms (SNPs) in exons and introns. By loading the relevant sequences, computing the L2 distances between embeddings of wild-type and mutated sequences, and visualizing the results, we quantified and compared the impact of SNPs in these regions. Our analysis revealed that SNPs in exons generally have a more significant impact than those in introns, a finding supported by statistical tests. This tutorial underscores the power of pre-trained models in bioinformatics, offering valuable insights into the functional consequences of genetic variations and paving the way for further research in genomics.
