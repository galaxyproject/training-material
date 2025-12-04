---
layout: tutorial_hands_on
title: "Practical deep learning with PyTorch"
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
questions:
- to do
objectives:
- Understanding tensors
- Handling data using Dataset and DataLoader
- Learning the concepts of a loss function and an optimizer
- Discovering the basic training loop
- Learning the concepts of an activation function
- Training a simple fully-connected neural network
time_estimation: 1H
redirect_from: 
- /topics/statistics/tutorials/deep-learning-without-gai-with-python/tutorial.md
key_points:
- Tensors are like NumPy arrays but with GPU acceleration and automatic differentiation.
- PyTorch provides a Dataset class to handle data and a DataLoader class to load data in batches.
- The train loop consists of forward pass, loss calculation, backward pass, and optimization.
- Activation functions introduce non-linearity into the model.
contributions:
  authorship:
  - ralfg
tags:
- elixir
- ai-ml
priority: 4
notebook:
  language: python
  pyolite: true
---
> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Neural networks with PyTorch

PyTorch is compiled in different versions for different systems. When working locally, check out the PyTorch [Get Started Guide](https://pytorch.org/get-started/locally/) to install PyTorch with the appropriate CUDA version for your system.


```python
import torch

_ = torch.manual_seed(42)
```

## Working with tensors

In their most basic form, tensors are just multi-dimensional arrays.


```python
t = torch.tensor([1, 2, 3])
```

Tensors can be stored on different devices, such as CPU or GPU. PyTorch provides a simple way to move tensors between devices using the `.to()` method. This is useful for leveraging the computational power of GPUs for deep learning tasks.


```python
t.device
```


```python
t.to(device='cuda:0')
```

As in Numpy, tensors support many data types, including `float`, `int`, and `bool`. A full list of data types can be found in the [PyTorch documentation](https://pytorch.org/docs/stable/tensors.html#torch-tensor-dtypes). When working with floats, 32-bit floats are most common.

To save memory, 16-bit floats can also be used, in what is sometimes termed "lower precision" or "mixed precision". This is especially useful when training large models on GPUs, as it can significantly reduce memory usage and speed up training.

While less common in deep learning tasks, 64-bit floats can also be used for high-precision calculations.


```python
t.dtype
```


```python
t.to(dtype=torch.float32)
```

What really sets apart tensors from Numpy arrays is that they can automatically record every operation performed on them. In other words, they build a *computation graph*, which is a directed acyclic graph tracing how each tensor value was computed. This graph is what enables *automatic differentiation*, the core mechanism behind *back-propagation* in neural-network training.

Every time you combine or transform tensors, PyTorch notes that operation and links the inputs and outputs in its graph. Once you’ve computed some final output (for instance, a loss), you call `.backward()`, and PyTorch traverses the graph in reverse to compute gradients - the derivatives of that output with respect to each tensor marked for tracking. Those gradients tell you how to adjust your model’s parameters to reduce the loss. Without a recorded graph, you’d have to derive and implement each derivative by hand.

Back propagation will be covered later in this tutorial. For now, it is important to know that:

1. Tensors can keep track of all operations that lead to the resulting values.
2. This graph is used to compute the gradients of the loss function with respect to the model parameters.
3. Gradients are essentially a combination of derivatives, which are used to update the model parameters during training.

We can enable gradient tracking by setting the `requires_grad` attribute to `True`. Lets perform a simple operation:


```python
t = torch.tensor([1.0, 2.0, 3.0], requires_grad=True)
out = t.pow(2).sum()
out
```

The grad_fn attribute of the tensor shows the function that created it. In this case, it is a `SumBackward` function, which indicates that the tensor was created by summing two other tensors.

Calling `.backward()` on the tensor will compute the gradients of all tensors that have `requires_grad=True` and are part of the computation graph leading to this tensor. The gradients will be stored in the `.grad` attribute of those tensors.


```python
out.backward()
t.grad
```

Indeed, calling `.backward()` on the `out` tensor computes the gradients of `out` with respect to `t`, and stores them in `t.grad`. The gradient (or derivatives) of the sum of squares with respect to each element is `2 * t[i]`, so for `t = [1.0, 2.0, 3.0]`, the gradients will be `[2.0, 4.0, 6.0]`.

> <hands-on-title></hands-on-title>
> 1. Create a tensor with `requires_grad=True`
> 2. Perform a sequence of operations on it. 
> 3. Call `backward()`
> 4. Check the gradient.
{: .hands-on}


```python
# Add your own code here
```

More information on tensors can be found in the [PyTorch documentation](https://docs.pytorch.org/docs/stable/tensors.html#initializing-and-basic-operations), including a full list of tensor operations. You can also check the [PyTorch tutorial on tensor basics](https://pytorch.org/tutorials/beginner/basics/tensorqs_tutorial.html) for more examples and explanations.

## Handling data sets in PyTorch


### Exploring the breast cancer dataset

The breast cancer dataset is a classic dataset for binary classification tasks. It contains features extracted from images of breast cancer tumors, along with labels indicating whether the tumor is malignant or benign. As the features have already been extracted and processed into a numerical format, the dataset is ideal for a simple neural network classification task.

The dataset is available in the [`sklearn.datasets`](https://scikit-learn.org/stable/modules/generated/sklearn.datasets.load_breast_cancer.html) module, and can be loaded using the [`load_breast_cancer()`](https://scikit-learn.org/stable/modules/generated/sklearn.datasets.load_breast_cancer.html#sklearn.datasets.load_breast_cancer) function.


```python
import pandas as pd
import seaborn as sns

from sklearn.datasets import load_breast_cancer
```


```python
cancer_data = load_breast_cancer(as_frame=True)
features = cancer_data['data']
targets = cancer_data['target']
target_names = cancer_data['target_names']
feature_names = cancer_data['feature_names']
```

Let's check out the target classes:


```python
for i, name in enumerate(target_names):
    print(i, name)
```

And lets see how many data points we have for each class:


```python
targets.value_counts()
```

*Question: Is the dataset balanced? What does that imply for the training task?*

Let's check out the features:


```python
features.head()
```

And let's visualize some of them:


```python
data = pd.concat([features, targets.astype(bool).map({False: "malignant", True: "benign"})], axis=1)

sns.pairplot(
    data[["worst radius", "worst texture", "worst symmetry", "worst fractal dimension", "target"]],
    hue="target",
    diag_kind='kde',
    plot_kws={"alpha": 0.2, "s": 10},
)
```

We can see that these features show some separation between the two classes, but not all of them are equally useful. For example, "worst fractal dimension" is not as informative as "worst radius". The neural network will allow us to combine these features in a way that maximizes the separation between the two classes.

### Creating a data set object

Data handling for deep learning is a bit different from regular machine learning. In deep learning, we often work with large datasets that are too big to fit into memory all at once. Instead, we load the data in batches, which allows us to train both save memory and speed up the training process. The latter is possible applying gradient descent per batch instead of for the entire dataset at once.

PyTorch provides to convenient helper classes for handling data: [`Dataset`](https://pytorch.org/docs/stable/data.html#torch.utils.data.Dataset) and [`DataLoader`](https://pytorch.org/docs/stable/data.html#torch.utils.data.DataLoader). The former is an abstract class that represents a single dataset, while the latter class is used to load data from a dataset in batches.


```python
from torch.utils.data import Dataset
```

As we already loaded the data into a Pandas DataFrame, the initialization of the dataset is straightforward. We just need to convert the DataFrame into a PyTorch tensor with the correct data types. The `__getitem__` method is used to retrieve a single sample from the dataset, while the `__len__` method returns the total number of samples in the dataset.


```python
class CancerDataset(Dataset):
    def __init__(self, features: pd.DataFrame, targets: pd.Series):
        self.features = torch.tensor(features.values, dtype=torch.float32)
        self.targets = torch.tensor(targets.values, dtype=torch.float32)

    def __len__(self):
        return len(self.features)

    def __getitem__(self, idx):
        x = self.features[idx]
        y = self.targets[idx]
        return x, y

cancer_dataset = CancerDataset(features, targets)
```

Now the dataset can be accessed by index, just like a list. Each item is a tuple containing the features and targets as tensors.

> <hands-on-title></hands-on-title>
> Try to access the first 10 samples of the dataset.
> 
> > <question-title></question-title>
> > What do you see?
> {: .question}
{: .hands-on}


```python
# Add your own code here
```

Learn more about Datasets on the [PyTorch documentation](https://pytorch.org/docs/stable/data.html).

### Creating a data loader

PyTorch provides a convenient way to load data in batches using the [`DataLoader`](https://pytorch.org/docs/stable/data.html#torch.utils.data.DataLoader) class. It takes a `Dataset` object and provides an iterable over the dataset, allowing you to easily iterate over the data in batches. Training data can be immediately shuffled by the `DataLoader`, which is useful for training models. The `DataLoader` class also provides options for parallel data loading, which can speed up the training process if data parsing is computationally expensive.


```python
from torch.utils.data import DataLoader
cancer_dataloader = DataLoader(cancer_dataset, batch_size=32, shuffle=True)
```

To test the dataloader, we can iterate over the data. Each iteration yields a batch of data, which is a tuple containing the features and labels.


```python
for batch_x, batch_y in cancer_dataloader:
    print(batch_x.shape, batch_y.shape)
    break
```

Learn more about `DataLoader` on the [PyTorch documentation](https://pytorch.org/docs/stable/data.html#torch.utils.data.DataLoader) and check out the [PyTorch tutorial on data loading and processing](https://pytorch.org/tutorials/beginner/data_loading_tutorial.html).

*Question: Go to the documentation and look up the difference between an iterable and map-style dataset. When would you use one over the other? What are the advantages and disadvantages of each?*

## Building a simple linear classifier

### Creating a Pytorch model object to define the model architecture

The [`nn.Module`](https://pytorch.org/docs/stable/generated/torch.nn.Module.html#torch.nn.Module) class is the base class for all neural network modules in PyTorch. It provides a convenient way to define and manage the architecture and all parameters of your model, as well as to define the forward pass of your model.
The `__init__` method is where you define the layers of your model. In this case, we are using a single linear layer with 30 input features and 2 output features. The `forward` method defines the forward pass of your model, which takes in the input data and passes it through the layers defined in the `__init__` method.

A `backward` method is not needed, as PyTorch automatically computes the gradients for you when you call `.backward()` on the output of your model. The backward steps are inferred using the computation graph.


```python
from torch import nn

class BreastCancerClassifier(nn.Module):
    def __init__(self):
        super().__init__()
        self.layer1 = nn.Linear(30, 1)

    def forward(self, x):
        x = self.layer1(x)
        return x
```

*Question: What are the numbers 30 and 1 in the [`nn.Linear`](https://pytorch.org/docs/stable/generated/torch.nn.Linear.html#torch.nn.Linear) constructor? Why were they chosen?*


```python
model = BreastCancerClassifier()
model.forward(batch_x).shape
```

### The loss function

The loss function is simply a measure of how well the model is performing. PyTorch provides a variety of loss functions in the [`torch.nn`](https://pytorch.org/docs/stable/nn.html#loss-functions) module, including mean squared error (MSE), cross-entropy loss, and binary cross-entropy loss. The choice of loss function depends on the type of problem you are trying to solve. For example, if you are working on a regression problem, you would typically use MSE, while for a classification problem, you would use cross-entropy loss.

In this case, we are using binary cross-entropy loss. Its values are between 0 and 1, where 0 means the model is perfect and 1 means the model is completely wrong. The closer the value is to 0, the better the model is performing. The following blog post provides a great visual explanation of the binary cross-entropy loss function: [Understanding binary cross-entropy / log loss: a visual explanation | Towards Data Science](https://towardsdatascience.com/understanding-binary-cross-entropy-log-loss-a-visual-explanation-a3ac6025181a).


```python
loss_function = nn.BCEWithLogitsLoss()
```

We can also calculate the ROC-AUC score on the entire validation set using the [`roc_auc_score`](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.roc_auc_score.html#sklearn.metrics.roc_auc_score) function from `sklearn.metrics`.


```python
from sklearn.metrics import roc_auc_score
```

### The optimizer

The optimizer is the algorithm that adjusts the parameters of a model to minimize the loss function. It takes the gradients computed during the backward pass and
determines how to use them to update the model parameters. The most important hyperparameter for the optimizer is the learning rate, which controls how much to adjust the parameters at each step. A smaller learning rate means smaller updates, while a larger learning rate means larger updates.

In PyTorch, optimizers are implemented as classes in the [`torch.optim`](https://pytorch.org/docs/stable/optim.html#algorithms) module. The most common optimizer is *stochastic gradient descent (SGD)*, but there are many other optimizers available, such as Adam, RMSprop, and Adagrad. Each has its own advantages and disadvantages, and the choice of optimizer can have a significant impact on the performance of your model.

When setting up the optimizer, you need to pass in the parameters of your model that you want to optimize. In this case, we are passing in the parameters of the `model` object. The optimizer will then update these parameters during training.


```python
optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)
```

### The training loop
The training loop is the core of the training process. It goes through the following steps:

1. Get a batch of data from the `DataLoader`.
2. Zero the gradients of the model parameters using `optimizer.zero_grad()`. Otherwise, the gradients would accumulate over multiple iterations, which is not what we want.
3. Perform a forward pass through the model to get its current predictions.
4. Compute the loss.
5. Perform a backward pass to compute the gradients.
6. Update the model parameters using the optimizer.

Here, we implement a single training step:


```python
for batch_x, batch_y in cancer_dataloader:                  # 1. Get a batch of data
    optimizer.zero_grad()                                   # 2. Reset gradients
    y_pred = model(batch_x)                                 # 3. Forward pass
    loss_value = loss_function(y_pred.squeeze(), batch_y)   # 4. Compute loss
    loss_value.backward()                                   # 5. Backward pass, based on the loss
    optimizer.step()                                        # 6. Update weights

    break  # For demonstration, we only do one batch

print(f"{loss_value.item():.4f}")
```

*Note: [`squeeze()`](https://pytorch.org/docs/stable/generated/torch.squeeze.html#torch-squeeze) is used to remove the extra dimension from the output of the model. The model outputs a tensor of shape `(batch_size, 1)`, but we want it to be of shape `(batch_size,)` for the loss function.*

### Training a model

Of course, we want to train the model for multiple *epochs*. An epoch is a complete pass through the training data. In each epoch, we will go through all the batches of data in the `DataLoader` and perform the training steps described above. Additionally, we want to see how the model is performing on the *validation set*. To do this, we will perform a single forward step with the current model on a validation dataset and compute the loss after each epoch.

*Question: What is the difference between a validation set and a test set? Why do we need both? When do you use one over the other?*

To split our data into a training and validation set, we can use the [`random_split`](https://pytorch.org/docs/stable/data.html#torch.utils.data.random_split) function.


```python
data_train, data_validation = torch.utils.data.random_split(cancer_dataset, [0.8, 0.2])

dataloader_train = DataLoader(data_train, batch_size=32, shuffle=True)
dataloader_validation = DataLoader(data_validation, batch_size=32, shuffle=False)
```

Complete the code below to train the model. Add the train step and validation step to the training loop. Take into account that not all points from the training step are needed in the validation step.


```python
model = BreastCancerClassifier()
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
loss_function = nn.BCEWithLogitsLoss()

# We will train the model for 100 epochs
n_epochs = 10
train_losses = []
validation_losses = []

for epoch in range(1, n_epochs + 1):
    # Train the model
    model.train()
    epoch_loss_train = 0  # Accumulate the loss for this epoch
    for batch_x, batch_y in dataloader_train:
        # Add your code here for a single training step. Make sure the loss is placed in
        # the loss_value variable.
        # ---


        # ---
        epoch_loss_train += loss_value.item()

    avg_epoch_loss_train = epoch_loss_train / len(dataloader_train)  # Average loss for the epoch
    train_losses.append(avg_epoch_loss_train)

    # Evaluate the model
    model.eval()
    epoch_loss_validation = 0
    with torch.no_grad():
        for batch_x, batch_y in dataloader_validation:
            # Add your code here for a single validation step. Make sure the loss is placed in
            # the loss_value variable.
            # ---


            # ---
            epoch_loss_validation += loss_value.item()

    avg_epoch_loss_validation = epoch_loss_validation / len(dataloader_validation)
    validation_losses.append(avg_epoch_loss_validation)

    print(f"Epoch {epoch}: Train Loss: {avg_epoch_loss_train:.4f}, Validation Loss: {avg_epoch_loss_validation:.4f}")

# Calculate the AUC score on the full validation set
x = data_validation[:][0]
y = data_validation[:][1]
y_pred = model(x)
print(f"ROC-AUC score: {roc_auc_score(y.numpy(), y_pred.detach().numpy()):.4f}")

# Plot the learning curves
pd.DataFrame({
    "Train Loss": train_losses,
    "Validation Loss": validation_losses,
}).plot(title="Losses", xlabel="Epoch", ylabel="Loss")
```

If everything is working correctly, you should see the train loss decreasing over time. The model should also be able to predict the labels of the validation set with a some accuracy.

## Building a linear deep neural network

So far, we have built a simple linear classifier. However, in practice, we often want to build more complex models that can learn more intricate patterns in the data. This means that we will add multiple layers to our model, making it a *deep neural network* (or in this case, rather a *shallow* neural network).

### Activation functions

For the breast cancer dataset, we will build a model with three hidden layers, each with a decreasing number of neurons. In between each layer, we will add an activation function. The activation function is a non-linear function that is applied to the output of each neuron in the layer. It introduces non-linearity into the model, allowing it to learn more complex patterns in the data. The most common activation functions are **ReLU (Rectified Linear Unit)**, **sigmoid**, and **tanh**. ReLU is the most commonly used activation function in deep learning, as it is computationally efficient and helps to mitigate the vanishing gradient problem. The sigmoid function is often used in the output layer for binary classification problems, as it outputs a value between 0 and 1, which can be interpreted as a probability. The tanh function is similar to the sigmoid function, but it outputs values between -1 and 1.

PyTorch provides activations functions as part of the [`torch.nn`](https://pytorch.org/docs/stable/nn.html#torch-nn) module. The most common activation functions are implemented as classes, and can be used as layers in the model. For example, to use the ReLU activation function, you can simply add [`nn.ReLU()`](https://pytorch.org/docs/stable/generated/torch.nn.ReLU.html#torch.nn.ReLU) to your model.

### Extending the model architecture

> <hands-on-title></hands-on-title>
> 1. Modify the `BreastCancerClassifier` class to include three hidden layers, each with a decreasing number of neurons (30, 15, 1). 
> 2. Use the ReLU activation function between each layer. 
>
> The final output layer should have 1 neuron and use the sigmoid activation function.
{: .hands-on}


```python
class BreastCancerClassifier(nn.Module):
    def __init__(self):
        super().__init__()
        # Add your code here
        # ---


        # ---

    def forward(self, x):
        # Add your code here
        # ---


        # ---
        return x
```

### Training the model

> <hands-on-title></hands-on-title>
> Repeat the training process by reusing your code from 1.3.3. 
>
> > <question-title></question-title>
> > 1. How does the model perform?
> > 2. Is it better than the previous model?
> {: .question}
{: .hands-on}


```python
# Add your own code here

```

# Working faster and cleaner with PyTorch Lightning

## Introduction

In the previous tutorial, we learned how to build a simple neural network using PyTorch. In this tutorial, we will learn how to use PyTorch Lightning to make our code cleaner and more efficient.

PyTorch Lightning is a lightweight wrapper around PyTorch that helps you organize your code and decouple the science code from the engineering code. It provides a high-level interface for training and testing models, making it easier to write clean and maintainable code. It also helps you to scale your code to multiple GPUs and TPUs, and to run it on different platforms (e.g., cloud, local, etc.) without changing your code.

We will also briefly look at Weights & Biases (wandb), a tool for tracking experiments and visualizing results. It is not required to use PyTorch Lightning, but it is a great tool to have in your toolbox. We will use it to log our training and validation metrics, and to visualize our results.

## PyTorch Lightning

### Repeating the basics

First, lets get some code from the previous tutorial. We will use the same dataset and model, but we will refactor the code to use PyTorch Lightning.


```python
import torch
import torch.nn as nn
import pandas as pd
from sklearn.datasets import load_breast_cancer
from torch.utils.data import Dataset, DataLoader

_ = torch.manual_seed(42)

class CancerDataset(Dataset):
    def __init__(self, features: pd.DataFrame, targets: pd.Series):
        self.features = torch.tensor(features.values, dtype=torch.float32)
        self.targets = torch.tensor(targets.values, dtype=torch.float32)

    def __len__(self):
        return len(self.features)

    def __getitem__(self, idx):
        x = self.features[idx]
        y = self.targets[idx]
        return x, y


class BreastCancerClassifier(nn.Module):
    def __init__(self):
        super().__init__()
        self.layer1 = nn.Linear(30, 20)
        self.layer2 = nn.Linear(20, 10)
        self.layer3 = nn.Linear(10, 1)
        self.activation = nn.ReLU()

    def forward(self, x):
        x = self.layer1(x)
        x = self.activation(x)
        x = self.layer2(x)
        x = self.activation(x)
        x = self.layer3(x)
        return x


cancer_data = load_breast_cancer(as_frame=True)
cancer_dataset = CancerDataset(cancer_data['data'], cancer_data['target'])

data_train, data_validation = torch.utils.data.random_split(cancer_dataset, [0.8, 0.2])
dataloader_train = DataLoader(data_train, batch_size=32, shuffle=True)
dataloader_validation = DataLoader(data_validation, batch_size=32, shuffle=False)
```

### Setting up the Lightning module

Take a look again at the training loop code from the previous tutorial. It contained quite some boilerplate code that would be highly similar across different ML projects. We will now refactor this code to use PyTorch Lightning, which will help us to remove a lot of this boilerplate code.

To do this, we will create a new class that inherits from [`LightningModule`](https://pytorch-lightning.readthedocs.io/en/stable/common/lightning_module.html#lightningmodule). This class will contain all the code for training and validating our model. Similar to the `nn.Module` class, we will need to implement the `__init__` method to initialize our model and the `forward` method to define the forward pass. However, we will also need to implement some additional methods that are specific to PyTorch Lightning:

- `training_step`: This method will be called for each batch of training data. It will contain the code for the forward pass and the loss calculation.
- `validation_step`: This method will be called for each batch of validation data. It will contain the code for the forward pass and the loss calculation.
- `configure_optimizers`: This method will be called to configure the optimizer and the learning rate scheduler.



```python
import lightning as L

class BreastCancerModule(L.LightningModule):
    def __init__(self, learning_rate=0.001):
        super().__init__()
        self.save_hyperparameters()
        self.model = BreastCancerClassifier()
        self.loss_fn = nn.BCEWithLogitsLoss()
        self.learning_rate = learning_rate
        self.optimizer = torch.optim.Adam

    def forward(self, x):
        return self.model(x).squeeze()

    def training_step(self, batch, batch_idx):
        x, y = batch
        y_hat = self(x)
        loss = self.loss_fn(y_hat, y)
        self.log('train_loss', loss, prog_bar=True)
        return loss

    def validation_step(self, batch, batch_idx):
        x, y = batch
        y_hat = self(x)
        loss = self.loss_fn(y_hat, y)
        self.log('val_loss', loss, prog_bar=True)

        return loss

    def configure_optimizers(self):
        return self.optimizer(self.parameters(), lr=self.learning_rate)
```

Note that the resulting code is much cleaner and more organized than the training loop code we had before.

### 2.2.3 Setting up the Lightning trainer

All configuration of the training process is done in the [`Trainer`](https://pytorch-lightning.readthedocs.io/en/stable/common/trainer.html#trainer) class. This class will take care of the training and validation loops, as well as logging and checkpointing. We will create an instance of this class and pass it our model and the training and validation data loaders.

> <hands-on-title></hands-on-title>
> Browse through the documentation for the trainer class and try to understand the different arguments of the class. Pay special attention to the `accelerator` and `devices` arguments, which are some of the most useful features of PyTorch Lightning.
{: .hands-on}


```python
trainer = L.Trainer(
    max_epochs=10,
    accelerator="auto",  # Automatically selects GPU or CPU
)
```

Note that the trainer immediately tells us whether we are using a GPU or not. This is a great feature of PyTorch Lightning, as it allows us to write code that is agnostic to the hardware we are using. We can run the same code on a CPU, a single GPU, or multiple GPUs without changing anything in our code.

### Fitting the model

Now we can simply call the `fit` method of the trainer to start training our model. It takes the model, the training data loader, and the validation data loader as arguments.


```python
trainer.fit(
    BreastCancerModule(learning_rate=0.001),
    train_dataloaders=dataloader_train,
    val_dataloaders=dataloader_validation,
)
```

While we had to take care of logging the progress ourselves in the previous tutorial, PyTorch Lightning takes care of this for us.

## Logging with Weights & Biases

Weights & Biases (wandb) allows us to log and visualize training and validation metrics across different training runs. It also has a feature for hyperparameter tuning, called Sweeps, which allows us to automatically search for the best hyperparameters for our model.

To use wandb, first go to the [wandb website](https://wandb.ai/) and create an account. You can easily sign in with an existing GitHub, Google, or Microsoft account. After signing in, you will be taken to the dashboard, where you can create a new project. You can also create a new API key, which you will need to use wandb in your code.




```python
import wandb
wandb.init(project="breast-cancer-classification")
```

Now we can add the wandb logger to our PyTorch Lightning trainer with the `logger`argument:


```python
trainer = L.Trainer(
    max_epochs=10,
    accelerator="auto",
    logger=L.pytorch.loggers.WandbLogger(
        project="breast-cancer-classification",
        log_model=True,
    ))

trainer.fit(
    BreastCancerModule(learning_rate=0.001),
    train_dataloaders=dataloader_train,
    val_dataloaders=dataloader_validation,
)

wandb.finish()  # Finish the wandb run
```

Go to the run URL as logged by wandb and check out the results. You should see a new run with the same name as the one you used in the code. You can click on it to see the details of the run, including the training and validation metrics, the model checkpoints, and the hyperparameters used for the run.

## Performing a hyperparameter sweep

Now that we have set up wandb, we can use it to perform a hyperparameter sweep. This will allow us to automatically search for the best hyperparameters for our model, such as number of layers, number of neurons per layer, learning rate, etc.

To do this, we must first update the module to accept hyperparameters as arguments. To make the number of layers and the number of neurons per layer configurable, we will implement a loop in the `__init__` method that creates the layers based on the hyperparameters.


```python
class BreastCancerClassifier(nn.Module):
    def __init__(self, hidden_layers=1, hidden_neurons=10):
        super().__init__()
        self.layers = nn.ModuleList()
        self.layers.append(nn.Linear(30, hidden_neurons))
        for _ in range(hidden_layers - 1):
            self.layers.append(nn.Linear(hidden_neurons, hidden_neurons))
        self.layers.append(nn.Linear(hidden_neurons, 1))
        self.activation = nn.ReLU()

    def forward(self, x):
        # Iterate over all layers except the last one and add activations
        for layer in self.layers[:-1]:
            x = layer(x)
            x = self.activation(x)
        # Last layer without activation
        x = self.layers[-1](x)
        return x
```

We must also slightly modify the Lightning module to accept the hyperparameters as arguments and pass them on to the model and the optimizer.


```python
class BreastCancerModule(L.LightningModule):
    def __init__(self, learning_rate=0.001, hidden_layers=1, hidden_neurons=10):
        super().__init__()
        self.save_hyperparameters()
        self.model = BreastCancerClassifier(hidden_layers, hidden_neurons)
        self.loss_fn = nn.BCEWithLogitsLoss()
        self.learning_rate = learning_rate
        self.optimizer = torch.optim.Adam

    def forward(self, x):
        return self.model(x).squeeze()

    def training_step(self, batch, batch_idx):
        x, y = batch
        y_hat = self(x)
        loss = self.loss_fn(y_hat, y)
        self.log('train_loss', loss, prog_bar=True)
        return loss

    def validation_step(self, batch, batch_idx):
        x, y = batch
        y_hat = self(x)
        loss = self.loss_fn(y_hat, y)
        self.log('val_loss', loss, prog_bar=True)

        return loss

    def configure_optimizers(self):
        return self.optimizer(self.parameters(), lr=self.learning_rate)
```

The following function will take the hyperparameters as arguments and setup the training run:


```python
def sweep(*args, **kwargs):
    # Create the trainer
    trainer = L.Trainer(
        max_epochs=10,
        accelerator="auto",
        logger=L.pytorch.loggers.WandbLogger(
            project="breast-cancer-classification",
            log_model=True,
        ))

    trainer.fit(
        BreastCancerModule(*args, **kwargs),
        train_dataloaders=dataloader_train,
        val_dataloaders=dataloader_validation,
    )
    wandb.finish()
```

Now, we can create a sweep configuration. This configuration will define the hyperparameters we want to search over and the values we want to try. We will use the `wandb` library to create a sweep configuration. To configure sweeps, check out the [wandb documentation](https://docs.wandb.ai/guides/sweeps).


```python
sweep_configuration = {
    "method": "random",
    "metric": {"goal": "minimize", "name": "val_loss"},
    "parameters": {
        "hidden_layers": {
            "values": [0, 1, 2]
        },
        "hidden_neurons": {
            "values": [2, 4, 8, 16, 32]
        },
        "learning_rate": {
            "min": 0.0001,
            "max": 0.01
        },
    },
}
```

Here, we initialize the sweep with its configuration:


```python
sweep_id = wandb.sweep(sweep=sweep_configuration, project="breast-cancer-classification")
```

With the `wandb.agent` function, we can start the sweep. Note that the `count` argument defines how many runs we want to perform. If it would not be set, the sweep would run indefinitely and keep testing different hyperparameters until we stop it manually.


```python
wandb.agent(sweep_id, function=sweep, count=20)
```

> <hands-on-title></hands-on-title>
>
> Go to the wandb dashboard and check out the results of the sweep. 
> 
> You should see a new sweep with the same name as the one you used in the code. You can click on it to see the details of the sweep, including the training and validation metrics, the model checkpoints, and the hyperparameters used for each run.
{: .hands-on}
