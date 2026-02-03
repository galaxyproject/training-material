---
layout: faq-page
---

## What is the HAM10000 dataset?

The dataset used in this tutorial is a preprocessed, balanced subset of the original HAM10000 collection. It contains 1,400 dermoscopic images (200 per class) covering seven types of pigmented skin lesions: melanoma, melanocytic nevus, basal cell carcinoma, actinic keratosis, benign keratosis, dermatofibroma, and vascular lesions. The preprocessing followed the methodology from Shetty et al. (2022), including image selection, resizing to 96×96 pixels, and horizontal flip augmentation to create a balanced dataset suitable for deep learning.

## What is class imbalance and why does it matter?

Class imbalance occurs when the training dataset has a significantly different number of samples for different classes. In the original HAM10000 dataset, melanocytic nevi represented 67% of samples while dermatofibroma and vascular lesions represented less than 2% each. This imbalance can cause machine learning models to perform poorly on minority classes. To address this, our dataset was preprocessed to create a balanced set with 200 images per class (14.3% each), using horizontal flip augmentation on a selected subset of 100 images per class.

## Why use data augmentation?

Data augmentation artificially increases the size and diversity of the training dataset by applying transformations to existing images. In this tutorial, horizontal flip augmentation was applied during preprocessing (following Shetty et al. 2022) because:
- Lesion orientation is not clinically significant for diagnosis
- It doubled our dataset from 100 to 200 images per class
- It helps create a balanced dataset addressing class imbalance
- It improves model generalization without requiring additional data collection

## What is transfer learning?

Transfer learning is a technique where a pre-trained model (trained on a large dataset like ImageNet) is adapted for a new task. Instead of training from scratch, we use the learned features from the pre-trained model as a starting point, then fine-tune them for skin lesion classification. This approach:
- Requires fewer training samples
- Trains faster
- Often achieves better performance
- Is especially useful for medical imaging

## What does "fine-tune" mean?

Fine-tuning refers to adjusting all or most of the weights in a pre-trained model for a specific new task. In the Image Learner, setting "Fine Tune" to "True" means that all layers of the pre-trained model will be updated during training, not just the final classification layers. This typically leads to better performance on new tasks.

## What are the key performance metrics?

- **Accuracy**: Percentage of samples correctly classified. High for balanced datasets.
- **Precision**: Of the samples predicted as positive, what fraction actually are? Important when false positives are costly.
- **Recall**: Of the true positive samples, what fraction did the model find? Important when false negatives are costly.
- **F1-Score**: Harmonic mean of precision and recall, balancing both metrics.
- **ROC-AUC**: Area under the Receiver Operating Characteristic curve. Measures discrimination across all thresholds. Ranges from 0 to 1.

## What does ROC-AUC mean?

The ROC (Receiver Operating Characteristic) curve plots the True Positive Rate (sensitivity) against the False Positive Rate (1-specificity) at different classification thresholds. The AUC (Area Under the Curve) summarizes this in a single value between 0 and 1:
- **AUC = 1.0**: Perfect classification
- **AUC = 0.5**: Random guessing
- **AUC = 0.9880**: Excellent discrimination (our model)

## What is a confusion matrix?

A confusion matrix shows the breakdown of predictions versus actual values for each class. It displays:
- True Positives (correct positive predictions)
- True Negatives (correct negative predictions)
- False Positives (incorrect positive predictions)
- False Negatives (incorrect negative predictions)

For multi-class problems, it shows these values for each class pair, helping identify which classes are confused with each other.

## How do I interpret the confusion matrix?

- Look at the diagonal: high values mean good classification
- Look at off-diagonal values: high values indicate confusion between classes
- Compare rows to see if one class is consistently confused with another
- Consider clinical significance: some misclassifications may be more important than others

## Can I use this model on new images?

Yes! After training, the model can be saved and used to make predictions on new, unseen images. Galaxy allows you to download the trained model file and use it with the Image Learner prediction tool or export it for use in other platforms.

## What does the learning rate control?

The learning rate controls how much the model's weights change during training:
- **High learning rate**: Faster training but may overshoot optimal values and diverge
- **Low learning rate**: Slower training but more stable, less likely to diverge
- **Too low**: Training becomes very slow
- **Too high**: Training may fail to converge

A learning rate of 0.001 is conservative and works well for fine-tuning pre-trained models.

## What are epochs and batch size?

- **Epochs**: The number of times the model sees the entire training dataset. 30 epochs is typical for this task.
- **Batch size**: The number of images processed together before updating model weights. 32 is a common balance between memory usage and gradient stability.

## Why use a 70/10/20 train/validation/test split?

- **Training set (70%)**: Used to train the model
- **Validation set (10%)**: Used to tune hyperparameters and monitor for overfitting
- **Test set (20%)**: Used for final evaluation (never seen during training)

This split ensures fair evaluation of model performance on unseen data.

## How can I improve model performance?

Potential strategies include:
- Increasing the number of epochs (if not overfitting)
- Adjusting the learning rate
- Using different data augmentation strategies
- Trying different model architectures
- Collecting more training data
- Addressing class imbalance with different techniques
- Using class weighting to penalize misclassifications of minority classes

## What is overfitting?

Overfitting occurs when a model learns the training data too well, including its noise and peculiarities, and performs poorly on new data. Signs include:
- Training accuracy much higher than validation/test accuracy
- Continuously decreasing training loss but increasing validation loss

## Can I run this tutorial on Galaxy?

Yes! This tutorial is designed to be run on Galaxy. You can access Galaxy at:
- [Galaxy US Server](https://usegalaxy.org)
- [Cancer Galaxy](https://cancer.usegalaxy.org)



