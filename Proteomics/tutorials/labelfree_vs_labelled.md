---
layout: tutorial_hands_on
topic_name: Label-free versus Labelled - How to Choose Your Quantitation Method
tutorial_name: labelfree-vs-labelled
---


# Label-free versus Labelled - How to Choose Your Quantitation Method

## Purpose
A basic question when planning any quantitative proteomic experiment is the choice between label-free versus labelled quantitation. Here, we want to explain and discuss benefits and drawbacks of each method. 
Finally, we will give a guideline what points may be vital for the decision in our opinion. We will not discuss different labelling methods here in detail, but focus on our basic question.

A basic overview of different quantitation techniques is in the figure by Marc Vaudel given below:

![Overview quantitation techniques](../images/Vaudel_label_vs_labelfree.png)

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. [Drawbacks and Benefits of Labelled Quantitation](#pro-contra)
> 2. [Guideline: How to Choose Your Technique](#conclusion)

<a name="pro-contra"/></a>
## Drawbacks and Benefits of Labelled Quantitation
Label-free techniques are immensely popular in the proteomics community and may be the most straight-forward laboratory technique in the field. However, labelled approaches have some benefits, which make them the method-of-choice for some scientific questions. 

We summarized the pros and cons in the table below:

### Drawbacks and Benefits: Overview Table
category | label-free | labelled
:--|:--:|:--:
machine time | more | **less**
wet lab complexity & time | **little** | medium
comparability of samples | difficult | **easy**
data analysis | complex | complex
experimental design | **flexibel** | fixed

The **superior technique** in each line is marked in **bold**.

### Drawbacks and Benefits: Detailed Explanation
- **Machine time**: In label-free experiments, each sample is measured in a separate MS run. In labelled experiments, samples of each condition are combined prior to the MS run. This cuts down the machine time needed by the complexity of the labelling technique (usually between 2 and 8 times less machine time).
- **Wet lab complexity & time**: While label-free samples can be measured without much preparation, all labelling techniques need additional pretreatments in the wet lab. The samples have to be labelled either metabolically (e.g. by SILAC) or chemically (e.g. ICAT or iTRAQ) and the different conditions have to be combined. Thus, label-free techniques are less prone to wet lab errors than labelling techniques.
- **Comparability of samples**: A drawback of the label-free approaches is that exterior conditions (e.g. temperature, experimenter) may differ between samples. Such differences do not occur in labelled experiments, because all samples are measured in the very same MS run. Thus, label-free are more prone to errors introduced by the measurement conditions than labelled. 
Including a well-chosen standard in label-free experiments (e.g. a labelled control sample mixed to each sample prior to the MS run) may reduce this problem, but has to be carefully planned in beforehand.
- **Data analysis**: Data analysis of each type of experiment has it's special pitfalls. In our opinion, the benefits and drawbacks cancel each other out.
- **Experimental design**: Label-free approaches have the advantage of being very adaptable, even after having started the experiment. New samples may be included at any time. In contrast, labelled approaches need the same number (n) of each condition. Thus, new samples cannot be measured, if they cannot be matched to a control. 
Labelled techniques like the ["Super-SILAC"](https://www.ncbi.nlm.nih.gov/pubmed/20364148) approach do reduce this problem, but also need to be carefully planned in beforehand.

<a name="conclusion"/></a>
## Guideline: How to Choose Your Technique
Considering the points mentioned above, the following guideline may help you in finding the right method for your scientific question:

1. **Experimental design** & **Comparability of samples**: If you cannot predict the availability of samples in beforehand, flexibility is highly important. Therefore, we advise you to use a label-free approach for animal or clinical samples.
2. **Wet lab complexity** & **Data analysis**: Your experience matters! If you see only weak advantages of one approach, stick to the technique you are already used to. It will save time and reduce the likelihood of errors, both in the wet lab and *in silico*.
3. **Machine time**: If you have only limited access to the mass spectrometer, choose a labelling approach.
4. **Wet lab complexity & time**: If you have limited lab resources (e.g. time, man-power, chemicals), choose a label-free approach.