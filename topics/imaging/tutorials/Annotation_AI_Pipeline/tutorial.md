---
layout: tutorial_hands_on

title: AI pipeline for annotating marine species (Project Moorev - Marine)
zenodo_link: https://doi.org/10.5281/zenodo.20639741
questions:
- How can we use artificial intelligence to annotate images of marine species?
- How can we check and correct these annotations?
objectives:
- Use SAM3 to automatically detect marine species from a simple keyword
- Correct generated annotations with Edit COCO Annotation
time_estimation: 30m
key_points:
- AI does most of the annotation work, but a human eye is still needed to check and correct it.
- Once set up, this pipeline can be reused for other species or other projects.
tags:
  - Marine ecosystems
  - biodiversity
  - AI segmentation
  - YOLO
  - deep learning
  - object detection
contributions:
  authorship:
    - TuturBaba
    - yvanlebras
    - NadineLeBris
  testing: 
    - kpayet26
  funding:
    - moorev
    - sorbonneuniv
    - ISYEB
    - mnhn
    - pndb
lang: en
translations:
    - fr
---

In this tutorial, we will explore a complete pipeline for annotating images and videos of marine species. As an example, we will use a video from the [Moorev](https://moorev.fr/) project, cut into a short clip and available on Zenodo.

This pipeline has two main steps:
1. **Automatic annotation** using a text prompt with {% tool [SAM3 Semantic Segmentation](toolshed.g2.bx.psu.edu/repos/ecology/sam3_semantic_segmentation/sam3_semantic_segmentation/1.0.1+galaxy6) %}
2. **Correction and validation** of the annotations with {% tool [Edit COCO Annotation](toolshed.g2.bx.psu.edu/repos/bgruening/edit_coco_annotation/edit_coco_annotation/1.0.0+galaxy0) %}

> <details-title>Learn more about the MOOREV project</details-title>
>
> **MOOREV**: Microclimates and observation tools for the responses of marine life on the seafloor
>
> **Observing interactions between marine species facing microclimatic gradients on the shore**
>
> The MOOREV project is led by Nadine Le Bris (Sorbonne University, Concarneau Marine Station),
> with support from the National Museum of Natural History, the CNRS, and the Fondation de France.
> It started in 2022. Its goal is to better understand, and help others understand, the effects of
> climate disturbance on coastal biodiversity.
>
> The project uses underwater imaging methods to observe benthic species at the individual level,
> across different shore habitats. Groups of school students repeat data collection at their study
> sites, which are labelled by the Educational Marine Areas program of the French Biodiversity
> Office, across tide cycles, seasons, and years.
>
> The project brings together researchers and environmental education professionals. It is
> co-built with school classes and their teachers, using the shore as a natural laboratory. In the
> long term, it aims to support protection and conservation measures, taking into account the links
> between climate change and marine socio-ecosystems.
>
{: .details}

> <agenda-title>In this tutorial, we will cover:</agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Loading the data

The video used in this tutorial is a short clip taken directly from the Moorev project, available on Zenodo. This clip was chosen because it shows most of the features of the tools presented here.

> <hands-on-title>Load the video into Galaxy</hands-on-title>
>
> 1. Create a new history for this tutorial (for example: "Annotation Moorev")
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
> 2. Import the file from [Zenodo]({{ page.zenodo_link }}):
>
>    ```
>    https://zenodo.org/records/20639741/files/Lieujaune-PorzBreign-GOPR5167_Edited.mp4?download=1
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Check that the file is shown in **green** in the history before you continue
>
> 4. Check the datatype
>
>    > <tip-title>Check the datatype</tip-title>
>    >
>    > Galaxy assigns the datatype automatically. Since most tools use it to filter their inputs, it is important to make sure it is correct. In our case, check that the type is `mp4`. If it is not, the import may have had an error, and you can change the type manually.
>    >
>    > ![Check the datatype](../../images/Annotation_AI_Pipeline/Check_data.png){: style="width:75%; display:block; margin:auto;"}
>    >
>    > {% snippet faqs/galaxy/datasets_change_datatype.md datatype="mp4" %}
>    >
>    {: .tip}
{: .hands_on}

# Automatic annotation with SAM3

SAM3 is a text-guided segmentation model. It automatically detects and segments objects that match your prompt, with no need for previous annotations.

> <tip-title>A dedicated SAM3 tutorial exists</tip-title>
>
> If you want to learn more about SAM3 and its parameters, check out the dedicated tutorial:
> [SAM3 Tutorial]({% link topics/ecology/tutorials/SAM3/tutorial.md %})
>
{: .tip}

> <hands-on-title>Setting up SAM3</hands-on-title>
>
> 1. {% tool [SAM3 Semantic Segmentation](toolshed.g2.bx.psu.edu/repos/ecology/sam3_semantic_segmentation/sam3_semantic_segmentation/1.0.1+galaxy6) %} with these parameters:
>
>    - {% icon param-file %} *"Model data"*: `Segment Anything Model 3 (SAM 3)` (default)
>    - {% icon param-select %} *"Input type"*: `One video`
>        - {% icon param-file %} *"Input video file"*: `1: Lieujaune-PorzBreign-GOPR5167_Edited.mp4`
>        - {% icon param-select %} *"Video quality"*: `Original quality (lossless)` (default)
>        - {% icon param-select %} *"Additional output formats"*: `Empty` (default)
>        - {% icon param-select %} *"Export individual frames"*: `Both annoted and raw frames`
>    - {% icon param-text %} *"Text prompt"*: `fish`
>    - {% icon version %} *"Confidence threshold"*: `0.5`
>    - {% icon version %} *"Video frame stride"*: `5`
>    - {% icon param-toggle %} *"Show bounding boxes on annotated output"*: `Yes`
>    - {% icon param-toggle %}  *"Inference image size"*: `644 (default, balanced)` (default)
>    - {% icon param-toggle %} *"Normalize outputs?"*: `No` (default)
>    - {% icon param-toggle %} *"Normalize outputs?"*: `No` (default)
>
> 2. Click **Run Tool**
>
>    > <comment-title>Processing time</comment-title>
>    >
>    > Processing can take several minutes, depending on the length of the video and the server resources. Wait until the outputs turn green in the history.
>    >
>    {: .comment}
>
> 3. Once it is finished, you will have these items in your history:
>    - **5: Annotation COCO**: `annotations.json` — the segmentation masks in COCO format
>    - **4: Annotated Frames**: the collection of extracted images (annotated)
>    - **3: Raw Frames**: the collection of extracted images (not annotated)
>    - **2: COCO Extracted Frames**: the annotated video, for visual checking
>
> 4. Check the results visually by clicking {% icon galaxy-eye %} on **Annotated Outputs**
>
>    ![SAM3 segmentation result on the video](../../images/Annotation_AI_Pipeline/Visualize_output_SAM3-ezgif.com-video-to-gif-converter.gif){: style="width:75%; display:block; margin:auto;"}
>
>    > <comment-title>Limits of SAM3</comment-title>
>    >
>    > As you can see, not all annotations are perfect. False positives are the most common problem. This is why the next step, correction and validation, is essential.
>    >
>    {: .comment}
>
{: .hands_on}

# Correcting and validating the annotations

We will now see how to correct the annotations with **Edit COCO Annotation**, and then check the result with **COCO Annotation Visualizer**.

## Correcting with Edit COCO Annotation

The **Edit COCO Annotation** tool lets you modify COCO annotations without opening the images. It has three modes:

- **Keep**: keep only the listed IDs (remove all the others) — useful when only a few IDs need to be kept; it can also rename them
- **Remove**: remove only the listed IDs
- **Rename**: rename the listed IDs without removing any others

> <hands-on-title>Keep, remove, or rename annotations</hands-on-title>
>
> 1. {% tool [Edit COCO Annotation](toolshed.g2.bx.psu.edu/repos/bgruening/edit_coco_annotation/edit_coco_annotation/1.0.0+galaxy0) %} with these parameters:
>
>    - {% icon param-file %} *"COCO annotation file"*: `5: Annotation COCO `
>    - {% icon param-select %} *"Mode"* : `Keep - keep only listed IDs (remove all others)`
>        - {% icon param-repeat %} *"1: Track to keep"*
>            - {% icon param-text %} *"Track IDs"*: `0,4-6`
>            - {% icon param-text %} *"Rename"*: `Gobiusculus flavescens`
>            - {% icon param-text %} *"Frame min"*: `Empty` (default)
>            - {% icon param-text %} *"Frame max"*: `Empty` (default)
>        - {% icon param-repeat %} *"2: Track to keep"*
>            - {% icon param-text %} *"Track IDs"*: `2`
>            - {% icon param-text %} *"Rename"*: `Crenilabre`
>            - {% icon param-text %} *"Frame min"*: `Empty` (default)
>            - {% icon param-text %} *"Frame max"*: `Empty` (default)
>        - {% icon param-text %} *"Frame max"*: `Empty` (default)
>
>    > <comment-title>The same result in two steps</comment-title>
>    >
>    > You can also get the same result in two steps: use **Remove** to remove IDs 1 and 3, then use **Rename** to rename the remaining groups.
>    >
>    > Note: SAM3 can sometimes give the same Track ID to two different objects at different moments in the video. This is why the **Frame min** and **Frame max** parameters were added.
>    >
>    {: .comment}
>
> 2. Once it is finished, you will have this item in your history:
>    - **54: Edited COCO annotations**: the JSON file in COCO format, modified with the parameters set above
>
{: .hands_on}

### Visualizing the annotations with COCO Annotation Visualizer

This tool lets you easily visualize a JSON file in COCO format, shown on top of a video or images. Here we use it to check our annotations after correction.

> <hands-on-title>Visualize the COCO annotations</hands-on-title>
>
> 1. {% tool [COCO Annotation Visualizer](toolshed.g2.bx.psu.edu/repos/bgruening/coco_annotation_visualizer/coco_annotation_visualizer/1.0.0) %} with these parameters:
>
>      - {% icon param-file %} *"Input video file"*: `1: Lieujaune-PorzBreign-GOPR5167_Edited.mp4`
>      - {% icon version %} *"Frame stride"*: `5`
>    - {% icon param-file %} *"COCO annotation file"*: `54: Edited COCO annotations`
>    - {% icon param-text %} *"Filter categories"*: `Empty` (default)
>    - In *"Display options"*:
>        - {% icon param-toggle %} *"Show bounding boxes"*: `No`
>        - {% icon param-toggle %} *"Show segmentation masks"*: `Yes` (default)
>        - {% icon param-toggle %} *"Show category labels"*: `Yes` (default)
>        - {% icon param-toggle %} *"Show annotation count"*: `No` (default)
>        - {% icon version %} *"Mask opacity"*: `0.4` (default)
>        - {% icon version %} *"Bounding box thickness"*: `2` (default)
>        - {% icon version %} *"Label font scale"*: `0.6` (default)
>        - {% icon param-select %} *"Color mode"*: `Per instance (different color for each annotation)`
>        - {% icon param-select %} *"Output image format"*: `PNG (lossless)` (default)
>    - {% icon param-select %} *"Output mode"* : `Both (frames and video)`
>        - {% icon version %} *"Video frame rate (FPS)"*: `5.0`
>    - {% icon param-toggle %} *"Annotated frames only"*: `Yes` (default)
>
> 2. Once it is finished, you will have these items in your history:
>    - **56: Annoted video**: the annotated video
>    - **55: Annotated images**: each annotated frame
>
>    ![Annotated video with the COCO segmentation masks](../../images/Annotation_AI_Pipeline/Visualize_output_Coco_Visualizer.gif "56: Annoted video "){: style="width:65%; display:block; margin:auto;"}
{: .hands_on}

# Conclusion

You now know how to use this complete pipeline to annotate videos of marine species:

- **SAM3** automatically generates annotations from a simple text prompt
- **Edit COCO** lets you quickly clean up annotations by keeping, removing, or renaming track IDs
- **COCO Annotation Visualizer** lets you visually check the results at each step

This pipeline can be reused for other species or other marine imaging projects — you just need to adapt the prompt and the SAM3 confidence settings.
