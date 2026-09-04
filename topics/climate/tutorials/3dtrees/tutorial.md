---
layout: tutorial_hands_on

title: '3Dtrees: Point cloud processing'
zenodo_link: ''
questions:
- What are LAS data ?
- How can I process and visualise point cloud data ? 
objectives:
- Learn to do point cloud processing
- 
time_estimation: 3H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
  autorship:
    - kgerb
    - Marie59

---

TO DO
=> Add description general of what does the workflow and how it can be usefull (small intro of 3Dtrees here (in a more details box))

=> Introduce the Data used and for what (small intro for data terra in a more details box)


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Title for your first section

Intro on the forest and its structure ? of IGN data. 

## Get data

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Recherche Data Gouv](https://entrepot.recherche.data.gouv.fr/) more specifically the Data register underthis [DOI](https://doi.org/10.15454/DMYWPB)
> 3. On the left pannel go on **Upload** select **{% icon pref-cloud %} Choose Remote Files**
> 
> ![Image of how the upload pannel looks like](../images/3dtrees/upload_pannel.png)
>
> 4. In the search bar write down "recherche", you should see the repo named **Recherche Data Gouv** appear, click on it.
> 5. In the search bar write down "pic", you should see the repo named **Données lidar acquises par drone au pied du Pic St Loup (34)** appear, click on it.
>
> ![Image of how the remote repo to select](../images/3dtrees/rdg_repo.png)
>
> 6. You should see 4 datasets appear, select them all and then in the bottom right click on **start**
>
> ![Image of how the data to upload](../images/3dtrees/data_upload.png)
>
> 7. The 4 datasets should appear into your history (in yellow when uploading) in green
>
> ![Image of how the data to uploading in the history](../images/3dtrees/data_inhistory.png)
>
> > <tip-title>las format</tip-title>
> > Check that the datatype is **las format**
> >
> >  {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
> {: .tip}
> 
> 8. Build a collection 
>
> {% snippet faqs/galaxy/collections_autobuild_list.md %}
{: .hands_on}


# Standardization

***TODO***: *Consider adding a detail box to expand the theory and add some info here*

> <details-title> More details about the theory </details-title>
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}


## 3DTrees: LAS/LAZ Standardization

> <hands-on-title> Standardize individual files in collection </hands-on-title>
>
> 1. {% tool [3DTrees: LAS/LAZ Standardization](toolshed.g2.bx.psu.edu/repos/bgruening/3dtrees_standardization/3dtrees_standardization/1.1.0+galaxy0) %} with the following parameters:
>    - *"Mode"*: `Single File Standardization`
>        - {% icon param-collection %} *"Point Cloud File"*: `output` (Input dataset collection)
>
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

In the same time you can do this 2nd step even if 1st one has not finished yet.

> <hands-on-title> 3DTrees Collection Check  </hands-on-title>
>
> Collection-level validation and consistency check using standardization tool in collection mode
>
> 1. {% tool [3DTrees: LAS/LAZ Standardization](toolshed.g2.bx.psu.edu/repos/bgruening/3dtrees_standardization/3dtrees_standardization/1.1.0+galaxy0) %} with the following parameters:
>    - *"Mode"*: `Collection Validation`
>        - {% icon param-collection %} *"Point Cloud Collection"*: `output` (Input dataset collection)
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

You should have 4 new objects in your history running yellow soon turning green.

![Standardization running](../images/3dtrees/standard_run.png)

# Thumbnail generation 

> <hands-on-title> 3D Trees Overview Generator </hands-on-title>
>
> Generate overview images from standardized point clouds
>
> 1. {% tool [3Dtrees: Overviews](toolshed.g2.bx.psu.edu/repos/bgruening/3dtrees_overviews/3dtrees_overviews/1.2.0+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Point Cloud Collection"*: `pc_standardized` (output of **3DTrees: LAS/LAZ Standardization** {% icon tool %})
>
> You should see to types of section images in your history (4 images in total have been generated)
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

![Example of the overview image on the section ns](../images/3dtrees/sectionns.png)

# Instance Segmentation by SegmentAnyTree

> <hands-on-title> 3DTrees: SmartTile </hands-on-title>
>
> 1. {% tool [3DTrees: SmartTile](toolshed.g2.bx.psu.edu/repos/bgruening/3dtrees_smart_tile/3dtrees_smart_tile/2.0.1+galaxy0) %} with the following parameters:
>    - *"Operation"*: `Tile (COPC normalize + tile + subsample)`
>        - {% icon param-file %} *"Input point clouds (LAZ/LAS)"*: `pc_standardized` (output of **3DTrees: LAS/LAZ Standardization** {% icon tool %})
>
> 5 new objects should appear in your history.
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

> <hands-on-title> 3Dtrees: SegmentAnyTree </hands-on-title>
>
> Deep learning tree segmentation
> 
> 1. {% tool [3Dtrees: SegmentAnyTree](toolshed.g2.bx.psu.edu/repos/bgruening/3dtrees_segmentanytree/3dtrees_segmentanytree/1.2.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input Dataset"*: `output_subsampled_res2` (output of **3DTrees: SmartTile** {% icon tool %})
>    - *"Predicted instance dimension name"*: `PredInstance`
>    - *"Predicted semantic dimension name"*: `PredSemantic`
>
> Your tool should look like that before running. It can take a few minutes to run.
> ![Image of the inputs within this tool](../images/3dtrees/smartile.png)
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}


> <hands-on-title> 3DTrees: SmartTile </hands-on-title>
>
> 1. {% tool [3DTrees: SmartTile](toolshed.g2.bx.psu.edu/repos/bgruening/3dtrees_smart_tile/3dtrees_smart_tile/2.0.1+galaxy0) %} with the following parameters:
>    - *"Operation"*: `Filter (remove overlap-only instances)`
>        - {% icon param-file %} *"Segmented tiles (LAZ/LAS)"*: `output` (output of **3Dtrees: SegmentAnyTree** {% icon tool %})
>        - {% icon param-file %} *"Tile layout JSON (tile_bounds_tindex.json)"*: `output_tile_bounds_json` (output of **3DTrees: SmartTile** {% icon tool %})
>        - *"Border zone width"*: `Derive from tile layout JSON (recommended)`
>        - *"Small cluster reassignment"*: `Enabled`
>        - *"After filtering"*: `Also remap filtered dimensions to another collection`
>            - *"Remap target files"*: `Originals and subsampled targets`
>                - {% icon param-file %} *"Original files (LAZ/LAS)"*: `output_original_copc` (output of **3DTrees: SmartTile** {% icon tool %})
>                - {% icon param-file %} *"Subsampled target files (LAZ/LAS)"*: `output_subsample_res1` (output of **3DTrees: SmartTile** {% icon tool %})
>            - *"Also write one merged output (LAZ)"*: `Yes`
>                - *"Enrich merged output with original attributes"*: `Yes`
>                - {% icon param-file %} *"Standardization summary (optional)"*: `collection_summary` (output of **3DTrees: LAS/LAZ Standardization** {% icon tool %})
>
>
>    > <comment-title> See exactly how the tols is filled in </comment-title>
>    >
>    > ![Image of how to fill in the tool form step 1](../images/3dtrees/form1.png)
>    > ![Image of how to fill in the tool form step 2](../images/3dtrees/form2.png)
>    > ![Image of how to fill in the tool form step 3](../images/3dtrees/form3.png)
>    {: .comment}
>
{: .hands_on}

# Point cloud visualization by Potree

> <hands-on-title> 3Dtrees: Potree Converter </hands-on-title>
>
> Convert segmented collection to Potree octree format
>
> 1. {% tool [3Dtrees: Potree Converter](toolshed.g2.bx.psu.edu/repos/bgruening/3dtrees_potree/3dtrees_potree/1.0.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Point Cloud Files"*: `output_merged_with_originals` (output of **3DTrees: SmartTile** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

# Conclusion

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.