---
layout: tutorial_hands_on
title: FAIR Bioimage Metadata

zenodo_link: ''

questions:
- What are the commonly used repositories for bioimaging data?
- Which repositories are suitable for my data?
- What are the requirements for submitting?

objectives:
- Locate bioimage data repositories
- Compare repositories to find which are suitable for your data
- Find out what the requirements are for submitting

time_estimation: "15m"

key_points:
- Data repositories such as BioImage Archive, Electron Microscopy Public Image Archive (EMPIAR) and Image Data Repository (IDR) are available to help make bioimaging data FAIR.
- Find out what are the repository's requirements to help decide which is suitable for your data.
- All repositories require some metadata, so collecting the metadata alongside data acquisition will make this process easier.

tags:
- fair
- data management
- bioimaging
  
priority: 4

contributions:
  authorship:
    - wee-snufkin
    - Laura190
    - kkamieniecka
    - poterlowicz-lab
  funding:
      - ELIXIR-UK-DaSH
subtopic: fair-data

requirements:
  - type: "internal"
    topic_name: fair
    tutorials:
      - fair-intro
      - data-management
   

follow_up_training:
  -
    type: "internal"
    topic_name: fair
    tutorials:
        - bioimage-REMBI
  -
    type: "internal"
    topic_name: imaging
        
---

# FAIR Bioimaging

Submitting your data to a repository is a good way to make the data FAIR. This will make it:
- **F**indable, as the data will be given specific identifiers
- **A**ccessible, as the data will be available online, open and free where possible
-	**I**nteroperable, as the repository will often enforce the use of formalised, consistent language
-	**R**eusable, as the data will be released under a license with detailed provenance

# Examples of bioimage data repositories

But the question remains: where can I submit my data? Currently the main repositories where you can submit your images are: BioImage Archive, Electron Microscopy Public Image Archive (EMPIAR) and Image Data Repository (IDR). Let's have a look at the questions below to explore those repositories more in depth!

> <question-title></question-title>
>
> Listed below are three examples of Bioimage Data Repositories:
> - [IDR: Image Data Repository](https://idr.openmicroscopy.org/)
> - [EMPIAR: Electron Microscopy Public Image Archive](https://www.ebi.ac.uk/empiar/)
> - [BioImage Archive](https://www.ebi.ac.uk/bioimage-archive/)
> 
> Visit their websites and find out what their scope is or what sorts of datasets they accept.
>
> > <solution-title></solution-title>
> >
> > - **IDR: Image Data Repository**: Curated datasets of cell and tissue microscopy images
> > - **EMPIAR: Electron Microscopy Public Image Archive**: Cryo-EM, Scanning Electron Microscopy, Soft X-ray tomography
> > - **BioImage Archive**: Everything else and some overlap with IDR and EMPIAR
> {: .solution}
>
{: .question}

> <tip-title>Repositories everywhere</tip-title>
>
> As well as these repositories, your Institute may have their own repository. For example, at the Warwick University, there is also [OMERO](https://warwick.ac.uk/fac/sci/med/research/biomedical/facilities/camdu/training/omero-warwick-guide_2.pdf) and [WRAP](https://wrap.warwick.ac.uk/).
> 
{: .tip}


> <comment-title></comment-title>
> The repositories we are looking at in this course are for bioimage data, not medical data. There are other specialist repositories available if you have medical data.
{: .comment}

# Things to consider when choosing a repository

Now we know what repositories are available, but how to decide which one is best given the files we want submit? Try to work through the below Question box and find the answer!

> <question-title></question-title>
>
> Choose one repository from above and look through its documentation. Try to find:
> 1. What data formats are accepted?
> 2. What license is recommended to publish the data?
> 3. Are there specific instructions for large datasets?
>
> > <solution-title></solution-title>
> >
> > - **IDR: Image Data Repository**:
> >   1. The IDR uses the Bio-Formats library for reading imaging data. Bio-Formats supports over 150 proprietary and open file formats (see the [full list](https://bio-formats.readthedocs.io/en/stable/supported-formats.html)).
> >   2. It is strongly recommended that submitters make their datasets available under [CC-BY](https://creativecommons.org/licenses/by/4.0/) license.
> >   3. As specified on the [IDR website](https://idr.openmicroscopy.org/about/submission.html), dataset size is typically not an issue, but for sizes significantly larger than 1000 GB special planning may be needed.
> > - **EMPIAR: Electron Microscopy Public Image Archive**:
> >   1. Provide image data in the formats in which they are uploaded, but recommended is the use of common formats in the field including MRC, MRCS, TIFF, DM4, IMAGIC, SPIDER, MRC FEI, RAW FEI and BIG DATA VIEWER HDF5. 
> >   2. All data in EMPIAR is freely and publicly available to the global community under the [CC0](https://creativecommons.org/share-your-work/public-domain/cc0/) license.
> >   3. As specified on the [EMPIAR page](https://www.ebi.ac.uk/empiar/deposition/manual/#manIntro), typically having more than 4000 files in a directory has a tendency to slow down access considerably. It is recommended in this case to sub-divide the directory into subdirectories with no more than 4000 files each. If you have a single file larger than 1 TB, contact EMPAIR in advance.
> >      To find out more, check the [FAQ page](https://www.ebi.ac.uk/empiar/faq).
> > - **BioImage Archive**:
> >   1. The BioImage Archive accepts all image data formats, although formats readable by [Bio-Formats library](https://bio-formats.readthedocs.io/en/stable/supported-formats.html), are preferable.
> >   2. According to [BioImage Archive Policies](https://www.ebi.ac.uk/bioimage-archive/help-policies/), all new data directly submitted to the BioImage Archive will be made available under a [CC0](https://creativecommons.org/share-your-work/public-domain/cc0/) licence, datasets brokered/imported from other resources may have other licenses though.
> >   3. There are different submission methods depending on data size:
> >      - Less than 50 GB total size, less than 20GB per file – use submission tool
> >      - Up to 1TB total size – use FTP
> >      - Anything larger – use Aspera
> > 
> >        To find out more, check the [FAQ page](https://www.ebi.ac.uk/bioimage-archive/help-faq/).
> {: .solution}
>
{: .question}


# What metadata to collect

Whichever repository you choose, you will be required to upload some metadata along with your data. In an ideal world, you would remember everything about your data when you submit it. In reality, this is unlikely, the data could have been collected over a long time period or by different people. To overcome these challenges, it is best to collect metadata alongside imaging experiments, don’t leave it all to the end!
However, this raises further challenges. At the time of data acquisition, you probably won’t know which repository you will submit it to, what the study results will be, or even who will be the target audience for the data. So what metadata do you need to collect?

Currently, there is no standard for bioimages, so here is the general outline how to proceed:
- If you have chosen a repository, use their template/guidelines.
- Otherwise, use [REMBI](https://www.nature.com/articles/s41592-021-01166-8). These are published guidelines which are explained in detail in the [REMBI tutorial]({% link topics/fair/tutorials/bioimage-REMBI/tutorial.md %}). REMBI is useful as it should cover most of the metadata requirements of the repositories, even if you haven’t decided which one you want to use yet.
- For medical images, see the [DICOM](https://www.dicomstandard.org/) standard.

# How to store metadata

Metadata should be stored somewhere that it can be viewed and edited by collaborators. This helps everyone to stay on the same page with regard to what data should and has been collected. This is also useful if different people contribute to different aspects of the image acquisition, e.g. maybe one person prepared the sample, another imaged it under the microscope, and a third person did the post-processing. Each person can then update the metadata related to their part of the study.

Recommendations for storing metadata:
- If a repository has been chosen, use their template (if provided).
- Use a Delimited text file format, e.g. .csv, .tsv. You can use spreadsheet software and save to this format. Try to use a plain format, e.g. avoid merged/split cells.
-  Use data management software suitable for your data, e.g. [OMERO](https://www.openmicroscopy.org/omero/).

# Further steps
{% icon congratulations %} Congratulations on successfully completing this tutorial! If you want to know more about FAIR data management, we provide training that you can find in the [FAIR Data, Workflows, and Research]({% link topics/fair %}) section on GTN. For more detials on FAIR data management in bioimaging, see the [REMBI tutorial]({% link topics/fair/tutorials/bioimage-REMBI/tutorial.md %}), and if you want to dive into the imaging analysis straight away, feel free to choose one of the Galaxy [tutorials]({% link topics/imaging %})! Additionally, [Global BioImaging](https://globalbioimaging.org/) offers [training courses](https://globalbioimaging.org/international-training-courses/repository/image-data) dedicated to image data management, sharing, reuse, and image data repositories that were mentioned in this tutorial.
