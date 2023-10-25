---
layout: tutorial_hands_on
title: REMBI - Recommended Metadata for Biological Images – metadata guidelines for bioimaging data

zenodo_link: ''

questions:
- What is REMBI and why should I use it?
- What information should be included when collecting bioimage data?

objectives:
- Organise bioimage metadata
- Find out what REMBI is and why it is useful
- Categorise what metadata belongs to each of the submodules of REMBI
- Gather the metadata for an example bioimage dataset
  
time_estimation: "15m"

key_points:
- REMBI describes useful guidelines for bioimaging that can help unification and FARIfication of the data.

tags:
- fair
- data management
- bioimaging
  
priority: 5

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
      - bioimage-metadata

follow_up_training:
  - type: "internal"
    topic_name: imaging
---

# Metadata guidelines for bioimaging data

REMBI (Recommended Metadata for Biological Images) was proposed as a draft metadata guidelines to begin addressing the needs of diverse communities within light and electron microscopy. Currently, these guidelines are in draft form to encourage discussion within the community, but they provide a useful guide as to what metadata should be gathered to make your image data FAIR. They divide the metadata requirements into eight modules which further split into attributes - that seems to be a daunting task, doesn't it? But at the same time it's exciting news for the community! To find out more, have a look at the [REMBI article](https://www.nature.com/articles/s41592-021-01166-8).


> <question-title></question-title>
>
> In the [REMBI paper](https://www.nature.com/articles/s41592-021-01166-8), the authors consider three potential user groups who require different metadata. Find out what are these three groups and their metadata requirements.
>
> > <solution-title></solution-title>
> > The identified three user groups are: Biologists, Imaging scientists, Computer-vision researchers. 
> > - A research biologist may be interested in the biological sample that has been imaged to compare it to similar samples that they are working with.
> > - An imaging scientist may be interested in how the image was acquired so they can improve upon current image acquisition techniques.
> > - A computer vision researcher may be interested in annotated ground-truth segmentations, that can be obtained from the image, so they can develop faster and more accurate  algorithms.
> {: .solution}
>
{: .question}

> <tip-title>Instructor Note</tip-title>
>
> If you're an instructor leading this training, you might ask people to work in small groups for this exercise and encourage the discussion. Ask group members to share which of the user groups they identify as and what metadata they would want.
> 
{: .tip}

# Categories of metadata 
REMBI covers different categories of metadata, such as: 
- study
- study component
- biosample
- specimen
- image acquisition
- image data
- image correlation
- analyzed data

Within each module, there are attributes that should be included to make the published data FAIR. We will explore all the modules and attributes suggested by REMBI and we'll show some examples as well. 

## Study 
The first module of REMBI metadata describes the Study and should include:
- Study type
- Study description
- General dataset information

### Study type

Ideally, the study type will be part of an ontology. You can look up the main subject of your study using a tool like [OLS](https://www.ebi.ac.uk/ols/index) to find a suitable ontology. This will help others to see where your study sits within the wider research area.

> <comment-title>Example</comment-title>
>
> <table>
>   <tr>
>     <th>Study type</th>
>     <td>Regulation of mitotic cell division</td>
>   </tr>
> </table>
>
{: .comment}


### Study description

A brief description of the project. The Study Description should include the title of the study, a brief description and any related publication details such as authors, title and  DOI. If you are gathering metadata prepublication, you can fill in the publication details later or enter a draft title or the journal name you plan to submit to. It’s still a good idea to include the category, so you don’t forget.

> <comment-title>Example</comment-title>
> <table>
>   <tr>
>     <th colspan="2">Study description</th>
>   </tr>
>   <tr>
>     <th>Title</th>
>     <td>Imaging mitotic cells</td>
>   </tr>
>   <tr>
>     <th>Description</th>
>     <td>Visualising HeLa cells using confocal microscopy</td>
>   </tr>
>   <tr>
>     <th>Publication details</th>
>     <td>TBC</td>
>   </tr>
> </table>
{: .comment}


### General dataset information

This should include all the information that relates to all the data in the project. This can include the names of contributors and the repository where the data is or will be stored. State the licence under which you intend to make the data available, the repository you intend to submit to and if you are using a schema for structuring your metadata. This helps to keep all collaborators on the same page. Any other general information with respect to the study can be included here, but try to keep this broad as more detailed information should be included in other sections of the metadata.

> <comment-title>Example</comment-title>
> <table>
>   <tr>
>     <th colspan="2">General Dataset Information</th>
>   </tr>
>   <tr>
>     <th>Contributors</th>
>     <td>Alica and Bob</td>
>   </tr>
>   <tr>
>     <th>Repository</th>
>     <td>Bioimage Archive</td>
>   </tr>
>   <tr>
>     <th>Licenses</th>
>     <td>CC-BY</td>
>   </tr>
>   <tr>
>     <th>Schemas</th>
>     <td>Datacite Metadata</td>
>   </tr>
> </table>
{: .comment}

## Study component

A study component can be thought of as an experiment, both the physical experiment and subsequent data analysis, or a series of experiments that have been conducted with the same aim in mind.

The associated metadata should describe the imaging method used and include a description of the image dataset. The REMBI guidelines store high-level metadata in the study component and then divide the more detailed metadata into other modules. 

Within the Study component we include the Imaging Method which should describe the techniques used to acquire the raw data. This could be one or multiple methods, which should be part of a relevant ontology. For Confocal Microscopy data, we can use the Biological Imaging Methods Ontology, although it is also present in a number of other ontologies.

The description of the study component should include an overview of what was imaged as well as any processed data that is created during analysis. 

> <comment-title>Example</comment-title>
> <table>
>   <tr>
>     <th>Imaging Method</th>
>     <td>Confocal Microscopy</td>
>   </tr>
>   <tr>
>     <th>Study Component Description</th>
>     <td>Images of cells and segmented binary masks</td>
>   </tr>
> </table>
{: .comment}


> <tip-title>Storing metadata</tip-title>
>
> You could either choose to store the metadata in the same file as your study data or have a new file for each study component. This could be stored in the same place as your study metadata, or you could create a subdirectory structure. 
{: .tip}

## Biosample

The first thing you need for the biosample metadata is an Identity. This is a code that you assign to each sample you are describing, which will link this metadata to the physical sample. Then, state what the biological entity is, which should come from a relevant ontology. Use a taxonomy to name the organism. Next, describe the variables in your experiment. The REMBI guidelines split the variables into three types:
- intrinsic - describe an innate trait of the biosample, such as a genetic alteration
- extrinsic - describe something you added to the sample, for example, a reagent
- experimental - things that you intentionally vary, like time

You can leave out some of the variables if they are not part of your experiment.

> <comment-title>Example</comment-title>
> <table>
>   <tr>
>     <th>Identity</th>
>     <td>CM001</td>
>   </tr>
>   <tr>
>     <th>Biological entity</th>
>     <td>JURKAT E-6.1 cell</td>
>   </tr>
>   <tr>
>     <th>Organism</th>
>     <td>Homo sapiens</td>
>   </tr>
>   <tr>
>     <th>Intrinsic variable</th>
>     <td>Jurkat E6.1 transfected with emerald-VAMP7</td>
>   </tr>
>   <tr>
>     <th>Extrinsic variable</th>
>     <td>Aspirin</td>
>   </tr>
>   <tr>
>     <th>Experimental variables</th>
>     <td>Dose response of aspirin</td>
>   </tr>
> </table>
{: .comment}

## Specimen

 The specimen metadata should include:
- the experimental status (control or test)
- the location within the biosample, such as a coordinate or a particular well in a plate
- how the sample was prepared
- how the signal is being generated
- the content and biological entities of different channels.

 Include enough information so that someone with experience in the field could reproduce a sample by following the information you provided. Assume they would know typical techniques and name them using terms from an ontology if possible. Only include lots of detail if you are describing a novel technique.

> <comment-title>Example</comment-title>
> <table>
>   <tr>
>     <th>Experimental status</th>
>     <td>Control</td>
>   </tr>
>   <tr>
>     <th>Location within biosample</th>
>     <td>Plasma membrane within 100 nm of coverslip (TIRF)</td>
>   </tr>
>   <tr>
>     <th>Preparation method</th>
>     <td>Cos-7 cells cultured in DMEM medium, and then plated on #1 coverslips and imaged live in L-15 medium</td>
>   </tr>
>   <tr>
>     <th>Signal/contrast mechanism</th>
>     <td>fluorescent proteins</td>
>   </tr>
>   <tr>
>     <th>Channel – content</th>
>     <td>Green: eGFP, Red: mCherry</td>
>   </tr>
>   <tr>
>     <th>Channel – biological entity</th>
>     <td>Green: EGFR, Red: Src</td>
>   </tr>
> </table>
{: .comment}

## Image acquisition
Here you should include all the information about the instrument you used and how it was set up. Like with the specimen metadata, describe this information as though you are speaking to someone who already knows how to use a similar instrument. What would they need to know to produce the same image data?
 
Check with your facility manager if they have any guidelines for what details need to be recorded for your particular instrument. Make sure that the parameters you record can actually be used by someone else if they don’t have exactly the same instrument or setup. For example, don’t say that you used a certain percentage of laser power, as this doesn’t tell you how much power was used unless you also provide the total power of the laser. If the instrument software has automatically generated a metadata file, remember to save this. Depending on its content, this may be sufficient.
 
Start with the details of the equipment for the Instrument Attributes. If this is commercial equipment, include the make and model, a short description of what type of instrument it is and details about its configuration. If the instrument is bespoke, you will need to include more details. Next, you should include image acquisition parameters. These relate to how the instrument was set up for the particular experiment. Some of these may be captured automatically by the instrument’s software, so make things easy for yourself and check if a file is generated and what’s in it. If a file is generated, then you only need to manually record anything that is missing from the file.

> <comment-title>Example</comment-title>
> <table>
>   <tr>
>     <th>Instrument attributes</th>
>     <td>Olympus FV3000, laser point scanning confocal, 500-550 nm filter, 37-degree chamber.</td>
>   </tr>
>   <tr>
>     <th>Image acquisition parameters</th>
>     <td> </td>
>   </tr>
>   <tr>
>     <th>Objective</th>
>     <td>Cos-7 cells cultured in DMEM medium, and then plated on #1 coverslips and imaged live in L-15 medium</td>
>   </tr>
>   <tr>
>     <th>Excitation Wavelength</th>
>     <td>488 nm</td>
>   </tr>
>   <tr>
>     <th>PMT gain</th>
>     <td>500 V</td>
>   </tr>
>   <tr>
>     <th>Pixel dwell time</th>
>     <td>2 𝜇s</td>
>   </tr>
>   <tr>
>     <th>Confocal aperture</th>
>     <td>200 𝜇m</td>
>   </tr>
> </table>
{: .comment}

> <tip-title>Helpful resources</tip-title>
>
> To help you collect the information for your own data, you might have a look at the local resources from your institution or universities. For example, at Warwick University, there are [webpages](https://warwick.ac.uk/fac/sci/med/research/biomedical/facilities/camdu/methodsreporting/) describing the metadata that needs to be collected for some of the microscopes. 
> 
{: .tip}

## Image data

In this section, you record all the information related to all the images you have. Not only the primary or raw images, but also any processed images,  perhaps such as binary files showing the resulting segmentation.

You need to say what format the images are in and if they have undergone any compression, the dimensions of the images, and what the physical size of the pixel or voxel is, including the units. Most of this information you should be able to get from the metadata or header of the image files.

Next, you need to state the physical size of the image or magnification, calculated from the pixel or voxel size and the dimension extents. Give any information related to how the channels are represented. For processed images, you need to provide the methods used for processing.

Finally, say you have used contrast inversion, do the bright features in the image correspond to areas of high signal, or is it the other way around?

> <comment-title>Example</comment-title>
> <table>
>   <tr>
>     <th>Type</th>
>     <td>Primary Image, Segmentation</td>
>   </tr>
>   <tr>
>     <th>Format and compression</th>
>     <td>Primary: .oir (Olympus), Segmentation: .tiff</td>
>   </tr>
>   <tr>
>     <th>Dimension extents</th>
>     <td>x: 512, y: 512, z: 25</td>
>   </tr>
>   <tr>
>     <th>Size description</th>
>     <td>153.6 x 153.6 x 25 𝜇m</td>
>   </tr>
>   <tr>
>     <th>Pixel/Voxel size description</th>
>     <td>0.3 x 0.3 x 0.1 𝜇m</td>
>   </tr>
>   <tr>
>     <th>Image processing method</th>
>     <td>Fiji: Median filter (3 pixel kernel), Otsu threshold</td>
>   </tr>
>   <tr>
>     <th>Contrast inversion</th>
>     <td>No</td>
>   </tr>
> </table>
{: .comment}


## Image correlation

If you have used different imaging modalities with the same sample, this part of the metadata should describe how the images relate to one another. You could use this section to describe generally the relationship between images. In the example below, images from different modalities have been aligned.

> <comment-title>Example</comment-title>
> <table>
>   <tr>
>     <th>Spatial and temporal alignment</th>
>     <td>Manual</td>
>   </tr>
>   <tr>
>     <th>Fiducials used</th>
>     <td>Soil grains</td>
>   </tr>
>   <tr>
>     <th>Transformation matrix</th>
>     <td>See file: Transforms.csv</td>
>   </tr>
>   <tr>
>     <th>Size description</th>
>     <td>153.6 x 153.6 x 25 𝜇m</td>
>   </tr>
>   <tr>
>     <th>Related images and relationship</th>
>     <td>Primary XCT: Data/XCT
>       Primary XRF: Data/XRF
>       Processed XRF: Data/Transformed_XRF</td>
>   </tr>
> </table>
{: .comment}

## Analysed data

This section should not include metadata for any image data, including processed images, as that should have been covered in the Image Data section. Instead, it should describe the analysis results you have, such as measurements. Have you done some numerical analysis or some phenotyping or something else? There is no need to describe the methods in great detail if they are already described in the relevant publication.

> <comment-title>Example</comment-title>
> <table>
>   <tr>
>     <th>Analysis results type</th>
>     <td>Speed of cell division</td>
>   </tr>
>   <tr>
>     <th>Data used for analysis</th>
>     <td>Preprocessed images, Cell tracks</td>
>   </tr>
>   <tr>
>     <th>Analysis method</th>
>     <td>Track cell lineage: BayesianTracker (btrack) with configuration track_config.json
>       Measure speed: Numerical analysis in Python</td>
>   </tr>
> </table>
{: .comment}

# Final notes

For more examples, check out REMBI Supplementary Information - either in [pdf](https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-021-01166-8/MediaObjects/41592_2021_1166_MOESM1_ESM.pdf) or [spreadsheet](https://docs.google.com/spreadsheets/d/1Ck1NeLp-ZN4eMGdNYo2nV6KLEdSfN6oQBKnnWU6Npeo/edit#gid=1023506919).


At first glance, it might seem to be quite a stretch to collect all that metadata! But don’t get discouraged - following those guidelines will ensure better communication between the scientists and will make your research FAIR: Findable, Accessible, Interoperable, Reusable. During big data era when we are surrounded by so much resources, it’s crucial to get good data management habits, share them with others and hence contribute to the development of Science toghether. 
