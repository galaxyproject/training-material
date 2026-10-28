---
layout: tutorial_hands_on

title: '3Dtrees: From LiDAR point clouds to individual trees'
level: Introductory
zenodo_link: ''
questions:
- How can LiDAR describe the three-dimensional structure of a forest?
- How can individual trees be identified in a point cloud?
- How can segmentation results be inspected in Galaxy?
objectives:
- Explain what a LiDAR point cloud represents
- Standardize and prepare one LAS file for tree segmentation
- Distinguish tree/non-tree classes from individual-tree identifiers
- Inspect tree predictions and identify obvious visual problems
time_estimation: 3H
subtopic: analysing
tags:
- earth-system
- land
- lidar
key_points:
- LiDAR records the 3D structure of vegetation and terrain
- SegmentAnyTree predicts tree and non-tree points and assigns an ID to each predicted tree
- Visual inspection can reveal obvious problems but is not a quantitative accuracy assessment
contributions:
  authorship:
  - kgerb
  - Marie59
requirements:
- type: internal
  topic_name: introduction
  tutorials:
  - galaxy-intro-101-everyone

---

Forests are three-dimensional. A view from above shows the canopy, but tells us less about the
vegetation layers below it. LiDAR (Light Detection and Ranging) uses laser pulses to measure
positions on vegetation and the ground. These measurements form a *point cloud*: a set of points
with horizontal coordinates and an elevation.

LiDAR collected from an aircraft or drone can help ecologists study tree height, crown size,
canopy gaps, and changes in forest structure. To study individual trees, we first need to group
the points that belong to each tree. This is called *individual-tree segmentation*.

This tutorial is for interested beginners, ecology students, and ecologists. You will use the
3Dtrees tools ({% cite Gerberding2026 %}) in Galaxy to prepare one point cloud and predict individual trees with
SegmentAnyTree (SAT) ({% cite Wielgosz2024 %}). At the end, you will have a point cloud with tree
labels and an interactive Potree view for inspecting the results in your browser.

> <tip-title>New to Galaxy?</tip-title>
>
> If Galaxy is new to you, first try [Galaxy Basics for everyone]({% link topics/introduction/tutorials/galaxy-intro-101-everyone/tutorial.md %}).
>
{: .tip}

> <comment-title>Which Galaxy server?</comment-title>
>
> This tutorial was prepared on [Galaxy Europe](https://usegalaxy.eu), where the 3Dtrees tools and the
> Recherche Data Gouv remote-file source used below are available. The screenshots come from its
> [Earth System](https://earth-system.usegalaxy.eu) view, which offers the same tools.
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

# From point clouds to individual trees

We will import one file, standardize it, and make quick overview images. Next, we will prepare
smaller processing units, predict trees, and transfer the predictions to a denser point cloud.
Finally, we will inspect the result in 3D. Each step prepares an output for the next one.

## Get data

The data were collected by drone over a Mediterranean oak forest near Pic Saint-Loup in southern
France ({% cite Alleaume2021 %}). We use only the smallest file, `Strip-20210804-150144_D.las`
(about 477 MB), to keep the tutorial manageable. LAS is a common format for LiDAR points; LAZ is
its compressed form. Both can store extra information for each point, such as color or a predicted
tree identifier.

> <details-title>About the Pic Saint-Loup dataset</details-title>
>
> The data were collected in June 2021 with a YellowScan Surveyor scanner on a DJI Matrice 600 Pro
> drone, as part of the SENTHYMED/MEDOAK experiment ({% cite Adeline2024 %}). The dates in the file
> names record when post-processing finished, not the flight date. The points use the RGF93 /
> Lambert-93 coordinate reference system. The dataset has about 200--250 points per square metre and
> is available under the Etalab Open License 2.0.
>
> The four files total approximately 4.1 GB. The associated paper states that they match four flight
> areas. A and B are spatially separate, while the published bounding boxes of C and D overlap;
> consequently, these files should not be treated as four non-overlapping tiles of one continuous
> survey.
>
> See the dataset record ({% cite Alleaume2021 %}) and its documentation ({% cite Adeline2024 %})
> for provenance and reuse information.
>
{: .details}

> <hands-on-title>Import one flight-area file</hands-on-title>
>
> 1. Create a new history for this tutorial.
>
>    {% snippet faqs/galaxy/histories_create_new.md box_type="none" %}
>
> 2. Import the file from its persistent download address:
>
>    ```
>    https://entrepot.recherche.data.gouv.fr/api/access/datafile/121030
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md format="laz" %}
>
>    > <comment-title>Why the laz datatype for a LAS file?</comment-title>
>    >
>    > The file is a LAS file, but the Standardization tool currently only offers datasets of type
>    > `laz` in its form. It reads LAS content correctly, and it writes a real LAZ file that all
>    > later steps use.
>    {: .comment}
>
>    > <tip-title>Browsing the repository instead</tip-title>
>    >
>    > Galaxy Europe can also browse the dataset: in **Upload**, select **Choose from repository**,
>    > open **Recherche Data Gouv**, search for `pic` and open
>    > **Données lidar acquises par drone au pied du Pic St Loup (34)**. In our test this import
>    > stopped early and produced a file of only 187 MB, without an error. If you use it, check that
>    > the dataset is 476.9 MB before you continue.
>    >
>    > ![Repository browser showing the four LAS files of the Pic Saint-Loup dataset](../../images/3dtrees/rdg_files.png "The four flight-area files in Recherche Data Gouv. Only file D is needed")
>    {: .tip}
>
> 3. Wait for the dataset to turn green in the history, then expand it and check that its size is
>    **476.9 MB**. A smaller file means the download stopped early; delete it and import it again.
>
>    ![Galaxy history with the dataset Strip-20210804-150144_D.las, 476.9 MB, format laz](../../images/3dtrees/data_inhistory.png "The imported file D in the history"){: width="35%"}
>
> 4. Rename the dataset to `Strip-20210804-150144_D.las` and check that its datatype is `laz`.
>    Do not build a four-file collection for this tutorial; we use only file D.
>
>    {% snippet faqs/galaxy/datasets_rename.md name="Strip-20210804-150144_D.las" box_type="none" %}
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="laz" %}
{: .hands_on}


# Standardize the point cloud

Files from different scanners can store point information in different ways. Standardization
checks the file and prepares a consistent format for the next tools. It checks the file header,
spatial bounds, point density, and coordinate reference system, which describes where the points
are located on Earth.

The tool writes a standardized LAZ file, a metadata summary, and a footprint showing the area
covered by the survey. This helps prevent file-format problems later; it does not identify trees.

> <hands-on-title>Standardize the example LAS file</hands-on-title>
>
> 1. {% tool [3DTrees: LAS/LAZ Standardization](toolshed.g2.bx.psu.edu/repos/bgruening/3dtrees_standardization/3dtrees_standardization/1.2.1+galaxy0) %} with the following parameters:
>    - *"Mode"*: `Single File Standardization`
>        - {% icon param-file %} *"Point Cloud File"*: `Strip-20210804-150144_D.las`
>
>    ![Standardization form with Mode set to Single File Standardization and Strip-20210804-150144_D.las as the point cloud file](../../images/3dtrees/standardization_form.png "The Standardization form for file D")
>
>    Three outputs appear in the history: `standardized` (LAZ), `metadata` (JSON), and
>    `convex_hull` (GeoJSON).
>
>    > <comment-title>Choose the output for the next step</comment-title>
>    >
>    > Use `standardized`, the standardized point cloud, in the next tools. Keep the metadata and
>    > footprint (`convex_hull`) as supporting information about the input.
>    {: .comment}
>
{: .hands_on}

{% include _includes/quiz.html id="standardization.yaml" %}

> <details-title>Explanation: Why standardize a point cloud?</details-title>
>
> Standardization checks the file and writes a consistent LAZ file, a metadata summary and a
> footprint. Tree identifiers are created later by the segmentation model. Standardization cannot
> prove that those predictions are accurate, and it does not turn separate flight areas into one
> continuous survey.
>
{: .details}

# Get a first impression with overview images

Before predicting trees, take a quick look at the input. The Overviews tool creates images you
can open directly in Galaxy, without a 3D viewer. Top views show the survey shape and canopy gaps.
Side views help you recognize the ground, taller vegetation, and layers within the forest.

> <hands-on-title>Create and inspect overview images</hands-on-title>
>
> 1. {% tool [3Dtrees: Overviews](toolshed.g2.bx.psu.edu/repos/bgruening/3dtrees_overviews/3dtrees_overviews/1.2.0+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Point Cloud Collection"*: `standardized` (output of **3DTrees: LAS/LAZ Standardization** {% icon tool %})
>
> 2. Open the **Top View Images** and **Section View Images** collections in the history and view
>    their images. With the default settings, expect four images in total.
>
>    > <comment-title>Reading the overview</comment-title>
>    >
>    > The input field is called **Point Cloud Collection**, but select only `standardized` for
>    > this tutorial. The output collections simply group the overview images.
>    {: .comment}
>
{: .hands_on}

![Galaxy preview of the section_ns image, where the sloping transect occupies only a thin strip in the middle of the frame](../../images/3dtrees/sectionns.png "The north--south section view as Galaxy shows it")

In the Galaxy preview, the transect fills only a thin strip of the image: the survey is about 270 m
long but the forest is less than 15 m tall. Open the image at full size, or download it, and zoom in
to see the canopy. The images below show what to look for. They use the same idea as the Overviews
tool, a north--south slice through the middle of file D, drawn at true scale and coloured by
height above the ground.

![Side view of a 10 m wide north--south slice through file D. The ground rises from about 378 m to 430 m along 270 m, and two boxes labelled A and B mark the close-ups](../../images/3dtrees/section_ns_full.png "The whole north--south transect. Boxes A and B mark the two close-ups")

![Close-up A, 40 m long and 4 m wide: several separate crowns 6 to 9 m tall with gaps between them, sparse understory points, and a clear ground line](../../images/3dtrees/section_ns_cutout_a.png "Close-up A: taller, separate crowns with gaps between them")

![Close-up B, 40 m long and 4 m wide: a lower, continuous canopy of 3 to 6 m on a steep slope, with crowns touching each other and few points between canopy and ground](../../images/3dtrees/section_ns_cutout_b.png "Close-up B: lower, continuous canopy on steeper ground")

> <details-title>How the close-ups were made</details-title>
>
> The close-ups were drawn outside Galaxy from `Strip-20210804-150144_D.las`, to show detail that is
> hard to see in the default overview. The full transect is 10 m wide, like the default
> **Section Width** of the Overviews tool. The close-ups use a 4 m wide slice, so fewer crowns lie
> behind each other. Height above ground is estimated from the lowest points along the transect and
> is only used for colour.
>
{: .details}

{% include _includes/quiz.html id="overview.yaml" %}

> <details-title>Explanation: What can the overview tell you?</details-title>
>
> The top view shows gaps in the canopy. A section view shows how the points are arranged above
> the ground. In close-up B, the lower crowns touch each other and form an almost continuous canopy,
> so there are few gaps to show where one tree ends and the next begins. Keep both areas in mind
> when you inspect the predictions in Potree.
>
{: .details}

# Prepare the data for segmentation

Large point clouds can be too big to process at once. SmartTile can divide them into smaller
spatial *tiles*. When tiles are created, an overlapping margin, or *buffer*, gives the model more
context around trees near each edge.

SmartTile also reduces the number of points, a step called *subsampling*. Points are grouped in
small cubes, called *voxels*, and each voxel is replaced by one point. SegmentAnyTree works on a
10 cm voxel size, which reduces processing work but also removes fine detail. The predictions are
transferred back to every point of the original file later.

Keep the default settings for this tutorial. The **Tiling threshold (MB)** (default 10000 MB)
means that a single file smaller than this, such as file D, stays as one processing unit.
Subsampling still takes place; you do not need multiple tiles to continue.

> <hands-on-title>Prepare processing units and subsampled points</hands-on-title>
>
> 1. {% tool [3DTrees: SmartTile](toolshed.g2.bx.psu.edu/repos/bgruening/3dtrees_smart_tile/3dtrees_smart_tile/2.4.1+galaxy0) %} with the following parameters:
>    - *"Task"*: `Tile and subsample`
>        - {% icon param-file %} *"Input LAZ/LAS files"*: `standardized` (output of **3DTrees: LAS/LAZ Standardization** {% icon tool %})
>        - *"Tiling mode"*: `Buffered spatial tiles`
>        - In *"Output"*:
>            - {% icon param-repeat %} *"1: Output"*
>                - *"Resolution"*: `Subsampled`
>                    - *"Voxel size (cm)"*: `10`
>                - *"Format"*: `LAZ`
>
>    ![SmartTile form with Task Tile and subsample, input standardized, buffered spatial tiles, tile length 300 m, buffer 20 m, and one subsampled output at 10 cm in LAZ format](../../images/3dtrees/smarttile_tile_form.png "SmartTile set up to tile and subsample file D at 10 cm")
>
>    Three outputs appear in the history:
>    - `Tiles, output 1`: the point cloud subsampled at 10 cm, one file per processing unit
>    - `Tiling preview`: a map of the processing units
>    - `Tile bounds JSON`: the processing boundaries, needed to put the predictions back together
>
>    > <comment-title>Keep the outputs for later steps</comment-title>
>    >
>    > Use `Tiles, output 1` for SegmentAnyTree and keep `Tile bounds JSON` for the merge step.
>    > The collection comes from your single input file; its elements are not the four original
>    > flight areas.
>    {: .comment}
>
{: .hands_on}

{% include _includes/quiz.html id="tiles.yaml" %}

> <details-title>Explanation: Why use overlapping processing tiles?</details-title>
>
> A buffer gives the model points from beyond the tile edge, which reduces the chance of cutting
> off a tree's surroundings, but it cannot guarantee that every tree is complete. A wider buffer
> gives more context but repeats more points, so processing takes longer and the same tree may be
> predicted in neighbouring tiles. The merge step removes those duplicates.
>
{: .details}

# Identify individual trees with SegmentAnyTree

SegmentAnyTree is a deep-learning model trained on laser scans from many forest types, sensors, and
platforms ({% cite Wielgosz2024 %}). It adds two predictions to each point.
*Semantic segmentation* assigns a class: `PredSemantic_SAT` distinguishes tree from non-tree points.
*Instance segmentation* separates individual trees: `PredInstance_SAT` gives the points belonging to
each predicted tree a shared identifier.

For example, two neighboring trees can have the same semantic class but different instance IDs.
Ground points belong to the non-tree class; SAT does not provide a detailed vegetation-type map.

> <hands-on-title>Predict trees on the coarse point cloud</hands-on-title>
>
> 1. {% tool [3Dtrees: SegmentAnyTree](toolshed.g2.bx.psu.edu/repos/bgruening/3dtrees_segmentanytree/3dtrees_segmentanytree/1.2.3+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input Dataset"*: `Tiles, output 1` (output of **3DTrees: SmartTile** {% icon tool %})
>    - *"Predicted instance dimension name"*: `PredInstance_SAT`
>    - *"Predicted semantic dimension name"*: `PredSemantic_SAT`
>
>    Select the collection with the {% icon param-collection %} **Dataset collection** button so
>    Galaxy runs SAT on each processing unit. Keep the default dimension names: the merge and
>    Potree steps look for `PredInstance_SAT`. Processing time depends on the data and on server
>    availability. In our run on Galaxy Europe, the job waited about an hour for a GPU and then ran
>    for about 20 minutes.
>
>    ![SegmentAnyTree form with the collection Tiles, output 1 selected, mapped over the tool, and the dimension names PredInstance_SAT and PredSemantic_SAT](../../images/3dtrees/segmentanytree_form.png "SegmentAnyTree runs once per processing unit")
>
>    > <comment-title>What SAT adds to the point cloud</comment-title>
>    >
>    > The `Processed Files` output collection contains point clouds with `PredSemantic_SAT` and
>    > `PredInstance_SAT` added. These labels are model predictions, not measurements made by the LiDAR
>    > sensor.
>    {: .comment}
>
{: .hands_on}

{% include _includes/quiz.html id="predictions.yaml" %}

> <details-title>Explanation: Instance and semantic predictions</details-title>
>
> `PredSemantic_SAT` describes the class (tree or non-tree). `PredInstance_SAT` identifies the
> individual tree. Both are predictions that SegmentAnyTree adds to the point cloud; the LiDAR
> sensor does not measure them.
>
{: .details}

# Combine predictions and transfer them to the original points

If processing tiles overlap, the same tree can appear in more than one prediction. SmartTile uses
the tile bounds to keep each tree only in the tile that owns it, and joins trees cut by tile
borders. This step addresses processing-tile overlap, not overlap between separate drone flights.
When file D stays as one processing unit, there is no tile overlap to remove, but this step is still
needed to build the final point cloud.

The **Remap target** setting then transfers the predictions from the 10 cm points back to every
point of the original file. The result keeps all original point attributes, such as intensity, and
adds `PredSemantic_SAT` and `PredInstance_SAT`. Transferring labels adds detail to the displayed
cloud, but does not independently predict or validate every added point.

> <hands-on-title>Merge the predictions and transfer them to the original file</hands-on-title>
>
> 1. {% tool [3DTrees: SmartTile](toolshed.g2.bx.psu.edu/repos/bgruening/3dtrees_smart_tile/3dtrees_smart_tile/2.4.1+galaxy0) %} with the following parameters:
>    - *"Task"*: `Filter / merge segmented tiles`
>        - {% icon param-collection %} *"Segmented tiles"*: `Processed Files` (output of **3Dtrees: SegmentAnyTree** {% icon tool %})
>        - *"Voxel size of the segmented tiles (cm)"*: `10`
>        - *"Model"*: `SegmentAnyTree, ForestMamba or other point predictions`
>            - *"Instance dimension"*: `PredInstance_SAT`
>        - {% icon param-file %} *"Tile bounds JSON"*: `Tile bounds JSON` (output of **3DTrees: SmartTile** {% icon tool %})
>        - *"Remap target"*: `Original files`
>            - {% icon param-file %} *"Original LAZ/LAS files"*: `standardized` (output of **3DTrees: LAS/LAZ Standardization** {% icon tool %})
>
>    ![SmartTile form with Task Filter / merge segmented tiles, Processed Files as segmented tiles, voxel size 10 cm, instance dimension PredInstance_SAT, the tile bounds JSON, and Remap target Original files with standardized selected](../../images/3dtrees/smarttile_merge_form.png "SmartTile set up to merge the predictions and transfer them to file D")
>
>    Two collections appear in the history: `Merged tiles`, the filtered predictions at 10 cm, and
>    `Originals with predictions`, file D at full resolution with the predictions added. In our run,
>    this step took about 50 minutes.
>
>    > <comment-title>Instance dimension</comment-title>
>    >
>    > The default **Instance dimension** is `PredInstance`. Change it to `PredInstance_SAT`, the
>    > name SegmentAnyTree used in the previous step, otherwise SmartTile cannot find the tree IDs.
>    {: .comment}
>
{: .hands_on}

Use `Originals with predictions` for visualization. Its single element is file D with the
predictions added; nothing from the other three flight areas is included.

# Inspect the result in Potree

The overview images gave us a quick first look. Potree now lets us rotate and zoom into the final
point cloud in a browser. We can inspect whether predicted tree boundaries look plausible and
look for missed trees, trees joined together, or one tree split into several instances.

> <hands-on-title>Create and open the interactive view</hands-on-title>
>
> 1. {% tool [3Dtrees: Potree Converter](toolshed.g2.bx.psu.edu/repos/bgruening/3dtrees_potree/3dtrees_potree/1.2.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Point Cloud Files"*: `Originals with predictions` (output of **3DTrees: SmartTile** {% icon tool %})
>    - *"Generate interactive Viewer"*: `Yes`
>    - *"Special Coloring"*: `Enabled`
>        - *"Special Coloring Instance Attributes"*: `PredInstance_SAT`
>        - *"Write Special Coloring JSON"*: `Yes`
>
> 2. When the tool finishes, open the **Potree Visualization** HTML output using its
>    {% icon galaxy-eye %} eye icon. Loading the viewer can take a minute.
> 3. The viewer opens with the points coloured by `coloring_id_sat`: each predicted tree has its own
>    colour, and points that belong to no tree are grey.
>
>    ![Potree view of file D from above, with predicted trees in different colours](../../images/3dtrees/potree_instances.png "Predicted trees coloured by coloring_id_sat")
>
> 4. Inspect a few neighbouring trees from above and from the side, for example the areas of
>    close-ups A and B. To see the tree/non-tree classes, open **Scene** in the menu on the left,
>    click the point cloud `potree-vis`, and set **Attribute** to `PredSemantic_SAT`.
>
>    ![Potree view of file D coloured by PredSemantic_SAT, with tree points and non-tree points in two colours](../../images/3dtrees/potree_semantic.png "Tree and non-tree points from PredSemantic_SAT")
>
>    > <comment-title>Special colouring</comment-title>
>    >
>    > Tree IDs are just numbers, so colouring `PredInstance_SAT` directly gives most neighbouring
>    > trees almost the same colour. **Special Coloring** adds `coloring_id_sat`, which gives nearby
>    > trees different colours from a small palette, and **Write Special Coloring JSON** tells the
>    > viewer which colour to use for each value. The colours represent identifiers, not species or
>    > confidence.
>    {: .comment}
>
{: .hands_on}

{% include _includes/quiz.html id="potree.yaml" %}

> <details-title>Explanation: Does the segmentation look plausible?</details-title>
>
> Colours that follow single crowns, different colours on neighbouring crowns, and ground points in
> the non-tree class all support a plausible result. Straight cuts along processing boundaries and
> one colour across unrelated crowns suggest errors. Potree helps you find visible problems, but it
> does not calculate segmentation accuracy: that needs independently labelled reference trees.
>
{: .details}

> <warning-title>Results depend on the input</warning-title>
>
> Results may change with forest type, point density, and survey method. A result that looks good
> may still be wrong. Measuring accuracy requires labelled reference trees.
>
{: .warning}

# Optional next steps

Once the predictions have been checked, individual-tree labels can support studies of tree
height, crown size, and tree density. The wider point cloud also provides information about canopy
gaps and forest structure.

Other 3Dtrees models available through Galaxy extend the workflow:

- {% tool [3Dtrees: ForestMamba](toolshed.g2.bx.psu.edu/repos/bgruening/3dtrees_forestmamba/3dtrees_forestmamba/1.0.0+galaxy0) %}
  offers an alternative for individual-tree segmentation and can also separate leaf and wood points.
- {% tool [3Dtrees: DetailView](toolshed.g2.bx.psu.edu/repos/bgruening/3dtrees_detailview/3dtrees_detailview/1.1.1+galaxy0) %}
  can predict tree species from segmented trees, within the classes supported by its trained model.

These optional steps could be covered in future Galaxy tutorials.

> <details-title>Working with multiple point-cloud files</details-title>
>
> The appropriate workflow depends on what the files represent:
>
> - Process independent scans or flight areas separately.
> - Files that are true tiles of one survey may be processed together after checking their coordinate
>   systems, shared boundaries, and attributes.
> - Overlapping drone flights need their own alignment and overlap checks. SmartTile only handles
>   overlap between the tiles it creates.
>
> The Pic Saint-Loup files are four flight areas, not ready-made tiles. A and B are separate, while
> C and D overlap. This tutorial therefore uses only file D.
>
{: .details}

# Conclusion

You turned one drone LiDAR file into a point cloud with predicted tree labels and an interactive
view. Along the way, you checked the input, reduced processing work, predicted trees, and prepared
the result for inspection.

`PredSemantic_SAT` tells you the predicted class; `PredInstance_SAT` tells you which predicted tree a point
belongs to. These labels provide a starting point for ecological analysis. Your visual checks can
reveal obvious errors; measuring accuracy requires independently labelled reference trees.
