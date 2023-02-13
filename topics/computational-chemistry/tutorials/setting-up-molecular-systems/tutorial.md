---
layout: tutorial_hands_on

title: Setting up molecular systems
zenodo_link: 'https://zenodo.org/record/2600690'
level: Intermediate
questions:
- How to get started modelling a protein and a ligand?
objectives:
- learn about the Protein Data Bank
- learn how to set up up a model protein and ligand system (with CHARMM-GUI)
- learn how to upload the system to Galaxy
requirements:
  -
    type: "internal"
    topic_name: computational-chemistry
    tutorials:
      - setting-up-molecular-systems
time_estimation: 2H
key_points:
  - "The PDB is a key resource for finding protein structures."
  - "Using CHARMM-GUI is one way to prepare a protein and ligand system."
  - "To get data into Galaxy you can upload a file from your computer or paste in a web address."
follow_up_training:
  -
    type: "internal"
    topic_name: computational-chemistry
    tutorials:
      - md-simulation-namd
contributors:
  - chrisbarnettster
  - simonbray
  - nagoue

---

> <comment-title>Audience</comment-title>
> This tutorial is intended for those who are new to the computational chemistry tools in Galaxy.
{: .comment}

# Introduction


In this tutorial, we'll cover the basics of molecular modelling by setting up a protein in complex with a ligand and uploading the structure to Galaxy. This tutorial will make use of CHARMM-GUI. Please note that the follow-up to this tutorial (located [here]({% link topics/computational-chemistry/tutorials/md-simulation-namd/tutorial.md %})) requires access to NAMD Galaxy tools, which can be accessed using the [Docker container](https://github.com/scientificomputing/BRIDGE) but are currently not available on any public Galaxy server.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Cellulase and cellulose

To start we'll look at the PDB and find the entry for a fungal enzyme that cleaves cellulose. The enzyme is 7CEL, a hydrolase as seen in [the figure.](#figure-1)

![Snapshot of 7CEL pdb with octaose ligand]({% link topics/computational-chemistry/images/enzyme.jpg %} "7CEL Cellulase with a short chain cellulose (octaose) ligand")

In this section we'll access the PDB, download the correct structure, import it and view in Galaxy.

> <details-title>Background: What is the PDB (Protein Data Bank) and format?</details-title>
>
> The Protein Data Bank (PDB) format contains atomic coordinates of biomolecules and provides a standard representation for macromolecular structure data derived from X-ray diffraction and NMR studies. Each structure is stored under a four-letter accession code. For example, the PDB file we will use is assigned the code [7CEL](https://www.rcsb.org/pdb/explore/explore.do?structureId=7CEL)).
>
> More resources:
>
>  -  Multiple structures are stored and can be queried at [https://www.rcsb.org/](https://www.rcsb.org/)
>  - Documentation describing the PDB file format is available from the wwPDB at [http://www.wwpdb.org/documentation/file-format.php](http://www.wwpdb.org/documentation/file-format.php).
{: .details}


> <details-title>Background: Why choose a cellulase?</details-title>
>
> Using enzymes to break down abundant cellulose into disaccharide units (cellobiose) is a method to optimise the
> biofuel process. {% cite barnett11 %}
>
> More resources:
>
  - [https://en.wikipedia.org/wiki/Cellulase](https://en.wikipedia.org/wiki/Cellulase)
  - [https://en.wikipedia.org/wiki/Biofuel](https://en.wikipedia.org/wiki/Biofuel)
  - [Fungal Cellulases](https://pubs.acs.org/doi/full/10.1021/cr500351c)
  - [Cellobiohydrolase I Induced Conformational Stability and Glycosidic Bond Polarization ](https://pubs.acs.org/doi/10.1021/ja103766w)
{: .details}

## Get data

The 7CEL [PDB](https://files.rcsb.org/download/7CEL.pdb) does not include a complete 8 unit substrate and some modelling is required. The correctly modelled substrate is provided for this tutorial.
> <details-title>More details about the modelling done</details-title>
> - VMD (visualisation software) was used for atomic placement and CHARMM was used for energy minimisation.
> - The PDB structure contains a mutation at position 217 (glumatate to glutamine). Our structure reverses this.
> - The ligand was modelled separately and inserted into the binding site.
>
{: .details}

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the files from the Zenodo link provided.
>
>    ```
>    https://zenodo.org/record/2600690/files/7cel_modeled.pdb?download=1
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets.
> 4. Check that the datatype is correct. The file should have the PDB datatype.
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="pdb" %}
>
{: .hands_on}

# Modelling with CHARMM-GUI
It is convenient to set up the molecular system outside Galaxy using a tool such as CHARMM-GUI. Alternative methods are possible - see the [GROMACS tutorial]({% link topics/computational-chemistry/tutorials/md-simulation-gromacs/tutorial.md %}) for an example. {% cite jo17 %}

> <tip-title>Viewing figures</tip-title>
> * Some of the figures are screenshots and it may be difficult to make out details
> * Right-click on the image and choose 'Open image in new tab' to view
> * Zoom in and out as needed to see the content
{: .tip}

![CHARMM-GUI interface]({% link topics/computational-chemistry/images/charmmgui.png %} "The CHARMM-GUI interface")

Go to the correct section depending on which MD engine you will be using.

## CHARMM

### Upload the PDB to CHARMM-GUI
[Navigate to CHARMM-GUI](http://www.charmm-gui.org/?doc=input/pdbreader) and use the Input Generator, specifically the PDB Reader tool and upload the Cellulase PDB file. Press 'Next Step: Select Model/Chain' in the bottom right corner.

> <hands-on-title>Upload the PDB to CHARMM-GUI</hands-on-title>
>
> 1. Retrieve the modelled PDB structure from [Zenodo](https://doi.org/10.5281/zenodo.2600690).
> 2. Upload the PDB and choose CHARMM format.
> ![Snapshot of CHARMM-GUI PDB reader section]({% link topics/computational-chemistry/images/charmmgui-reader.png %} "The CHARMM-GUI PDB Reader tool")
{: .hands_on}

### Select both protein and ligand models

> <hands-on-title>Generate PDB file</hands-on-title>
> Two model chains are presented for selection: the protein (PROA) and the hetero residue, which is the ligand or glycan in this case (HETA). Select both, and press 'Next Step: Generate PDB' in the bottom right corner.
> ![Snapshot of CHARMM-GUI model section]({% link topics/computational-chemistry/images/charmmgui-modelchain.png %} "Select both ligand and protein models in CHARMM-GUI")
{: .hands_on}

### Manipulate the system
> <hands-on-title>Make necessary modifications</hands-on-title>
> Rename the hetero chain to BGLC and add ten disulfide bonds to the protein, as shown in the figure. Then press 'Next Step: Generate PDB' in the bottom right corner.
> ![Snapshot of CHARMM-GUI renaming section]({% link topics/computational-chemistry/images/charmmgui-manipulate.png %} "Rename the chains in CHARMM-GUI")
{: .hands_on}

### Download the output
> <hands-on-title>Download CHARMM output</hands-on-title>
> The output is a .tgz file (a tarball or zipped tarball). Inside the archive you will see all inputs and outputs from CHARMM-GUI.
> ![Snapshot of CHARMM-GUI CHARMM output section]({% link topics/computational-chemistry/images/charmmgui-charmmoutput.png %} "CHARMM output from CHARMM-GUI")
{: .hands_on}

> <tip-title>What is a .tgz file?</tip-title>
>
> This is a compressed file which contains all the output files created by the CHARMM-GUI. To access them, the .tgz file needs to be decompressed. There should be a tool available on your operating system for this. If you prefer to use the command line, tar will work fine on Linux or Mac `tar -zxvf example.tgz`. On Windows use [7zip](https://www.7-zip.org/download.html), or download Git for windows and use Git Bash.
{: .tip}


### Upload to Galaxy
> <hands-on-title>Upload files to Galaxy</hands-on-title>
> Upload the step1_pdbreader.psf and step1_pdbreader.crd files to your Galaxy instance and run the system setup tool.
{: .hands_on}

## NAMD

### Upload the PDB to CHARMM-GUI

> <hands-on-title>Upload the PDB to CHARMM-GUI</hands-on-title>
> Retrieve the modelled PDB structure from [Zenodo](https://doi.org/10.5281/zenodo.2600690).
> [Navigate to CHARMM-GUI](http://www.charmm-gui.org/?doc=input/mdsetup) and use the Input Generator, specifically the Solution Builder tool. Upload the PDB file, selecting 'CHARMM' as the file format. Press 'Next Step: Select Model/Chain' in the bottom right corner.
> ![Snapshot of CHARMM-GUI Solution Builder tool ]({% link topics/computational-chemistry/images/charmmgui-mdsimulator-solution-builder.png %} "The CHARMM-GUI Solution Builder tool")
{: .hands_on}


### Select both protein and ligand models
> <hands-on-title>Generate PDB file</hands-on-title>
> Two model chains are presented for selection: the protein (PROA) and the hetero residue, which is the ligand or glycan in this case (HETA). Select both, and press 'Next Step: Generate PDB' in the bottom right corner.
> ![Snapshot of CHARMM-GUI model section]({% link topics/computational-chemistry/images/charmmgui-modelchain.png %} "Select both ligand and protein models in CHARMM-GUI")
{: .hands_on}


### Manipulate the system
> <hands-on-title>Make necessary modifications</hands-on-title>
> Rename the hetero chain to BGLC and add disulfide bonds.Press 'Next Step: Generate PDB' in the bottom right corner.
> ![Snapshot of CHARMM-GUI renaming section]({% link topics/computational-chemistry/images/charmmgui-manipulate.png %} "Rename the chains in CHARMM-GUI")
{: .hands_on}


### Set up the waterbox and add ions
> <hands-on-title>Solvate the protein</hands-on-title>
> Set up a waterbox. Use a size of 10 angstroms and choose a cubic box ('rectangular' option).
> ![Snapshot of CHARMM-GUI waterbox section]({% link topics/computational-chemistry/images/charmm-gui-waterbox-from-solution-builder.png %} "Setting up a waterbox in CHARMM-GUI")
{: .hands_on}

> <question-title></question-title>
>
> Why is 10 angstrom a fair choice for the buffer?  Why choose 0.15M NaCl?
>
> > <solution-title></solution-title>
> > Under periodic boundary conditions, we need to ensure the protein can never interact with its periodic image, otherwise artefacts are introduced. Allowing 10 angstroms between the protein and the box edge ensures the two images will always be at minimum 20 angstroms apart, which is sufficient.
> >
> > Some of the residues on the protein surface are charged and counter-ions need to be present nearby to neutralise them. Failure to explicitly model salt ions may destabilise the protein.
> {: .solution}
{: .question}


### Generate the FFT automatically

> <hands-on-title>Generate the FFT</hands-on-title>
> Particle Mesh Ewald (PME) summation is the method being used to calculate long-range interactions in this system. To improve the computational time a Fast Fourier Transform (FFT) is used. A detailed discussion of FFT will not be presented here; there are many articles on the subject. Try [Wikipedia](https://en.wikipedia.org/wiki/Ewald_summation) and [Ewald summation techniques in perspective: a survey](https://doi.org/10.1016/0010-4655(96)00016-1).
> ![Snapshot of CHARMM-GUI FFT section]({% link topics/computational-chemistry/images/charmmgui-fft.png %} "Setting up a FFT in CHARMM-GUI")
{: .hands_on}

### Download the output
> <hands-on-title>Solvate the protein</hands-on-title>
> The output is a .tgz file (a tarball or zipped tarball). Inside the archive you will see all inputs and outputs from CHARMM-GUI.
> ![Snapshot of CHARMM-GUI NAMD output section]({% link topics/computational-chemistry/images/charmmgui-namdoutput.png %} "NAMD output from CHARMM-GUI")
{: .hands_on}

> <tip-title>What is a .tgz file?</tip-title>
>
> This is a compressed file and needs to be uncompressed using the correct tool.
> On Linux or Mac: tar will work fine `tar -zxvf example.tgz`.
> On Windows use [7zip](https://www.7-zip.org/download.html) or download Git for windows and use Git Bash.
{: .tip}


### Upload to Galaxy
> <hands-on-title>Upload files to Galaxy</hands-on-title>
Upload the following files to your Galaxy instance and ensure the correct datatype is selected:
 - step3_pbcsetup.psf -> xplor psf input (psf format)
 - step3_pbcsetup.pdb -> pdb input (pdb format)
 - Checkfft.str -> PME grid specs (txt format)
- step2.1_waterbox.prm -> waterbox prm input (txt format)
{: .hands_on}

You are now ready to run the NAMD workflow, which is discussed in another [tutorial]({% link topics/computational-chemistry/tutorials/md-simulation-namd/tutorial.md %}).


# Conclusion


{% icon trophy %} Well done! You have started modelling a cellulase protein and uploaded it into Galaxy. The next step is running molecular dynamics simulations ([tutorial]({% link topics/computational-chemistry/tutorials/md-simulation-namd/tutorial.md %}))
