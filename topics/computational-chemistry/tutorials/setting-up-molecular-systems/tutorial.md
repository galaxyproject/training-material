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

---

> ### {% icon comment %} Audience
> This tutorial is intended for those who are new to the computational chemistry tools in Galaxy.
{: .comment}

# Introduction
{:.no_toc}

In this tutorial, we'll cover the basics of molecular modelling by setting up a protein in complex with a ligand and uploading the structure to Galaxy. This tutorial will make use of CHARMM-GUI. Please note that the follow-up to this tutorial (located [here]({% link topics/computational-chemistry/tutorials/md-simulation-namd/tutorial.md %})) requires access to NAMD Galaxy tools, which can be accessed using the [Docker container](https://github.com/scientificomputing/BRIDGE) but are currently not available on any public Galaxy server.

> ### Agenda
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

> ### {% icon details %} Background: What is the PDB (Protein Data Bank) and format?
>
> The Protein Data Bank (PDB) format contains atomic coordinates of biomolecules and provides a standard representation for macromolecular structure data derived from X-ray diffraction and NMR studies. Each structure is stored under a four-letter accession code. For example, the PDB file we will use is assigned the code [7CEL](https://www.rcsb.org/pdb/explore/explore.do?structureId=7CEL)).
>
> More resources:
>
>  -  Multiple structures are stored and can be queried at [https://www.rcsb.org/](https://www.rcsb.org/)
>  - Documentation describing the PDB file format is available from the wwPDB at [http://www.wwpdb.org/documentation/file-format.php](http://www.wwpdb.org/documentation/file-format.php).
{: .details}


> ### {% icon details %} Background: Why choose a cellulase?
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
> ### {% icon details %} More details about the modelling done
> - VMD (visualisation software) was used for atomic placement and CHARMM was used for energy minimisation.
> - The PDB structure contains a mutation at position 217 (glumatate to glutamine). Our structure reverses this.
> - The ligand was modelled separately and inserted into the binding site.
>
{: .details}

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial.
>
>    {% include snippets/create_new_history.md %}
>
> 2. Import the files from the Zenodo link provided.
>
>    ```
>    https://zenodo.org/record/2600690
>    ```
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Rename the datasets.
> 4. Check that the datatype is correct. The file should have the PDB datatype.
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
{: .hands_on}

# Modelling with CHARMM-GUI
It is convenient to set up the molecular system outside Galaxy using a tool such as CHARMM-GUI. Alternative methods are possible - see the [GROMACS tutorial]({% link topics/computational-chemistry/tutorials/md-simulation-gromacs/tutorial.md %}) for an example. {% cite jo17 %}

> ### {% icon tip %} Tip: Viewing figures
> * Some of the figures are screenshots and it may be difficult to make out details
> * Right-click on the image and choose 'Open image in new tab' to view
> * Zoom in and out as needed to see the content
{: .tip}

![CHARMM-GUI interface]({% link topics/computational-chemistry/images/charmmgui.png %} "The CHARMM-GUI interface")

Go to the correct section depending on which MD engine you will be using.

## CHARMM

### Upload the PDB to CHARMM-GUI
[Navigate to CHARMM-GUI](http://www.charmm-gui.org/?doc=input/pdbreader) and use the Input Generator, specifically the PDB Reader tool and upload the Cellulase PDB file. Press 'Next Step: Select Model/Chain' in the bottom right corner.

> ### {% icon hands_on %} Hands-on: Upload the PDB to CHARMM-GUI
>
> 1. Retrieve the modelled PDB structure from [Zenodo](https://doi.org/10.5281/zenodo.2600690).
> 2. Upload the PDB and choose CHARMM format.
> ![Snapshot of CHARMM-GUI PDB reader section]({% link topics/computational-chemistry/images/charmmgui-reader.png %} "The CHARMM-GUI PDB Reader tool")
{: .hands_on}

### Select both protein and ligand models

> ### {% icon hands_on %} Hands-on: Generate PDB file
> Two model chains are presented for selection: the protein (PROA) and the hetero residue, which is the ligand or glycan in this case (HETA). Select both, and press 'Next Step: Generate PDB' in the bottom right corner.
> ![Snapshot of CHARMM-GUI model section]({% link topics/computational-chemistry/images/charmmgui-modelchain.png %} "Select both ligand and protein models in CHARMM-GUI")
{: .hands_on}

### Manipulate the system
> ### {% icon hands_on %} Hands-on: Make necessary modifications
> Rename the hetero chain to BGLC and add ten disulfide bonds to the protein, as shown in the figure. Then press 'Next Step: Manipulate PDB' in the bottom right corner.
> ![Snapshot of CHARMM-GUI renaming section]({% link topics/computational-chemistry/images/charmmgui-manipulate.png %} "Rename the chains in CHARMM-GUI")
{: .hands_on}

### Download the output
> ### {% icon hands_on %} Hands-on: Download CHARMM output
> The output is a .tgz file (a tarball or zipped tarball). Inside the archive you will see all inputs and outputs from CHARMM-GUI.
> ![Snapshot of CHARMM-GUI CHARMM output section]({% link topics/computational-chemistry/images/charmmgui-charmmoutput.png %} "CHARMM output from CHARMM-GUI")
{: .hands_on}

> ### {% icon tip %} What is a .tgz file?
>
> This is a compressed file which contains all the output files created by the CHARMM-GUI. To access them, the .tgz file needs to be decompressed. There should be a tool available on your operating system for this. If you prefer to use the command line, tar will work fine on Linux or Mac `tar -zxvf example.tgz`. On Windows use [7zip](https://www.7-zip.org/download.html), or download Git for windows and use Git Bash.
{: .tip}


### Upload to Galaxy
> ### {% icon hands_on %} Hands-on: Upload files to Galaxy
> Upload the step1_pdbreader.psf and step1_pdbreader.crd files to your BRIDGE instance and run the system setup tool.
{: .hands_on}

## NAMD

### Upload the PDB to CHARMM-GUI

> ### {% icon hands_on %} Hands-on: Upload the PDB to CHARMM-GUI
> Retrieve the modelled PDB structure from [Zenodo](https://doi.org/10.5281/zenodo.2600690).
[Navigate to CHARMM-GUI](http://www.charmm-gui.org/?doc=input/mdsetup) and use the Input Generator, specifically the Quick MD Simulator tool. Upload the PDB file, selecting 'CHARMM' as the file format. Press 'Next Step: Select Model/Chain' in the bottom right corner.
> ![Snapshot of CHARMM-GUI Quick MD Simulator tool ]({% link topics/computational-chemistry/images/charmmgui-mdsimulator.png %} "The CHARMM-GUI Quick MD Simulator tool")
{: .hands_on}

### Select both protein and ligand models
> ### {% icon hands_on %} Hands-on: Generate PDB file
> Two model chains are presented for selection: the protein (PROA) and the hetero residue, which is the ligand or glycan in this case (HETA). Select both, and press 'Next Step: Generate PDB' in the bottom right corner.
> ![Snapshot of CHARMM-GUI model section]({% link topics/computational-chemistry/images/charmmgui-modelchain.png %} "Select both ligand and protein models in CHARMM-GUI")
{: .hands_on}

### Manipulate the system
> ### {% icon hands_on %} Hands-on: Make necessary modifications
> Rename the hetero chain to BGLC and add disulfide bonds.
> ![Snapshot of CHARMM-GUI renaming section]({% link topics/computational-chemistry/images/charmmgui-manipulate.png %} "Rename the chains in CHARMM-GUI")
{: .hands_on}

### Set up the waterbox and add ions
> ### {% icon hands_on %} Hands-on: Solvate the protein
> Set up a waterbox. Use a size of 10 angstroms and choose a cubic box ('rectangular' option).
> ![Snapshot of CHARMM-GUI waterbox section]({% link topics/computational-chemistry/images/charmmgui-waterbox.png %} "Setting up a waterbox in CHARMM-GUI")
{: .hands_on}

> ### {% icon question %} Question
>
> Why is 10 angstrom a fair choice for the buffer?  Why choose 0.15M NaCl?
>
> > ### {% icon solution %} Solution
> > Under periodic boundary conditions, we need to ensure the protein can never interact with its periodic image, otherwise artefacts are introduced. Allowing 10 angstroms between the protein and the box edge ensures the two images will always be at minimum 20 angstroms apart, which is sufficient.
> >
> > Some of the residues on the protein surface are charged and counter-ions need to be present nearby to neutralise them. Failure to explicitly model salt ions may destabilise the protein.
> {: .solution}
{: .question}


### Generate the FFT automatically

> ### {% icon hands_on %} Hands-on: Generate the FFT
> Particle Mesh Ewald (PME) summation is the method being used to calculate long-range interactions in this system. To improve the computational time a Fast Fourier Transform (FFT) is used. A detailed discussion of FFT will not be presented here; there are many articles on the subject. Try [Wikipedia](https://en.wikipedia.org/wiki/Ewald_summation) and [Ewald summation techniques in perspective: a survey](https://doi.org/10.1016/0010-4655(96)00016-1).
> ![Snapshot of CHARMM-GUI FFT section]({% link topics/computational-chemistry/images/charmmgui-fft.png %} "Setting up a FFT in CHARMM-GUI")
{: .hands_on}

### Download the output
> ### {% icon hands_on %} Hands-on: Solvate the protein
> The output is a .tgz file (a tarball or zipped tarball). Inside the archive you will see all inputs and outputs from CHARMM-GUI.
> ![Snapshot of CHARMM-GUI NAMD output section]({% link topics/computational-chemistry/images/charmmgui-namdoutput.png %} "NAMD output from CHARMM-GUI")
{: .hands_on}

> ### {% icon tip %} What is a .tgz file?
>
> This is a compressed file and needs to be uncompressed using the correct tool.
> On Linux or Mac: tar will work fine `tar -zxvf example.tgz`.
> On Windows use [7zip](https://www.7-zip.org/download.html) or download Git for windows and use Git Bash.
{: .tip}


### Upload to Galaxy
> ### {% icon hands_on %} Hands-on: Upload files to Galaxy
Upload the following files to your BRIDGE instance and ensure the correct datatype is selected:
 - step3_pbcsetup.xplor.ext.psf -> xplor psf input (psf format)
 - step3_pbcsetup.pdb -> pdb input (pdb format)
 - Checkfft.str -> PME grid specs (txt format)
- step2.1_waterbox.prm -> waterbox prm input (txt format)
{: .hands_on}

You are now ready to run the NAMD workflow, which is discussed in another [tutorial]({% link topics/computational-chemistry/tutorials/md-simulation-namd/tutorial.md %}).


# Conclusion
{:.no_toc}

{% icon trophy %} Well done! You have started modelling a cellulase protein and uploaded it into Galaxy. The next step is running molecular dynamics simulations ([tutorial]({% link topics/computational-chemistry/tutorials/md-simulation-namd/tutorial.md %}))
