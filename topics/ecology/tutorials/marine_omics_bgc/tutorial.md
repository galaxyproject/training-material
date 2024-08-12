---
layout: tutorial_hands_on

title: Marine Omics identifying biosynthetic gene clusters
zenodo_link: ''
questions:
- Which biological questions are addressed by the tutorial?
- Which bioinformatics techniques are important to know for this type of data?
objectives:
- Follow a marine omics analysis
- Learn to conduct a secondary metabolite biosynthetic gene cluster (SMBGC) Annotation
- Discover Sanntis a tool for identifying BGCs in genomic & metagenomic data
- Manage fasta files
time_estimation: 3H
key_points:
- SMBGC annotation
- Marine omics analysis
tags:
  - earth-system
  - ocean
  - marine omics
contributions:
  authorship:
    - Marie59
  funding:
    - fairease
    - eurosciencegateway

---

# Introduction
Through this tutorial, you will learn in the first part how to produce a protein fasta file from a nucleotide fasta file using Prodigal (Fast, reliable protein-coding gene prediction for prokaryotic genomes).

Then, you'll be using InterProscan to create a tabular. Interproscan is a batch tool to query the InterPro database. It provides annotations based on multiple searches of profile and other functional databases.

And finally, you will discover Sanntis both to build genbank and especially to conduct annotation using Neural Networks Trained on Interpro Signatures.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get data

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial and give it a name (for example “Marine Omics: SMBGC annotation”) for you to find it again later if needed.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the files
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>
> 3. Rename the datasets
>
>
>
{: .hands_on}

{% include _includes/cyoa-choices.html option1="Tools" option2="Workflow"
       text="Do you want to run the workflow or to discover the tools one by one ?" %}
       
<div class="Workflow" markdown="1">

# Import and launch the workflow
> <hands-on-title>Import the workflow</hands-on-title>
>   - Click on Workflow on the top menu bar of Galaxy. You will see a list of all your workflows.
>   - Option 1: use the URL
>   	- Click on {% icon galaxy-upload %} Import at the top-right of the screen
>   	- Paste the URL of the workflow into the box labelled “Archived Workflow URL” `https://earth-system.usegalaxy.eu/u/marie.josse/w/marine-omics-identifying-biosynthetic-gene-clusters`
>   - Option 2: use the workflow name
>   	- Click on **Public workflows** at the top-right of the screen
>       - Search for **Marine Omics identifying biosynthetic gene clusters** 
>       - In the workflow preview box click on {% icon galaxy-upload %} Import
>   - Click the Import workflow button
{: .hands_on}

> <hands-on-title>Run the workflow</hands-on-title>
>    - Click on Workflow on the top menu bar of Galaxy. You will see a list of all your workflows.
>    - Click on the {% icon workflow-run %} (Run workflow) button next to your workflow
>    - /!\ Select **Yes** for **Workflow semi automatic** 
>    - Configure the workflow as needed with the 2 datasets you uploaded right before (**BGC0001472.fna**)
>    - Click the Run Workflow button at the top-right of the screen
>    - You may have to refresh your history to see the queued jobs
{: .hands_on}

Now you don't have to do anything else. You should see all the different steps of the worklow appear in your history.
When the workflow is fully runned you should have the folowing history.

![Image of the history with all the steps of the workflow](../../images/marineomics/history.png "History")

> <tip-title>Extract a RO-Crate</tip-title>
> Workflows are a powerful Galaxy feature that allows you to scale up your analysis by performing an end-to-end analysis with a single click of a button. In order to keep provenance of the workflow invocation (an invocation of a workflow means one run or execution of the workflow) it can be exported from Galaxy in the form of a Workflow Run Crate RO-Crate profile.
> 
>
> After the workflow has completed, we can export the RO-Crate. The crate does not appear in your history, but can be accessed from the User -> Workflow Invocations menu on the top bar.
> 
> > <hands-on-title> Extract a RO-Crate </hands-on-title>
> >
> > 1. In the top right of your history, go to {% icon galaxy-history-options %} -> Show Invocations
> >
> > ![Image of the history options](../../images/marineomics/showinvoc.png "History options")
> >
> > Our latest workflow run should be listed at the top.
> >
> > 2. Click on it to expand it:
> >
> > ![Image of the workflow invocation](../../images/marineomics/workflowinvoc.png "Workflow invocation")
> > 
> > 3. Click on the Export tab in the expanded view of the workflow invocation.
> > 	You should see a page that contains three download options:
> >		- Research Object Crate (RO-Crate) 
> >		- BioCompute Object
> >		- File
> > 4. Click on the **Generate** {% icon galaxy-download %} option of the RO-Crate box (1st box)
> > 
> > ![Image of the RO-Crate download](../../images/marineomics/rocrate.png "RO-Crate")
> {: .hands-on}
> 
> 
> Great work! You have created a Workflow Run Crate. This makes it easy to track the provenance of the executed workflow.
{: .tip}

</div>


<div class="Tools" markdown="1">
# Prodigal Gene Predictor: generate a protein fasta file



## Prodigal Gene Predictor

> <hands-on-title> Run prodigal </hands-on-title>
>
> 1. {% tool [Prodigal Gene Predictor](toolshed.g2.bx.psu.edu/repos/iuc/prodigal/prodigal/2.6.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Specify input file"*: `BGC0001472.fna` (Input the nucleotide fasta file)
>    - *"Specify mode"*: `Meta : Anonymous sequences, analyze using preset training files, ideal for metagenomic data or single short sequences`
>
>    You don't need to change any other parameters leave them on the default input.
>
> 2. Click on **Run Tool**
>
{: .hands_on}

You should have 4 new outputs appearing in your history. In these outputs you should have **Prodigal Gene Predictor on data 1 : protein translations file**.
You can click on it and then click on the {% icon galaxy-eye %} (eye).

![Image of the protein fasta file from prodigal and of the new history with all the new outputs](../../images/marineomics/prodigal.png "Prodigal outputs")

You can notice here that at each end of the sequence there's a *. Later on we will need to remove this star. But, first we are going to use this protein file to build the Genbank that Sanntis need to make a SMBGC annotation. 

## Sanntis for building a Genbank file

> <hands-on-title> Build Genbank </hands-on-title>
>
> 1. {% tool [Sanntis biosynthetic gene clusters](toolshed.g2.bx.psu.edu/repos/ecology/sanntis_marine/sanntis_marine/0.9.3.5+galaxy1) %} with the following parameters:
>    - *"Do you want to build a genbank or to make a SMBGC Annotation?"*: `Build genbank`
>        - {% icon param-file %} *"Input a nucleotide fasta file"*: `BGC0001472.fna` (Input the nucleotide fasta file)
>        - {% icon param-file %} *"Input a protein fasta file"*: `Prodigal Gene Predictor on data 1 : protein translations file` (output of **Prodigal Gene Predictor** {% icon tool %})
> 2. Click on **Run Tool**
>
{: .hands_on}

![Image of the Genbank file produces by Sanntis](../../images/marineomics/genbank.png "Genbank file")

## Regex Find And Replace
Remember earlier we noticed the star * in the protein fasta file ? 

Now is the time to remove it ! This is to simplify the next step with InterProScan.

> <hands-on-title> Remove * </hands-on-title>
>
> 1. {% tool [Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regex1/1.0.3) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `Prodigal Gene Predictor on data 1 : protein translations file` (output of **Prodigal Gene Predictor** {% icon tool %})
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `\*$`
>            - *"Replacement"*: `` (leave an empty box there)
>
> 2. Click on **Run Tool**
>
{: .hands_on}

Check if the * were well removed. 

# InterProScan

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [InterProScan](toolshed.g2.bx.psu.edu/repos/bgruening/interproscan/interproscan/5.59-91.0+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Protein FASTA File"*: `Regex Find And Replace on data **` (output of **Regex Find And Replace** {% icon tool %})
>    - *"Use applications with restricted license, only for non-commercial use?"*: `No`
> You can leave all the other parameters on the default input.
>
> 2. Click on **Run Tool**
{: .hands_on}

![Image of InterProScan output in tabular (tsv) in the history](../../images/marineomics/interproscan.png "InterProScan output")

# Sanntis for annotating biosynthetic gene clusters

> <hands-on-title> Identify biosynthetic gene clusters </hands-on-title>
>
> 1. {% tool [Sanntis biosynthetic gene clusters](toolshed.g2.bx.psu.edu/repos/ecology/sanntis_marine/sanntis_marine/0.9.3.5+galaxy1) %} with the following parameters:
>    - *"Do you want to build a genbank or to make a SMBGC Annotation?"*: `Run sanntis`
>        - {% icon param-file %} *"Input the tabular file from InterProScan"*: `InterProScan on data **` (output of **InterProScan** {% icon tool %})
>        - {% icon param-file %} *"Input a Genbank file"*: `Sanntis output data genbank` (output of **Sanntis biosynthetic gene clusters** {% icon tool %})
>
> 2. Click on **Run Tool**
{: .hands_on}

</div>

Finally, you should have one gff3 file in your history under **Sanntis output data**
![Image of Sanntis output](../../images/marineomics/sanntis.png "Sanntis output")


# Conclusion
Here you now know how to conduct a Marine Omics analysis

# Extra information

Coming up soon even more tutorials on and other Earth-System related trainings. Keep an {% icon galaxy-eye %} open if you are interested!
