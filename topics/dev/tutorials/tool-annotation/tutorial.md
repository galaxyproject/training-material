---
layout: tutorial_hands_on
title: Annotation of a tool with EDAM ontology terms using bio.tools
level: Introductory
subtopic: tooldev
questions:
- How are Galaxy tools linked to EDAM ontology?
- How to connect Galaxy tools to bio.tools?
objectives:
- Identify Galaxy tools without bio.tools entry
- Create a bio.tools entry
- Update a bio.tools entry 
- Add EDAM ontology terms to a bio.tools entry
- Link a Galaxy tool to its corresponding bio.tools entry
time_estimation: 1H
key_points:
- Galaxy tools can get EDAM ontology terms from bio.tools
- bio.tools entry can be created and modified to provide the best EDAM annotations
- bio.tools entry can easily be added to a Galaxy tool
contributions:
  authorship:
    - bebatut

---

Galaxy offers thousands of tools. A vast majority of these tools come with poor metadata. 

This prevents filtering for all tools in a specific research community or domain, and makes it all but impossible to employ advanced filtering with ontology terms like ones from EDAM or to group tools given ontology to improve the Galaxy tool panel.

[EDAM](https://edamontology.org/page) ({% cite black2021edam %}) is a comprehensive ontology of well-established, familiar concepts that are prevalent within bioscientific data analysis and data management. It includes 4 main sections of concepts (sub-ontologies):

- **Topic**:  A category denoting a rather broad domain or field of interest, of study, application, work, data, or technology. Topics have no clearly defined borders between each other
- **Operation**> A function that processes a set of inputs and results in a set of outputs, or associates arguments (inputs) with values (outputs)
- **Data**: Information, represented in an information artefact (data record) that is "understandable" by dedicated computational tools that can use the data as input or produce it as output
- **Format**: A defined way or layout of representing and structuring data in a computer file, blob, string, message, or elsewhere.

![Simplified data flow diagram in EDAM architecture: boxes for concepts, lines for relations. Streamlined data management.](./images/EDAMrelations.png "EDAM architecture is simple. Boxes indicate top-level concepts (sections, sub-ontologies), and lines indicate types of relations. Source: <a href="https://edamontology.org/page">EDAM website</a>")

The ontology can be navigated using [EDAM browser](https://edamontology.github.io/edam-browser/):

<iframe id="edam" src="https://edamontology.github.io/edam-browser/#operation_0291" frameBorder="0" width="80%" height="600px"> ![Krona at bacteria level](./images/edam_browser.png) </iframe>

A tool or software can be then characterized by different EDAM terms:
- A topic term, *e.g.*  "Proteomics",
- An operation (a specific scientific thing that a tool does) term, *e.g.* "Peptide identification",
- A data term for the type of biological data, *e.g.* "Mass spectrometry spectra",
- A format term, *e.g.* "Thermo RAW".

The annotation of tools can be done on [bio.tools](https://bio.tools/). bio.tools ({% cite ison2016tools %}) is a global portal for bioinformatics resources that helps researchers to find, understand, compare, and select resources suitable for their work. It relies on the EDAM ontology for standardazing the annotations.

In Galaxy, tools can be annotated with EDAM terms, either by putting them directly in the wrapper or extracting them from their corresponding bio.tools entry by linking to it. The advantage of the second approach is that only one place (bio.tools entry) needs to store and keep updated EDAM terms.

The aim is this tutorial is to improve the annotation of a Galaxy tool by linking it to a bio.tools after its creation or update with proper EDAM terms.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Choose a tool without bio.tool identifier

To start, we need to select a tool without bio.tool identifier.

> <hands-on-title>Choose a tool without bio.tool identifier</hands-on-title>
>
> 1. Open [the list of Galaxy tool](https://galaxyproject.github.io/galaxy_tool_metadata_extractor/)
> 2. Click on **Add Condition**
> 3. Select *bio.tool id* in **Data** drop-down
> 4. Select *Empty* in **Condition** drop-down
> 5. Select a tool in the list
>
{: .hands_on}

# Identify the tool bio.tools entry

Let's now search for our selected tool on bio.tools.

> <hands-on-title>Search a tool in bio.tools</hands-on-title>
>
> 1. Open [bio.tools](https://bio.tools/)
> 2. Type the name of your tool in the "Search bio.tools" bar on the top
>
{: .hands_on}

{% include _includes/cyoa-choices.html option1="No existing entry" option2="Existing entry" default="No bio.tool entry" text="Have you found the tool in bio.tools?" disambiguation="biotool"%}

<div class="No-existing-entry" markdown="1">

# Create a bio.tools entry for a tool

When the tool is not on bio.tools, we need to create an entry for it and populate it with metadata

> <hands-on-title> Create a bio.tools entry with minimum metadata </hands-on-title>
>
> 1. Sign up for bio.tools
> 2. Select **Add a tool** from the drop-down **Menu**
> 3. Fill in general information
>    1. Fill in **Tool name** 
>    2. Fill in **Description**
>    3. Fill in **Homepage URL**
> 4. Add EDAM operation term and EDAM data term for inputs and outputs
>    
>    {% include topics/dev/tutorials/tool-annotation/add_edam_function.md %}
> 
> 5. Add EDAM topic term, tool type, language and others
>    
>    {% include topics/dev/tutorials/tool-annotation/add_edam_topic.md %}
> 
> 6. Click on **Validate** on the top
> 7. Click on **Save** to create the bio.tools entry
> 8. Copy the bio.tools id
{: .hands_on}

</div>

<div class="Existing-entry" markdown="1">

# Review and update the EDAM terms for a bio.tools entry

Before linking the Galaxy tools with its corresponding bio.tools entry, we need to check if the tool is correctly annotated with EDAM terms.

> <hands-on-title>Check EDAM terms in a bio.tools entry</hands-on-title>
>
> 1. Open the bio.tool entry for the tool
> 2. Check the EDAM topic terms by looking at the green boxes (if they exist) below the tool name, URL, and available versions
> 3. Check the EDAM operation terms by looking at the blue boxes (if they exist) below the tool description
{: .hands_on}

{% include _includes/cyoa-choices.html option1="Terms to be modified" option2="Correct terms" default="Terms to be modified" text="What do you think about the EDAM terms in the bio.tools entry?" disambiguation="edamUpdate"%}

<div class="Terms-to-be-modified" markdown="1">

To modify EDAM terms in a bio.tools entry, we need to request editing rights and then modify this entry.

> <hands-on-title>Modify EDAM terms in a bio.tools entry</hands-on-title>
>
> 1. Sign up for bio.tools
> 2. Click on **Request editing rights** on the bottom of bio.tool entry page
> 3. Wait for the request to be approved
> 4. Click on **Update this record**
> 5. Update or add EDAM operation term and EDAM data term for inputs and outputs
>    
>    {% include topics/dev/tutorials/tool-annotation/add_edam_function.md %}
> 
> 6. Update or add EDAM topic term, tool type, language and others
>    
>    {% include topics/dev/tutorials/tool-annotation/add_edam_topic.md %}
> 7. Click on **Validate** on the top
> 8. Click on **Save** to create the bio.tools entry
> 9. Copy the bio.tools id
{: .hands_on}

</div>

</div>

# Link the Galaxy tool and a bio.tools entry

To link the Galaxy tool to its corresponding bio.tools entry, we need to first find the source of the wrapper.

> <hands-on-title>Find the Galaxy wrapper</hands-on-title>
>
> 1. Go to the tool on any Galaxy server
> 2. Click on the drop-down menu next to the **Run tool** button
> 3. Select **See in Tool Shed**
> 4. Click on the link to the development repository
> 5. Fork the repository
{: .hands_on}

Once we have the wrapper, we can add the bio.tools entry.

> <hands-on-title>Add bio.tools entry to the Galaxy wrapper</hands-on-title>
> 1. Open the Galaxy tool XML file
> 2. Add the xref snippet below `macros` and before `requirements`
>
>    ```
>    <xrefs>
>        <xref type="bio.tools">biotool-id</xref>
>    </xrefs>
>    ```
>
> 3. Replace `biotool-id` above by the bio.tool id of your tool
> 3. Commit the change on a new branch
> 4. Make a pull request (PR) against the original repository
> 5. Wait patiently for the PR to be merged, at which point the new bio.tools reference will be added to the Galaxy tool wrapper. Make sure to respond to any feedback from the owner of the wrapper.
>
{: .hands_on}

# Conclusion