---
layout: tutorial_hands_on

title: Annotate, prepare tests and publish Galaxy workflows in workflow registries
tags:
- workflows
- FAIR
questions:
- How can a Galaxy workflow be annotated to improve reusability?
- What are the best practices for testing a Galaxy workflow?
- How can a Galaxy workflow be submitted to the Intergalactic Workflow Commission (IWC)?
objectives:
- Annotate a Galaxy workflows with essential metadata
- Annotate and apply best practices to data analysis Galaxy workflows for consistency and reusability
- Implement robust tests to ensure workflow reliability and accuracy
- Successfully integrate key data analysis Galaxy workflows into IWC, improving accessibility and usability
requirements:
- type: internal
  topic_name: galaxy-interface
  tutorials:
      - workflow-editor
time_estimation: 2H
key_points:
- Workflow annotation and best practices improve clarity and reusability
- Standardized metadata enhances workflow documentation
- Testing ensures workflow reliability and accuracy
- IWC submission enables version control and usability by the Galaxy community
contributions:
  authorship:
  - clsiguret
  - nagoue
  - bebatut
  funding:
  - abromics
  - ifb
level: Intermediate
subtopic: workflows
---

Research data is accumulating at an unprecedented rate, presenting significant challenges for achieving fully reproducible science. As a result, implementing high-quality management of scientific data has become a global priority. One key aspect of this effort is the use of computational workflows, which describe the complex, multi-step methods used for data analysis. In Galaxy, workflows are a powerful feature that allows researchers to link multiple steps of complex analyses seamlessly. To maximize their impact, these workflows should adhere to best practices that make them **FAIR: Findable, Accessible, Interoperable, and Reusable**.

The FAIR principles —-Findable, Accessible, Interoperable, and Reusable—- provide practical guidelines for enhancing the value of research data:

- **Findable**: Easy to locate through rich metadata and unique identifiers.
- **Accessible**: Stored in a way that allows them to be retrieved by those who need them.
- **Interoperable**: Usable across different systems and platforms without extensive adaptation.
- **Reusable**: Well-documented and clearly licensed to enable reuse in different contexts.
> <comment-title></comment-title>
> You can learn more about FAIR Data, Workflows, and Research in our [dedicated topic with numerous tutorials]({% link topics/fair/index.md %})
{: .comment}
Applying these principles to workflows is equally important for good data management:

- **Enhanced Discoverability**: Well-annotated and documented workflows are easier to find, making them more likely to be used and cited.
- **Improved Reproducibility**: Standardized and tested workflows ensure that analyses can be reproduced, validating research findings.
- **Community Collaboration**: Sharing workflows through centralized registries fosters collaboration and innovation within the bioinformatics community.
- **Sustainability**: Regular updates and versioning ensure that workflows remain compatible with the latest tools, extending their lifespan and utility.

While making a workflow FAIR might seem complicated ({% cite wilkinson2025applying %}, {% cite goble2020fair %}), publications like "Ten quick tips for building FAIR workflows" ({% cite de2023ten %}) provide practical guidelines to simplify the process:

![This image provides ten quick tips for building FAIR workflows, focusing on findability (registering the workflow and describing it with rich metadata), accessibility (making source code available in a public repository and providing example input data and results), interoperability (adhering to file format standards and making the workflow portable), and reusability (providing a reproducible computational environment, adding a configuration file with defaults, modularizing the workflow, and offering clear documentation).](./images/journal.pcbi.1011369.g001.PNG "Ten quick tips for building FAIR workflows. Source: {% cite de2023ten %}")

Using Galaxy as a workflow management system, many of these tips (all related to **Interoperability** and **Reusability**) are already fulfilled:
- Tip 5 (**Interoperability**): *"The tools integrated in a workflow should adhere to file format standards."*
- Tip 6 (**Interoperability**): *"Make the workflow portable."*
- Tip 7 (**Reusability**): *"Provide a reproducible computational environment to run the workflow."*
- Tip 8 (**Reusability**): *"Add a configuration file with defaults."*
- Tip 9 (**Reusability**): *"Modularize the workflow."*

In this tutorial, we will demonstrate how to fulfill the remaining tips:
- Tip 1 (**Findability**): *"Register the workflow."*
- Tip 2 (**Findability**): *"Describe the workflow with rich metadata."*
- Tip 3 (**Accessibility**): *"Make source code available in a public code repository."*
- Tip 4 (**Accessibility**): *"Provide example input data and results along with the workflow."*
- Tip 10 (**Reusability**): *"Provide clear and concise workflow documentation."*

To illustrate the process, we will use a simple workflow with 2 steps as an example.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Prepare the workflow

Here, we will use a workflow running Falco (FastQC alternative, {% cite de2021falco %}) and MultiQC ({% cite ewels2016multiqc %}) but you can use your own workflow.

> <hands-on-title>Create the workflow into Galaxy</hands-on-title>
>
> 1. Go to the workflow page
> 2. Create a new workflow
>
>    {% snippet faqs/galaxy/workflows_create_new.md %}
>
> 3. Add a single data input
>
> 4. Add {% tool [Falco](toolshed.g2.bx.psu.edu/repos/iuc/falco/falco/1.2.4+galaxy0) %} and {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.27+galaxy3) %}
>
> 5. Connect the input to **Falco**
> 6. Define **MultiQC** parameter
>    - **Which tool was used generate logs?**: `FastQC`
> 
> 7. Connect `text_file` **Falco** to **MultiQC** input
{: .hands_on}

# Annotate the Galaxy workflow with essential metadata

The first step to FAIRify a Galaxy workflow is to fulfill Tip 2 (**Findability**): *"Describe the workflow with rich metadata"*. Describing the workflow with rich metadata helps both humans and machines understand its purpose, facilitating discovery by search engines. Galaxy allows to associate workflows with metadata directly in its interface.

> <hands-on-title>Open workflow attributes</hands-on-title>
>
> 1. Click on **Attributes** on the left side of the workflow editor (or right side for older Galaxy version for clicking on {% icon galaxy-pencil %} **Edit Attributes**)
{: .hands_on}

Some workflow metadata should appear.

> <question-title></question-title>
>
> 1. What is the name of the workflow?
>
> > <solution-title></solution-title>
> > 1. "Unnamed Workflow"
> {: .solution}
{: .question}

Metadata should include details about the workflow's components and clearly outline its purpose, scope, and limitations. This ensures the workflow is easily findable and understandable. 

> <question-title></question-title>
>
> Which metadata is supported in Galaxy workflow interface?
>
> > <solution-title></solution-title>
> > Galaxy workflow interface supports:
> > - Name
> > - Version
> > - Annotation: These notes will be visible when this workflow is viewed.
> > - License
> > - Creator: Either a person or an organization. 
> > - Tags: Apply tags to make it easy to search for and find items with the same tag. 
> {: .solution}
{: .question}

Galaxy workflow interface supports some metadata. Are they enough to fulfill Tip 2 and also Tip 1 *"Register the workflow."*? WorkflowHub, a workflow registry we will explain more later, supports the [following metadata](https://about.workflowhub.eu/docs/metadata-list/):

Name | Description | Mandatory
--- | --- |
**Title** | This field is mandatory and is with some workflow types pre-filled with the title of the workflow. | Yes
**Description** | If a CWL (abstract) file is given, the description will be parsed automatically out of the doc attribute. In any other case this field can be used to write some documentation that will be shown on the workflow page. | No
**Source** | If the workflow came from an external repository (i.e. GitHub), you can include its original URL here. | No
**Maturity** | This field can be used to specify in which maturity state the workflow is. The two available options are: `work-in-progress`, `stable` | No
**Teams** | Every workflow registration is linked to one or more teams. | Yes
**Licence** | The standard licence is Apache Software Licence 2.0. If you did not make the workflow yourself, be sure that the licence corresponds to the licence where you took the workflow from (for example github licences). | No
**Sharing** | Specify who can view the summary, get access to the content, and edit the Workflow. This is possibly already filled in according to the selected project. | No
**Tags** | Choose an appropriate tag for your workflow. Please check if your tag is already available and use the existing one if so. If you make a new tag, keep it simple without capitals or spaces. For example all new covid-19 workflows need to be tagged with covid-19. | No
**Creators** | This is an important section where all the people that were involved in making / publishing this workflow are listed. | No

> <question-title></question-title>
>
> 1. Which WorkflowHub metadata is supported in Galaxy workflow interface?
> 2. Which values should we put for our workflow?
>
> > <solution-title></solution-title>
> >
> > Workflow metadata | Galaxy metadata | Mapping workflow value
> > --- | --- | ---
> > Title | Name | Quality control and mapping
> > Description | Annotation | Workflow runs quality control and mapping of paired-end short-reads data.
> > Source | *Not supported in Galaxy workflow interface* | 
> > Teams | *Not supported in Galaxy workflow interface* | 
> > License | License | MIT License
> > Sharing | *Supported in Galaxy workflow interface but not in the same way as WorkflowHub* | 
> > Tags | Tags | `sequence-analysis`, `mapping`
> > Creators | Creator | 
> >
> {: .solution}
{: .question}

Let's annotate the workflow

> <hands-on-title>Add metadata to the workflow</hands-on-title>
>
> 1. Add a proper title in **Name**
> 2. Add a workflow description in **Annotation**
> 3. Select a license in **License**
>
>    > <comment-title></comment-title>
>    > The most appropriate **license** are **MIT** or **GNU Affero** (if you want a [copyleft](https://en.wikipedia.org/wiki/Copyleft)).
>    {: .comment}
>
> 4. Add creators in **Creator** with their name and a unique identifier (typically an orcid.org ID)
>
>    > <comment-title></comment-title>
>    > Do not forget to click on **Save**
>    {: .comment}
>
> 5. Add tags in **Tags**
> 6. Save the workflow
{: .hands_on}

Tip 2 is now fulfilled. We can now move toward the other recommendations:

- Tip 1 (**Findability**): *"Register the workflow."*
- Tip 3 (**Accessibility**): *"Make source code available in a public code repository."*
- Tip 4 (**Accessibility**): *"Provide example input data and results along with the workflow."*
- Tip 10 (**Reusability**): *"Provide clear and concise workflow documentation."*

   The workflow annotation provides a short and concise description of the workflow. So Tip 10 is partially fullfilled. In addition to the workflow annotation, each step could annotated and [commented](https://galaxyproject.org/news/2024-04-26-workflows-workflows-workflows/#-workflow-comments) so users could have an idea about the purpose of the steps

   > <details-title>Add step descriptions to fulfill Tip 10</details-title>
   > 
   > Let's add label and annotation to each step and comments to the workflow.
   >
   > > <hands-on-title>Add tool label and annotation and comments</hands-on-title>
   > >
   > > 1. Add tool label and step annotation to **Falco**
   > >    1. Click on **Falco**
   > >    2. Fill in **Label** with `Quality control` in the right panel
   > >    3. Fill in **Step Annotation** with `This step uses Falco (FastQC alternative) to generate statistics of raw reads quality including basic statistics, per base sequence quality, per sequence quality scores, adapter content, etc.`
   > >    4. Save the workflow
   > > 2. Add tool label and step annotation to **MultiQC**
   > > 3. Add comments and annotations to workflow to create an high resolution image of the Galaxy Workflow by following the [dedicated tutorial]({% link topics/galaxy-interface/tutorials/workflow-posters/tutorial.md %})
   > {: .hands_on}
   >
   {: .details}

Before we care about Tip 3 and 4, lets use the build-in Galaxy Wizard to check if our workflow adheres to best practices of Galaxy workflows.

# Adhere to best practices for Galaxy workflows

Galaxy has a build-in Wizard to check against [community developed best practice workflows](https://planemo.readthedocs.io/en/latest/best_practices_workflows.html). It will recommend to use a License, add Authors etc but also ensures that workflows are easy to test, that they are usable as subworkflows and invocation reports, and can be consumed easily via an API. The reusablility will be greatly enhanced if you stick to those recommendations.

> <hands-on-title>Check alignment with best practices</hands-on-title>
>
> 1. Click on {% icon galaxy-wf-best-practices %} **Best Practices** on the right
>
>    {% snippet faqs/galaxy/workflows_best_practices.md %}
>
{: .hands_on}

A side panel will open showing a review of the best practices with
- Green checks indicating that certain elements are already correctly defined (e.g., creator information, license).
- Orange warnings highlighting potential issues such as missing annotations, unconnected inputs, or other workflow problems.

By following the steps with orange warnings and checking the side panel, we will ensure our **workflow adheres to best practices** and is **ready for public use**. 

> <question-title></question-title>
>
> How many orange warnings do we have?
>
> > <solution-title></solution-title>
> >
> > 2
> >
> {: .solution}
{: .question}

The first warning is `Some workflow inputs are missing labels and/or annotations`. To follow best practices, all inputs should be explicit (with labelled input nodes) and tool steps should not have disconnected data inputs (even though the GUI can handle this) or consume workflow parameters. Older style "runtime parameters" should only be used for post job actions and newer type workflow parameter inputs should be used to manipulate tool logic.

> <hands-on-title>Add label to input</hands-on-title>
>
> 1. Click on **Input dataset: Missing a label and annotation**
> 2. Fill in **Label** in the `1: Input dataset` side panel that opened on the right
> 3. Save the workflow
> 4. Check the best practices
{: .hands_on}

> <question-title></question-title>
>
> Is it okay for the input?
>
> > <solution-title></solution-title>
> >
> > The input is missing an annotation 
> >
> {: .solution}
{: .question}

Let's add annotation to input.

> <hands-on-title>Add annotation to input</hands-on-title>
>
> 1. Click on **Input dataset: Missing a label and annotation**
> 2. Fill in **Annotation** in the `1: Input dataset` side panel that opened on the right
> 3. Save the workflow
> 4. Check the best practices
{: .hands_on}

We now have a green check for the input. The second orange warning mentions `This workflow has no labeled outputs, please select and label at least one output.` 

As for inputs, workflows should define explicit, labeled outputs. While Galaxy does not require this, declaring and labeling outputs offers significant advantages. A workflow with clearly defined outputs provides an explicit interface, making it easier to generate reports, test the workflow, and document them.

> <hands-on-title>Add labels to outputs</hands-on-title>
>
> 1. Rename labels for all outputs to add tool name (e.g. Falco `text_file` becomes `falco_text_file`)
> 2. Check boxes on the left of HTML outputs of both Falco and MultiQC in the middle pannel
> 3. Save the workflow
> 4. Check the best practices
{: .hands_on}

> <question-title></question-title>
>
> Is the workflow following all best practices?
>
> > <solution-title></solution-title>
> >
> > Yes!
> >
> {: .solution}
{: .question}

After confirming that all the best practices are applied, **make the workflow public**. This allows other users to access, reuse, and share the workflow in the Galaxy community.

> <hands-on-title> Make the workflow public</hands-on-title>
> 1. Make the workflow public
>
>    {% snippet faqs/galaxy/workflows_publish.md %}
{: .hands_on}

# Prepare example input data and results

To fulfill Tip 4 (*"Provide example input data and results along with the workflow."*), we will **create test data** that can be used **to test the workflow**. The data should be small and simple, so it can interact with the workflow and generate results without excessive processing time. 

> <hands-on-title>Create an example Galaxy History</hands-on-title>
>
> 1. **Create a new Galaxy history** for this analysis.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. **Rename the history**.
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
> 3. {% icon galaxy-upload %} **Upload** test data to the newly created history from the following link:
>
>    ```
>    https://zenodo.org/record/3977236/files/female_oral2.fastq-4143.gz
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 4. {% icon galaxy-eye %} **Preview the datasets** to verify that the uploaded data looks correct and is ready for analysis.
>
> 5. Run the workflow using the uploaded test data.
>
>    {% snippet faqs/galaxy/workflows_run.md %}
>
> 4. Make the **history public**
>
>    {% snippet faqs/galaxy/histories_sharing.md %}
>
{: .hands_on}

With this example history we fulfilled Tip 4. We will now fulfill Tip 1 (*"Register the workflow."*) and Tip 3 (*"Make source code available in a public code repository."*) at the same time

# Register the workflow by making the workflow available in the public IWC GitHub repository

To make your workflow FAIR, it's essential to start by making it findable. Let's now register the workflow in a public registry that enables systematic scientific annotations and supports multiple workflow languages. We recommend using specialized workflow registries such as [WorkflowHub](https://workflowhub.eu/) ({% cite gustafsson2024workflowhub %}) or [Dockstore](https://dockstore.org/) ({% cite yuen2021dockstore %}), which cater to workflows written in different languages and provide unique features like digital object identifiers (DOIs) for easy citation and version tracking.

A Galaxy workflow can be registered by anyone on the WorkflowHub and Dockstore by following their documentation. In this section, we will see how to register a workflow on both WorkflowHub and Dockstore (Tip 1 - *"Register the workflow."*) by making it available in a public GitHub repository (Tip 3 - *"Make source code available in a public code repository."*) maintained by the Intergalactic Workflow Commission (IWC), a Galaxy community effort.

The IWC maintains high-quality Galaxy Workflows via a [Galaxy Workflows Library](https://iwc.galaxyproject.org/). All workflows are reviewed and tested before publication and with every new Galaxy release. Deposited workflows follow best practices and are versioned using GitHub releases. Workflows also contain important metadata. Additionally the IWC collects further best practices, tips and tricks, FAQs and assists the community in designing high-quality Galaxy workflows.

IWC offer [guidelines for adding workflow](https://github.com/galaxyproject/iwc/blob/main/workflows/README.md#adding-workflows) that will go through:

1. Check for workflow eligibility
2. Ensure workflows follow best-practices (already done)
3. Generate tests

## Check for workflow eligibility

IWC collects production workflows targeted at users that want to analyze their own data. As such, the workflow should be sufficiently generic that users can provide their own data.

IWC encourage, but do not require, links to related Galaxy Training Network Tutorials. Importantly, each workflow should be described in a way that a user can run the workflow on their own data without modifying the workflow. If we want to deposit a workflow that accompanies a tutorial we have to make sure that the workflow does not refer to datasets that only make sense in the context of the tutorial.

By fulfilling the first Tips as we did before, we guarantee for workflow eligibility.

## Generate tests

This is usually the most difficult part and we encourage all new contributors to IWC to propose their workflows even if they did not managed to generate tests. However, the publication of these workflow will be speed up if tests are already present. 

### Find input datasets

To test the workflow, we need input datasets. By fulfilling Tip 4 (*"Provide example input data and results along with the workflow."*), we generated a toy dataset. We now need to publish it to [**Zenodo**](https://zenodo.org/) to have a permanent URL, also allowing others to **easily retrieve and reuse the data** when running or validating the workflow.

> <hands-on-title>Upload test data to Zenodo</hands-on-title>  
>
> 1. **Create an account** on [Zenodo](https://zenodo.org/) or **log in** using your **GitHub** account.  
> 2. Click on the **"+"** button at the top of the page and select **"New upload"**.  
> 3. **Drag and drop** the test dataset files into the upload area.  
> 4. **Fill in the metadata**:  
>    - **Resource type**: `Dataset`.  
>    - **Title**: `Dataset for "..." workflow ` 
>    - **Creator**: Add the relevant author(s).  
>    - **Description**: `This dataset is associated with the Galaxy workflow "..."`  
>    - **License**: Choose `GNU General Public License v3.0 or later`
> 5. Click **"Publish"** to make the dataset publicly available.  
>
{: .hands_on}

Uploading test data to Zenodo ensures that it has a **permanent DOI**, making it easy to reference in workflow documentation, publications, and testing pipelines.  

### Generate test from a workflow invocation

To generate tests, we can either write test cases by hand, or use a workflow invocation to generate a test case. We will do the second case. Let's first prepare the folder to store the tests in the IWC git repository.

> <hands-on-title>Prepare the folder for the tests</hands-on-title>  
>
> 1. **Fork** the [IWC GitHub repository](https://github.com/galaxyproject/iwc)to your GitHub account.  
>
> 2. **Clone your fork locally**:  
>
>    ```bash
>    git clone https://github.com/yourusername/iwc.git
>    cd iwc
>    ```
>
> 3. **Create a new branch** for adding the workflow:  
>
>    ```bash
>    git checkout -b add-new-workflow
>    ```
>
> 4. **Create a new directory** under one of the directories that represent categories
>
>    > <comment-title></comment-title>
>    > If no category is suitable, we can create a new category directory. We should name the directory that contains our workflow(s) appropriately, as it will become the name of the repository deployed to [iwc-workflows GitHub organization](https://github.com/iwc-workflows) and only use lower-case and `-` in names of categories and repositories.
>    {: .comment}
>
> 5. Move the newly created directory
>
{: .hands-on}

We have earlier executed the workflow on our example dataset. We can use this **workflow invocation** to generate the tests using Planemo, a software development kit for tools and workflows ({% cite bray2023planemo %})

> <hands-on-title>Extract the tests from the example history workflow invocation</hands-on-title>  
>
> 1. **Install Planemo** (if not already installed) by following the [Planemo installation guide](https://planemo.readthedocs.io/en/stable/installation.html).  
>
> 3. Get your **Galaxy API key** from the Galaxy server.
>
>    {% snippet faqs/galaxy/preferences_admin_api_key.md %} 
>
> 4. Get the workflow invocation.
>
>    {% snippet faqs/galaxy/workflows_get_invocation_id.md %} 
>
> 5. Extract the workflow with test data using:  
>
>    ```bash
>    planemo workflow_test_init \
>       --from_invocation <workflow_invocation_id> \
>       --galaxy_url <galaxy_server_url> \
>       --galaxy_user_key <your_api_key>
>    ```
>
{: .hands-on}

This will place in the current working directory:
- a `<workflow_name>.ga` with the workflow
- a `<workflow_name-test>.yml` to describe the tests
- a `test-data` folder with input file and selected outputs

> <question-title></question-title>
>
> In the `<workflow_name-test>.yml`,
> 1. How is the input file provided?
> 2. How are the workflow outputs checked for validity?
> 3. How big are the files in the `test-data` file?
>
> > <solution-title></solution-title>
> >
> > 1. The input file is given by:
> >
> >    ```
> >    job:
> >      FastQ:
> >        class: File
> >        path: test-data/FastQ.fastqsanger.gz
> >        filetype: fastqsanger.gz
> >    ```
> >
> > 2. The workflow outputs are checked for validity by comparing the generated outputs to the files stored in the `test-data` folder.
> > 
> >    ```
> >    outputs:
> >      multiqc_html:
> >        path: test-data/multiqc_html.html
> >      falco_html:
> >        path: test-data/falco_html.html
> >    ```
> >
> > 3. The files are quite large:
> >    - `FastQ.fastqsanger.gz`: 132 KB
> >    - `falco_html.html`: 169 KB
> >    - `multiqc_html.html`: 5 MB
> {: .solution}
{: .question}

The HTML files in the `test-data` folder are quite large, above the 1 MB limit of IWC. To limit the size in the IWC GitHub repository, we will remove test files in `test-data` folder. To do that, we will:
1. use the file stored on Zenodo for the input.
2. edit the test comparisons to use assertions testing the output content, rather than comparing the entire output file with test data.

   > <details-title>Assertion documentation</details-title>
   > The description of assertion can be find in the [Galaxy XML documentation](https://docs.galaxyproject.org/en/latest/dev/schema.html#tool-tests-test-output-assert-contents).
   >
   > Do not hesitate to look at different test files in the [IWC GitHub](https://github.com/galaxyproject/iwc/tree/main/workflows) for examples.
   {: .details}
   
> <hands-on-title>Edit the test description</hands-on-title>
>
> 1. Remove the `test-data` folder
> 2. Open the `<workflow_name-test>.yml` file
> 
> 3. Replace the `path` of input by the **Zenodo URL**:
>
>      ```yaml
>      https://zenodo.org/record/3977236/files/female_oral2.fastq-4143.gz
>      ```
>
> 4. Ensure all test **output names** are included.
>
> 5. For Falco step 
>    
>    1. Open the Falco HTML output on Galaxy
>    2. Search for a text specific to this data (e.g. `female_oral2_fastq-4143_gz.gz` here)
>    3. Replace in the `<workflow_name-test>.yml` file the `path` line of `falco_html` output by:
> 
>       ```
>       asserts:
>          has_text:
>            text: "female_oral2_fastq-4143_gz.gz"
>       ```
>
> 6. Do the same for MultiQC
> 7. Save the file
>
{: .hands-on}

> <comment-title></comment-title>
> If the workflow is using build-in indexes, we should use the available indexes on [CernVM-FS (CVMFS)](https://datacache.galaxyproject.org/), a distributed filesystem perfectly designed for sharing readonly data across the globe. IWC continuous integration uses it for the tests.
{: .comment}

### Lint the workflow

To make sure the workflow and its test are syntactically correct, we can now run `planemo workflow_lint` :

> <hands-on-title>Lint the workflow</hands-on-title>  
> 1. **Run the tests** on the extracted workflow using
>
>     ```bash
>     planemo workflow_lint <workflow_file.ga>
>     ```
{: .hands-on}

### Test the workflow against an instance which have all tools installed.

Before submitting the workflow to IWC, it might be interesting to run the tests against a Galaxy instance using Planemo so we can easily see what is failing and what are the differences between our expectations and the output we get:

> <hands-on-title>Run the tests</hands-on-title>  
> 1. **Run the tests** on the extracted workflow using
>
>     ```bash
>     planemo test \
>        --galaxy_url <galaxy_server_url> \
>        --galaxy_user_key <your_api_key> \
>        <workflow_file.ga>
>     ```
{: .hands-on}

> <details-title>Using Planemo to check that the test is valid against the same workflow invocation</details-title>
>
> If the tests are not passing because of an error into the test file, we can modify the test file and use Planemo to check that the test is valid against the same invocation.
>
> ```bash
> planemo workflow_test_on_invocation \
>    --galaxy_url <galaxy_server_url> \
>    --galaxy_user_key <your_api_key> \
>    <workflow-tests.yml> \
>    <workflow_invocation_id>
> ```
{: .details}


### Add required metadata

Once the workflow tests has been created and tested, we are almost ready to submit it to IWC. Few last steps are missing to add required metadata


The first step is to generate a `.dockstore.yml` file that contains metadata needed for Dockstore.

> <hands-on-title>Generate a `.dockstore.yml` file</hands-on-title>  
>
> 1. Run:  
>    
>    ```bash
>    planemo dockstore_init
>    ```
> 
> 2. Open the generated `docker.yml` and verify that author names and details are correct.  
> 
{: .hands-on}

We now need 
1. a `README.md` file that briefly describes the workflow
2. a `CHANGELOG.md` file that lists changes, additions and fixed. 

   For that, we need to follow the formatting and principles proposed on [keepachangelog.com](https://keepachangelog.com/en/1.0.0/).

> <hands-on-title>Create README and CHANGELOG files</hands-on-title>  
>
> 1. Create a `README.md` file inside the workflow folder ([example README](https://github.com/galaxyproject/iwc/blob/main/workflows/microbiome/pathogen-identification/nanopore-pre-processing/README.md)).  
>
> 2. Create a `CHANGELOG.md` file inside the workflow folder given the following template:
>
>    ```md
>    # Changelog
>    
>    ## [0.1] yyyy-mm-dd
>    
>    First release.
>    ```
> 
{: .hands-on}

Finally, there is currently no user interface within Galaxy to define release versions, so we have to manually set a `release: "0.1"` key value pair in the `.ga` file. 

> <hands-on-title>Edit the release version in the workflow</hands-on-title>  
>
> 1. Edit the `workflow.ga` file to include the **release number**
>
>    ```bash
>    ],
>    "format-version": "0.1",
>    "license": "GPL-3.0-or-later",
>    "release": "0.1",
>    "name": "NameOfTheWorkflow",
>    "steps": {
>    "0": {
>    ```
>
{: .hands_on} 

## Submit the workflow to IWC

> <hands-on-title>Submit the workflow to IWC</hands-on-title>  
>
> 1. Commit the changes to the createw branch
> 2. Push the branch to your fork
> 3. **Open a Pull Request (PR)** on the IWC GitHub repository to submit your workflow.  
>
{: .hands_on}

Once the PR is reviewed and merged, the workflow will be available in **Galaxy's IWC GitHub repository** (Tip 3 - *"Make source code available in a public code repository."*) and indexed in both [**WorkflowHub**](https://workflowhub.eu/workflows/) and [**Dockstore**](https://dockstore.org/workflows/github.com/iwc-workflows/) registries (Tip 1 - *"Register the workflow."*), making it **discoverable and reusable** by the community.

> <comment-title></comment-title>
> The [**Best practices for workflows in GitHub repositories**]({% link topics/fair/tutorials/ro-crate-galaxy-best-practices/tutorial.md %}) GTN provides additional recommendations for structuring and maintaining workflows in Github.
{: .comment}

# Conclusion

By following this tutorial, you have ensured that your Galaxy workflow is following the **"Ten quick tips for building FAIR workflows"** ({% cite de2023ten %})  and adheres to **best practices**. From creating and refining the workflow to creating test data and submitting it to the **IWC**, each step contributes to making workflows more **reliable, shareable, and reusable**.  

For further workflow optimization and maintenance, the **Galaxy community** provides a [comprehensive guide on best practices](https://planemo.readthedocs.io/en/latest/best_practices_workflows.html). This guide extends beyond the **best practices panel in Galaxy** and includes additional recommendations, such as:  

- **Adding automated tests** to validate workflow functionality.  
- **Publishing workflows** on platforms like GitHub, GitLab, or other public version-controlled repositories.  
- **Registering workflows** in well-known registries like [WorkflowHub](https://workflowhub.eu/workflows/) and [Dockstore](https://dockstore.org/workflows/github.com/iwc-workflows/), ensuring wider accessibility and discoverability.  

By continuously following these practices, you contribute to a **stronger, more open, and collaborative** bioinformatics community, where workflows are easily shared, improved, and adapted to new challenges.