---
layout: tutorial_hands_on

title: "Creating a new tutorial - Defining the technical infrastructure"
questions:
  - "How can we define the technical infrastructure for a tutorial?"
  - "How to define the tools needed for a tutorial?"
  - "How to add the needed data directly in an instance?"
  - "How to add the workflows related to a tutorial?"
  - "How can we check the technical infrastructure is working?"
  - "How can we make an existing Galaxy instance able to run a tutorial?"
objectives:
  - "Extracting the technical description for a tutorial"
  - "Populating an existing instance with the needed tools, data and workflows for a tutorial"
  - "Creating a Galaxy Docker flavor with the needed tools, data and workflows for a tutorial"
  - "Testing the Galaxy Docker flavor of a tutorial"
time_estimation: "30m"
key_points:
  - "Tools, data and workflows can be easily integrated in a Docker flavor to have a useful technical support for a tutorial"
  - "A Galaxy Docker flavor is a great support for training"
  - "A Galaxy Docker flavor can be deployed 'anywhere' and is scalable"
contributors:
  - bebatut
  - bgruening
  - shiltemann
  - erasche
---

# Building a Galaxy instance specifically for your training
{:.no_toc}

To be able to run the tutorial, we need a Galaxy instance where all of the needed tools and data are available. Thus we need to describe the needed technical infrastructure.

This files we define in this tutorial will be used to automatically build a Docker Galaxy flavour, and also to test if a public Galaxy instance is able to run the tool.

In this tutorial, you will learn how to create a virtualised Galaxy instance, based on Docker, to run your training - either on normal computers or cloud environments.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Extracting workflows

Once the tutorial is ready, we need to develop a workflow that represents the steps taken in the tutorial, and then extract these workflow(s) and add them to the `workflows` directory in the tutorial. Additionally we will need to add some explanation about the workflow(s) in a `README.md` file

> ### {% icon hands_on %} Hands-on: Extract the workflow
>
> 1. Download the workflow for the tutorial
> 2. Save it in the `workflow` directory of the tutorial
{: .hands_on}

## Creating the `tools.yaml` (recommended)

The first file to fill out is the `tools.yaml` file which contains the list of the required tools that could be installed from the ToolShed.

This file looks like:

```yaml
---
api_key: admin
galaxy_instance: http://localhost:8080
tools:
- name: tool1
  owner: owner
  tool_panel_section_label: "Section1"
- name: tool2
  owner: owner
  tool_panel_section_label: "Section2"
```

with:

- `name`: the name of the wrapper of the tool in the ToolShed
- `owner`: the owner of the wrapper of the tool in the ToolShed
- `tool_panel_section_label`: section where to put the tool (in the left panel in the Galaxy instance)

This list of tools can be automatically extracted from the workflow using [Ephemeris](https://ephemeris.readthedocs.io/en/latest/index.html) (which should be in the conda environment):

```console
$ workflow-to-tools -w path/to/workflow -o path/to/tools.yaml
```

After the extraction, some formatting is needed:

1. Add at the beginning:

    ```yaml
    ---
    api_key: admin
    galaxy_instance: http://localhost:8080
    ```

2. Change the `tool_panel_section_label` to something more informative

> ### {% icon hands_on %} Hands-on: Creating the `tools.yaml` from your workflow
>
> 1. Create the `tools.yaml` file using your workflow and Ephemeris
> 2. Correct the formatting of the `tools.yaml` file
{: .hands_on}

## Testing the workflow (recommended)

Workflow testing is a great way to get feedback that your tutorial can be run successfully on a given server. When you're giving a training this can provide peace of mind, not only are the tools installed (as is indicated by the badges we provide) but they also work.

Given the workflow you created above and have included in the tutorial folder, you'll need to create a corresponding `-test.yml` file.

> ### {% icon hands_on %} Hands-on: Creating the `-test.yml` file for your workflow
>
> 1. Find the correct name for the file; if your workflow was `unicycler.ga`, then your test file should be `unicycler-test.yml`, they need to share the same prefix.
>
> 2. Create the following structure:
>
>    ```yaml
>    ---
>    - doc: Test sample data for the workflow
>      job:
>        an_input_file:
>          class: File
>          location: https://....
>          filetype: fasta
>      outputs:
>        ffn:
>          asserts:
>            has_text:
>              text: ">A"
>            has_text:
>              text: ">B"
>    ```
>
>
{: .hands_on}

You'll need to edit the `job` and `outputs` sections according to your workflow's inputs and outputs. Additionally you will need to edit the steps of your workflow `.ga` file appropriately.

### Inputs

Your workflow **must** use "Data Inputs" for each input dataset. For each of these input step in the `.ga` file, you'll need to do the following:

1. Edit the `label`
2. Edit the `name`
3. Edit the `inputs[0].name`
4. Edit the `tool_state`

In a normal workflow you have exported from Galaxy, you'll see something like

```json
{
    "id": 0,
    "input_connections": {},
    "inputs": [
        {
            "description": "",
            "name": "patient1_ChIP_ER_good_outcome.bam"
        }
    ],
    "label": null,
    "name": "Input dataset",
    "outputs": [],
    "position": {
        "left": 10,
        "top": 10
    },
    "tool_id": null,
    "tool_state": "{\"name\": \"patient1_ChIP_ER_good_outcome.bam\"}",
    "tool_version": null
}
```

You should synchronize the aforementioned fields so it looks like this:

```json
{
    "id": 0,
    "input_connections": {},
    "inputs": [
        {
            "description": "",
            "name": "good_outcome"
        }
    ],
    "label": "good_outcome",
    "name": "good_outcome",
    "outputs": [],
    "position": {
        "left": 10,
        "top": 10
    },
    "tool_id": null,
    "tool_state": "{\"name\": \"good_outcome\"}",
    "tool_version": null
}
```

This will allow you to specify `good_outcome` in your job to load a file:

```
- doc: ...
  job:
    good_outcome:
      class: File
      location: ...
      filetype: ...
```

The filetype should be the Galaxy datatype of your file, for example `fastqsanger`, `tabular`, `bam`.

### Outputs

For the outputs the process is somewhat simpler:

1. Identify a step, the outputs of which you would like to test
2. Convert the relevant `outputs` to `workflow_outputs`

   In a normal workflow you see

   ```json
   {
       "outputs": [
           {
               "type": "txt",
               "name": "ofile"
           },
           {
               "type": "txt",
               "name": "ofile2"
           }
       ],
       "workflow_outputs": []
   }
   ```

   If you want to test the contents of `ofile`, you should change it to

   ```json
   {
       "outputs": [
           {
               "type": "txt",
               "name": "ofile"
           },
           {
               "type": "txt",
               "name": "ofile2"
           }
       ],
       "workflow_outputs": [
           {"output_name": "ofile", "label": "my_output"}
       ]
   }
   ```

3. You can now use the label you chose (here `my_output`) in your test case:

   ```yaml
   - doc:
     job: ...
     outputs:
       my_output:
         asserts:
           has_text:
             text: 'some-string'
   ```

### Running the Tests

You can test the file you've written with the following command and a recent version (>=0.56.0) of planemo:

```console
planemo test \
	--galaxy_url "$GALAXY_URL" \
	--galaxy_user_key "$GALAXY_USER_KEY" \
	--no_shed_install \
	--engine external_galaxy \
	workflow.ga
```

Planemo will autodetect that the `workflow-test.yml` file and load that for the testing.

# Creating the `data-library.yaml` (recommended)

The datasets needed for a tutorial can also be integrated in the Galaxy instance inside of data libraries. These allow the datasets to be easily shared with all users of a Galaxy instance. Additionally it lets trainees avoid each re-downloading the input data.

These datasets are described in the `data-library.yaml` files:

```yaml
---
destination:
  type: library
  name: GTN - Material
  description: Galaxy Training Network Material
  synopsis: Galaxy Training Network Material. See https://training.galaxyproject.org
items:
- name: Title of the topic
  description: Summary of the topic
  items:
  - name: Title of the tutorial
    items:
    - name: 'DOI: 10.5281/zenodo....'
      description: latest
      items:
      - info: https://doi.org/10.5281/zenodo....
        url: https://zenodo.org/api/files/URL/to/the/input/file
        ext: galaxy-datatype
        src: url
```

> ### {% icon hands_on %} Hands-on: Creating the `data-library.yaml`
>
> 1. Copy the Zenodo link
> 2. Generate the `data-library.yaml` file and update the tutorial metadata with the link:
>
>    ```
>    $ planemo training_fill_data_library \
>             --topic_name "my-topic" \
>             --tutorial_name "my-new-tutorial" \
>             --zenodo_link "URL to the Zenodo record"
>    ```
>
> 3. Check that the `data-library.yaml` has been generated (or updated)
> 4. Check tha the Zenodo link is in the metadata at the top of the `tutorial.md`
{: .hands_on}

# Creating the `data-manager.yaml` (optional)

Some of the tools may require specific databases, specifically prepared for the tool. In this case, some Galaxy tools come with "data managers" to simplify this process.

If you need such data managers for your training, then you should describe how to run them in the `data-manager.yaml` file:

```yaml
data_managers:
    - id: url to data manager on ToolShed
      params:
        - 'param1': '{{ item }}'
        - 'param2': 'value'
      # Items refer to a list of variables you want to run this data manager. You can use them inside the param field with {{ item }}
      # In case of genome for example you can run this DM with multiple genomes, or you could give multiple URLs.
      items:
        - item1
        - item2
      # Name of the data-tables you want to reload after your DM are finished. This can be important for subsequent data managers
      data_table_reload:
        - all_fasta
        - __dbkeys__
```

# Creating the Galaxy Interactive Tour (optional)

A Galaxy Interactive Tour is a way to go through an entire analysis, step by step inside Galaxy in an interactive and explorative way.
It is a great way to help users run the tutorial directly inside Galaxy. To learn more about creating a Galaxy tour please have a look at our [dedicated tour training]({{site.baseurl}}/topics/contributing/tutorials/create-new-tutorial-tours/tutorial.html).

# Testing the technical infrastructure

Once we have defined all the requirements for running the tutorial, we can test these requirements, either in a locally running Galaxy or in a Docker container. Please see our tutorial about [Setting up Galaxy for Training]({{site.baseurl}}/topics/instructors/tutorials/setup-galaxy-for-training/tutorial.html) about how to test your tutorial requirements.


# Conclusion
{:.no_toc}
