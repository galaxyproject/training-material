---
layout: tutorial_hands_on
title: "RO-Crate in Python"
time_estimation: "30M"
questions:
- What data is contained within an RO-Crate
- How can I create an RO-Crate myself?
objectives:
- Create a custom, annotated RO-Crate
- Use ORCIDs and other linked data to annotate datasets contained within the crate
key_points:
- RO-Crates can be created by hand with essentially arbitrary data, using the rocrate python module
- However the rocrate command line tool adds several commands to make it easier to automatically generate crates based on existing folder structures.
tags:
- ro-crate
priority: 2
contributions:
  authorship:
    - simleo
    - kinow
  editing:
    - hexylena
    - elichad
  testing:  
    - elichad
  funding:
    - by-covid
license: Apache-2.0
subtopic: ro-crate
notebook:
  language: python
---

This tutorial will show you how to manipulate [RO-Crates](https://w3id.org/ro/crate/) in Python using the [ro-crate-py](https://github.com/ResearchObject/ro-crate-py) package. It is based on the [ro-crate-py documentation](https://github.com/ResearchObject/ro-crate-py/blob/e1218fbca595f4c33059cfe15849ee2ae9e6896b/README.md).

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

Let's start by installing the library via [pip](https://docs.python.org/3/installing/). Note that the name of the package is `rocrate`.

```bash
pip install rocrate
```

# Creating an RO-Crate

In its simplest form, an RO-Crate is a directory tree with an `ro-crate-metadata.json` file at the top level. This file contains metadata about the other files and directories, represented by [data entities](https://www.researchobject.org/ro-crate/1.1/data-entities.html). These metadata consist both of properties of the data entities themselves and of other, non-digital entities called [contextual entities](https://www.researchobject.org/ro-crate/1.1/contextual-entities.html). A contextual entity can represent, for instance, a person, an organization or an event.

Suppose Alice and Bob worked on a research project together, and then wrote a paper about it; additionally, Alice prepared a spreadsheet containing experimental data, which Bob then used to generate a diagram. For the purpose of this tutorial, you can just create placeholder files for the documents:

```python
import os

data_dir = "exp"
os.mkdir(data_dir)

for filename in ["paper.pdf", "results.csv", "diagram.svg"]:
    with open(os.path.join(data_dir, filename), "w") as file:
        pass
```

Let's make an RO-Crate to represent this information:

```python
from rocrate.rocrate import ROCrate

crate = ROCrate()
paper = crate.add_file("exp/paper.pdf", properties={
    "name": "manuscript",
    "encodingFormat": "application/pdf"
})
table = crate.add_file("exp/results.csv", properties={
    "name": "experimental data",
    "encodingFormat": "text/csv"
})
diagram = crate.add_file("exp/diagram.svg", dest_path="images/figure.svg", properties={
    "name": "bar chart",
    "encodingFormat": "image/svg+xml"
})
```

We've started by adding the data entities. Now we add contextual entities representing Alice and Bob:

```python
from rocrate.model.person import Person

alice_id = "https://orcid.org/0000-0000-0000-0000"
bob_id = "https://orcid.org/0000-0000-0000-0001"
alice = crate.add(Person(crate, alice_id, properties={
    "name": "Alice Doe",
    "affiliation": "University of Flatland"
}))
bob = crate.add(Person(crate, bob_id, properties={
    "name": "Bob Doe",
    "affiliation": "University of Flatland"
}))
```

At this point, we have a representation of the various entities. Now we need to express the relationships between them. This is done by adding properties that reference other entities:

```python
paper["author"] = [alice, bob]
table["author"] = alice
diagram["author"] = bob
```

You can also add whole directories together with their contents. In an RO-Crate, a directory is represented by the `Dataset` entity:

```python
logs_dir = os.path.join(data_dir, "logs")
os.mkdir(logs_dir)

for filename in ["log1.txt", "log2.txt"]:
    with open(os.path.join(logs_dir, filename), "w") as file:
        pass

logs = crate.add_dataset("exp/logs")
```

Finally, we serialize the crate to disk:

```python
crate.write("exp_crate")
```

This should generate an `exp_crate` directory containing copies of all the files we added and an `ro-crate-metadata.json` file containing a [JSON-LD](https://json-ld.org) representation of the metadata. Note that we have chosen a different destination path for the diagram, while the paper and the spreadsheet have been placed at the top level with their names unchanged (the default).

Some applications and services support RO-Crates stored as archives. To save the crate in zip format, you can use `write_zip`:

```python
crate.write_zip("exp_crate.zip")
```

> <comment-title>How `rocrate` handles the contents of `exp/logs`</comment-title>
>
> Exploring the `exp_crate` directory, we see that all files and directories contained in `exp/logs` have been added recursively to the crate. However, in the `ro-crate-metadata.json` file, only the top level Dataset with `@id` `"exp/logs"` is listed. This is because we used `crate.add_dataset("exp/logs")` rather than adding every file individually. There is no requirement to represent every file and folder within the crate in the `ro-crate-metadata.json` file - in fact, if there were many files in the crate it would be impractical to do so.
>
> If you do want to add files and directories recursively to the metadata, use `crate.add_tree` instead of `crate.add_dataset` (but note that it only works on local directory trees).
{: .comment}

## Appending elements to property values

What ro-crate-py entities actually store is their JSON representation:

```python
paper.properties()
```

```json
{
  "@id": "paper.pdf",
  "@type": "File",
  "name": "manuscript",
  "encodingFormat": "application/pdf",
  "author": [
    {"@id": "https://orcid.org/0000-0000-0000-0000"},
    {"@id": "https://orcid.org/0000-0000-0000-0001"},
  ]
}
```

When `paper["author"]` is accessed, a new list containing the `alice` and `bob` entities is generated on the fly. For this reason, calling `append` on `paper["author"]` won't actually modify the `paper` entity in any way. To add an author, use the `append_to` method instead:

```python
donald = crate.add(Person(crate, "https://en.wikipedia.org/wiki/Donald_Duck", properties={
  "name": "Donald Duck"
}))
paper.append_to("author", donald)
```

Note that `append_to` also works if the property to be updated is missing or has only one value:

```python
for n in "Mickey_Mouse", "Scrooge_McDuck":
    p = crate.add(Person(crate, f"https://en.wikipedia.org/wiki/{n}"))
    donald.append_to("follows", p)
```

## Adding remote entities

Data entities can also be remote:

```python
input_data = crate.add_file("http://example.org/exp_data.zip")
```

By default the file won't be downloaded, and will be referenced by its URI in `ro-crate-metadata.json`:

```json
{
  "@id": "http://example.org/exp_data.zip",
  "@type": "File"
}
```

If you add `fetch_remote=True` to the `add_file` call, however, the library (when `crate.write` is called) will try to download the file and include it in the output crate.

Another option that influences the behavior when dealing with remote entities is `validate_url`, also `False` by default: if it's set to `True`, when the crate is serialized, the library will try to open the URL to add / update metadata such as the content's length and format.

## Adding entities with an arbitrary type

An entity can be of any type listed in the [RO-Crate context](https://www.researchobject.org/ro-crate/1.1/context.jsonld). However, only a few of them have a counterpart (e.g., `File`) in the library's class hierarchy, either because they are very common or because they are associated with specific functionality that can be conveniently embedded in the class implementation. In other cases, you can explicitly pass the type via the `properties` argument:

```python
from rocrate.model.contextentity import ContextEntity

hackathon = crate.add(ContextEntity(crate, "#bh2021", properties={
    "@type": "Hackathon",
    "name": "Biohackathon 2021",
    "location": "Barcelona, Spain",
    "startDate": "2021-11-08",
    "endDate": "2021-11-12"
}))
```

Note that entities can have multiple types, e.g.:

```
    "@type" = ["File", "SoftwareSourceCode"]
```

# Consuming an RO-Crate

An existing RO-Crate package can be loaded from a directory or zip file:

```python
crate = ROCrate('exp_crate')  # or ROCrate('exp_crate.zip')
for e in crate.get_entities():
    print(e.id, e.type)
```

```
./ Dataset
ro-crate-metadata.json CreativeWork
paper.pdf File
results.csv File
images/figure.svg File
https://orcid.org/0000-0000-0000-0000 Person
https://orcid.org/0000-0000-0000-0001 Person
...
```

The first two entities shown in the output are the [root data entity](https://www.researchobject.org/ro-crate/1.1/root-data-entity.html) and the [metadata file descriptor](https://www.researchobject.org/ro-crate/1.1/metadata.html), respectively. The former represents the whole crate, while the latter represents the metadata file. These are special entities managed by the `ROCrate` object, and are always present. The other entities are the ones we added in the [section on RO-Crate creation](#creating-an-ro-crate). As shown above, `get_entities` allows to iterate over all entities in the crate. You can also access only data entities with `crate.data_entities` and only contextual entities with `crate.contextual_entities`. For instance:

```python
for e in crate.data_entities:
    author = e.get("author")
    if not author:
        continue
    elif isinstance(author, list):
        print(e.id, [p.get("name") for p in author])
    else:
        print(e.id, repr(author.get("name")))
```

```
paper.pdf ['Alice Doe', 'Bob Doe']
results.csv 'Alice Doe'
images/figure.svg 'Bob Doe'
```

You can fetch an entity by its `@id` as follows:

```python
article = crate.dereference("paper.pdf")  # or crate.get("paper.pdf")
```

# Command Line Interface

> <comment-title>Jupyter Notebook users: switch to a terminal</comment-title>
> The code cells in this section use Unix shell commands, which can't be run within a notebook. Open a Unix/Linux terminal to follow along.
{: .comment}

`ro-crate-py` includes a hierarchical command line interface: the `rocrate` tool. `rocrate` is the top-level command, while specific functionalities are provided via sub-commands. Currently, the tool allows to initialize a directory tree as an RO-Crate (`rocrate init`) and to modify the metadata of an existing RO-Crate (`rocrate add`).

```console
$ rocrate --help
Usage: rocrate [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  add
  init
  write-zip
```

## Crate initialization

The `rocrate init` command explores a directory tree and generates an RO-Crate metadata file (`ro-crate-metadata.json`) listing all files and directories as `File` and `Dataset` entities, respectively.

```console
$ rocrate init --help
Usage: rocrate init [OPTIONS]

Options:
  --gen-preview
  -e, --exclude CSV
  -c, --crate-dir PATH
  --help                Show this message and exit.
```

The command acts on the current directory, unless the `-c` option is specified. The metadata file is added (overwritten if present) to the directory at the top level, turning it into an RO-Crate.

## Adding items to the crate

The `rocrate add` command allows to add files, datasets (directories), workflows, and other entity types (currently [testing-related metadata](https://crs4.github.io/life_monitor/workflow_testing_ro_crate)) to an RO-Crate:

```console
$ rocrate add --help
Usage: rocrate add [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  dataset
  file
  test-definition
  test-instance
  test-suite
  workflow
```

Note that data entities (e.g., workflows) must already be present in the directory tree: the effect of the command is to register them in the metadata file.

## Example

To run the following commands, we need a copy of the ro-crate-py repository:

```bash
git clone https://github.com/ResearchObject/ro-crate-py
```

Navigate to the following directory in the repository we just cloned:

```bash
cd ro-crate-py/test/test-data/ro-crate-galaxy-sortchangecase
```

This directory is already an RO-Crate. Delete the metadata file to get a plain directory tree:

```bash
rm ro-crate-metadata.json
```

Now the directory tree contains several files and directories, including a Galaxy workflow and a Planemo test file, but it's not an RO-Crate anymore, since there is no metadata file. Initialize the crate:

```bash
rocrate init
```

This creates an `ro-crate-metadata.json` file that lists files and directories rooted at the current directory. Note that the Galaxy workflow is listed as a plain `File`:

```json
{
  "@id": "sort-and-change-case.ga",
  "@type": "File"
}
```

To register the workflow as a `ComputationalWorkflow`, run the following:

```bash
rocrate add workflow -l galaxy sort-and-change-case.ga
```

Now the workflow has a type of `["File", "SoftwareSourceCode", "ComputationalWorkflow"]` and points to a `ComputerLanguage` entity that represents the Galaxy workflow language. Also, the workflow is listed as the crate's `mainEntity` (this is required by the [Workflow RO-Crate profile](https://w3id.org/workflowhub/workflow-ro-crate/1.0), a subtype of RO-Crate which provides extra specifications for workflow metadata).

To add files or directories after crate initialization:

```bash
cp ../sample_file.txt .
rocrate add file sample_file.txt -P name=sample -P description="Sample file"
cp -r ../test_add_dir .
rocrate add dataset test_add_dir
```

The above example also shows how to set arbitrary properties for the entity with -P. This is supported by most `rocrate add` subcommands.
