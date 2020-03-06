---
layout: tutorial_hands_on

title: Data Libraries
questions:
- How do data libraries work?
- What are they good for?
- How can I use them?
- How can I setup permissions for them?
objectives:
- Setup a data library
- Manage permissions
- Import data from disk
time_estimation: "30m"
key_points:
- Data libraries are a great way to share data with groups of users
contributors:
  - hexylena
subtopic: features
tags:
  - storage
requirements:
 - type: "internal"
   topic_name: admin
   tutorials:
     - ansible
     - ansible-galaxy
---


> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Adding local data

Before we can import local data, we need to configure Galaxy to permit this:

> ### {% icon hands_on %} Hands-on: Setting up Grafana
>
> 1. We will add a pre-task to clone [a data repository](https://github.com/galaxyproject/galaxy-test-data) into your machine. We will use this as the source for a library dataset.
>
>
>
>        - name: Create the second storage directory
>          file:
>            owner: galaxy
>            group: galaxy
>            path: /libraries/
>            state: directory
>            mode: '0755'
>        - name: Create the second storage directory
>          file:
>            owner: galaxy
>            group: galaxy
>            path: /libraries/user/
>            state: directory
>            mode: '0755'
>        - git:
>            repo: 'https://github.com/galaxyproject/galaxy-test-data'
>            dest: /libraries/admin
>
> 4. Edit the file `group_vars/galaxyservers.yml` and set the following variables:
>
>    ```yaml
>    galaxy_config:
>      galaxy:
>        library_import_dir: /libraries/admin
>        user_library_import_dir: /libraries/user
>        user_library_import_dir_auto_creation: true
>    ```
>
> 5. Run the playbook:
>
>    ```
>    ansible-playbook galaxy.yml
>    ```
>
{: .hands_on}

data is set up

- walk through import interface
- unprivileged user import?


# Installing remote data

if your data is accessible via URL, you can write a yaml file to import and setup the data library automatically:

```
---
destination:
  type: library
  name: Mouse sequencing project
  description: some data
  synopsis: samples collected from somewhere
items:
- url: https://zenodo.org/api/files/287966da-5411-4f79-8cfb-0ffa84d0d6cc/wildtype.fna
  src: url
  ext: fasta
  info: https://doi.org/10.5281/zenodo.582600
- name: A directory
  description: Exome sequencing means that all protein-coding genes in a genome are
  items:
  - url: https://zenodo.org/api/files/287966da-5411-4f79-8cfb-0ffa84d0d6cc/mutant_R1.fastq
    src: url
    ext: fastqsanger
    info: https://doi.org/10.5281/zenodo.582600
  - url: https://zenodo.org/api/files/287966da-5411-4f79-8cfb-0ffa84d0d6cc/mutant_R2.fastq
    src: url
    ext: fastqsanger
    info: https://doi.org/10.5281/zenodo.582600
```

- Write to a file
- run ephemeris:

    setup-data-libraries -g https://gat-0.student.galaxy.training -u admin@example.org -p galaxy --training -i data-library.yaml --legacy


# GTN Data

Join the GTN shared data repository:

https://github.com/usegalaxy-eu/shared-data

We will keep the GTN data on your server, updated.
