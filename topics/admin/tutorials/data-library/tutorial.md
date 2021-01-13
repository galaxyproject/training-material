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

TODO: overview.

- save space for users, files don't count against quota
- great for e.g. sharing sequencing run with a bunch of people.

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Adding local data

Before we can import local data, we need to configure Galaxy to permit this. Additionally we will setup an example data library which we can use for demonstrative purposes.

> ### {% icon hands_on %} Hands-on: Setting up Grafana
>
> 1. We will add a pre-task to clone [a data repository](https://github.com/galaxyproject/galaxy-test-data) into your machine. We will use this as the source for a library dataset.
>
>    ```diff
>    --- a/galaxy.yml
>    +++ b/galaxy.yml
>    @@ -5,6 +5,9 @@
>         - name: Install Dependencies
>           package:
>             name: ['git', 'make', 'python3-psycopg2', 'virtualenv', 'tar', 'bzip2']
>    +    - git:
>    +        repo: 'https://github.com/usegalaxy-eu/libraries-training-repo'
>    +        dest: /libraries/
>       handlers:
>         - name: Restart Galaxy
>           systemd:
>    ```
>
> 4. Edit the file `group_vars/galaxyservers.yml` and set the following variables:
>
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -28,6 +28,9 @@ miniconda_manage_dependencies: false
>
>     galaxy_config:
>       galaxy:
>    +    library_import_dir: /libraries/admin
>    +    user_library_import_dir: /libraries/user
>         job_resource_params_file: "{{ galaxy_config_dir }}/job_resource_params_conf.xml"
>         tool_destinations_config_file: "{{ galaxy_config_dir }}/tool_destinations.yml"
>    ```
>
> 5. Run the playbook:
>
>    > ### {% icon code-in %} Input: Bash
>    > ```
>    > ansible-playbook galaxy.yml
>    > ```
>    {: .code-in}
>
{: .hands_on}

# Importing Data

There are multiple options for importing data from your server, we'll go through all of your choices below. But first, let's take a quick look at the example library structure we've provided.

> > ### {% icon code-in %} Input: Bash
> > ```bash
> > tree /libraries
> > ```
> {: .code-in}
>
> > ### {% icon code-out %} Output: Bash
> > ```
> > /libraries/
> > ├── admin
> > │   └── admin-wildtype.fna
> > ├── example-library.yaml
> > ├── README.md
> > └── user
> >     ├── admin@example.com
> >     │   └── user-wildtype.fna
> >     └── admin@example.org
> >         └── user-wildtype.fna
> >
> > 4 directories, 5 files
> > ```
> {: .code-out}
{: .code-2col}

We have a directory named `admin`, which will be available to all admin users (we set `library_import_dir: /libraries/admin` earlier.)

Additionally we have a `user` directory, below the user directory are more directories with the user's email as they directory key. Data can be placed in here, and it will become accessible to those users (we set `user_library_import_dir: /libraries/user` for this.)

![An add datasets dropdown menu in galaxy showing the options from history, from user directory, and under admins only, from import directory](../../images/data/import-menu.png)

## from History

This is by far the easiest and most convenient option for small datasets, or datasets that are just already in a history

![A select box is shown listing files in the user's history](../../images/data/import-admin.png)

You can easily select multiple files and quickly import them.

## from User Directory

If user directories are configured, as we did at the beginning of this tutorial, then users will be able to import any files under their personal directory. This can be used for a wide variety of setups, e.g. providing the output of sequencing machines to users. This can point to the same directory structure that's used by the FTP service, if you want your users to be able to import files directly from FTP.

![Import popup with a list of files with one file, user-wildtype.fna, and buttons for configuring import behaviour.](../../images/data/import-user.png)

## from import Directory (Admins only)


![Same as the previous image, import popup listing options and one file, admin-wildtype.fna](../../images/data/import-admin.png)



select-from-lib0.png
select-from-lib1.png
select-from-lib2.png


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
