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
  - shiltemann
subtopic: data
tags:
  - ansible
  - storage
  - git-gat
requirements:
 - type: "internal"
   topic_name: admin
   tutorials:
     - ansible
     - ansible-galaxy
     - tool-management
---

Data libraries are a great way to provide structured repositories of data to
everyone on your server, or to a subset of your users. Datasets in libraries do
not count against user quotas, so they are commonly used to provide the results
of a sequencing run or similar project to a group of users on the servers.

![screenshot of data libraries]({{site.baseurl}}/topics/admin/images/data/data-library-screenshot.png)


> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="data-library" %}

# Setup

Before we can import local data, we need to configure Galaxy to permit this. Additionally we will setup an example data library which we can use for demonstrative purposes.

> <hands-on-title>Setting up Data Libraries</hands-on-title>
>
> 1. We will add a pre-task to clone [a data repository](https://github.com/usegalaxy-eu/libraries-training-repo) into your machine. We will use this as the source for a library dataset.
>
>    {% raw %}
>    ```diff
>    --- a/galaxy.yml
>    +++ b/galaxy.yml
>    @@ -8,6 +8,9 @@
>         - name: Install Dependencies
>           package:
>             name: ['acl', 'bzip2', 'git', 'make', 'python3-psycopg2', 'tar', 'virtualenv']
>    +    - git:
>    +        repo: 'https://github.com/usegalaxy-eu/libraries-training-repo'
>    +        dest: /libraries/
>       roles:
>         - galaxyproject.postgresql
>         - role: galaxyproject.postgresql_objects
>    {% endraw %}
>    ```
>    {: data-commit="Add the git repository to the pre-tasks"}
>
>    {% snippet topics/admin/faqs/diffs.md %}
>
> 2. Take a minute to explore the [folders in our sample library.](https://github.com/usegalaxy-eu/libraries-training-repo). These will be important when we start loading data.
>
> 4. Edit the file `group_vars/galaxyservers.yml` and set the following variables:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -29,6 +29,8 @@ miniconda_manage_dependencies: false
>     
>     galaxy_config:
>       galaxy:
>    +    library_import_dir: /libraries/admin
>    +    user_library_import_dir: /libraries/user
>         dependency_resolvers_config_file: "{{ galaxy_config_dir }}/dependency_resolvers_conf.xml"
>         containers_resolvers_config_file: "{{ galaxy_config_dir }}/container_resolvers_conf.xml"
>         tool_data_table_config_path: /cvmfs/data.galaxyproject.org/byhand/location/tool_data_table_conf.xml,/cvmfs/data.galaxyproject.org/managed/location/tool_data_table_conf.xml
>    {% endraw %}
>    ```
>    {: data-commit="Configure the library import directories"}
>
>    Note that the `/libraries/admin` will refer to a folder within the libraries-training-repo that we cloned in the pre-task, and likewise for `/libraries/user`
>
> 5. Run the playbook:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
{: .hands_on}

# Importing Data

There are multiple options for importing data from your server, we'll go through all of your choices below. But first, let's take a quick look at the example library structure we've provided.

> > <code-in-title>Bash</code-in-title>
> > ```bash
> > tree /libraries
> > ```
> {: .code-in}
>
> > <code-out-title>Bash</code-out-title>
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

> <tip-title>admin@example.org</tip-title>
> Note that in the user directories, admin@example.com and example.org are used, if you've used a different email address for *your* admin user, you'll need to copy one of these directories.
{: .tip}

> ```bash
> 1.sh
> ```
> {: data-test="true"}
{: .hidden}

We have a directory named `admin`, which will be available to all admin users (we set `library_import_dir: /libraries/admin` earlier.)

Additionally we have a `user` directory, below the user directory are more directories with the user's email as they directory key. Data can be placed in here, and it will become accessible to those users (we set `user_library_import_dir: /libraries/user` for this.)

![An add datasets dropdown menu in galaxy showing the options from history, from user directory, and under admins only, from import directory](../../images/data/import-menu.png)

## from History

This is by far the easiest and most convenient option for small datasets, or datasets that are just already in a history

![A select box is shown listing files in the user's history](../../images/data/import-history.png)

You can easily select multiple files and quickly import them.

## from User Directory

If user directories are configured, as we did at the beginning of this tutorial, then users will be able to import any files under their personal directory. This can be used for a wide variety of setups, e.g. providing the output of sequencing machines to users. This can point to the same directory structure that's used by the FTP service, if you want your users to be able to import files directly from FTP.

![Import popup with a list of files with one file, user-wildtype.fna, and buttons for configuring import behaviour.](../../images/data/import-user.png)

This will enable the option for everyone, any unprivileged user with a folder in the `user_library_import_dir` directory *and* permissions on a data library can import things from their import directory.

> <tip-title>Do you not see a directory here?</tip-title>
> We pre-created a directory for someone registered as `admin@example.org`, but if you are logged in with an account registered with a different email, you'll see nothing.
>
> You can fix this by going into /libraries/user and `cp -Rv admin@example.org whatever-email-you-used`
> Then it should appear in Galaxy.
{: .tip}

## from import Directory (Admins only)

Similarly to the user import directories, there is the option to have an admin only import area. If one admin is responsible for creating the data libraries and importing data this can be a good option.

![Same as the previous image, import popup listing options and one file, admin-wildtype.fna](../../images/data/import-admin.png)

An important feature of data libraries importing is the ability to "link files instead of copying". If you know your data will not be deleted (e.g. sequencing data sitting on an NFS), then you can use this option to further decrease the amount of storage your Galaxy server requires. The datasets can be used as if they were normal datasets, imported by users and analysed, but not imported into the Galaxy data storage area.

# Automatically Populating a Data Library

If your data is accessible via URL, you can write a yaml file to import and setup the data library automatically.
We've included this file in the example data libraries repository we cloned at the beginning:

> <code-in-title>Bash</code-in-title>
> ```bash
> cat /libraries/example-library.yaml
> ```
{: .code-in}

> <code-out-title>Bash</code-out-title>
> ```yaml
> ---
> destination:
>   type: library
>   name: Mouse sequencing project
>   description: some data
>   synopsis: samples collected from somewhere
> items:
> - url: https://zenodo.org/api/files/287966da-5411-4f79-8cfb-0ffa84d0d6cc/wildtype.fna
>   src: url
>   ext: fasta
>   info: https://doi.org/10.5281/zenodo.582600
> - name: A directory
>   description: Exome sequencing means that all protein-coding genes in a genome are
>   items:
>   - url: https://zenodo.org/api/files/287966da-5411-4f79-8cfb-0ffa84d0d6cc/mutant_R1.fastq
>     src: url
>     ext: fastqsanger
>     info: https://doi.org/10.5281/zenodo.582600
>   - url: https://zenodo.org/api/files/287966da-5411-4f79-8cfb-0ffa84d0d6cc/mutant_R2.fastq
>     src: url
>     ext: fastqsanger
>     info: https://doi.org/10.5281/zenodo.582600
> ```
{: .code-out}

Let's try setting that up in our Galaxy!

{% snippet faqs/galaxy/preferences_admin_api_key.md admin=true %}

> <hands-on-title>hands-on: importing data library yaml.</hands-on-title>
>
> 1. re-activate the virtualenv you created for the [ephemeris tool management tutorial]({% link topics/admin/tutorials/tool-management/tutorial.md %}).
>
>    > <code-in-title>input: bash</code-in-title>
>    > ```
>    > . ~/ephemeris_venv/bin/activate
>    > ```
>    {: .code-in}
>
>    > <tip-title>missing?</tip-title>
>    >
>    > then you might need to re-run the steps:
>    >
>    > ```bash
>    > virtualenv -p python3 ~/ephemeris_venv
>    > . ~/ephemeris_venv/bin/activate
>    > pip install ephemeris
>    > ```
>    {: .tip}
>
> 2. we'll use the `setup-data-libraries` command to install the data in this yaml file into a library in our galaxy.
>
>    > <code-in-title>input: bash</code-in-title>
>    > ```bash
>    > setup-data-libraries -g https://your-galaxy -a <api-key> --training -i /libraries/example-library.yaml --legacy
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
>    > <code-out-title>output</code-out-title>
>    > ```
>    > library name: mouse sequencing project
>    > ```
>    {: .code-out}
>
> 3. this command is safe to re-run. for example if you update the yaml, it will simply report that the library exists before ensuring all files exist in their appropriate location:
>
>    ```bash
>    library name: mouse sequencing project
>    library already exists! id: f597429621d6eb2b
>    ```
>
{: .hands_on}

![Screenshot of data libraries, we're in a library folder named "Mouse sequencing project" and a directory and single file are shown. The file has an ugly URL as its name](../../images/data/imported.png)

That's it! You should be able to see your newly created data library in your Galaxy.

> ```bash
> 2.sh
> ```
> {: data-test="true"}
{: .hidden}

{% snippet topics/admin/faqs/missed-something.md step=6 %}

Note that we've used some special flags here, `--training` and `--legacy`. Training sets some defaults that make sense for the GTN (mostly around library descriptions / etc.)

### `--legacy`

This enables the use of legacy APIs even for newer Galaxies that should have a batch upload API enabled. Unfortunately the new batch upload API is not able to update existing data libraries in place and will always create new libraries with the same name. So currently `--legacy` is quite useful for maintaining a YAML file, and running `setup-data-libraries` regularly whenever that file has updates.

But this comes with the downside that the entire URL is used for the filename.

# GTN Data

This process scales quite well. Galaxy Europe, in their mission to support all of the GTN materials on their server setup [a shared-data repository](https://github.com/usegalaxy-eu/shared-data), a single giant YAML file with all of the data from the GTN tutorials. This was then expanded to other UseGalaxy.* servers to provide synced training libraries for all of the major servers.

Do you want the training data synced on your server? If so join the shared data repository! If you provide us with a non-admin API key and a data library to upload data into, we can sync this automatically.

# Using Data from Libraries

Users can now conveniently use datasets in libraries when they are running analyses, without having to import them first.

The tool form provides a button on the right of every dataset selector which allows users to "Browse datasets"

![tool form with the folder icon to the right of the dataset selector shown, the button with tooltip "Browse datasets" is highlighted](../../images/data/select-from-lib0.png)

Users can then choose from datasets in their history, or browse through the data libraries

![Popup with "Data Libraries", datasets from histories, and an upload button](../../images/data/select-from-lib1.png)

Here users can see every file accessible to them in the data library to begin analysing.

![Same popup as previous image, but the file listing now shows the contents of the data library we created during the hands-on step with Mouse data](../../images/data/select-from-lib2.png)
