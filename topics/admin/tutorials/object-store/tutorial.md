---
layout: tutorial_hands_on

title: "Distributed Object Storage"
questions:
- How does Galaxy locate data?
- How can I have Galaxy use multiple storage locations?
objectives:
- Setup Galaxy with both the Hierarachical and Distributed Object Storages
time_estimation: "30m"
key_points:
- The distributed object store configuration allows you to easily expand that storage that is attached to your Galaxy.
- You can move data around without affecting users.
contributors:
  - natefoo
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

# Expanding Storage

{:.no_toc}

You may find that your Galaxy files directory has run out of space, but you don't want to move all of the files from one filesystem to another. One solution to this problem is to use Galaxy's hierarchical object store to add an additional file space for Galaxy.

Alternatively, you may wish to write new datasets to more than one filesystem. For this, you can use Galaxy's distributed object store.

This tutorial assumes you have done the "Ansible for installing Galaxy" tutorial, it references the base configuration set up in that tutorial in numerous places.

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Hierarchical Object Store

First, note that your Galaxy datasets have been created thus far in the directory `/data`, due to `galaxy_config: galaxy: file_path`. In some cases, we run out of storage in a particular location. Galaxy allows us to add additional storage locations where it will create new datasets, while still looking in the old locations for old datasets. You will not have to migrate any of your datasets, and can just "plug and play" with new storage pools.


> ### {% icon hands_on %} Hands-on: Adding Hierarchical Storage
>
> 1. Open your group variables file and set the `object_store_config_file` variable:
>
>    ```yaml
>    galaxy_config:
>      galaxy:
>        object_store_config_file: {% raw %}"{{ galaxy_config_dir }}/object_store_conf.xml"{% endraw %}
>    ```
>
> 2. In your group variables file, add it to the `galaxy_config_files` section:
>
>    ```yaml
>    galaxy_config_files:
>      - src: files/galaxy/config/object_store_conf.xml
>        dest: {% raw %}"{{ galaxy_config.galaxy.object_store_config_file }}"{% endraw %}
>    ```
>
> 3. Create and edit `files/galaxy/config/object_store_conf.xml` with the following contents:
>
>    ```xml
>    <?xml version="1.0"?>
>    <object_store type="hierarchical">
>        <backends>
>            <backend id="newdata" type="disk" order="0">
>                <files_dir path="/data2" />
>                <extra_dir type="job_work" path="/data2/job_work_dir" />
>            </backend>
>            <backend id="olddata" type="disk" order="1">
>                <files_dir path="/data" />
>                <extra_dir type="job_work" path="/data/job_work_dir" />
>            </backend>
>        </backends>
>    </object_store>
>    ```
>
> 4. Add a `pre_task` to create the `/data2` folder [using the file module](https://docs.ansible.com/ansible/latest/modules/file_module.html), exactly like for the `/data` folder.
>
> 5. Run the playbook and restart Galaxy
>
> 6. Run a couple of jobs after Galaxy has restarted, run a couple of jobs.
>
>    > ### {% icon question %} Question
>    >
>    > Where is the data now stored?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > You should see `/data2` in the `Full Path`, if not, something went wrong. Check that your "order" is correct
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
{: .hands_on}

# Distributed Object Store

Rather than searching a hierarchy of object stores until the dataset is found, Galaxy can store the ID (in the database) of the object store in which a dataset is located when the dataset is created. This allows Galaxy to write to more than one object store for new datasets.


> ### {% icon hands_on %} Hands-on: Distributed Object Store
>
> 1. Edit your `files/galaxy/config/object_store_conf.xml` file and replace the contents with:
>
>    ```xml
>    <?xml version="1.0"?>
>    <object_store type="distributed">
>        <backends>
>            <backend id="newdata" type="disk" weight="1">
>                <files_dir path="/data2"/>
>                <extra_dir type="job_work" path="/data2/job_work_dir"/>
>            </backend>
>            <backend id="olddata" type="disk" weight="1">
>                <files_dir path="/data"/>
>                <extra_dir type="job_work" path="/data/job_work_dir"/>
>            </backend>
>        </backends>
>    </object_store>
>    ```
>
> 2. Run the playbook, restart Galaxy
>
> 3. Run 4 or so jobs, and check where the output appear. You should see that they are split relatively evenly between the two data directories.
>
{: .hands_on}

Sites like UseGalaxy.eu use the distributed object store in order to balance dataset storage across 10 different storage pools.

> ### {% icon details %} More documentation
>
> More information can be found in the [sample file](https://github.com/galaxyproject/galaxy/blob/dev/config/object_store_conf.xml.sample).
>
{: .details}

> ### {% icon warning %} Warning: switching object store types will cause issues
> We have switched between two different object stores here, but this is not supported. If you need to do this, you will need to update datasets in Galaxy's database. Any datasets that were created as hierarchical will lack the `object_store_id`, and you will need to supply the correct one. Do not just blindly copy these instructions, please understand what they do before running them and talk to us on [Gitter](http://gitter.im/galaxyproject/Lobby) for more help
>
> 1. Move the datasets to their new location: `sudo -u galaxy rsync -avr /hierarchical/000/ /distributed/000/`
>
> 2. Update the database: `sudo -Hu galaxy psql galaxy -c "UPDATE dataset SET object_store_id='data';" `
>
> 3. Restart your Galaxy
>
{: .warning}
