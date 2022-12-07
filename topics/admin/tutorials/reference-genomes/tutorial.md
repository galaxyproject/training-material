---
layout: tutorial_hands_on

title: "Reference Data with Data Managers"
zenodo_link: ""
questions:
objectives:
  - Have an understanding of the way in which Galaxy stores and uses reference data
  - Be able to download and use data managers to add a reference genome and its pre-calculated indices into the Galaxy reference data system
  - Use an Ansible playbook for all of the above
time_estimation: "1h"
key_points:
  - Understand how Galaxy stores and uses its reference data
  - Understand how to manually add a reference genome and tool indices if required
  - Understand and how to use data managers to make all of this much much easier
contributors:
  - Slugger70
  - afgane
  - natefoo
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
subtopic: data
---

# Overview

**The problem**

The Galaxy server administrator needed to know how to update each type of reference data, how to run the indexers, where to get the data from, and how to update Galaxy's *Tool Data Table* and *location* configuration files.

**Data managers to the rescue**

Data Managers are a special class of Galaxy tool which allows for the download and/or creation of data that is stored within Tool Data Tables and their underlying flat location (e.g. .loc) files. These tools handle the creation of indices and the addition of entries/lines to the data table / .loc file via the Galaxy admin interface.

Data Managers can be defined locally in the data manager and tool data table configuration files or installed through the Tool Shed. When Data Managers are installed from the Tool Shed, their configuration is added to the shed versions of the data manager configuration and tool data table configuration files.

They are a flexible framework for adding reference data to Galaxy (not just genomic data). They are workflow compatible and can run via the Galaxy API.


- Data Manager configuration (e.g. `data_manager_conf.xml`)
- Data Manager tool

Data managers automatically update the appropriate location files when new data are installed.

More detailed background information on data managers can be found [here](https://galaxyproject.org/admin/tools/data-managers/), and details on how to define a data manager for a tool can be found [here](https://galaxyproject.org/admin/tools/data-managers/how-to/define/).

Data Managers are composed of two components:
> <comment-title>Pre-built data are available</comment-title>
>
> The usegalaxy.\* servers and Galaxy Community have a large amount of reference data online and available for use by your Galaxy server. For instructions on how to access and use these data, see [the Reference Data with CVMFS tutorial]({% link topics/admin/tutorials/cvmfs/tutorial.md %}).
>
> If your data are not available as part of the CVMFS repository, Galaxy Data Managers can be used to locally install and build reference data.
>
{: .comment}

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data Manager Configuration

**FIXME: is any of this necessary?**

We need to tell Galaxy where to find the Data Managers and their configuration.

> <hands-on-title>Configuring Galaxy to use Data Managers</hands-on-title>
>
> 1. In your working directory, edit the `group_vars/galaxyservers.yml` file.
>
>    If the terms "Ansible", "role" and "playbook" mean nothing to you, please checkout [the Ansible introduction slides]({% link topics/admin/tutorials/ansible/slides.html %}) and [the Ansible introduction tutorial]({% link topics/admin/tutorials/ansible/tutorial.md %})
> 
> 1. **Comment out existing definition of `tool_data_table_config_path`**
>
> 2. Add:
>
>    ```yaml
>    galaxy_config:
>      galaxy:
>        enable_data_manager_user_view: "True"
>        galaxy_data_manager_data_path: "{{ galaxy_mutable_data_dir }}/tool-data"
>    ```
{: .hands_on}

Where:
- *enable_data_manager_user_view* allows non-admin users to view the available data that has been managed.
- *galaxy_data_manager_data_path* defines the location to use for storing the files created by Data Managers. When not configured it defaults to the value of `tool_data_path`.

# Install a Data Manager from the Tool Shed

We will install a data manager that can fetch the various genome sequences from multiple sources.

> <hands-on-title>Install the "Fetch Genome" Data Manager </hands-on-title>
>
> > <tip-title>Intall with Ephemeris</tip-title>
> > This hands-on exercise installs a tool through the Galaxy UI, but you are encouraged to install tools in a deterministic, recordable way through the use of Ephemeris, which is descirbed in [the Galaxy Tool Management with Ephemeristutorial]({ link topics/admin/tutorials/tool-management/tutorial.md %}).
> {: .tip}
>
> 1. Access the **Admin** menu from the top bar (you need to be logged-in with an email specified in the `admin_users` setting)
> 2. Click **Install and Uninstall**, which can be found on the left, under **Tool Management**
> 3. Enter `fetch_genome` in the search interface
> 4. Click on the first hit, having `devteam` as owner
> 5. Click the **Install** button for the latest revision
{: .hands_on}

View in the file system where the various elements land. Have a look in the configuration files located in config directory.

> <question-title></question-title>
>
> What did this tool installation change?
>
> > <solution-title></solution-title>
> >
> > - The data manager and its data tables are added to the Galaxy-managed "shed" versions of the data manager config (`/srv/galaxy/var/config/shed_data_manager_conf.xml`) and data table config (`/srv/galaxy/var/config/shed_tool_data_table_conf.xml`)
> > - The data manager tool is installed along side other Galaxy tools in the shed tools directory
> >
> > > <code-in-title>Bash</code-in-title>
> > > Let's investigate the data manager config file.
> > > ```bash
> > > cat /srv/galaxy/var/config/shed_data_manager_conf.xml
> > > ```
> > {: .code-in}
> >
> > > <code-out-title>Bash</code-out-title>
> > > ``` xml
> > > <?xml version="1.0"?>
> > > <data_managers>
> > >     <data_manager guid="toolshed.g2.bx.psu.edu/repos/devteam/data_manager_fetch_genome_dbkeys_all_fasta/data_manager/fetch_genome_all_fasta_dbkeys/0.0.1" id="fetch_genome_all_fasta_dbkeys" shed_conf_file="./config/shed_tool_conf.xml">
> > >         <tool file="toolshed.g2.bx.psu.edu/repos/devteam/data_manager_fetch_genome_dbkeys_all_fasta/b1bc53e9bbc5/data_manager_fetch_genome_dbkeys_all_fasta/data_manager/data_manager_fetch_genome_all_fasta_dbkeys.xml" guid="toolshed.g2.bx.psu.edu/repos/devteam/data_manager_fetch_genome_dbkeys_all_fasta/data_manager_fetch_genome_all_fasta_dbkey/0.0.2"><tool_shed>toolshed.g2.bx.psu.edu</tool_shed><repository_name>data_manager_fetch_genome_dbkeys_all_fasta</repository_name><repository_owner>devteam</repository_owner><installed_changeset_revision>b1bc53e9bbc5</installed_changeset_revision><id>toolshed.g2.bx.psu.edu/repos/devteam/data_manager_fetch_genome_dbkeys_all_fasta/data_manager_fetch_genome_all_fasta_dbkey/0.0.2</id><version>0.0.2</version></tool><data_table name="all_fasta">
> > >             <output>
> > >                 <column name="value" />
> > >                 <column name="dbkey" />
> > >                 <column name="name" />
> > >                 <column name="path" output_ref="out_file">
> > >                     <move type="file">
> > >                         <source>${path}</source>
> > >                         <target base="${GALAXY_DATA_MANAGER_DATA_PATH}">${dbkey}/seq/${path}</target>
> > >                     </move>
> > >                     <value_translation>${GALAXY_DATA_MANAGER_DATA_PATH}/${dbkey}/seq/${path}</value_translation>
> > >                     <value_translation type="function">abspath</value_translation>
> > >                 </column>
> > >             </output>
> > >         </data_table>
> > >         <data_table name="__dbkeys__">
> > >             <output>
> > >                 <column name="value" />
> > >                 <column name="name" />
> > >                 <column name="len_path" output_ref="out_file">
> > >                     <move type="file">
> > >                         <source>${len_path}</source>
> > >                         <target base="${GALAXY_DATA_MANAGER_DATA_PATH}">${value}/len/${len_path}</target>
> > >                     </move>
> > >                     <value_translation>${GALAXY_DATA_MANAGER_DATA_PATH}/${value}/len/${len_path}</value_translation>
> > >                     <value_translation type="function">abspath</value_translation>
> > >                 </column>
> > >             </output>
> > >         </data_table>
> > >     </data_manager>
> > > </data_managers>
> > > ```
> > {: .code-out.code-max-300}
> >
> > > <code-in-title>Bash</code-in-title>
> > > Let's also investigate the tool data table config file.
> > > ```bash
> > > cat /srv/galaxy/var/config/shed_tool_data_table_conf.xml
> > > ```
> > {: .code-in}
> >
> > > <code-out-title>Bash</code-out-title>
> > > ``` xml
> > > <?xml version="1.0"?>
> > > <tables>
> > > <table comment_char="#" name="all_fasta">
> > >         <columns>value, dbkey, name, path</columns>
> > >         <file path="/srv/galaxy/tool-data/toolshed.g2.bx.psu.edu/repos/devteam/data_manager_fetch_genome_dbkeys_all_fasta/b1bc53e9bbc5/all_fasta.loc" />
> > >         <tool_shed_repository>
> > >             <tool_shed>toolshed.g2.bx.psu.edu</tool_shed>
> > >             <repository_name>data_manager_fetch_genome_dbkeys_all_fasta</repository_name>
> > >             <repository_owner>devteam</repository_owner>
> > >             <installed_changeset_revision>b1bc53e9bbc5</installed_changeset_revision>
> > >             </tool_shed_repository>
> > >     </table>
> > > <table comment_char="#" name="__dbkeys__">
> > >         <columns>value, name, len_path</columns>
> > >         <file path="/srv/galaxy/tool-data/toolshed.g2.bx.psu.edu/repos/devteam/data_manager_fetch_genome_dbkeys_all_fasta/b1bc53e9bbc5/dbkeys.loc" />
> > >         <tool_shed_repository>
> > >             <tool_shed>toolshed.g2.bx.psu.edu</tool_shed>
> > >             <repository_name>data_manager_fetch_genome_dbkeys_all_fasta</repository_name>
> > >             <repository_owner>devteam</repository_owner>
> > >             <installed_changeset_revision>b1bc53e9bbc5</installed_changeset_revision>
> > >             </tool_shed_repository>
> > >     </table>
> > > </tables>
> > > ```
> > {: .code-out.code-max-300}
> >
> {: .solution }
>
{: .question}

# Download and install a reference genome sequence

Next, we will install some reference data. Specifically, we will grab sacCer2 (version 2 of the Saccharomyces cerevisiae genome).

> <hands-on-title>Download and install sacCer2</hands-on-title>
>
> 1. Access the Admin menu from the top bar
> 2. Click **Local Data**, which can be found on the left, under **Server**
>
>    You should see something like this:
>
>    **TODO: recapture and rehost** ![nearly empty data manager tool list](https://github.com/galaxyproject/admin-training/raw/2018-oslo/docs/05-reference-genomes/images/nearly_empty_data_manager_tool_list.png)
>
> 3. Click **all_fasta** under **View Tool Data Table Entries**
>
>    You should see the current contents of `tool-data/all_fasta.loc`, which will be empty.
>
> 4. Under **Run Data Manager Tools**, click {% tool [Create DBKey and Reference Genome - fetching]({{data_manager_fetch_genome_all_fasta_dbkey}}) %}. NOTE: If you receive the error "Uncaught exception in exposed API method:", you will need to restart Galaxy first.
>    - {% icon param-select %} *"DBKEY to assign to data"*: `sacCer2`
>    - {% icon param-text %} *"Name of sequence"*: `S. cerevisiae June 2008 (SGD/sacCer2) (sacCer2)``
>    - {% icon param-text %} *"ID for sequence"*: leave this field blank
>
>    **TODO: what do these fields mean and how determine their values?**
>
> 5. Click **Execute**. In your history, you will see a new dataset for the data manager run. When the job has finished, go back to the Data Manager view on the Galaxy Admin page (Click **Local Data**).
> 6. Click **all_fasta** under *View Tool Data Table Entries*
>
>    You should see that sacCer2 has been added to all_fasta.
>
>    **TODO: recapture and rehost** ![all_fasta.png](https://github.com/galaxyproject/admin-training/raw/2018-oslo/docs/05-reference-genomes/images/all_fasta.png)
>
> **TODO: identify loc files and inspect changes**
{: .hands_on}

# Download and install the BWA data manager

Having the genome is a prerequisite for our ultimate goal, which is to use the sacCer2 genome as a *reference genome* for the BWA tool. BWA, like many tools, needs an *index* of the reference genome, and has its own format for that index. Thankfully, the BWA/BWA-MEM data manager will build that index for us.

In this part we will repeat the same process as when we installed the Fetch Genome data manager, except that we will install the BWA/BWA-MEM data manager this time.

> <hands-on-title>Install the "Fetch Genome" Data Manager </hands-on-title>
>
> 1. Access the **Admin** menu from the top bar
> 2. Click **Install and Uninstall**, which can be found on the left, under **Tool Management**
> 3. Enter `bwa_mem_index` in the search interface
> 4. Click on the first hit, having `devteam` as owner
> 5. Click the **Install** button for the latest revision
{: .hands_on}

# Build the BWA index for sacCer2

In this part we will actually build the BWA index for sacCer2. It will automatically be added to our list of available reference genomes in the BWA tool.

* From the Galaxy Admin page, click **Local data**
* Click on **BWA-MEM index - builder** under *Run Data Manager Tools*
  * Select *sacCer2* for Source Fasta Sequence
  * Put sacCer2 into the other two blank fields.
  * Click **Execute**. NOTE: If you receive the error "Parameter all_fasta_source requires a value, but has no legal values defined.", you will need to restart Galaxy first.

> <hands-on-title>Download and install sacCer2</hands-on-title>
>
> 1. Access the Admin menu from the top bar
> 2. Click **Local Data**, which can be found on the left, under **Server**
> 3. Under **Run Data Manager Tools**, click {% tool [BWA-MEM index - builder]({{data_manager_bwa_mem_index_builder}}) %}. NOTE: If you receive the error "Uncaught exception in exposed API method:", you will need to restart Galaxy first.
>    - {% icon param-select %} *"Source Fasta Sequence"*: `sacCer2`
>    - {% icon param-text %} *"FIXME: unknown field 1"*: `sacCer2`
>    - {% icon param-text %} *"FIXME: unknown field 2*: `sacCer2`
>
>    **TODO: what do these fields mean and how determine their values?**
> 4. Click **Execute**. NOTE: If you receive the error "Parameter all_fasta_source requires a value, but has no legal values defined.", you will need to restart Galaxy first.
> 5. Verify that the new BWA index for sacCer2 has been built and the .loc file has been filled in. From the **Local Data** page in the Admin section, click on **bwa mem indexes** under *View Tool Data Table Entries*
>
>    S. cerevisiae sacCer2 should now appear in the list!
{: .hands_on}

# Run BWA with the new reference data

Finally, we will run the BWA tool and check to see if the reference data is listed and the tool works with it.

- **TODO: hands-on**: Run the BWA tool on the 2 fast files we loaded earlier, using sacCer2 as the reference.

How cool is that? No editing `.loc` files, no making sure you've got TABs instead of spaces. Fully auto!

# Addendum: Installing and Running a DM with ephemeris

Create a config file for `run-data-managers` named `fetch-sacCer3.yml`:

```yaml
data_managers:
    # Data manager ID
    - id: toolshed.g2.bx.psu.edu/repos/devteam/data_manager_fetch_genome_dbkeys_all_fasta/data_manager_fetch_genome_all_fasta_dbkey/0.0.4
      # tool parameters, nested parameters should be specified using a pipe (|)
      params:
        - 'dbkey_source|dbkey': '{{ item }}'
        - 'reference_source|reference_source_selector': 'ucsc'
        - 'reference_source|requested_dbkey': '{{ item }}'
      # Items refere to a list of variables you want to run this data manager. You can use them inside the param field with {{ item }}
      # In case of genome for example you can run this DM with multiple genomes, or you could give multiple URLs.
      items:
        - sacCer3
        #- dm3
        #- mm9
      # Name of the data-tables you want to reload after your DM are finished. This can be important for subsequent data managers
      data_table_reload:
        - all_fasta
        - __dbkeys__
```

Install the fetch DM:

```console
$ shed-tools install -g https://gat-0.student.galaxy.training --name data_manager_fetch_genome_dbkeys_all_fasta --owner devteam -a abbacadabba
Storing log file in: /tmp/ephemeris_2qpg_hrq
(1/1) repository data_manager_fetch_genome_dbkeys_all_fasta already installed at revision 4d3eff1bc421. Skipping.
Installed repositories (0): []
Skipped repositories (1): [('data_manager_fetch_genome_dbkeys_all_fasta', '4d3eff1bc421')]
Errored repositories (0): []
All repositories have been installed.
Total run time: 0:00:00.770248
```

Run the fetch DM:

> <code-in-title></code-in-title>
> ```
> run-data-managers --config fetch-sacCer3.yml -g https://gat-0.student.galaxy.training -a abbacadabba
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> Storing log file in: /tmp/ephemeris_kfsmjk2a
> Running data managers that populate the following source data tables: ['all_fasta']
> Dispatched job 1. Running DM: "toolshed.g2.bx.psu.edu/repos/devteam/data_manager_fetch_genome_dbkeys_all_fasta/data_manager_fetch_genome_all_fasta_dbkey/0.0.4" with parameters: {'dbkey_source|dbkey': 'sacCer3', 'reference_source|reference_source_selector': 'ucsc', 'reference_source|requested_dbkey': 'sacCer3'}
> Job 1 finished with state ok.
> Running data managers that index sequences.
> Finished running data managers. Results:
> Successful jobs: 1
> Skipped jobs: 0
> Failed jobs: 0
> ```
