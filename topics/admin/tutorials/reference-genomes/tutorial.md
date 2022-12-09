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
  authorship:
    - slugger70
    - afgane
  editing:
    - hexylena
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

They are a flexible framework for adding reference data to Galaxy (not just genomic data). They are workflow compatible and can run via the Galaxy API. Data managers automatically update the appropriate location files when new data are installed.

More detailed background information on data managers can be found [here](https://galaxyproject.org/admin/tools/data-managers/), and details on how to define a data manager for a tool can be found [here](https://galaxyproject.org/admin/tools/data-managers/how-to/define/). There is also discussion of this topic in [the Reference Genomes in Galaxy slides]({% link topics/admin/tutorials/reference-genomes/slides.html %}).

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

# Galaxy Reference Data Components

As discussed in the overview, Galaxy Data Managers handle the population of numerous Galaxy components. From the highest level to the lowest:

- Addition of tool data table entries into the shed version of the tool data table configuration file, `shed_tool_data_table_conf.xml`. Multiple data tables with the same `name` attribute are automatically merged by Galaxy.
- Population of the corresponding location (`.loc`) file when new reference data are installed. The `.loc` file used corresponds to the Data Manager and its specific version that was executed to install the data. The path to the location file can be found in the Data Manager's entry in the tool data table configuration file (above).
- Generation of the reference data on disk, underneath the directory specified by `tool_data_path` in `galaxy.yml`.

# Install a Data Manager from the Tool Shed

> <comment-title>Galaxy Configuration</comment-title>
>
> No special configuration for Galaxy beyond creating an admin user is necessary in order to use Data Managers. However, you may want to be aware of the `enable_data_manager_user_view` option, which enables unprivileged user to browse data installed via Data Managers. This option is not required, and users will be able to *use* data installed via Data Managers regardless.
>
> Data are installed in the path specified by `tool_data_path` in `galaxy.yml`. If you prefer to keep hand-managed reference data separate from DM-managed reference data, you can set `galaxy_data_manager_data_path` instead.
>
{: .comment}

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
> > > ```xml
> > > <data_managers><data_manager id="fetch_genome_all_fasta_dbkeys" guid="toolshed.g2.bx.psu.edu/repos/devteam/data_manager_fetch_genome_dbkeys_all_fasta/data_manager/fetch_genome_all_fasta_dbkeys/0.0.1" shed_conf_file="/srv/galaxy/var/config/shed_tool_conf.xml">
> > >         <tool file="toolshed.g2.bx.psu.edu/repos/devteam/data_manager_fetch_genome_dbkeys_all_fasta/4d3eff1bc421/data_manager_fetch_genome_dbkeys_all_fasta/data_manager/data_manager_fetch_genome_all_fasta_dbkeys.xml" guid="toolshed.g2.bx.psu.edu/repos/devteam/data_manager_fetch_genome_dbkeys_all_fasta/data_manager_fetch_genome_all_fasta_dbkey/0.0.4"><tool_shed>toolshed.g2.bx.psu.edu</tool_shed><repository_name>data_manager_fetch_genome_dbkeys_all_fasta</repository_name><repository_owner>devteam</repository_owner><installed_changeset_revision>4d3eff1bc421</installed_changeset_revision><id>toolshed.g2.bx.psu.edu/repos/devteam/data_manager_fetch_genome_dbkeys_all_fasta/data_manager_fetch_genome_all_fasta_dbkey/0.0.4</id><version>0.0.4</version></tool><data_table name="all_fasta">
> > >             <output>
> > >                 <column name="value"/>
> > >                 <column name="dbkey"/>
> > >                 <column name="name"/>
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
> > >                 <column name="value"/>
> > >                 <column name="name"/>
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
> > > ```xml
> > > <?xml version="1.0" ?>
> > > <tables>
> > >     <table name="all_fasta" comment_char="#">
> > >         <columns>value, dbkey, name, path</columns>
> > >         <file path="/srv/galaxy/var/tool-data/toolshed.g2.bx.psu.edu/repos/devteam/data_manager_fetch_genome_dbkeys_all_fasta/4d3eff1bc421/all_fasta.loc"/>
> > >         <tool_shed_repository>
> > >             <tool_shed>toolshed.g2.bx.psu.edu</tool_shed>
> > >             <repository_name>data_manager_fetch_genome_dbkeys_all_fasta</repository_name>
> > >             <repository_owner>devteam</repository_owner>
> > >             <installed_changeset_revision>4d3eff1bc421</installed_changeset_revision>
> > >         </tool_shed_repository>
> > >     </table>
> > >     <table name="__dbkeys__" comment_char="#">
> > >         <columns>value, name, len_path</columns>
> > >         <file path="/srv/galaxy/var/tool-data/toolshed.g2.bx.psu.edu/repos/devteam/data_manager_fetch_genome_dbkeys_all_fasta/4d3eff1bc421/dbkeys.loc"/>
> > >         <tool_shed_repository>
> > >             <tool_shed>toolshed.g2.bx.psu.edu</tool_shed>
> > >             <repository_name>data_manager_fetch_genome_dbkeys_all_fasta</repository_name>
> > >             <repository_owner>devteam</repository_owner>
> > >             <installed_changeset_revision>4d3eff1bc421</installed_changeset_revision>
> > >         </tool_shed_repository>
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
>    ![nearly empty data manager tool list](../../images/dm-admin-local-data.png)
>
>    If you instead see the message "You do not currently have any Data Managers installed," you may need to restart your Galaxy server.
>
> 3. Click **all_fasta** under **View Tool Data Table Entries**
>
>    You should see the current contents of `tool-data/all_fasta.loc`, which will be empty.
>
> 4. Click **Local Data** to return to the previous page
>
> 5. Click {% tool [Create DBKey and Reference Genome - fetching]({{data_manager_fetch_genome_all_fasta_dbkey}}) %} under **Installed Data Managers**.
>    - {% icon param-select %} *"Use existing dbkey or create a new one"*: `Existing`
>    - {% icon param-select %} *"DBKEY to assign to data"*: `sacCer2`
>    - {% icon param-text %} *"Name of sequence"*: `S. cerevisiae June 2008 (SGD/sacCer2)`
>    - {% icon param-select %} *"Choose the source for the reference genome"*: `UCSC`
>
>    > <details-title>What do these fields mean?</details-title>
>    >
>    > The data manager tool form fields correspond to the columns in the tool data table, as referenced in `shed_tool_data_table_conf.xml`.
>    >
>    > Tool Field                | Tool Data Table Entry         | Purpose
>    > ---                       | -----                         | ---
>    > *ID for sequence*         | `value` column in `all_fasta` | Identifier of specific build of the genome, if it differs from the DBKEY. Used for variant builds such as the `hg19female` build of `hg19`.
>    > *DBKEY to assign to data* | `dbkey` column in `all_fasta` | Identifier of the genome, the UCSC build ID for UCSC builds like `hg19`.
>    > *Name of sequence*        | `name` column in `all_fasta`  | Display name of the genome, displayed in Galaxy tool dropdowns.
>    {: .details}
>
> 5. Click **Execute**. In your history, you will see a new dataset for the data manager run. When the job has finished, go back to the Data Manager view on the Galaxy Admin page (Click **Local Data**).
> 6. Click **all_fasta** under *View Tool Data Table Entries*
>
>    You should see that *sacCer2* has been added to all_fasta.
>
>    ![populated all_fasta data table](../../images/dm-all-fasta-populated.png)
>
{: .hands_on}

View in the file system where the changes to the location file and the reference genome fasta file.

> <question-title></question-title>
>
> What did this data manager execution change?
>
> > <solution-title></solution-title>
> >
> > - An entry for the *sacCer2* genome has been added to the `all_fasta.loc` file for the `data_manager_fetch_genome_all_fasta_dbkey` data manager. The path to the location file can be found in "shed" versions of the data manager config (`/srv/galaxy/var/config/shed_data_manager_conf.xml`) and data table config (`/srv/galaxy/var/config/shed_tool_data_table_conf.xml`)
> > - The data manager tool is installed along side other Galaxy tools in the shed tools directory
> >
> > > <code-in-title>Bash</code-in-title>
> > > Get the path to the `all_fasta.loc` file:
> > > ```bash
> > > grep 'data_manager_fetch_genome_dbkeys_all_fasta.*all_fasta.loc' /srv/galaxy/var/config/shed_tool_data_table_conf.xml
> > > ```
> > {: .code-in}
> >
> > > <code-out-title>Bash</code-out-title>
> > > ```xml
> > >         <file path="/srv/galaxy/var/tool-data/toolshed.g2.bx.psu.edu/repos/devteam/data_manager_fetch_genome_dbkeys_all_fasta/4d3eff1bc421/all_fasta.loc"/>
> > > ```
> > {: .code-out}
> >
> > > <code-in-title>Bash</code-in-title>
> > > Let's investigate the contents of the location file identified above (the path to yours may be different if you installed a different version of the data manager):
> > > ```bash
> > > cat /srv/galaxy/var/tool-data/toolshed.g2.bx.psu.edu/repos/devteam/data_manager_fetch_genome_dbkeys_all_fasta/4d3eff1bc421/all_fasta.loc
> > > ```
> > {: .code-in}
> >
> > > <code-out-title>Bash</code-out-title>
> > > ```bash
> > > #This file lists the locations and dbkeys of all the fasta files
> > > #under the "genome" directory (a directory that contains a directory
> > > #for each build). The script extract_fasta.py will generate the file
> > > #all_fasta.loc. This file has the format (white space characters are
> > > #TAB characters):
> > > #
> > > #<unique_build_id>  <dbkey>     <display_name>  <file_path>
> > > #
> > > #So, all_fasta.loc could look something like this:
> > > #
> > > #apiMel3    apiMel3 Honeybee (Apis mellifera): apiMel3      /path/to/genome/apiMel3/apiMel3.fa
> > > #hg19canon  hg19        Human (Homo sapiens): hg19 Canonical        /path/to/genome/hg19/hg19canon.fa
> > > #hg19full   hg19        Human (Homo sapiens): hg19 Full         /path/to/genome/hg19/hg19full.fa
> > > #
> > > #Your all_fasta.loc file should contain an entry for each individual
> > > #fasta file. So there will be multiple fasta files for each build,
> > > #such as with hg19 above.
> > > #
> > > sacCer2 sacCer2 S. cerevisiae June 2008 (SGD/sacCer2)   /srv/galaxy/var/tool-data/sacCer2/seq/sacCer2.fa
> > > ```
> > {: .code-out}
> >
> > > <code-in-title>Bash</code-in-title>
> > > Finally, verify that we have some sequences in the genome fasta file:
> > > ```bash
> > > grep '^>' /srv/galaxy/var/tool-data/sacCer2/seq/sacCer2.fa
> > > ```
> > {: .code-in}
> >
> > > <code-out-title>Bash</code-out-title>
> > > ```
> > > >2micron
> > > >chrI
> > > >chrII
> > > >chrIII
> > > >chrIV
> > > >chrIX
> > > >chrM
> > > >chrV
> > > >chrVI
> > > >chrVII
> > > >chrVIII
> > > >chrX
> > > >chrXI
> > > >chrXII
> > > >chrXIII
> > > >chrXIV
> > > >chrXV
> > > >chrXVI
> > > ```
> > {: .code-out}
> >
> {: .solution }
>
{: .question}

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

> <hands-on-title>Download and install sacCer2</hands-on-title>
>
> 1. Access the Admin menu from the top bar
> 2. Click **Local Data**, which can be found on the left, under **Server**
> 3. Click **{% tool [BWA-MEM index - builder]({{data_manager_bwa_mem_index_builder}}) %}** under **Installed Data Managers**.
>    - {% icon param-select %} *"Source Fasta Sequence"*: `S. cerevisiae June 2008 (SGD/sacCer2)`
> 4. Click **Execute**.
> 5. Verify that the new BWA index for sacCer2 has been built and the .loc file has been filled in. From the **Local Data** page in the Admin section, click on **bwa_mem_indexes** under *View Tool Data Table Entries*
>
>    S. cerevisiae sacCer2 should now appear in the list!
{: .hands_on}

> <question-title></question-title>
>
> What changes were made by the **BWA-MEM index builder**?
>
> > <solution-title></solution-title>
> >
> > - An entry for the *sacCer2* BWA/BWA-MEM index has been added to the `bwa_mem_index.loc` file for the `data_manager_bwa_mem_index_builder` data manager.
> > - A symlink to the sacCer2 genome fasta file installed by the Fetch Genome data manager above was created in the sacCer2 BWA index directory, `/srv/galaxy/var/tool-data/sacCer2/bwa_mem_index/sacCer2/`.
> > - The BWA/BWA-MEM index for the sacCer2 genome has been built.
> >
> > > <code-in-title>Bash</code-in-title>
> > >
> > > ```bash
> > > ls -l /srv/galaxy/var/tool-data/sacCer2/bwa_mem_index/sacCer2
> > > ```
> > {: .code-in}
> >
> > > <code-out-title>Bash</code-out-title>
> > > ```console
> > > total 20800
> > > lrwxrwxrwx 1 galaxy galaxy       20 Dec  9 18:29 sacCer2.fa -> ../../seq/sacCer2.fa
> > > -rw-r--r-- 1 galaxy galaxy       14 Dec  9 18:29 sacCer2.fa.amb
> > > -rw-r--r-- 1 galaxy galaxy      591 Dec  9 18:29 sacCer2.fa.ann
> > > -rw-r--r-- 1 galaxy galaxy 12163076 Dec  9 18:29 sacCer2.fa.bwt
> > > -rw-r--r-- 1 galaxy galaxy  3040750 Dec  9 18:29 sacCer2.fa.pac
> > > -rw-r--r-- 1 galaxy galaxy  6081552 Dec  9 18:29 sacCer2.fa.sa
> > > ```
> > {: .code-out}
> >
> {: .solution }
>
{: .question}

# Verify that BWA can access the new reference data

Finally, we will verify that the BWA tool can see the new sacCer2 indexes.

> <hands-on-title>Configure Galaxy</hands-on-title>
>
> 1. Install the BWA and BWA-MEM tools, if needed.
>
>    {% snippet topics/admin/faqs/install_tool.md query="bwa" name="Map with BWA" section="Mapping" %}
>
> 2. Click the "Home" icon at the top to return to the Galaxy analysis interface
>
> 3. Load the {% tool [Map with BWA](bwa) %} tool and verify that the sacCer2 genome appears in the {% icon param-select %} *"Using reference genome"* option.
{: .hands_on}

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
