---
layout: tutorial_hands_on

title: Understanding Galaxy Data Managers
level: Introductory
subtopic: tooldev
questions:
- What are Data Managers and why are they needed?
- How do you write a Data Manager Tool?
- How do you test a Data Manager Tool?
objectives:
- Understand the idea behind and the concept of Galaxy Data Managers
- Understand what components are needed to write your own Data Manager Tool
- Know how to test a Data Manager and be aware of the limitations of the current test framework with respect to Data Managers
time_estimation: 2H
key_points:
- Data Managers are tools to be run by admins of a Galaxy instance.
- They automate reference data collection and preparation and they write Data Table (.loc file) records.
- In addition to a regular tool wrapper xml file, Data Managers require several xml config files that define the interaction of the Data Manager with the Galaxy framework and with the Data Tables they are supposed to populate.
- Currently, Galaxy only allows for partial automation of Data Manager testing. Some manual testing is required.
contributions:
  authorship:
    - wm75

---

# What are Data Managers and why are they needed?

Many tools run with two kinds of input data: some experimental data specific to
the tool run (like, e.g., sequencing data) and other data, which stays the
same across a range of different tool runs (e.g. a reference genome, or some
genome annotations that are the same across runs for the same organism).

Forcing users to provide that second type of data for every tool run is undesirable
because:

- some of that data is complicated to gather from public sources or needs some preprocessing
- it leads to unnecessary copies of data (of often very significant size) that would better be reused across user accounts

One possible solution (which was actually used in the early days of Galaxy) is to have
Galaxy server admins collect and prepare commonly used data, store it on the server and
record the location along with other metadata in a simple tab-separated so-called .loc file.
Tools can then declare select parameters that are populated with the records of specific .loc files,
can access the paths stored in them and, at tool run time, retrieve the data cached on the server.

While this is user-friendly, it shifts the burden of making cached data available to admins.
With more and more tools requiring data of very different formats,
this approach becomes increasingly unmanageable because, for each tool,
an admin has to research where to obtain the data, or how to calculate it,
and whether it needs some reformatting or other pre-processing before being usable by tools.

![The issue of maintaining different kinds of data with the .loc files approach illustrated through three examples](../../images/loc-files-approach.png "The issue of maintaining different kinds of data with the .loc files approach")

In principle, admins could automate some of this work through scripts, but it would be nice to not have each admin reinvent the wheel.

## The idea behind Data Managers

Just like admin-installed data found via .loc files frees the *user* from having to care about all the details of how to obtain the data,
special tools, ideally written by people who know how to obtain and set up certain types of data,
should automate data collection and preparation, and the writing of .loc file records, and make life easier for *admins*.
In other words, admins become users themselves. They don't have to know all the details, but just run a tool that "knows" how to install a certain type of data on the server.

> <comment-title>Reference document</comment-title>
>
> This tutorial is an attempt to describe the different parts and functions of Data Managers in
> a way that is structured as logically as possible.
>
> An in-depth, technical explanation of this matter is provided at
> <https://docs.galaxyproject.org/en/latest/dev/data_managers.html>
> and, when in doubt, that material should be considered the reference document
> for Data Managers.
>
{: .comment}

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Components of a Data Manager

Here you see a tree view of the files that together constitute the widely used [bowtie2 Data Manager](https://github.com/galaxyproject/tools-iuc/tree/main/data_managers/data_manager_bowtie2_index_builder):

![Terminal tree view of the folder structure of the bowtie2 Data Manager](../../images/data-manager-file-layout.png "Layout of a typical Data Manager folder"){:max-width="70%"}

Lets look at these components one-by-one:

1. There's a `data_manager` subfolder with the actual Data Manager **Tool** defined through a familiar tool xml file.

   This is what the admin is interacting with when installing new data (bowtie2 indices in this case).

2. In the root folder of the Data Manager, there is a **Data Manager Configuration** file.

   This file is always named `data_manager_conf.xml`.
   It declares, which .loc files the Data Manager will write to,
   defines how the output of the Data Manager Tool is to be translated into .loc file records,
   and where exactly Galaxy should store the data downloaded by the Data Manager Tool.

   When talking about Data Managers, their .loc files are called Data Tables.

3. A `.shed.yml` file that serves the same purpose (of declaring metadata for the Galaxy toolshed) as for regular tools.

   We will not discuss the contents of this file any further here.

4. A `test-data` folder with, you may guess it, data for testing the Data Manager.

   We will discuss testing Data Managers at the end of the tutorial.

5. A `tool-data` folder

   All tools, including regular ones, that use Data Tables must include a `tool-data` folder.
   In it, there needs to be one, typically empty or comments-only .loc.sample file for every used Data Table.

   Data Manager Tools are special in that they write to Data Tables instead of just reading from them, but the rule of
   "one .loc.sample file for each Data Table used" is independent of whether the Data Tables are read from or written to.
   In fact, Data Manager Tools may both write to Data Tables and use other Data Tables as a source for populating
   select boxes in their tool interface.
   The bowtie2 Data Manager, for example, lets the admin select the reference genome to build an index for
   from the list of installed genomes read from the `all_fasta` Data Table.
   The index it builds is useful for both bowtie2 and tophat2 so the path to the installed index will get recorded in the corresponding two Data Tables.
   The `tool-data` folder of this Data Manager, therefor, has three .loc.sample files -
   `all_fasta.loc.sample`, `bowtie2_indices.loc.sample` and `tophat2_indices.loc.sample`.
   Each of these files consists only of comment lines describing the .loc file purpose and expected structure.

   Comment lines in .loc.sample files are optional and just that: comments that are ignored by Galaxy.
   The only place that *Galaxy* reads .loc file structure from is the Data Table configuration file discussed next.

   > <details-title>Content of .loc.sample files</details-title>
   >
   > Tools can use .loc.sample files with actual record lines to point to data directly *shipping* with the tool.
   > That data would then also be stored in the `tool-data` folder (which explains its name).
   > Most tools and Data Managers, however, do not ship data directly and provide empty .loc.sample files with optional
   > comment lines like the bowtie2 Data Manager example.
   >
   {: .details}

6. A Data Table configuration file

   This file is always named `tool_data_table_conf.xml.sample` and provides the *layout* information for all Data Tables the Data Manager operates on and uses.

   For the bowtie2 Data Manager, this file has the following content:

   ```
   <tables>
       <!-- Locations of all fasta files under genome directory -->
       <table name="all_fasta" comment_char="#">
           <columns>value, dbkey, name, path</columns>
           <file path="tool-data/all_fasta.loc" />
       </table>
       <!-- Locations of indexes in the Bowtie2 mapper format -->
       <table name="bowtie2_indexes" comment_char="#">
           <columns>value, dbkey, name, path</columns>
           <file path="tool-data/bowtie2_indices.loc" />
       </table>
       <!-- Locations of indexes in the Bowtie2 mapper format for TopHat2 to use -->
       <table name="tophat2_indexes" comment_char="#">
           <columns>value, dbkey, name, path</columns>
           <file path="tool-data/tophat2_indices.loc" />
       </table>
   </tables>
   ```

   The Data Tables mentioned here are the same as the ones the `tool-data` folder has .loc.sample files for,
   but here we find the metadata for these tables, i.e. their name and .loc file name (which *could*,
   but probably shouldn't be different) and the identifiers for the different columns in each table.

   > <warning-title>Released Data Table layout information can never be changed again!</warning-title>
   >
   > Every version of every regular tool and of every Data Manager Tool using a given Data Table needs to declare its expected layout in a `tool_data_table_conf.xml.sample`.
   > An instance of Galaxy that gets conflicting layout information about the same Data Table from different tools or tool versions will refuse to load the Data Table
   > (or even refuse to start at all) until an admin resolves the conflict.
   >
   > This means that once you're publicly releasing a Data Manager Tool or regular tool with layout info for a new Data Table,
   > e.g. by uploading the tool to a Galaxy Toolshed, you must *never* change its layout again in a later version of that tool, or in a new tool reusing the Data Table!
   > No renaming of columns, no reshuffling, and no addition of new columns! You need to keep that layout frozen!
   >
   > It is, therefor, important to consider a new Data Table's columns very carefully.
   > If you are unsure whether the data it references will be versioned at some point in the future, better assume it will be and add that version column now!
   >
   > The only way to amend an inadequate table layout later is to have new tool versions declare a *new* Data Table. Try to avoid that confusing scenario if you can!
   >
   {: .warning}

7. A Data Table configuration test file

   This file is always named `tool_data_table_conf.xml.test`, is very similar to the `tool_data_table_conf.xml.sample`,
   but exists only for testing purposes, which, again, will be discussed at the end of the tutorial.

# How are Data Managers different from regular tools?

1. (Normally) only admins can run them through the Admin interface of Galaxy

2. They do create an output dataset in your history,
   but their **"side-effects"** are what really matters. As such side-effects
   Data Managers will:

   - typically, download or compute some data.
   - always instruct Galaxy to write at least one line of tab-separated info to
     at least one Data Table file.
   - if data was downloaded or computed, instruct Galaxy to

     - move the data to a permanent storage location
     - record the path to that location in the newly created Data Table record

# How does a Data Manager communicate with Galaxy?

1. It declares itself a Data Manager through the `tool_type="manage_data"` attribute.

   `<tool id="example_dm" name="An example Data Manager" version="1.0" tool_type="manage_data" profile="23.0">`

   This has several consequences:

   1. The tool will only appear in the Admin user interface.
   2. Galaxy will expect the output of the tool to be of `format="data_manager_json"`
      and its content to describe which columns should be added to new lines in which Data Tables.
   3. Galaxy will expect any data downloaded or computed by the tool to live in that output's
      `extra_files_path`.

      For example, if the output section of the Data Manager looks like this:

      ```
      <outputs>
          <data name="out_file" format="data_manager_json" label="${tool.name}"/>
      </outputs>
      ```
      it should make sure that its command section deposits data to be stored
      by Galaxy in `'$out_file.extra_files_path'`.
   4. The `data_manager_json` output file that the wrapper declares will exist
      **before** the command section runs and will contain a mapping of the input
      parameters,
      [among other things](https://docs.galaxyproject.org/en/latest/dev/data_managers.html#example-json-input-to-tool).
      You do not have to read the file if you don't want to, but you will have to overwrite it!

   > <comment-title>Minimal profile version for Data Managers</comment-title>
   >
   > Data Managers were executed in Galaxy's main environment until release 18.09!
   >
   > This means:
   >
   > - if you want to use `requirements` in a Data Manager Tool, you should set
   >   `profile="18.09"` or higher
   > - if you are bumping the profile version of an existing Data Manager to
   >   beyond 18.09, you may have to add requirements to it that bring in things
   >   the old version happened to find in Galaxy's environment.
   >
   {: .comment}

2. In its command section (or in a helper script called from there), the Data Manager Tool

   1. **overwrites** the already existing output file with a json of the items that should be
      added to one or more Data Tables.

   2. **creates** the folder at `output.extra_files_path` and deposits any data
      there that Galaxy should move then to a permanent storage location

3. The Data Manager Tool ships with a `data_manager_conf.xml` file,
   which forms the bridge between the `data_manager_json` file that it produces
   as output and the Data Table files Galaxy is supposed to add lines to.

   An example config file:

    ```
    <?xml version="1.0"?>
    <data_managers>
      <data_manager tool_file="data_manager/data_manager_cat.xml" id="data_manager_cat" >
        <data_table name="cat_database">  <!-- Defines a Data Table to be modified. -->
          <output> <!-- Handle the output of the Data Manager Tool -->
            <column name="value" /> <!-- columns that are going to be specified by the Data Manager Tool -->
            <column name="name" />  <!-- columns that are going to be specified by the Data Manager Tool -->
            <column name="database_folder" output_ref="out_file" >
              <move type="directory" relativize_symlinks="True">
                <source>${database_folder}</source>
                <target base="${GALAXY_DATA_MANAGER_DATA_PATH}">CAT/${database_folder}</target>
              </move>
              <value_translation>${GALAXY_DATA_MANAGER_DATA_PATH}/CAT/${database_folder}</value_translation>
              <value_translation type="function">abspath</value_translation>
            </column>
            <column name="taxonomy_folder" output_ref="out_file" >
              <move type="directory" relativize_symlinks="True">
                <source>${taxonomy_folder}</source>
                <target base="${GALAXY_DATA_MANAGER_DATA_PATH}">CAT/${taxonomy_folder}</target>
              </move>
              <value_translation>${GALAXY_DATA_MANAGER_DATA_PATH}/CAT/${taxonomy_folder}</value_translation>
              <value_translation type="function">abspath</value_translation>
            </column>
          </output>
        </data_table>
      </data_manager>
    </data_managers>
    ```

   This file declares column names for a single Data Table (`cat_database`) that
   Galaxy should add lines to based on the `data_manager_json` file returned by
   the `data_manager_cat` tool, and which might look like this:

   ```
   {'data_tables': {
       'cat_database': [
           {
               'database_folder': '<extra_files_path>/a_CAT_database',
               'name': '<extra_files_path>',
               'taxonomy_folder': '<extra_files_path>/a_taxonomy',
               'value': '<extra_files_path>'
           }
       ]
   }}
   ```

   Here, each innermost dictionary corresponds to one line that Galaxy should
   add to the Data Table `cat_database` and the keys in it match the column names
   declared in the `data_manager_conf.xml` file so Galaxy knows which dict value
   it should write into which column of the Data Table.

   The example tool downloads data, then extracts it into two folders, `a_CAT_database`
   and `a_taxonomy` under its output's `extra_files_path` folder. The tool wants
   Galaxy to record the `extra_files_path` folder name both in the value column and in the
   name column of the `cat_database` Data Table.

   It also wants to store the paths to the extracted `database_folder` and `taxonomy_folder` so that tools that later want to use that data can discover
   it from the corresponding columns of the Data Table.

   However, here's the issue:
   The Data Manager Tool at run time knows only the `extra_files_path`, but not the
   ultimate location that Galaxy will move the data to. This is where the more
   complicated parts of the above `data_manager_conf.xml` file enter the scene:

   ```
   <column name="database_folder" output_ref="out_file" >
       <move type="directory" relativize_symlinks="True">
           <source>${database_folder}</source>
           <target base="${GALAXY_DATA_MANAGER_DATA_PATH}">CAT/${database_folder}</target>
       </move>
       <value_translation>${GALAXY_DATA_MANAGER_DATA_PATH}/CAT/${database_folder}</value_translation>
       <value_translation type="function">abspath</value_translation>
   </column>
   <column name="taxonomy_folder" output_ref="out_file" >
       <move type="directory" relativize_symlinks="True">
           <source>${taxonomy_folder}</source>
           <target base="${GALAXY_DATA_MANAGER_DATA_PATH}">CAT/${taxonomy_folder}</target>
       </move>
       <value_translation>${GALAXY_DATA_MANAGER_DATA_PATH}/CAT/${taxonomy_folder}</value_translation>
       <value_translation type="function">abspath</value_translation>
   </column>
   ```

   The definitions of the `database_folder` column hold two types of instructions for Galaxy:

   1. The `<move>` element says that Galaxy should take (see the `<source>` element) the data that lives where the `${database_folder}` item of the `data_manager_json` output says it lives and move it to a destination `CAT/${database_folder}` under the base path indicated by `${GALAXY_DATA_MANAGER_DATA_PATH}` (which itself is the configured cached data storage path of the Galaxy instance).

   2. The first `<value_translation>` element says that Galaxy should not write the Data Manager Tool-provided value for `database_folder` directly, but instead first translate it to `${GALAXY_DATA_MANAGER_DATA_PATH}/CAT/${database_folder`. If you compare the resulting string with the `<move>` instructions, you will see that it will now be the same as the ultimate path to the folder after Galaxy has moved it.

      The second `<value_translation>` element simply says that Galaxy should turn the result of the first translation into an absolute path on the system. The result is then the value that will get written into the `database_folder` column of the `cat_database` table.

   The same logic is then used again to move the `taxonomy` folder to its final destination and to obtain the value to write to the corresponding Data Table column.

# How to test Data Manager Tools?

Unfortunately, the only automated tests you can run on a Data Manager Tool are the ones available for regular tools, too.

This means that you can test assumptions about the tool's `data_manager_json` output, about the command line formed and the stdout and stderr generated, but you can **not** verify that the Data Manager framework detects any data in the `extra_files_path` and moves it to the intended location.

For this reason, `planemo serve` is a very important command to use during any work on Data Managers!

> <comment-title>Getting planemo to work on a Data Managers</comment-title>
>
> `planemo test` and `planemo serve` work just fine for Data Managers,
> if you keep in mind that a Data Manager is more than just the tool xml file.
> If you're following the standard layout of Data Managers with the tool xml
> file in a subfolder, you need to run planemo from *outside* that subfolder to have
> it discover all required files beyond the tool xml, but point it to the tool
> xml to test or serve.
>
> - `planemo serve data_manager/bowtie2_index_builder.xml` run from the parent directory of the `data_manager` subfolder and
> - `planemo serve data_managers/data_manager_bowtie2_index_builder/data_manager/bowtie2_index_builder.xml` run from the root folder of the tools-iuc repo
>
> would both work fine.
>
> With `planemo serve` specifically, don't be surprised if the tool doesn't show
> up in the tools panel - it's not supposed to, but it's accessible from the
> Admin interface under *"Local Data"*.
>
{: .comment}

With a correctly written `tool_data_table_conf.xml.test` the Data Manager, during testing, will read from and write to the .loc files in its `test-data` folder.
This data is persistent across planemo runs as is the actual installed data, so after testing with `planemo serve`, you can inspect the .loc file records that have been written and check the path recorded there to see if the data has been installed the way you intended.

Before committing the test-data folder you may want to consider clearing the Data Tables you may have populated in it.

## Testing Data Managers and their client tools in combination

As said above, `planemo serve` and `planemo test` need to be run from the root folder of the Data Manager, but it's possible to point planemo to multiple xml files to test and you can use this to test both the Data Manager Tool and a tool using its Data Table in one session through, e.g.:

- `planemo serve data_manager/bowtie2_index_builder.xml ../../tools/bowtie2/bowtie2_wrapper.xml` run from the parent directory of the `data_manager` subfolder, or
- `planemo serve data_managers/data_manager_bowtie2_index_builder/data_manager/bowtie2_index_builder.xml tools/bowtie2/bowtie2_wrapper.xml` run from the root folder of the tools-iuc repo.

If you have previously served your Data Manager in isolation and installed some data, then, because this brings the Data Manager's test-data folder back into scope, that data will be immediately usable by the client tool.

# Data Manager Checklist

To sum up the key points of the above discussion, here is a practical checklist for Data Manager Tools.
In a few places this list goes beyond Galaxy's requirements for Data Manager Tools, but recommends standard approaches.

> <hands-on-title>Checklist for Data Manager Tools</hands-on-title>
>
> 1. Quick check of required files and key content
>
>    - Data Manager Tool XML
>
>      - {% icon point-right %} Tool XML (and helper scripts, if any) in `data_manager` subfolder
>      - {% icon galaxy-pencil %} Tool XML declares `tool_type="manage_data"` and a recent `profile`
>
>    - Declaration of Data Tables
>
>      - {% icon point-right %} `tool_data_table_conf.xml.sample` present in root folder
>      - {% icon galaxy-pencil %} file declares the layout of every Data Table touched by the tool and
>      - {% icon galaxy-pencil %} lists the columns of each table, including a *version* column if versioning the managed data might ever make sense and
>      - {% icon galaxy-pencil %} references a .loc file for every declared Data Table via a `<file path="tool-data/[Table name].loc" />` line
>      - {% icon point-right %} `tool-data` subfolder exists and has a `[Table name].loc.sample` file for every declared Data Table
>
>    - Data Manager framework integration
>
>      - {% icon point-right %} `data_manger_conf.xml` file present in root folder
>
>    - Toolshed readiness
>
>      - {% icon point-right %} `.shed.yml` file present in root folder
>
>    - Data Manager tests
>
>      - {% icon point-right %} `test-data` folder exists
>      - {% icon point-right %} `tool_data_table_conf.xml.test` present in root folder
>
> 2. Logic checks
>
>    - {% icon point-right %} Data Manager (or helper script) writes JSON to primary output
>    - {% icon point-right %} Data to be managed gets written to `extra_files_path`
>    - {% icon point-right %} `<move>` logic in `data_manger_conf.xml` handles transfer of all relevant content under `extra_files_path` to `<target base="${GALAXY_DATA_MANAGER_DATA_PATH}">` destination
>    - {% icon point-right %} `<value_translation>` logic in `data_manger_conf.xml` handles rewrite of all paths to managed data to point to target destinations
>    - {% icon point-right %} Data Manager (or helper script) deletes any irrelevant data left behind under `extra_files_path`
>
{: .hands_on}
