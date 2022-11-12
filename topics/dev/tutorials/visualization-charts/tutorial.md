---
layout: tutorial_hands_on

title: "JavaScript plugins"
questions:
  - "How can I make a custom visualization plugin for Galaxy?"
objectives:
  - "Learn how to add JavaScript plugins to Galaxy using the Charts visualization framework"
time_estimation: "1h"
key_points:
  - "Demonstrating a pluggable extension system for JavaScript visualizations"
  - "With three primary files we can integrate any JavaScript visualization into Galaxy"
subtopic: viz
contributors:
  - shiltemann
  - yhoogstrate
  - bgruening
  - guerler
  - dannon
---

## Introduction


In this tutorial we are going to demonstrate how to add a third party
JavaScript-based visualization to Galaxy, and we'll talk about what the benefits
are. The plugin we've selected for this exercise is the [*PV-Javascript Protein
Viewer*](https://biasmv.github.io/pv/). It is an open source protein structure
viewer for `PDB`-files. There are many other popular protein structure viewers
available for the visualization of `PDB`-files such as e.g.
[NGL](https://arose.github.io/ngl/) (also available in Galaxy) and
[JSMol](https://chemapps.stolaf.edu/jmol/jsmol/jsmol.htm).

> <details-title>Background: What is the PDB (Protein Data Bank) file format?</details-title>
>
> The `PDB`-file format contains atomic coordinates of biomolecules derived through a range of experimental and computational methods. Most commonly the file contains a spatial crystallographic snapshot of a protein. There are 100s of thousands of protein structures publicly available at the Protein Data Bank ([https://www.rcsb.org](https://www.rcsb.org)). Proteins are usually labeled by a four-letter code.
> Here is an example of a `PDB`-file for a hydrolase bond to its inhibitor (PDB: [1ACB](https://www.rcsb.org/pdb/explore/explore.do?structureId=1acb)):
>
> ```
> HEADER    HYDROLASE/HYDROLASE INHIBITOR           08-NOV-91   1ACB
> TITLE     CRYSTAL AND MOLECULAR STRUCTURE OF THE BOVINE ALPHA-CHYMOTRYPSIN-EGLIN
> TITLE    2 C COMPLEX AT 2.0 ANGSTROMS RESOLUTION
> ...
> KEYWDS    SERINE PROTEASE, HYDROLASE-HYDROLASE INHIBITOR COMPLEX
> AUTHOR    R.Z.KRAMER,L.VITAGLIANO,J.BELLA,R.BERISIO,L.MAZZARELLA,
> AUTHOR   2 B.BRODSKY,A.ZAGARI,H.M.BERMAN
> ...
> REMARK   1 REFERENCE 1
> REMARK   1  AUTH   M.BOLOGNESI,L.PUGLIESE,G.GATTI,F.FRIGERIO,A.CODA,L.ANTOLINI,
> REMARK   1  AUTH 2 H.P.SCHNEBLI,E.MENEGATTI,G.AMICONI,P.ASCENZI
> REMARK   1  TITL   X-RAY CRYSTAL STRUCTURE OF THE BOVINE
> REMARK   1  TITL 2 ALPHA-CHYMOTRYPSIN(SLASH)EGLIN C COMPLEX AT 2.6 ANGSTROMS
> ...
> SEQRES   1 E  245  CYS GLY VAL PRO ALA ILE GLN PRO VAL LEU SER GLY LEU
> SEQRES   2 E  245  SER ARG ILE VAL ASN GLY GLU GLU ALA VAL PRO GLY SER
> SEQRES   3 E  245  TRP PRO TRP GLN VAL SER LEU GLN ASP LYS THR GLY PHE
> ...
> ATOM      1  N   CYS E   1       2.323 -16.405  18.812  1.00 43.48           N
> ATOM      2  CA  CYS E   1       3.017 -15.136  18.786  1.00 35.11           C
> ATOM      3  C   CYS E   1       4.134 -15.068  19.799  1.00 32.90           C
> ATOM      4  O   CYS E   1       4.173 -15.810  20.772  1.00 41.38           O
> ATOM      5  CB  CYS E   1       2.052 -13.969  19.139  1.00 31.14           C
> ATOM      6  SG  CYS E   1       1.246 -14.085  20.788  1.00 34.72           S
> ATOM      7  N   GLY E   2       4.993 -14.081  19.607  1.00 21.94           N
> ...
> HETATM 2292  O   HOH E 406      12.343   1.842  12.901  0.86 18.70           O
> HETATM 2293  O   HOH E 407      -4.767  17.237  10.630  1.00 59.78           O
> HETATM 2294  O   HOH E 408      11.489  -6.278  18.740  0.96 20.00           O
> ...
> ```
>
> More resources on this file format:
>
>   - [https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format) ](https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format))
>   - [https://www.wwpdb.org/documentation/file-format ](http://www.wwpdb.org/documentation/file-format)
{: .details}

As mentioned above we will be focusing on the *PV-Javascript Protein Viewer* in
this tutorial. Now that we have learned about the underlying file format, let us
continue by visiting the protein viewer developer site at
[https://biasmv.github.io/pv/](https://biasmv.github.io/pv/) to get familiar
with the plugin.

> <hands-on-title>Hands-on</hands-on-title>
>
> 1. View the plugin in action, rotate the molecule and change its style.
>
> 2. Under which license is this plugin distributed?
>
> 3. Can you find the minified code file of this plugin?
>
{: .hands_on}

> <agenda-title></agenda-title>
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}


## Section 1 - Basic plugin setup

### Directory and plugin preparations

In this section we will download the viewer and add it to a local *Galaxy*
instance. All development takes place within the *Galaxy* codebase. The first
thing we are going to do is to clone a *Galaxy* instance and prepare the
directory structure for the new visualization plugin.

> <hands-on-title>Hands-on</hands-on-title>
>
> 1. Clone an instance of *Galaxy* in a path, further referred to as `$GALAXY_ROOT`:
>    ```bash
>    $ git clone https://github.com/galaxyproject/galaxy
>    ```
>
> 2. Navigate to the visualization plugin directory:
>    ```bash
>    $ cd $GALAXY_ROOT/config/plugins/visualizations
>    ```
>
> 3. Copy the existing example into a new directory:
>    ```bash
>    $ cp -r example myviz
>    ```
>
{: .hands_on}

This example visualization provides a great place to start with most of the
basics already in place.  Now that the directory structure is in place, let us
review the example visualization. Each visualization contains a set of <b>3</b>
files:

- Logo (`static/logo.png`) which will appear in *Chart*'s plugin selection interface.
- Configuration (`config/example.xml`) describing input parameters and options.
- Wrapper (`src/script.js`) which serves as a bridge between *Galaxy* and our 3rd-party plugin.

In the following sections we are going to discuss these files in more detail and
modify them to incorporate a new visualization. Let's start with the logo for
our visualization.

### Your visualization needs a logo

Each visualization is represented by a logo in the Galaxy interface. This makes
it easier for users to find and configure their visualization. The logo should
be in the `png`-file format. It will appear with a width of 120 pixels.

Here's an example [logo](../../files/charts-plugins/pdb/logo.png):

![Logo](../../files/charts-plugins/pdb/logo.png)

> <hands-on-title>Hands-on</hands-on-title>
>
> 1. Find an arbitrary image in `PNG`-file format. Possibly using *Google*'s [image search](https://images.google.com).
>
> 2. Copy it to the `myviz/static` directory and name it `logo.png`.
{: .hands_on}

### Configure the visualization

Each visualization has a configuration file. In this case it is named
`example.xml`. This file has conceptual similarities with a Tool's XML-file. It
allows developers to specify a variety of attributes and input parameters for
their visualization. Throughout this tutorial we are going to gradually augment
this file but for now we keep it simple.

> <hands-on-title>Hands-on</hands-on-title>
>
> 1. Rename the file to `config/myviz.xml`
>
> 2. Edit the file named `config/myviz.xml` and change the name and description.
>
> 3. Go to the Galaxy root directory to install dependencies, activate the virtual environment, and invoke the Galaxy client build.
>    ```bash
>    $ cd $GALAXY_ROOT
>    $ make setup-venv
>    $ source .venv/bin/activate
>    $ make client
>    ```
>
> 4. Run Galaxy
>    ```bash
>    $ ./run.sh
>    ```
>
{: .hands_on}

Your visualization should now be loaded.  You can verify that now by clicking
on `Visualize > Create Visualization` in the top menu bar of Galaxy and finding
your plugin with its new logo in the list there.

### Assign a new datatype to your visualization

> <hands-on-title>Hands-on</hands-on-title>
>
> 1. Open the file named `config/myviz.xml` and find the `<data_source>` section.
>
> 2. Replace the existing two `<test>` sections with:
>
>    `<test type="isinstance" test_attr="datatype" result_type="datatype">molecules.PDB</test>`
>
> 3. Remove the `<settings>` and `<groups>` sections.
>
{: .hands_on}

This links the plugins to the `PDB`-file format, which means that for any
history item of this file type the plugin will automatically be available.

### Modifying the wrapper

Now let's take a look at the wrapper which connects our visualization with
Galaxy. The wrapper consists of a module written in *JavaScript* and is
available at `src/script.js`:

The wrapper receives an `options` dictionary with the following <b>four</b> components:
 - *charts*: The model of the current visualization with attributes, settings etc.
 - *process*: A [jQuery.Deferred()](https://api.jquery.com/jquery.deferred/) object to allow asynchronous data requests within the wrapper
 - *dataset*: Details on the selected datasets such as url, ids etc. which can be used to access the dataset
 - *targets*: The DOM ids of the container elements to draw into

In this tutorial we will implement the *PV-Viewer* plugin. In order to execute a
3rd-party plugin we need to figure out how it works. This can be done by finding
a working example or documentation. Fortunately the *PV-Viewer* comes with both.
Let's take a look at the [documentation](https://pv.readthedocs.io/).

> <hands-on-title>Exploring the PV-Viewer</hands-on-title>
>
> 1. Identify the parameter which is needed to initialize the plugin when calling [*pv.Viewer()*](https://pv.readthedocs.io/en/v1.8.1/viewer.html#pv.Viewer).
>
> 2. Which of the wrapper option components represents this parameter?
>
> 3. Can you identify which `mode` settings are valid to render the structure with [*pv.Viewer.renderAs()*](https://pv.readthedocs.io/en/v1.8.1/viewer.html#pv.Viewer.renderAs)?
>
{: .hands_on}

Now that we have learned the basics on how the viewer plugin works, we can edit it and adjust  `script.js`.

> <hands-on-title>Hands-on</hands-on-title>
>
> 1. Access your visualization's `myviz/src` directory.
>    ```bash
>    $ cd $GALAXY_ROOT/config/plugins/visualizations/myviz
>    ```
>
> 2. Install the package for the *PV-Viewer*:
>    ```bash
>    yarn add bio-pv
>    ```
>
> 3. Modify your plugin's entrypoint script `src/script.js` to import the new `bio-pv` library we just added, add this at the top of the file:
>    ```js
>    import * as pv from 'bio-pv';
>    ```
>
> 4. Lastly, replace the `load` function contents in the same entrypoint file with the following code that will set up the protein viewer:
>     ```js
>     var viewer = pv.Viewer( document.getElementById( options.targets[ 0 ] ), {
>         width       : 'auto',
>         height      : 'auto',
>         antialias   : true,
>         outline     : true
>     });
>     $.ajax( {
>         url     : options.dataset.download_url,
>         success : function( response ) {
>             var structure = pv.io.pdb( response );
>             viewer.clear();
>             viewer.renderAs( 'protein', structure, 'cartoon', {} );
>             viewer.centerOn( structure );
>             viewer.autoZoom();
>             options.chart.state('ok', 'Chart drawn.');
>             options.process.resolve();
>         }
>     });
>     ```
{: .hands_on}

### Build the package

Now that we have completed the basic plugin definition, it is time to build the
scripts and libraries into a single bundle that Galaxy can use.  Galaxy
visualizations typically use [*Parcel*](https://parceljs.org) for this due to
its simplicity, and this is how the example project is already configured.  If
you're interested in more details, take a look at the `package.json` file in
your `myviz` plugin directory.

After it has been built and staged the plugin will be accessible through
*Galaxy*'s user interface. This process not require restarting your *Galaxy*
instance, just make sure to properly refresh your browser.

> <hands-on-title>Hands-on</hands-on-title>
>
> 1. Navigate to your visualization's root directory:
>    ```bash
>    $ cd $GALAXY_ROOT/config/plugins/visualizations/myviz
>    ```
>
> 2. Install any necessary JavaScript dependencies:
>    ```bash
>    $ yarn install
>    ```
>
> 3. Run the plugin build process using yarn:
>    ```bash
>    $ yarn run build
>    ```
>
> 4. Stage the scripts and run Galaxy:
>    ```bash
>    $ cd $GALAXY_ROOT
>    $ make client
>    $ ./run.sh
>    ```
>
{: .hands_on}

Lets test this.

### Test the visualization

In this section we will select a `PDB`-file from the Protein Data Bank and visualize it with our new plugin.

> <hands-on-title>Hands-on</hands-on-title>
>
> 1. Visit [https://www.rcsb.org ](http://www.rcsb.org) and select a protein structure e.g. [1ACB](http://www.rcsb.org/pdb/explore/explore.do?structureId=1acb)
>
> 2. Copy the link to the raw `PDB`-file e.g.
>
>    ```bash
>    https://files.rcsb.org/view/1ACB.pdb
>    ```
>
> 3. Start your Galaxy instance
>
>    ```bash
>    $ cd $GALAXY_ROOT
>    $ ./run.sh
>    ```
>
> 4. Open the upload dialog, paste the above link and click on *Start*.
>
> 5. Close the upload dialog, and select the file from the history panel on the right.
>
> 6. Click on the *diagram* icon. (You must be logged in)
>
> 7. Find your visualization and click on its logo.
>
{: .hands_on}

![First view](../../images/pv_1.png)

## Section 2 - Allow different rendering modes

In this section we are going to augment the visualization such that users may
select different rendering modes. Similar to a Tool's XML file, developers may
specify input parameters which then will be presented to the user. The
definition of Tool and Visualization input parameters are similar, however the
latter is provided in *JavaScript* and not as XML.

More information on parameters can be found in the [wiki](https://docs.galaxyproject.org/en/latest/dev/schema.html).

> <hands-on-title>Hands-on</hands-on-title>
>
> 1. Add the following block into the `myviz.xml` file:
>    ```xml
>    <settings>
>       <input>
>           <name>mode</name>
>           <label>Display mode</label>
>           <type>select</type>
>           <display>radio</display>
>           <value>cartoon</value>
>           <help>Select the rendering mode.</help>
>           <data>
>               <data>
>                   <value>cartoon</value>
>                   <label>Cartoon</label>
>               </data>
>               <data>
>                   <value>lines</value>
>                   <label>Lines</label>
>               </data>
>               <data>
>                    <value>points</value>
>                    <label>Points</label>
>               </data>
>           </data>
>       </input>
>    </settings>
>    ```
>
> 2. Change the following line in `script.js`:
>    ```js
>    viewer.renderAs( 'protein', structure, 'cartoon', {} );
>    ```
>    to
>
>    ```js
>    var settings = options.chart.settings;
>    viewer.renderAs( 'protein', structure, settings.get( 'mode' ), {} );
>    ```
>
> 3. Rebuild the plugin
>
>    ```bash
>    $ yarn run build
>    $ cd $GALAXY_ROOT
>    $ make client
>    ```
>
> 4. Refresh your browser.
>
> 5. Load your visualization and test different rendering modes in the *Customization* tab of your visualization.
>
{: .hands_on}

## Conclusion


First of all, thank you for completing this tutorial. We have learned how to add
JavaScript visualizations to Galaxy utilizing the Charts framework.
