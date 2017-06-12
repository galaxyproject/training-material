---
layout: tutorial_hands_on
topic_name: Dev-Corner
tutorial_name: visualizations
---

# Introduction

In this tutorial we are going to demonstrate how to add a 3rd-party visualization to *Charts* and what the benefits are. The plugin we select for this purpose is the *PV-Javascript Protein Viewer*. It is an open source, protein structure viewer for `PDB`-files.

> ### What is the PDB-file format
> 
> ```bash
HEADER    EXTRACELLULAR MATRIX                    22-JAN-98   1A3I
TITLE     X-RAY CRYSTALLOGRAPHIC DETERMINATION OF A COLLAGEN-LIKE
TITLE    2 PEPTIDE WITH THE REPEATING SEQUENCE (PRO-PRO-GLY)
...
EXPDTA    X-RAY DIFFRACTION
AUTHOR    R.Z.KRAMER,L.VITAGLIANO,J.BELLA,R.BERISIO,L.MAZZARELLA,
AUTHOR   2 B.BRODSKY,A.ZAGARI,H.M.BERMAN
...
REMARK 350 BIOMOLECULE: 1
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, B, C
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
...
SEQRES   1 A    9  PRO PRO GLY PRO PRO GLY PRO PRO GLY
SEQRES   1 B    6  PRO PRO GLY PRO PRO GLY
SEQRES   1 C    6  PRO PRO GLY PRO PRO GLY
...
ATOM      1  N   PRO A   1       8.316  21.206  21.530  1.00 17.44           N
ATOM      2  CA  PRO A   1       7.608  20.729  20.336  1.00 17.44           C
ATOM      3  C   PRO A   1       8.487  20.707  19.092  1.00 17.44           C
ATOM      4  O   PRO A   1       9.466  21.457  19.005  1.00 17.44           O
ATOM      5  CB  PRO A   1       6.460  21.723  20.211  1.00 22.26           C
...
HETATM  130  C   ACY   401       3.682  22.541  11.236  1.00 21.19           C
HETATM  131  O   ACY   401       2.807  23.097  10.553  1.00 21.19           O
HETATM  132  OXT ACY   401       4.306  23.101  12.291  1.00 21.19           O
...
> ```

[More resources on PDB] (https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format))

Let's start by visiting the viewers developer site at https://biasmv.github.io/pv/ to get familiar with this plugin.

# Part 1

In the following sections we will download this viewer and add it to a local Galaxy instance. This development takes place within the Galaxy codebase. The first thing we are going to do is to clone a Galaxy instance and prepare the directory structure for the new visualization.

> ### Directory and plugin preparations
>
> 1. Clone an instance of Galaxy in a path, further referred to as `$GALAXY_ROOT`:
>    ```bash
>    $ git clone https://github.com/galaxyproject/galaxy 
>    ```
>
> 2. Navigate to the *Charts* repository root:
>
>    ```bash
>    $ cd $GALAXY_ROOT/config/plugins/visualizations/charts/static/repository
>    ```
>
> 3. Register your visualization by adding a line to the `registry.json`:
>
>    ```bash
>    "myviz" : [ "pdb" ]
>    ```
>
> 4. Create a new directory:
>
>    ```bash
>    $ mkdir -p visualizations/myviz/pdb
>    $ cd visualizations/myviz/pdb
>    ```
>
> 5. Download the PV-Viewer from Github:
>
>    ```bash
>    curl https://raw.githubusercontent.com/biasmv/pv/master/bio-pv.min.js -o plugin.js
>    ```
>

Each visualization requires three files. A logo (`logo.png`), a configuration file (`config.js`) and a wrapper (`wrapper.js`). We are about to create and discuss these files now.

## Your visualization needs a logo

Each visualization is represented by a logo in the *Charts* interface. This makes it easier for users to find and configure their visualization. The logo file should be in `png` format. Find and download a png-file, name it `logo.png` and copy it into the the `myviz/pdb` directory.

![Specific region](/Dev-Corner/images/exercise.png)

## Annotate the visualization

Each visualization has a configuration file named `config.js`. This file has conceptual similarities with a Tool's XML-file. It allows developers to specify a variety attributes and input parameters for their visualization. In the following sections we are going to gradually augment this file. For now, we keep it simple. Create a `config.js` file with the following contents:

```js
define( [], function() {
    return {
        title       : 'A PDB viewer',
        library     : 'My Visualization',
        description : 'Displays Protein Structures.',
        datatypes   : [ 'pdb' ],
        keywords    : []
    }
});
```

This configures the plugin's name and a description which will appear on the Charts selection interface. It also links the plugin to the PDB-file format, which means that for any history item of these file type the plugin will automatically become available. Keywords are optional and can help to improve the annotation.

http://pv.readthedocs.io/en/v1.8.1/viewer.html#pv.Viewer
http://pv.readthedocs.io/en/v1.8.1/viewer.html#pv.Viewer.renderAs

## Add a wrapper to connect Galaxy with the PV-Viewer plugin

```js
define( [ 'visualizations/myviz/pdb/plugin' ], function( pv ) {
    return Backbone.Model.extend({
        initialize: function( options ) {
            var viewer = pv.Viewer( document.getElementById( options.targets[ 0 ] ), {
                width       : 'auto',
                height      : 'auto',
                antialias   : true,
                outline     : true
            });
            $.ajax( {
                url     : options.dataset.download_url,
                success : function( response ) {
                    var structure = pv.io.pdb( response );
                    viewer.clear();
                    viewer.renderAs( 'protein', structure, 'cartoon', {} );
                    viewer.centerOn( structure );
                    viewer.autoZoom();
                    options.process.resolve();
                }
            });
        }
    });
});

```


Lets pack the visualization:

```bash
$ cd $GALAXY_ROOT/config/plugins/visualizations/charts
$ npm install
$ webpack
```

Lets go and test it.

http://www.rcsb.org/pdb/explore/explore.do?structureId=1acb
https://files.rcsb.org/view/1ACB.pdb

> ### :pencil2: Hands-on: Data upload
>
> 1. Run your Galaxy instance
> 2. Open a new tab and go to the http://www.rcsb.org
> 3. Download a PDB-file or copy the link to it
> 4. Run you Galaxy instance
> 5. Upload the PDB file to your instance
> 6. Select the file in the history
> 7. Start Charts 

# Part 2

Lets build a form. TODO: Links https://docs.galaxyproject.org/en/latest/dev/schema.html. Explain how XML-parameters map to this json format.

## Add input parameters

```js
define( [], function() {
    return {
        title       : 'A PDB viewer',
        library     : 'My Visualization',
        datatypes   : [ 'pdb' ],
        keywords    : 'pdb protein structure',
        description : 'Galaxy tutorial.',
        settings    : {
            mode : {
                label   : 'Render as:',
                help    : 'Select the rendering mode.',
                type    : 'select',
                display : 'radio',
                value   : 'cartoon',
                data    : [ { label : 'Cartoon',        value : 'cartoon' },
                            { label : 'Lines',          value : 'lines' },
                            { label : 'Points',         value : 'points' },
                            { label : 'Spheres',        value : 'spheres' },
                            { label : 'Trace',          value : 'trace' },
                            { label : 'Trace (line)',   value : 'lineTrace' },
                            { label : 'Trace (smooth)', value : 'sline' },
                            { label : 'Tube',           value : 'tube' } ]
            }
        }
    }
});

```

## Update the wrapper

Change the following line in `wrapper.js`
```js
    viewer.renderAs( 'protein', structure, 'cartoon', {} );
```
to
```js
    var settings = options.chart.settings;
    viewer.renderAs( 'protein', structure, settings.get( 'mode' ), settings.attributes );
```

## Rebuild the plugin

```bash
$ rm static/repository/build/myviz_pdb.js
$ webpack
```

# Part 3

Add more options to configure.


## Add 'pointSize', 'lineWidth' and 'radius' to `config.js`

```js
define( [], function() {
    return {
        title       : 'A PDB viewer',
        library     : 'My Visualization',
        datatypes   : [ 'pdb' ],
        keywords    : 'pdb protein structure',
        description : 'Galaxy tutorial.',
        settings    : {
            mode : {
                label   : 'Render as:',
                help    : 'Select the rendering mode.',
                type    : 'select',
                display : 'radio',
                value   : 'cartoon',
                data    : [ { label : 'Cartoon',        value : 'cartoon' },
                            { label : 'Lines',          value : 'lines' },
                            { label : 'Points',         value : 'points' },
                            { label : 'Spheres',        value : 'spheres' },
                            { label : 'Trace',          value : 'trace' },
                            { label : 'Trace (line)',   value : 'lineTrace' },
                            { label : 'Trace (smooth)', value : 'sline' },
                            { label : 'Tube',           value : 'tube' } ]
            },
            pointSize: {
                label : 'Point size',
                help  : 'Specify the point size.',
                type  : 'float',
                min   : 0.1,
                max   : 10,
                value : 1
            },
            lineWidth : {
                label : 'Line width',
                help  : 'Specify the line width.',
                type  : 'float',
                min   : 0.1,
                max   : 10,
                value : 4
            },
            radius : {
                label : 'Radius',
                help  : 'Radius of tube profile. Also influences the profile thickness for helix and strand profiles.',
                type  : 'float',
                min   : 0.1,
                max   : 3,
                value : 0.3
            }
        }
    }
});

```

```bash
$ rm static/repository/build/myviz_pdb.js
$ webpack
```