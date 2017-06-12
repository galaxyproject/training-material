---
layout: tutorial_hands_on
topic_name: Dev-Corner
tutorial_name: visualizations
---

# Introduction


# Part 1

In this tutorial we are going to add a 3rd-party visualization to the *Charts* framework.

We will use the *PV-Javascript Protein Viewer* for this purpose. It is an open source, protein structure viewer for the PDB-file format.
See: https://biasmv.github.io/pv/.

The development of a Galaxy visualization takes place within the Galaxy codebase.

> ### :pencil2: Hands-on: Directory and plugin preparations
>
> 1. Clone an instance of Galaxy in a path, further referred to as `$GALAXY_ROOT`
> 2. Explore the plugin directory as follows:
>
>    ```bash
>    $ cd $GALAXY_ROOT/config/plugins/visualizations/charts
>    $ cd static/repository/visualizations/
>    ```
>
> 3. Add the visualization to `registry.json`:
>
>    ```bash
>    "myviz" : [ "pdb" ]
>    ```
>
> 4. Create a new directory for our new plugin project
>
>    ```bash
>    $ mkdir myviz
>    $ mkdir myviz/pdb
>    $ cd myviz/pdb
>    ```
>
> 5. Download the 3rd-party code:
>
>    ```bash
>    curl https://raw.githubusercontent.com/biasmv/pv/master/bio-pv.min.js -o plugin.js
>    ```
>

Each visualization requires three files. A logo (`logo.png`), a configuration file (`config.js`) and a wrapper (`wrapper.js`). We are about to create and discuss these files now.

## Your visualization needs a logo

Each visualization is represented by a logo in the *Charts* interface. This makes it easier for users to find and configure their visualization. The logo file should be in `png` format. Find and download a png-file, name it `logo.png` and copy it into the the `myviz/pdb` directory.

![Specific region](/Dev-Corner/images/exercise.png)

## Visualization annotation

Each visualization requires a configuration file named `config.js`. 

> ### :pencil2: Hands-on: Data upload
>
> Create the file  `config.js` with the following contents:
>
> ```js
define( [], function() {
    return {
        title       : 'A PDB viewer',
        library     : 'My Visualization',
        description : 'Displays Protein Structures.',
        datatypes   : [ 'pdb' ],
        keywords    : []
    }
});
> ```

TODO: Explain available options.
This configures the plugin's name and a description which will appear on the Charts selection interface. It also links the plugin to the PDB-file format, which means that for any history item of these file type the plugin will automatically become available. Keywords are optional and can help to improve the annotation.

http://pv.readthedocs.io/en/v1.8.1/viewer.html#pv.Viewer
http://pv.readthedocs.io/en/v1.8.1/viewer.html#pv.Viewer.renderAs

## Creating the visualization

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

```python
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

## Creating the visualization

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

```js
define( [ 'visualizations/myviz/pdb/plugin' ], function( pv ) {
    return Backbone.Model.extend({
        initialize: function( options ) {
            var settings = options.chart.settings;
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
                    viewer.renderAs( 'protein', structure, settings.get( 'mode' ), settings.attributes );
                    viewer.centerOn( structure );
                    viewer.autoZoom();
                    options.process.resolve();
                }
            });
        }
    });
});

```

```bash
$ rm static/repository/build/myviz_pdb.js
$ webpack
```

# Part 3

Lets build a form.


## Creating the visualization

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