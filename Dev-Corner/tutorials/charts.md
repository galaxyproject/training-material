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

> ### :pencil2: Hands-on: Data upload
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
>
>    $ cd myviz/pdb
>    ```
>
> 5. Download the 3rd-party code:
>
>    ```bash
>    curl https://raw.githubusercontent.com/biasmv/pv/master/bio-pv.min.js -o plugin.js
>    ```
>

## Find a logo

Find a png-file, name it `logo.png` and copy it into the the `myviz/pdb` directory.

## Linking the plugin with Galaxy

To create a bridge between our not-yet-written plugin and Galaxy, we need to write a
configuration in XML format.

> ### :pencil2: Hands-on: Data upload
>
> Create the file  `config.js` with the following contents:
>
> ```js
define( [], function() {
    return {
        title       : 'A PDB viewer',
        library     : 'My Visualization',
        datatypes   : [ 'pdb' ],
        keywords    : 'pdb protein structure',
        description : 'Galaxy tutorial.'
    }
});
> ```

This configures the plugin's name, which shall appear on pressing the visualization button in
the history menu. It also links the plugin to two file formats: BAM and SAM, which means that
for any history item of these file formats the plugin will automatically become available.

It also includes a reference to a mako template file (HTML + Python syntax), to be found in the
`templates` directory (we will create this file in the next section). The `var_name_in_template`
parameter is set to the value `hda`, which will be the name of the variable in the mako template
corresponding to the dataset to be visualized.

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