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