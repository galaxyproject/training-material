Visualizations in Galaxy
============================================

:grey_question: ***Questions***

- *How can visualization plugins benefit science?*

:dart: ***Objectives***

- *Implement a first Galaxy visualization*
- *Understand the client side vs. server side principle*

:heavy_check_mark: ***Requirements***

- *Galaxy introduction*
- *Javascript knowledge*

:hourglass: ***Time estimation*** *1.5h*

# Introduction

Visualizations may be very helpful in understanding data better. There is a whole
range of visualizations, from rather simple scatter and barplots up to projections
of high dimensional data or even entire genomes. Many of these visualizations often
require a lot of tweaking and changes in settings like zooming in and assigning colors, etc.
Therefore, visualizations are ideally interactive, and changing settings is often
an initial step in exploring data. For this reason it may be inconvenient to make use
of static galaxy tools because it lacks these interactive features. For these situations Galaxy
offers the option to create *visualizations plugins*, file format specific javascripts
that integrate with the history menu, without making redundant copies of data.

In this tutorial we shall go through how this system works and create a simple visualization
plugin. The tool will create a visualization of the number of aligned reads per
chromosome of a BAM file, and we will discuss possible optimizations and advantages
and disadvantages of the proposed implementation.

If you want to make visualizations ready for production, it is essential to have a good
understanding of HTML5 and JavaScript as these are the basic languages in which they are written.
However, for this tutorial we will keep it basic.

Additional documentation about Galaxy visualizations can be found here:

- https://wiki.galaxyproject.org/VisualizationsRegistry
- https://wiki.galaxyproject.org/VisualizationsRegistry/Cookbook
- https://wiki.galaxyproject.org/VisualizationsRegistry/Code
- https://wiki.galaxyproject.org/DataProviders
- https://wiki.galaxyproject.org/DataProviders/Cookbook
- https://wiki.galaxyproject.org/Develop/Visualizations

# Part 1

The visualization we would like to write is a tool that shows the number of aligned
reads per chromosome, of a BAM file. The first thing we need to do is to come up with a name.
Let's call it *alignment_rname_boxplot*. Note that the reference sequences (usually chromosomes)
to which we align are named `RNAME` in the BAM/SAM specification.

The development of a Galaxy visualization takes place within the Galaxy codebase.
After cloning an instance of Galaxy in a path, further referred to as `$GALAXY_ROOT`,
we find the plugin directory as follows:

```bash
$ cd $GALAXY_ROOT/config/plugins/visualizations
```

In here we need to make a new directory for our new plugin project:

```bash
$ mkdir alignment_rname_boxplot
$ cd alignment_rname_boxplot
```

We need to make three (sub-)directories to complete the structure of the project:

```bash
$ mkdir config
$ mkdir static
$ mkdir templates
```

## Linking the plugin with Galaxy

To create a bridge between the not yet written plugin and Galaxy, we need to write a configuration in XML format.
Create the following file:  `config/alignmnet_rname_boxplot.xml` and give it the following data:

```xml
  <?xml version="1.0" encoding="UTF-8"?>
  <!DOCTYPE visualization SYSTEM "../../visualization.dtd">
  <visualization name="alignment_rname_boxplot">
      <data_sources>
          <data_source>
              <model_class>HistoryDatasetAssociation</model_class>

              <test type="isinstance" test_attr="datatype" result_type="datatype">binary.Bam</test>
              <test type="isinstance" test_attr="datatype" result_type="datatype">tabular.Sam</test>

              <to_param param_attr="id">dataset_id</to_param>
          </data_source>
      </data_sources>
      <params>
          <param type="dataset" var_name_in_template="hda" required="true">dataset_id</param>
      </params>
      <template>alignment_rname_boxplot.mako</template>
  </visualization>
```

This configures the projects' name, which shall appear on pressing the visualization button in
the history menu. It also links the plugin to two file formats: BAM and SAM, which means that
for any history item of these file formats the plugin will automatically become available.
It also includes a configuration of the template file, which is a mako template
(HTML + Python syntax) found in the `./templates` directory. The `var_name_in_template` is set to `hda`, which will be the name of the variable in the mako template corresponding to the visualized dataset.

## Creating the visualization

We have linked Galaxy to a mako file (which we did not yet create).
This file is a blueprint for the visualization.
This means that for every invocation of the visualization, the mako file will be compiled to render a HTML file.
It is trivial to understand that compilation needs to take place in a few ms because otherwise loading becomes too slow, so computational intensive tasks can not be done prior to loading the visualization.
A bit of server side rendering at itself is not a problem, but the visualizations written in HTML and/or JS are supposed to to the actual calculations and conversions at the client side (in the browser).
Therefore, unlike static galaxy tools, parsing files does not take place at the server, but instead data will be downloaded by the client via an exposed Galaxy URL prior to client side rendering.

The most trivial part of the mako file are the variables used for further web development, given below:

```python
  <!DOCTYPE HTML>
  <%
      ## Generates hash (hdadict['id']) of history item
      hdadict = trans.security.encode_dict_ids( hda.to_dict() )

      ## Finds the parent directory of galaxy (/, /galaxy, etc.)
      root     = h.url_for( '/' )

      ## Determines the exposed URL of the ./static directory
      app_root = root + 'plugins/visualizations/'+visualization_name+'/static/'

      ## Actual file URL:
      file_url = root + 'datasets/' + hdadict['id'] + "/display?to_ext=" + hda.ext;
  %>
```

The `hdadict` is a variable that contains a file identifier that has been encoded to it's exposable uid.
Here `root` is the root of Galaxy on the webserver (e.g. /, /galaxy/, /galaxy-pub/) and `app_root` is the exposed url of the static files.
The `file_url` is the exposed url of the dataset selected for visualization.

We could obtain the BAM file client side by downloading the BAM/SAM file via *file_url* and a javascript.
However, BAM files can become really large and it is not desired to pump such data over the network.
It is also inconvenient to parse the BAM file via JS just to count the number of reads.

Fortunately, BAM files (should) have indices.
These indices are brief summaries describing the number of entries per chromosome, to access them more quickly.
In the python part of the mako template we can access an index as *hda.metadata.bam_index*.
Please remark that this is a file path on the server and not an exposed URL.
It is probably not that difficult to obtain the exposed Galaxy URL of the index file, but this file is binary and compressed and therefore requires more advanced programming to extract the information, which is not the scope of this practical.

Samtools idxstats is a program in samtools able to do this advanced binary decompression for you.
However, since visualizations do not have dependency management, it is very tricky to let the mako template do a system call to samtools.
Fortunately, the galaxy ecosystem ships with a built-in pysam dependency, a library that can do any native samtools command within python.

The `.metadata.bam_index` is a special kind of file in the Galaxy ecosystem.
It is actually an invisible file in Galaxy, linked to another history item, but does have a unique filename.
So, for the BAM file ./database/files/000/dataset_001.dat, our BAI file is **NOT** ./database/files/000/dataset_001.dat.bai but could be ./database/files/000/dataset_002.dat or ./database/files/000/dataset_003.dat.
To solve this we can create a symlink - if it does not already exists - to make sure they have the same prefix to make them compatible with samtools and pysam.
This trick is rather similar to how it is applied in static Galaxy tools.
We can make this symlink as follows:

```python
  ## Ensure BAI index is symlinked
  bai_target = hda.file_name+'.bai'
  import os

  if not os.path.isfile(bai_target):
      os.symlink(hda.metadata.bam_index.file_name, bai_target)
```

Now the BAM file has a BAI file with the same name-prefix, and we can extract the idxstats as follows:

```python
  import pysam
  data = pysam.idxstats(hda.file_name)
```

We could decide to remove the symlink after using it, but this may have deeper implications, so until we study the consequences I would recommend to keep it.
With the lines of python code above, the idxstats data is parsed into the RAM of python during compilation at the server, but is not yet exported into the HTML page nor parsed by JS.
To dump the idxstats into a HTML file, please create the mako file: ./templates/alignment_rname_boxplot.mako and fill it with the following code:

```python
  <!DOCTYPE HTML>
  <%
      ## Generates hash (hdadict['id']) of history item
      hdadict = trans.security.encode_dict_ids( hda.to_dict() )

      ## Finds the parent directory of galaxy (/, /galaxy, etc.)
      root     = h.url_for( '/' )

      ## Determines the exposed URL of the ./static directory
      app_root = root + 'plugins/visualizations/'+visualization_name+'/static/'

      ## Actual file URL:
      file_url = root + 'datasets/' + hdadict['id'] + "/display?to_ext=" + hda.ext;


      ## Ensure BAI index is symlinked
      bai_target = hda.file_name+'.bai'
      import os

      if not os.path.isfile(bai_target):
          os.symlink(hda.metadata.bam_index.file_name, bai_target)


      ## Extract idxstats
      import pysam
      bam_idxstats_data = pysam.idxstats(hda.file_name)
  %>
  <html>
      <head>
         <title>${hda.name | h} | ${visualization_name | h}</title>
      </head>
      <body>
          ${bam_idxstats_data | h}
      </body>
  </html>
```

In here you see `${bam_idxstats_data | h}` which prints the python variable into the HTML page and does also HTML escaping by giving the ` | h`-flag (for security reasons).
We can now test this basic visualization and we need a (small) BAM file for it.
You can find small examples in the test files of IUC tools.
Please download: https://github.com/galaxyproject/tools-iuc/raw/master/tools/hisat2/test-data/hisat_output_1.bam
Go the galaxy root directory and start Galaxy:

```bash
$ cd $GALAXY_ROOT
$ ./run.sh
```

Upload the BAM file to the history and take a look at the history item.
If everything went well the visualization has appeared in the history item and proceed with the visualization.
Once you open it, it will show the contents of idxstats, compiled to HTML:
```
['phiX174\t5386\t19\t1\n', '*\t0\t0\t0\n']
```

It contains two entries, one to `phiX` and one to `*`, of which the former is a chromosome and the latter the unmapped reads.
Entries are tab delimited (`\t`) and for the phiX entry it indicates that the length of the RNAME is 5386 bases and 19 reads are aligned to it.

The trick we need to apply is to make the data a bit more usable for Javascript by converting it into a simple dictionary of the following syntax:

```python
  {'phiX': 19, '*': 0}
```

Although it is possible to do this in python I recommend to do this in JS.
The var dump provided by python/pysam is actually a valid syntax for Javascript too.
So getting the raw data into javascript is rather easy:

```html
  <script>
      bam_idxstats_data = ${bam_idxstats_data};
  </script>
```

Converting the data is not the scope of the tutorial, so here we provide such a function:

```html
<script>
  bam_idxstats_data = ${bam_idxstats_data};
  function parse_data(bam_idxstats_data) {
    var output = {};

    for(var i = 0; i < data.length ; i++) {
      var line = data[i];
      var chunks = line.split("\t");

      if(chunks[0].split("_").length == 1) { // only if it does not contain underscore
        output[chunks[0]] = parseInt(chunks[2]);
      }
    }

   return output;
  }
  </script>
```

The great thing about the mako system is that it does not require to restart galaxy in order to make functional changes to the mako files.
So, please open our visualization in a new browser tab, so we can press [F5] to test changes.
Change the mako file to the following:

```python
    <!DOCTYPE HTML>
    <%
        ## Generates hash (hdadict['id']) of history item
        hdadict = trans.security.encode_dict_ids( hda.to_dict() )

        ## Finds the parent directory of galaxy (/, /galaxy, etc.)
        root     = h.url_for( '/' )

        ## Determines the exposed URL of the ./static directory
        app_root = root + 'plugins/visualizations/'+visualization_name+'/static/'

        ## Actual file URL:
        file_url = root + 'datasets/' + hdadict['id'] + "/display?to_ext=" + hda.ext;


        ## Ensure BAI index is symlinked
        bai_target = hda.file_name+'.bai'
        import os

        if not os.path.isfile(bai_target):
            os.symlink(hda.metadata.bam_index.file_name, bai_target)


        ## Extract idxstats
        import pysam
        bam_idxstats_data = pysam.idxstats(hda.file_name)


        ## Cleanup ? os.remove(bai_target)
    %>
    <html>
        <head>
           <title>${hda.name | h} | ${visualization_name | h}</title>
            <script>
               bam_idxstats_data = ${bam_idxstats_data};
                function parse_data(data) {

                    var output = {};

                    for(var i = 0; i < data.length ; i++) {
                        var line = data[i];
                        var chunks = line.split("\t");

                        if(chunks[0].split("_").length == 1) { // only if it does not contain underscore
                            output[chunks[0]] = parseInt(chunks[2]);
                        }
                    }

                    return output;
                }
            </script>
        </head>
        <body onload="bam_idxstats = parse_data(bam_idxstats_data);">
            ${bam_idxstats_data | h}
        </body>
    </html>
```

Refresh the browser tab containing the visualization with [F5] and open the developers console in the browser.
Please ask for help if you don't know how to do this, this setting is different for every browser.
Within the console, type: `idx_stats` and press [Enter].
This should give the parsed contents as a dictionary, which can directly be used in javascript.

From this point I would like you to continue on your own to see if you are able to create a simple visualization of this dictionary.
You may think of tables, DIVs or even more complicated solutions :).

If you want to check out a simple solution I made, please have a look at:
[../plugins/alignment_rname_boxplot/](alignment_rname_boxplot)

In the given example, the `RNAME` queries containing an underscore were removed.
This is because there are many alternative chromosomes, making the list very large for certain reference genomes.
However, for certain studies it might be desired to just look at those.
Also, the current plot is a box plot, but I can imagine that a pie-chars can be convenient too.
All of those additional settings can be implemented for interactive behaviour, contributing to quicker understanding of the data which is generally not so convenient using static Galaxy tools.

### Tip: static files

In the example we included Javascript and CSS into the HTML website.
Remember that for every new invocation of the visualization the entire CSS en JS are copied and transferred as well.
This is a waste of (redundant) bandwidth as we could save the files in the static directory and refer to them within the HTML.
The browser shall check it's cache for the presence of libs and style sheets and only update them if they have changed.

### Improvements
Another thing you may realize is that we still do the calculation (pysam.idxstats) server side.
Although this is a marginal calculation, it is causing delay and the bigger the files the worse.
The underlying problem is that we did not obtain and parse the BAI file via javascript.
This would be a more elaborate solution, but requires more development time as we need to develop a function able to parse the binary file.

Another thing we could do is create a static `samtools idxstats` tool, that creates a file of datatype `tabular.Idxstats` and include that datatype into Galaxy.
We then make the visualization specific for that datatype, just plotting the results of the `idxstats`.

A fundamental and more complicated problem is that BAM files are simply too big to transfer for these kind of applications.
It would be ideal to have web server integration that allows to query specific locations or metadata within or from a BAM file where indexing operations are taken care of at the server side.
This is what has been done in Trackster.

# Conclusion

We have just created a visualization plugin in Galaxy, that is able to visualize the number of alignments per `RNAME` in a BAM file.

:grey_exclamation: ***Key Points***

- *Visualizations require a different way of thinking*
    * server and client side
    * downloading files rather than system level access
- *Interactivity is what makes visualizations different from static tools*
- *Requires understanding of both the Galaxy ecosystem as well as HTML5/JS*
- *Performance is more important than for static Galaxy tools*

# :clap: Thank you
