Visualizations in Galaxy
============================================

:grey_question: ***Questions***

- *Major question that would be addressed in this tutorial (mostly general biological questions)*
- *Second question*
- *Third question*
- *...*

:dart: ***Objectives***

- *First objective of this tutorial (It is a single sentence describing what a learner will be able to do once they have sat through the lesson. The objectives must be technical, but also theoretical, objectives. You can check [SWC lessons](http://swcarpentry.github.io/instructor-training/19-lessons/) to help you writing learning objectives.)*
- *Second objective*
- *Third objective*
- *...*

:heavy_check_mark: ***Requirements***

- *Galaxy introduction*
- *Javascript Knowledge*
- *...*

:hourglass: ***Time estimation*** *1d/3h/6h*

# Introduction

Visualizations may be very helpful in understanding data better.
These can be simple barplots but also projections of high dimensional data or even genomes.
What is characteristical about visualizations is that they often require endless tweaking and change in settings like zoom and colors.
Most often visualizations are ideally interactive, where changing the settings is effectively the browsing through the data.
For this reason it may be inconvenient to make use of a static galaxy tool because it misses interactive features.
Therefore galaxy offers the option to create *visualizations plugins*, file format specific javascripts that intergrate with the history menu.

In this tutorial we shall go through how this works and create a simple tool ourselves.
The tool will create a visualization of the number of aligned reads per chromosome in a BAM file, and we will discuss possible optimizations.

Because visualizations are written in HTML5 and JavaScript it is essential to have a good understand if you want to make visualizations ready for production.
However, for this tutorial we will keep it very basic.


Documentation about Galaxy visualizations can be found here. It can be used as base:

- https://wiki.galaxyproject.org/VisualizationsRegistry
- https://wiki.galaxyproject.org/VisualizationsRegistry/Cookbook
- https://wiki.galaxyproject.org/VisualizationsRegistry/Code
- https://wiki.galaxyproject.org/DataProviders
- https://wiki.galaxyproject.org/DataProviders/Cookbook
- https://wiki.galaxyproject.org/Develop/Visualizations

# Part 1

The visualization we want to write is a tool that shows the number of aligned reads per chromosome, of a BAM file.
The first thing we need to do is to come up with a name.
Let's call it *alignment_rname_boxplot*.
Remark that the chromosomes to which we align are named RNAME in the BAM/SAM specification.

When we start our development, we do that within the Galaxy codebase.
After cloning an instance of galaxy in $GALAXY_ROOT, we go into the plugin directory:

  cd $GALAXY_ROOT/config/plugins/visualizations ;

Here we need to make a new directory for our new plugin:

  mkdir alignment_rname_boxplot ;
  cd alignment_rname_boxplot ;

We need to make three directories to complete the structure of the project:

  mkdir config ;
  mkdir static ;
  mkdir templates ;

## Linking the plugin with Galaxy

To create a bridge between the not yet written plugin and Galaxy, we need to write a configuration in XML format.
Create the file:  `config/alignmnet_rname_boxplot.xml` and give it the following syntax:

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

This configures the name, which shall appear on hovering the visualization button in the history menu.
It also links the plugin to two file formats: BAM and SAM.
This means that for any history item in Galaxy of those file formats, the plugin will become available once they're in a history.
The *var_name_in_template* is set to hda, which will be the name of the variable corresponding to the visualized dataset.
The last configuration is the template, which is a mako template (HTML + Python syntax) found in the ./templates directory.

## Creating the visualization

We have just linked Galaxy to a mako file, which does not yet exist.
Please create a mako file: ./templates/alignment_rname_boxplot.mako
This file is a blueprint for the visualization. 
This means that for every invocation of the visualization, the mako file will be compiled to render a HTML file.
It is trivial to understand that compiliation needs to take place in a few ms because otherwise the visualization becomes slow.
This at itself is not a problem, as the visualizations are supposed to be written in HTML and/or JS and that most of the calculation takes place at the client side (in the browser).
Another thing, that is fundamentally different from developing static Galaxy tools, is that access to the data is not via python anymore.
Instead, the data will be downloaded / queried by the client via a URL and parsing happens in the browser via JS libs too.
This may sound pretty scary but we will explain later on.

We start a mako file with the variables we need:

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

This allows us to query the BAM/SAM file by adding the *file_url* into a HTML section.
However, BAM files can become really large and undesired to pump over the network.
Additionally, it is also inconvenient to parse the BAM file via JS to count the number of reads.

Fortunately, BAM files (should) have indices.
These indices are brief summaries describing the number of entries per chromosome, to access them more quickly.
In the python part of the mako temolate we can access an index as *hda.metadata.bam_index*.
Please remark that this is a path in the server and not an exposed URL.
It is probably not that difficult to extract the URL of the index.
However, this file is binary compressed and requires more advanced programming to extract the information from it which is not the scope of this practical.
Samtools idxstats is a program in samtools that is able to do this advanced binary decompression for you.
However, since visualizations do not have dependency management, it is very tricky to let the mako template do a system call to samtools.
Fortunately, the galaxy ecosystem has a builtin pysam dependency, a library that can do any native samtools command within python.

The *bam_index* is a special kind of file. 
It is actually an invisible file in Galaxy, but does have its own unique filename.
So, for the BAM file ./database/files/000/dataset_001.dat, our BAI file *is NOT* ./database/files/000/dataset_001.dat.bai but could be ./database/files/000/dataset_002.dat or ./database/files/000/dataset_003.dat.
To solve this in a safeway, we can create a symlink - if it does not already exists - to make sure they have the same prefix  to make them compatible with pysam.
This trick is rather similar to how it is applied in static Galaxy tools.
We can make this symlink as follows:

    ## Ensure BAI index is symlinked
    bai_target = hda.file_name+'.bai'
    import os
    
    if not os.path.isfile(bai_target):
        os.symlink(hda.metadata.bam_index.file_name, bai_target)

Now the BAM file has a BAI file with the same name prefix, we can extract the idxstats as follows:

    import pysam
    data = pysam.idxstats(hda.file_name)


We could decide to remove the symlink after using it, but may have deeper implications than I am aware off, but it's definitely worth considering it.
Currently the idxstats data is present in the virtual memory of python during compilation, but is not exported to HTML or JS.
To figure out what the data looks like, let's proceed with the following code:

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
        </head>
        <body>
            ${bam_idxstats_data | h}
        </body>
    </html>
    

We are now going to test it, but we need a (small) BAM file for it.
You can find small examples in the test files of the aligners in tools IUC (hisat2, rna-star).
Let's proceed with the following file:
https://github.com/galaxyproject/tools-iuc/raw/master/tools/hisat2/test-data/hisat_output_1.bam
Go the the galaxy root directoryt (cd $GALAXY_ROOT) and run galaxy: ./run.sh .
Upload the BAM file and take a look at the visualization.
Once you open the visualization you will see the contents of idxstats.
It is python dump of a list: ['phiX174\t5386\t19\t1\n', '*\t0\t0\t0\n']
It contains two entries, one to phiX and one to *, of which the former is a chromosome and the latter are unmapped.
The entry of phiX is tab delimited and indicates that the length of the RNAME is 5386 bases and 19 reads are aligned to it.

The trick we need to do now is to make this data ready to be read by Javascript and to parse it to a  simple dictionary of the following syntax:
    {'phiX': 19, '*': 0}
This could be done in python but I would recommend to do this in JS.
The var dump in python is actually a valid syntax for Javascript too.
So getting the raw data into javascript if rather easy:  `bam_idxstats_data = ${bam_idxstats_data};`
Parsing is not the scope of the tutorial, so here it is:

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

The great thing about the mako system is that it does not require to restart galaxy in order to test changes to the mako files.
So, please open our visualization in a new tab, so we can press [F5] to test changes.
Change the file to the following:


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

Reload the page and open the developers console in the browser.
Please ask help if you don't know how to do this, this setting is different for every browser.
Within the console, type: `idx_stats` and press enter.
This should give the parsed contents as a dictionary, which can directly be used in javascript.

From this point I would like you to continue on your own to see if you are able to create a simple visualization of this dictionary.
You may think of tables, DIVs or even more complicated solutions :)

If you want to check out the solution I made, please have a look at:
[../plugins/alignment_rname_boxplot](alignment_rname_boxplot)

In the given example, the RNAME queries containing an underscore were removed.
This is because there are many alterative chromosomes, making the list very large for certain reference genomes.
However, for certain studies it might be desired to just look at those.
Also,  the current plot is a boxplot, but I can imagine that a pie-chars can be convient too.
All of those can be implemented for interactive change, contributing to quicker understanding of the data which is generally not feasible using static Galaxy tools.

In the example we had included the Javascript and CSS into the HTML website.
This at itself is not a problem, but remember that for every new invocation of the visualization the entire CSS en JS are copied and transferred as well.
This is a waste of bandwith as we could save the files in the static directory and refer to them within the HTML.
The browser shall check it's cache for the presence of libs and stylesheets and only update them if they have changed.

Another thing that you may realize is that we still do the calculation (pysam.idxstats) server side.
Although this is a marginal calculation, it is causing delay and the bigger the files get the worse.
The underlying problem here is that we did not obtain and parse the BAI file via javascript.
This would be a more elaborate solution, but requires more development time.

Another thing we could do is create a separate `samtools idxstats` tool, that creates a file of datatype `txt_idxstats` and include that datatype into Galaxy.
Then we make the visualization specific for that datatype, just plotting the results of the idxstats.

A more complicated problem is that BAM files are simply too big to transfer for this kind of applications.
It would be ideal to have some webserver integration that allows to query specific locations within a BAM file while indexing is done server side.
This is exaxtly what has been done in Trackster.

# Conclusion

We have created a visulization plugin in Galaxy, that is able to visualize the number of alignments per RNAME in a BAM file.

:grey_exclamation: ***Key Points***

- *Visualizations require a different way of thinking*
- *Interactivity is what makes visualizations different from static tools*
- *Requires understanding of both the Galaxy ecosystem as well as HTML5/JS*
- *Performance is much more important than for static Galaxy tools*

# :clap: Thank you
