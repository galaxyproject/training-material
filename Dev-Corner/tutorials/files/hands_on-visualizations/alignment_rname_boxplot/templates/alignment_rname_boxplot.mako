<!DOCTYPE HTML>
<%
    import os
    import tempfile
    import pysam

    ## Generates hash (hdadict['id']) of history item
    hdadict = trans.security.encode_dict_ids( hda.to_dict() )

    ## Finds the parent directory of galaxy (/, /galaxy, etc.)
    root     = h.url_for( '/' )

    ## Determines the exposed URL of the ./static directory
    app_root = root + 'plugins/visualizations/'+visualization_name+'/static/'

    ## Actual file URL:
    file_url = os.path.join(root, 'datasets', hdadict['id'], "display?to_ext="+hda.ext)

    # File path in tmp dir:
    temp_data_dir = tempfile.gettempdir()+'/plugins/visualizations/'+visualization_name
    alignment_file = temp_data_dir + '/' + hdadict['id'] + '_' + os.path.basename(hda.file_name)
    bai_target = alignment_file + '.bai'

    ## Ensure temp dir exists
    if not os.path.exists(temp_data_dir):
        os.makedirs(temp_data_dir)

    ## Make symlink to datafile in tmp dir
    if not os.path.isfile(alignment_file):
        os.symlink(hda.file_name , alignment_file)

    ## Make symlink to datafile in temp dir (to keep Galaxy storage dir clean)
    if not os.path.isfile(bai_target):
        if not os.path.isfile(hda.metadata.bam_index.file_name):
            pysam.index(alignment_file)
        else:
            os.symlink(hda.metadata.bam_index.file_name, bai_target)

    ## Extract idxstats
    bam_idxstats_data = pysam.idxstats(alignment_file)
%>
<html>
     <head>
        <title>${hda.name | h} | ${visualization_name | h}</title>

        <style>
             .chart div {
               font: 10px sans-serif;
               background-color: steelblue;
               text-align: right;
               padding: 3px;
               margin: 1px;
               color: white;
             }
        </style>
        <script>
            bam_idxstats_data = ${bam_idxstats_data};
            function parse_data(data) {
                 /*
                  Data comes in as tuple of unsplit lines:
                  ["chr1\t1000\t0\t0", "chr2\t2500\t0\t0"]

                  We need to split it up, and ideally only keep reference names without an underscore
                 */
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

            function calc_max(parsed) {
                output = 0;
                for (var key in parsed) {
                   if (parsed[key] > output){
                     output = parsed[key];
                 }
                }
                return output;
            }

            function plot_data(parsed) {
                var max = calc_max(parsed);

                for (var key in parsed) {
                     var value = parsed[key];
                     if (max > 0) {
                        var ratio = 100.0 * value / max;
                     }
                     else {
                        var ratio = 100.0; // fixes null division
                     }

                     var div = document.createElement("div");
                     div.innerHTML = '<nobr>'+key+'</nobr>';
                     document.getElementById("chart_names").appendChild(div);

                     var div = document.createElement("div");
                     div.innerHTML = '<nobr>'+value+" ("+Math.round(ratio*100)/100+"%)</nobr>";
                     div.title = key+': '+value+" ("+Math.round(ratio*100)/100+"%)";
                     div.style.width =  ratio+'%';
                     document.getElementById("chart").appendChild(div);
                }
            }
        </script>
     </head>
     <body onload="plot_data(parse_data(bam_idxstats_data));">
         <center>
             <h1>Bam contents of ${hda.name | h}</h1>

             <table border="0" borderpadding="0" borderpanning=")" style="width: 500px;">
                 <tr>
                     <td style="width:50px;">
                         <div id="chart_names" class="chart" style="width: 100%; border: 1px dashed gray;text-align: left;" />
                     </td>
                     <td style="width:450px;">
                         <div id="chart" class="chart" style="width:100%; border: 1px dashed gray;text-align: left;" />
                     </td>
                 </tr>
             </table>
         </center>
     </body>
</html>
