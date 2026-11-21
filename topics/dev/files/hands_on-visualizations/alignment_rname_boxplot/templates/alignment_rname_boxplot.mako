<!DOCTYPE HTML>
<%
    import os

    ## Generates hash (hdadict['id']) of history item
    hdadict = trans.security.encode_dict_ids( hda.to_dict() )

    ## Finds the parent directory of galaxy (/, /galaxy, etc.)
    root     = h.url_for( '/' )

    ## Determines the exposed URL of the ./static directory
    app_root = root + 'plugins/visualizations/'+visualization_name+'/static/'

    ## Actual file URL:
    file_url = os.path.join(root, 'datasets', hdadict['id'], "display?to_ext="+hda.ext)

    ## Ensure BAI index is symlinked
    bai_target = hda.file_name+'.bai'

    if not os.path.isfile(bai_target):
        os.symlink(hda.metadata.bam_index.file_name, bai_target)

    ## Extract idxstats
    import pysam
    bam_idxstats_data = pysam.idxstats(hda.file_name)
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

            function calc_stats(parsed) {
                max = 0;
                sum = 0;
                for (var key in parsed) {
                    if (parsed[key] > max){
                        max = parsed[key];
                    }
                    sum += parsed[key]
                }
                return [max, sum];
            }

            function plot_data(parsed) {
                var max = calc_stats(parsed)[0];
                var sum = calc_stats(parsed)[1];

                for (var key in parsed) {
                     var value = parsed[key];
                     var ratio = 100.0 * value / sum;
                     var ratio2 = 100.0 * value / max;

                     var div = document.createElement("div");
                     div.innerHTML = '<nobr>'+key+'</nobr>';
                     document.getElementById("chart_names").appendChild(div);

                     var div = document.createElement("div");
                     div.innerHTML = '<nobr>'+value+" ("+Math.round(ratio*100)/100+"%)</nobr>";
                     div.title = key+': '+value+" ("+Math.round(ratio*100)/100+"%)";
                     div.style.width =  ratio2+'%';
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
