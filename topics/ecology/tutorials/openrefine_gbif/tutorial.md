---
layout: tutorial_hands_on

title: Cleaning GBIF data using OpenRefine
zenodo_link: ''
questions:
- How can I use OpenRefine to clean data?
- How do I check and clean biodiversity data using OpenRefine?

objectives:
- Use OpenRefine faceting functionalities to apply mass editing and manage duplicates
- Use OpenRefine clustering and filtering functionalities to edit, transform data notably using regular expression
- Use OpenRefine to apply API services on your data
time_estimation: 0h45m
key_points:
- OpenRefine is a powerful tool, with many functionalities to check, clean and enrich your data
contributions:
    authorship:
        - yvanlebras
        - sophiepamerlon
        - laurarussell
    reviewing:
        - Sch-Da
    funding:
        - gbif
        - pndb
tags:
  - biodiversity
subtopic: ecologymetadatamgt
---


In this tutorial we will use OpenRefine tool to clean occurrence records retrieved from GBIF (Global Biodiversity Information Facility).

This tutorial is based on the GBIF {% cite GBIFtraining2021 %} "Biodiversity data mobilization course" accessible on [GBIF website](https://docs.gbif.org/course-data-mobilization/en/) on Use case 1 ["Herbarium Specimens"](https://docs.gbif.org/course-data-mobilization/en/exercise-3a-c.html) using file from Exercise 3C, at the bottom of the page, reachable through this direct [URL](https://docs.gbif.org/course-data-mobilization/exercise-data/UC1-3c-open-refine.csv).

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Retrieve data from GBIF "Biodiversity data mobilization course"


> <hands-on-title> Data import and project creation </hands-on-title>
>
> 1. Create a new history for this tutorial and give it a name (example: “OpenRefine tutorial”) for you to find it again later if needed.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
> 2. Import data file from GBIF {% cite GBIFtraining2021 %} "Biodiversity data mobilization course" accessible on [GBIF website](https://docs.gbif.org/course-data-mobilization/en/) on Use case 1 ["Herbarium Specimens"](https://docs.gbif.org/course-data-mobilization/en/exercise-3a-c.html) using file from Exercice 3C, at the bottom of the page
>
>    "Herbarium Specimens" data file as a csv file named "UC1-3c-open-refine.csv"
>    ```
>    https://docs.gbif.org/course-data-mobilization/exercise-data/UC1-3c-open-refine.csv
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>
> 3. {% tool [Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regex1/1.0.3) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `UC1-3c-open-refine.csv`
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `;`
>            - *"Replacement"*: `\t`
>    - click button "run tool"   
> 4. Change the datatype from "txt" to "tabular"
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tabular" %}
>
> 6. Check your output. All columns must be separated by a tabulation now. See example below.
>
>    ![Output of Regex Find And Replace example](../../images/openrefine/1_openrefine.png "Output of Regex Find And Replace example"){: style="width:50%"}
>
>
> 7. Rename the datasets if needed, for example with "UC1-3c-open-refine-tabular" to keep orginal name and specify the tabular format.
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
>
>
{: .hands_on}


# Deploy an OpenRefine instance and push your data in


> <hands-on-title> Deploy an OpenRefine instance and import data from your Galaxy history </hands-on-title>
>
> 1. {% tool [OpenRefine tool](interactive_tool_openrefine) %} with the following parameters:
>    - *"Input file in tabular format"*: `UC1-3c-open-refine-tabular`
>    - click the button "run tool"
> You will then be redirected. A notification "Interactive Tools" shows up on the left-hand side, indicating you have an interactive tool running.
> Click on the notification, click on the empty box in front of the text saying "Openrefine".
> Then, click on the text "Openrefine" to open the tool in a new tab. It might take a moment until your cursor changes and the text becomes a hyperlink.
>
> 2. Open the OpenRefine instance
>    - go to "Open Project" space
>    - select the project named "Galaxy file"
> 3. You will now have your project with first 10 lines of your file displayed
> ![OpenRefine project start](../../images/openrefine/openrefine_start.PNG)
>
{: .hands_on}


# Faceting

Faceting is a feature that will allow us to get a big picture overview of the data, and to filter down to just the subset of rows
that we want to change or view in bulk. It facilitates the use and analysis of data and can be done with cells containing any
kind of text, numbers and dates.

> <hands-on-title> Faceting and mass editing </hands-on-title>
>
> 1. Go to column "kingdom", and then click on the column menu and follow the route to `Text facet`.
>    - On the left a window with the name of the column will appear, that is the facet.
>    ![OpenRefine facet kingdom](../../images/openrefine/openrefine_facet_kingdom.PNG)
> 2. Click on count to sort by count, then click on name to sort alphabetically
> 3. Fix the spelling mistakes (Plante -> Plantae). Place the cursor over the text in the window and click on `edit`, then fix the error in the text box, and to save click on `apply`.
>    - All the values will be fixed automatically.
>
{: .hands_on}

> <question-title></question-title>
>
> 1. How many Plantae occurences you now have?
>
> > <solution-title></solution-title>
> >
> > 1. 736
> >
> {: .solution}
>
{: .question}

> <hands-on-title> Faceting and white space 1 </hands-on-title>
>
> 1. Go to "Country col." column and click on column menu and perform a `Text Facet`
>    - On a quick view, the country appears to be spelled correctly, but the facet shows three different values due to the extra spaces at the end of the text.
> 2. Fix the error from the column menu on "Country col." column, following the route by clicking the dropdown menu, `Edit Cells > Common transforms > Trim leading and trailing whitespace`.
>    - You will see a notification message "Text transform on 38 cells in column Country col.: value.trim()" and the possibility to come back on previous state clicking on the `undo` hyperlink.
> 3. Now check the facet window; only one value will remain.
>
{: .hands_on}

> <hands-on-title> {% icon hands_on %} Faceting and white space 2 </hands-on-title>
>
> 1. Go to column "Full name" and click on the dropdown menu, then go to `Text facet`. Then click on `count`.
>    - *Guzmania lingulata* is the first item in the list with 25 specimens, but it is also present in the 3rd position with 20 specimens.
> 2. Fix the error from the "Full name" column menu, `Edit Cells > Common transforms > Collapse consecutive whitespaces`.
>
{: .hands_on}

> <question-title></question-title>
>
> 1. How many *Guzmania lingulata* records only appear in the list finally?
>
> > <solution-title></solution-title>
> >
> > 1. 45
> >
> {: .solution}
>
{: .question}

> <hands-on-title> {% icon hands_on %} Faceting and duplicates </hands-on-title>
>
> 1. Go to column catalog in "Cat. Numb", and follow the route `Facet > Customized facets > Duplicates facet`.
>    - The facet will show 4 duplicates
> 2. Click on true, and you'll see the values in the main window
> 3. After a check with the specimens labels, fix the values clicking edit directly on the cell with the correct catalogue numbers
>    - UWP:122470 for Vargas P
>    - UWP:122471 for Vargas I
>    - UWP:157351 for Betancur H
>    - UWP:157339 for Betancur J
>    
>     click `apply` after you finish each edit
{: .hands_on}

# Filtering

> <hands-on-title> {% icon hands_on %} Basic filter </hands-on-title>
>
> 1. Go again to "Full name" column menu and perform a `Text facet` to visualize the values
> 2. Then go again to the column menu and click on `Text filter`
> 3. Perform the following filters and fix them as described below:
>      - search for "sp1" entries -> Then remove it and obtain "Cyperus", clicking on `Edit` directly in the cell
>      - search for "SP2" entries, check `case sensitive` -> To remove it and obtain "Cyperus", you can `Edit` directly in the cell
>      - search for "spp".
>         - Go to "Full name" column menu, then click `Edit cells > Transform`
>         - In the text box paste the formula `value.replace(" spp.", "")` and click `ok`
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Do you think you can use also `Edit cells > Transform` formula to apply the changes you made manually for the two first terms searched ("sp1" and "SP2")?
> 2. Why can you easily edit these two first entries directly but not the third one?
>
> > <solution-title></solution-title>
> >
> > 1. Yes, using respectively `value.replace(" sp1", "")` and `value.replace(" SP2", "")` formulas
> > 2. Because the first two searches only concern one entry for each when the third one returns several.
> >
> {: .solution}
>
{: .question}

> <hands-on-title> {% icon hands_on %} Advance filter 1 </hands-on-title>
>
> 1. Go to column "genus" and perform a `Text filter`.
> 2. Check `regular expression` and `case sensitive` , then paste the expression "^[a-z]"
>    - This regular expression filters the strings in which the first letter is lowercase.
> 3. Perform a correction since the genus should be capitalized.
>
>
>    > <comment-title></comment-title>
>    >
>    > If you want to know more about regular expressions, you can consult the [Openrefine regular expressions documentation](https://openrefine.org/docs/manual/expressions#regular-expressions) and specifically the [documentation for GREL function](https://openrefine.org/docs/manual/grelfunctions)
>    {: .comment}
>
>    {% snippet faqs/galaxy/analysis_regular_expressions.md %}
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Try to find a regular expression to apply the same kind of `Text filter` operation as in the previous exercise. Which regular expression can work?
>
> > <solution-title></solution-title>
> >
> > 1. On my side, I found this formla `value.replace(value.substring(0,1),toUppercase(value.substring(0,1)))`
> >
> {: .solution}
>
{: .question}

> <hands-on-title> {% icon hands_on %} Advance filter 2 </hands-on-title>
>
> 1. Go to column "Full name" and perform a `Text filter`.
> 2. Check `regular expression` and `case sensitive`, then paste the expression `^[A-Z].*\s[A-Z]`
>    - This regular expression filters the strings that start with a capital letter followed by 0 or more characters, then a space, then a capital letter.
> 3. Perform a correction since the second word of the name should be lowercase.
>
{: .hands_on}


# Clustering

> <hands-on-title> Basic clustering </hands-on-title>
>
> 1. Go to "County" and perform a `Text facet`.
>    - Keep in mind that the correct counties are: "Flores", "La Libertad", "Melchor de Mencos", "San Andres" and "San Jose".
> 2. On the top right of the facet window click on `Cluster`, a new window will appear.
> 3. Click on the `Cluster` button from this new window.
> 4. Now you can see information about the clusters:
>    - **"Cluster size":** the number of different versions that the clustering algorithm believes to be the same.
>    - **"Row count":** the number of records with any of the cluster values.
>    - **"Values in cluster":** the actual values that the algorithm believes to be the same. There is also the number of records with each particular value, and the possibility to browse the contents of the cluster in a different tab.
>    - **"Merge?":** check if values are to be merged into a single standard value.
>    - **"New cell value":** the value to be applied to every record in the cluster. By default, it is the value with most records. You can also click on any value to apply that to the New cell value.
>
>
>    > <comment-title></comment-title>
>    >
>    > If you want to know more about clustering have a look at the [manual](https://openrefine.org/docs/manual/cellediting#cluster-and-edit)
>    {: .comment}
>
> 5. Click on `Select All` and then on `Merge Selected & close`
>    - You will see a notification message "Mass edit 119 cells in column County".
> 6. To fix the remaining counties go again to Cluster in the facet window of Count.
> 7. In the *Cluster and edit* window, go to `Keying Function`, then select `ngram-fingerprint`
>    - set "1" as the value in `n-Gram Size`.
> 8. Press the `Cluster` button, you normally see a cluster about "San Andres" of size "4".
> 9. Click on `Select All` and then on `Merge Selected & close`, you will see a notification message "Mass edit 360 cells in column County".
>    - Your counties are now fixed! Congratulation!
{: .hands_on}


# Exporting

> <hands-on-title> {% icon hands_on %} Exporting cleaned file into your Galaxy history </hands-on-title>
>
> 1. On the upper right corner click on `Export` and select `Galaxy exporter`.
>    - A notification message as "Dataset has been exported to Galaxy, please close this tab" is displayed.
>    - You normally have your resulting data file exported on your Galaxy history as "openrefine-Galaxt file.tsv" dataset.
>
>    - When you want to close OpenRefine, go back to the Galaxy window. Go to the middle window, "Active Interactive Tools". Below `Openrefine` click `stop`.
>    - The interactive tool is then stopped and the respective task shows red in your history. This does not affect your data, which is safely stored in the "openrefine-Galaxt file.tsv" dataset.
>    - You can continue working with it in Galaxy for the following tasks.
>
>    > <comment-title></comment-title>
>    >
> > You can also download the file directly on your computer, choosing `Export` and `Custom tabular exporter`, allowing you to select notable columns and rows you want to download.
>    {: .comment}
>
{: .hands_on}


# API use

Reconciliation matches the information in one of your columns to an outside database. This is particularly helpful when it
comes to name validation, as it proves the name you have exists somewhere else. This is a really useful service but can be
time-consuming. In this case, we will go through the process with only three records using the API from GBIF.
To open your dataset again from Galaxy in OpenRefine, follow the earlier steps from this tutorial "Deploy an OpenRefine instance and push your data in".

> <hands-on-title> {% icon hands_on %} Higher taxonomy </hands-on-title>
>
> 1. In your opened dataset in OpenRefine, go to the "Collector" column, then make a `Text facet`. Select the collector "Elsa P".
> 2. Under "Full name", click on column menu and then `Edit column > Add column by fetching URLs…`
>    - call the new column "Api_name"
> 3. Change the Throttle Delay to 250 and paste the expression `"http://api.gbif.org/v1/species/match?verbose=true&name="+escape(value,'url')`
> 4. Click ok and wait, this might take some time depending on internet connection and the number of taxa.
> 5. Go to "Api_name", click on column menu and then `Edit column > Add column based on this column...`.
>    - Call the new column "higherClassification" and paste the expression:
>
>    ```
>    value.parseJson().get("kingdom")+
>    ", "+value.parseJson().get("phylum")+
>    ", "+value.parseJson().get("class")+
>    ", "+value.parseJson().get("order")+
>    ", "+value.parseJson().get("family")
>    ```
>   click `ok`. 
>   You will see the Kingdom, Phylum, Class, Order and family of each taxon.
> 6. Under "higherClassification" follow the route `Edit column > Split into several columns…`, leave the initial settings and click `ok`.
> 7. Now you know how to obtain the taxonomic categories of a given taxon if this is available in the GBIF API. Column names can be edited in `Edit column > Rename this column`.
> 8. For the purpose of the original GBIF workshop, the columns created in this exercise (Higher taxonomy) must be deleted.
>    -  Under All, which is the first column, go to `Edit columns > Re-order / remove columns…`.
>    - Remove columns "Api_name", "higherClassification 1", "higherClassification 2", "higherClassification 3", "higherClassification 4" and "higherClassification 5".
> 9. No need to export this file as it normally comes back to the previous version you already exported.
>
{: .hands_on}



# Conclusion

Here you learned how to use OpenRefine tool from Galaxy platform to clean Biodiversity data. This tutorial notably allowed you to apply some basic but powerful functionalities of OpenRefine to clean your data.
