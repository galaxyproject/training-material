---
title: Using IGV with Galaxy
area: visualisation
box_type: tip
layout: faq
contributors: [shiltemann]
---

You can send data from your Galaxy history to IGV for viewing:

1. Install IGV on your computer ([IGV download page](https://software.broadinstitute.org/software/igv/download))
2. Start IGV
3. In recent verions of IGV, you will have to enable the port:
   - In IGV, go to `View > Preferences > Advanced`
   - Check the box `Enable Port`

3. In Galaxy, expand the dataset you would like to view in IGV
   - Make sure you have set a reference genome/database correctly (dbkey) ([instructions]({% link faqs/galaxy/datasets_change_dbkey.md %}))
   -  Under `display in IGV`, click on `local`
