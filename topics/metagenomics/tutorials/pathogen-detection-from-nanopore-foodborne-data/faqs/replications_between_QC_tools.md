---
title: Between Nanoplot and FAStQC or FastP, is there a recalibration done?
box_type: question
layout: faq
contributors: [EngyNasr]
---

we can reproduce this discrepancy. we observe there are 2 tables output by Nanoplot: stats vs stats post filtering, we suspect something is happening there. But even that does not seem to fully explain it, because fastqc quality distribution pre-processing is:

|Quality| Count |
|-------|-------|	
|11     |11.0   |
|12     |93.0   |
|13     |580.0  |
|14     |1713.0 |
|15     |3188.0 |
|16     |4893.0 |
|17     |6953.0 |
|18     |9377.0 |
|19     |12743.0|
|20     |15741.0|
|21     |17226.0|
|22     |15567.0|
|23     |11749.0|
|24     |7809.0 |
|25     |4183.0 |
|26     |1810.0 |
|27     |660.0  |
|28     |261.0  |
|29     |97.0   |
|30     |66.0   |
|31     |46.0   |
|32     |30.0   |
|33     |28.0   |
|34     |24.0   |
|35     |16.0   |
|36     |18.0   |
|37     |14.0   |
|38     |13.0   |
|39     |14.0   |
|40     |15.0   |
|41     |6.0    |
|42     |8.0    |
|43     |4.0    |
|44     |14.0   |
|45     |5.0    |
|46     |1.0    |
|47     |4.0    |
|48     |2.0    |
|49     |0.0    |
|50     |3.0    |
|51     |0.0    |
|52     |0.0    |
|53     |0.0    |
|54     |1.0    |

