---
title: What advantages does a Chromatogram Library have over a DDA-generated library or predicted spectral library?
box_type: question
layout: faq
contributors: [emmaleith]
---

While generating a Chromatogram Library is the most time consuming step of the EncyclopeDIA workflow, it is beneficial to DIA data analysis. DIA is a novel technique and methods for DIA data analysis are still being developed. One method commonly used includes searching DIA data against DDA-generated libraries. However, there are limitations in this method. Firstly, DDA-generated libraries are not always an accurate representation of DIA data: differences in the methods of data collection play an important role in the efficacy of the library. Secondly, DDA-generated libraries often require labs to run completely separate DDA experiments to simply generate a library with which to analyze their DIA data.
Chromatogram Libraries mitigate some of the previous shortcomings mentioned. DIA data is incorporated into the generation of the Chromatogram Library and therefore provides context to the DIA data being analyzed. Secondly, the ELIB format of the Chromatogram Library allows for extra data to be included in the analysis of the DIA data, including intensity, m/z ratio, and retention time compared to the use of a DDA-generated DLIB library. Lastly, a Chromatogram Library can be generated without the use of a spectral library (as mentioned in the last question). Therefore, it is possible to forgo DDA data collection as the DLIB DDA-generated library is not strictly needed for Chromatogram Library generation and to run the EncyclopeDIA workflow (saving time and resources).

