---
title: Downloading datasets using command line
area: datasets
box_type: tip
layout: faq
contributors: [jennaj, beachyesh]
---

From the terminal window on your computer, you can use **wget** or **curl**.

1. Make sure you have **wget** or **curl** installed.  
2. Click on the Dataset name, then click on the copy link icon {% icon galaxy-link %}. This is the direct-downloadable dataset link. 
3. Once you have the link, use any of the following commands:   
   - For **wget**
   >**``wget '<link>'``   
   > `` wget -O '<link>'``  
   > ``wget -O --no-check-certificate '<link>'  # ignore SSL certificate warnings``   
   > ``wget -c '<link>' # continue an interrupted download``**   
   - For **curl**
   > **``curl -o outfile '<link>' ``  
   > ``curl -o outfile --insecure '<link>'     # ignore SSL certificate warnings``  
   > ``curl -C - -o outfile '<link>'           # continue an interrupted download``**
4. For dataset collections and datasets within collections you have to supply your API key with the request
   - Sample commands for **wget** and **curl** respectively are:
   > 
   > **``wget https://usegalaxy.org/api/dataset_collections/d20ad3e1ccd4595de/download?key=MYSECRETAPIKEY``**
   > 
   >**``curl -o myfile.txt https://usegalaxy.org/api/dataset_collections/d20ad3e1ccd4595de/download?key=MYSECRETAPIKEY``**
 
