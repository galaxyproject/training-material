---
username: tklingstrom
layout: training_philosophy
---

I have been working in three very different training environments.

1. Lecturing with a biology students in an applied bioinformatics course.
2. Classroom training with biobank experts working to create a more comprehensive service for users.
3. One of one with system administrators expanding into the technical part of bioinformatics support.

I agree with everything said above so I will try to add some more specialised flavour. I would also like to add that I think the Galaxy Training resources are absolutely fantastic and if I get the opportunity to hold a more comprehensive course (ie more than half a day) I will structure it so that we have one quick joint introduction and then conduct a number of interactive tutorials using the training portal.

**Case 1 students**
- For students context is king. I think almost everyone how many open resources there are available in the world and this is especially true for students. Knowing what is available for download and/or through the "Get data" tab of Galaxy is therefore a valuable exercise in itself, likewise informing the students about key resources such as the GTN, Biostars and Stack Overflow is very important.
- Students often spend so much time writing/clicking commands that they forget about the biological context. As a teacher I think it therefore is very important to help students take a step back and ask themselves why they do a step. For example, after doing the steps to assembly RNA-seq reads it is good to break the work, show an assembly of an mRNA and help them remember what we have done and then explain how the next steps are to map it to the genome so that we can identify it in an unambiguous manner before doing the statistics.

**Case 2 experts**
- Catching the attention of experts (in another field than bioinformatics) is often a bit of a challenge. Not due to a lack of interest but because their focus is often on a narrow but very well defined area. Meaning that it is often very valuable to pre-plan exercises to be directly applicable within their own research. Ideally I would like to work give a training in such a situation that even the "get data step" is based on their local environment but that is unlikely to be possible in most cases.

**Case 3 system administrators**
During our set up in Gambia I had a very positive experience regarding using the Dockerized Galaxy set up. Configuring environmental variables and then just running the docker run command gives a very good insight in how to adapt the system for a local environment without running the risk of ruining your environments. At SLU I am now working to take this a step further with a Nextcloud implementation that provide us with a back end for handling FTP data to and from the docker FTP server.
