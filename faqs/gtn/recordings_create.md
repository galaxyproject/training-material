---
title: Recording a video tutorial
area: contributors
layout: faq
box_type: tip
contributors: [shiltemann, hexylena]
---

This FAQ describes some **general guidelines** for recording your video

**Anybody is welcome to record one of the GTN tutorials**, even if another recording already exists!
Both the GTN tutorial and Galaxy itself change significantly over time, and having regular and/or multiple recordings of tutorials is great!

**Done with your recording?** Check out the instructions for adding it to the GTN:

[Submitting Recordings to the GTN]({% link faqs/gtn/recordings_add.md %}){: .btn.btn-info}


#### Video content

1. **Start of video**
  - **Introduce yourself**
  - Discuss the **questions and learning objectives** of the tutorial
  - Give a **basic introducion about the topic**, many participants will be novices

2. **Guide the learners through the tutorial step by step**
   - Explain the scientific background of the analysis
   - Explain where you are clicking in Galaxy
   - Explain what tool parameters mean
   - Explain what the tool does
   - Discuss the output files
   - Discuss how to interpret the results
   - Discuss question boxes from the tutorial

3. **Speak slowly and clearly**
  - Take your time, we are not in a hurry
  - It is often a lot of new information for participants, give them a chance to process all of it
  - Speaking slowly and clearly will improve the quality of the auto-generated captions, and will be less work for you to fix captions.

4. **If things go wrong that is OK!**
  - It's a great teaching moment!
  - Explain the steps you are taking to determine what went wrong, and how you are fixing it.
  - It makes participants feel less bad if things go wrong for them

5. **If your tutorial is long**
  - Indicate good places for people to take a break
  - e.g. when a tool takes a while to run

6. **End of video**
  - Go over some of the **take-home messages (key-points)** of the tutorial
  - Remind viewers about the **feedback form** embedded at the end of the tutorial
  - Share your recommendations for **follow-up tutorials**
  - Share any other tips for where to learn more about the topic
  - Share **how to connect with the community** (e.g. Matrix, Help Forum, social media, etc)

7. If you are doing both a lecture and a hands-on training, please create 2 separate videos


#### Technical Guidelines

1. Start a [Zoom](https://zoom.us/) call with yourself, record that.
   - For Mac users, QuickTime Player is also a nice option.
   - Have another preference like OBS? Totally OK too!
   - We recommend zoom to folks new to video production as it is the easiest to get started and produces quite small file sizes.

2. Do a short **test recording** first
   - Is the **audio quality** good enough?
     - Wearing a headset often improves the audio quality.
   - **Screen sharing:** is your screen readable?
     - Make sure you **zoom in** enough for it to be clearly visible what you are doing in Galaxy.
     - Test watching the video in a non-maximised window. Is it still legible?
     - If the participant is using 50% of their screen for the video, 50% for Galaxy, will it be legible?

3. Need to edit your video after recording?
  - For example to merge multiple videos together?
  - Software like [KDEnlive](https://kdenlive.org/en/) can help here.
  - Feel free to ask us for help if you need!


#### Standards

1. **Zoom in**, in every interface you're covering! Many people will be watching the video while they're doing the activity, and won't have significant monitor space. Which video below would you rather be trying to follow?

   Bad | Good 😍
   --- | ---
   ![default size screenshot of usegalaxy.eu]({% link faqs/gtn/images/bad.png %}) | ![zoomed in screenshot of usegalaxy.eu, now much more legible]({% link faqs/gtn/images/good.png %})

   Bad | Good 🤩
   --- | ---
   ![green text on black background console with tiny font]({% link faqs/gtn/images/bad-console.png %}) | ![zoomed in screenshot of a console with high contrast black and white content]({% link faqs/gtn/images/good-console.png %})

2. (Especially for introductory videos!) Clearly call out what you're doing, especially on the first occurrence

   Bad | Good
   --- | ---
   "Re-run the job" | "We need to re-run the job which we can do by first clicking to expand the dataset, and then using the re-run job button which looks like a refresh icon."

   Bad | Good
   --- | ---
   "As you can see here the report says X" | "I'm going to view the output of this tool, click on the eyeball icon, and as you can see the report says X."

   But the same goes for terminal aliases, please disable all of your favourite terminal aliases and quick shortcuts that you're used to using, disable your bashrc, etc. These are all things students will try and type, and will fail in doing so. We need to be very clear and explicit because people will type exactly what is on the screen, and their environment should at minimum match yours.

   Bad | Good
   --- | ---
   `lg file`| `ls -al | grep file`
   `z galaxy`| `cd path/to/the/galaxy`

3. Consider using a pointer that is more visually highlighted.

   ![mouse pointer with circle around it that follows it around]({% link faqs/gtn/images/mouse.png %})

   There are themes available for your mouse pointer that you can temporarily use while recording that can make it easier for watchers to see what you're doing.

   - [Windows](https://www.microsoft.com/en-us/p/mouse-pointer-highlight/9p7sb9s4rq7z?activetab=pivot:overviewtab)
   - [Linux](https://askubuntu.com/questions/777896/how-do-i-highlight-my-mouse-pointer-while-screen-recording/917587#917587)


