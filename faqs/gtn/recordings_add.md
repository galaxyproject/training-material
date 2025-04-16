---
title: Adding your recording to a tutorial or slide deck
area: contributors
layout: faq
box_type: tip
contributors: [shiltemann]
---

**We welcome anybody to submit their recordings!** Your videos can be used in (online) training events, or for self-study by learners on the GTN.

For some tips and tricks about recording the video itself, please ensure your recording conforms to our recommendations:

[Recording Tips & Tricks]({% link faqs/gtn/recordings_create.md %}){: .btn.btn-info}
[Submit a Recording](https://forms.gle/qNG8FkTN1yRZPNZY6){: .btn.btn-info}

#### Submission process

The process of adding recordings to the GTN is as follows:

1. **Instructor:** Record video ([tips & tricks]({% link faqs/gtn/recordings_create.md %}))
2. **Instructor:** Submit your video using this [Google Form](https://forms.gle/qNG8FkTN1yRZPNZY6)
3. **GTN:** A [GTN GitHub pull request (PR)](https://github.com/galaxyproject/training-material/pulls) will be made by our bot based on the form.
4. **GTN:**: We will upload your video to the [GalaxyProject YouTube channel](https://www.youtube.com/c/galaxyproject)
5. **GTN:**: We will put the auto-generated captions from YouTube into a Google Doc
6. **Instructor:**: Check and fix the auto-generated captions
7. **GTN:** Upload the fixed captions to YouTube
8. **GTN:** Merge the Pull Request on GitHub
9. Done! Your recording will now show up on the tutorial for anybody to use and re-use


**Note:** If you are submitting a video to use in an event, please submit your recording **2 weeks before the start of your course** to allow ample time to complete the submission process.

#### Recordings Metadata

Our bot will add some metadata about your recording to the tutorial or slide deck in question, and looks as follows:

```
recordings:
  - speakers:       # speakers must be defined in the CONTRIBUTORS.yaml file
    - shiltemann
    - hexylena
    captioners:     # captioners must also be present in the CONTRIBUTORS.yaml file
    - bebatut
    type:           # optional, will default to Tutorial or Lecture, but if you do something different, set it here (e.g. Demo, Lecture & Tutorial, Background, Webinar)
    date: '2024-06-12'         # date on which you recorded the video
    galaxy_version: '24.0'     # version of Galaxy you used during the recording, can be found under 'Help->About' in Galaxy
    length: 1H17M              # length of your video, in format: 17M or 2H34M  etc
    youtube_id: "dQw4w9WgXcQ"  # the bit of the YouTube URL after youtube.com/watch?v=

  - speakers:
    - shiltemann
    captioners:
    - hexylena
    - bebatut
    date: '2020-06-12'
    galaxy_version: '20.05'
    length: 51M
    youtube_id: "oAVjF_7ensg"

```

#### Misc

**Note:** If your videos are already uploaded to YouTube, for example as part of a different project's account, you can add this metadata to the tutorial or slides manually, without using our submission form.
Note that we do require all videos to have good-quality English captions, and we will not be able to help you configure these on other YouTube accounts.
