---
title: "New Feature: Automatic Slides-to-video conversion"
contributors: [hexylena,delphine-l]
tags: [new feature, videos, gtn, pandemic, remote-teaching]
layout: news
---

This past year has been a struggle for all of us, but one task has become apparent as one of the most time consuming portions of course preparation during a pandemic: **Video Production**.

We have found ourselves spending countless hours recording video lectures, producing videos of us giving slide decks, and videos of us giving trainings on how to run tools in Galaxy. Worse, whenever changes were made we would have to re-record lessons!

Inspired by the efforts of {% include _includes/contributor-badge.html id="delphine-l" %} in her [video-lectures](https://github.com/galaxyproject/video-lectures/) series, {% include _includes/contributor-badge.html id="hexylena" %} has implemented a similar feature of automatic Text-to-Speech (TTS) in the [Galaxy Training Network]({% link topics/contributing/tutorials/slides-with-video/tutorial.md %}).

The GTN had existing infrastructure for producing slide decks from easy-to-maintain markdown documents, and, key to the success, existing syntax for writing "Speaker Notes" on slides. These were intended to be little reminders to instructors giving slide decks about what they should say on each slide.

Thanks to this existing groundwork, after we implemented <abbr title="Text to Speech">TTS</abbr>, contributors only had to add `video: true` to the top of their slides, and they would opt-in to automatic video production. Now our contributors could easily produce slide deck recordings with very little overhead. Best of all, anytime a contributor updated their slide decks, the videos could be re-built as well. As we know the text that would be spoken we could also trivially add accurate subtitles, ensuring accessibility for our deaf/Deaf/<abbr title="Hard of Hearing">HoH</abbr> community members.

All of this process is completely automated, from building the videos, to uploading them to our S3 bucket and including them in the GTN. Since implementation in November, we have added 19 videos to the GTN! See the full list in the [GTN video page]({% link videos/index.md %})

<video controls="" preload="metadata" width="600" height="447" aria-label="a video produced by the GTN text-to-slides implementation. It is very inaccessible for blind and vision-impaired users, we recommend just reading the slides and speaker notes. Images there have alt-text, and speaker notes are 100% of the spoken content.">
	<source src="https://galaxy-training.s3.amazonaws.com/videos/topics/introduction/tutorials/galaxy-intro-short/slides.mp4" type="video/mp4">
	<track label="English" kind="captions" srclang="en" src="/assets/slides.en.vtt" default="">
</video>
