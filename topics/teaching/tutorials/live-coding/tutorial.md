---
layout: tutorial_hands_on

title: Live Coding is a Skill
subtopic: practises
enable: false
time_estimation: 1h
questions:
  - Why do we teach programming using participatory live coding?
objectives:
  - Explain the advantages and limitations of participatory live coding.
  - Summarize the key dos and do nots of participatory live coding.
  - Demonstrate participatory live coding.
key_points:
  - Live coding forces the instructor to slow down.
  - Coding-along gives learners continuous practice and feedback.\
  - Mistakes made during participatory live coding are valuable learning opportunities.
contributors:
  - bebatut
  - fpsom
  - carpentries
---

# Introduction


---

One of the most effective ways to teach new computational skills is through live coding {% cite Rubi2013 %} {% cite Haar2017 %} {% cite Raj2018 %}: *instructors do not use slides to teach the material*, but work through the lesson material, typing in the code or going through the instructions, with the workshop participants following along. This tutorial explains how it works, why we use it, and gives general tips for an effective participatory live coding presentation. We will finish this tutorial by practicing ourselves and providing feedback for each other.


## Why Participatory Live Coding?

We do not use slides in our lessons. Instead, instructors plug their laptop into the projector and work through the lesson, typing in the code, reformatting data, and talking as we go. This is called ["live coding"](https://en.wikipedia.org/wiki/Live_coding).

However, the instructor is not live coding in a vacuum. Importantly, learners are strongly encouraged to "code-along" with the instructor. We refer to the practice of having the instructor live code and the learners code along as "participatory live coding" or, less formally, 'code-along sessions'.

> <hands-on-title>Up and Down (5 min)</hands-on-title>
>
> List some advantages and challenges of participatory live coding from both a learner's and an instructor's point of view.
>
> > <solution-title></solution-title>
> > Some advantages are:
> > * Watching a program being written is more compelling than watching someone page through slides that present bits and pieces of the same code.
> > * It enables instructors to be more responsive to "what if?" questions. Where a slide deck is like a railway track, participatory live coding allows instructors to go off-road and follow their learners' interests.
> > * Lateral knowledge transfer: participatory live coding facilitates the transfer of [tacit knowledge](https://jonudell.net/udell/2006-09-19-screencasting-of-tacit-knowledge.html) -- people learn more than we realized we were teaching by watching *how* instructors do things.
> > * It slows the instructor down: if she has to type in the program as she goes along, she can only go twice as fast as her learners, rather than ten-fold faster as she could with slides.
> > * Learners get to see instructors' mistakes *and how to diagnose and correct them*. Novices are going to spend most of their time doing this, but it is left out of most textbooks.
> >
> > Some challenges are:
> > * It requires instructors to be able to improvise when things go wrong or when learners have questions not directly addressed in the text of the lesson.
> > It can be hard for learners to listen and type at the same time, due to the *split-attention effect* that was discussed in the [learning principles tutorial]({% link topics/contributing/tutorials/learning-principles/tutorial.md %}). This is why it is very important that instructors first explain what they are going to do, then say what they are typing as they type it, and then explain what they did again afterwards.
> > It may take a bit of practice for instructors to get used to thinking aloud while coding in front of an audience.
> >
> {: .solution}
{: .hands_on}

Live coding fits well into the practice-feedback model we have been discussing - by providing learners with continuous opportunities for practice (every time they type in a line of code) and continuous feedback (their code either works or fails with an error message). It is important to keep in mind, however, that feedback is not helpful if you cannot understand it.
Many error messages are obscure and not written with novices in mind. Continue to use the strategies for error framing that
[was discussed in the motivation tutorial]({% link topics/teaching/tutorials/learner_participation_engagement/tutorial.md %}) to make sure this feedback is useful to learners.

> <hands-on-title>Compare and Contrast (15 min)</hands-on-title>
>
> Watch this first participatory live coding demo video: [https://youtu.be/bXxBeNkKmJE](https://youtu.be/bXxBeNkKmJE) and this second demo video: [https://youtu.be/SkPmwe_WjeY](https://youtu.be/SkPmwe_WjeY) as a group and then summarize your feedback on both. Use the 2x2 rubric for feedback that was [discussed in the feedback tutorial]({% link topics/teaching/tutorials/assessment/tutorial.md %}).
>
> In the videos, the bash shell `for` loop is taught, and it is assumed learners are familiar with how to use a variable, the `head` command and the content of the `basilisk.dat unicorn.dat` files.
>
> Note: Sometime sounds in the room can be poor. Turning on closed captioning by pressing the cc button will improve the accessibility of these videos.List some advantages and challenges of participatory live coding from both a learner's and an instructor's point of view.
>
> > <solution-title></solution-title>
> > The main idea is for the Instructor to lead a discussion about the videos and your feedback on them, making sure that the points of the Top Ten Tips below have been made.
> > If this a self-taught session, please do reflect on your feedback by contrasting against the Top Ten Tips below.
> {: .solution}
>
{: .hands_on}

## Top Ten Tips for Participatory Live Coding in a Workshop
1. **Stand up and move around the room if possible.** This makes the experience more interactive and less monotonous. Use a microphone if one is available to make it easier for people with hearing difficulties to hear you.
2. **Go slowly.** For every command you type, every word of code you write, every menu item or website button you click, say out loud what you are doing while you do it.  Then point to the command and its output on the screen and go through it a second time.  This slows you down and allows learners to copy what you do, or to catch up.  Do not copy-paste code.
3. **Mirror your learner's environment.** Try to create an environment that is as similar as possible to what your learners have to reduce cognitive load. Avoid using keyboard shortcuts.
4. **Use your screen wisely.** Use a big font, and maximize the window.  A black font on a white background works better than a light font on a dark background.  When the bottom of the projector screen is at the same height, or below, the heads of the learners, people in the back will not be able to see the lower parts.  Draw up the bottom of your window(s) to compensate. Pay attention to the lighting (not too dark, no lights directly on/above the presenter's screen) and if needed, re-position the tables so all learners can see the screen, and helpers can easily reach all learners.
5. **Use illustrations** to help learners understand and organize the material. You can also generate the illustrations on the board as you progress through the material.  This allows you to build up diagrams, making them increasingly complex in parallel with the material you are teaching.  It helps learners understand the material, makes for a more lively workshop and gathers the learners' attention to you as well.
6. **Turn off notifications** on your laptop and phone.
7. **Stick to the lesson material.** The core Carpentries lessons are developed collaboratively by many instructors and tried and tested at many workshops.  This means they are very streamlined - which is great when you start teaching them for the first time.  It may be tempting to deviate from the material because you would like to show a neat trick, or demonstrate some alternative way of doing something.  Do not do this, since there is a fair chance you will run into something unexpected that you then have to explain.  If you really want to use something outside of the material, try it out thoroughly before the workshop: run through the lesson as you would during the actual teaching and test the effect of your modification.
Some instructors use printouts of the lesson material during teaching. Others use a second device (tablet or laptop) when teaching, on which they can view their notes and the Etherpad session.  This seems to be more reliable than displaying one virtual desktop while flipping back and forth to another.
8. **Leave no learner behind.** Use sticky notes, see below, to gauge learners' progress and understanding.
9. **Embrace mistakes.** No matter how well prepared you are, you will make mistakes. This is OK! Use these opportunities to do error framing ([as discussed in the motivation tutorial]({% link topics/teaching/tutorials/learner_participation_engagement/tutorial.md %})) and to help your learners learn the art of troubleshooting.
10. **Have fun!** It is OK to use humor and improvisation to liven up the workshop. This becomes easier when you are more familiar with the material, and more relaxed. Start small, even just saying 'that was fun' after something worked well is a good start.

Read more in [Ten quick tips for teaching with participatory live-coding] {% cite Nederbragt2020 %}, as well as the excellent section on "Live Coding" from the book "Teaching Tech Together" by Greg Wilson {% cite Wilson2019 %}.


> <hands-on-title>Practice Teaching (25 min)</hands-on-title>
>
> 1. Split into groups of three.
> 1. Assign roles, which will rotate: presenter, timekeeper, note-taker.
> 2. Have each group member teach 3 minutes of your chosen lesson episode using live coding. For this exercise, your peers will not "code-along." Before you begin, briefly describe what you will be teaching and what has been learned previously. Do not record this exercise.
> 3. After each person finishes, each group member should share feedback (starting with themselves) using the same 2x2 rubric as yesterday. The timekeeper should keep feedback discussion to about 1 minute per person; this may leave some time at the end for general discussion. The note-taker should record feedback in the Etherpad.
> 4. Trade off roles.
>
> > <solution-title></solution-title>
> > If this a self-taught session, one way to approach this exercise is by recording your teaching demonstration (e.g. through a cell-phone), and attempting to provide an objective self-feedback based on the 2x2 rubric.
> {: .solution}
{: .hands_on}




# Recap

##
