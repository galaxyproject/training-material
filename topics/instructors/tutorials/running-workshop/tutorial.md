---
layout: tutorial_hands_on

title: Running a workshop as instructor
questions:
-
objectives:
-
time_estimation: 3H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- bebatut
- bgruening
- shiltemann
- erasche
---


# Introduction
{:.no_toc}

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Tips for instructors

## Using the Training Material

1. Decide (individually or together) if you will:

     1. Run through the tutorial with them step by step, waiting for everyone at each step
     2. OR allow trainees to follow along at their own pace

   Each of these options has pros and cons:

   - Step by step

      - Pro: everyone is on the same step, helpers can see similar errors at similar times that they can point out to the group at once
      - Con: Some trainees will be faster, and some will be slower. If you wait for everyone it may cause you to go at a slower pace than expected, depending on their skill levels

   - At their own pace:

      - Pro: Those who want to explore features have time to do so. If they are intrigued by the workflow editor, they can play around with it for some time before moving on.
      - Con: People will finish at wildly different times. Be prepared for some users to be done and maybe encourage them with a list of things they can do after they've completed the tutorial, or additional reading material, or just encourage them to try different parameters.

2. Will you present slides? Make sure these are prepared ahead of time. The Training Material offers introductory slides for many topics and many tutorials. See if you can contribute your slides back, if the topic does not currently have any.

## (Un)expected situations

### Students asking extremely specific questions

Occasionally students will either want to show off their knowledge on specific topics, or will just want special attention with their particular research question.
If you can answer the question briefly then this can be good. Otherwise just ask them to bring it up with you during the next break.

### Interruptions

There are many kinds of interruptions that can occur:

Fire alarm / evacuate the classroom
:    There is not much you can do in this situation. If they bring their stuff with them (they shouldn't if it's a fire alarm), then consider continuing outside. However this can be unpleasant for participants due to weather, or too much sun on their laptops. Make a break, point out nearby attractions, and make a plan to meat back in 30 minutes.

Internet / Server Outage
:    Can you estimate how long this will last? If it is going to be a short interruption, this can be a good time to recap briefly and then ask questions of the students, or suggest students should ask their questions. If it will be a long outage, you can switch to your backup server, OR you can walk the students through the rest of the steps. Be sure to go slowly and be detailed with what would have happened during the steps that they will miss.

(Any others?)

### Unsocial participants

Not everyone is always outwardly enthusiastic about the course for any number of reasons; some students may be very introverted, maybe haven't had their coffee yet, were told to be there by their boss, or just not willing to be social. While you cannot force them everyone, you can encourage try and make activities cover a wide range of groups, from small to large. Some people will be more willing to join smaller groups, or will be more than willing to interact with you one-on-one. Ice breakers that involve one-on-one conversations can be good, a nice one is to have a 'bingo' chart with various aspects about researchers (has PhD, works with viruses, is from local city, is from outside of the country, etc.) and ask them to find people matching these attributes. This allows them to have more personal conversations.

# Checklists

## Before the workshop

1. Read the [teaching recommendations](https://carpentries.github.io/instructor-training/22-practices/) of Software Carpentry
2. Read the [organiser recommendations]({{ site.baseurl }}/topics/instructors/tutorials/organize-workshop/tutorial.html) and see if there is anything you can help with
2. Decide on the order of the lessons and who will teach what with the organizers
3. Assist the host in recruiting helpers (people to walk around the classroom and help people if they get stuck)
4. Assist the host in ensuring the workshop location is accessible and everything is ready for the students
5. Share emergency contact information with host in case of last minute changes
6. Practice teaching the material
7. Remind the host of any special equipment you might need (e.g. Mac connector)
8. Confirm criteria for reimbursement (per diem or save receipts)
9. Create a recap with the different formats and workflows taught during the workshop and print them

### Prepare your Galaxy

1. Identify where you will run the training:

   1. A Public Galaxy?

      - UseGalaxy.org? UseGalaxy.eu? .org.au? (etc.)
      - If you run on UseGalaxy.eu, they provide free [Training Infrastructure]({{ site.baseurl }}/topics/instructors/tutorials/setup-tiaas-for-training/tutorial.html)
      - Work with the Galaxy administrator to ensure that the tools and datasets required are available

   2. A Private Galaxy?

      Check our tutorial on [how to set up a Galaxy for Training]({{ site.baseurl}}/topics/instructors/tutorials/setup-galaxy-for-training/tutorial.html)

      - Technical infrastructure
      - Install the needed tools and data
      - Check the scalability for the number of participants
      - Test the tutorials on the infrastructure

   3. **No matter what**: Have a backup plan!

2. Run the tutorials several times on your server and your backup server

3. Identify the time / memory limiting steps and plan accordingly

## Before the participants arrive

1. Check accessibility

    - Set up the computer and go to the back of the room. Is the font used in your slides big enough? Is the font for showing the training materials?
    - If your room has a mic, use it, and have a helper stand in the back of the room to confirm that you can be heard

## During the workshop

1. Review Code of Conduct with participants
2. Support host in collecting attendee names and emails
3. Remind learners to use sticky notes to give feedback

   Ask "Who is not yet finished?" and remind them to place the post-it notes so you can see quickly whether or not you can continue (assuming you have chosen step-by-step instructionm). If you are using the at-your-own-pace instruction, you can ask them to place postit notes if they have questions (better than their arms getting tired.)

4. Make regular breaks (every 1-2 hours for 15-30 minutes)
5. Engage and interact with the participants

    You can use gamification and interactive examples to help participants learn. For example: use LEGO bricks to demonstrate how mapping works (TODO: Link to how to do this demo.)

    Find someone you haven't talked to yet today, and talk to them. Make sure all the helpers have the same goal.

6. Regularly ask questions of the trainees:

  - What data formats work here?
  - Where is the data coming from?
  - What data transformation is being made?
  - What were the data formats used today?
  - What were the main steps of today's workflow?

7. Regularly remind trainees of things you've taught them

   Repeat again and again, in different ways, to reinforce concepts

   Explain and re-explain file formats and group them for their uses (e.g. SAM/BAM/CRAM are sequence alignment formats, fasta is for genomic sequences, FastQ variants all store sequencing reads)

8. Draw the workflow, maybe interactively ask students to remind you of next steps
9. Do recap at the end of each day about what happened during the day
10. Do recap each morning about what happened yesterday
11. Provide the students with a recap of the formats and workflows as a take-home

## After the workshop

1. Submit receipts for reimbursement
2. Report any training errors or suggest improvements on github, gitter or by email to the GTN
3. Update this document with anything that you did that wasn't listed.

# Conclusion
{:.no_toc}
