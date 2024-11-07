---
title: Creating a GTN Event
area: contributors
layout: faq
box_type: tip
contributors: [shiltemann]
---


To add your event to the GTN, you will need to supply your course information (dates, location, program, etc). You will then get [an event page like this]({% link events/2024-04-01-example-event.md %}) which you can use during your training. This page includes a course overview, course handbook (full program with links to tutorials) and setup instructions for participants.

Your event will also be shown on the [GTN event horizon]({% link events/index.md %}) and on the homepage. We are also happy to advertise your event on social media and Matrix channels.


**Already have your own event page?** No problem! You can add your event as and external event (see below) and we will simply link to your page!

To add your event to the GTN:

1. Create a page in the `events/` folder of the [GTN repository](https://github.com/galaxyproject/training-material)
2. Have a look at example event definitions in this folder:
   - [2024-04-01-example-event.md](https://github.com/galaxyproject/training-material/blob/main/events/2024-04-01-example-event.md)
   - or [2024-04-01-example-event-external.md](https://github.com/galaxyproject/training-material/blob/main/events/2024-04-01-example-event-external.md) if you already have an event page elsewhere
3. Adapt one of these example pages to fit your event
4. Create a pull request on the GTN

**We are also happy to help you** to add your event, please [contact us on Matrix](https://matrix.to/#/#Galaxy-Training-Network_Lobby:gitter.im) to discuss the details of your course with us.

For a full list of metadata fields for events, please have a look at our [schema documentation page](https://github.com/galaxyproject/training-material/blob/main/_layouts/event.html)

Please also feel free to contact us with ideas for improvements! We know that training comes in many different forms, so if something in your event is not yet supported, let us know and we are happy to add it!


### External events

Already have a course webpage? Great! In this case, you only have to provide the most basic information about your course (title, desciption, dates, location).

The easiest method is to fill in our Google Form:

[Events Google Form!](https://forms.gle/4KjCKKrZ6kamg81o7){: .btn.btn-success}

Or you can create the event file manually. See also [2024-04-01-example-event-external.md](https://github.com/galaxyproject/training-material/blob/main/events/2024-04-01-example-event-external.md) for an example definition.

```
---
layout: event-external
title: My External Training Event Title

external: "https://galaxyproject.org/events/"
description:

date_start:
date_end:  # optional, for multi-day events

location:
  name:
  city:
  country:

contributions:
  organisers:
    - name1
    - name2
```


