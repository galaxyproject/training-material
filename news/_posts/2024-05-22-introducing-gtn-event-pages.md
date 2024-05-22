---
title: Introducing GTN Event Pages
layout: news

tags:
- gtn infrastructure
- new feature
- events
contributions:
  authorship:
    - shiltemann
    - hexylena
---

Are you organizing a training event using GTN materials? You can now add your event directly to the GTN!

View all upcoming and past events on the brand-new [GTN Events page]({% link events/index.md %}).

## What you will get

- Full event webpage complete with
  - Course overview
  - Course handbook with your full program including links to all materials participants will need on the day
  - Setup instructions for participants
  - See [this example event]({% link events/2024-04-01-example-event.md %})

- Your event is advertised on the GTN Event page
- We are happy to help announce your event on social media
- Already have a course webpage? No problem! you can still add your event to the GTN and we will simply to your event webpage.


## Adding your event

GTN events are defined in similar was as learning pathways, just with some additional practical information (dates/times, location, registration etc)

Have a look at [this FAQ]({% link faqs/gtn/gtn_events_create.md %}) for detailed instructions for adding an event.


## Adding an external event
Already have a course webpage? No problem! you can still add your event to the GTN and we will simply to your event webpage. In this case, you only have to provide the basic information about your course (title, desciption, dates, location)


```
---
layout: event-external
title: My External Training Event Title

external: "https://galaxyproject.org/events/"
description:

date_start:
data_end

location:
  name:
  city:
  country:

contributions:
  organisers:
    - name1
    - name2
```



