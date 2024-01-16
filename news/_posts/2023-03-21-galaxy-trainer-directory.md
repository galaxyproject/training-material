---
title: "New Feature: Trainer Directory! (Add yourself today!)"
tags: [new feature, community building, capacity building]
contributions:
    authorship: [hexylena, lldelisle]
    funding: [gallantries]
layout: news
---

Numerous other training groups keep "Trainer Directories" as a way to help students meet trainers in their area, especially for hosting workshops for their local community. There have been a handful of attempts at this in the past, both in [galaxyproject/galaxy-maps](https://github.com/galaxyproject/galaxy-maps) and in [a Google Map](https://www.google.com/maps/d/u/0/viewer?mid=1r8LJy6la-JeIrg23aZjpwVNjJDE) on the GTN homepage.

We will now store this data within the GTN itself, and trainers can manage this data themselves. If you are a trainer that offers training courses in your local area, or are interested in helping others teach in that region, please feel free to update your Contributors profile with the new metadata:

```yaml
hexylena:
    name: Helena Rasche
    orcid: 0000-0001-9760-8992
    ...
    contact_for_training: true
    location:
      country: NL
      lat: 51.91
      lon: 4.46
```

And we'll add a nice map to the GTN, and potentially the Galaxy Hub as well! The `contact_for_training` field is a true or false value which lets you indicate you're open to being contacted regarding giving training. Other trainers in your region might use this if they're looking for help giving a course or want experts in a field for a course they plan to host. We do not have plans currently expose this to "end users" looking for training, we will probably attempt to push them first to the existing events list.

[Check out the map]({% link hall-of-fame.md %}){: .btn .btn-primary}
