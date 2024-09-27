---
title: How to Wrap Your Tool
description: There are many different ways to wrap a tool and some of the choices you'll want to make change depending on what you have available.
area: tool-dev
box_type: tip
layout: faq
contributors: [hexylena]
---

## Dependencies

{% include _includes/cyoa-choices.html title="How are your dependencies currently available?" disambiguation="1" option1="Conda" option2="Docker" option3="BioContainers" option4="None of the above" default="__hidden__" %}


<div class="Conda" markdown=1>
Great! You're good to go then! Galaxy natively uses Conda based dependencies.

in your tool XML you'll just need to specify packages and their versions, based on your conda environment.

```xml
<requirements>
  <requirement type="package" version="1.0.0">my-package</requirement>
</requirements>
```

</div>

<div class="Docker" markdown=1>

Great! You're good to go then! Galaxy can natively use Docker based dependencies. In your tool XML you'll just need to specify the container and it's version

```xml
<requirements>
  <container type="docker">quay.io/galaxy/homer-galaxy:4.11</container>
</requirements>
```

TODO: how to inputs/outputs work? container commands?

</div>

<div class="BioContainers" markdown=1>

If your dependency is a [BioContainers](https://biocontainers.pro/) container, there's a chance it's available in Conda as well, as biocontainers are generated from conda packages automatically. You can check the conda repositories for the upstream package, and if it's available, you can treat it like a conda package:

```xml
<requirements>
  <requirement type="package" version="1.0.0">my-package</requirement>
</requirements>
```

</div>
