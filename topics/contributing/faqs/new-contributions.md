---
title: Using the new Contributions Annotation framework
area: gtn
box_type: tip
layout: faq
contributors: [hexylena]
---

If you are writing a tutorial or slides, there are two ways to annotate contributions:

The old way, which doesn't accurately track roles

```
contributors: [hexylena, shiltemann]
```

And the new way which lets you annotate who has helped build a tutorial in a much richer way:

```yaml
contributions:
    authorship:
        - shiltemann
        - bebatut
    editing:
        - hexylena
        - bebatut
        - natefoo
    testing:
        - bebatut
    infrastructure:
        - natefoo
    translation:
        - shiltemann
    funding:
        - gallantries
```

This is especially important if you want to track funding or infrastructure contributions. The old way doesn't allow for this, and thus we would *strongly* recommend you use the new format!
