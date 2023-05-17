---
title: Pick the right Concatenate tool
description: Most Galaxy servers will have two Concatenate tools installed - know which one to pick!
area: analysis
box_type: tip
layout: faq
contributors: [wm75]
---

On most Galaxy servers you will find two {% icon tool %} **Concatenate datasets** tools installed:

1. **Concatenate datasets** tail-to-head
2. **Concatenate datasets** tail-to-head (cat)

{% if include.toolspec %}
{% icon tip %} The tool you want to use at this step is `{{ include.toolspec }}`, but if you are interested in the details, here is why:
{% endif %}

The two tools have nearly identical interfaces, but behave differently in certain situations, specifically:
- The second tool, the one with "(cat)" in its name, simply concatenates everything you give to it into a single output dataset.

  Whether you give it multiple datasets or a collection as the first parameter, or some datasets as the first and some others as the second parameter, it will always concatenate them all.
  In fact, the only reason for having multiple parameters for this tool is that by providing inputs through multiple parameters, you can make sure they are concatenated in the order you pass them in.
- The first tool, on the other hand, will only ever concatenate inputs provided through *different* parameters.

  This tool allows you to specify an arbitrary number of {% icon param-file %} single datasets, but if you also want to use {% icon param-files %} multiple datasets or {% icon param-collection %} a collection for some of the Dataset parameters, then all of these need to be of the same type (multiple datasets *or* collections) and have the same number of inputs.

  Now depending on the inputs, one of the following behaviors will occur:
  - If all the different inputs are {% icon param-file %} single datasets, the tool will concatenate them all and produce a single output dataset.
  - If all the different inputs are specified *either* as {% icon param-files %} multiple datasets *or* as {% icon param-collection %}, and all have the same number of datasets, then the tool will concatenate the first datasets of each input parameter, the second datasets of each input parameter, the third, etc., and produce an output collection with as many elements as there are inputs per Dataset parameter.
  - In extension of the above, if *some* additional inputs are provided as {% icon param-file %} single datasets, the content of these will be recycled and be reused in the concatenation of all the nth elements of the other parameters.
