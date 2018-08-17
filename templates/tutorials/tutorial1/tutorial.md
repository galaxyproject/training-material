---
layout: tutorial_hands_on
topic_name: templates
tutorial_name: tutorial1
---

# Introduction
{:.no_toc}

<!-- This is a comment. -->

General introduction about the topic and then an introduction of the tutorial (the questions and the objectives). It is nice also to have a scheme to sum up the pipeline used during the tutorial. The idea is to give to trainees insight into the content of the tutorial and the (theoretical and technical) key concepts they will learn.

**Please follow our [tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Part 1

Introduction about this part

## Subpart 1

Short introduction about this subpart.

<!--
{% icon hands_on %} will render the hands_on icon as specified in
_config.yml in the root of this repository.
-->

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Step1
> 2. Step2
>
>    > ### {% icon comment %} Comments
>    > A comment
>    {: .comment}
>
>    > ### {% icon tip %} Tip: A tip
>    >
>    > * Step1
>    > * Step2
>    {: .tip}
{: .hands_on}

## Subpart 2

Short introduction about this subpart.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Step1
> 2. Step2
>
>    > ### {% icon question %} Question
>    >
>    > Question?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > Answer to question
>    > >
>    > {: .solution}
>    >
>    {: .question}
{: .hands_on}

Some blabla
> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Step1
> 2. **My Tool** {% icon tool %} with the following parameters
>   - *"param1"*: the file `myfile`
>   - *"param2"*: `42`
>   - *"param3"*: `Yes`
>
> 3. **My Tool** {% icon tool %} with the following parameters
>   - {% icon param-text %} *"My text parameter"*: `my value`
>   - {% icon param-file %} *"My input file"*: `my file`
>   - {% icon param-files %} *"My multiple file input or collection"*: `my collection`
>   - {% icon param-select %} *"My select menu"*: `my choice`
>   - {% icon param-check %} *"My check box"*: `yes`
>
>    > ### {% icon question %} Questions
>    >
>    > 1. Question1?
>    > 2. Question2?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > 1. Answer for question1
>    > > 2. Answer for question2
>    > >
>    > {: .solution}
>    >
>    {: .question}
>
> 3. Step3
{: .hands_on}

> ### {% icon warning %} Warning: Be careful about ...
>
> Add more details in Markdown.
>
{: .warning-box}

# Part 2

Short introduction about this subpart.

> ### {% icon comment %} Comment
>
> Do you want to learn more about the principles behind mapping? Follow our [training](../../NGS-mapping)
{: .comment}



> ### {% icon details %} Background: More details on the ....
>
> Add more details in Markdown. By default the box is collapsed. And is expanded when clicked
>
{: .details}


# Conclusion
{:.no_toc}

Conclusion about the technical key points. And then relation between the techniques and the biological question to end with a global view.
