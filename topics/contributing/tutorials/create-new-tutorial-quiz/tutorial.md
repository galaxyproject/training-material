---
layout: tutorial_hands_on

title: "Adding Quizzes to your Tutorial"
questions:
  - "How to make a quiz?"
objectives:
  - "Create a quiz"
time_estimation: "15m"
key_points:
  - "Quizzes are helpful for both self-directed learning, and ensuring that in synchronous classes, students are all following the material"
contributors:
  - hexylena
---

Interactive quizzes can be used either alone, or with a classroom of students, to check student's knowledge.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## How it Works

We developed a small Kahoot-like interface, where a teacher can initiate a quiz, students can join the teacher's session (using PeerJS), and do a collaborative quiz.

## Quiz Format

Each quiz starts a lot like the tutorials!

```yaml
title: SQL Advanced Recap
contributors:
- hexylena
```

Within the questions section we have a list of question

```yaml
questions:
- title: How do you find the number of records in a query result?
  answers:
    - count(name)
    - count(*)
    - sum(id)
    - max(id)
  correct:
    - count(name)
    - count(*)
  timeout: 20
  type: choose-many
```

There are a few different types of questions:

- `choose-1`, one correct answer
- `choose-many`, potentially multiple correct answers
- `poll`, there's no right answer, just ask the students how they're feeling!


You can include images in the main area, if you need some context for a question

```yaml
- title: Which of these joins is NOT valid
  image: /training-material/topics/data-science/images/carpentries-sql/sql-join-structure.svg
  answers:
    - Select * From P as P1 Join P as P2 on P1.id = P2.id
    - SELECT * From P Join Q Join V on P.id = Q.person and Q.taken = V.id
    - SELECT * From S join Q on S.name=V.dated and V.site = Q.quant
    - SELECT * From S Join V on S.name = V.site
  correct: SELECT * From S join Q on S.name=V.dated and V.site = Q.quant
  timeout: 60
  type: choose-1
```

Poll's let you check in with students, how they're feeling, or take their opinion on what they think.

```yaml
- title: How are you feeling?
  answers:
    - Great
    - Horrid
  timeout: 20
  type: poll
  live: true
```

The `live` key there allows polls to show the results "live" and let students change answers while the time runs.

## Folder Structure

You place quizzes in a `quiz` subdirectory of your tutorial

```console
.
├── quiz
│   └── a.yaml
└── tutorial.md
```

## Inserting quizzes throughout a tutorial

If you want to show a quiz at a specific place within a tutorial

```markdown
{% raw %}
{% include _includes/quiz.html id="a.yaml" %}
{% endraw %}
```

{% include _includes/quiz.html id="a.yaml" %}

And then it's available!

## Default Location

By default all quizzes are found at the bottom of the tutorial.
