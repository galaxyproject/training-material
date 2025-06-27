---
layout: tutorial_hands_on

title: An introduction to Jupyter-based teaching
subtopic: introduction
draft: true
time_estimation: 15m
questions:
  - What are the general advantages and limitations of Jupyter-based teaching?
objectives:
  - An overview is provided to familiarize teachers with the advantages and limitations of Jupyter-based teaching.
key_points:
  - Jupyter-based teaching offers an easy entrypoint for students, but deprives them of the need to learn basic setup steps themselves.
contributors:
  - mittler-works
---

Jupyter Notebooks are a promising way to provide participants of a course or workshop with an easy-to-use, webbased, interactive and reproducible environment to learn and work with.

Jupyter Notebooks makes it possible to Code in various programming languages, writing documentation and other text, doing math, printing plots and much more in one single environment. All relevant resources are ready to be used.

Jupyter Notebooks can be run on an individual machine as well as centrally managed by a JupyterHub server. If provided via JupyterHub, the webbased access truly allows participants to bring their own device, as their devices only need to satisfy common web standards and do not need to satisfy the requirements of the software running inside the notebook, as this will be run on the server. No installation on the participants device required.

This allows for an low-level access especially for beginner courses or workshops.

As a teacher you may provide a predefined learning environment, i.e. a predefined software stack, data, configuration and more. This allows you to cover all the necessary requirements in advance, so your participants do not have to carry out further setup steps during the course, as these tend to be frustrating in many szenarios.

If you use an infrastructure provider for your course, you only have to deal with your JupyterLab customizations in order to adapt JupyterLab to your course. The deployment of these JupyterLabs in order to provide them to your participants is the infrastructure providers' job. With this in mind it makes little difference if you plan your JupyterLab for three or three hundred participants, as it is horizontally scalable.

# Jupyter Terminology

To get started, we have to elaborate on some terminology. Especially the term `Jupyter Notebook` is heavily overloaded and may be used in many szanarios while having very different meanings. Over all, in the context of Jupyter, there are some components with similar terminology involved. It is important to know how to keep them apart.

> <question-title>Which Jupyter-related terms do you now?</question-title>
> 
> Just take a moment and think about what Jupyter-related terms you have already heard about. For each term that comes to your mind, can you explain how it fits into the Jupyter landscape?
{: .question}


Let's start this tutorial with a list of a few very common and important Jupyter terms.

### Jupyter
Jupyter itself is the name of the umbrella open source project that provides software components like Jupyter Notebook, JupyterLab, JupyterHub and many more projects. Jupyter itself is not a software nor a file format.

### Jupyter Notebook
Jupyter Notebook is a computer program that serves a webbased interactive tool for writing, documenting, and sharing code. It is commonly used on a local machine. Note, this term is heavily overloaded and it is common that the `.ipnb` files themselves are reffered to as Jupyter Notebooks too.

Please always mind the context to distinguish between the Jupyter Notebook program and the Jupyter Notebook file.

### JupyterLab
JupyterLab is a computer program that serves an extensible webbased environment for interactive coding, that combines Jupyter Notebook with other components like text editors, terminals, kernels, plots and more.

### JupyterHub
JupyterHub is a centralized organizational platform for managing JupyterLab instances on a server. It allows to provide many users with JupyterLab instances and handles authentication and authorization.

# Advantages for teaching

* User management

* Integratable with LMS (Moodle etc)

* Interactive teaching

* Code, Data, Documentation, Results --> All in one place

* Customizable for your course needs

* Predefined Environment

* Web-based, BYD --> No onboarding time, instant teaching, Easy to use

* Localization

* Reproducability

# Disadvantages for teaching

* Important first steps are not need to be taught

* Bad Code Practices can easily creep in

* Less debugging possibilities

* Dependecies

* Habituation effects

* Jail effects

# (WIP)
