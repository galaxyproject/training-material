---
layout: tutorial_hands_on

title: "Galaxy Development"
questions:
  - "How do I debug Galaxy?"
objectives:
  - "Fix a broken branch"
  - "Interpret the results of failed tests on GitHub"
  - "Run individual tests locally"
  - "Fix simple errors identified by failing tests"
  - "Write a simple test exposing the identified bug"
time_estimation: "2h"
key_points:
  - "Unit Tests"
  - "API or Integration Tests"
  - "Client Linting Failure"
  - "Client Tests"
  - "Selenium Tests"
  - "Runtime Error"
requirements:
  -
    type: "internal"
    topic_name: dev
    tutorials:
      - architecture
  - Familiarity with basic Git commands
  - Basic knowledge of Python and Javascript
  - Mac OS or Linux that can run Galaxy & your favorite IDE or editor
contributors:
  - assuntad23
  - ic4f
  - jmchilton
---

## Introduction

In this tutorial we are going to demonstrate how to find common types of bugs or errors that you may encounter as a contributor to Galaxy.We aim to step you through the process of finding and fixing a bug - from identifying them on GitHub, to finding and developing a solution, to committing your code change.

With the skills from this tutorial, it is our hope that you will feel more prepared to develop solutions for Galaxy and more confidently navigate the obstacles along the way. 


## Some Basics

For almost every type of error you encounter, there will be a relatively similar process:

1. View failing test/error on GitHub workflows
2. Run failing test locally
3. Search for issue causing failure
4. Analyze the code, develop and implement a solution
5. Ensure test passes locally
6. Push fixed branch to local fork

Note: This process will not be the same for the runtime error, since there will be no failed test to review initially.

## Getting Started: Pulling the Buggy Fork

Good news! All of the failing tests are on the same fork of Galaxy! That means you only need to download and work off of one branch. The downside is that this is one buggy branch with issues all over the place. 

We have no doubts that you can solve them though!

The fork is located at : TODO [insert link] 


### Directions for Cloning/Basing off this fork

TODO

> ###  {% icon comment %} Galaxy Repository
> Note: This tutorial uses a fork to house the issues created for the purposes of this tutorial. In the future, when you contribute to the Galaxy project, you'll want to use one of the repositories located at [https://github.com/galaxyproject]. Most likely you'll be using the main repo, simply titled "galaxy".
{: .comment}

## Unit Test Failure

We're going to start with our failing unit test. 
TODO definition/purpose of tests? 

### Finding the failing test on GitHub

If you navigate to TODO [insert link] you'll see that there are quite a few failing checks. One of them says TODO (check text >)"Unit Test". Feel free to click on that to see a log of all the tests. 

Try to identify which test is failing based on the logs.

TODO hidden solution test/unit/util/test_utils.py::test_strip_control_characters

### Running the test locally

TODO

### Finding and Fixing the Issue

TODO

If you think you've got it, try running the test locally again to be sure. 

### A good commit message

TODO

Congratulations - you found and fixed a Unit Test failure! 

## API or Integration Test Failure

Next up is the API test. 
TODO definition/purpose of tests?  

### Finding the failing test on GitHub

One of the failing tests on GitHub says TODO (check text >)"API Test". Feel free to click on that to see a log of all the tests. 

Try to identify which test is failing based on the logs.

TODO hidden solution 

### Running the test locally

TODO

### Finding and Fixing the Issue

TODO

If you think you've got it, try running the test locally again to be sure. 

### A good commit message

TODO

Woo-hoo - you fixed a Failing API Test! 

## Client Linting Failure

Client Linting is next on deck. 
TODO definition/purpose of tests?  

### Finding the failing test on GitHub

One of the failing tests on GitHub says TODO (check text >) "Client Linting". Feel free to click on that to see a log of all the tests. 

Try to identify which test is failing based on the logs.

TODO hidden solution 

### Running the test locally

TODO

### Finding and Fixing the Issue

TODO

If you think you've got it, try running the test locally again to be sure. 

### A good commit message

TODO

Awesome work - you now understand how to solve client linting issues!

## Client Test Failure

A failing Client Test walks into a bar... 
TODO definition/purpose of tests?  

### Finding the failing test on GitHub

One of the failing tests on GitHub says TODO (check text >) "Client Test". Feel free to click on that to see a log of all the tests. 

Try to identify which test is failing based on the logs.

TODO hidden solution 

### Running the test locally

TODO

### Finding and Fixing the Issue

TODO

If you think you've got it, try running the test locally again to be sure. 

### A good commit message

TODO

Nice Job - you now have a handle on client tests!

## Selenium Test Failure

Selenium tests can be difficult to fix, since they cover an expanse of code, but I'm sure you can do it.
TODO definition/purpose of tests?  

### Finding the failing test on GitHub

One of the failing tests on GitHub says TODO (check text >) "Selenium Test (3.1, 3.0? --- what does this even mean?)". Feel free to click on that to see a log of all the tests. 

Try to identify which test is failing based on the logs.

TODO hidden solution 

### Running the test locally

TODO

### Finding and Fixing the Issue

TODO

If you think you've got it, try running the test locally again to be sure. 

### A good commit message

TODO

Excellent! Selenium tests are in the bag!

## Runtime Error 

TODO

## Finally --- Back in the Green



## Conclusion

First of all, thank you for completing this tutorial. We hope you feel more confident debugging your code on Galaxy. Of course, we are always available to help answer questions and support you!

