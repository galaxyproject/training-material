---
layout: tutorial_hands_on

title: "Debugging with Galaxy"
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
contributors:
  - assuntad23
  - ic4f
  - jmchilton
requirements:
  - type: "none"
    title: "Familiarity with basic Git commands"
  - type: "none"
    title: "Basic knowledge of Python and Javascript"
  - type: "none"
    title: "Mac OS or Linux that can run Galaxy & your favorite IDE or editor"
  - type: "internal"
    topic_name: dev
    tutorials:
      - architecture
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

The fork is located at : TODO

### Directions for Cloning/Basing off this fork

TODO

> ### {% icon comment %} Galaxy Repository
>
> Note: This tutorial uses a fork to house the issues created for the purposes of this tutorial. In the future, when you contribute to the Galaxy project, you'll want to use one of the repositories located at https://github.com/galaxyproject. Most likely you'll be using the main repo, simply titled "galaxy".
{: .comment}

## Unit Test Failure

TODO

### Finding the failing test on GitHub

TODO

### Running the test locally

TODO

### Finding and Fixing the Issue

TODO

If you think you've got it, try running the test locally again to be sure.

### A good commit message

TODO

Congratulations - you found and fixed a Unit Test failure!

## API or Integration Test Failure

TODO

### Finding the failing test on GitHub

TODO

### Running the test locally

TODO

### Finding and Fixing the Issue

TODO

If you think you've got it, try running the test locally again to be sure.

### A good commit message

TODO

Woo-hoo - you fixed a Failing API Test!

## Client Linting Failure

Client Linting is important to reduce errors, improve quality of the codebase, and enforce common rules. These rules are defined by the eslinter defined with Vue best practices, primarily. On Github, the command `eslint` is run to start the process for linting code on the client side.

> ### {% icon comment %} Python Linting
>
> For backend consistency, we also have a python linter. While we won't be running you through that exercise today, you may encounter a python linting issue in the future. The way to diagnose and treat these are very similar to the way to do it for the front end.
{: .comment}

### Finding the failing test on GitHub

One of the failing tests on GitHub says "get_code_and_test". Clicking on the word "Details" assocated with this failure will open a brief view of the CircleCI Checks, with a list of tests that are run as part of that command. You'll see that the one that is failing is js_lint. And if you click that, you'll be navigated to the CicleCI app, where you can see much more information about the failures.

Here you can see the exact command that ran to lint our codebase; it's `cd client && yarn run eslint`.

We can also see that it has found 1 error, where the error is and a description of the problem.

### Running the linter locally

Alright, so we've seen the failure, but we want to verify it locally now. For this, you'll need to open up a terminal.

The first thing you'll need to do is activate your virtual environment:

> `. /path_to_galaxy_on_your_machine/.venv/bin/activate`

Then you've got two options. You can run the same exact command that CircleCI uses `cd client && yarn run eslint` or `make client-lint`, which does `cd client && yarn run eslint && yarn run prettier-check`. This will also make sure your code conforms to formatting standards.

After running that, you should see the same linting error from CircleCI; you should know where it is and how to fix it too!

### Finding and Fixing the Issue

> ### {% icon solution %} Solution
>
> In the file `client/src/components/RuleBuilder/RuleComponent`, there is an unused import. Simple delete the line that says `import Vue from "vue";` 
{: .solution }

If you think you've got it, try running the test locally again to be sure.

### Finishing up

Whoo-hoo! If you're here, you've found the problem, ran it locally, fixed the problem, ran it locally again without errors, and are ready to commit your change. Write a descriptive commit message like `fixed a client linting issue in RuleComponent` or whatever you like, as long as it's clear what you've done.

Awesome work - you now understand how to solve client linting issues!

## Client Test Failure

Client tests are tests written by developers to test front-end code. In the case of the test we're going to walk through right now, we're looking at a Vue test failure, which means the component that it's testing is definitely a Vue component.

### Finding the failing test on GitHub

One of the failing tests on GitHub says "Client Unit Testing / jest". Clicking on Details beside that failure, will open up a the terminal output from that test. Here you should be able to see what test is failing.

> ### {% icon solution %} Solution
>
> The failing test file is `src/components/RuleBuilder/SavedRulesSelector.test.js`, and the failing test in that file is `SavedRulesSelector › disables history icon if there is no history`. 
{: .solution }

### Running the test locally

The first thing you'll need to do is activate your virtual environment:

> `. /path_to_galaxy_on_your_machine/.venv/bin/activate`

To run the client tests locally, you must be in the client directory; if you're just in the galaxy directory, `cd client` should take you to the right place.

Next, type `yarn run test` or `yarn jest` to run all the client tests.

The output should match what you found on Github.

### Finding and Fixing the Issue

Next up, finding the error. Here are some tips:

- Take a look at the failing test and the imports to find out which Vue component the test correlates to.
- Look at the test to see what it does. Is there an error in the test?
- Look at the Vue component where the failure seems to be taking place. Is the error located there?

The solution is hidden below, but try your hand at fixing it first:

> ### {% icon solution %} Solution
>
> The failure is only happening on the disabled test case, so take a look in the template of `SavedRuleSelector` where the conditions for adding disabled are defined.  
> The way it's worded in the test (and the way one could think about it logically), is that if there is no history, the button should be disabled, but the code says `:class="{ disabled: numOfSavedRules == 1 }"` -- whoops! A very silly off by one error. On that very same line, change the 1 to a 0, and run the test again. 
{: .solution }

If you think you've got it, try running the test locally again to be sure.

> ### {% icon comment %} Running fewer client tests
>
> If you don't want to run the whole suite of client tests, you can add keywords that match the test path/name. For example both `yarn run jest rule` and `yarn jest selector` will work, but the former will run other tests with rule in the path/module name.  
{: .comment}

### Finishing up

Nice Job! If you're here, you've found the problem, ran it locally, fixed the problem, ran it locally again without errors, and are ready to commit your change. Write a descriptive commit message like `fixed a logic error in SavedRuleSelector` or whatever you like, as long as it's clear what you've done.

Awesome work - you now understand how to solve client linting issues!- you now have a handle on client tests!

## Selenium Test Failure

Selenium is the end-to-end or integration testing framework that we use; so these tests often are designed to test both frontend and backend code. This makes fixing them a little difficult, but I'm sure you can do it!

### Finding the failing test on GitHub

There are actually a few failing Selenium tests, but the good new is, they are the result of only one little error. One of the failing tests on GitHub says "Selenium Test 3.7, 0" and "Selenium Test 3.7, 1". When you click on the details of this one, it's not always very clear which test is causing the failure, but if you download the artifacts, it's very clear which tests have failed.

Try to identify which tests are failing based on artifacts you downloaded.

> ### {% icon solution %} Solution
>
> The failing tests are:
>
> - test_build_list_and_show_items
> - test_data_column_input
> - test_upload_pair_specify_extension
{: .solution }

### Running the tests locally

We do not want to run the full suite of Selenium tests here - that would take a lot of time and we already know that most of them passed.

Instead, we can build our commands to only run a subset of the tests. Here's an example: `./run_tests.sh -selenium lib/galaxy_test/selenium/test_sign_out.py`

A command like that will still run all of the tests defined in that python file, so we need to find out where the tests that failed are defined. Try seaching the baseline for our three failing Selenium tests; what files are they in?

> ### {% icon solution %} Solution
>
> The failing tests are located:
>
> - test_collection_builders.py
> - test_workflow_editor.py
> - test_uploads.py
{: .solution }

Now that we have the files, we can run the selenium tests. The command to do this should be very similar to `./run_tests.sh -selenium lib/galaxy_test/selenium/test_sign_out.py`. The only difference should be the python file where the tests are defined. 

Running any of these three tests will show us that the tests are still failing on our local machine. So, you should see a failing test that matches one of the three we identified from the downloaded artifacts. Onwards to fixing them!  

### Finding and Fixing the Issue

Good news - there's an easy way to figure out what's going wrong with a selenium test! In your terminal type `export GALAXY_TEST_SELENIUM_HEADLESS=false`. This will allow you to see the test running in real time. Re-run the selenium test to watch the headless browser step through the test. Looking at the test run in the headless browser, alongside the definition for the test should help you understand what might be going wrong. Try and figure out where the problem is on your own, before looking at the solution.

> ### {% icon solution %} Solution
>
> The problem is that the tests are expecting a different HID than what is being made available to them.  Since all the failing tests have to do with Lists, you might want to check out the ListCollectionCreator.
> 
> For comparisons sake, open one of the other CollectionCreators, as well and see if you can spot anything that is amiss.
>
> Pay careful attention in the template to the <collection-creator> props.
> 
> Here you should notice that in the PairCollectionCreator or in the PairedListCollectionCreator, `:hide-source-items="hideSourceItems"` instead of `:hide-source-items="!hideSourceItems"`. It's a little change, just flipping the boolean, but it makes a big difference! 
{: .solution }

Awesome! If you're here, you've found and fixed the problem, at least, in theory.  

If you think you've got it, try running the tests locally again to be sure. In this simulated case, you can get by with just running one of the tests, but in reality, should probably run all three Selenium tests to be sure that it's totally fixed. 

### Finishing up

Excellent! If you're here, you've found the problem, ran it locally, fixed the problem, ran it locally again without errors, and are ready to commit your change. Write a descriptive commit message like `fixed a logic error in ListCollectionCreatory` or whatever you like, as long as it's clear what you've done.

Selenium tests are in the bag!

## Runtime Error

TODO

## Finally --- Back in the Green

## Conclusion

First of all, thank you for completing this tutorial. We hope you feel more confident debugging your code on Galaxy and understanding how to navigate the testing frameworks we utilize. Of course, we are always available to help answer questions and support you!
