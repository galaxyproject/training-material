---
layout: tutorial_hands_on

title: "Debugging Galaxy"
questions:
  - "How do I debug Galaxy?"
objectives:
  - "Fix a broken branch"
  - "Interpret the results of failed tests on GitHub"
  - "Run individual tests locally"
  - "Fix errors identified by failing tests"
  - "Debug simple runtime errors using the Python debugger"
  - "Write tests exposing the identified bug"
time_estimation: "2h"
key_points:
  - "TODO Unit Tests"
  - "TODO API Tests"
  - "TODO Client Linting Failures"
  - "TODO Client Tests"
  - "TODO Selenium Tests"
  - "TODO Runtime Errors"
contributors:
  - assuntad23
  - ic4f
  - jmchilton
requirements:
  - type: none
    title: "Familiarity with basic Git commands"
  - type: none
    title: "Basic knowledge of Python and JavaScript"
  - type: none
    title: "Mac OS or Linux that can run Galaxy & your favorite IDE or editor"
  - type: "internal"
    topic_name: dev
    tutorials:
      - architecture
subtopic: core
---

# Introduction
{:.no_toc}


In this tutorial we will demonstrate how to find and fix common types of bugs you may encounter as a contributor to Galaxy. We will step you through the process of finding and fixing a bug - from locating specific errors in the logs of Galaxy's GitHub workflows, to identifying their cause, developing a solution and committing your edits.

With the skills from this tutorial, it is our hope that you will feel more prepared to develop solutions for Galaxy and more confidently navigate any obstacles along the way.

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Local development environment setup

For almost every type of error you encounter, there will be a relatively similar process:

1. Locate failing test or error on GitHub workflows
2. Run failing test locally
3. Identify issue causing failure
4. Develop a solution
5. Ensure test passes locally
6. Commit and push fixed branch to local fork

The process for dealing with runtime errors will be somewhat different, since there will be no failed test to review initially. However, in either case, the first step is to setup your local development environment.

## Contributing to Galaxy

To contribute to galaxy, a GitHub account is required. Changes are proposed via a [pull request](https://docs.github.com/en/github/collaborating-with-pull-requests). This allows the project maintainers to review the changes and suggest improvements.

The general steps are as follows:
1. Fork the Galaxy repository
2. Clone your fork
3. Make changes in a new branch
4. Commit your changes, push branch to your fork
5. Open a pull request for this branch in the upstream Galaxy repository

In this tutorial, instead of cloning Galaxy's default branch, we will clone a branch that has been created for this tutorial and contains several bugs that you'll need to identify and fix. Also, we won't be opening a pull request in the upstream Galaxy repository.

## Getting started: Pulling the buggy branch

Good news! All of the failing tests are on the same Galaxy branch! That means you only need to download and work off of one branch. The downside is that this is one buggy branch with issues all over the place. We have no doubts that you can solve them though!

> ### {% icon hands_on %} Hands-on: Setup your local Galaxy instance
> 
> 1. Use GitHub UI to fork Galaxy's repository at `galaxyproject/galaxy`. 
> 2. Clone your forked repository to a local path, further referred to as GALAXY_ROOT and `cd` into GALAXY_ROOT. Note that we specify the tutorial branch with the `-b` option:
> ```shell
> $ git clone -b training_dev_debugging https://github.com/galaxyproject/galaxy GALAXY_ROOT
> $ cd GALAXY_ROOT
> ```
> Note: In the future, when you contribute to Galaxy, you'll need to clone the `dev` branch, which is the default, so you don't need to specify the `-b` option in the `git clone` command).
> 
> 3. Before we can use Galaxy, we need to create a virtual environment and install the required dependencies. First, let's create a virtual environment:
> ```shell
> $ virtualenv .venv
> ```
> Make sure your Python version is at least 3.6 (you can check your Python version with `python --version`). If your system uses an older version, you may specify an alternative Python interpreter using the `-p` option:
> ```shell
> $ virtualenv -p PATH-TO-PYTHON-INTERPETER .venv
> ```
> 4. Activate your new virtual environment:
> ```shell
> . .venv/bin/activate
> ```
> Once activated, you'll see the name of the virtual environment prepended to your shell prompt: `(.venv)$`.
> 5. Now you are ready to install the dependencies. Since you wil be doing development, you need to install the dependencies listed in `lib/galaxy/dependencies/dev-requirements.txt`:
> ```shell
> $ pip install -r lib/galaxy/dependencies/dev-requirements.txt
> ```
> 6. Finally, let's create a new branch for your edits:
> ```shell
> $ git checkout -b my-training
> ```
> Now when you run `git branch` you'll see that your new branch is activated:
> ```shell
>   training_dev_debugging
> * my-training
> ```
> Note: `my-training` is just an example; you can call your new branch anything you like.
{: .hands_on}

# Unit test failure

Most of Galaxy unit tests are designed to test a separate component or function, are very fast, and do not require a running Galaxy instance. Finding the issue causing a unit test to fail is relatively straightforward. 

> ### {% icon hands_on %} Hands-on: Fixing a failing unit test
> 
> 1. __Finding the failing test on GitHub__
>    > 
>    > Head over to the [Galaxy repository on GitHub](https://github.com/galaxyproject/galaxy) and in the branch drop-down select the "training_dev_debugging" branch. Click the red "x" next to the last commit, you'll get a drop-down listing of workflow results; find "Unit tests / Test (3..." and click "Details". This will display a detailed report that includes information about the failed tests. Try to find the failed tests.
>    > 
>    >  > ### {% icon solution %} Solution
>    >  > 
>    >  > The log lists 2 failed tests:
>    >  > ```
>    >  > - FAILED lib/galaxy/model/tags.py::galaxy.model.tags.TagHandler.parse_tags
>    >  > - FAILED test/unit/util/test_utils.py::test_strip_control_characters
>    >  > ```
>    >  > The first is a doctest. While it will provide us with hints on what may have caused the error, the unit test is generally more helpful. 
>    >  {: .solution }
>    >  
> 2. __Running the test locally__
>    > 
>    > An essential requirement for productive development is a fast feedback loop: we want to make a change, run a test, and get immediate feedback. If we push our edits to the remote fork (or submit a pull request to the upstream repository) and wait for the the tests to run remotely, it could take hours before we get any feedback, which is very inneficient (and detrimental to sustaining the joy of programming).
>    > 
>    > Instead, we will run the tests locally. Furthermore, we will run a specific test that we know is failing: this will give us the instant feedback we need.
>    > 
>    > Make sure your virtual environment is activated. If not, activate it: `. .venv/bin/activate`
>    > 
>    > The simplest way to run the failing test locally is using pytest directly:
>    > ```shell
>    > $ pytest test/unit/util/test_uril.py
>    > ```
>    > Do you see the error message?
>    > 
>    >  > ### {% icon solution %} Solution
>    >  > 
>    >  > The test gives us a very specific error message:
>    >  > ```
>    >  > AssertionError: assert 'b l a' == 'bla'
>    >  > - bla
>    >  > + b l a
>    >  > ?  + +
>    >  > ```
>    >  > This means that instead of the expected `bla` we got `b l a` (with spaces). 
>    >  {: .solution }
>    >  
> 3. __Finding and fixing the issue__
>    > 
>    > Our next step is to open the test file and see what part of the codebase (i.e., what module, class + method or function) is failing to provide the expected result. 
>    > 
>    > The failing test is calling the `strip_control_characters` function in the `util` module. That function, apparently, is not producing the expected output. Let's head over to `lib/galaxy/utl/__init__.py`, where we'll find the function definition:
>    > 
>    > ```python
>    > def strip_control_characters(s):
>    >     """Strip unicode control characters from a string."""
>    >     return " ".join(c for c in unicodify(s) if unicodedata.category(c) != "Cc")
>    > ```
>    > The function takes a string as input, and, in theory, returns a copy of that string with certain characters stripped away (it doesn't really matter what characters are stripped away). However, from the test we see that the string is modified in a different way: spaces have been added between the characters. Can you see what's causing this behavior?
>    > 
>    >  > ### {% icon solution %} Solution
>    >  > 
>    >  > The genexp is evaluated and subsequently joined using a space as a delimiter. The space is, obviously, a typo: we don't want to *add* anything to the string, so the delimiter should be the empty string instead:
>    >  > 
>    >  > ```python
>    >  > return "".join(c for c in unicodify(s) if unicodedata.category(c) != "Cc")
>    >  > ```
>    > {: .solution }
>    > 
>    > Make the change, save the file, run the same test. Now it should pass.
>    >
> 4. __Finishing up__
>    > 
>    > Now you are ready to commit your edits and push the branch to your fork. So, `git add [your modified files]` and `git commit -m "[insert your commit message here]"`. Please consult [Galaxy's contributor's guide](https://github.com/galaxyproject/galaxy/blob/dev/CONTRIBUTING.md#how-to-contribute) for details on how to write a useful and properly formatted commit message, and why it matters.
>    > 
>    > Congratulations - you found and fixed a unit test failure!
{: .hands_on}

# API test failure

API test various aspects of the Galaxy API, as well as general backend aspects of Galaxy using the API. These tests require a Galaxy instance.

> ### {% icon hands_on %} Hands-on: Fixing a failing API test
> 
> 1. __Finding the failing test on GitHub__
>    > 
>    > We'll follow the same process we used in the unit test section. On GitHub, click the same red "x" and select the "API tests / Test (3..." workflow and click "Details". Try to find the failed tests.
>    > 
>    >  > ### {% icon solution %} Solution
>    >  > 
>    >  > The log lists 2 relevant failed tests:
>    >  > ```
>    >  > FAILED lib/galaxy_test/api/test_licenses.py::LicensesApiTestCase::test_404_on_unknown_license
>    >  > FAILED lib/galaxy_test/api/test_licenses.py::LicensesApiTestCase::test_get_license
>    >  > ```
>    >  > Note: occasionally, there could be more failing tests; however, for the purposes of this exercise we are concerned with the licenses tests.
>    > {: .solution }
>    > 
> 2. __Running the test locally__
>    > 
>    > Again, we will run the individual test locally. While it is possible to run the test directly with pytest, in this case we will use Galaxy's script `run_tests.sh`. The script optimizes test runs by sharing the same instance of Galaxy across multiple tests; without the script a new Galaxy will be started for each TestCase class. There are other reasons to run the script; for more information check the documentation block at the top of the script.
>    > 
>    > Let's run the failed test:
>    > ```shell
>    > ./run_tests.sh -api lib/galaxy_test/test_licenses.py
>    > ```
>    > Two tests are failing. Do you see the error messages?
>    > 
>    >  > ### {% icon solution %} Solution
>    >  > 
>    >  > We see the following error messages:
>    >  > ```
>    >  > FAILED lib/galaxy_test/api/test_licenses.py::LicensesApiTestCase::test_404_on_unknown_license - AssertionError: Request status code (500) was not expected value 404. Body was INVALID JSON RESPONSE <Internal Server Error>
>    >  > FAILED lib/galaxy_test/api/test_licenses.py::LicensesApiTestCase::test_get_license - AssertionError: Request status code (404) was not expected value 200. Body was {'err_msg': "License 'Apache-2.0' not found", 'err_code': 404001}
>    >  > ```
>    >  > This means that instead of HTTP status code 200 (OK/success) we got a 500 (server error) and 404 (resource not found). Sometimes multiple test failures may be caused by the same error; this is the case here.
>    > {: .solution }
>    > 
> 3. __Finding and fixing the issue__
>    > 
>    > We can use either one of the failing tests to start debugging. Let's use `test_get_license`. As with the unit test, first let's take a look at the test function definition (located in `lib/galaxy_test/api/test_licenses.py`):
>    > 
>    > ```python
>    > def test_get_license(self):
>    >     licenseId = "Apache-2.0"
>    >     response = self._get(f"licenses/{licenseId}")
>    >     self._assert_status_code_is(response, 200)
>    >     self._assert_matches_license_id(licenseId, response.json())
>    > ```
>    > The important line here is line 3: it tells us what API endpoint is called: `licenses/{licenseId}`.
>    > 
>    > Our next step is to examine the controller for this endpoint. For that, we head over to `lib/galaxy/webapps/galaxy/api/licenses.py`.
>    > 
>    > There are 2 controllers (a legacy controller and a FastAPI controller) that have identical functionality; we can use either one. If you look at the FastAPI controller, you'll see the `get` method that is responsible for fetching the data in response to a `GET` HTTP request received at the licenses/{licenseId}` endpoint:
>    > 
>    > ```python
>    >     @router.get('/api/licenses/{id}',
>    >         summary="Gets the SPDX license metadata associated with the short identifier",
>    >         response_description="SPDX license metadata")
>    >     async def get(self, id=LicenseIdPath) -> LicenseMetadataModel:
>    >         """Returns the license metadata associated with the given
>    >         [SPDX license short ID](https://spdx.github.io/spdx-spec/appendix-I-SPDX-license-list/)."""
>    >         return self.licenses_manager.get_license_by_id(id)
>    > ```
>    > 
>    > The only thing in this method that matters to us is the last line: `return self.licenses_manager.get_license_by_id(id)` (Why? Because despite all this syntax, the method simply redirects the request to `self.licenses_manager`). If we read the code, we'll quickly see that this variable holds an instance of `LicencesManager`, which (look at the import statements at the top of the file) is located in `lib/galaxy/managers/licenses`. Which is where we should look next.
>    > 
>    > If we carefully trace the code in `LicensesManager`, we'll find the error. Do you see it?
>    > 
>    >  > ### {% icon solution %} Solution
>    >  > 
>    >  > The `get_licenses_by_id()` method raises an `ObjectNotFound` error (which results in the 404 status code we see in the test result). This error is a direct result of the `get` method called on the previous line - take a look at that method: 
>    >  > 
>    >  > ```python
>    >  >     def get(self, uri):
>    >  >         if uri not in self._by_index:
>    >  >             return self._by_index[uri]
>    >  >             ...
>    >  > ```
>    >  > Let's rewrite the first clause in pseudo code to expose the problem:
>    >  > ```
>    >  > if key not in dictionary: 
>    >  >     return dictionary[key]
>    >  > ```
>    >  > Well, that makes no sense whatsoever! It should be the opposite: if key **is** in dictionary, then get the item with that key. The `not` is clearly a typo (or, more accurately, a wicked edit made for the purpose of setting up a learning experience). Remove this negation and rerun the test.  
>    >  > 
>    > {: .solution }
>    > 
> 4. __Finishing up__
>    > 
>    > Well, that was fun. Congratulations! Now commit your edit (don't forget to supply a useful commit message), and move on to client errors.
{: .hands_on}

# Client linting error

Client Linting is important to reduce errors, improve quality of the codebase, and enforce common rules. These rules are defined by the eslinter defined with Vue best practices, primarily. On Github, the command `eslint` is run to start the process for linting code on the client side.

> ### {% icon comment %} Python Linting
> For backend consistency, we also have a Python linter. While we won't be running you through that exercise today, you may encounter a Python linting issue in the future. The way to diagnose and treat these is very similar to the way to do it for the front end.
{: .comment}

> ### {% icon hands_on %} Hands-on: Fixing client linting error
> 
> 1. __Finding the failing test on GitHub__
> 
>    > One of the failing tests on GitHub says "get_code_and_test". Clicking on the word "Details" assocated with this failure will open a brief view of the CircleCI Checks, with a list of tests that are run as part of that command. You'll see that the one that is failing is js_lint. And if you click that, you'll be navigated to the CicleCI app, where you can see much more information about the failures.
>    > 
>    > Here you can see the exact command that ran to lint our codebase; it's `cd client && yarn run eslint`.
>    > 
>    > We can also see that it has found 1 error, where that error is located, and a description of the problem.
>    > 
> 2. __Running the linter locally__
>    > 
>    > Alright, so we've seen the failure, but we want to verify it locally now. For this, you'll need to open up a terminal.
>    > 
>    > Make sure your virtual environment is activated. If not, activate it: `. .venv/bin/activate`
>    > 
>    > Then you've got two options. You can run the same exact command that CircleCI uses `cd client && yarn run eslint` or `make client-lint`, which does `cd client && yarn run eslint && yarn run prettier-check`. This will also make sure your code conforms to formatting standards.
>    > 
>    > After running that, you should see the same linting error from CircleCI; you should know where it is and how to fix it too!
>    > 
> 3. __Finding and fixing the issue__
>    > 
>    > > ### {% icon solution %} Solution
>    > > In the file `client/src/components/RuleBuilder/RuleComponent`, there is an unused import. Simply delete the line that says `import Vue from "vue";` 
>    > {: .solution }
>    > 
>    > If you think you've got it, try running the test locally again to be sure.
>    > 
> 4. __Finishing up__
>    > 
>    > Whoo-hoo! If you're here, you've identified the error, ran the appropriate command locally, found the problem, fixed it, confirmed that the command runs without errors, and are ready to commit your change. Write a useful commit message like `Fix a client linting issue in RuleComponent` or something similar, making sure it is clear what you've done.
>    > (See [Galaxy's contributor's guide](https://github.com/galaxyproject/galaxy/blob/dev/CONTRIBUTING.md#how-to-contribute) for details on how to write a useful and properly formatted commit message.)
>    > 
>    > Awesome work - you now understand how to solve client linting issues!
{: .hands_on}

# Client unit test failure

Client tests are tests written by developers to test front-end code. In the case of the test we're going to walk through right now, we're looking at a Vue test failure, which means the component that it's testing is definitely a Vue component.

> ### {% icon hands_on %} Hands-on: Fixing a failing client unit test
> 
> 1. __Finding the failing test on GitHub__
>    >  
>    >  One of the failing tests on GitHub says "Client Unit Testing / jest". Clicking on Details beside that failure, will open up a the terminal output from that test. Here you should be able to see what test is failing.
>    >  
>    >  > ### {% icon solution %} Solution
>    >  >
>    >  > The failing test file is `src/components/RuleBuilder/SavedRulesSelector.test.js`, and the failing test in that file is `SavedRulesSelector › disables history icon if there is no history`. 
>    >  {: .solution }
>    >  
> 2. __Running the test locally__
>    >  
>    > Make sure your virtual environment is activated. If not, activate it: `. .venv/bin/activate`
>    >  
>    >  To run the client tests locally, you must be in the client directory; if you're just in the galaxy directory, `cd client` should take you to the right place.
>    >  
>    >  Next, type `yarn run test` or `yarn jest` to run all the client tests.
>    >  
>    >  The output should match what you found on Github.
>    >  
> 3. __Finding and fixing the issue__
>    >  
>    >  Next up, finding the error. Here are some tips:
>    >  
>    >  - Take a look at the failing test and the imports to find out which Vue component the test correlates to.
>    >  - Look at the test to see what it does. Is there an error in the test?
>    >  - Look at the Vue component where the failure seems to be taking place. Is the error located there?
>    >  
>    >  The solution is hidden below, but try your hand at fixing it first:
>    >  
>    >  > ### {% icon solution %} Solution
>    >  >
>    >  > The failure is only happening on the disabled test case, so take a look in the template of `SavedRuleSelector` where the conditions for adding disabled are defined.  
>    >  > The way it's worded in the test (and the way one could think about it logically), is that if there is no history, the button should be disabled, but the code says 
>    >  > ```javascript
>    >  > :class="{ disabled: numOfSavedRules == 1 }"
>    >  > ``` 
>    >  > -- whoops! A very silly off by one error. On that very same line, change the 1 to a 0, and run the test again. 
>    >  {: .solution }
>    >  
>    >  If you think you've got it, try running the test locally again to be sure.
>    >  
>    >  > ### {% icon comment %} Running fewer client tests
>    >  >
>    >  > If you don't want to run the whole suite of client tests, you can add keywords that match the test path/name. For example both `yarn run jest rule` and `yarn jest selector` will work, but the former will run other tests with rule in the path/module name.  
>    >  {: .comment}
>    >  
> 4. __Finishing up__
>    >  
>    >  Nice Job! If you're here, you've found the problem, ran it locally, fixed the problem, ran it locally again without errors, and are ready to commit your change. Write a useful commit message like `Fix logic error in SavedRuleSelector` or something similar, making sure it is clear what you've done.
>    > (See [Galaxy's contributor's guide](https://github.com/galaxyproject/galaxy/blob/dev/CONTRIBUTING.md#how-to-contribute) for details on how to write a useful and properly formatted commit message.)
>    >  
>    >  Awesome work - you now understand how to solve client linting issues!- you now have a handle on client tests!
{: .hands_on}

# Selenium test failure

Selenium is the end-to-end or integration testing framework that we use; so these tests often are designed to test both frontend and backend code. This makes fixing them a little difficult, but we are sure you can do it!

> ### {% icon hands_on %} Hands-on: Fixing a failing Selenium test
> 
> 1. __Finding the failing test on GitHub__
>    >  
>    >  There are actually a few failing Selenium tests, but the good new is, they are the result of only one little error. One of the failing tests on GitHub says "Selenium Test 3.7, 0" and "Selenium Test 3.7, 1". When you click on the details of this one, it's not always very clear which test is causing the failure, but if you download the artifacts, it's very clear which tests have failed.
>    >  
>    >  Try to identify which tests are failing based on artifacts you downloaded.
>    >  
>    >  > ### {% icon solution %} Solution
>    >  >
>    >  > The failing tests are:
>    >  >
>    >  > - `test_build_list_and_show_items`
>    >  > - `test_data_column_input`
>    >  > - `test_upload_pair_specify_extension`
>    >  {: .solution }
>    >  
> 2. __Running the tests locally__
>    >  
>    >  We do not want to run the full suite of Selenium tests here - that would take a lot of time and we already know that most of them passed.
>    >  
>    >  Instead, we can build our commands to only run a subset of the tests. Here's an example: `./run_tests.sh -selenium lib/galaxy_test/selenium/test_sign_out.py`
>    >  
>    >  A command like that will still run all of the tests defined in that python file, so we need to find out where the tests that failed are defined. Try seaching the baseline for our three failing Selenium tests; what files are they in?
>    >  
>    >  > ### {% icon solution %} Solution
>    >  >
>    >  > The failing tests are defined in these files:
>    >  >
>    >  > - `lib/galaxy_test/selenium/test_collection_builders.py`
>    >  > - `lib/galaxy_test/selenium/test_workflow_editor.py`
>    >  > - `lib/galaxy_test/selenium/test_uploads.py`
>    >  {: .solution }
>    >  
>    >  Now that we have the files, we can run the selenium tests. The command to do this should be very similar to `./run_tests.sh -selenium lib/galaxy_test/selenium/test_sign_out.py`. The only difference should be the python file where the tests are defined. 
>    >  
>    >  Running any of these three tests will show us that the tests are still failing on our local machine. So, you should see a failing test that matches one of the three we identified from the downloaded artifacts. Onward to fixing them!
>    >  
> 3. __Finding and fixing the issue__
>    >  
>    >  Good news - there's an easy way to figure out what's going wrong with a selenium test! In your terminal type `export GALAXY_TEST_SELENIUM_HEADLESS=false`. This will allow you to see the test running in real time. Re-run the selenium test to watch the headless browser step through the test. Looking at the test run in the headless browser, alongside the definition for the test should help you understand what might be going wrong. Try and figure out where the problem is on your own, before looking at the solution.
>    >  
>    >  > ### {% icon solution %} Solution
>    >  >
>    >  > The problem is that the tests are expecting a different HID than what is being made available to them.  Since all the failing tests have to do with Lists, you might want to check out the ListCollectionCreator.
>    >  > 
>    >  > For comparisons sake, open one of the other CollectionCreators, as well and see if you can spot anything that is amiss.
>    >  >
>    >  > Pay careful attention in the template to the <collection-creator> props.
>    >  > 
>    >  > Here you should notice that in the PairCollectionCreator or in the PairedListCollectionCreator, `:hide-source-items="hideSourceItems"` instead of `:hide-source-items="!hideSourceItems"`. It's a little change, just flipping the boolean, but it makes a big difference! 
>    >  {: .solution }
>    >  
>    >  Awesome! If you're here, you've found and fixed the problem, at least, in theory.  
>    >  
>    >  If you think you've got it, try running the tests locally again to be sure. In this simulated case, you can get by with just running one of the tests, but in reality, you should probably run all three Selenium tests to be sure that it's totally fixed. 
>    >  
> 4. __Finishing up__
>    >  
>    >  Excellent! If you're here, you've found the problem, ran it locally, fixed the problem, ran it locally again without errors, and are ready to commit your change. Write a useful commit message like `Fix logic error in ListCollectionCreatory` or something similar, making sure it is clear what you've done.
>    > (See [Galaxy's contributor's guide](https://github.com/galaxyproject/galaxy/blob/dev/CONTRIBUTING.md#how-to-contribute) for details on how to write a useful and properly formatted commit message.)
>    >  
>    >  Selenium tests are in the bag!
{: .hands_on}

# Runtime error

TODO

# Conclusion

Finally, we're back in the Green - all tests pass!

First of all, thank you for completing this tutorial. We hope you feel more confident debugging Galaxy code and understanding how to navigate the testing frameworks we utilize. Of course, we are always available to help answer questions and support you!
