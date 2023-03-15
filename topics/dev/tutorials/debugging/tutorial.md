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
time_estimation: "4h"
key_points:
  - "Always run individual failing tests locally to minimize the feedback loop."
  - "To fix a bug, start by looking for the precise location in the code where the unexpected result is being produced, starting from the failing test or the error message in the log. Once the location is identified, finding and correcting the error is often trivial."
  - "When debugging an API test failure, first determine what endpoint is being utilized; then look in `lib/galaxy/webapps/galaxy/api/` to locate the controller that handles that enpoint. The controller will point you to a manager module, which is where you're most likely to find the bug."
  - "To locate the probable cause of a Selenium test failure, run Selenium with `GALAXY_TEST_SELENIUM_HEADLESS=false` to watch the browser step through the test."
  - "When debugging a runtime error, use pdb or any other interactive Python debugger to step through the code to identify the location of the error. Remember to use `GALAXY_CONFIG_DEBUG=1` to enable debugging."
  - "If you find a bug in a code block that is not covered by a test, before you fix the bug, write a test exposing the bug; run the test to ensure it fails; then fix the bug; then run the test again to ensure it passes."
  - "Write useful and property formatted commit messages."
  - "Pdb is your friend."
contributors:
  - assuntad23
  - jdavcs
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


In this tutorial we will demonstrate how to find and fix common types of bugs you may encounter as a contributor to Galaxy. We will step you through the process of finding and fixing a bug - from locating specific errors in the logs of Galaxy's GitHub Actions, to identifying their cause, developing a solution and committing your edits

With the skills from this tutorial, it is our hope that you will feel more prepared to develop solutions for Galaxy and more confidently navigate any obstacles along the way.

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Local development environment setup

For almost every type of error you encounter, there will be a relatively similar process:

1. Locate failing test or error on GitHub Actions
2. Run failing test locally
3. Identify issue causing failure
4. Develop a solution
5. Ensure test passes locally
6. Commit and push fixed branch to local fork

The process for dealing with runtime errors will be somewhat different, since there will be no failed test to review initially. However, in either case, the first step is to setup your local development environment.

## Contributing to Galaxy

{% snippet topics/dev/faqs/contributing.md %}

In this tutorial, instead of cloning Galaxy's default branch, we will clone a branch that has been created for this tutorial and contains several bugs that you'll need to identify and fix. Also, we won't be opening a pull request in the upstream Galaxy repository.

## Getting started: Pulling the buggy branch

Good news! All of the failing tests are on the same Galaxy branch! That means you only need to download and work off of one branch. The downside is that this is one buggy branch with issues all over the place. We have no doubts that you can solve them though!

> <hands-on-title>Setup your local Galaxy instance</hands-on-title>
>
> 1. Use GitHub UI to fork Galaxy's repository at `galaxyproject/galaxy`.
> 2. Clone your forked repository to a local path, further referred to as `GALAXY_ROOT` and `cd` into `GALAXY_ROOT`. Note that we specify the tutorial branch with the `-b` option:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > git clone -b training2022 https://github.com/<your-username>/galaxy GALAXY_ROOT
>    > cd GALAXY_ROOT
>    > ```
>    {: .code-in}
>
>    Note: In the future, when you contribute to Galaxy, you'll need to clone the `dev` branch, which is the default, so you don't need to specify the `-b` option in the `git clone` command).
>
> 3. Before we can use Galaxy, we need to create a virtual environment and install the required dependencies. This is generally done with the `common_startup.sh` script:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > bash scripts/common_startup.sh --dev-wheels
>    > ```
>    {: .code-in}
>
>    Make sure your Python version is at least 3.6 (you can check your Python version with `python --version`). If your system uses an older version, you may specify an alternative Python interpreter using the `GALAXY_PYTHON` environment variable (`GALAXY_PYTHON=/path/to/alt/python bash scripts/common_startup.sh --dev-wheels`).
>
> 4. Activate your new virtual environment:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > . .venv/bin/activate
>    > ```
>    {: .code-in}
>
>    Once activated, you'll see the name of the virtual environment prepended to your shell prompt: `(.venv)$`.
>
> 5. Finally, let's create a new branch for your edits:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > git checkout -b my-training
>    > ```
>    {: .code-in}
>
>    Now when you run `git branch` you'll see that your new branch is activated:
>
>    > > <code-in-title>Bash</code-in-title>
>    > > ```bash
>    > > git branch
>    > > ```
>    > {: .code-in}
>    >
>    > > <code-out-title></code-out-title>
>    > > ```bash
>    > >   training2022
>    > > * my-training
>    > > ```
>    > {: .code-out}
>    {: .code-2col}
>
>    Note: `my-training` is just an example; you can call your new branch anything you like.
{: .hands_on}

# Unit test failure

Most of Galaxy unit tests are designed to test a separate component or function, are very fast, and do not require a running Galaxy instance. Finding the issue causing a unit test to fail is relatively straightforward.

> <hands-on-title>Fixing a failing unit test</hands-on-title>
>
> 1. **Finding the failing test on GitHub**
>
>    Head over to the [Galaxy repository on GitHub](https://github.com/galaxyproject/galaxy) and in the branch drop-down select the "training" branch. Click the red "x" next to the last commit, you'll get a drop-down listing of workflow results; find "Unit tests / Test (3..." and click "Details". This will display a detailed report that includes information about the failed tests. Try to find the failed tests.
>
>    > <solution-title></solution-title>
>    > The log lists 2 failed tests:
>    > ```
>    > - FAILED lib/galaxy/model/tags.py::galaxy.model.tags.TagHandler.parse_tags
>    > - FAILED test/unit/util/test_utils.py::test_strip_control_characters
>    > ```
>    > The first is a doctest. While it will provide us with hints on what may have caused the error, the unit test is generally more helpful.
>    {: .solution }
>
> 2. **Running the test locally**
>
>    An essential requirement for productive development is a fast feedback loop: we want to make a change, run a test, and get immediate feedback. If we push our edits to the remote fork (or submit a pull request to the upstream repository) and wait for the the tests to run remotely, it could take hours before we get any feedback, which is very inefficient (and detrimental to sustaining the joy of programming).
>
>    Instead, we will run the tests locally. Furthermore, we will run a specific test that we know is failing: this will give us the instant feedback we need.
>
>    > <code-in-title>Bash</code-in-title>
>    > Make sure you are in `GALAXY_ROOT` and your virtual environment is activated. If not, activate it:
>    > ```bash
>    > . .venv/bin/activate
>    > ```
>    {: .code-in}
>
>    The simplest way to run the failing test locally is using pytest directly:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > pytest test/unit/util/test_utils.py
>    > ```
>    {: .code-in}
>
>    Do you see the error message?
>
>    > <solution-title></solution-title>
>    > The test gives us a very specific error message:
>    > ```
>    > AssertionError: assert 'b l a' == 'bla'
>    > - bla
>    > + b l a
>    > ?  + +
>    > ```
>    > This means that instead of the expected `bla` we got `b l a` (with spaces).
>    {: .solution }
>
> 3. **Finding and fixing the issue**
>
>    Our next step is to open the test file and see what part of the codebase (i.e., what module, class + method or function) is failing to provide the expected result.
>
>    The failing test is calling the `strip_control_characters` function in the `util` module. That function, apparently, is not producing the expected output. Let's head over to `lib/galaxy/util/__init__.py`, where we'll find the function definition:
>
>    ```python
>    def strip_control_characters(s):
>        """Strip unicode control characters from a string."""
>        return " ".join(c for c in unicodify(s) if unicodedata.category(c) != "Cc")
>    ```
>
>    The function takes a string as input, and, in theory, returns a copy of that string with certain characters stripped away (it doesn't really matter what characters are stripped away). However, from the test we see that the string is modified in a different way: spaces have been added between the characters.
>
>    Can you see what's causing this behavior?
>
>    > <solution-title></solution-title>
>    > The genexp is evaluated and subsequently joined using a space as a delimiter. The space is, obviously, a typo: we don't want to _add_ anything to the string, so the delimiter should be the empty string instead:
>    > ```python
>    > return "".join(c for c in unicodify(s) if unicodedata.category(c) != "Cc")
>    > ```
>    {: .solution }
>
>    Make the change, save the file, run the same test. Now it should pass.
>
> 4. **Finishing up**
>
>    Now you are ready to commit your edits and push the branch to your fork. So, `git add [your modified files]` and `git commit -m "[insert your commit message here]"`. Please consult [Galaxy's contributor's guide](https://github.com/galaxyproject/galaxy/blob/dev/CONTRIBUTING.md#how-to-contribute) for details on how to write a useful and properly formatted commit message, and why it matters.
>
>    Congratulations - you found and fixed a unit test failure!
{: .hands_on}

# API test failure

API tests test various aspects of the Galaxy API, as well as general backend aspects of Galaxy using the API. These tests require a Galaxy instance.

> <hands-on-title>Fixing a failing API test</hands-on-title>
>
> 1. **Finding the failing test on GitHub**
>
>    We'll follow the same process we used in the unit test section. On GitHub, click the same red "x" and select the "API tests / Test (3..." workflow and click "Details". Try to find the failed tests.
>
>    > <solution-title></solution-title>
>    > The log lists multiple failed tests. Some of them were caused by bugs you'll be fixing in later sections of this tutorial. For this exercise, we care about these two failed tests:
>    > ```
>    > FAILED lib/galaxy_test/api/test_licenses.py::LicensesApiTestCase::test_404_on_unknown_license
>    > FAILED lib/galaxy_test/api/test_licenses.py::LicensesApiTestCase::test_get_license
>    > ```
>    > Note: occasionally, there could be more failing tests; however, for the purposes of this exercise we are concerned with the licenses tests.
>    {: .solution }
>
> 2. **Running the test locally**
>
>    Again, we will run the individual test locally. While it is possible to run the test directly with pytest, in this case we will use Galaxy's script `run_tests.sh`. The script optimizes test runs by sharing the same instance of Galaxy across multiple tests; without the script a new Galaxy will be started for each TestCase class. There are other reasons to run the script; for more information check the documentation block at the top of the script.
>
>    Make sure you are in `GALAXY_ROOT` (if your virtual environment is not activated, the script will do it for you).
>
>    Let's run the failed test:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ./run_tests.sh -api lib/galaxy_test/api/test_licenses.py
>    > ```
>    {: .code-in}
>
>    Two tests are failing. Do you see the error messages?
>
>    > <solution-title></solution-title>
>    > We see the following error messages:
>    > ```
>    > FAILED lib/galaxy_test/api/test_licenses.py::LicensesApiTestCase::test_404_on_unknown_license - AssertionError: Request status code (500) was not expected value 404. Body was INVALID JSON RESPONSE <Internal Server Error>
>    > FAILED lib/galaxy_test/api/test_licenses.py::LicensesApiTestCase::test_get_license - AssertionError: Request status code (404) was not expected value 200. Body was {'err_msg': "License 'Apache-2.0' not found", 'err_code': 404001}
>    > ```
>    > This means that instead of HTTP status code 200 (OK/success) we got a 500 (server error) and 404 (resource not found). Sometimes multiple test failures may be caused by the same error; this is the case here.
>    {: .solution }
>
> 3. **Finding and fixing the issue**
>
>    We can use either one of the failing tests to start debugging. Let's use `test_get_license`. As with the unit test, first let's take a look at the test function definition (located in `lib/galaxy_test/api/test_licenses.py`):
>
>    ```python
>    def test_get_license(self):
>        licenseId = "Apache-2.0"
>        response = self._get(f"licenses/{licenseId}")
>        self._assert_status_code_is(response, 200)
>        self._assert_matches_license_id(licenseId, response.json())
>    ```
>
>    The important line here is line 3: it tells us what API endpoint is called: `licenses/{licenseId}`.
>
>    Our next step is to examine the controller for this endpoint. For that, we head over to `lib/galaxy/webapps/galaxy/api/licenses.py`.
>
>    If you look at the FastAPI controller, you'll see the `get` method that is responsible for fetching the data in response to a `GET` HTTP request received at the `licenses/{licenseId}` endpoint:
>
>    ```python
>        @router.get('/api/licenses/{id}',
>            summary="Gets the SPDX license metadata associated with the short identifier",
>            response_description="SPDX license metadata")
>        async def get(self, id=LicenseIdPath) -> LicenseMetadataModel:
>            """Returns the license metadata associated with the given
>            [SPDX license short ID](https://spdx.github.io/spdx-spec/appendix-I-SPDX-license-list/)."""
>            return self.licenses_manager.get_license_by_id(id)
>    ```
>
>    The only thing in this method that matters to us is the last line:
>
>    ```python
>    return self.licenses_manager.get_license_by_id(id)
>    ```
>
>    Why? Because despite all this syntax, the method simply redirects the request to `self.licenses_manager`. If we read the code, we'll quickly see that this variable holds an instance of `LicencesManager`, which (look at the import statements at the top of the file) is located in `lib/galaxy/managers/licenses.py`. Which is where we should look next.
>
>    If we carefully trace the code in `LicensesManager`, we'll find the error. Do you see it?
>
>    > <solution-title></solution-title>
>    > The `get_license_by_id()` method raises an `ObjectNotFound` error (which results in the 404 status code we see in the test result). This error is a direct result of the LicensesManager `get` method called on the previous line (`license = self.get(id)`) - take a look at that method:
>    > ```python
>    >     def get(self, uri):
>    >         if uri not in self._by_index:
>    >             return self._by_index[uri]
>    >             ...
>    > ```
>    > Let's rewrite the first clause in pseudo code to expose the problem:
>    > ```
>    > if key not in dictionary:
>    >     return dictionary[key]
>    > ```
>    > Well, that makes no sense whatsoever! It should be the opposite: if key **is** in dictionary, then get the item with that key. The `not` is clearly a typo (or, more accurately, a wicked edit made for the purpose of setting up a learning experience). Remove this negation and rerun the test.
>    {: .solution }
>
> 4. **Finishing up**
>
>    Well, that was fun. Congratulations! Now commit your edit (don't forget to supply a useful commit message), and move on to client errors.
{: .hands_on}

# Client linting error

Client Linting is important to reduce errors, improve quality of the codebase, and enforce common rules. These rules are defined by the eslinter defined with Vue best practices, primarily. On Github, the command `eslint` is run to start the process for linting code on the client side.

> <comment-title>Python Linting</comment-title>
> For backend consistency, we also have a Python linter. While we won't be running you through that exercise today, you may encounter a Python linting issue in the future. The way to diagnose and treat these is very similar to the way to do it for the front end.
{: .comment}

> <hands-on-title>Fixing client linting error</hands-on-title>
>
> 1. **Finding the failing test on GitHub**
>
>    With this error, one of the failing test workflows on GitHub would normally be “get_code_and_test”. However, since we are working on a branch that contains multiple bugs by design, some of the test failures overlap, as a result of which the “get_code_and_test” workflow has been skipped in the last commit.
>
>    Therefore, in this section, for the remote results, we'll have to describe what would have been displayed.
>
>    One of the failing tests on GitHub says "get_code_and_test". Clicking on the word "Details" assocated with this failure will open a brief view of the CircleCI Checks, with a list of tests that are run as part of that command. You'll see that the one that is failing is js_lint. And if you click that, you'll be navigated to the CicleCI app, where you can see much more information about the failures.
>
>    Here you can see the exact command that ran to lint our codebase:
>
>    ```
>    cd client && yarn run eslint
>    ```
>
>    We can also see that it has found 1 error, where that error is located, and a description of the problem.
>
> 2. **Running the linter locally**
>
>    Alright, so we've seen the failure, but we want to verify it locally now. For this, you'll need to open up a terminal.
>
>    > <code-in-title>Bash</code-in-title>
>    > Make sure you are in `GALAXY_ROOT` and your virtual environment is activated. If not, activate it:
>    > ```bash
>    > . .venv/bin/activate
>    > ```
>    {: .code-in}
>
>    Then you've got two options. You can run the same exact command that CircleCI uses, `cd client && yarn run eslint` or `make client-lint`, which does `cd client && yarn run eslint && yarn run prettier-check`. This will also make sure your code conforms to formatting standards.
>
>    After running that, you should see the same linting error from CircleCI; you should know where it is and how to fix it too!
>
> 3. **Finding and fixing the issue**
>
>    > <solution-title></solution-title>
>    >
>    > You are currently in `GALAXY_ROOT/client`.
>    > In the file `src/components/RuleBuilder/RuleComponent.vue`, there is an unused import. Simply delete the line that says `import Vue from "vue";`
>    >
>    > Alternatively, you can run `yarn eslint --fix` or `yarn run prettier`.
>    {: .solution }
>
>    If you think you've got it, try running the test locally again to be sure.
>
> 4. **Finishing up**
>
>    Whoo-hoo! If you're here, you've identified the error, ran the appropriate command locally, found the problem, fixed it, confirmed that the command runs without errors, and are ready to commit your change. Write a useful commit message like `Fix a client linting issue in RuleComponent` or something similar, making sure it is clear what you've done.
>    (See [Galaxy's contributor's guide](https://github.com/galaxyproject/galaxy/blob/dev/CONTRIBUTING.md#how-to-contribute) for details on how to write a useful and properly formatted commit message.)
>
>    Awesome work - you now understand how to solve client linting issues!
{: .hands_on}

# Client unit test failure

Client tests are tests written by developers to test front-end code. In the case of the test we're going to walk through right now, we're looking at a Vue test failure, which means the component that it's testing is definitely a Vue component.

> <hands-on-title>Fixing a failing client unit test</hands-on-title>
>
> 1. **Finding the failing test on GitHub**
>    In this section we, again, can't demonstrate the remote output, so we descibe it instead. 
>    One of the failing tests on GitHub says "Client Unit Testing / jest". Clicking on Details beside that failure, will open up a the terminal output from that test. Here you should be able to see what test is failing.
>
>    > <solution-title></solution-title>
>    >
>    > The failing test file is `src/components/RuleBuilder/SavedRulesSelector.test.js`, and the failing test in that file is `SavedRulesSelector › disables history icon if there is no history`.
>    {: .solution }
>
> 2. **Running the test locally**
>
>    > <code-in-title>Bash</code-in-title>
>    > Make sure you are in `GALAXY_ROOT` and your virtual environment is activated. If not, activate it:
>    > ```bash
>    > . .venv/bin/activate
>    > ```
>    {: .code-in}
>
>    To run the client tests locally, you must be in the client directory; if you're just in the galaxy directory, `cd client` should take you to the right place.
>
>    Next, type `yarn run test` or `yarn jest` to run all the client tests.
>
>    The output should match what you found on Github.
>
> 3. **Finding and fixing the issue**
>
>    Next up, finding the error. Here are some tips:
>
>    - Take a look at the failing test and the imports to find out which Vue component the test correlates to.
>    - Look at the test to see what it does. Is there an error in the test?
>    - Look at the Vue component where the failure seems to be taking place. Is the error located there?
>
>    The solution is hidden below, but try your hand at fixing it first:
>
>    > <solution-title></solution-title>
>    > The failure is only happening on the disabled test case, so take a look in the template of `SavedRulesSelector` where the conditions for adding disabled are defined.
>    >
>    > The way it's worded in the test (and the way one could think about it logically), is that if there is no history, the button should be disabled, but the code says
>    >
>    > ```javascript
>    > :class="{ disabled: numOfSavedRules == 1 }"
>    > ```
>    >
>    > -- whoops! A very silly off by one error. On that very same line, change the 1 to a 0, and run the test again.
>    {: .solution }
>
>    If you think you've got it, try running the test locally again to be sure.
>
>    > <comment-title>Running fewer client tests</comment-title>
>    > If you don't want to run the whole suite of client tests, you can add keywords that match the test path/name. For example both `yarn run jest rule` and `yarn jest selector` will work, but the former will run other tests with rule in the path/module name.
>    {: .comment}
>
> 4. **Finishing up**
>
>    Nice Job! If you're here, you've found the problem, ran it locally, fixed the problem, ran it locally again without errors, and are ready to commit your change. Write a useful commit message like `Fix logic error in SavedRuleSelector` or something similar, making sure it is clear what you've done.
>
>    (See [Galaxy's contributor's guide](https://github.com/galaxyproject/galaxy/blob/dev/CONTRIBUTING.md#how-to-contribute) for details on how to write a useful and properly formatted commit message.)
>
>    Awesome work - you now understand how to solve client linting issues!- you now have a handle on client tests!
{: .hands_on}

# Selenium test failure

Selenium is the end-to-end or integration testing framework that we use; so these tests often are designed to test both frontend and backend code. This makes fixing them a little difficult, but we are sure you can do it!

> <hands-on-title>Fixing a failing Selenium test</hands-on-title>
>
> 1. **Finding the failing test on GitHub**
>
>    There may be actually a few failing Selenium tests, but the good news is, they are the result of only one little error. The failing tests on GitHub are in the "Selenium Test 3.7,0" and "Selenium Test 3.7,1" workflows. When you click on the details of this one, it's not always very clear which test is causing the failure, but if you download the artifacts, it's very clear which tests have failed.
>
>    Try to identify which tests are failing based on artifacts you downloaded.
>
>    > <solution-title></solution-title>
>    >
>    > The failing tests are:
>    >
>    > - `test_build_list_and_show_items`
>    > - `test_upload_pair_specify_extension`
>    {: .solution }
>
> 2. **Running the tests locally**
>
>    Make sure you are in `GALAXY_ROOT` (if your virtual environment is not activated, the script will do it for you).
>
>    We do not want to run the full suite of Selenium tests here - that would take a lot of time and we already know that most of them passed.
>
>    Instead, we can build our commands to only run a subset of the tests. Here's an example: `./run_tests.sh -selenium lib/galaxy_test/selenium/test_sign_out.py`
>
>    A command like that will still run all of the tests defined in that python file, so we need to find out where the tests that failed are defined. Try seaching the baseline for our three failing Selenium tests; what files are they in?
>
>    > <solution-title></solution-title>
>    >
>    > The failing tests are defined in these files:
>    >
>    > - `lib/galaxy_test/selenium/test_collection_builders.py`
>    > - `lib/galaxy_test/selenium/test_uploads.py`
>    {: .solution }
>
>    Now that we have the files, we can run the selenium tests. The command to do this should be very similar to `./run_tests.sh -selenium lib/galaxy_test/selenium/test_sign_out.py`. The only difference should be the python file where the tests are defined. Try it yourself, before checking out the solution.
>
>    > <solution-title></solution-title>
>    >
>    > Either of these two commands will work:
>    >
>    > > <code-in-title>Bash</code-in-title>
>    > > ```bash
>    > > ./run_tests.sh -selenium lib/galaxy_test/selenium/test_collection_builders.py
>    > > ```
>    > > or
>    > > ```bash
>    > > ./run_tests.sh -selenium lib/galaxy_test/selenium/test_uploads.py
>    > > ```
>    > {: .code-in}
>    {: .solution }
>
>    Running any of these tests will show us that the tests are still failing on our local machine. So, you should see a failing test that matches one of the three we identified from the downloaded artifacts. Onward to fixing them!
>
> 3. **Finding and fixing the issue**
>
>    Good news - there's an easy way to figure out what's going wrong with a selenium test!
>
>    In your terminal type `export GALAXY_TEST_SELENIUM_HEADLESS=false`. This will allow you to see the test running in real time.
>
>    Re-run the selenium test to watch the headless browser step through the test.
>
>    Looking at the test run in the headless browser, alongside the definition for the test should help you understand what might be going wrong. Try and figure out where the problem is on your own, before looking at the solution.
>
>    > <solution-title></solution-title>
>    >
>    > The problem is that the tests are expecting a different HID (history identifier) than what is being made available to them. Since all the failing tests have to do with Lists, you might want to check out the ListCollectionCreator.
>    >
>    > For comparisons sake, open one of the other CollectionCreators, as well and see if you can spot anything that is amiss.
>    >
>    > Pay careful attention in the template to the <collection-creator> props.
>    >
>    > Here you should notice that in the PairCollectionCreator or in the PairedListCollectionCreator, `:hide-source-items="hideSourceItems"` instead of `:hide-source-items="!hideSourceItems"`. It's a little change, just flipping the boolean, but it makes a big difference!
>    {: .solution }
>
>    Awesome! If you're here, you've found and fixed the problem, at least, in theory.
>
>    If you think you've got it, try running the tests locally again to be sure. In this simulated case, you can get by with just running one of the tests, but in reality, you should probably run all three Selenium tests to be sure that it's totally fixed.
>
> 4. **Finishing up**
>
>    Excellent! If you're here, you've found the problem, ran it locally, fixed the problem, ran it locally again without errors, and are ready to commit your change. Write a useful commit message like `Fix logic error in ListCollectionCreator` or something similar, making sure it is clear what you've done.
>    (See [Galaxy's contributor's guide](https://github.com/galaxyproject/galaxy/blob/dev/CONTRIBUTING.md#how-to-contribute) for details on how to write a useful and properly formatted commit message.)
>
>    Selenium tests are in the bag!
{: .hands_on}

# Runtime error

Our last error happens at runtime, which means we don't have a failing test; instead we have a bug report describing the error. Our goal is to fix the error. Our steps are as follows:

1. Reproduce the error
2. Locate the error
3. Identify the problem
4. Write a test exposing the bug. The test must fail.
5. Fix the error.
6. Ensure the error no longer occurs.

> <hands-on-title>Fixing a runtime error</hands-on-title>
>
> 1. **Reproduce the error**
>
>    Here's the bug report:
>
>    ```In the User menu, clicking the Datasets option causes an error message to be displayed on the page: "Uncaught exception in exposed API method".```
>
>    Make sure you are in `GALAXY_ROOT`. Then start your local Galaxy using the `uvicorn` command*:
>
>    > > <code-in-title>Bash</code-in-title>
>    > > ```bash
>    > > GALAXY_CONFIG_FILE=config/galaxy.yml uvicorn --app-dir lib --factory galaxy.webapps.galaxy.fast_factory:factory
>    > > ```
>    > {: .code-in}
>    >
>    > > <code-out-title></code-out-title>
>    > > The first startup may take a few minutes. Eventually, you'll see the following output (the PID will be different):
>    > >
>    > > ```
>    > > Uvicorn running on http://127.0.0.1:8000 (Press CTRL+C to quit)
>    > > ```
>    > {: .code-out}
>    {: .code-2col}
>
>    *NOTE: We are not using `run.sh` because that would prevent us from using the Python interactive debugger, which we will be utilizing in this section.
>
>    Now you can access your Galaxy instance from the browser at [`http://127.0.0.1:8000`](http://127.0.0.1:8000).
>
>    To reproduce this error, we need to access the User menu. For this, we need to be logged in to Galaxy:
>
>    1. Click on "Login or Register:
>    2. At the bottom of the form, click "Register here"
>    3. Fill out the form. You may use any data, as long as the email you provide is a valid email format and the password is at least 6 characters long. For this exercise, it doesn't matter what email, name or password you use.
>    4. Click "Create"
>
>    Now the User menu is available to us. Click "User", then select "Datasets" from the dropdown menu. You should be able to see the error message displayed on the page: "Uncaught exception in exposed API method".
>
>    Being able to reproduce a bug is good: this means you're on your way to fixing it, and you know what to look for! In a real scenario, it is not uncommon that a bug is not easily reproduceable, which makes fixing it a much more complicated task.
>
> 2. **Locate the problem**
>
>    To figure out what's happening, our best bet is to look at the Galaxy log (look at the terminal window from which you launched Galaxy). Can you see the error message?
>
>    > <solution-title></solution-title>
>    >
>    > You should see something like this:
>    >
>    > ```
>    > galaxy.web.framework.decorators ERROR 2021-06-08 17:52:08,350 [pN:main.web.1,p:1583948,w:1,m:0,tN:uWSGIWorker1Core1] Uncaught exception in exposed API method:
>    > Traceback (most recent call last):
>    >   File "lib/galaxy/web/framework/decorators.py", line 312, in decorator
>    >     rval = func(self, trans, *args, **kwargs)
>    >   File "lib/galaxy/webapps/galaxy/api/datasets.py", line 100, in index
>    >     raise Exception('This should not happen!')
>    > Exception: This should not happen!
>    > ```
>    >
>    > From this error log, we can easily tell that the error happens in the `lib/galaxy/webapps/galaxy/api/datasets.py` file on line 100. The specific error message doesn't tell us much - which is often what happens in the wild.
>    {: .solution }
>
>    Let's head over to the `datasets'py` file, line 100:
>
>    ```python
>     97         if str_as_bool(trans.app.config.get("show_datasets", "True")):
>     98             return self.service.index(trans, history_id, serialization_params, filter_query_params)
>     99         else:
>    100             raise Exception("This should not happen!")
>    ```
>
>    The good news is that the error happens inside an `if/else` block, which narrows down our search to line 97: that line evaluates to `False`, causing the `else` clause to execute. The bad news is that we have no idea what causes that line to evaluate to `False`.
>
> 3. **Identify the problem**
>
>    We will investigate using pdb - our trusted Python debugger!
>
>    pdb offers a wide range of functionality and is exceptionally useful for debugging at runtime. We encourage you to read its [documentation](https://docs.python.org/3/library/pdb.html); however for this exercise we will only use the built-in [breakpoint()](https://docs.python.org/3/library/functions.html#breakpoint) function, that drops us into pdb.
>
>    Add a `breakpoint()` statement to your code, right above the `if/else` block.
>
>    > <solution-title></solution-title>
>    >
>    > ```python
>    > 126         breakpoint()
>    >  97         if str_as_bool(trans.app.config.get("show_datasets", "True")):
>    >  98             return self.service.index(trans, history_id, serialization_params, filter_query_params)
>    >  99         else:
>    > 100             raise Exception("This should not happen!")
>    > ```
>    >
>    {: .solution }
>
>    Now we need to restart Galaxy. However, to use pdb, we need to enable Galaxy's debug configuration option. In general, it may be a good idea to enable this in your `galaxy.yml` configuration file (just copy or rename `config/galaxy.yml.sample` and uncomment the option you want to set). However, for the purposes of this exercise, it's enough to set the `GALAXY_CONFIG_DEBUG` environment variable when running Galaxy:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > GALAXY_CONFIG_DEBUG=1 GALAXY_CONFIG_FILE=config/galaxy.yml uvicorn --app-dir lib --factory galaxy.webapps.galaxy.fast_factory:factory
>    > ```
>    {: .code-in}
>
>    > <comment-title>Why we set the debug option</comment-title>
>    >
>    > The reason for setting the `debug` option for releases prior to 22.05 was quite esoteric: among other things, it prevents uWSGI (no longer used as of 22.05) from remapping `stdin` to `dev/null` which would prevent tools like pdb from running.
>    {: .comment}
>
>    Now your Galaxy is running in debug mode. Repeat your steps that reproduced the error. When you click "Datasets", head over to your terminal. At the bottom of the log you'll see something like this:
>
>    ```
>    > /home/sergey/2sandbox/galaxy/dev_training/lib/galaxy/webapps/galaxy/api/datasets.py(97)index()
>    -> if str_as_bool(trans.app.config.get('show_datasets', 'True')):
>    (Pdb)
>    ```
>
>    You're at the pdb prompt. From here you can execute Python code interactively. The logs printed to `stdout` can be distracting, so try to ignore them. Let's explore what's happening in the next line:
>
>    > <code-in-title>Pdb</code-in-title>
>    > ```
>    > (Pdb) trans
>    > <galaxy.webapps.base.webapp.GalaxyWebTransaction object at 0x7f3664eee2b0>
>    > (Pdb) trans.app
>    > <galaxy.app.UniverseApplication object at 0x7f36ccc18790>
>    > (Pdb) trans.app.config
>    > <galaxy.config.GalaxyAppConfiguration object at 0x7f369addcac0>
>    > (Pdb) trans.app.config.get
>    > <bound method CommonConfigurationMixin.get of <galaxy.config.GalaxyAppConfiguration object at 0x7f369addcac0>>
>    > (Pdb)
>    > ```
>    {: .code-in}
>
>    We could have just typed `trans.app.config.get`, of course; the steps above give us some additional context.
>
>    Now let's take a look at the definition of the `get` method. We need to find `CommonConfigurationMixin` which is used by `galaxy.config.GalaxyAppConfiguration`.
>
>    Head over to `lib/galaxy/config/__init__.py` and look for the class `CommonConfigurationMixin`. Once you find it, look for the definition of its `get` method. When you find it, look at the method signature. What is the meaning of the argument `'True'` that is passed to this method in `datasets.py` on line 97?
>
>    > <solution-title></solution-title>
>    >
>    > ```python
>    > 584     def get(self, key, default=None):
>    > 585         # Warning: the value of self.config_dict['foo'] may be different from self.foo
>    > 586         return self.config_dict.get(key, default)
>    > ```
>    >
>    > The argument `True` corresponds to the parameter `default`. In the method's body, we see that this value will be passed to `self.config_dict.get()`; `self.config_dict` is a Python dictionary, so the second argument to `get()` is the default value that will be retuned if `key` is not in `self.config_dict`.
>    {: .solution }
>
>    Back to `datasets.py`, line 97. Now we know that `True` is the default value that will be returned if the key `show_datasets` does not exist in the `config_dict` dictionary. Moving on:
>
>    > <code-in-title>Pdb</code-in-title>
>    > ```
>    > (Pdb) trans.app.config.get('show_datasets')
>    > (Pdb) trans.app.config.get('show_datasets', 'True')
>    > 'True'
>    > ```
>    {: .code-in}
>
>    So, we see that `show_datasets` does not exist and we get the default value, which is the string `'True'`.
>
>    Next step: `str_to_bool`. First, let's pass it the value `'True'`:
>
>    > <code-in-title>Pdb</code-in-title>
>    > ```
>    (Pdb) str_as_bool('True')
>    False
>    ```
>    {: .code-in}
>
>    That doesn't look right! Exit the debugger (type in `q`, then `Enter`, then `CTRL-C ` to exit Galaxy.
>    Head over to the function's definition (you can tell from the import statement at the top of the file that it's in `lib/galaxy/util/__init__.py`) and try to figure out if there is an error.
>
>    > <solution-title></solution-title>
>    >
>    > ```python
>    > 117 def str_as_bool(string):
>    > 118     """This is for training only."""
>    > 119     if str(string) in ("true", "yes", "on", "1"):
>    > 120         return True
>    > 121     else:
>    > 122         return False
>    > ```
>    >
>    > At first glance, everything seems right. But think of what value we are passing: `'True'` - it is _not_ in the `('true', 'yes', 'on', '1')` tuple! Aparently, the developer did not account for upper vs. lower case!
>    >
>    > Do NOT fix the error just yet.
>    {: .solution }
>
> 4. **Write a test exposing the bug**
>
>    Before you fix code that was not covered by a test, write the test to expose the bug! In this case we have a simple, completely isolated function - a simple unit test should do. Your test should go into this module: `test/unit/util/test_utils.py`.
>
>    Start with a blank test function. Then call the funtion under test (`str_as_bool`) and assert that when you pass it the value `'True'`, it returns `True`. Can you think of more cases to test?
>
>    > <solution-title></solution-title>
>    >
>    > Your first step could be something like this:
>    >
>    > ```python
>    > def test_str_as_bool():
>    >     assert util.str_as_bool('True')
>    > ```
>    >
>    > However, keep in mind that there may be other cases of mixed or uppercase values. You should also test for lowercase values - to make sure that while fixing the bug you don't break the current functionality. And remember to test for the `False` cases as well.
>    >
>    > There are many ways to implement this, feel free to write your own. Here's one option:
>    >
>    > ```python
>    > def test_str_as_bool():
>    >     for value in ('true', 'yes', 'on', '1', 'True', 'Yes', 'On','TRUE', 'YES', 'ON'):
>    >         assert util.str_as_bool(value)
>    >     for value in ('false', '0', 'something else'):
>    >         assert not util.str_as_bool(value)
>    > ```
>    >
>    {: .solution }
>
>    Run your test. Does it fail? Good! Now it's time to fix the error.
>
> 5. **Fix the error**
>
>    Can you think of a simple way to fix the function so that our test passes?
>
>    > <solution-title></solution-title>
>    >
>    > One way to do it is to simply convert the input value to lowercase:
>    >
>    > ```python
>    > 117 def str_as_bool(string):
>    > 118     """This is for training only."""
>    > 119     if str(string).lower() in ("true", "yes", "on", "1"):
>    > 120         return True
>    > 121     else:
>    > 122         return False
>    > ```
>    >
>    > In fact, that's exactly what the real function does! Check line 1017 - you'll see almost the exact same function, named `string_as_bool`: we didn't modify it because that might have broken other parts of Galaxy which would have been a distraction from the core aspects of the exercise.
>    {: .solution }
>
>    Run the test - now it should pass!
>
> 6. **Ensure the error no longer occurs**
>
>    Now that we have the `str_as_bool` function covered by test, as a bonus, we can safely do some minor refactoring: the test will prevent us from breaking things. Can you simplify the function's code?
>
>    > <solution-title></solution-title>
>    >
>    > ```python
>    > 117 def str_as_bool(string):
>    > 118     """ This is for training only."""
>    > 119     return str(string).lower() in ('true', 'yes', 'on', '1')
>    > ```
>    >
>    > Much better!
>    {: .solution }
>
>    Finally, head back to your local Galaxy and verify that the runtime error no longer occurs:
>
>    1. Remove the `breakpoint()` statement you added to `lib/galaxy/webapps/galaxy/api/datasets.py`.
>    2. Start Galaxy and repeat the steps from the bug report (select "User" > "Datasets"). You should see the Datasets page now - no more error message!
>
>    Congratulations - you've completed the last and, possibly, most challenging section of this tutorial!
{: .hands_on}

# Conclusion

Finally, we're back in the Green - all tests pass!

First of all, thank you for completing this tutorial. We hope you feel more confident debugging Galaxy code and understanding how to navigate the testing frameworks we utilize. Of course, we are always available to help answer questions and support you!
