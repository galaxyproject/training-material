---
layout: tutorial_hands_on

title: "Writing Automated Tests for Galaxy"
questions:
  - TODO
objectives:
  - "Learn about the different types of automated tests in Galaxy"
  - "Learn to write unit tests"
  - "Learn to write API tests"
  - "Learn to write end-to-end tests"
time_estimation: "3h"
key_points:
  - TODO
contributors:
  - jdavcs
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

## Testing in Galaxy
{:.no_toc}

TODO

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}


## Local development environment setup

{% snippet topics/dev/faqs/contributing.md %}

> ### {% icon hands_on %} Hands-on: Setup your local Galaxy instance
>
> 1. Use GitHub UI to fork Galaxy's repository at `galaxyproject/galaxy`.
> 2. Clone your forked repository to a local path, further referred to as `GALAXY_ROOT` and `cd` into `GALAXY_ROOT`. Note that we specify the tutorial branch with the `-b` option:
>
>    > ### {% icon code-in %} Input: Bash
>    > ```bash
>    > git clone -b training https://github.com/<your-username>/galaxy GALAXY_ROOT
>    > cd GALAXY_ROOT
>    > ```
>    {: .code-in}
>
>    Note: In the future, when you contribute to Galaxy, you'll need to clone the `dev` branch, which is the default, so you don't need to specify the `-b` option in the `git clone` command).
>
> 3. Before we can use Galaxy, we need to create a virtual environment and install the required dependencies. This is generally done with the `common_startup.sh` script:
>
>    > ### {% icon code-in %} Input: Bash
>    > ```bash
>    > bash scripts/common_startup.sh --dev-wheels
>    > ```
>    {: .code-in}
>
>    Make sure your Python version is at least 3.6 (you can check your Python version with `python --version`). If your system uses an older version, you may specify an alternative Python interpreter using the `GALAXY_PYTHON` environment variable (`GALAXY_PYTHON=/path/to/alt/python bash scripts/common_startup.sh --dev-wheels`).
>
> 4. Activate your new virtual environment:
>
>    > ### {% icon code-in %} Input: Bash
>    > ```bash
>    > . .venv/bin/activate
>    > ```
>    {: .code-in}
>
>    Once activated, you'll see the name of the virtual environment prepended to your shell prompt: `(.venv)$`.
>
> 5. Finally, let's create a new branch for your edits:
>
>    > ### {% icon code-in %} Input: Bash
>    > ```bash
>    > git checkout -b my-training
>    > ```
>    {: .code-in}
>
>    Now when you run `git branch` you'll see that your new branch is activated:
>
>    > > ### {% icon code-in %} Input: Bash
>    > > ```bash
>    > > git branch
>    > > ```
>    > {: .code-in}
>    >
>    > > ### {% icon code-out %} Output
>    > > ```bash
>    > >   training
>    > > * my-training
>    > > ```
>    > {: .code-out}
>    {: .code-2col}
>
>    Note: `my-training` is just an example; you can call your new branch anything you like.
{: .hands_on}




## Unit tests

There are numerous unit tests in the Galaxy code code base. However, existing unit tests do not cover the entire code, nor do they verify all the relevant logic. In this tutorial, we will be writing unit tests for such code.

#### Testing simple functions

We'll start by looking at the module `lib/galaxy/util/bytesize.py`. This module is not covered by unit tests. There are 2 obvious targets for unit tests: the `parse bytesize` function and the `to_unit` method of the `ByteSize` class. Our goal is to come up with a reasonable selection of unit tests that verify essential aspects of this functionality.

We'll start with the `parse_bytesize` function. Look at the code of the function and identify the logic that, you think, needs testing. Essentially, you want to ensure that any combination of valid input produces the expected output. Then add a few tests. (see https://docs.pytest.org/en/7.1.x/how-to/assert.html for examples)

{% include topics/dev/tutorials/writing_tests/unit1.md %}

Run the tests to verify they pass. You can use the `pytest` command:

> ### {% icon code-in %} Input: Bash
> ```bash
> pytest test/unit/util/test_bytesize.py
> ```
{: .code-in}

#### Verifying that an error is raised

Now let's add a test that verifies that invalid input raises an error, as expected.

{% include topics/dev/tutorials/writing_tests/unit2.md %}

Run the tests to verify they pass.

#### Using fixtures to reduce code duplication

Let's move on to the `to_unit` method of the `ByteSize`class. Just like you did for the previous function, add a few tests that verify that valid input produces expected output, and invalid input raises an error.

{% include topics/dev/tutorials/writing_tests/unit3.md %}

Run the tests to verify they pass.

Our tests are looking good for the most part, although code duplication is starting to creep in. It's time to refactor. 

Let's start by factoring out the ByteSize object creation into a fixture (see https://docs.pytest.org/en/7.1.x/how-to/fixtures.html). In this particular case moving this code into a fixture might be unnecessary, however, it provides an example of an approach that is very useful in more complex scenarios and is heavily used in Galaxy testing code.

{% include topics/dev/tutorials/writing_tests/unit4.md %}

Run the tests to verify they pass.

#### Parametrization of test functions

Finally, our `test_bytesize_to_unit` test has a lot of assert statements of the same form. We can do better! Let's use pytest's test parametrization feature (see https://docs.pytest.org/en/7.1.x/how-to/parametrize.html) to eliminate this redundancy. Again, as with the fixture example, this particular case would be fine without this feature. However, sometime you want to run the same test function on hundreds of different input combinations, in which case this feature is invaluable (e.g. Galaxy's tool tests, or integration tests for configuration settings).

{% include topics/dev/tutorials/writing_tests/unit5.md %}

Run the tests to verify they pass.

#### Dealing with less trivial code

So far we've been dealing with simple clean functions that had no external dependencies. However, quite often you will encounter lengthy blocks of code implementing nontrivial logic, that depend on objects that would be hard to instantiate in the context of a test. There are many excellent resources on this topic (the book Working Effectively with Legacy Code by Michael Feathers is but one example). In this tutorial, we'll demonstrate a few basic approaches to handling such cases.

Your target module for the following exercises is `lib/galaxy/security/validate_user_input.py`. The unit tests for this module are located at `test/unit/data/security/test_validate_user_input.py`. You will be adding your code to this testing module.

If you look at the tests for this module, you'll note that most of its functionality is directly or indirectly covered by tests. However, the `validate_email` function does not have tests - and that's a shame, for despite being relatively small, it does contain logic that is not immediately obvious (consider how the function behaves under different combinations of the allowlist and blocklist). That's the function we'll be testing.

The `validate_email` function is not fun to test. It depends on a transaction and a user object - one does not simply ~~walk into Mordor~~ instantiate those in the context of a unit test. It has multiple paths of execution, to determine which it accesses the app object (which is a reference to a running instance of Galaxy) to access configuration values, the User object to check the value of its email attribute, and, of course, it calls the database. None of these dependencies require testing, at least not in the context of a unit test. However, we do want to verify the function's core logic: the cases when no validation is required, when the provided email already exists, and the different combinations of the allowlist and blocklist.

#### Mocking a dependency

First, let's test the case when the function returns early - when the email argument is the empty string or when its value is the same as the value of the `email` attribute of the provided `user` object. The first case is trivial, but the next one requires a user object. However, if we think about it, we only care about that user having a given string as its `email` attribute. So let's use a "test double" - a stand-in empty object, and give it an `email` attribute. We'll call it `MockUser` (strictly speaking, it's a stub, and "mocks aren't stubs"](https://martinfowler.com/articles/mocksArentStubs.html), but we'll follow the naming convention often used in Galaxy). We can pass `None` as the transaction object in both cases.

{% include topics/dev/tutorials/writing_tests/unit6.md %}

Now run the tests for this module:

> ### {% icon code-in %} Input: Bash
> ```bash
> pytest test/unit/data/security/test_validate_user_input.py
> ```
{: .code-in}

As a word of caution, mocking (or any kind of patching to accommodate testing) may lead to brittle tests: by relying too much on *how* the logic is implemented, we'll be forced to adjust the test each time that implementation changes. It's a tradeoff between having narrowly-scoped tests that will point to the exact location of the problem but may lock the code under test into a given state, and higher-level integration tests that test the end result without "looking under the hood", but may be less helpful when pinpointing the cause of a failed test. So, use sparingly and keep it simple.

#### Refactoring for testability

Next, we'd like to test the case when the provided email already exists, or doesn't exist. To determine the existence of an email the function calls the database - we don't want to do that. Instead, we can simulate both conditions, like we simulated the user email in the previous exercise.

For this, we need to factor out the code that calls the database into its own function - then we can override that function for our test. So let's modify `lib/galaxy/security/validate_user_input.py`.

{% include topics/dev/tutorials/writing_tests/unit7.md %}

#### "Monkeypatching" the module under test

Now we can use pytest's monkeypatch built-in fixture to replace that function with a stub method and then use it to return `True` or `False` to test both cases. (see https://docs.pytest.org/en/7.1.x/how-to/monkeypatch.html?highlight=monkeypatch for more details).

{% include topics/dev/tutorials/writing_tests/unit8.md %}

#### More refactoring

Run the tests to verify they pass after our latest modifications. But one test fails! The problem is that our validation function is not yet ready for testing: it still relies on the transaction object to access the running instance of Galaxy to retrieve 2 configuration settings: `email_domain_allowlist_content` and `email_domain_blocklist_content`.

How would you fix this problem?

{% include topics/dev/tutorials/writing_tests/unit9.md %}

#### More "monkeypatching"

Run the tests. They still fail. That's because we haven't replaced the new functions with a mock for our test. Try to do that - it's the same approach as we used for monkeypatching the `check_for_existing_email` function. Note that both lists should be `None` to be ignored.

{% include topics/dev/tutorials/writing_tests/unit10.md %}

Run the tests - they should pass now.

#### Utilizing the mocks to test the logic

Now that our tests pass, it's time to test the last bit of logic: the allowlist and the blocklist. Our first step is to add a simple test to verify that if the allow-list is not empty, only emails with allowed domain names will be validated.

{% include topics/dev/tutorials/writing_tests/unit11.md %}

Run the tests to verify we haven't broken anything.

However, again, we are getting a lot of duplication from our monkeypatching code, and there are more tests to come. Time to move some of it into fixtures; **leave only the code that is specific to the test**.

{% include topics/dev/tutorials/writing_tests/unit12.md %}

Again, run the tests to verify we haven't broken anything.

#### Reformatting for improved readability

Before we add more tests, let's reformat our code to group all the tests that target the `validate_email` function in one class. Keep in mind that each test gets its own instance of the class - so you cannot share instance state across tests. However, grouping related tests makes the testing module easier to navigate (there are more benefits; see https://docs.pytest.org/en/7.1.x/getting-started.html#group-multiple-tests-in-a-class).

{% include topics/dev/tutorials/writing_tests/unit13.md %}

Again, run the tests to verify we haven't broken anything.

#### Adding the remaining tests

We have all the infrastructure in place. We need to add 2 last tests. First, you may notice that if the allowlist is not empty, the blocklist is ignored. Let's write a test to verify this.

{% include topics/dev/tutorials/writing_tests/unit14.md %}

Run the test - they should pass.

Finally, let's test the blocklist.

{% include topics/dev/tutorials/writing_tests/unit15.md %}

Run the test one last time - they should pass.

And we are done! Here's the final version of the code we added to this testing module.

{% include topics/dev/tutorials/writing_tests/unit16.md %}

## API tests

TODO

## End-to-end tests

TODO

## Conclusion

TODO
