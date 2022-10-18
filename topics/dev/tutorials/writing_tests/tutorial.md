---
layout: tutorial_hands_on

title: "Writing Automated Tests for Galaxy"
objectives:
  - "Learn about the different types of automated tests in Galaxy"
  - "Learn to write API tests"
  - "Learn to write unit tests"
time_estimation: "3h"
key_points:
  - Read the Writing Tests for Galaxy documentation article
  - Check for existing examples of similar tests before implementing your own
  - Prefer Galaxy's abstractions to writing low level code whenever possible
  - Write tests that do not depend on each other or the state of the database
  - When testing the API, verify the return status code before checking the response data
  - Refactor the code under test if that's needed to make it testable
  - Move redundant setup code into fixtures
  - Use test doubles (mocks, stubs, etc.) sparingly
contributors:
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

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Testing in Galaxy

The Galaxy code base contains thousands of tests that include tests of different types (unit vs. functional vs. end-to-end; client vs. backend, etc.) that are supported by a variety of testing frameworks and libraries. In this tutorial, we will offer a small, yet representative sample of the types of tests you might write, as well as the concepts and issues you may need to be familiar with when writing tests for Galaxy code, whether as part of a new feature you are implementing, or as a standalone contribution to Galaxy's testing code.

A good way to start learning about Galaxy's testing infrastructure and how to use it is to read the documentation article on the different types of tests that are present in the code base as well as how to determine which type is most appropriate for a given scenario (see [Writing Tests for Galaxy](https://docs.galaxyproject.org/en/master/dev/writing_tests.html)).

In addition to reading the documentation, it is essential to have a reasonable understanding of Galaxy's code organization. The best documentation resource for that is the [Galaxy Code Architecture slides](https://training.galaxyproject.org/training-material/topics/dev/tutorials/architecture/slides.html#1). In addition, we recommend the [Contributing a New Feature to Galaxy Core](https://training.galaxyproject.org/training-material/topics/dev/tutorials/core-contributing/tutorial.html) tutorial, which, among other things, combines some of the concepts from the architecture slides with the material covered in this tutorial.

The most detailed and up-to-date documentation on running Galaxy tests is located at the top of the `run_tests.sh` script in Galaxy's root directory.

To run Galaxy tests, you may also use the `pytest` command directly (except for client tests, which provide their own infrastructure and scripts described in ``client/README.md``. Using pytest directly is most convenient for running individual tests (although you might have to manually set the required environment variables for all but unit tests, as per documentation in `run_tests.sh`). However, keep in mind that the `run_scripts.sh` script optimizes the tests by reusing the same Galaxy instance and database for API tests, so running API tests in batch with pytest would be very inefficient.

Another useful resource is the [Debugging Galaxy](https://training.galaxyproject.org/training-material/topics/dev/tutorials/debugging/tutorial.html) tutorial which contains a lot of useful information on how to debug test failures that occur both locally and remotely.

Finally, nothing can substitute studying Galaxy's test code - we encourage you to always look for examples of similar tests and testing scenarios before you write your own.

# Local development environment setup

{% snippet topics/dev/faqs/contributing.md %}

> <hands-on-title>Setup your local Galaxy instance</hands-on-title>
>
> 1. Use GitHub UI to fork Galaxy's repository at `galaxyproject/galaxy`.
> 2. Clone your forked repository to a local path, further referred to as `GALAXY_ROOT` and `cd` into `GALAXY_ROOT`. Note that we specify the tutorial branch with the `-b` option:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > git clone https://github.com/<your-username>/galaxy GALAXY_ROOT
>    > cd GALAXY_ROOT
>    > ```
>    {: .code-in}
>
> 3. Before we can use Galaxy, we need to create a virtual environment and install the required dependencies. This is generally done with the `common_startup.sh` script:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > bash scripts/common_startup.sh --dev-wheels
>    > ```
>    {: .code-in}
>
>    Make sure your Python version is at least 3.7 (you can check your Python version with `python --version`). If your system uses an older version, you may specify an alternative Python interpreter using the `GALAXY_PYTHON` environment variable (`GALAXY_PYTHON=/path/to/alt/python bash scripts/common_startup.sh --dev-wheels`).
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
>    > >   training
>    > > * my-training
>    > > ```
>    > {: .code-out}
>    {: .code-2col}
>
>    Note: `my-training` is just an example; you can call your new branch anything you like.
{: .hands_on}

# API tests

We turn to API tests when we need to test some feature that requires a running Galaxy instance and, in many cases, the Galaxy database. These tests use the Galaxy API to drive the test. 

In a way, these tests may be the simplest to write. While the required setup of an API test is, certainly, more involved than that of a basic unit test, Galaxy's testing infrastructure takes care of all the heavy lifting, such as creating a test database, configuring and starting up a Galaxy instance, and tearing down the setup upon test completion. The testing infrastructure also provides a wealth of convenient abstractions that simplify pre-populating the database with the necessary state for each test, interacting with the API, as well as expressing expectations about the outcomes of a test. Thus, writing an API test boils down to calling the appropriate API endpoint and verifying the result.

We are not developing any new features in this tutorial (check out the [Contributing a New Feature to Galaxy Core](https://training.galaxyproject.org/training-material/topics/dev/tutorials/core-contributing/tutorial.html) tutorial). However, we need to exercise the API, so we will be using existing Galaxy functionality as our testing context. Much of that functionality is already covered by tests. Existing tests may include advanced concepts or details on Galaxy's internals that are beyond the scope of a tutorial on testing, so they are not optimal as examples for training. Therefore, the approach we will use is as follows:

1. Pick a controller `foo` in the `lib/galaxy/webapps/galaxy/api` directory.
2. Pick one of the existing endpoints from the selected controller. If the endpoint represents a `GET` request (e.g. `@router.get("/api/foos`), you may be able to view the results on a running instance at `[host]/api/foos`
3. Write a new test to verify the specified functionality.

Ideally, we'd like to follow the process of test-driven development: write a test, run it and watch it fail, implement the functionality under test, rerun the test to verify it passes. Although we are not developing any new functionality here, we still would like to verify that the test we have written is a legitimate test of the correctness of the underlying code. In other words, we want it to be deterministic, and fail if the input is invalid. For this, you may want to intentionally break the test on its first run to verify it doesn't pass under any condition. Then fix the test and rerun. Although this doesn't guarantee correctness, it helps prevent obvious errors in the testing code (such as forgetting to edit the endpoint after cutting and pasting code from a similar test function.

## Basic API test

Let's start by writing a basic test for a very simple API endpoint: `api/version`. The controller for this endpoint is located at `lib/galaxy/webapps/galaxy/api/configuration.py`. You can check the format of the data at `https://usegalaxy.org/api/version`.

First, you need to create a new file at `lib/galaxy_test/api/test_mytutorial.py`. For simplicity, we'll place all new tests in this module. Next, add a class definition for `MyTutorialApiTestCase` that should be a subclass of `ApiTestCase`. Then add a test method, `test_version_is_current`, where you will (1) call the API via the ``_get`` method that returns a response object; and (2) verify that the response contains the "22.09" version number (assuming you have cloned the "dev" branch; otherwise your version may be different).

If your test fails, one way to debug it is to insert a `breakpoint()` statement into the body of the test (right after the call to ``_get`` would be a logical spot), and then use [pdb](https://docs.python.org/3/library/pdb.html), Python's interactive debugger, to explore the response at runtime with the test paused (see [Debugging Galaxy](https://training.galaxyproject.org/training-material/topics/dev/tutorials/debugging/tutorial.html) for more details on using pdb to debug Galaxy).

{% include topics/dev/tutorials/writing_tests/api1.md %}

You may notice that in the provided solution we also call ``response.raise_for_status()``: this is a requests library function that verifies that the request was successful (of course, we can also explicitly check the status code). It is helpful to include this check, because otherwise, if an error is raised, the error message in the error log will be not as helpful. Try using a nonextistant endpoint and running the test with and without the `raise_for_status` call. Compare these error messages:

With the status check:
- `requests.exceptions.HTTPError: 404 Client Error: Not Found for url: http://127.0.0.1:9584/api/versionx`

Without the status check:
- `requests.exceptions.JSONDecodeError: Expecting value: line 1 column 1 (char 0)`

## Test creating a simple object

Now let's write an API test for adding a simple object. By "simple", we mean that this object does not depend on any other objects, so, for example, we don't need to provide a reference to a History or User object to create it. We will be creating a role, using the `api/roles` `POST` endpoint located at `lib/galaxy/webapps/galaxy/api/roles.py`.

To create a role, we need to call the ``_post`` method, passing it 4 arguments:
- the string value of the endpoint
- the payload containing the data for the new Role object
- a boolean argument setting the format of the payload data: `json=True`
- a boolean argument setting the test user's admin status: `admin=True`

How do we know what data to include in the payload? Look at the method signature in the API:
```
self, trans: ProvidesUserContext = DependsOnTrans, role_definition_model: RoleDefinitionModel = Body(...)
```
`RoleDefinitionModel` is the Pydantic model that describes the structure of the data passed with this argument. Looking at its definition (check the import statement at the top of the file to find its location), we see the following:

```python
class RoleDefinitionModel(BaseModel):
    name: str = RoleNameField
    description: str = RoleDescriptionField
    user_ids: Optional[List[EncodedDatabaseIdField]] = Field(title="User IDs", default=[])
    group_ids: Optional[List[EncodedDatabaseIdField]] = Field(title="Group IDs", default=[])
```

Thus, we need to provide a name and a description, both of which are string values (the other members are optional), formatted as JSON. (For more details on how the API works, see [FastAPI's excellent documentation](https://fastapi.tiangolo.com/tutorial/).)

The response will return the Role object, so testing it is straightforward.

{% include topics/dev/tutorials/writing_tests/api2.md %}

## A better test for creating a simple object

Can we do better? We can verify that the new object has been added by using another endpoint to retrieve that object from the database. Here's a first take:

```python
def test_create_role_WRONG(self):
    ROLE_COUNT = 1
    # verify no roles except one built-in role for test user
    response = self._get("roles")
    response.raise_for_status()
    data = response.json()
    assert len(data) == ROLE_COUNT

    # create role
    name = 'my cool name'
    description = 'description of this cool role'
    payload = {
        "name": name,
        "description": description,
    }
    response = self._post("roles", payload, admin=True, json=True)
    response.raise_for_status()
    role = response.json()
    assert role["name"] == name
    assert role["description"] == description

    # verify role has been added
    response = self._get("roles")
    response.raise_for_status()
    data = response.json()
    assert len(data) == ROLE_COUNT + 1
```

Everything is straightforward: we verify that there is only n=1 role currently present (a built-in for the test user), we add the role, then retrieve all roles and verify that the new total is n + 1. Unfortunately, ***this is the wrong approach***. Try running this together with the previous version of the test - you'll get an assertion error: there's an extra role in the database! 

```
FAILED lib/galaxy_test/api/test_mytutorial.py::MyTutorialApiTestCase::test_create_role_WRONG - AssertionError: assert 2 == 1
```

The reason for that is that we are using the same database for both tests, so whatever artifacts are created in one test will affect the following test if we make any kind of assumptions about the state of the database.

How do we rewrite this test without making assumptions about database state? We can't retrieve the roles and use their initial quantity as our baseline: that's a perfect recipe for a race condition (if some other test creates a role after we have counted the existing roles but before we have created our new role, our test will fail). Instead, we can use the role's name (which is unique across all roles) to verify that our object has been added. 

We could generate our own unique name, but we can also simply use an existing method used across Galaxy's tests: ``get_random_name()`` (we need to add an import and override the ``setUp`` method for that: see the solution code). Now we have a test we can safely run together with our old test (or with any number of other tests that create role objects).

{% include topics/dev/tutorials/writing_tests/api3.md %}

Can't we destroy and create or re-populate the database for each test? We can, and if you absolutely need to do that, there is infrastructure for that in `test/unit/data/model/testing_utils` - you can see examples of its usage in the mapping and migrations unit tests (check these subdirectories). However, setting up and initializing the database is an expensive operation, which will be noticeable for a single test; Galaxy has hundreds of tests that rely on the database, so providing a fresh copy of the database for each test function is infeasible. 

## Use Galaxy's populators for setting up database state

As our last example, we'll look at how to take advantage of some of the many abstractions provided by Galaxy's testing infrastructure. Supposing you'd like to test the "datasets-filter-by-history" functionality exposed via the ``api/datasets`` endpoint (see `lib/galaxy/webapps/galaxy/api/datasets.py`). The test design is straightforward: 
- Create 2 history objects h1 and h2.
- Create n and m datasets associated with h1 and h2 respectively.
- Retrieve datasets filtering by h1. Verify there are n results.
- Retrieve datasets filtering by h2. Verify there are m results.

The question is, how do we create the datasets and the histories (and all their respective dependencies?). One way would be to use the API for each such object. However, that would result in a barely readable test function which would mostly consist of setup code. Instead, we want our tests to be simple and understandable at a glance, and contain as little infrastructure code as possible. Enter dataset populators! This is one of many abstractions available to you that make it much easier to write tests that otherwise would have required nontrivial setup.

{% include topics/dev/tutorials/writing_tests/api4.md %}

Populators are used extensively throughout API tests as well as integration and Selenium tests to both setup database state, as well as access information from the Galaxy server. For more details, see the populators module at ``lib/galaxy_test/base/populators.py``.

And we are done! Here's the final version of the code we added to this testing module.

{% include topics/dev/tutorials/writing_tests/api5.md %}

# Unit tests

Unit tests are best suited for testing well-defined, isolated functionality and, with rare exceptions, should not require a running Galaxy instance, a database, or any other external infrastructure. Unit tests are located in the `test/unit` directory.

There are numerous unit tests in the Galaxy code base. However, existing unit tests do not cover the entire code base and do not verify all the relevant logic (so there is plenty of space for improvement, and new contributions are always welcome!) In this tutorial, we will be writing unit tests for such code.

## Testing simple functions

We'll start by looking at the module `lib/galaxy/util/bytesize.py`. This module is not covered by unit tests. There are 2 obvious targets for unit tests: the `parse bytesize` function and the `to_unit` method of the `ByteSize` class. Our goal is to come up with a reasonable selection of unit tests that verify essential aspects of this functionality.

We'll start with the `parse_bytesize` function. Look at the code of the function and identify the logic that, you think, needs testing. Essentially, you want to ensure that any combination of valid input produces the expected output. Next, add a few tests. (see [https://docs.pytest.org/en/7.1.x/how-to/assert.html](https://docs.pytest.org/en/7.1.x/how-to/assert.html) for examples)

{% include topics/dev/tutorials/writing_tests/unit1.md %}

Run the tests to verify they pass. You can use the `pytest` command:

> <code-in-title>Bash</code-in-title>
> ```bash
> pytest test/unit/util/test_bytesize.py
> ```
{: .code-in}

## Verifying that an error is raised

Now let's add a test that verifies that invalid input raises an error, as expected.

{% include topics/dev/tutorials/writing_tests/unit2.md %}

Run the tests to verify they pass.

## Using fixtures to reduce code duplication

Let's move on to the `to_unit` method of the `ByteSize`class. Just like you did for the previous function, add a few tests that verify that valid input produces expected output, and invalid input raises an error.

{% include topics/dev/tutorials/writing_tests/unit3.md %}

Run the tests to verify they pass.

Our tests are looking good for the most part, although code duplication is starting to creep in. It's time to refactor! 

Let's start by factoring out the ByteSize object creation into a fixture (see [https://docs.pytest.org/en/7.1.x/how-to/fixtures.html](https://docs.pytest.org/en/7.1.x/how-to/fixtures.html)). In this particular case moving this code into a fixture might be unnecessary, however, it provides an example of an approach that is very useful in more complex scenarios and is heavily used in Galaxy's testing code.

{% include topics/dev/tutorials/writing_tests/unit4.md %}

Run the tests to verify they pass.

## Parametrization of test functions

Finally, our `test_bytesize_to_unit` test has a lot of assert statements of the same form. We can do better! Let's use pytest's test parametrization feature (see [https://docs.pytest.org/en/7.1.x/how-to/parametrize.html](https://docs.pytest.org/en/7.1.x/how-to/parametrize.html)) to eliminate this redundancy. Again, as with the fixture example, this particular case would be fine without this feature. However, sometime you want to run the same test function on hundreds of different input combinations, in which case this feature is invaluable (e.g. Galaxy's tool tests, or integration tests for configuration settings).

{% include topics/dev/tutorials/writing_tests/unit5.md %}

Run the tests to verify they pass.

## Dealing with less trivial code

So far we've been dealing with simple clean functions that had no external dependencies. However, quite often you will encounter lengthy blocks of code implementing nontrivial logic, that depend on objects that would be hard to instantiate in the context of a test. There are many excellent resources on this topic (the book [Working Effectively with Legacy Code](https://www.oreilly.com/library/view/working-effectively-with/0131177052/) by Michael Feathers is but one example). In this tutorial, we'll demonstrate a few basic approaches to handling such cases.

Your target module for the following exercises is `lib/galaxy/security/validate_user_input.py`. The unit tests for this module are located at `test/unit/data/security/test_validate_user_input.py`. You will be adding your code to this testing module.

If you look at the tests for this module, you'll note that most of its functionality is directly or indirectly covered by tests. However, the `validate_email` function does not have tests - and that's a shame for, despite being relatively small, it contains logic that is not immediately obvious (consider how the function behaves under different combinations of the allowlist and blocklist). That's the function we'll be testing.

The `validate_email` function is not fun to test. It depends on a Transaction and User objects - one does not simply ~~walk into Mordor~~ instantiate those in the context of a unit test. It has multiple paths of execution, to determine which it accesses the GalaxyApplication object (which is a reference to a running instance of Galaxy) to access configuration values, the User object to check the value of its email attribute, and in addition to all that, of course, it calls the database. None of these dependencies require testing, at least not in the context of a unit test. However, we do want to verify the function's core logic: the cases when no validation is required, when the provided email already exists, and the different combinations of the allowlist and blocklist.

## Mocking a dependency

First, let's test the case when the function returns early - when the email argument is the empty string or when its value is the same as the value of the `email` attribute of the provided User object. The first case is trivial, but the next one requires a User object. However, if we think about it, we only care about that user having a given string as its `email` attribute. So let's use a "test double" - a stand-in empty object, and give it an `email` attribute. We'll call it `MockUser` (strictly speaking, it's a stub, and ["mocks aren't stubs"](https://martinfowler.com/articles/mocksArentStubs.html), but we'll follow the naming convention often used in Galaxy). We can pass `None` as the transaction object in both cases.

{% include topics/dev/tutorials/writing_tests/unit6.md %}

Now run the tests for this module:

> <code-in-title>Bash</code-in-title>
> ```bash
> pytest test/unit/data/security/test_validate_user_input.py
> ```
{: .code-in}

> <warning-title>Warning</warning-title>
> As a word of caution, mocking (or any kind of patching to accommodate testing) may lead to brittle tests: by relying too much on *how* the logic is implemented, we'll be forced to adjust the test each time that implementation changes. It's a tradeoff between having narrowly-scoped tests that will point to the exact location of the problem but may lock the code under test into a given state, and higher-level integration-type tests that test the end result without "looking under the hood", but may be less helpful when pinpointing the exact cause of a failed test. So, use sparingly and keep it simple.
{: .warning}

## Refactoring for testability

Next, we'd like to test the case when the provided email already exists, or doesn't exist. To determine the existence of an email the function calls the database - we don't want to do that. Instead, we can simulate both conditions, like we simulated the user email in the previous exercise.

For this, we need to factor out the code that calls the database into its own function - then we can override that function for our test. So let's modify `lib/galaxy/security/validate_user_input.py`.

{% include topics/dev/tutorials/writing_tests/unit7.md %}

## "Monkeypatching" the module under test

Now we can use pytest's "monkeypatch" built-in fixture to replace that function with a stub method and then use it to return `True` or `False` to test both cases. (see [https://docs.pytest.org/en/7.1.x/how-to/monkeypatch.html?highlight=monkeypatch](https://docs.pytest.org/en/7.1.x/how-to/monkeypatch.html?highlight=monkeypatch) for more details).

{% include topics/dev/tutorials/writing_tests/unit8.md %}

## More refactoring

Run the tests to verify they pass after our latest modifications. But one test fails! The problem is that our validation function is not yet ready for testing: it still relies on the transaction object to access the running instance of Galaxy to retrieve 2 configuration settings: `email_domain_allowlist_content` and `email_domain_blocklist_content`.

How would you fix this problem?

{% include topics/dev/tutorials/writing_tests/unit9.md %}

## More "monkeypatching"

Run the tests. They still fail. That's because we haven't replaced the new functions with a mock for our test. Try to do that - it's the same approach as we used for monkeypatching the `check_for_existing_email` function. Note that both lists should be `None` to be ignored.

{% include topics/dev/tutorials/writing_tests/unit10.md %}

Run the tests - they should pass now.

## Utilizing the mocks to test the logic

Now that our tests pass, it's time to test the last bit of logic: the allowlist and the blocklist. Our first step is to add a simple test to verify that if the allow-list is not empty, only emails with allowed domain names will be validated.

{% include topics/dev/tutorials/writing_tests/unit11.md %}

Run the tests to verify we haven't broken anything.

However, again, we are getting a lot of duplication from our monkeypatching code, and there are more tests to come. Time to move some of it into fixtures; **leave only the code that is specific to the test**.

{% include topics/dev/tutorials/writing_tests/unit12.md %}

Again, run the tests to verify we haven't broken anything.

## Reformatting for improved readability

Before we add more tests, let's reformat our code to group all the tests that target the `validate_email` function in one class. Keep in mind that each test gets its own instance of the class - so you cannot share instance state across tests. However, grouping related tests makes the testing module easier to navigate (there are more benefits; see [https://docs.pytest.org/en/7.1.x/getting-started.html#group-multiple-tests-in-a-class](https://docs.pytest.org/en/7.1.x/getting-started.html#group-multiple-tests-in-a-class)).

{% include topics/dev/tutorials/writing_tests/unit13.md %}

Again, run the tests to verify we haven't broken anything.

## Adding the remaining tests

Now we have all the infrastructure in place. We need to add 2 more tests. First, you may notice that if the allowlist is not empty, the blocklist is ignored. Let's write a test to verify this.

{% include topics/dev/tutorials/writing_tests/unit14.md %}

Run the tests - they should pass.

Finally, let's test the blocklist.

{% include topics/dev/tutorials/writing_tests/unit15.md %}

Run the tests one last time - they should pass.

And we are done! Here's the final version of the code we've added to this testing module.

{% include topics/dev/tutorials/writing_tests/unit16.md %}
