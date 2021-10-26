---
layout: tutorial_hands_on

title: "Contributing a New Feature to Galaxy Core"
questions:
 - How can I develop extensions to Galaxy data model?
 - How can I implement new API functionality within Galaxy?
 - How can I extend the Galaxy user interface with VueJS components?
objectives:
time_estimation: "3H"
contributors:
 - jmchilton
key_points:
 - Galaxy database interactions are mitigated via SQL Alchemy code in lib/galaxy/model.
 - Galaxy API endpoints are implemented in lib/galaxy/webapps/galaxy, but generally defer to application logic in lib/galaxy/managers.
 - Galaxy client code should do its best to separate API interaction logic from display components.
subtopic: core
requirements:
  -
    type: "internal"
    topic_name: contributing
    tutorials:
        - github-command-line-contribution
  -
    type: "internal"
    topic_name: dev
    tutorials:
        - architecture
---

# Contributing a New Feature to Galaxy Core

This tutorial walks you through an extension to Galaxy and how to contribute back to the core project.

To setup the proposed extension imagine you're running a specialized Galaxy server and each of your users only use a few of Galaxy datatypes. You'd like to tailor the UI experience by allowing users of Galaxy to select their favorite extensions for additional filtering in downstream applications, UI extensions, etc..

Like many extensions to Galaxy, the proposed change requires persistent state. Galaxy stores most persistent state in a relational database. The Python layer that defines Galaxy data model is setup by defining SQLAlchemy models.

The proposed extension could be implemented in several different ways on Galaxy's backend. We will choose one for this example for its simplicity, not for its correctness or cleverness, because our purpose here is to demonstrate modifying and extending various layers of Galaxy.

With simplicity in mind, we will implement our proposed extension to Galaxy by adding a single new table to Galaxy's data model called ``user_favorite_extension``. The concept of a favorite extension will be represented by a one-to-many relationship from the table that stores Galaxy's user records to this new table. The extension itself that will be favorited will be stored as a ``Text`` field in this new table. This table will also need to include an integer primary key named ``id`` to follow the example set by the rest of the Galaxy data model.

## Forking Galaxy

{% snippet topics/dev/faqs/contributing.md %}

> ### {% icon hands_on %} Hands-on: Setup your local Galaxy instance
>
> 1. Use GitHub UI to fork Galaxy's repository at `galaxyproject/galaxy`.
> 2. Clone your forked repository to a local path, further referred to as `GALAXY_ROOT` and `cd` into `GALAXY_ROOT`. Note that we specify the tutorial branch with the `-b` option:
>
>    > ### {% icon code-in %} Input: Bash
>    > ```bash
>    > git clone https://github.com/<your-username>/galaxy GALAXY_ROOT
>    > cd GALAXY_ROOT
>    > ```
>    {: .code-in}
>
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
>    > git checkout -b my-feature
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
>    > >   dev
>    > > * my-feature
>    > > ```
>    > {: .code-out}
>    {: .code-2col}
>
>    Note: `my-feature` is just an example; you can call your new branch anything you like.
{: .hands_on}

## Models

The relational database tables consumed by Galaxy are defined in ``lib/galaxy/model/mapping.py``.

> ### {% icon question %} Questions about Mapping
>
> 1. What should be the SQLAlchemy model named corresponding to the table ``user_favorite_extension`` based on [other examples in the mapping file](https://github.com/galaxyproject/galaxy/blob/dev/lib/galaxy/model/mapping.py)?
> 2. What table stores Galaxy's user records?
> 3. What is another simple table with a relationship with the Galaxy's user table?
>
> > ### {% icon solution %} Solution
> > 1. ``UserFavoriteExtension``
> > 2. ``galaxy_user``
> > 3. An example table might be the ``user_preference`` table.
> {: .solution }
{: .question }

Implement the required changes in ``mapping.py`` to add a mapping for
the proposed ``user_favorite_extension`` table.

{% include topics/dev/tutorials/core-contributing/mapping.py_diff.md %}

The Python model objects used by Galaxy corresponding to these tables are
defined in ``lib/galaxy/model/__init__.py``.

Modify ``lib/galaxy/model/__init__.py`` to add a model class
called ``UserFavoriteExtension`` as described above.

{% include topics/dev/tutorials/core-contributing/__init__.py_diff.md %}

## Migrations

There is one last database issue to consider before moving on to considering the API.
Each successive release of Galaxy requires recipes for how to migrate old database schemes
to updated ones. These recipes are called versions and are currently implemented using SQLAlchemy
 Migrate. These versions are stored in ``lib/galaxy/model/migrate``.

Each of these versions is prefixed with a 4-digit number (e.g. ``0125``) to specify the linear
order these migrations should be applied in.

We've developed a lot of abstractions around the underlying migration library we use, so it
is probably better to find existing examples inside of Galaxy for writing migrations than
consulting the SQLAlchemy Migrate documentation.

A good example for the table we need to construct for this example is again based on the existing
user preferences concept in Galaxy.

> ### {% icon question %} Questions about Migrations
>
> What existing Galaxy migration added the concept of user preferences to the Galaxy codebase?
>
> > ### {% icon solution %} Solution
> > [``lib/galaxy/model/migrate/versions/0021_user_prefs.py``](https://github.com/galaxyproject/galaxy/blob/dev/lib/galaxy/model/migrate/versions/0021_user_prefs.py)
> {: .solution }
{: .question }

Add a new file to ``lib/galaxy/model/migrate/versions/`` prefixed appropriately.

{% include topics/dev/tutorials/core-contributing/0177_add_user_favorite_extensions.py_diff.md %}

## Test Driven Development

With the database model in place, we need to start adding the rest of
the Python plumbing required to implement this feature. We will do
this with a test-driven approach and start by implementing an API test
that exercises operations we would like to have available for favorite
extensions.

We will stick a test case for user extensions in
``lib/galaxy_test/api/test_users.py`` which is a relatively
straightforward file that contains tests for other user API
endpoints.

Various user-centered operations have endpoints under
``api/user/<user_id>`` and ``api/user/current`` is sometimes
substituable as the current user.

We will keep things very simple and only implement this functionality
for the current user.

We will implement three simple API endpoints.

Method   | Route                                                                  | Definition
------   | -----                                                                  | ---
`GET`    | `<galaxy_root_url>/api/users/current/favorites/extensions`             | This should return a list of favorited extensions for the current user.
`POST`   | `<galaxy_root_url>/api/users/current/favorites/extensions/<extension>` | This should mark an extension as a favorite for the current user.
`DELETE` | `<galaxy_root_url>/api/users/current/favorites/extensions/<extension>` | This should unmark an extension as a favorite for the current user.

Please review ``test_users.py`` and attempt to write a test case that:

- Verifies the test user's initially favorited extensions is an empty list.
- Verifies that a ``POST`` to ``<galaxy_root_url>/api/users/current/favorites/extensions/fasta`` returns a 200 status code indicating success.
- Verifies that after this ``POST`` the list of user favorited extensions contains ``fasta`` and is of size 1.
- Verifies that a ``DELETE`` to ``<galaxy_root_url>/api/users/current/favorites/extensions/fasta`` succeeds.
- Verifies that after this ``DELETE`` the favorited extensions list is again empty.

{% include topics/dev/tutorials/core-contributing/test_users.py_diff.md %}

## Run the Tests

Verify this test fails when running stand-alone.

> ### {% icon code-in %} Input: Bash
> ```bash
> ./run_tests.sh -api lib/galaxy_test/api/test_users.py::UsersApiTestCase::test_favorite_extensions
> ```
{: .code-in}

## Implementing the API

Add a new API implementation file to
``lib/galaxy/webapps/galaxy/api/`` called ``user_favorites.py`` with
an API implementation of the endpoints we just outlined.

To implement the API itself, add three methods to the user manager in
``lib/galaxy/managers/users.py``.

- ``get_favorite_extensions(user)``
- ``add_favorite_extension(user, extension)``
- ``delete_favorite_extension(user, extension)``

{% include topics/dev/tutorials/core-contributing/user_favorites.py_diff.md %}

{% include topics/dev/tutorials/core-contributing/users.py_diff.md %}

This part is relatively challenging and takes time to really become an
expert at - it requires greping around the backend to find similar examples,
lots of trial and error, debugging the test case and the implementation in unison,
etc..

Ideally, you'd start at the top of test case - make sure it fails on the first API request,
implement ``get_favorite_extensions`` on the manager and the API code to wire it up, and continue
with ``add_favorite_extension`` before finishing with ``delete_favorite_extension``.

> ### {% icon code-in %} Input: Bash
> ```bash
> ./run_tests.sh -api lib/galaxy_test/api/test_users.py::UsersApiTestCase::test_favorite_extensions
> ```
{: .code-in}

## Building the UI

Once the API test is done, it is time to build a user interface for this addition to Galaxy. Let's get
some of the plumbing out of the way right away. We'd like to have a URL for viewing the current user's
favorite extensions in the UI. This URL needs to be registered as a client route in ``lib/galaxy/webapps/galaxy/buildapp.py``.

Add ``/user/favorite/extensions`` as a route for the client in ``buildapp.py``.

{% include topics/dev/tutorials/core-contributing/buildapp.py_diff.md %}

Let's add the ability to navigate to this URL and future component to
the "User" menu in the Galaxy masthead. The file ``client/src/layout/menu.js``
contains the "model" data describing the masthead. Add a link to the route
we described above to ``menu.js`` with an entry titled "Favorite Extensions".

{% include topics/dev/tutorials/core-contributing/menu.js_diff.md %}

The next piece of this plumbing is to respond to this route in the
analysis router. The analysis router maps URLs to UI components to
render.

Assume a VueJS component will be available called
``FavoriteExtensions`` in the file
``components/User/FavoriteExtensions/index.js``.  In
``client/src/entry/analysis/AnalysisRouter.js`` respond to the route
added above in ``buildapp.py`` and render the fictious VueJS component
``FavoriteExtensions``.

{% include topics/dev/tutorials/core-contributing/AnalysisRouter.js_diff.md %}

There are many ways to perform the next steps, but like the API entry-point lets
start with a test case describing the UI component we want to write. Below is a
Jest unit test for a VueJS component that mocks out some API calls to ``/api/datatypes``
and the API entry points we implemented earlier and renders an editable list of
extensions based on it.

{% include topics/dev/tutorials/core-contributing/List.test.js_diff.md %}

Sketching out this unit test would take a lot of practice, this is a step
that might be best done just by copying the file over. Make sure the
individual components make sense before continuing though.

Next implement a VueJS component in that same directory called
``List.vue`` that fullfills the contract described by the unit test.

{% include topics/dev/tutorials/core-contributing/List.vue_diff.md %}

Finally, we added a level of indirection when we utilized this
component from the analysis router above by importing it from
``index.js``. Let's setup that file and import the component from
``List.vue`` and export as a component called
``FavoriteExtensions``.

{% include topics/dev/tutorials/core-contributing/index.js_diff.md %}
