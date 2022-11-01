---
layout: tutorial_hands_on

title: "Contributing a New Feature to Galaxy Core"
questions:
 - How can I add a new feature to Galaxy that involves modifications to the database, the API, and the UI?
objectives:
 - Learn to develop extensions to the Galaxy data model
 - Learn to implement new API functionality within Galaxy
 - Learn to extend the Galaxy user interface with VueJS components
time_estimation: "3H"
contributors:
 - jmchilton
 - jdavcs
key_points:
 - Galaxy database interactions are mitigated via SQLAlchemy code in lib/galaxy/model.
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

# Introduction


This tutorial walks you through developing an extension to Galaxy, and how to contribute back to the core project.

To setup the proposed extension imagine you're running a specialized Galaxy server and each of your users only use a few of Galaxy datatypes. You'd like to tailor the UI experience by allowing users of Galaxy to select their favorite extensions for additional filtering in downstream applications, UI extensions, etc..

Like many extensions to Galaxy, the proposed change requires persistent state. Galaxy stores most persistent state in a relational database. The Python layer that defines Galaxy's data model is setup by defining SQLAlchemy models.

The proposed extension could be implemented in several different ways on Galaxy's backend. We will choose one for this example for its simplicity, not for its correctness or cleverness, because our purpose here is to demonstrate modifying and extending various layers of Galaxy.

With simplicity in mind, we will implement our proposed extension to Galaxy by adding a single new table to Galaxy's data model called ``user_favorite_extension``. The concept of a favorite extension will be represented by a one-to-many relationship from the table that stores Galaxy's user records to this new table. The extension itself that will be favorited will be stored as a ``Text`` field in this new table. This table will also need to include an integer primary key named ``id`` to follow the example set by the rest of the Galaxy data model.

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Forking Galaxy

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
>    > git checkout -b my-feature
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
>    > >   dev
>    > > * my-feature
>    > > ```
>    > {: .code-out}
>    {: .code-2col}
>
>    Note: `my-feature` is just an example; you can call your new branch anything you like.
>
> 6. As one last step, you need to initialize your database. This only applies if you are working on
>    a clean clone and have not started Galaxy (starting Galaxy will initialize the database). 
>    Initializing the database is necessary because you will be making changes to the database
>    schema, which cannot be applied to a database that has not been initialized. 
>
>    To initialize the database, you can either start Galaxy (might take some time when executing
>    for the first time):
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > sh run.sh
>    > ```
>    {: .code-in}
>
>    or you may run the following script (faster):
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > sh create_db.sh
>    > ```
>    {: .code-in}
>
{: .hands_on}

# Models

Galaxy uses a relational database to persist objects and object relationships. Galaxy's data model represents the object view of this data. To map objects and their relationships onto tables and rows in the database, Galaxy relies on SQLAlchemy, which is a SQL toolkit and object-relational mapper.

The mapping between objects and the database is defined in `lib/galaxy/model/__init__.py` via "declarative mapping", which means that models are defined as Python classes together with the database metadata that describes the database table corresponding to each class. For example, the definition of the JobParameter class includes the database table name:

```python
__tablename__ = "job_parameter"
```

and four `Column` attributes that correspond to table columns with the same names:

```python
id = Column(Integer, primary_key=True)
job_id = Column(Integer, ForeignKey("job.id"), index=True)
name = Column(String(255))
value = Column(TEXT)
```

Associations between objects are usually defined with the `relationship` construct. For example, the `UserAddress` model has an association with the `User` model and is defined with the `relationship` construct as the `user` attribute.

> <question-title>about Mapping</question-title>
>
> 1. What should be the SQLAlchemy model named corresponding to the table ``user_favorite_extension`` based on [other examples](https://github.com/galaxyproject/galaxy/blob/dev/lib/galaxy/model/__init__.py)?
> 2. What table stores Galaxy's user records?
> 3. What is another simple table with a relationship with the Galaxy's user table?
>
> > <solution-title></solution-title>
> > 1. ``UserFavoriteExtension``
> > 2. ``galaxy_user``
> > 3. An example table might be the ``user_preference`` table.
> {: .solution }
{: .question }

To implement the required changes to add the new model, you need to create a new class in `lib/galaxy/model/__init__.py` with appropriate database metadata:

1. Add a class definition for your new class
2. Your class should be a subclass of `Base`
3. Add a `__tablename__` attribute
4. Add a `Column` attribute for the primary key (should be named `id`)
5. Add a `Column` attribute to store the extension
6. Add a `Column` attribute that will serve as a foreign key to the `User` model
7. Use the `relationship` function to define an association between your new model and the `User` model.

To define a regular `Column` attribute you must include the datatype:

```python
foo = Column(Integer)
```

To define an attribute that is the primary key, you need to include the `primary_key` argument:
(a primary key is one or more columns that uniquely identify a row in a database table)

```python
id = Column(Integer, primary_key=True)
```

To define an attribute that is a foreign key, you need to reference the associated table + its primary key column using the `ForeignKey` construct (the datatype will be derived from that column, so you don't have to include it):

```python
bar_id = Column(ForeignKey("bar.id"))
```

To define a relationship between tables, you need to set it on *both* tables:

```python
class Order(Base):
    …
    items = relationship("Item", back_populates="order")

class Item(Base):
    …
    order_id = Column(ForeignKey("order.id")
    awesome_order = relationship("Order", back_populates="items")
    # relationship not named "order" to avoid confusion: this is NOT the table name for Order
```

Now modify `lib/galaxy/model/__init__.py` to add a model class called `UserFavoriteExtension` as described above.

{% include topics/dev/tutorials/core-contributing/__init__.py_diff.md %}

# Migrations

There is one last database issue to consider before moving on to considering the API. Each successive release of Galaxy requires recipes for how to migrate old database schemas to updated ones. These recipes are called versions, or revisions, and are implemented using Alembic. 

Galaxy's data model is split into the galaxy model and the install model. These models are persisted in one combined database or two separate databases and are represented by two migration branches: "gxy" (the galaxy branch) and "tsi" (the tool shed install branch). Schema changes for these branches are defined in these revision modules:
- `lib/galaxy/model/migrations/alembic/versions_gxy` (galaxy model)
- `lib/galaxy/model/migrations/alembic/versions_tsi` (install model)

We encourage you to read [Galaxy's documentation on migrations](https://github.com/galaxyproject/galaxy/tree/dev/lib/galaxy/model/migrations), as well as relevant [Alembic documentation](https://alembic.sqlalchemy.org/en/latest/tutorial.html#create-a-migration-script).

For this tutorial, you'll need to do the following:
1. Create a revision template
2. Edit the revision template, filling in the body of the upgrade and downgrade functions. 
3. Run the migration.

> <question-title> about generating a revision template</question-title>
>
> What command should you run to generate a revision template?
>
> > <solution-title></solution-title>
> > ```bash
> > sh run_alembic.sh revision --head=gxy@head -m "Add user_favorite_extentions table"
> > ```
> > The title of the revision is an example only.
> {: .solution }
{: .question }

To fill in the revision template, you need to populate the body of the `upgrade` and `downgrade`
functions. The `upgrade` function is executed during a schema upgrade, so it should create your table. 
The `downgrade` function is executed during a schema downgrade, so it should drop your table.

Note that although the table creation command looks similar to the one we used to define the model,
it is not the same. Here, the `Column` definitions are arguments to the `create_table` function.
Also, while you didn't have to specify the datatype of the `user_id` column in the model, you must
do that here.

{% include topics/dev/tutorials/core-contributing/revisionid_add_user_favorite_extensions.py_diff.md %}

> <question-title>about running the migration</question-title>
>
> What command should you run to upgrade your database to include the new table? 
>
> > <solution-title></solution-title>
> > ```bash
> > sh manage_db.sh upgrade
> > ```
> {: .solution }
{: .question }


To verify that the table has been added to your database, you may use the SQLite CLI tool. First,
you login to your database; then you display the schema of the new table; and, finally, you verify
that the database version has been updated (the first record stored in the `alembic_version` table
is the revision identifier that corresponds to the revision identifier in the revision file you
added in a previous step.

{% include topics/dev/tutorials/core-contributing/sqlite_cli.md %}

# Test Driven Development

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

# Run the Tests

Verify this test fails when running stand-alone.

> <code-in-title>Bash</code-in-title>
> ```bash
> ./run_tests.sh -api lib/galaxy_test/api/test_users.py::UsersApiTestCase::test_favorite_extensions
> ```
{: .code-in}

# Implementing the API

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

Ideally, you'd start at the top of the test case - make sure it fails on the first API request,
implement ``get_favorite_extensions`` on the manager and the API code to wire it up, and continue
with ``add_favorite_extension`` before finishing with ``delete_favorite_extension``.

> <code-in-title>Bash</code-in-title>
> ```bash
> ./run_tests.sh -api lib/galaxy_test/api/test_users.py::UsersApiTestCase::test_favorite_extensions
> ```
{: .code-in}

# Building the UI

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
``client/src/entry/analysis/router.js`` respond to the route
added above in ``buildapp.py`` and render the fictitious VueJS component
``FavoriteExtensions``.

{% include topics/dev/tutorials/core-contributing/router.js_diff.md %}

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
