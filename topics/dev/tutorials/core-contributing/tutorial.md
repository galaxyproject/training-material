
---
layout: tutorial_hands_on

title: "Contributing to Galaxy Core"
questions:
 - How can I develop extensions to Galaxy data model?
 - How can I implement new API functionality within Galaxy?
 - How can I extend the Galaxy user interface with VueJS components?
objectives:
time_estimation: "180M"
contributors:
 - jmchilton
key_points:
 - Galaxy database interactions are mitigated via SQL Alchemy code in lib/galaxy/model.
 - Galaxy API endpoints are implemented in lib/galaxy/webapps/galaxy, but generally defer to application logic in lib/galaxy/managers. 
 - Galaxy client code should do its best to separate API interaction logic from display components.
---

# Contributing to Galaxy Core 

This tutorial walks you through an extension to Galaxy and how to contribute back to the core project.

To setup the proposed extension imagine you're running a specialized Galaxy server and each of your users only use a few of Galaxy datatypes. You'd like to tailor the UI experience by allowing user's of Galaxy to select their favorite extensions for additional filtering in downstream applications, UI extensions, etc..

Like many extensions to Galaxy, the proposed change requires persistent state. Galaxy stores most persistent state in its relational database. The Python layer that defines Galaxy data model is setup by defining SQL Alchemy models.

The proposed extension could be implemented several different ways on Galaxy's backend and we will choose one for this example for its simplicity not for its correctness or cleverness, because our purpose here to demonstrate modifying and extending various layers of Galaxy.

With simplicity in mind, we will implement our proposed extension to Galaxy by adding a single new table to Galaxy's data model called ``user_favorite_extension``. The concept of a favorite extension will be represented by a one-to-many relationship from the table that stores Galaxy's user records to this new table. The extension itself that will be favorited will be stored as a ``Text`` field in this new table. This table will also need to include a integer primary key named ``id`` to follow the example set by the rest of the Galaxy data model.

The relational database tables consumed by Galaxy are defined in ``lib/galaxy/model/mapping.py``.

> ### {% icon question %} Questions about Mapping
>
> 1. What should the SQL Alchemy model named corresponding to the table ``user_favorite_extension`` based on other examples in the file.
> 2. What table stores Galaxy's user records?
> 3. What is another simple table with a relationship with the Galaxy's user table?
>
> > ### {% icon solution %} Solution
> > 1. ``UserFavoriteExtension`` 
> > 2. ``galaxy_user``
> > 3. An example table might be the ``user_preference`` table.
> {: .solution }
{: .question }

Implement the required changes ``mapping.py`` to add a mapping for
the proposed ``user_favorite_extension`` table.

{% include dev/tutorials/core-contributing/mapping.py_diff.md %}

The Python model objects used by Galaxy corresponding to these tables are
defined in ``lib/galaxy/model/__init__.py``.

Modify ``lib/galaxy/model/__init__.py`` to add a model class
called ``UserFavoriteExtension`` as described above.

{% include topics/dev/tutorials/core-contributing/__init__.py_diff.md %}

There is one last database issue to consider before moving on to considering the API.
Each successive release of Galaxy requires recipes for how to migrate old database schemes
to updated ones. These recipes are called versions and currently implemented using SQL
Alchemy migrate. These versions are stored in ``lib/galaxy/model/migrate``.

Each of these versions is prefixed with a number 4 digit (e.g. ``0125``) to specify the linear
order these migrations should be applied in.

We've devloped a lot of abstractions around the underlying migration library we use, so it
is probably better to find existing examples inside of Galaxy for writing migrations than
consulting the SQL Alchemy Migrate documentation.

A good example for the table we need to construct for this example is again based on the existing
user preferences concept in Galaxy.

> ### {% icon question %} Questions about Migrations
>
> 1. What existing Galaxy migration added the concept of user preferences to the Galaxy codebase?
>
> > ### {% icon solution %} Solution
> > 1. ``lib/galaxy/model/migrate/versions/0021_user_prefs.py``
> {: .solution }
{: .question }


Add a new file to ``lib/galaxy/model/migrate/versions/`` prefixed appropriately.

{% include topics/dev/tutorials/core-contributing/0175_add_user_favorite_extensions.py_diff.md %}

With the database model in place, we need to start adding the rest of
the Python plubming required to implement this feature. We will do
this with a test-driven approach and start by implementing an API test
that exercises operations we would like to have available for favorite
extensions.

We will stick a test case for user extensions in
``lib/galaxy_test/api/test_users.py`` which is a relatively
straight-forward file that contains tests for other user API
endpoints.

Various user centered operations have endpoints under
``api/user/<user_id>`` and ``api/user/current`` is sometimes
substituable as the current user.

We will keep things very simple and only implement this functionality
for the current user.

We will implement three very API endpoints.

- ``GET <galaxy_root_url>/api/users/current/favorites/extensions``. This should return a list of favorited extensions for the current user.
- ``POST <galaxy_root_url>/api/users/current/favorites/extensions/<extension>``. This should mark an extension as a favorite for the current user.
- ``DELETE <galaxy_root_url>/api/users/current/favorites/extensions/<extension>``. This should unmark an extension as a favorite for the current user.

Please review ``test_users.py`` and attempt to write a test case that:

- Verifies the test user's initially favorited extensions is an empty list.
- Verifies that a ``POST`` to ``<galaxy_root_url>/api/users/current/favorites/extensions/fasta`` returns a 200 status code indicating success.
- Verifies that after this ``POST`` the list of user favorited extensions contains ``fasta`` and is of size 1.
- Verifies that a ``DELETE`` to ``<galaxy_root_url>/api/users/current/favorites/extensions/fasta`` succeeds.
- Verifies that after this ``DELETE`` the favorited extensions list is again empty.

{% include topics/dev/tutorials/core-contributing/test_users.py_diff.md %}

Verify this test fails when running stand-alone.

```
./run_tests.sh -api lib/galaxy_test/api/test_users.py::UsersApiTestCase::test_favorite_extensions
```

Add a new API implementation file to
``lib/galaxy/webapps/galaxy/api/`` called ``user_favorites.py`` with
an API implementation of the endpoints we just outlined.


get_favorite_extensions, add_favorite_extension, delete_favorite_extension.
