---
layout: tutorial_hands_on

title: "Data source integration"
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

# Contributed to Galaxy Core 

This tutorial walks you through an extension to Galaxy and how to contribute back to the core project.

To setup the proposed extension imagine you're running a specialized Galaxy server and each of your users only use a few of Galaxy datatypes. You'd like to tailor the UI experience by allowing user's of Galaxy to select their favorite extensions for additional filtering in downstream applications, UI extensions, etc..

Like many extensions to Galaxy, the proposed change requires persistent state. Galaxy stores most persistent state in its relational database. The Python layer that defines Galaxy data model is setup by defining SQL Alchemy models.

The proposed extension could be implemented several different ways on Galaxy's backend and we will choose one for this example for its simplicity not for its correctness or cleverness, because our purpose here to demonstrate modifying and extending various layers of Galaxy.

With simplicity in mind, we will implement our proposed extension to Galaxy by adding a single new table to Galaxy's data model called ``user_favorite_extension``. The concept of a favorite extension will be represented by a one-to-many relationship from the table that stores Galaxy's user records to this new table. The extension itself that will be favorited will be stored as a ``Text`` field in this new table. This table will also need to include a integer primary key named ``id`` to follow the example set by the rest of the Galaxy data model.

The relational database tables consumed by Galaxy are defined in ``lib/galaxy/model/mapping.py``.

> ### {% icon question %} Questions
>
> 1. What should the SQL Alchemy model named corresponding to the table ``user_favorite_extension`` based on other examples in the file.
> 2. What table stores Galaxy's user records?
> 2. What is another simple table with a relationship with the Galaxy's user table?
>
> > ### {% icon solution %} Solution
> > 1. ``UserFavoriteExtension`` 
> > 2. ``galaxy_user``
> > 3. An example table might be the ``user_preference`` table.
> >  {: .solution }
{: .question}


> ### {% icon hands_on %} Hands-on: Implement the SQL Alchemy Mapping
>    > Implement the required changes ``mapping.py`` to add a mapping for
>    > the proposed ``user_favorite_extension`` table.
>    >
>    > ### {% icon solution %} Solution
>    >
>    > {% include dev/tutorials/core-contributing/mapping.py_diff.md %}
>    >
>    > {: .solution }
{: .hands_on}


The Python model objects used by Galaxy corresponding to these tables are
defined in ``lib/galaxy/model/__init__.py``.

> ### {% icon hands_on %} Hands-on: Implement a Model Object
>    > Modify ``lib/galaxy/model/__init__.py`` to add a model class
>    > called ``UserFavoriteExtension`` as described above.
>    >
>    > ### {% icon solution %} Solution
>    >
>    > {% include dev/tutorials/core-contributing/__init__.py_diff.md %}
>    >
>    > {: .solution }
{: .hands_on}

There is one last database issue to consider before moving on to considering the API.
Each successive release of Galaxy requires recipes for how to migrate old database schemes
to updated ones. In Galaxy these a
