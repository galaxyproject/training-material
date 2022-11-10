---
layout: tutorial_hands_on

title: "Contributing to BioBlend as a developer"
questions:
    - "How to get started making contributions to BioBlend?"
objectives:
    - "Learn the basics behind BioBlend and the Galaxy API from a developer perspective."
    - "Learn how to implement a simple method in BioBlend."
    - "Learn how to run the BioBlend test suite."
requirements:
    -
        title: "Basic command line knowledge"
        type: none
    -
        title: "Familiarity with git and access to a GitHub account"
        type: none
    -
        title: "Familiarity with the Python programming language"
        type: none
    -
        title: "Familiarity with HTTP methods and JavaScript Object Notation (JSON)"
        type: none
time_estimation: "3h"
key_points:
    - "BioBlend is a Python library that provides methods for easy interaction with the Galaxy API."
    - "Implementing BioBlend methods is generally quite an easy process, making it well suited to beginners and a viable stepping stone into Galaxy development."
subtopic: api
contributors:
    - rikeshi
    - simonbray

abbreviations:
    LDA: Library Dataset Association
    HDA: History Dataset Association
    LDCA: Library Dataset Collection Association
    HDCA: History Dataset Collection Association
    API: Application Programming Interface
    PR: Pull Request
---

# Introduction

BioBlend ({% cite Sloggett2013 %}) is a Python library to enable simple interaction with Galaxy ({% cite Afgan2018 %}) via the command line or scripts.

Galaxy is a data analysis platform for accessible, reproducible and transparent computational research. It includes a web interface through which users can design and perform tasks in a visual and interactive manner. The Galaxy server also exposes this functionality through its REST-based {API}.

Computer programs can communicate with the Galaxy server through this {API} and perform similar tasks as can be achieved manually via the web interface. The program sends a network request to a URL of an {API} endpoint. The server then computes the result for this request and sends back a response to the program.

BioBlend provides classes and methods that handle the specific details of this communication for Python programs. Similar libraries also exist for interacting with Galaxy via other programming languages:

* blend4j (Java): [https://github.com/galaxyproject/blend4j](https://github.com/galaxyproject/blend4j)
* blend4php (PHP): [https://github.com/galaxyproject/blend4php](https://github.com/galaxyproject/blend4php)
* clj-blend (Clojure): [https://github.com/chapmanb/clj-blend](https://github.com/chapmanb/clj-blend)

> <agenda-title></agenda-title>
>
> In this tutorial, we will learn the basics of BioBlend development:
>
> 1. TOC
> {:toc}
>
{:.agenda}

# Development on GitHub

Galaxy and BioBlend are developed as open source projects on GitHub and contributions are welcome!

The GitHub repositories can be found at:
* Galaxy: [https://github.com/galaxyproject/galaxy](https://github.com/galaxyproject/galaxy)
* BioBlend: [https://github.com/galaxyproject/bioblend](https://github.com/galaxyproject/bioblend)

##  Contributing on GitHub

To contribute to BioBlend, a GitHub account is required.

Changes are proposed via a {PR}. This allows the project maintainers to review the changes and suggest improvements.

The general steps are as follows:

1. Fork the BioBlend repository
2. Make changes in a new branch
3. Open a pull request for this branch in the upstream BioBlend repository

It is generally a good idea to enable the “Allow edits and access to secrets by maintainers” option for the {PR}. Enabling this option gives maintainers more freedom to help out.

# Downloading Galaxy and BioBlend

Now we are ready to set up our development environment! Since BioBlend communicates with the Galaxy {API}, we must also set up a local Galaxy server in order to test BioBlend functionality.
We use the [git](https://git-scm.com) versioning tool to download the repositories.

For this we require a public SSH key associated with our GitHub account. It makes sense to set this up now, since pushing changes to GitHub without a public key will prompt for credentials every time. See the [GitHub Docs](https://docs.github.com/en/github-ae@latest/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent) for information on setting up an SSH key.

> <hands-on-title>Downloading repositories</hands-on-title>
> Download Galaxy and set it up as a local repository:
>
> ```shell
> git clone git@github.com:galaxyproject/galaxy.git
> ```
>
> Likewise, download BioBlend:
>
> ```shell
> git clone git@github.com:galaxyproject/bioblend.git
> ```
>
> If you have already forked these repositories, you probably want to replace the links in the commands above to point to the forks.
{:.hands_on}

# Interacting with the API

## Structure of the Galaxy API

The source code for the {API} endpoints is contained in various files under [lib/galaxy/webapps/galaxy/api/](https://github.com/galaxyproject/galaxy/tree/dev/lib/galaxy/webapps/galaxy/api) in the Galaxy source code. Each of these files contains a controller which exposes functionality for a specific entity. For example, the dataset-related functionality is contained in the `DatasetsController` class in the [datasets.py](https://github.com/galaxyproject/galaxy/blob/v21.01/lib/galaxy/webapps/galaxy/api/datasets.py#L33) file.

Additionally, the [buildapp.py](https://github.com/galaxyproject/galaxy/blob/dev/lib/galaxy/webapps/galaxy/buildapp.py) file contains a complete listing of all the endpoints.

> <comment-title>Note</comment-title>
> Different versions of Galaxy differ in the functionality that is implemented. The development version is named `dev` and includes the most recent changes. Release versions are named according to their release date. These Galaxy releases can be accessed via git branches, e.g. `release_20.01`. At the time of writing BioBlend supports Galaxy releases 17.09 and later, although some functionality is not supported by all of these versions.
{:.comment}


> <tip-title>Going deeper into the Galaxy back-end code</tip-title>
> The various Galaxy {API} methods that are exposed to the outside are contained in controller classes. These contain the endpoint methods corresponding to the URLs of the Galaxy {API}. These endpoint methods handle the incoming requests from BioBlend.
>
> A focus for development of the Galaxy {API} is to make these controllers as “terse” as possible. This basically means that any logic not strictly required by the {API} endpoint method is moved to a corresponding manager class. This approach separates the {API} more cleanly from internal functionality. The manager classes should contain as much of the functionality as possible. An example is the [`DatasetManager`](https://github.com/galaxyproject/galaxy/blob/v21.01/lib/galaxy/managers/datasets.py#L26), referenced to regularly by the `DatasetsController` which is exposed directly to the {API}.
>
> At the time of writing this transition is ongoing. Therefore, source code of methods for one controller might look and function differently than those of another controller. This might be encountered when debugging a failing request in BioBlend.
>
{: .tip}



## Core concepts of the Galaxy interface

This section provides brief descriptions of the core concepts of the Galaxy interface, from the perspective of a BioBlend/Galaxy developer, rather than the typical Galaxy user. This list is by no means complete, but it should help those not already familiar with Galaxy.

BioBlend provides separate clients for interacting with each of these entities. For example, the BioBlend client for jobs is the `JobsClient` class. The source code for each `Client` can be found under [bioblend/galaxy/](https://github.com/galaxyproject/bioblend/tree/main/bioblend/galaxy).

<!-- This is the syntax for a definition list in markdown, it's semantically the most correct option here -->

Tool
: A Galaxy tool packages a particular piece of scientific software, keeping track of its required dependencies. It also contains additional metadata such as the BibTex references to the associated research paper(s). Tools generally take datasets and/or dataset collections as inputs and produce new datasets and/or collections as outputs.

Dataset
: A dataset represents digital data that can serve as input for a job. It can be associated with (contained in) a history or a library. A dataset can be shared between users. Each time a dataset is copied between histories or libraries, a new {HDA} or {LDA} is created. This naming scheme is used internally by the Galaxy relational database, but the terms might be encountered during BioBlend development as well.

Dataset Collection
: A dataset collection represents a group of related datasets. Like a dataset, it can be associated with a history or a library. Each time it is copied or shared, a new {HDCA} or a {LDCA} is created.<br><br>There are two kinds of dataset collections, `list` and `paired`. The former is a simple list of elements. The latter has a particular application in bioinformatics; it contains only two elements, corresponding to the `forward` and `reverse` strands of a piece of DNA.<br><br> Collection elements can be other collections; in other words, collections can be nested. These are designated (for example) as `list:list` or `list:paired`.

Job
: Jobs are associated with the execution of a tool on the Galaxy server. Each dataset is created by running a tool and therefore has an associated job. The computation of each step in a workflow has an associated job as well.<br><br>Jobs run as background processes on the Galaxy server. This is relevant because reading results before the appropriate jobs have finished can lead to missing data in the output. BioBlend methods must handle this by waiting for relevant jobs to finish before returning a response to the user.<br><br>Job creation can be different when dealing with dataset collections. If a tool input receives a collection instead of a dataset, it will create a separate job for each element, resulting in an output collection with the same structure as the input. This is referred to as 'mapping over' a collection.

History
: A history represents a container that keeps track of actions performed on its contents. Datasets are associated with a particular history; tools and workflows also need to be run in a particular history, which keeps track of the outputs produced by each of the resulting jobs or invocations. While the history is a very important concept to users of the graphical user interface, as it provides a constantly visible workspace to organize individual analyses or projects, it is slightly less relevant for BioBlend users.

Workflow
: A workflow represents a computation pipeline, made up of multiple tools chained together. It is created and contained in a history. A workflow is comprised of steps, which can be either inputs or tools. Its structure can be imagined as a directed graph of nodes that process inputs and pass on their outputs to the next nodes. The initial input to a workflow is commonly one or multiple datasets or dataset collections, as well as input parameters.

Invocation
: Just like a job represents a tool executed on certain inputs, an invocation represents the execution of a workflow on specified inputs and with certain specified parameters. Every invocation corresponds to only one workflow, and every step of the workflow has an associated invocation step. Each invocation step has jobs associated with it.

Library
: A library represents a container for datasets and dataset collections. Datasets that are contained in a library can be shared between multiple Galaxy users.

## Structure of the BioBlend library

BioBlend methods for the Galaxy API are associated with the various `Clients` already mentioned, stored under [bioblend/galaxy/](https://github.com/galaxyproject/bioblend/tree/main/bioblend/galaxy) in the BioBlend GitHub repository. The functionality for each Galaxy entity is in a separate folder. For example, the methods for interacting with workflows are stored under [bioblend/galaxy/workflows](https://github.com/galaxyproject/bioblend/tree/main/bioblend/galaxy/workflows) in the `WorkflowsClient` class.

The tests for these methods are stored under [tests/](https://github.com/galaxyproject/bioblend/tree/main/bioblend/_tests). For example, the tests related to workflows are in the file [tests/TestGalaxyWorkflows.py](https://github.com/galaxyproject/bioblend/blob/main/bioblend/_tests/TestGalaxyWorkflows.py).

The classic approach for accessing the Galaxy API is using the various clients and their methods. Each of these methods corresponds to one of the API endpoints. They send a request to the Galaxy server and return a Python object representing the parsed JSON response. Alternatively, there is also BioBlend.objects ({% cite Leo2014 %}), which provides an object-oriented API. Interaction with the API occurs via wrapper classes. These wrappers represent entity instances and provide methods to interact with and manipulate them. For example, a wrapper could represent a single workflow.


> <question-title>What is the difference between BioBlend and BioBlend.objects?</question-title>
>
> Look at the following code:
> ```python
> # classic BioBlend
> from bioblend.galaxy import GalaxyInstance
> gi = GalaxyInstance(url=galaxy_url, key=galaxy_key)
> workflows = gi.workflows.get_workflows()
> for wf_dict in workflows:
>     print(f"{wf_dict['name']} : {wf_dict['id']}"))
>
> # BioBlend.objects
> from bioblend.galaxy.objects import GalaxyInstance
> gi = GalaxyInstance(galaxy_url, galaxy_key)
> workflows = gi.workflows.get_previews()
> for wf in workflows:
>     print(f"{wf.name} : {wf.id}")
> ```
> It does the same thing (prints all available workflows) using both BioBlend and BioBlend.objects. What differences can you see?
> > <solution-title></solution-title>
> >
> > The output of the classic BioBlend code is provided as JSON, which the user can continue to manipulate themselves. The BioBlend.objects code provides the same information, but in the format of a Python object. This has various properties which the user can access, using for example `name` or `id`.
> >
> > You can view [numerous examples of both here](https://github.com/galaxyproject/bioblend/tree/main/docs/examples).
> {: .solution }
{: .question}


> <comment-title>Note</comment-title>
> BioBlend also provides methods for interacting with the Galaxy [ToolShed](https://galaxyproject.org/toolshed/) and [Cloudman](https://galaxyproject.org/cloudman/). However, this tutorial focuses on the part of BioBlend which provides access to the Galaxy API.
{:.comment}

# Communicating with the Galaxy API

We will make some manual requests to get a better feel for the Galaxy API. But first we need a valid API key in order to communicate with the Galaxy server through the API.

{% snippet faqs/galaxy/preferences_admin_api_key.md dev=true %}

In the following examples we make use of a simple networking tool named [curl](https://curl.se). Each command is listed together with the corresponding BioBlend and BioBlend.objects code.

We assume basic familiarity with HTTP request methods. For an overview, see the Mozilla Developer Network [documentation](https://developer.mozilla.org/en-US/docs/Web/HTTP/Methods).

Now, open a terminal window and let's make a few manual requests! In each case, the corresponding BioBlend and BioBlend.objects code is also provided for comparison.

> <hands-on-title>Making requests to the API</hands-on-title>
> 1. __GET request__
>
>    This is a simplest type of request. We request data from a specific URL, but we do not provide any further information apart from our API key. In this example we request a listing of all our invocations.
>
>    > <code-in-title>curl</code-in-title>
>    > ```shell
>    > curl -X GET -H 'x-api-key: <API_KEY>' 'localhost:8080/api/invocations'
>    > ```
>    >
>    > The `-X` flag specifies the HTTP method of our request. In this case we are making a GET request. GET is the default method, so we could also simply omit it here.
>    >
>    > The `-H` flag is used to specify HTTP headers. We authenticate our requests by specifying our API key in the `x-api-key` header.
>    {:.code-in}
>
>    > <code-out-title>Example JSON response</code-out-title>
>    > ```javascript
>    > [{
>    >     'history_id': '2f94e8ae9edff68a',
>    >     'id': 'df7a1f0c02a5b08e',
>    >     'model_class': 'WorkflowInvocation',
>    >     'state': 'new',
>    >     'update_time': '2015-10-31T22:00:22',
>    >     'uuid': 'c8aa2b1c-801a-11e5-a9e5-8ca98228593c',
>    >     'workflow_id': '03501d7626bd192f',
>    > }]
>    > ```
>    >
>    > The above response would mean that our user has only invoked a workflow once. The workflow's encoded ID is listed (`workflow_id`), as well as the encoded ID of the history containing this workflow (`history_id`), and the ID of the invocation itself (`id`). The state `new` indicates that the request to invoke this workflow is still pending and that the invocation has not been scheduled yet.
>    {:.code-out}
>
>    BioBlend code:
>
>    ```python
>    from bioblend.galaxy import GalaxyInstance
>    gi = GalaxyInstance('http://localhost:8080', key=<API_KEY>)
>    invs = gi.invocations.get_invocations()
>    print(invs)
>    ```
>
>    BioBlend.objects code:
>
>    ```python
>    from bioblend.galaxy.objects import GalaxyInstance
>    obj_gi = GalaxyInstance('http://localhost:8080', key=<API_KEY>)
>    previews = obj_gi.invocations.get_previews()
>    for preview in previews:
>        print(preview.id)
>    ```
>
> 2. __GET request with query parameters__
>
>    > <code-in-title>curl</code-in-title>
>    > In this example we include two query parameters in our request, `limit` and `include_terminal`.
>    >
>    > ```shell
>    > curl -X GET -H 'x-api-key: <API_KEY>' \
>    >      'localhost:8080/api/invocations?limit=5&include_terminal=False'
>    > ```
>    >
>    > With `limit=5` we limit the result to at most 5 invocations. With `include_terminal=False` we indicate that invocations in a terminal state should be excluded from the results.
>    {:.code-in}
>
>    > <code-out-title>Example JSON response</code-out-title>
>    > ```javascript
>    > []
>    > ```
>    >
>    > The above response would mean that there are no workflows in a non-terminal state - any workflows which have been invoked have been successfully scheduled.
>    {:.code-out}
>
>    BioBlend code:
>    ```python
>    invs = gi.invocations.get_invocations(limit=5, include_terminal=False)
>    print(invs)
>    ```
>
>    BioBlend.objects code:
>    ```python
>    invs = obj_gi.invocations.list(limit=5, include_terminal=False)
>    for invocation in invs:
>        print(invocation.id)
>    ```
>
> 3. __POST request with a payload__
>
>    > <code-in-title>curl</code-in-title>
>    > ```shell
>    > curl -X POST \
>    >     -H 'x-api-key: <API_KEY>' \
>    >     -H 'Content-Type: application/json' \
>    >     -d '{"name":"MyNewHistory"}' \
>    >     'localhost:8080/api/histories'
>    > ```
>    >
>    > This creates a new history with the name ‘MyNewHistory’.
>    {:.code-in}
>
>    > <code-out-title>Example JSON response</code-out-title>
>    > ```javascript
>    > {'update_time': '2021-05-24T13:04:05.393590',
>    > 'user_id': 'e9376175b7f708de',
>    > 'name': 'MyNewHistory',
>    > 'genome_build': None,
>    > 'username_and_slug': None,
>    > 'deleted': False,
>    > 'published': False,
>    > 'slug': None,
>    > 'model_class': 'History',
>    > 'state': 'new',
>    > 'id': '464411b81613694d',
>    > 'empty': True,
>    > 'annotation': None,
>    > 'importable': False,
>    > 'state_ids': {'new': [],
>    >  'upload': [],
>    >  'queued': [],
>    >  'running': [],
>    >  'ok': [],
>    >  'empty': [],
>    >  'error': [],
>    >  'discarded': [],
>    >  'paused': [],
>    >  'setting_metadata': [],
>    >  'failed_metadata': []},
>    > 'create_time': '2021-05-24T13:04:05.393579',
>    > 'purged': False,
>    > 'tags': [],
>    > 'contents_url': '/api/histories/464411b81613694d/contents',
>    > 'size': 0,
>    > 'state_details': {'new': 0,
>    >  'upload': 0,
>    >  'queued': 0,
>    >  'running': 0,
>    >  'ok': 0,
>    >  'empty': 0,
>    >  'error': 0,
>    >  'discarded': 0,
>    >  'paused': 0,
>    >  'setting_metadata': 0,
>    >  'failed_metadata': 0},
>    > 'url': '/api/histories/464411b81613694d'}
>    > ```
>    >
>    > The response contains some basic information about the newly created history, including the `history_id`, which we would need for further API calls.
>    {:.code-out}
>
>    BioBlend code:
>    ```python
>    history = gi.histories.create_history('MyNewHistory')
>    print(history)
>    ```
>
>    BioBlend.objects code:
>
>    ```python
>    history = obj_gi.histories.create('MyNewHistory')
>    print(history.name)
>    ```
{:.hands_on}

> <tip-title>Encoded and non-encoded IDs in Galaxy</tip-title>
> While making the API calls now, we encountered various IDs - for example, we created a history with ID `464411b81613694d`. How do Galaxy IDs work?
>
> The Galaxy server assigns a numeric (integer) ID to every entity, which serves as its unique identifier. (The exception are tools, which have a predefined string ID specified in the [XML wrapper](https://docs.galaxyproject.org/en/latest/dev/schema.html#attributes) defining the tool.) For example, each workflow, invocation, history and dataset is assigned a unique ID. If you access the Galaxy database directly, using a tool such as [gxadmin](https://github.com/usegalaxy-eu/gxadmin), these integers are the IDs that you will see.
>
> However, these IDs are not (or should not be!) exposed through the API. The Galaxy server encodes the numeric IDs into hexadecimal strings before sending them to external API clients. These are referred to as encoded IDs. BioBlend only deals with these encoded IDs.
>
{: .tip}

> <tip-title>Params versus payload: what's the difference?</tip-title>
>
> In BioBlend, the `Client._get()` helper method is used for GET requests. It takes a `params` argument. The `Client._post()` and `Client._put()` helper methods make POST and PUT requests respectively. They expect a `payload` argument. In this section we briefly explain the difference.
>
> An HTTP request is comprised of (1) a request line, (2) headers with metadata, (3) an empty line and (4) an optional message body. The HTTP GET method uses query parameters to specify which resource should be retrieved from the server. These parameters are included in the request line as part of the URL. The HTTP POST and PUT methods create or update resources. This data is specified in the payload, which is included in the message body of the request.
>
> We can see the difference using [curl](https://curl.se) and [ncat](https://nmap.org/ncat/).
>
> Open a terminal and run the `ncat` command in listen mode on port 8080. The loop will print out any incoming requests; each request is separated by a line containing `---`.
>
> ```shell
> while ncat --listen 8080; do printf "\n---\n"; done
> ```
>
> > <comment-title>Note</comment-title>
> > The `ncat` command is a low-level networking tool that does not respond to any incoming connections. Since HTTP is done over [TCP](https://en.wikipedia.org/wiki/Transmission_Control_Protocol), the `curl` command will hang expecting a response. Simply press <kbd>CTRL</kbd> + <kbd>C</kbd> to terminate the `curl` command.
> {:.comment}
>
> #### Parameters
>
> Open a second terminal and send a GET request with a query parameter:
>
> > <code-in-title>curl</code-in-title>
> > ```shell
> > curl -X GET localhost:8080/?message=Hello
> > ```
> {:.code-in}
>
> > <code-out-title>ncat</code-out-title>
> > ```
> > GET /?message=Hello HTTP/1.1
> > Host: localhost:8080
> > User-Agent: curl/7.76.1
> > Accept: */*
> >
> >
> > ```
> >
> > The two blank lines indicate that the payload is empty.
> {:.code-out}
>
> #### Payload
>
> Send a POST request with the payload ‘My name is curl’.
>
> > <code-in-title>curl</code-in-title>
> > ```shell
> > curl -X POST -d 'My name is curl' localhost:8080
> > ```
> {:.code-in}
>
> > <code-out-title>ncat</code-out-title>
> > ```
> > POST / HTTP/1.1
> > Host: localhost:8080
> > User-Agent: curl/7.76.1
> > Accept: */*
> > Content-Length: 15
> > Content-Type: application/x-www-form-urlencoded
> >
> > My name is curl
> > ```
> >
> > The message body now contains our payload.
> {:.code-out}
{: .tip}

# Creating a simple BioBlend method

As an introduction to BioBlend development, we will extend the BioBlend `ToolsClient` class by adding a new method named `uninstall_dependencies()`. This method will send a request to the Galaxy server to uninstall any dependencies for the specified tool.

> <comment-title>Note</comment-title>
> This example is based on the pull request “[Improving Tools API coverage](https://github.com/galaxyproject/bioblend/pull/390)” in the BioBlend repository.
{:.comment}

## Identifying the relevant Galaxy endpoint

We should check the Galaxy controller to verify whether this functionality is already implemented in the Galaxy API. Let's check the `dev` version of Galaxy. Our method interacts with the tools API, so let's check the [api/tools.py](https://github.com/galaxyproject/galaxy/blob/dev/lib/galaxy/webapps/galaxy/api/tools.py) file.

As it turns out, a corresponding [`uninstall_dependencies()`](https://github.com/galaxyproject/galaxy/blob/v21.01/lib/galaxy/webapps/galaxy/api/tools.py#L293) endpoint method is already available in the Galaxy `dev` branch. Lucky us! If this was not the case, then we would have first had to implement the endpoint in Galaxy, before writing a BioBlend method to interact with it.

> <hands-on-title>Creating the method and a simple test</hands-on-title>
> 1. __Constructing the correct method signature__
>
>    Let's look at the signature of the endpoint method:
>
>    __File [api/tools.py](https://github.com/galaxyproject/galaxy/blob/v21.01/lib/galaxy/webapps/galaxy/api/tools.py#L293), line 291__:
>    ```python
>    def uninstall_dependencies(self, trans: GalaxyWebTransaction, id, **kwds):
>    ```
>    The `GalaxyWebTransaction` is an object which represents our API request. It is internal to the Galaxy server. We do not specify it as a parameter in BioBlend.
>
>    We can see that it also expects an `id` parameter, which is the tool ID. The `kwds` argument indicates that the endpoint might accept additional keyword arguments as well.
>
>    Let's check the method documentation:
>    ```python
>    """
>    DELETE /api/tools/{tool_id}/dependencies
>    Attempts to uninstall requirements via the dependency resolver
>
>    parameters:
>
>        index:
>            index of dependency resolver to use when installing dependency.
>            Defaults to using the highest ranking resolver
>
>        resolver_type:
>            Use the dependency resolver of this resolver_type to install
>            dependency
>    """
>    ```
>    The endpoint expects a DELETE request at URL `/api/tools/{tool_id}/dependencies`. The `id` parameter in the method signature corresponds to the variable `{tool_id}` present in the endpoint URL.
>
>    There are two additional parameters: `index` and `resolver_type`. These parameters specify how the server should resolve the tool dependencies. For an average API user these options are probably too specific as they require knowledge of the resolvers used by the Galaxy server back end. We will skip these parameters for this example.
>
>    Lastly, we need to find out what the endpoint returns. Unfortunately the return type is not listed in the method documentation. It is a reality that various parts of the Galaxy API lack complete documentation, although this is being worked on.
>
>    Let's look at the return statement of the endpoint:
>    ```python
>    return tool.tool_requirements_status
>    ```
>    Here the tool API controller calls the `tool_requirements_status` property of the `Tool` class instance, which represents the tool corresponding to the `id` parameter. Let's also check the definition of this property.
>
>    __File [tools/\_\_init\_\_.py](https://github.com/galaxyproject/galaxy/blob/v21.01/lib/galaxy/tools/__init__.py#L1833), line 1833__:
>    ```python
>    def tool_requirements_status(self):
>        """
>        Return a list of dictionaries for all tool dependencies with their
>        associated status
>        """
>        return self._view.get_requirements_status(
>            {self.id: self.tool_requirements},
>            self.installed_tool_dependencies
>        )
>    ```
>    The documentation of this method mentions that it returns a list of dictionaries for all tool dependencies with their associated status. Now we know the return type of the endpoint method!
>
>    Our BioBlend method will receive this list of dictionaries encoded as JSON in the server response. It will parse the JSON and return the list. Thus, the signature of our new `uninstall_dependencies()` method will be as follows:
>
>    ```python
>    def uninstall_dependencies(self, tool_id: str) -> List[dict]:
>    ```
>
>    The `self` parameter corresponds to the `ToolsClient` instance. The parameter `tool_id` corresponds to the encoded ID for a given tool. This ID will be substituted in place of `{tool_id}` in the request URL.
>
>
> 2. __Adding the method body__
>
>    Our method is part of the `ToolsClient` class. This means that our method has access to all its helper methods.
>    Let's recall that the endpoint expects a DELETE request at URL `/api/tools/{tool_id}/dependencies`.
>
>    To translate this into BioBlend source code, we do the following:
>
>    - We use the `self._make_url()` helper method to construct a valid URL with the supplied `tool_id` parameter. Finally, we need to append `'/dependencies'` for our endpoint.
>
>        ```python
>        url = self._make_url(tool_id) + '/dependencies'
>        ```
>
>    - BioBlend clients also include helper methods for making requests of various types. In this case we want to make a DELETE request, so we use the `self._delete()` helper method with our constructed URL as argument. These helper methods are inherited from the `Client` class located at [bioblend/galaxy/client.py](https://github.com/galaxyproject/bioblend/blob/v0.16.0/bioblend/galaxy/client.py#L19).
>
>
>    The `_delete()` method also expects a `payload` argument. Our payload will be empty (`{}`), since we use `tool_id` only to construct the URL and there are no additional parameters to pass. This helper method automatically parses the JSON response into a Python object, so we can directly return its result.
>
>    ```python
>    return self._delete(payload={}, url=url)
>    ```
>
> 3. __Adding the docstring__
>
>    We should document the method and its parameter for the user.
>    ```python
>    """
>    Uninstall dependencies for a given tool via a resolver.
>
>    :type tool_id: str
>    :param tool_id: id of the requested tool
>
>    :rtype: list of dicts
>    :return: Tool requirement statuses
>    """
>    ```
>
>    This is the docstring format used in BioBlend. At the top we provide a short general description of the method. Then we document each parameter with its type and a short description. The return value of the method is documented last. After the PR is merged, the docstring will be rendered as part of the [BioBlend documentation](https://bioblend.readthedocs.io/en/latest/api_docs/galaxy/all.html).
>
> 4. __Checking Galaxy version compatibility__
>
>    In case certain versions of Galaxy do not support our method's functionality, we should append a note to the method's docstring.
>
>    __Example__:
>
>    ```
>    .. note::
>        This method is only supported by Galaxy 19.09 or later.
>    ```
>
>    Previously we made sure that the `dev` version of Galaxy supported our method. Let's now check if older versions of Galaxy also support it. At the time of writing, Galaxy 17.09 is the earliest version supported by BioBlend, so it makes sense to check it first.
>
>    It turns out that Galaxy 17.09 also [supports](https://github.com/galaxyproject/galaxy/blob/v17.09/lib/galaxy/webapps/galaxy/api/tools.py#L159) our method! It's safe to assume that all later versions do as well, and thus we do not need to add a note.
>
>     > <tip-title>Determining the Galaxy version from within BioBlend</tip-title>
>     > In certain situations it might be useful to determine the version of the Galaxy server from BioBlend. We can determine the version using the `self.gi.config.get_version()` method.
>     >
>     > __Example__:
>     >
>     > ```python
>     > if self.gi.config.get_version()['version_major'] >= '21.01':
>     >     print('Newer version of Galaxy!')
>     > else:
>     >     print('Older version of Galaxy!')
>     > ```
>     > An example where this is useful is the [`download_dataset_collection()`](https://github.com/galaxyproject/bioblend/blob/v0.16.0/bioblend/galaxy/dataset_collections/__init__.py#L125) method. It downloads a dataset collection as an archive file. Early versions of Galaxy return a ‘tar.gz’ archive. Later versions return a ‘zip’ archive. The BioBlend method uses the Galaxy version to determine the archive type.
>     {: .tip}
>
> 5. __Putting it all together__
>
>    In the end, our new method should look something like this:
>
>    ```python
>    def uninstall_dependencies(self, tool_id: str) -> List[dict]:
>        """
>        Uninstall dependencies for a given tool via a resolver.
>
>        :type tool_id: str
>        :param tool_id: id of the requested tool
>
>        :rtype: list of dicts
>        :return: Tool requirement statuses
>        """
>        url = self._make_url(tool_id) + '/dependencies'
>        return self._delete(payload={}, url=url)
>    ```
>
> 6. __Writing a test for our new method__
>
>    It is best to keep BioBlend tests simple. We need to check that we get an expected response, but testing minute details is unnecessary, since those are the responsibility of the Galaxy server itself.<br>
>
>    Therefore, let's check that, when run on a given tool ID, the returned status indicates that the dependencies are null. We can reuse a tool ID from an existent test method. Let's take `CONVERTER_fasta_to_bowtie_color_index`, which has only one dependency.
>
>    ```python
>    def test_tool_dependency_uninstall(self):
>        statuses = self.gi.tools.uninstall_dependencies(
>             'CONVERTER_fasta_to_bowtie_color_index'
>        )
>        self.assertEqual(statuses[0]['model_class'], 'NullDependency')
>    ```
>
>     > <tip-title>Skip tests with the `@test_util.skip_unless_galaxy` decorator</tip-title>
>     > On our local machine we normally test our changes against a single version of Galaxy, for example `dev`. However, GitHub Actions continuous integration runs the BioBlend tests against all supported versions of Galaxy.
>     >
>     > At the time of writing BioBlend supports Galaxy 17.09 and later. However, certain functionality might not be available in all these versions. For example, `download_dataset_collection` is only supported by Galaxy 18.01 and later. We don't want to run tests specific to this functionality on earlier versions of Galaxy. They would simply fail. We can skip a test for specific versions of Galaxy by adding the `@test_util.skip_unless_galaxy` decorator above the test method in question.
>     >
>     > __Example__:
>     > ```python
>     > @test_util.skip_unless_galaxy('release_19.09')
>     > def test_some_functionality(self):
>     > ```
>     > This test will only be run against Galaxy 19.09 and later.
>     {: .tip}
{:.hands_on}


> <tip-title>Galaxy API parameter formats</tip-title>
>
> Implementing a BioBlend method often requires knowledge about the parameters which are accepted by the corresponding Galaxy API endpoint. At the time of writing there are a few different formats in which endpoints might accept their input parameters.
>
> #### Direct parameters
> A parameter can be specified explicitly in the signature of the endpoint method.
>
> __Example__:
>
> In the histories API [`delete()`](https://github.com/galaxyproject/galaxy/blob/v21.01/lib/galaxy/webapps/galaxy/api/histories.py#L384) endpoint method, the `id` parameter is specified explicitly.
>
> ```python
> def delete(self, trans, id, **kwd):
> ```
>
> #### Keyword parameters
>
> Other parameters might get passed to the endpoint as an element of `kwd`. If the documentation of the endpoint is lacking, these keyword parameters can be difficult to identify.
>
> #### The `q` and `qv` format
>
> Filtering parameters might also be expected as pairs of `q` and `qv` values. These are then parsed into filters by the Galaxy server.
>
> The general format is `q={filter}-{operation}` and `qv={value}`. The accepted values for these parameters must be identified in order to construct valid `q` and `qv` pairs.
>
> __Example__:
>
> The Datasets API [`index()` method](https://github.com/galaxyproject/galaxy/blob/v21.01/lib/galaxy/webapps/galaxy/api/datasets.py#L54) uses this format. The `q` and `qv` parameters are passed to the endpoint as part of `kwd`. They are not used in the endpoint itself, but rather implicitly passed into the `self.parse_filter_params()` helper method as part of `kwd`. See [the corresponding method](https://github.com/galaxyproject/bioblend/blob/v0.16.0/bioblend/galaxy/datasets/__init__.py#L146) of the `DatasetClient` to get a feeling how this can be handled in BioBlend.
>
{: .tip}

# Running the BioBlend tests

Finally, we should run our test to make sure that it passes. This would indicate that our method works correctly. After testing we will be ready to push our changes to GitHub and open a pull request!

This section outlines two approaches for running BioBlend tests: A simple approach that should suffice for basic testing, and a more involved approach that speeds up testing.

It makes sense to test our method against the most recent changes of Galaxy, so let's use Galaxy `dev`.

> <comment-title>Note</comment-title>
> On our local machine we generally only test against Galaxy `dev`, so it might happen that tests on GitHub still fail because there tests are run against all supported Galaxy versions. When this happens we simply need to update our changes to fix those errors as well, and we might need to check out older Galaxy versions for debugging locally.
{:.comment}

> <hands-on-title>Using the `run_bioblend_tests.sh` script</hands-on-title>
>
> This script is provided in the repository. It will lint the Python code and run the BioBlend tests. It works well for running all BioBlend tests at once, or for small tests that require little to no debugging or experimentation.
>
> Since BioBlend requires a Galaxy server, we must specify the path to the Galaxy server directory:
>
> -  `-g <galaxy_path>`
>
> Two useful optional parameters:
>
> -  `-e <python_version>`
> -  `-t <path_to_test>`
>
>     We can specify a subset of tests to run by supplying this argument. This is specified in [pytest](https://pytest.org/) format. See the [documentation](https://docs.pytest.org/en/latest/how-to/usage.html#specifying-which-tests-to-run) for more information.
>
> __Examples__:
>
> Let's assume that our galaxy directory is located at `../galaxy`.
>
> Then this command runs all BioBlend tests:
> ```shell
> ./run_bioblend_tests -g ../galaxy -e py39
> ```
> And this command runs only the test named `test_get_jobs`:
> ```shell
> ./run_bioblend_tests.sh \
>     -g ../galaxy \
>     -e py39 \
>     -t tests/TestGalaxyJobs.py::TestGalaxyJobs::test_get_jobs
> ```
> __Downside of the `run_bioblend_tests.sh` script__:
>
> The script starts and stops a new instance of the Galaxy server every time it is run. This makes it robust. However, it also means we have to wait for the server to start up every time. This makes it painful to debug more complicated tests.
{:.hands_on}

## Running the tests directly

We can speed up our test environment by controlling the Galaxy server and the BioBlend tests directly. The initial setup will be more involved, but this way we can keep the server running in the background while we do our testing.

> <warning-title>Caveat</warning-title>
>
> When the Galaxy server is not reset between test runs, tests that implicitly depend on a clean server state might start failing after the first run.
>
> For example, let's assume there is only a single test that creates a new history and checks that the total number of histories is equal to 1. This test will pass for a new Galaxy server because no histories exist yet. However, the next run there would be two histories and thus the test would fail.
>
> Although this is generally not a big issue, it is something to watch out for. There are a few BioBlend tests that will fail in this manner. Most of these tests measure an array length property that increases with every run and therefore fail after the first run.
>
> It is preferable to write tests that are robust in this regard!
{:.warning}

> <hands-on-title>Starting the Galaxy server with a custom configuration</hands-on-title>
>
> Open a terminal in the Galaxy base directory and execute the following commands:
>
> 1. __Declare an API key of our choosing to use for testing__
>
>    We will include this key in the Galaxy server configuration. The server will then accept incoming requests that use this key.
>
>    ```shell
>    API_KEY=LetMeIn
>    ```
>
> 2. __Create a temporary directory__
>
>    ```shell
>    TEMP_DIR=$(mktemp -d)
>    echo "Galaxy directory: $TEMP_DIR"
>    ```
>
>    The temporary directory will contain useful information for debugging tests, such as the ‘main.log' file and the ‘universe.sqlite' database.
>
> 3. __Export environment variables required by the Galaxy server__
>
>    ```shell
>    export SKIP_GALAXY_CLIENT_BUILD=1
>    export GALAXY_LOG_FILE=$TEMP_DIR/main.log
>    export GALAXY_CONFIG_FILE=$TEMP_DIR/test_galaxy.ini
>    ```
>
> 4. __Copy a sample tool configuration__
>
>    ```shell
>    cp config/tool_conf.xml.sample $TEMP_DIR/tool_conf.xml.sample
>    ```
>
> 5. __Declare the desired logging level__
>
>    ```shell
>    LOG_LEVEL=DEBUG
>    ```
>
>    The available logging levels can be found in the [documentation](https://docs.python.org/3/library/logging.html#logging-levels) of the Python [logging](https://docs.python.org/3/library/logging.html) library.
>
> 6. __Create the Galaxy configuration file__
>
>    ```shell
>    cat << EOF > $GALAXY_CONFIG_FILE
>    [server:main]
>
>    use = egg:Paste#http
>    port = 8080
>
>    [app:main]
>
>    log_level = $LOG_LEVEL
>    paste.app_factory = galaxy.web.buildapp:app_factory
>    database_connection = sqlite:///$TEMP_DIR/universe.sqlite?isolation_level=IMMEDIATE
>    file_path = $TEMP_DIR/files
>    new_file_path = $TEMP_DIR/tmp
>    tool_config_file = $TEMP_DIR/tool_conf.xml.sample
>    conda_auto_init = True
>    job_working_directory = $TEMP_DIR/jobs_directory
>    allow_library_path_paste = True
>    admin_users = test@test.test
>    allow_user_deletion = True
>    allow_user_dataset_purge = True
>    enable_beta_workflow_modules = True
>    master_api_key = $API_KEY
>    enable_quotas = True
>    cleanup_job = onsuccess
>    EOF
>    ```
>
> 7. __Start the Galaxy server__
>
>    ```shell
>    ./run.sh
>    ```
{:.hands_on}

Now we have a Galaxy server running with our specified configuration. We can run tests against this server in a separate terminal.

> <hands-on-title>Running the tests</hands-on-title>
>
> Navigate to the BioBlend base directory and execute the following commands:
>
> 1. __Declare the same API key we used for the Galaxy server__
>
>    ```shell
>    API_KEY=LetMeIn
>    ```
>
> 2. __Export environment variables required by BioBlend__
>
>    ```shell
>    export BIOBLEND_GALAXY_API_KEY=$API_KEY
>    export BIOBLEND_GALAXY_URL=http://localhost:8080/
>    ```
>
> 3. __Activate the virtual environment__
>
>    ```shell
>    source <GALAXY_DIRECTORY>/.venv/bin/activate
>    ```
>
> 4. __Declare the desired test logging level__
>
>    ```shell
>    LOG_LEVEL=WARNING
>    ```
>
> 5. __Run the desired BioBlend test(s)__
>
>    ```shell
>    # selection format: tests/<file>::<module>::<test>
>    pytest --override-ini log_cli_level=$LOG_LEVEL \
>        tests/TestGalaxyWorkflows.py::TestGalaxyWorkflows::test_show_versions
>    ```
>
> 6. __Lint our code to detect bad practices and potential errors__
>
>    ```shell
>    tox -e lint
>    ```
>
> 7. __Leave the virtual environment when we are done testing__
>
>    ```shell
>    deactivate
>    ```
{:.hands_on}

You may find your tests do not pass straight away. In this case, you need to identify the issue and modify the code to fix it. This might take several iterations.

Once tests are passing, the final step would be to [open a PR](https://github.com/galaxyproject/bioblend/compare) at the BioBlend GitHub repository! This will trigger automatic tests and after some time your new code will be reviewed by the BioBlend maintainers. After any concerns are addressed and all tests are passing, your PR will be merged into the `main` branch and incorporated into the next BioBlend release.

# Conclusion

We covered the basics of BioBlend development. We also touched briefly on some key concepts relating to the Galaxy API and the structure of the Galaxy back-end code.

Compared to Galaxy, BioBlend is a smaller project with limited complexity. Changes are generally not very hard to implement, which makes it a good candidate for contributions from those not yet familiar with Galaxy in its entirety.

Good luck!
