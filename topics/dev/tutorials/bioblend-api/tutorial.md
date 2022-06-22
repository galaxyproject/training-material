---
layout: tutorial_hands_on

title: "Scripting Galaxy using the API and BioBlend"
level: Introductory
requirements: []
follow_up_training: []
questions:
  - What is a REST API?
  - How to interact with Galaxy programmatically?
  - Why and when should I use BioBlend?

objectives:
  - Interact with Galaxy via BioBlend.

time_estimation:  2h
key_points:
  - The API allows you to use Galaxy's capabilities programmatically.
  - BioBlend makes using the Galaxy API from Python easier.
  - BioBlend objects is an object-oriented interface for interacting with Galaxy.

subtopic: api
contributions:
  authorship:
  - nsoranzo
  editing:
  - claresloggett
  - nturaga

notebook:
  language: python
  pyolite: true
---



> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Interacting with histories in Galaxy API

We are going to use the [requests](http://python-requests.org/) Python library to communicate via HTTP with the Galaxy server. To start, let's define the connection parameters.

**You need to insert the API key for your Galaxy server in the cell below**:
1. Open the Galaxy server in another browser tab
2. Click on "User" on the top menu, then "Preferences"
3. Click on "Manage API key"
4. Generate an API key if needed, then copy the alphanumeric string and paste it as the value of the `api_key` variable below.

```python
from pprint import pprint

import json

import requests
from six.moves.urllib.parse import urljoin

server = 'https://usegalaxy.eu/'
api_key = ''
base_url = urljoin(server, 'api')
base_url
```

We now make a GET request to retrieve all histories owned by a user:

```python
params = {'key': api_key}
r = requests.get(base_url + '/histories', params)
print(r.text)
hists = r.json()
pprint(hists)
```

As you can see, GET requests in Galaxy API return JSON strings, which need to be **deserialized** into Python data structures. In particular, GETting a resource collection returns a list of dictionaries.

Each dictionary returned when GETting a resource collection gives basic info about a resource, e.g. for a history you have:
- `id`: the unique **identifier** of the history, needed for all specific requests about this resource
- `name`: the name of this history as given by the user
- `deleted`: whether the history has been deleted
- `url`: the relative URL to get all info about this resource.

There is no readily-available filtering capability, but it's not difficult to filter histories **by name**:

```python
pprint([_ for _ in hists if _['name'] == 'Unnamed history'])
```

If you are interested in more **details** about a given resource, you just need to append its `id` to the previous collection request, e.g. to the get more info for a history:

```python
hist0_id = hists[0]['id']
print(hist0_id)
params = {'key': api_key}
r = requests.get(base_url + '/histories/' + hist0_id, params)
pprint(r.json())
```

As you can see, there are much more entries in the returned dictionary, e.g.:
- `create_time`
- `size`: total disk space used by the history
- `state_ids`: ids of history datasets for each possible state.

To get the list of **datasets contained** in a history, simply append `/contents` to the previous resource request.

```python
params = {'key': api_key}
r = requests.get(base_url + '/histories/' + hist0_id + '/contents', params)
hdas = r.json()
pprint(hdas)
```

The dictionaries returned when GETting the history content give basic info about each dataset, e.g.: `id`, `name`, `deleted`, `state`, `url`...

To get the details about a specific dataset, you can use the `datasets` controller:

```python
hda0_id = hdas[0]['id']
print(hda0_id)
params = {'key': api_key}
r = requests.get(base_url + '/datasets/' + hda0_id, params)
pprint(r.json())
```

Some of the interesting additional dictionary entries are:
- `create_time`
- `creating job`: id of the job which created this dataset
- `download_url`: URL to download the dataset
- `file_ext`: the Galaxy data type of this dataset
- `file_size`
- `genome_build`: the genome build (dbkey) associated to this dataset.

**New resources** are created with POST requests. The uploaded **data needs to be serialized** in a JSON string. For example, to create a new history:

```python
params = {'key': api_key}
data = {'name': 'New history'}
r = requests.post(base_url + '/histories', data=json.dumps(data), params=params, headers={'Content-Type': 'application/json'})
new_hist = r.json()
pprint(new_hist)
```

The return value of a POST request is a dictionary with detailed info about the created resource.

To **update** a resource, make a PUT request, e.g. to change the history name:

```python
params = {'key': api_key}
data = {'name': 'Updated history'}
r = requests.put(base_url + '/histories/' + new_hist['id'], json.dumps(data), params=params, headers={'Content-Type': 'application/json'})
print(r.status_code)
pprint(r.json())
```

The return value of a PUT request is usually a dictionary with detailed info about the updated resource.

Finally to **delete** a resource, make a DELETE request, e.g.:

```python
params = {'key': api_key}
r = requests.delete(base_url + '/histories/' + new_hist['id'], params=params)
print(r.status_code)
```


## Exercise: Galaxy API

**Goal**: Upload a file to a new history, import a workflow and run it on the uploaded dataset.

> ### {% icon question %} Question: Initialise
> First, define the connection parameters. What variables do you need?
> > ### {% icon solution %} Solution
> > ```python
> > from pprint import pprint
> >
> > import json
> >
> > import requests
> > from six.moves.urllib.parse import urljoin
> >
> > server = 'https://usegalaxy.eu/'
> > api_key = ''
> > base_url = urljoin(server, 'api')
> > ```
> >
> {: .solution}
{: .question}

```python
# Try it out here!

```


> ### {% icon question %} Question: New History
> Next, create a new Galaxy history via POST to the correct API.
> > ### {% icon solution %} Solution
> > ```python
> > params = {'key': api_key}
> > data = {'name': 'New history'}
> > r = requests.post(base_url + '/histories', json.dumps(data), params=params, headers={'Content-Type': 'application/json'})
> > new_hist = r.json()
> > pprint(new_hist)
> > ```
> >
> {: .solution}
{: .question}

```python
# Try it out here!

```


> ### {% icon question %} Question: Upload a dataset
> **Upload** the local file "test-data/1.txt" to a new dataset in the created history. You need to run the special 'upload1' tool by making a `POST` request to `/api/tools`. You don't need to pass any inputs to it apart from attaching the file as 'files_0|file_data'.
>
> You can obtain this file from the following URL, you'll need to download it first.
>
> ```
> https://raw.githubusercontent.com/nsoranzo/bioblend-tutorial/main/test-data/1.txt
> ```
>
> > ### {% icon solution %} Solution
> > ```python
> > params = {'key': api_key}
> > data = {
> >     'history_id': new_hist['id'],
> >     'tool_id': 'upload1'}
> > with open("test-data/1.txt", 'rb') as f:
> >     files = {'files_0|file_data': f}
> >     r = requests.post(base_url + '/tools', data, params=params, files=files)
> > ret = r.json()
> > pprint(ret)
> > ```
> >
> {: .solution}
{: .question}

```python
# Try it out here!

```


> ### {% icon question %} Question: Find the dataset in your history
> Find the new uploaded dataset, either from the dict returned by the POST request or from the history contents.
> > ### {% icon solution %} Solution
> > ```python
> > hda = ret['outputs'][0]
> > pprint(hda)
> > ```
> >
> {: .solution}
{: .question}

```python
# Try it out here!

```

> ### {% icon question %} Question: Import a workflow
> **Import a workflow** from the local file "test-data/convert_to_tab.ga" by making a `POST` request to `/api/workflows`. The only needed data is 'workflow', which must be a deserialized JSON representation of the workflow.
>
> You can obtain this file from the following URL, you'll need to download it first.
>
> ```
> https://raw.githubusercontent.com/nsoranzo/bioblend-tutorial/main/test-data/convert_to_tab.ga
> ```
> > ### {% icon solution %} Solution
> >
> > ```python
> > params = {'key': api_key}
> > with open('test-data/convert_to_tab.ga', 'r') as f:
> >     workflow_json = json.load(f)
> > data = {'workflow': workflow_json}
> > r = requests.post(base_url + '/workflows', json.dumps(data), params=params, headers={'Content-Type': 'application/json'})
> > wf = r.json()
> > pprint(wf)
> > ```
> {: .solution}
{: .question}

```python
# Try it out here!

```

> ### {% icon question %} Question: View the WF details
> View the details of the imported workflow by making a GET request to `/api/workflows`.
> > ### {% icon solution %} Solution
> > ```python
> > params = {'key': api_key}
> > r = requests.get(base_url + '/workflows/' + wf['id'], params)
> > wf = r.json()
> > pprint(wf)
> > ```
> >
> {: .solution}
{: .question}

```python
# Try it out here!

```


> ### {% icon question %} Question: Run the WF
> **Run** the imported workflow on the uploaded dataset **inside the same history** by making a `POST` request to `/api/workflows/<id>/invocations`. The only needed data are 'history' and 'inputs'.
> > ### {% icon solution %} Solution
> >
> > ```python
> > params = {'key': api_key}
> > inputs = {0: {'id': hda['id'], 'src': 'hda'}}
> > data = {
> >     'history': 'hist_id=' + new_hist['id'],
> >     'inputs': inputs}
> > r = requests.post(base_url + '/workflows/' + wf['id'] + '/invocations', json.dumps(data), params=params, headers={'Content-Type': 'application/json'})
> > pprint(r.json())
> > ```
> {: .solution}
{: .question}


> ### {% icon question %} Question: View the results
> View the results on the Galaxy server with your web browser. Were you successful? Did it run?
{: .question}

# Interacting with histories in BioBlend

**You need to insert the API key for your Galaxy server in the cell below**:
1. Open the Galaxy server in another browser tab
2. Click on "User" on the top menu, then "Preferences"
3. Click on "Manage API key"
4. Generate an API key if needed, then copy the alphanumeric string and paste it as the value of the `api_key` variable below.

The user interacts with a Galaxy server through a `GalaxyInstance` object:

```python
from pprint import pprint

import bioblend.galaxy

server = 'https://usegalaxy.eu/'
api_key = ''
gi = bioblend.galaxy.GalaxyInstance(url=server, key=api_key)
```

The `GalaxyInstance` object gives you access to the various controllers, i.e. the resources you are dealing with, like `histories`, `tools` and `workflows`.
Therefore, method calls will have the format `gi.controller.method()`. For example, the call to retrieve all histories owned by the current user is:

```python
pprint(gi.histories.get_histories())
```

As you can see, methods in BioBlend do not return JSON strings, but **deserialize** them into Python data structures. In particular, `get_` methods return a list of dictionaries.

Each dictionary gives basic info about a resource, e.g. for a history you have:
- `id`: the unique **identifier** of the history, needed for all specific requests about this resource
- `name`: the name of this history as given by the user
- `deleted`: whether the history has been deleted
- `url`: the relative URL to get all info about this resource.

**New resources** are created with `create_` methods, e.g. the call to create a new history is:

```python
new_hist = gi.histories.create_history(name='BioBlend test')
pprint(new_hist)
```

As you can see, to make POST requests in BioBlend it is **not necessary to serialize data**, you just pass them explicitly as parameters. The return value is a dictionary with detailed info about the created resource.

`get_` methods usually have **filtering** capabilities, e.g. it is possible to filter histories **by name**:

```python
pprint(gi.histories.get_histories(name='BioBlend test'))
```

It is also possible to specify the unique **id** of the resource to retrieve, e.g. to get back the history we created before:

```python
hist_id = new_hist['id']
pprint(gi.histories.get_histories(history_id=hist_id))
```

Please note that independently of which parameters are passed to the `get_` method, it always returns a list.


To **upload** files to the new history, run the special upload tool by calling the `upload_file` method of the `tools` controller:

```python
pprint(gi.tools.upload_file('test-data/1.txt', hist_id))
```

If you are interested in more **details** about a given resource, you can use the corresponding `show_` method. For example, to the get more info for the history we have just populated:

```python
pprint(gi.histories.show_history(history_id=hist_id))
```

As you can see, there are much more entries in the returned dictionary, e.g.:
- `create_time`
- `size`: total disk space used by the history
- `state_ids`: ids of history datasets for each possible state.

To get the list of **datasets contained** in a history, simply add `contents=True` to the previous call.

```python
hdas = gi.histories.show_history(history_id=hist_id, contents=True)
pprint(hdas)
```

The dictionaries returned when showing the history content give basic info about each dataset, e.g.: `id`, `name`, `deleted`, `state`, `url`...

To get the details about a specific dataset, you can use the `datasets` controller:

```python
hda0_id = hdas[0]['id']
print(hda0_id)
pprint(gi.datasets.show_dataset(hda0_id))
```

Some of the interesting additional dictionary entries are:
- `create_time`
- `creating job`: id of the job which created this dataset
- `download_url`: URL to download the dataset
- `file_ext`: the Galaxy data type of this dataset
- `file_size`
- `genome_build`: the genome build (dbkey) associated to this dataset.

To **update** a resource, use the `update_` method, e.g. to change the name of the new history:

```python
pprint(gi.histories.update_history(new_hist['id'], name='Updated history'))
```

The return value of `update_` methods is usually a dictionary with detailed info about the updated resource.

Finally to **delete** a resource, use the `delete_` method, e.g.:

```python
pprint(gi.histories.delete_history(new_hist['id']))
```

## Exercise: BioBlend

**Goal**: Upload a file to a new history, import a workflow and run it on the uploaded dataset.


> ### {% icon question %} Question: Initialise
> Create a `GalaxyInstance` object.
> > ### {% icon solution %} Solution
> > ```python
> > from pprint import pprint
> >
> > import bioblend.galaxy
> >
> > server = 'https://usegalaxy.eu/'
> > api_key = ''
> > gi = bioblend.galaxy.GalaxyInstance(url=server, key=api_key)
> > ```
> >
> {: .solution}
{: .question}

```python
# Try it out here!

```


> ### {% icon question %} Question: New History
> Create a new Galaxy history.
> > ### {% icon solution %} Solution
> > ```python
> > new_hist = gi.histories.create_history(name='New history')
> > pprint(new_hist)
> > ```
> {: .solution}
{: .question}

```python
# Try it out here!

```


> ### {% icon question %} Question: Upload a dataset
> **Upload** the local file "test-data/1.txt" to a new dataset in the created history using `tools.upload_file()` .
>
> You can obtain this file from the following URL, you'll need to download it first.
>
> ```
> https://raw.githubusercontent.com/nsoranzo/bioblend-tutorial/main/test-data/1.txt
> ```
> > ### {% icon solution %} Solution
> > ```python
> > ret = gi.tools.upload_file("test-data/1.txt", new_hist['id'])
> > pprint(ret)
> > ```
> {: .solution}
{: .question}

```python
# Try it out here!

```


> ### {% icon question %} Question: Find the dataset in your history
> Find the new uploaded dataset, either from the dict returned by `tools.upload_file()` or from the history contents.
> > ### {% icon solution %} Solution
> > ```python
> > hda = ret['outputs'][0]
> > pprint(hda)
> > ```
> {: .solution}
{: .question}



> ### {% icon question %} Question: Import a workflow
> **Import a workflow** from the local file "test-data/convert_to_tab.ga" using `workflows.import_workflow_from_local_path()` .
>
> You can obtain this file from the following URL, you'll need to download it first.
>
> ```
> https://raw.githubusercontent.com/nsoranzo/bioblend-tutorial/main/test-data/convert_to_tab.ga
> ```
> > ### {% icon solution %} Solution
> > ```python
> > wf = gi.workflows.import_workflow_from_local_path('test-data/convert_to_tab.ga')
> > pprint(wf)
> > ```
> {: .solution}
{: .question}

```python
# Try it out here!

```

> ### {% icon question %} Question: View the WF details
> View the details of the imported workflow using `workflows.show_workflow()`
> > ### {% icon solution %} Solution
> > ```python
> > wf = gi.workflows.show_workflow(wf['id'])
> > pprint(wf)
> > ```
> {: .solution}
{: .question}

```python
# Try it out here!

```

> ### {% icon question %} Question: Run the WF
> **Run** the imported workflow on the uploaded dataset **inside the same history** using `workflows.invoke_workflow()` .
> > ### {% icon solution %} Solution
> > ```python
> > inputs = {0: {'id': hda['id'], 'src': 'hda'}}
> > ret = gi.workflows.invoke_workflow(wf['id'], inputs=inputs, history_id=new_hist['id'])
> > pprint(ret)
> > ```
> {: .solution}
{: .question}

```python
# Try it out here!

```

> ### {% icon question %} Question: View the results
> View the results on the Galaxy server with your web browser. Were you successful? Did it run?
{: .question}

# Interacting with histories in BioBlend.objects

**You need to insert the API key for your Galaxy server in the cell below**:

1. Open the Galaxy server in another browser tab
2. Click on "User" on the top menu, then "Preferences"
3. Click on "Manage API key"
4. Generate an API key if needed, then copy the alphanumeric string and paste it as the value of the `api_key` variable below.

The user interacts with a Galaxy server through a `GalaxyInstance` object:

```python
from pprint import pprint

import bioblend.galaxy.objects

server = 'https://usegalaxy.eu/'
api_key = ''
gi = bioblend.galaxy.objects.GalaxyInstance(url=server, api_key=api_key)
```

All `GalaxyInstance` method calls have the `client.method()` format, where `client` is the name of the resources you dealing with. There are 2 methods to get the list of resources:

- `get_previews()`: lightweight (one GET request), retrieves basic resources' info, returns a list of **preview** objects
- `list()`: one GET request for each resource, retrieves full resources' info, returns a list of **full** objects.

For example, the call to retrieve previews of all histories owned by the current user is:

```python
pprint(gi.histories.get_previews())
```

**New resources** are created with `create()` methods, e.g. to create a new history:

```python
new_hist = gi.histories.create(name='BioBlend test')
new_hist
```

As you can see, the `create()` methods in BioBlend.objects returns an object, not a dictionary.

Both `get_previews()` and `list()` methods usually have **filtering** capabilities, e.g. it is possible to filter histories **by name**:

```python
hists = gi.histories.list(name='BioBlend test')
pprint(hists)
```

Both `HistoryPreview` and `History` objects have many of their properties available as **attributes**, e.g. the id.

If you need to specify the unique **id** of the resource to retrieve, you can use the `get()` method, e.g. to get back the history we created before:

```python
hist0_id = hists[0].id
print(hist0_id)
h = gi.histories.get(hist0_id)
h
```

To **upload** files to the new history, run the special upload tool by calling the `upload_file` method of the `History` object:

```python
h.upload_file('test-data/1.txt')
```

To get the list of **datasets contained** in a history, simply look at the `content_infos` attribute of the `History` object.

```python
hdas = h.content_infos
pprint(hdas)
```

To get the details about one dataset, you can use the `get_dataset()` method of the `History` object:

```python
hda0_id = hdas[0].id
print(hda0_id)
h.get_dataset(hda0_id)
```

You can also filter history datasets by name using the `get_datasets()` method of `History` objects.


To **update** a resource, use the `update()` method of its object, e.g. to change the history name:

```python
new_hist.update(name='Updated history')
```

The return value of `update()` methods is the updated object.

Finally to **delete** a resource, you can use the `delete()` method of the object, e.g.:

```python
new_hist.delete()
```


## Exercise: BioBlend.objects

**Goal**: Upload a file to a new history, import a workflow and run it on the uploaded dataset.

> ### {% icon question %} Question: Initialise
> Create a `GalaxyInstance` object.
> > ### {% icon solution %} Solution
> > ```python
> > from pprint import pprint
> >
> > import bioblend.galaxy
> >
> > server = 'https://usegalaxy.eu/'
> > api_key = ''
> > gi = bioblend.galaxy.objects.GalaxyInstance(url=server, key=api_key)
> > ```
> >
> {: .solution}
{: .question}

```python
# Try it out here!

```

> ### {% icon question %} Question: New History
> Create a new Galaxy history.
> > ### {% icon solution %} Solution
> > ```python
> > new_hist = gi.histories.create(name='New history')
> > new_hist
> > ```
> {: .solution}
{: .question}

```python
# Try it out here!

```

> ### {% icon question %} Question: Upload a dataset
> **Upload** the local file "test-data/1.txt" to a new dataset in the created history using the `upload_file()` method of `History` objects.
>
> You can obtain this file from the following URL, you'll need to download it first.
>
> ```
> https://raw.githubusercontent.com/nsoranzo/bioblend-tutorial/main/test-data/1.txt
> ```
> > ### {% icon solution %} Solution
> > ```python
> > hda = new_hist.upload_file("test-data/1.txt")
> > hda
> > ```
> {: .solution}
{: .question}

```python
# Try it out here!

```

Please note that with BioBlend.objects there is no need to find the upload dataset, since `upload_file()` already returns a `HistoryDatasetAssociation` object.

> ### {% icon question %} Question: Import a workflow
> **Import a workflow** from the local file "test-data/convert_to_tab.ga" using `workflows.import_new()`
>
> You can obtain this file from the following URL, you'll need to download it first.
>
> ```
> https://raw.githubusercontent.com/nsoranzo/bioblend-tutorial/main/test-data/convert_to_tab.ga
> ```
> > ### {% icon solution %} Solution
> > ```python
> > with open('test-data/convert_to_tab.ga', 'r') as f:
> >     wf_string = f.read()
> > wf = gi.workflows.import_new(wf_string)
> > wf
> > ```
> {: .solution}
{: .question}

```python
# Try it out here!

```

> ### {% icon question %} Question: View the WF inputs
> > ### {% icon solution %} Solution
> > ```python
> > wf.inputs
> > ```
> {: .solution}
{: .question}

```python
# Try it out here!

```

> ### {% icon question %} Question: Run the WF
> **Run** the imported workflow on the uploaded dataset **inside the same history** using the `invoke()` method of `Workflow` objects.
> > ### {% icon solution %} Solution
> > ```python
> > inputs = {'0': {'id': hda.id, 'src': hda.SRC}}
> > wf.invoke(inputs=inputs, history=new_hist)
> > ```
> {: .solution}
{: .question}

```python
# Try it out here!

```

> ### {% icon question %} Question: View the results
> View the results on the Galaxy server with your web browser. Were you successful? Did it run?
{: .question}


# Optional Extra Exercises

If you have completed the exercise, you can try to perform these extra tasks with the help of the online documentation:

1. Download the workflow result to your computer
2. Publish your history
