#!/usr/bin/env python
import yaml
import json
import os
import sys
import requests
import glob
import time
import multiprocessing

if len(sys.argv) == 2:
    crates = [sys.argv[1]]
else:
    crates = glob.glob(
        "_site/training-material/api/workflows/**/rocrate.zip", recursive=True
    )

def doUpload(crate_path):
    p = crate_path.split("/")
    (topic, tutorial, workflow) = p[4:7]

    payload = {
        "ro_crate": (crate_path, open(crate_path, "rb")),
        "workflow[project_ids][]": (None, 63),
    }
    headers = {"authorization": "Token " + os.environ["DEV_WFH_TOKEN"], 'User-Agent': 'GTN (github.com/galaxyproject/training-material@1.0)'}

    response = requests.post(
        "https://dev.workflowhub.eu/workflows/submit", files=payload, headers=headers
    )
    code = response.status_code
    if code != 200:
        print(f"Error {code} uploading {crate_path}")
        print(response.text)
        return None
    else:
        print(json.loads(response.text)['data']['links']['self'])

    data = response.json()
    wf_id = data["data"]["id"]
    return (topic, tutorial, workflow, wf_id)


with multiprocessing.Pool(4) as p:
    results = p.map(doUpload, crates)

results = [x for x in results if x is not None]

data = {}
for (topic, tutorial, workflow, wf_id) in results:
    data.setdefault(topic, {}).setdefault(tutorial, {})[workflow] = wf_id

with open('metadata/workflowhub.yml', 'w') as handle:
    yaml.dump(data, handle)
