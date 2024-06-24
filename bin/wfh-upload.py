#!/usr/bin/env python
import yaml
import json
import os
import sys
import requests
import glob
import time
import multiprocessing

GTN_PROJECT_ID = 63

if len(sys.argv) == 2:
    crates = [sys.argv[1]]
else:
    crates = glob.glob(
        "_site/training-material/api/workflows/**/rocrate.zip", recursive=True
    )

def doUpload(crate_path):
    p = crate_path.split("/")
    (topic, tutorial, workflow) = p[3:6]

    payload = {
        "ro_crate": (crate_path, open(crate_path, "rb")),
        "workflow[project_ids][]": (None, GTN_PROJECT_ID), # GTN's ID.
    }
    headers = {
        "authorization": "Token " + os.environ["DEV_WFH_TOKEN"], 
        'User-Agent': 'GTN (github.com/galaxyproject/training-material@1.0)',
    }

    response = requests.post(
        "https://dev.workflowhub.eu/workflows/submit", files=payload, headers=headers
    )
    code = response.status_code
    if code != 200:
        print(f"Error {code} uploading {crate_path}")
        print(response.text)
        return None

    wfid = response.json()['data']['id']
    permissions_update = {
      "data": {
        "id": wfid,
        "type": "workflows",
        "attributes": {
          "policy": {
            "access": "download",
            "permissions": [
              {
                "resource": {
                  "id": str(GTN_PROJECT_ID),
                  "type": "projects"
                },
                "access": "manage"
              }
            ]
          }
        }
      }
    }
    headers.update({
        'Content-type': 'application/json',
        'Accept': 'application/json',
    })
    response2 = requests.put(f"https://dev.workflowhub.eu/workflows/{wfid}", headers=headers, json=permissions_update)
    if response2.status_code != 200:
        print(f"Error {response2.status_code} updating permissions for {wfid}: {response2.text}")

    return (topic, tutorial, workflow, wfid)


with multiprocessing.Pool(4) as p:
    results = p.map(doUpload, crates)

results = [x for x in results if x is not None]

data = {}
for (topic, tutorial, workflow, wf_id) in results:
    data['/'.join([topic, tutorial, workflow])] = wf_id

with open('metadata/workflowhub.yml', 'w') as handle:
    yaml.dump(data, handle)
