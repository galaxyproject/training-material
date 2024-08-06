#!/usr/bin/env python
import yaml
import json
import os
import sys
import requests
import glob
import time
import multiprocessing
import argparse

# GTN_BOT_ID = https://workflowhub.eu/people/731

argparser = argparse.ArgumentParser()
argparser.add_argument("--prod", action="store_true", help="Upload to production")
argparser.add_argument("--crate", type=str, help="Upload a single crate (defaults to all)")
args = argparser.parse_args()

if args.prod:
    WORKFLOWHUB = "https://workflowhub.eu"
    GTN_PROJECT_ID = 12
else:
    WORKFLOWHUB = "https://dev.workflowhub.eu"
    GTN_PROJECT_ID = 63

if args.crate:
    crates = [args.crate]
else:
    crates = glob.glob(
        "_site/training-material/api/workflows/**/rocrate.zip", recursive=True
    )

def doUploadWrap(crate_path):
    try:
        return doUpload(crate_path)
    except (requests.exceptions.SSLError, requests.exceptions.ConnectionError) as e:
        try:
            time.sleep(5)
            return doUpload(crate_path)
        except (requests.exceptions.SSLError, requests.exceptions.ConnectionError) as e:
            print(f"Error uploading {crate_path}: {e}")
            return None

def doUpload(crate_path):
    payload = {
        "ro_crate": (crate_path, open(crate_path, "rb")),
        "workflow[project_ids][]": (None, GTN_PROJECT_ID), # GTN's ID.
    }
    headers = {
        "authorization": "Token " + os.environ["WFH_TOKEN"],
        'User-Agent': 'GTN (github.com/galaxyproject/training-material@1.0)',
        # 'Content-type': 'application/json',
        'Accept': 'application/json',
    }

    response = requests.post(
        f"{WORKFLOWHUB}/workflows/submit", files=payload, headers=headers
    )
    code = response.status_code
    if code != 200:
        print(f"Error {code} uploading {crate_path}")
        print(response.text)
        return None

    add_discussion_channel = False

    wfid = response.json()['data']['id']
    print(f"Uploaded {crate_path} as {wfid}")
    permissions_update = {
      "data": {
        "id": wfid,
        "type": "workflows",
        "attributes": {
        }
      }
    }
    current_policy = response.json()['data']['attributes']['policy']

    updated_policy = {}
    if current_policy['access'] != 'download':
        updated_policy['access'] = 'download'

    gtn_permission = [x for x in current_policy['permissions'] if x['resource']['id'] == str(GTN_PROJECT_ID)]
    if len(gtn_permission) != 1:
        updated_policy['permissions'] = current_policy['permissions'] + [
            {
                "resource": {
                    "id": str(GTN_PROJECT_ID),
                    "type": "projects"
                },
                "access": "manage"
            }
        ]
    else:
        if gtn_permission[0]['access'] != 'manage':
            gtn_permission[0]['access'] = 'manage'
            updated_policy['permissions'] = current_policy['permissions']

    # "policy": {
    #   "access": "download",
    #   "permissions": [
    #     {
    #       "resource": {
    #         "id": str(GTN_PROJECT_ID),
    #         "type": "projects"
    #       },
    #       "access": "manage"
    #     }
    #   ]
    # }

    # {'access': 'download',
    #  'permissions': [{'access': 'manage',
    #                   'resource': {'id': '63', 'type': 'projects'}}]}
    #
    push = False
    if updated_policy:
        push = True
        permissions_update['data']['attributes']['policy'] = updated_policy

    dls = response.json()['data']['attributes']['discussion_links']
    if dls is None or len(dls) == 0 or (not any(
        x['label'] == 'GTN Matrix' for x in response.json()['data']['attributes']['discussion_links']
    )):
        push = True
        permissions_update['data']['attributes']['discussion_links'] = [
            {
                "label": "GTN Matrix",
                "url": "https://matrix.to/#/%23galaxyproject_admin:gitter.im" # Must encode the second #
            }
        ]

    # https://github.com/seek4science/seek/issues/1957
    if push:
        print(f"Permissions update required for {wfid}: {permissions_update['data']['attributes']}")
        headers.update({
            'Content-type': 'application/json',
            'Accept': 'application/json',
        })
        response2 = requests.put(f"{WORKFLOWHUB}/workflows/{wfid}", headers=headers, json=permissions_update)
        if response2.status_code != 200:
            print(f"Error {response2.status_code} updating permissions for {wfid}: {response2.text}")

    p = crate_path.split("/")
    windex = p.index("workflows")
    (topic, tutorial, workflow) = p[windex + 1:windex + 4]
    return (topic, tutorial, workflow, wfid)



with multiprocessing.Pool(4) as p:
    results = p.map(doUploadWrap, crates)

results = [x for x in results if x is not None]

if os.path.exists('metadata/workflowhub.yml'):
    with open('metadata/workflowhub.yml') as handle:
        data = yaml.safe_load(handle)
else:
    data = {}

for (topic, tutorial, workflow, wf_id) in results:
    data['/'.join([topic, tutorial, workflow])] = wf_id

with open('metadata/workflowhub.yml', 'w') as handle:
    yaml.dump(data, handle)
