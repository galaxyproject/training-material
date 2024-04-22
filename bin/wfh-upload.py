#!/usr/bin/env python
import os
import sys
import requests
import glob
import time
import multiprocessing

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
    headers = {"authorization": "Token " + os.environ["DEV_WFH_TOKEN"]}

    response = requests.post(
        "https://dev.workflowhub.eu/workflows/submit", files=payload, headers=headers
    )
    code = response.status_code
    if code != 200:
        print(f"Error {code} uploading {crate_path}")
        # sys.stdout.write('!')
        return response.text
    # sys.stdout.write('.')

    data = response.json()
    wf_id = data["data"]["id"]
    return f"Uploaded {crate_path} as workflow {wf_id}"


with multiprocessing.Pool(8) as p:
    results = p.map(doUpload, crates)

print()
for result in results:
    print(result)
