#!/usr/bin/env python
import os
import requests
import glob
import time

crates = glob.glob(
    "_site/training-material/api/workflows/**/rocrate.zip", recursive=True
)

for crate_path in crates:
    # _site/training-material/api/workflows/proteomics/proteogenomics-novel-peptide-analysis/galaxy-workflow-mouse_novel_peptide_analysis/rocrate.zip
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
        print(response.text)
        continue

    data = response.json()
    wf_id = data["data"]["id"]
    print(f"Uploaded {crate_path} as workflow {wf_id}")
    # time.sleep(1)
