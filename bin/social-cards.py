import requests
import json
import os
import sys
import subprocess

GTN = 'https://training.galaxyproject.org/training-material'
# GTN = 'http://localhost:4000/training-material'

meta = requests.get(f"{GTN}/api/social-meta.json").json()

last_timestamp = int(sys.argv[1])

# meta looks like:
# {
# "/topics/proteomics/tutorials/protein-id-sg-ps/workflows/wf_proteinID_SG_PS.svg": 1710761924,
# "/topics/proteomics/tutorials/protein-id-sg-ps/workflows/wf_proteinID_SG_PS_multipleFiles.svg": 1710761924,
# "/topics/proteomics/tutorials/metaquantome-function/workflows/main_workflow.svg": 1710761924,
# }

count = 0
max_ts = 0
for path, time in meta.items():
    if time < last_timestamp:
        continue

    if '/tags/' in path:
        continue

    # Do something with the path
    print(f"New social card: {path}")
    # download the svg
    out = f'../gtn-social/{path}'
    os.makedirs(os.path.dirname(out), exist_ok=True)

    if not os.path.exists(out):
        subprocess.check_call(['wget', GTN + path, '-O', out])

    # Convert to png
    if not os.path.exists(out.replace('.svg', '.png')):
        subprocess.check_call(['convert', '-density', '100', out, out.replace('.svg', '.png')])
        subprocess.check_call(['optipng', out.replace('.svg', '.png')])
    max_ts = max(max_ts, time)

with open('../gtn-social/timestamp.txt', 'w') as f:
    f.write(str(max_ts))
