#!/usr/bin/env python

import glob
import yaml


def extend_dict(merged, a):
    if isinstance(merged, dict):
        for k, v in a.items():
            if k in merged:
                extend_dict(merged[k], v)
            else:
                merged[k] = v
    else:    
        if isinstance(merged, list):
            extend_list(merged, a)
        else:
            merged += a


def extend_list(merged, a):
    missing = []
    for itema in a:
        if not isinstance(itema, dict) or itema in missing:
            continue
        if 'name' in itema:
            match = next((i for i in merged if i["name"] == itema["name"]), False)
            if match:
                extend_dict(match, itema)
            else:
                missing += [itema, ]
    merged += missing


merged = {}

for filename in glob.iglob('./**/data-library.yaml'):
    a = yaml.load(open(filename))
    extend_dict(merged, a)


print(yaml.dump(merged, default_flow_style=False, default_style=''))
