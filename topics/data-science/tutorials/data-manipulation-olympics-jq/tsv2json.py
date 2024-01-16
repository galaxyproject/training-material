#!/usr/bin/env python
import json
import sys

path = sys.argv[1]
data = []
with open(path, 'r') as handle:
    header = next(handle).strip().split('\t')

    for row in handle:
        line = dict(zip(header, row.strip().split('\t')))
        for k in ['athlete_id', 'birth_year', 'height', 'weight', 'year', 'medal']:
            try:
                line[k] = int(line[k])
            except:
                if line[k] == "NA" or line[k] == "":
                    line[k] = None
                pass
        data.append(line)


print(json.dumps(data))
