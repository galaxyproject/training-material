#!/usr/bin/env python

import os
import glob

# get all Dockerfiles and compare the first part of the file paths with the changed files.
dockerfiles = glob.glob(r'./**/Dockerfile', recursive=True)

dockerdirs = set()

for line in open('./CHANGES.list'):
    line = line.strip()
    for dockerfile in dockerfiles:
        if os.path.dirname(dockerfile).replace('/docker', '')[2:] in line:
            dockerdirs.add( os.path.dirname(dockerfile) )

with open('DOCKER_BUILDS.list', 'w+') as handle:
    for entry in dockerdirs:
        print(entry)
        handle.write('%s\n' % entry)