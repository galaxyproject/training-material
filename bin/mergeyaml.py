#!/usr/bin/env python

import glob
import yaml

merged = {'libraries': list()}
merged_libs = dict()

def is_url_included(url, lib):
    for files in lib.values():
        #print f
        for _url, _file_type in files:
            #print x
            if _url == url:
                return True
    return False

for filename in glob.iglob('./**/data-library.yaml'):
    a = yaml.load(open(filename))
    #print a
    #print '##########'
    libs = a['libraries']
    for lib in libs:
        name = lib['name']
        if name in merged_libs:
            for f in lib['files']:
                if not is_url_included(f['url'], merged_libs):
                    merged_libs[name].append(f)
        else:
            merged_libs[name] = lib['files']

for name, files in merged_libs.items():
    merged['libraries'].extend([{'files': files, 'name': name}])

print yaml.dump( merged, default_flow_style=False, default_style='' )


