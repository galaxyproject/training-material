#!/bin/bash
# This script is used to postprocess the video files to add chapter markers based on the input file
import json
import argparse

argparser = argparse.ArgumentParser(description='Postprocess video files to add chapter markers')
argparser.add_argument('input', help='Input cast file', type=argparse.FileType('r'))
argparser.add_argument('output', help='Save to this file', type=argparse.FileType('w'))
args = argparser.parse_args()
import re

def load():
    for line in args.input.readlines():
        if 'Next step: ' in line:
            parsed_line = json.loads(line)
            # This is a chapter start
            m = re.findall(r'Next step: [a-z0-9/_-]*: ([A-Za-z-0-9 _-]*)', line)
            if len(m) > 0:
                yield json.dumps([
                    parsed_line[0],
                    "m",
                    m[0]
                ]) + '\n'
            yield line
        else:
            yield line

# line = '[0.269871, "o", "\u001b[1mNext step: dev/tool-from-scratch/0000: Planemo init\u001b(B\u001b[m\r\n"]'
# print(re.findall(r'Next step: [a-z0-9_/-]*: ([A-Za-z-0-9 _-]*)', line))

for line in load():
    args.output.write(line)
