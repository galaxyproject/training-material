#!/usr/bin/env python
import copy
import re
import sys

# read in a tutorial, and check the structure of it.
tuto = open(sys.argv[1], 'r')

boxes = r'^([\s>]*>[\s>]*)'
box_open = r'### {% icon (.*) %} (.*)'
box_close = r'{: *(\.[^} ]*)}'
whitespace = r'^(\s*)'

tag_stack = []
prev_depth = 0


def stripN(line, count):
    c = copy.copy(line)
    for i in range(count):
        c = c.lstrip()
        c = c.lstrip('>')
        c = c.lstrip()
    return c




for line, text in enumerate(tuto.read().split('\n')):

    m = re.match(boxes, text)
    if m:
        depth = m.group(1).count('>')
    else:
        depth = 0

    mw = re.match(whitespace, text).group(1)

    if depth > prev_depth:
        unprefixed = stripN(text, depth)
        m1 = re.match(box_open, unprefixed)
        if m1:
            tag = m1.group(1)
            # print(f'{mw}{">" * depth} {tag}')

            tag_stack.append({
                'tag': tag,
                'line': line,
                'text': text,
                'mw': mw,
            })
        else:
            pass
            # error?

    elif depth < prev_depth:
        unprefixed = stripN(text, depth)
        m1 = re.match(box_close, unprefixed)
        print(f'd={depth} prev={prev_depth} line={line} text={text} m={m1}')



    prev_depth = depth
