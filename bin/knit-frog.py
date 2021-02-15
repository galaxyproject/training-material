#!/usr/bin/env python
import copy
import re
import sys

# read in a tutorial, and check the structure of it.
tuto = open(sys.argv[1], 'r')

boxes = r'^([\s>]*>[\s>]*)'
box_open = r'\s*```diff'
box_close = r'\s*{: data-commit="([^"]*)"}'
whitespace = r'^(\s*)'

diffs = []

def stripN(line, count):
    c = copy.copy(line)
    for i in range(count):
        c = c.lstrip()
        c = c.lstrip('>')
    return c

current = None
for line, text in enumerate(tuto.read().split('\n')):
    m = re.match(boxes, text)
    if m:
        depth = m.group(1).count('>')
    else:
        depth = 0

    unprefixed = stripN(text, depth)
    m1 = re.match(box_open, unprefixed)
    m2 = re.match(box_close, unprefixed)
    # print(current is not None, text)

    if m1:
        if current is not None:
            diffs.append(current)
        current = []
    elif current is not None:
        current.append(unprefixed)

    if current and len(current) > 2 and current[-2].strip() == '```' and not m2:
        current = None

    if m2:
        diffs.append(current)
        current = None

postfix = [
    '--',
    '2.25.1',
    '',
    ''
]

# import sys; sys.exit()


for idx, diff in enumerate(diffs):
    commit_msg = re.match(box_close, diff[-1].strip()).group(1)
    safe_commit = re.sub('[^a-z0-9-]', '-', commit_msg.lower())
    prefix = [
        'From: The Galaxy Training Network <gtn@example.org>',
        'Date: Mon, 15 Feb 2021 14:06:56 +0100',
        f'Subject: PATCH [{idx + 1}/{len(diffs)}]: {commit_msg}',
        '',
        '',
    ]

    # Fix whitespace.
    whitespace_prefix = [
        len(re.match(whitespace, line).group(1))
        for line in diff
    ]
    # __import__('pprint').pprint(whitespace_prefix)
    # Remove zeros (blank lines have no whitespace, not even the standard of
    # the rest)
    whitespace_prefix = [x for x in whitespace_prefix if x != 0]
    # print(whitespace_prefix)
    shortest = min(whitespace_prefix)

    diff = [
        x[shortest:]
        for x in diff
    ]

    patch_id = diff[-1]
    # Remove patch id, ```
    diff = diff[0:-2]
    if diff[-1] == '{% endraw %}':
        diff = diff[0:-1]

    fn = f'commit-{idx:04d}-{safe_commit}.patch'
    with open(fn, 'w') as handle:
        print(fn)
        handle.write('\n'.join(prefix + diff + postfix))
