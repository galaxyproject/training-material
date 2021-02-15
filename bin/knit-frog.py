#!/usr/bin/env python
import re
import sys
import knittingneedles as knit

# read in a tutorial, and check the structure of it.
tuto = open(sys.argv[1], 'r')


diffs = []
current = None
for line, text in enumerate(tuto.read().split('\n')):
    m = re.match(knit.BOXES, text)
    if m:
        depth = m.group(1).count('>')
    else:
        depth = 0

    (unprefixed, prefix) = knit.stripN(text, depth)
    m1 = re.match(knit.BOX_OPEN, unprefixed)
    m2 = re.match(knit.BOX_CLOSE, unprefixed)
    # print(current is not None, '|', prefix, '|', text)

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
    commit_msg = knit.extractCommitMsg(diff)
    safe_commit = re.sub('[^a-z0-9-]', '-', commit_msg.lower())

    prefix = [
        'From: The Galaxy Training Network <gtn@example.org>',
        'Date: Mon, 15 Feb 2021 14:06:56 +0100',
        f'Subject: PATCH [{idx + 1}/{len(diffs)}]: {commit_msg}',
        '',
        '',
    ]

    (_, diff) = knit.removeWhitespacePrefix(diff)

    patch_id = diff[-1]
    # Remove patch id, ```
    diff = diff[0:-2]
    if diff[-1] == '{% endraw %}':
        diff = diff[0:-1]

    fn = f'commit-{idx:04d}-{safe_commit}.patch'
    with open(fn, 'w') as handle:
        print(fn)
        handle.write('\n'.join(prefix + diff + postfix))
