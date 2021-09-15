#!/usr/bin/env python
import argparse
import re
import sys
import knittingneedles as knit


parser = argparse.ArgumentParser(
    description="Extract specially formatted diffs from tutorials as a series of patches to apply"
)
parser.add_argument("tutorial", type=argparse.FileType("r+"), help="Input tutorial")
parser.add_argument(
    "--patches", nargs='+', help="The patches to knit together with the tutorial"
)
args = parser.parse_args()
tutorial_contents = args.tutorial.read().split("\n")
chunks = []


# load patches
diffs = []
for patch in args.patches:
    with open(patch, "r") as handle:
        data = handle.readlines()
        subject = [x for x in data[0:5] if x.startswith('Subject:')]
        if len(subject) != 1:
            commit = "MISSING"
        else:
            if subject[0].count(": ") >= 2:
                commit = subject[0].split(": ")[2]
            else:
                commit = "MISSING"

        indexindex = [idx for (idx, x) in enumerate(data) if x.startswith('index ')]
        if len(indexindex) == 0:
            # take a guess, it's our format?
            data = [x.rstrip("\n") for x in data[5:-3]]
        else:
            # If 'index ' then we start at the next line
            data = [x.rstrip("\n") for x in data[indexindex[0] + 1:-3]]

        diffs.append((commit, data))

current = None
diff_idx = 0
prev_line = None
for line, text in enumerate(tutorial_contents):
    m = re.match(knit.BOXES, text)
    if m:
        depth = m.group(1).count(">")
    else:
        depth = 0
    (unprefixed, prefix) = knit.stripN(text, depth)
    m0 = re.match(knit.BOX_PREFIX, unprefixed)
    m1 = re.match(knit.BOX_OPEN, unprefixed)
    m2 = re.match(knit.BOX_CLOSE, unprefixed)

    if m0:
        current = []
    if m1 and current is None:
        current = []

    if current is not None:
        current.append(unprefixed)
    else:
        chunks.append(text)

    if current and len(current) > 2 and (current[-2].strip() == "```") and not m2:
        chunks.extend([prefix + x for x in current])
        # chunks.append('REJECT)')
        current = None

    if m2 and current and len(current) > 0:
        (amount, _) = knit.removeWhitespacePrefix(current)
        prefix_text = prefix + (amount * " ")
        from_patch = diffs[diff_idx]

        # Compare the diff messages
        obs_msg = knit.extractCommitMsg(current).strip()
        exp_msg = from_patch[0].strip()
        obs_extra = m2.groups()[1]
        if obs_msg != exp_msg:
            # Most just got truncated
            if exp_msg not in obs_msg:
                print(f"WARNING: Diff messages do NOT line up: {obs_msg} != {exp_msg}")
            # Replace with current message
            exp_msg = obs_msg

        # Knit it anyway
        new_diff = from_patch[1]
        chunks.append(prefix_text + "{% raw %}")
        chunks.append(prefix_text + "```diff")
        chunks.extend([f'{prefix}{amount * " "}{line}' for line in new_diff])
        chunks.append(prefix_text + "{% endraw %}")
        chunks.append(prefix_text + "```")
        chunks.append(prefix_text + '{: data-commit="%s"%s}' % (exp_msg, obs_extra))
        diff_idx += 1
        current = None

    prev_line = text

# Overwrite the original tutorial
args.tutorial.seek(0)
args.tutorial.truncate(0)
# The last chunk is a newline for some reason.
for c in chunks[0:-1]:
    args.tutorial.write(c + "\n")
args.tutorial.close()
