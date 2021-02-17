#!/usr/bin/env python
import argparse
import re
import sys
import knittingneedles as knit


# read in a tutorial, and check the structure of it.
parser = argparse.ArgumentParser(
    description="Extract specially formatted diffs from tutorials as a series of patches to apply"
)
parser.add_argument("tutorial", type=argparse.FileType("r"), help="Input tutorial")
parser.add_argument(
    "prefix", help="Some prefix for the output commit files to keep them separate"
)
args = parser.parse_args()


fnparts = args.tutorial.name.split('/')
fn_topic = fnparts[fnparts.index('tutorials') - 1]
fn_tutorial = fnparts[fnparts.index('tutorials') + 1]

diffs = []
current = None
for line, text in enumerate(args.tutorial.read().split("\n")):
    m = re.match(knit.BOXES, text)
    if m:
        depth = m.group(1).count(">")
    else:
        depth = 0

    (unprefixed, prefix) = knit.stripN(text, depth)
    m1 = re.match(knit.BOX_OPEN, unprefixed)
    m2 = re.match(knit.BOX_CLOSE, unprefixed)

    if m1:
        if current is not None:
            diffs.append(current)
        current = []
    elif current is not None:
        current.append(unprefixed)

    if current and len(current) > 2 and current[-2].strip() == "```" and not m2:
        current = None

    if m2:
        diffs.append(current)
        current = None

postfix = ["--", "2.25.1", "", ""]

# import sys; sys.exit()


for idx, diff in enumerate(diffs):
    commit_msg = knit.extractCommitMsg(diff)
    safe_commit = re.sub("[^a-z0-9-]", "-", commit_msg.lower())

    prefix = [
        "From: The Galaxy Training Network <gtn@example.org>",
        "Date: Mon, 15 Feb 2021 14:06:56 +0100",
        f"Subject: {fn_topic}/{fn_tutorial}/{idx:04d}: {commit_msg}",
        "",
        "",
    ]

    (_, diff) = knit.removeWhitespacePrefix(diff)

    patch_id = diff[-1]
    # Remove patch id, ```
    diff = diff[0:-2]
    if diff[-1] == "{% endraw %}":
        diff = diff[0:-1]

    fn = f"{args.prefix}-commit-{idx:04d}-{safe_commit}.patch"
    with open(fn, "w") as handle:
        print(fn)
        handle.write("\n".join(prefix + diff + postfix))
