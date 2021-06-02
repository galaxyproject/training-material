#!/usr/bin/env python
import shutil
import os
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
currentCmd = None
currentTest = None
for line, text in enumerate(args.tutorial.read().split("\n")):
    m = re.match(knit.BOXES, text)
    if m:
        depth = m.group(1).count(">")
    else:
        depth = 0

    (unprefixed, prefix) = knit.stripN(text, depth)
    m1 = re.match(knit.BOX_OPEN, unprefixed)
    m2 = re.match(knit.BOX_CLOSE, unprefixed)
    c1 = re.match(knit.CMD_OPEN, unprefixed)
    c2 = re.match(knit.CMD_CLOSE, unprefixed)
    t1 = re.match(knit.TEST_OPEN, unprefixed)
    t2 = re.match(knit.TEST_CLOSE, unprefixed)

    if m1:
        if current is not None:
            diffs.append(current)
        current = []
    elif current is not None:
        current.append(unprefixed)

    if c1:
        if currentCmd is not None:
            diffs.append(currentCmd)
        currentCmd = []
    elif currentCmd is not None:
        currentCmd.append(unprefixed)

    if t1:
        if currentTest is not None:
            diffs.append(currentTest)
        currentTest = []
    elif currentTest is not None:
        currentTest.append(unprefixed)

    if current and len(current) > 2 and current[-2].strip() == "```" and not m2:
        current = None

    if currentCmd and len(currentCmd) > 2 and currentCmd[-2].strip() == "```" and not c2:
        currentCmd = None

    if currentTest and len(currentTest) > 2 and currentTest[-2].strip() == "```" and not t2:
        currentTest = None

    # Tests are known short, same with commands. If we haven't exited by now, it's an issue..
    if currentTest and len(currentTest) > 4:
        currentTest = None

    # interesting = (t1, t2, currentCmd, currentTest)
    # if any([x is not None for x in interesting]):
        # print(*interesting)

    if m2:
        diffs.append(current)
        current = None

    if c2:
        diffs.append(currentCmd)
        currentCmd = None

    if t2:
        diffs.append(currentTest)
        currentTest = None

postfix = ["--", "2.25.1", "", ""]

# import sys; sys.exit()


GITGAT = os.path.dirname(args.prefix)
BASE = os.path.basename(args.prefix)

cmdhandle = open(f"{GITGAT}/.scripts/{BASE}-run.sh", 'w')
cmdhandle.write("#!/bin/bash\n")
cmdhandle.write("set -e\n\n")

lastCommit = None
for idx, diff in enumerate(diffs):
    if 'data-commit' in diff[-1]:
        commit_msg = knit.extractCommitMsg(diff)
        safe_commit = re.sub("[^a-z0-9-]", "-", commit_msg.lower())

        prefix = [
            "From: The Galaxy Training Network <galaxytrainingnetwork@gmail.com>",
            "Date: Mon, 15 Feb 2021 14:06:56 +0100",
            f"Subject: {fn_topic}/{fn_tutorial}/{idx:04d}: {commit_msg}",
            "",
            "",
        ]

        lastCommit = f"{fn_topic}/{fn_tutorial}/{idx:04d}"

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
    elif 'data-cmd' in diff[-1]:
        cmdhandle.write("\n# CMD\n")
        if lastCommit is not None:
            cmdhandle.write(f'## Checkout\ngit checkout $(git log main --pretty=oneline | grep "{lastCommit}" | cut -c1-40)\n')
        for line in diff[0:-2]:
            cmdhandle.write("## Run command\n")
            if 'openssl rand' in line and 'vault-password' in line:
                cmdhandle.write("## The students should use a random password, we override with 'password' for reproducibility\n")
                cmdhandle.write("echo 'password' > .vault-password.txt;\n")
            else:
                cmdhandle.write(line.strip() + "\n")
    elif 'data-test' in diff[-1]:
        cmdhandle.write("\n# TEST\n")
        if lastCommit is not None:
            cmdhandle.write(f'## Checkout\ngit checkout $(git log main --pretty=oneline | grep "{lastCommit}" | cut -c1-40)\n')
        for line in diff[0:-2]:
            testdir = f"topics/{fn_topic}/tutorials/{fn_tutorial}/tests"
            gittestdir = f"{GITGAT}/.scripts/{BASE}-test"
            if os.path.exists(testdir) and not os.path.exists(gittestdir):
                shutil.copytree(testdir, gittestdir)

            cmdhandle.write("## Install test dependencies\n")
            if os.path.exists(f"{gittestdir}/requirements.txt"):
                cmdhandle.write(f"pip install -r {gittestdir}/requirements.txt\n")

            cmdhandle.write("## Run test case\n")
            cmdhandle.write(f"bash {gittestdir}/{line.strip()}\n")
    else:
        print("Unknown!")
