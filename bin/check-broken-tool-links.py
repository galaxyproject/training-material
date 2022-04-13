#!/usr/bin/env python3
import sys
import re

ALLOWED_SHORT_IDS = [
    'intermine',
    'Count1',
    'Cut1',
    'Filter1',
    'Grep1',
    'Grouping1',
    'Paste1',
    'Remove beginning1',
    'Show beginning1',
    'Summary_Statistics1',
    'addValue',
    'cat1',
    'comp1',
    'gene2exon1',
    'join1',
    'random_lines1',
    'ucsc_table_direct1',
    'upload1',
    'sort1',
    'wig_to_bigWig',
    'ChangeCase',
    'param_value_from_file',
    'Convert characters1'
]

RE_TOOL = '{% tool ([^[%}]*)\s*%}'
RE_TOOL_ID = '{% tool \[([^]]*)\]\(([^%]*)\)\s*%}'
RE_LINE = '(topics[^:]*):([0-9]+):(.*)'

WROTE = 0

for line in sys.stdin.readlines():
    fn, lineNo, lineContents = re.match(RE_LINE, line).groups()

    # There isn't anything to check here.
    # for m in re.findall(RE_TOOL, lineContents):
        # pass

    for m in re.finditer(RE_TOOL_ID, lineContents):
        (title, toolid) = m.groups()
        # TS tool
        if '/' in toolid:
            if toolid.count('/') < 5:
                print(f'::error file="{fn}",line="{lineNo}",col={m.start()}::Broken tool link, {toolid} has fewer than the expected number of slashes, meaning you have probably left off the version component.')
                WROTE = 1

            if 'testtoolshed' in toolid:
                print(f'::warning file="{fn}",line="{lineNo}",col={m.start()}::The GTN strongly avoids using testtoolshed tools in your tutorials or workflows')
                WROTE = 1
        else:
            if '+' in toolid:
                print(f'::error file="{fn}",line="{lineNo}",col={m.start()}::Broken tool link, {toolid} has an unnecessary +')
                WROTE = 1

            if toolid not in ALLOWED_SHORT_IDS and not toolid.startswith('interactive_tool_') and not re.match('__[A-Z_]+__', toolid):
                print(f'::error file="{fn}",line="{lineNo}",col={m.start()}::Unknown short tool ID {toolid}. If this is correct, please add to bin/check-broken-tool-links.py')
                WROTE = 1

        if '://' in toolid or '?' in toolid or '#' in toolid:
            print(f'::error file="{fn}",line="{lineNo}",col={m.start()}::Broken tool link, {toolid} seems to include part of a URL which is incorrect. It should just be the tool ID.')
            WROTE = 1

sys.exit(WROTE)
