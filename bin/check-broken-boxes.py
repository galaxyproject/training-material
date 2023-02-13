#!/usr/bin/env python
import copy
import re
import sys

# read in a tutorial, and check the structure of it.
file = sys.argv[1]
tuto = open(file, 'r')

boxes = r'^([\s>]*>[\s>]*)'
box_open = r'<([a-z-]*)-title>(.*)<\/[a-z-]*-title>'
box_close = r'{:\s*(\.[^}]*)\s*}'
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

def skippable_tag(tags, check):
    for tag in tags:
        if tag in check:
            return True
    return False

BASE_SKIPPABLE_TAGS = ['hidden', 'quote', 'spoken', 'code-2col']
BASE_EXPECTED_TAGS = [
    'hands_on',
    'tip',
    'details'
]

for line, text in enumerate(tuto.read().split('\n')):

    m = re.match(boxes, text)
    if m:
        depth = m.group(1).count('>')
    else:
        depth = 0
    # print(f'{depth} {line} {text}')

    mw = re.match(whitespace, text).group(1)

    if depth > prev_depth:
        unprefixed = stripN(text, depth)
        m1 = re.match(box_open, unprefixed)
        if m1:
            tag = m1.group(1).replace('-', '_')
            # print(f'{mw}{" " * (depth - 1)}<{tag}>')

            tag_stack.append({
                'tag': tag,
                'line': line,
                'text': text,
                'mw': mw,
            })
        else:
            # print(f'{mw}{" " * (depth - 1)}<NONE>')
            tag_stack.append({
                'tag': None,
                'line': line,
                'text': text,
                'mw': mw,
            })
            # print(f"{line} Potential Quote/Hidden/Spoken")
            # error?

    elif depth < prev_depth:
        unprefixed = stripN(text, depth)
        m1 = re.match(box_close, unprefixed)
        # print(f'prev={prev_depth} -> curr={depth} line={line} m={m1} ts={len(tag_stack)} tss={[x["tag"] for x in tag_stack]}')
        if m1 is None:
            # print(f'NONE {line} |{text}|')
            # print(tag_stack[-1])
            # print(f"{mw}{' ' * (depth)}</NONE>")
            tag_stack.pop()
        else:
            if any([skippable_tag(BASE_SKIPPABLE_TAGS, t) for t in m1.groups()]):
                # print(f"{mw}{' ' * (depth)}</NONE-skip tag={m1.groups()[0]}>")
                # print(f'NONE {line} |{text}|')
                # print(tag_stack[-1])
                if 'code-2col' in m1.groups()[0] and (len(tag_stack) == 0 or tag_stack[-1]['tag'] is not None):
                    pass
                    # This is a special case.
                    # Here we don't have tag_stack[-1]['tag'] is None
                    # Because there wasn't even a blank header line before the 2col started.
                else:
                    tag_stack.pop()
            else:
                closing_tags = m1.groups(1)[0].replace('-', '_').lstrip('.').split('.')
                closing_tag = closing_tags[0].strip()
                # print(tag_stack[-1])
                # print(f"{mw}{' ' * (depth)}</{closing_tag}>")

                if len(tag_stack) == 0:
                    message = f"Potential broken was closed with {closing_tag} on line {line}"
                    print(f"{file}:{line}: {message}")
                if tag_stack[-1]['tag'] == closing_tag:
                    p = tag_stack.pop()
                else:
                    # print(f'prev={prev_depth} -> curr={depth} line={line} m={m1} c={closing_tags}')
                    if not (tag_stack[-1]['tag'] is None and closing_tag not in BASE_EXPECTED_TAGS):
                        message = f"A {tag_stack[-1]['tag']} was opened, but closed with {closing_tag} on line {line}"
                        print(f"{file}:{line}: {message}")
                    else:
                        p = tag_stack.pop()

    prev_depth = depth
