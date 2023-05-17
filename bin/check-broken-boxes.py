#!/usr/bin/env python
import copy
import re
import sys
import logging

# read in a tutorial, and check the structure of it.
file = sys.argv[1]
if len(sys.argv) > 2 and sys.argv[2] == '--debug':
    logging.basicConfig(level=logging.DEBUG)

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
    logging.debug(f'{depth} {line} {text}')

    mw = re.match(whitespace, text).group(1)

    if depth > prev_depth:
        unprefixed = stripN(text, depth)
        m1 = re.match(box_open, unprefixed)
        if m1:
            tag = m1.group(1).replace('-', '_')
            logging.debug(f'line={line + 1} {mw}{" " * (depth - 1)}<{tag}>')

            tag_stack.append({
                'tag': tag,
                'line': line,
                'text': text,
                'mw': mw,
            })
        else:
            logging.debug(f'{mw}{" " * (depth - 1)}<NONE>')
            tag_stack.append({
                'tag': None,
                'line': line,
                'text': text,
                'mw': mw,
            })
            logging.debug(f"{line} Potential Quote/Hidden/Spoken")
            # error?

    elif depth < prev_depth:
        unprefixed = stripN(text, depth)
        m1 = re.match(box_close, unprefixed)
        logging.debug(f'prev={prev_depth} -> curr={depth} line={line + 1} m1={m1} ts={len(tag_stack)} tss={[x["tag"] for x in tag_stack]}')
        if m1 is None:
            message = f"Potential broken box. A {tag_stack[-1]['tag']} was opened, but not closed on line {line}"
            print(f"{file}:{line}: {message}")
            logging.debug(f'NONE {line} |{text}|')
            logging.debug(tag_stack[-1])
            logging.debug(f"{mw}{' ' * (depth)}</NONE>")
            tag_stack.pop()
        else:
            if any([skippable_tag(BASE_SKIPPABLE_TAGS, t) for t in m1.groups()]):
                logging.debug(f"{mw}{' ' * (depth)}</NONE-skip tag={m1.groups()[0]}>")
                logging.debug(f'NONE {line} |{text}|')
                logging.debug(tag_stack[-1])
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
                logging.debug(tag_stack[-1])
                logging.debug(f"{mw}{' ' * (depth)}</{closing_tag}>")

                if len(tag_stack) == 0:
                    message = f"Potential broken was closed with {closing_tag} on line {line}"
                    print(f"{file}:{line}: {message}")
                if tag_stack[-1]['tag'] == closing_tag:
                    p = tag_stack.pop()
                else:
                    logging.debug(f'prev={prev_depth} -> curr={depth} line={line} m={m1} c={closing_tags}')
                    if not (tag_stack[-1]['tag'] is None and closing_tag not in BASE_EXPECTED_TAGS):
                        message = f"A {tag_stack[-1]['tag']} was opened, but closed with {closing_tag} on line {line}"
                        print(f"{file}:{line}: {message}")
                    else:
                        p = tag_stack.pop()

    prev_depth = depth
