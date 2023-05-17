#!/usr/bin/env python
import copy
import re
import sys
import logging

# read in a tutorial, and check the structure of it.
file = sys.argv[1]
if len(sys.argv) > 2 and sys.argv[2] == '--debug':
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.ERROR)

tuto = open(file, 'r')

boxes = r'^([\s>]*>[\s>]*)'
box_open = r'<([a-z-]*)-title>(.*)<\/[a-z-]*-title>'
box_close = r'{:\s*(\.[^}]*)\s*}'
whitespace = r'^(\s*)'
pre = r'^```'
in_pre = False

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
    logging.debug(f'A|{"[pre]" if in_pre else ""} depth={depth} line={line} {text:120s} tss={[x["tag"] for x in tag_stack]}')

    mw = re.match(whitespace, text).group(1)

    # Handling pre toggling
    unprefixed = stripN(text, depth)
    pre_toggle = re.match(pre, unprefixed)

    if not in_pre and depth > prev_depth:
        unprefixed = stripN(text, depth)
        m1 = re.match(box_open, unprefixed)
        if m1:
            tag = m1.group(1).replace('-', '_')
            logging.debug(f'B|{"[pre]" if in_pre else ""} line={line} {mw}{" " * (depth - 1)}<{tag}>')

            tag_stack.append({
                'tag': tag,
                'line': line,
                'text': text,
                'mw': mw,
            })
        else:
            logging.debug(f'C|{"[pre]" if in_pre else ""} LINE={line} {mw}{" " * (depth - 1)}<NONE>')
            tag_stack.append({
                'tag': None,
                'line': line,
                'text': text,
                'mw': mw,
            })
            logging.debug(f"D|{'[pre]' if in_pre else ''} {line} Potential Quote/Hidden/Spoken")
            # error?

    elif not in_pre and depth < prev_depth:
        unprefixed = stripN(text, depth)
        m1 = re.match(box_close, unprefixed.strip())
        logging.debug(f'E|{"[pre]" if in_pre else ""} prev={prev_depth} -> curr={depth} line={line} m1={m1} ({box_close} =~ {unprefixed}) ts={len(tag_stack)} tss={[x["tag"] for x in tag_stack]}')
        if m1 is None:
            message = f"Potential broken box. A {tag_stack[-1]['tag']} was opened on {tag_stack[-1]['line']}, but not closed on line {line}"
            print(f"{file}:{line}: {message}")
            logging.warning(f"{file}:{line}: {message}")
            logging.debug(f'F|{"[pre]" if in_pre else ""} NONE {line} |{text}|')
            logging.debug(tag_stack[-1])
            logging.debug(f"{mw}{' ' * (depth)}</NONE>")
            tag_stack.pop()
        else:
            if any([skippable_tag(BASE_SKIPPABLE_TAGS, t) for t in m1.groups()]):
                logging.debug(f"G|{'[pre]' if in_pre else ''} {mw}{' ' * (depth)}</NONE-skip tag={m1.groups()[0]}>")
                logging.debug(f'H|NONE {line} |{text}|')
                if ('code-2col' in m1.groups()[0] or 'hidden' in m1.groups()[0]) and (len(tag_stack) == 0 or tag_stack[-1]['tag'] is not None):
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
                logging.debug(f"I|{mw}{' ' * (depth)}</{closing_tag}>")

                if len(tag_stack) == 0:
                    message = f"Potential broken was closed with {closing_tag} on line {line}"
                    print(f"{file}:{line}: {message}")
                    logging.warning(f"{file}:{line}: {message}")
                if tag_stack[-1]['tag'] == closing_tag:
                    p = tag_stack.pop()
                else:
                    logging.debug(f'J|{"[pre]" if in_pre else ""} prev={prev_depth} -> curr={depth} line={line} m={m1} c={closing_tags}')
                    if not (tag_stack[-1]['tag'] is None and closing_tag not in BASE_EXPECTED_TAGS):
                        message = f"A {tag_stack[-1]['tag']} was opened, but closed with {closing_tag} on line {line}"
                        print(f"{file}:{line}: {message}")
                        logging.warning(f"{file}:{line}: {message}")
                    else:
                        p = tag_stack.pop()
    else:
        pass
        #unprefixed = stripN(text, depth)
        #pre_toggle = re.match(pre, unprefixed)
        #if pre_toggle:
        #    in_pre = not in_pre
        #    logging.debug(f'{"[pre]" if in_pre else ""} line={line} PRE {"OPEN" if in_pre else "CLOSE"} {text}')
    if pre_toggle:
        in_pre = not in_pre


    prev_depth = depth
