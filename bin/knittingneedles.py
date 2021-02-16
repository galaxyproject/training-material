import copy
import re

BOXES = r"^([\s>]*>[\s>]*)"
BOX_PREFIX = r"\s*{% raw %}"
BOX_OPEN = r"\s*```diff"
BOX_CLOSE = r'\s*{: data-commit="([^"]*)"}'
WHITESPACE = r"^(\s*)"


def stripN(line, count):
    c = copy.copy(line)
    for i in range(count):
        c = c.lstrip()
        c = c.lstrip(">")
    removed = len(line) - len(c)
    return (c, line[0:removed])


def removeWhitespacePrefix(lines):
    # Obtain whitespace amount
    whitespace_prefix = [len(re.match(WHITESPACE, line).group(1)) for line in lines]
    # Remove zeros (blank lines have no whitespace, not even the standard of
    # the rest)
    whitespace_prefix = [x for x in whitespace_prefix if x != 0]
    amount = min(whitespace_prefix)
    diff = [x[amount:] for x in lines]
    return (amount, diff)


def extractCommitMsg(lines):
    return re.match(BOX_CLOSE, lines[-1].strip()).group(1)
